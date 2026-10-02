import argparse
import os
import subprocess
import sys
from unittest.mock import patch

import pytest

from sage.cli import main
from sage.cli.options import CliOptions
from sage.cli.run_file_cmd import RunFileCmd, _run_file


@pytest.fixture(autouse=True)
def _preserve_process_state():
    """Keep direct ``RunFileCmd`` tests isolated from process-wide state."""
    import multiprocessing.spawn

    argv = sys.argv
    path = sys.path.copy()
    get_preparation_data = multiprocessing.spawn.get_preparation_data
    main_modules = {
        name: (name in sys.modules, sys.modules.get(name))
        for name in ('__main__', '__mp_main__')
    }
    yield
    sys.argv = argv
    sys.path[:] = path
    multiprocessing.spawn.get_preparation_data = get_preparation_data
    for name, (was_present, module) in main_modules.items():
        if was_present:
            sys.modules[name] = module
        else:
            sys.modules.pop(name, None)


def _parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--quiet", action="store_true", default=False)
    RunFileCmd.extend_parser(parser)
    return parser


def test_run_file_parser_forwards_arguments_after_file():
    args = _parser().parse_args(["-q", "test.sage", "42", "-y", "7", "--quiet"])

    assert args.quiet is True
    assert args.file == ["test.sage", "42", "-y", "7", "--quiet"]


def test_run_file_parser_handles_explicit_terminator():
    args = _parser().parse_args(["--", "-test.sage", "-y"])

    assert args.file == ["-test.sage", "-y"]


def test_main_forwards_script_arguments(capsys, tmp_path):
    """Regression test for #40871 and #41908."""
    file = tmp_path / "test.sage"
    file.write_text(
        "import argparse\n"
        "import sys\n"
        "parser = argparse.ArgumentParser()\n"
        "parser.add_argument('x', type=ZZ)\n"
        "parser.add_argument('-y', '--why', type=ZZ)\n"
        "parser.add_argument('-c', '--command')\n"
        "args = parser.parse_args()\n"
        "print(sys.argv[1:])\n"
        "print(args.x + args.why)\n"
        "print(args.command)\n"
    )

    with patch.object(
        sys, 'argv', ["sage", "-q", str(file), "42", "-y", "7", "-c", "script-c"]
    ):
        assert main() == 0

    captured = capsys.readouterr()
    assert captured.out == "['42', '-y', '7', '-c', 'script-c']\n49\nscript-c\n"


def test_run_file_cmd(capsys, tmp_path):
    file = tmp_path / "test.sage"
    file.write_text("print(3^33)")
    options = CliOptions(file=[str(file)])
    run_file_cmd = RunFileCmd(options)
    path_entry = sys.path[0]

    run_file_cmd.run()
    captured = capsys.readouterr()
    assert captured.out == "5559060566555523\n"
    assert sys.path[0] == path_entry


def test_run_file_cmd_with_args(capsys, tmp_path):
    with patch.object(sys, 'argv', ["python3", "test.sage", "1", "1"]):
        file = tmp_path / "test.sage"
        file.write_text("import sys; print(int(sys.argv[1]) + int(sys.argv[2]))")
        options = CliOptions(file=[str(file), "1", "1"])
        run_file_cmd = RunFileCmd(options)

        run_file_cmd.run()
        captured = capsys.readouterr()
        assert captured.out == "2\n"


def test_run_file_cmd_argv_matches_python(capsys, tmp_path):
    """``sys.argv`` inside a script must look like ``[script, *args]``."""
    with patch.object(sys, 'argv', ["python3", "test.sage", "-y", "7"]):
        file = tmp_path / "test.sage"
        file.write_text("import sys; print(sys.argv)")
        options = CliOptions(file=[str(file), "-y", "7"])
        RunFileCmd(options).run()
        captured = capsys.readouterr()
        assert captured.out == f"[{str(file)!r}, '-y', '7']\n"


def test_run_file_cmd_rebuilds_argv_for_each_run(capsys, tmp_path):
    """Construction has no side effect and scripts cannot poison a rerun."""
    file = tmp_path / "test.py"
    file.write_text("import sys\nprint(sys.argv)\nsys.argv.append('changed')\n")

    with patch.object(sys, 'argv', ['launcher']):
        run_file_cmd = RunFileCmd(CliOptions(file=[str(file), 'argument']))
        assert sys.argv == ['launcher']
        run_file_cmd.run()
        run_file_cmd.run()

    captured = capsys.readouterr()
    expected = f"[{str(file)!r}, 'argument']\n"
    assert captured.out == expected * 2


def test_interleaved_run_file_cmd_instances_keep_their_argv(capsys, tmp_path):
    first = tmp_path / "first.py"
    first.write_text("import sys\nprint(sys.argv)\n")
    second = tmp_path / "second.py"
    second.write_text("import sys\nprint(sys.argv)\n")

    first_cmd = RunFileCmd(CliOptions(file=[str(first), 'first']))
    RunFileCmd(CliOptions(file=[str(second), 'second']))
    first_cmd.run()

    captured = capsys.readouterr()
    assert captured.out == f"[{str(first)!r}, 'first']\n"


def test_run_file_cmd_dunder_name_is_main(capsys, tmp_path):
    """Inside the script, ``__name__`` must equal ``'__main__'`` (#42159)."""
    file = tmp_path / "test.sage"
    file.write_text(
        "import __main__\n"
        "import sys\n"
        "print(__name__)\n"
        "print(__file__)\n"
        "print(__package__)\n"
        "print(__spec__)\n"
        "print(__loader__)\n"
        "print(__cached__)\n"
        "print(__main__.__name__)\n"
        "print(__main__.__file__)\n"
        "print(sys.modules['__mp_main__'] is __main__)\n"
        "try:\n"
        "    from . import options\n"
        "except ImportError:\n"
        "    print('relative import failed')\n"
    )
    options = CliOptions(file=[str(file)])
    RunFileCmd(options).run()
    captured = capsys.readouterr()
    assert captured.out == (
        "__main__\n"
        f"{file}\n"
        "None\n"
        "None\n"
        "None\n"
        "None\n"
        "__main__\n"
        f"{file}\n"
        "True\n"
        "relative import failed\n"
    )


def test_run_file_cmd_reset_preserves_dunder_name(capsys, tmp_path):
    """``reset()`` must not clobber any of the script module metadata."""
    file = tmp_path / "test.sage"
    file.write_text(
        "'''script docstring'''\n"
        "import __main__\n"
        "_metadata = (__name__, __file__, __package__, __loader__, "
        "__spec__, __cached__, __doc__)\n"
        "reset()\n"
        "print(_metadata == (__name__, __file__, __package__, __loader__, "
        "__spec__, __cached__, __doc__))\n"
        "print(__main__.__dict__ is globals())\n"
    )
    options = CliOptions(file=[str(file)])
    RunFileCmd(options).run()
    captured = capsys.readouterr()
    assert captured.out == "True\nTrue\n"


def test_run_file_cmd_relinks_lazy_imports_after_restore(capsys, monkeypatch, tmp_path):
    """Copied lazy imports must always resolve in the script namespace."""
    import sage.all as sage_all
    from sage.misc.lazy_import import LazyImport

    name = 'RunFileLazyImportProbe'
    proxy = LazyImport(
        'sage.rings.integer',
        'Integer',
        as_name=name,
        namespace=sage_all.__dict__,
    )
    monkeypatch.setitem(sage_all.__dict__, name, proxy)

    file = tmp_path / "test.py"
    file.write_text(
        "from sage.misc.lazy_import import attributes as _lazy_attributes\n"
        "print(_lazy_attributes(RunFileLazyImportProbe)['_namespace'] is globals())\n"
        "_resolved = RunFileLazyImportProbe._get_object()\n"
        "print(globals()['RunFileLazyImportProbe'] is _resolved)\n"
        "reset()\n"
        "print(_lazy_attributes(RunFileLazyImportProbe)['_namespace'] is globals())\n"
        "_resolved = RunFileLazyImportProbe._get_object()\n"
        "print(globals()['RunFileLazyImportProbe'] is _resolved)\n"
        "restore('RunFileLazyImportProbe')\n"
        "print(_lazy_attributes(RunFileLazyImportProbe)['_namespace'] is globals())\n"
        "_resolved = RunFileLazyImportProbe._get_object()\n"
        "print(globals()['RunFileLazyImportProbe'] is _resolved)\n"
    )
    RunFileCmd(CliOptions(file=[str(file)])).run()

    captured = capsys.readouterr()
    assert captured.out == "True\n" * 6


def test_cython_compiles_from_the_original_working_directory(monkeypatch, tmp_path):
    """Relative Cython build inputs use the cwd recorded by the bootstrap."""
    source_cwd = tmp_path / 'source'
    source_cwd.mkdir()
    runtime_cwd = tmp_path / 'runtime'
    runtime_cwd.mkdir()
    file = source_cwd / 'test.spyx'
    file.write_text('')
    observed = []

    def fake_load_cython(filename):
        observed.append((filename, os.getcwd()))
        return 'pass'

    monkeypatch.setattr('sage.cli.run_file_cmd.load_cython', fake_load_cython)
    monkeypatch.chdir(runtime_cwd)
    _run_file(str(file), {'__name__': '__mp_main__'}, 'bootstrap.py', str(source_cwd))

    assert observed == [(str(file), str(source_cwd))]
    assert os.getcwd() == str(runtime_cwd)


@pytest.mark.parametrize(
    ('suffix', 'source', 'expected'),
    [
        (
            '.py',
            "from multiprocessing import get_context\n"
            "value = ZZ(5)\n"
            "class Result:\n"
            "    def __init__(self, value):\n"
            "        self.value = value\n"
            "    def __repr__(self):\n"
            "        return f'Result({self.value})'\n"
            "def grandchild(queue):\n"
            "    queue.put(Result(int(value + 2)))\n"
            "def worker(queue):\n"
            "    context = get_context('spawn')\n"
            "    child_queue = context.Queue()\n"
            "    child = context.Process(target=grandchild, args=(child_queue,))\n"
            "    child.start()\n"
            "    child.join(30)\n"
            "    queue.put((child.exitcode, child_queue.get(timeout=10)))\n",
            '(0, Result(7))',
        ),
        (
            '.sage',
            "from multiprocessing import get_context\n"
            "R.<x> = QQ[]\n"
            "def worker(queue):\n"
            "    queue.put(str((x + 1)^2))\n",
            'x^2 + 2*x + 1',
        ),
    ],
)
def test_run_file_cmd_supports_multiprocessing_spawn(
    suffix, source, expected, tmp_path
):
    """Spawn must reload both Python and preparsed Sage input correctly."""
    child_dir = tmp_path / 'child'
    child_dir.mkdir()
    file = tmp_path / f'test{suffix}'
    file.write_text(
        source
        + "if __name__ == '__main__':\n"
        + "    import os\n"
        + "    os.chdir('child')\n"
        + "    context = get_context('spawn')\n"
        + "    queue = context.Queue()\n"
        + "    process = context.Process(target=worker, args=(queue,))\n"
        + "    process.start()\n"
        + "    process.join(30)\n"
        + "    print(process.exitcode)\n"
        + "    print(queue.get(timeout=10))\n"
    )

    result = subprocess.run(
        [sys.executable, '-m', 'sage.cli', '-q', file.name],
        cwd=tmp_path,
        capture_output=True,
        check=False,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == f"0\n{expected}\n"


def test_run_file_cmd_keeps_main_for_excepthook_and_atexit(tmp_path):
    """The script main module must live for the whole CLI process."""
    file = tmp_path / 'test.py'
    file.write_text(
        "import atexit\n"
        "import pickle\n"
        "import sys\n"
        "class C:\n"
        "    pass\n"
        "def check(label):\n"
        "    import __main__\n"
        "    restored = pickle.loads(pickle.dumps(C()))\n"
        "    print(label, __main__.__file__ == __file__, "
        "type(restored) is C, flush=True)\n"
        "sys.excepthook = lambda *_: check('hook')\n"
        "atexit.register(check, 'atexit')\n"
        "raise RuntimeError('boom')\n"
    )

    result = subprocess.run(
        [sys.executable, '-m', 'sage.cli', '-q', str(file)],
        capture_output=True,
        check=False,
        text=True,
        timeout=60,
    )
    assert result.returncode != 0
    assert result.stdout == "hook True True\natexit True True\n"
