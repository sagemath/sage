import argparse
import os
import sys
import types

from sage.all import sage_globals
from sage.cli.options import CliOptions
from sage.misc.lazy_import import clean_namespace
from sage.misc.sage_ostools import restore_cwd
from sage.misc.temporary_file import tmp_filename
from sage.repl.load import load_cython
from sage.repl.preparse import preparse_file_named


def _spawn_bootstrap(source_file, source_cwd):
    r"""Create a bootstrap script for ``multiprocessing``.

    The ``spawn`` and ``forkserver`` start methods reload ``__main__`` with
    :mod:`runpy`.  A Sage source file cannot be reloaded directly: ``.sage``
    files first need preparsing, and Python files need the Sage global
    namespace.  The small bootstrap created here sends child processes through
    :func:`_run_file` instead.
    """
    bootstrap_file = tmp_filename(name='_sage_run_file_', ext='.py')
    with open(bootstrap_file, 'w', encoding='ascii') as f:
        f.write(
            'from sage.cli.run_file_cmd import _run_file as _sage_run_file\n'
            f'_sage_run_file({source_file!a}, globals(), '
            f'{bootstrap_file!a}, {source_cwd!a})\n'
        )
    return bootstrap_file


def _install_multiprocessing_spawn_support():
    """Teach ``multiprocessing`` how to replay the current Sage source."""
    import multiprocessing.spawn

    original = multiprocessing.spawn.get_preparation_data
    if getattr(original, '_sage_run_file', False):
        return

    def get_preparation_data(name):
        data = original(name)
        main_module = sys.modules.get('__main__')
        bootstrap_file = getattr(main_module, '__sage_spawn_bootstrap__', None)
        if bootstrap_file is None:
            return data

        # Keep normal script semantics (``__spec__ is None``) while replacing
        # only multiprocessing's private replay instruction.
        data.pop('init_main_from_name', None)
        data.pop('init_main_from_path', None)
        data['init_main_from_path'] = bootstrap_file
        return data

    get_preparation_data._sage_run_file = True
    get_preparation_data.__wrapped__ = original
    multiprocessing.spawn.get_preparation_data = get_preparation_data


def _run_file(source_file, namespace, bootstrap_file, source_cwd):
    r"""Execute ``source_file`` in a Sage-populated module namespace."""
    module_name = namespace.get('__name__', '__main__')
    namespace.update(sage_globals())
    clean_namespace(namespace)
    namespace.pop('_sage_run_file', None)
    namespace['__name__'] = module_name
    namespace['__file__'] = source_file
    namespace['__package__'] = None
    namespace['__loader__'] = None
    namespace['__spec__'] = None
    namespace['__cached__'] = None
    namespace['__doc__'] = None
    namespace['__sage_spawn_bootstrap__'] = bootstrap_file
    _install_multiprocessing_spawn_support()

    # ``runpy.run_path`` temporarily points argv[0] at the bootstrap.  User
    # code should continue to see its own filename.
    sys.argv[0] = source_file

    input_file = source_file
    if input_file.endswith('.sage'):
        input_file = str(preparse_file_named(input_file))

    if input_file.endswith('.pyx') or input_file.endswith('.spyx'):
        # Keep the historical Cython extension-module semantics.  Its
        # top-level code executes in the generated extension module before
        # public names are imported into this command-line main namespace.
        with restore_cwd(source_cwd):
            source = load_cython(input_file)
        exec(compile(source, tmp_filename(), 'exec'), namespace)
    else:
        with open(input_file, 'rb') as f:
            source = f.read()
        exec(compile(source, input_file, 'exec'), namespace)


class _RunFileArgs(argparse.Action):
    """Collect the script filename and all following arguments."""

    def __call__(self, parser, namespace, values, option_string=None):
        if values and values[0] == "--":
            values = values[1:]
        setattr(namespace, self.dest, values)


class RunFileCmd:
    @staticmethod
    def extend_parser(parser: argparse.ArgumentParser):
        r"""
        Extend the parser with the "run file" command.

        INPUT:

        - ``parsers`` -- the parsers to extend.

        OUTPUT:

        - the extended parser.
        """
        parser.add_argument(
            "file",
            nargs=argparse.REMAINDER,
            action=_RunFileArgs,
            help="execute the given file as sage code",
        )

    def __init__(self, options: CliOptions):
        r"""
        Initialize the command.
        """
        self.options = options
        # Rebuild ``sys.argv`` so that the executed script sees the same
        # arguments it would receive under plain Python:
        # ``sys.argv == [<script>, *script_args]``.  ``script_args`` is the
        # tail captured by the run-file parser after the script filename.
        script = options.file[0] if options.file else sys.argv[0]
        forwarded = list(options.file[1:] if options.file else [])
        self.argv = [script, *forwarded]

    def run(self) -> int:
        r"""
        Execute the given command.
        """
        # Rebuild this for every invocation: constructing another command, or
        # a previous script mutating ``sys.argv``, must not change our input.
        sys.argv = self.argv.copy()
        source_file = self.argv[0]
        absolute_source_file = os.path.abspath(source_file)
        source_cwd = os.getcwd()
        bootstrap_file = _spawn_bootstrap(absolute_source_file, source_cwd)

        main_module = types.ModuleType('__main__')
        sys.modules['__main__'] = main_module
        # ``multiprocessing`` aliases the process main module under this name
        # so objects created while a spawned child replays the script can be
        # unpickled by the parent process.
        sys.modules['__mp_main__'] = main_module
        _run_file(source_file, main_module.__dict__, bootstrap_file, source_cwd)
        return 0
