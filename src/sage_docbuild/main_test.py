"""Executable package tests for :mod:`sage_docbuild.__main__`.

Sage deliberately skips doctests in ``__main__.py`` files, so command-line
parsing behavior belongs in pytest tests rather than inert docstring examples.
"""

import errno
import hashlib
import json
import os
import stat
import sys
from importlib.util import find_spec

import pytest

if find_spec('sage_docbuild') is None or find_spec('sphinx') is None:
    pytest.skip(
        'sage_docbuild and Sphinx are optional build dependencies',
        allow_module_level=True,
    )

from sage_docbuild import __main__ as docbuild_main
from sage_docbuild import builders
from sage_docbuild.__main__ import (
    _options_generation_lock,
    command_line_args,
    format_columns,
    setup_parser,
    source_dir_for_help,
)
from sage_docbuild.build_options import BuildOptions
from sage_docbuild.builders import (
    SingleFileBuilder,
    _single_file_output_owned,
    get_builder,
)


def test_format_columns():
    assert format_columns(['tutorial', 'website'], width=30) == (
        '    tutorial   website    \n'
    )
    assert format_columns([]) == '\n'


@pytest.mark.parametrize(
    ('options', 'expected'),
    [
        ('  --no-plot  ', ['--no-plot', 'reference', 'html']),
        (
            "--source '/a path/with spaces'",
            ['--source', '/a path/with spaces', 'reference', 'html'],
        ),
        (
            r'-S=-Dfoo=\usepackage{microtype}',
            [r'-S=-Dfoo=\usepackage{microtype}', 'reference', 'html'],
        ),
    ],
)
def test_command_line_args(monkeypatch, options, expected):
    monkeypatch.setattr(
        sys, 'argv', ['sage-docbuild', 'reference', 'html']
    )
    monkeypatch.setenv('SAGE_DOCBUILD_OPTS', options)
    assert command_line_args() == expected


def test_command_line_args_rejects_malformed_environment(monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['sage-docbuild'])
    monkeypatch.setenv('SAGE_DOCBUILD_OPTS', "--source 'unterminated")
    with pytest.raises(
        SystemExit, match='error: could not parse SAGE_DOCBUILD_OPTS:'
    ):
        command_line_args()


def test_source_dir_for_help(monkeypatch, tmp_path):
    source = tmp_path / 'documentation source'
    source.mkdir()
    expected = source.absolute()

    # Without an explicit prog, ArgumentParser names itself after sys.argv[0]
    # and raises IndexError when that list is empty.  This helper is handed
    # its arguments directly, so it must not depend on the process-global one.
    monkeypatch.setattr(sys, 'argv', [])
    assert source_dir_for_help(['-D', '--source', str(source)]) == expected
    assert source_dir_for_help(['-D', f'--source={source}']) == expected
    assert source_dir_for_help(['-D', f'--sour={source}']) == expected


def test_source_dir_for_help_accepts_short_attached_value(monkeypatch, tmp_path):
    source = tmp_path / 'documentation source'
    source.mkdir()
    argument = f'-s={source}'

    # argparse accepts an equals sign between a short option and its attached
    # value.  The help actions exit while parsing, so their pre-parser must
    # recognize the same spelling itself.
    parsed = setup_parser().parse_args([argument, 'reference', 'html'])
    assert parsed.source_dir == source
    assert source_dir_for_help(['-D', argument]) == source.absolute()


@pytest.mark.parametrize(
    'arguments',
    [
        ['-D', '--source'],
        ['-D', '--source', '--help'],
        ['-D', '--source', '--not-an-option'],
    ],
)
def test_source_dir_for_help_requires_value(arguments):
    with pytest.raises(
        SystemExit,
        match=r'error: argument -s/--source: expected one argument',
    ):
        source_dir_for_help(arguments)


def test_source_dir_for_help_rejects_missing_directory(tmp_path):
    missing = tmp_path / 'missing'
    with pytest.raises(
        SystemExit,
        match=r'error: source directory .*missing does not exist',
    ):
        source_dir_for_help([f'--source={missing}', '-D'])


def test_source_dir_for_help_honors_end_of_options(monkeypatch, tmp_path):
    source = tmp_path / 'default'
    source.mkdir()
    monkeypatch.setenv('SAGE_DOC_SRC', str(source))
    assert source_dir_for_help(
        ['-D', '--', '--source', str(tmp_path / 'missing')]
    ) == source.absolute()


def test_source_dir_for_help_treats_negative_value_as_path(
    monkeypatch, tmp_path
):
    source = tmp_path / '-3'
    source.mkdir()
    monkeypatch.chdir(tmp_path)
    assert source_dir_for_help(['-D', '--source', '-3']) == source.absolute()


def test_single_file_builder_creates_default_output_root(monkeypatch, tmp_path):
    import sage.env

    source = tmp_path / 'issue_42481.py'
    source.write_text('"""A module documented on its own."""\n')
    dot_sage = tmp_path / 'new-sage-home'
    monkeypatch.setattr(sage.env, 'DOT_SAGE', str(dot_sage))
    monkeypatch.setattr(builders, 'DOT_SAGE', str(dot_sage), raising=False)
    from sage.env import DOT_SAGE as resolved
    assert resolved == str(dot_sage)
    options = BuildOptions(
        source_dir=tmp_path,
        output_dir=tmp_path / 'unused-documentation-output',
        output_dir_given=False,
    )

    previous_umask = None
    if os.name != 'nt':
        previous_umask = os.umask(0o022)
    try:
        builder = get_builder(f'file={source}', options)
    finally:
        if previous_umask is not None:
            os.umask(previous_umask)
    expected = dot_sage / 'docbuild' / source.stem

    assert isinstance(builder, SingleFileBuilder)
    assert builder._options is options
    assert builder._single_file_base_dir == expected
    if os.name != 'nt':
        assert stat.S_IMODE(dot_sage.stat().st_mode) == 0o700
    assert _single_file_output_owned(expected)
    assert (expected / 'source' / 'index.rst').is_file()


def test_ambiguous_source_abbreviation_matches_argparse(monkeypatch, tmp_path):
    source = tmp_path / 'default'
    source.mkdir()
    monkeypatch.setenv('SAGE_DOC_SRC', str(source))
    # ``--s`` is ambiguous with --sphinx-opts, so the help pre-parser leaves
    # it alone and argparse rejects it rather than selecting a source tree.
    assert source_dir_for_help(['-D', f'--s={tmp_path}']) == source.absolute()
    with pytest.raises(SystemExit):
        setup_parser().parse_args(['-D', f'--s={tmp_path}'])


def test_needs_options_page_reads_files_instead_of_consulting_umask(
    monkeypatch, tmp_path
):
    options = tmp_path / 'en' / 'reference' / 'repl' / 'options.txt'
    options.parent.mkdir(parents=True)
    content = 'usage: sage [options]\n'
    options.write_text(content, encoding='utf-8')

    input_digest = 'cli-input-digest'
    manifest = tmp_path / 'options-manifest.json'
    manifest.write_text(
        json.dumps(
            {
                'version': docbuild_main._OPTIONS_MANIFEST_VERSION,
                'inputs': input_digest,
                'output': hashlib.sha256(content.encode('utf-8')).hexdigest(),
            }
        ),
        encoding='utf-8',
    )

    # Both files are readable by this process.  Their mode need not match the
    # current umask: ACLs and the process identity, not creation defaults,
    # determine whether a read succeeds.
    options.chmod(0o600)
    manifest.chmod(0o600)
    monkeypatch.setattr(
        docbuild_main, '_options_input_digest', lambda root: input_digest
    )
    monkeypatch.setattr(
        docbuild_main,
        '_file_mode',
        lambda: pytest.fail('_needs_options_page consulted the current umask'),
    )

    assert not docbuild_main._needs_options_page(
        options, tmp_path / 'source root', manifest
    )


# The lock around the generated options page only keeps parallel builds from
# repeating the same atomic write, so no way of failing to hold it may fail a
# build.  Each test below breaks one step and asks that the body still run.

def _manifest(tmp_path):
    return tmp_path / 'sub' / 'sage-cli-options-manifest.json'


# Permission bits do not stop root, so the two tests that revoke them cannot
# say anything when the suite runs as root, as it does in some containers.
skip_as_root = pytest.mark.skipif(
    hasattr(os, 'geteuid') and os.geteuid() == 0,
    reason='root ignores the permissions these tests revoke',
)


@skip_as_root
def test_options_lock_survives_unopenable_lock_file(tmp_path):
    # A lock file left behind unreadable by an earlier build, in a directory
    # that is otherwise perfectly writable.
    manifest = _manifest(tmp_path)
    manifest.parent.mkdir(parents=True)
    lock_path = manifest.with_suffix('.lock')
    lock_path.touch()
    lock_path.chmod(0o000)
    try:
        with _options_generation_lock(manifest):
            ran = True
    finally:
        lock_path.chmod(0o600)
    assert ran


@skip_as_root
def test_options_lock_survives_unwritable_directory(tmp_path):
    manifest = _manifest(tmp_path)
    manifest.parent.mkdir(parents=True)
    manifest.parent.chmod(0o500)
    try:
        with _options_generation_lock(manifest):
            ran = True
    finally:
        manifest.parent.chmod(0o700)
    assert ran


@pytest.mark.parametrize('failing_operation', ['LOCK_EX', 'LOCK_UN'])
def test_options_lock_survives_refused_locking(
    monkeypatch, tmp_path, failing_operation
):
    # A filesystem that answers ENOLCK, as a network home directory may.
    fcntl = pytest.importorskip('fcntl')
    refused = getattr(fcntl, failing_operation)
    real_flock = fcntl.flock

    def flock(file, operation):
        if operation == refused:
            raise OSError(errno.ENOLCK, 'No locks available')
        return real_flock(file, operation)

    monkeypatch.setattr(fcntl, 'flock', flock)
    manifest = _manifest(tmp_path)
    with _options_generation_lock(manifest):
        ran = True
    assert ran
