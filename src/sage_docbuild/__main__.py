# sage.doctest: needs sphinx
r"""
Sage docbuild main

This module defines the Sage documentation build command::

    sage --docbuild [OPTIONS] DOCUMENT (FORMAT | COMMAND)

If ``FORMAT`` is given, it builds ``DOCUMENT`` in ``FORMAT``. If ``COMMAND`` is
given, it returns information about ``DOCUMENT``.

Run ``sage --docbuild`` to get detailed explanations about
arguments and options.

Positional arguments::

  DOCUMENT              name of the document to build. It can be either one of
                        the documents listed by -D or 'file=/path/to/FILE' to
                        build documentation for this specific file.
  FORMAT or COMMAND     document output format (or command)

Standard options::

  -h, --help            show a help message and exit
  -H, --help-all        show an extended help message and exit
  -D, --documents       list all available DOCUMENTs
  -F, --formats         list all output FORMATs
  -C DOC, --commands DOC
                        list all COMMANDs for DOCUMENT DOC; use 'all' to list all
  -i, --inherited       include inherited members in reference manual; may be
                        slow, may fail for PDF output
  -u, --underscore      include variables prefixed with '_' in reference
                        manual; may be slow, may fail for PDF output
  -j, --mathjax, --jsmath
                        ignored for backwards compatibility
  --no-plot             do not include graphics auto-generated using the '.. plot' markup
  --include-tests-blocks
                        include TESTS blocks in the reference manual
  --no-pdf-links        do not include PDF links in DOCUMENT 'website';
                        FORMATs: html, json, pickle, web
  --live-doc            make Sage code blocks live for html FORMAT
  --warn-links          issue a warning whenever a link is not properly
                        resolved; passes '-n' (nitpicky) to sphinx, except
                        for first-pass (inventory) builds
  --check-nested        check picklability of nested classes in DOCUMENT 'reference'
  --no-prune-empty-dirs
                        do not prune empty directories in the documentation source
  --use-cdns            assume internet connection and use CDNs; in particular,
                        use MathJax CDN
  -N, --no-colors       do not color output; does not affect children
  -q, --quiet           work quietly; same as --verbose=0
  -v LEVEL, --verbose LEVEL
                        report progress at LEVEL=0 (quiet), 1 (normal), 2
                        (info), or 3 (debug); does not affect children
  -o DIR, --output DIR  if DOCUMENT is a single file ('file=...'), write output
                        to this directory

Advanced options::

  Use these options with care.

  -S OPTS, --sphinx-opts OPTS
                        pass comma-separated OPTS to sphinx-build; must precede
                        OPTS with '=', as in '-S=-q,-aE' or '-S="-q,-aE"'
  -U, --update-mtimes   before building reference manual, update modification
                        times for auto-generated reST files
  -k, --keep-going      Do not abort on errors but continue as much as possible
                        after an error
  --all-documents ARG   if ARG is 'reference', list all subdocuments of
                        en/reference. If ARG is 'all', list all main documents
"""

import argparse
import functools
import hashlib
import json
import logging
import os
import shlex
import stat
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING

try:
    import fcntl
except ModuleNotFoundError:
    # Windows has no fcntl, and Sage does not build there.  Listing the
    # documents of a source tree is part of configuring a build, though, and
    # meson runs it wherever it configures: this module has to import even
    # where the locks below cannot be taken.  See _lock_exclusive().
    fcntl = None

import sphinx.ext.intersphinx

from . import build_options
from .builders import (
    _output_formats,
    get_all_documents,
    get_all_reference_documents,
    get_builder,
    get_documents,
)

if TYPE_CHECKING:
    from .build_options import BuildOptions

logger = logging.getLogger(__name__)

_OPTIONS_MANIFEST_VERSION = 1
_OPTIONS_MANIFEST_NAME = 'sage-cli-options-manifest.json'
_CLI_HELP_COLUMNS = '80'


def _stat_identity(info) -> tuple[int, ...]:
    """Return metadata used to detect a file changing and changing back."""

    return (info.st_dev, info.st_ino, info.st_mode, info.st_nlink,
            info.st_uid, info.st_gid, info.st_size, info.st_mtime_ns,
            info.st_ctime_ns)


def _path_status(path: Path):
    """Return a followed status, suppressing absence but no other error."""

    try:
        path.lstat()
    except FileNotFoundError:
        return None
    try:
        return path.stat()
    except FileNotFoundError as error:
        raise OSError(f'broken or moving path: {path}') from error


def _path_is_directory(path: Path) -> bool:
    """Test for a directory without hiding permission or I/O errors."""

    info = _path_status(path)
    return info is not None and stat.S_ISDIR(info.st_mode)


def _path_is_regular_file(path: Path) -> bool:
    """Test for a regular file without hiding permission or I/O errors."""

    info = _path_status(path)
    return info is not None and stat.S_ISREG(info.st_mode)


@functools.cache
def _file_mode() -> int:
    """Return the ordinary creation mode selected by the process umask.

    Probing costs a temporary directory and half a dozen system calls, and
    every generated file asks.  Nothing that runs while the documentation
    sources are written changes the umask, so the first answer stands for
    the whole run.
    """

    directory = Path(tempfile.mkdtemp(prefix='sage-docbuild-mode-'))
    descriptor = None
    probe = directory / 'probe'
    try:
        descriptor = os.open(probe, os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                             0o666)
        return stat.S_IMODE(os.fstat(descriptor).st_mode)
    finally:
        if descriptor is not None:
            os.close(descriptor)
        probe.unlink(missing_ok=True)
        directory.rmdir()


def _file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as file:
        for block in iter(lambda: file.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _input_file_digest(path: Path) -> tuple[str, tuple[tuple[int, ...], ...]]:
    """Hash a regular input and return stable path/target identities."""

    path_before = path.lstat()
    target_before = path.stat()
    if not stat.S_ISREG(target_before.st_mode):
        raise OSError(f'command-line source is not a regular file: {path}')
    digest = hashlib.sha256()
    with path.open('rb') as file:
        opened_before = os.fstat(file.fileno())
        if _stat_identity(opened_before) != _stat_identity(target_before):
            raise OSError(f'command-line source changed before reading: {path}')
        for block in iter(lambda: file.read(1024 * 1024), b''):
            digest.update(block)
        opened_after = os.fstat(file.fileno())
    path_after = path.lstat()
    target_after = path.stat()
    if (_stat_identity(path_after) != _stat_identity(path_before)
            or _stat_identity(target_after) != _stat_identity(target_before)
            or _stat_identity(opened_after) != _stat_identity(target_before)):
        raise OSError(f'command-line source changed while reading: {path}')
    identities = (_stat_identity(path_before), _stat_identity(target_before))
    return digest.hexdigest(), identities


def _walk_cli_tree(root: Path):
    """Yield the CLI tree without suppressing traversal errors."""

    for path in sorted(root.iterdir()):
        info = path.lstat()
        if stat.S_ISDIR(info.st_mode) and path.name == '__pycache__':
            continue
        yield path
        if stat.S_ISDIR(info.st_mode):
            yield from _walk_cli_tree(path)


@contextmanager
def _lock_exclusive(file):
    """Hold an exclusive lock on ``file`` where the platform provides one.

    The lock is an optimization, not a correctness requirement: it keeps
    parallel builds from each running the help subprocess and rewriting the
    same page.  Racing writers would still publish identical content, since
    every one of them writes atomically and the inputs are checked before and
    after.  So a platform without :mod:`fcntl`, and a filesystem that refuses
    to lock - a network home directory answering ``ENOLCK``, say - run without
    the lock rather than failing a build that has nothing else wrong with it.
    """
    if fcntl is None:
        yield
        return
    try:
        fcntl.flock(file, fcntl.LOCK_EX)
    except OSError as error:
        logger.debug('Not locking %s: %s', file.name, error)
        yield
        return
    try:
        yield
    finally:
        # Closing the handle releases the lock in any case, so a refused
        # release must not escape a body that has already run to completion.
        try:
            fcntl.flock(file, fcntl.LOCK_UN)
        except OSError as error:
            logger.debug('Not unlocking %s: %s', file.name, error)


@contextmanager
def _options_generation_lock(manifest: Path):
    """Serialize writers of ``options.txt`` and its manifest.

    Opening the lock is as optional as taking it, and for the same reason:
    the manifest and the page it describes are published atomically, so a
    build that cannot open the lock still writes correctly.  The lock file
    persists between builds, and one left behind unreadable by a partial or a
    differently permissioned build must not stop an otherwise writable tree.
    """

    lock_path = manifest.with_suffix('.lock')
    try:
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        lock = lock_path.open('a+b')
    except OSError as error:
        logger.debug('Not opening %s: %s', lock_path, error)
        yield
        return
    with lock, _lock_exclusive(lock):
        yield


def _write_text_atomic(path: Path, content: str) -> None:
    """Publish ``content`` atomically, preserving an unchanged file's mtime."""

    encoded = content.encode('utf-8')
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = _file_mode()
    try:
        info = path.lstat()
    except FileNotFoundError:
        info = None
    if (info is not None
            and not (stat.S_ISREG(info.st_mode)
                     or stat.S_ISLNK(info.st_mode))):
        raise OSError(
            f'refusing to replace non-file generated path: {path}')
    if info is not None and stat.S_ISREG(info.st_mode):
        unchanged = path.read_bytes() == encoded
        if unchanged:
            if stat.S_IMODE(info.st_mode) != mode:
                path.chmod(mode)
            return
    handle, temporary = tempfile.mkstemp(dir=path.parent, prefix=path.name,
                                         suffix='.tmp')
    try:
        with os.fdopen(handle, 'w', encoding='utf-8') as file:
            file.write(content)
            os.fchmod(file.fileno(), mode)
        # os.replace() atomically replaces a regular file or a symlink.  A
        # directory or special file cannot be a generated artifact and is
        # preserved rather than recursively removed.
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def format_columns(lst, align='<', cols=None, indent=4, pad=3, width=80):
    """
    Utility function that formats a list as a simple table and returns
    a Unicode string representation.

    The number of columns is
    computed from the other options, unless it's passed as a keyword
    argument.  For help on Python's string formatter, see

    https://docs.python.org/library/string.html#format-string-syntax

    A source tree holding no document at all is listed as the empty table it
    is, rather than making the help fail.
    """
    # Can we generalize this (efficiently) to other / multiple inputs
    # and generators?
    if not lst:
        return "\n"
    size = max(map(len, lst)) + pad
    if cols is None:
        import math
        cols = math.trunc((width - indent) / size)
    s = " " * indent
    for i in range(len(lst)):
        if i != 0 and i % cols == 0:
            s += "\n" + " " * indent
        s += "{0:{1}{2}}".format(lst[i], align, size)
    s += "\n"
    return s


def help_usage(s="", compact=False):
    """
    Append and return a brief usage message for the Sage documentation builder.

    If 'compact' is False, the function adds a final newline character.
    """
    s += "sage --docbuild [OPTIONS] DOCUMENT (FORMAT | COMMAND)"
    if not compact:
        s += "\n"
    return s


def help_description(s="", compact=False):
    """
    Append and return a brief description of the Sage documentation builder.

    If 'compact' is ``False``, the function adds a final newline character.
    """
    s += "Build or return information about Sage documentation. "
    s += "A DOCUMENT and either a FORMAT or a COMMAND are required."
    if not compact:
        s += "\n"
    return s


def help_examples(s=""):
    """
    Append and return some usage examples for the Sage documentation builder.
    """
    s += "Examples:\n"
    s += "    sage --docbuild -C all\n"
    s += "    sage --docbuild constructions pdf\n"
    s += "    sage --docbuild reference html -jv3\n"
    s += "    sage --docbuild reference print_unincluded_modules\n"
    s += "    sage --docbuild developer html --sphinx-opts='-q,-aE' --verbose 2"
    return s


def command_line_args() -> list[str]:
    """
    Return the arguments that the builder is to parse.

    Those of :envvar:`SAGE_DOCBUILD_OPTS` come first, so that the ones written
    on the command line override them.

    """
    docbuild_opts = os.getenv("SAGE_DOCBUILD_OPTS", "")

    # This variable is command-line text, not merely a whitespace-separated
    # list: quoting is needed in particular for a source path containing a
    # space.
    try:
        lexer = shlex.shlex(docbuild_opts, posix=True)
        lexer.whitespace_split = True
        lexer.commenters = ''
        # Sphinx option values commonly contain LaTeX backslashes.  They have
        # already survived the user's shell before reaching the environment,
        # so treating a backslash as another layer of escaping would eat it.
        lexer.escape = ''
        options = list(lexer)
    except ValueError as error:
        raise SystemExit(f"error: could not parse SAGE_DOCBUILD_OPTS: {error}")
    return options + sys.argv[1:]


def source_dir_for_help(argv=None) -> Path:
    """
    Return the documentation sources that the lists printed by the help are
    read from.

    The help is printed while the command line is still being parsed, so a
    ``--source`` that follows the option asking for the list has not reached
    the parser yet: it is read from ``argv`` here, whatever its place, and the
    default of :func:`main` is used when there is none.  ``argv`` defaults to
    what :func:`main` parses, :envvar:`SAGE_DOCBUILD_OPTS` included.

    A directory named explicitly is not second-guessed: the lists have to
    describe the tree that was asked for, or none at all.

    An option left without its value is reported rather than passed over, as
    the parser proper would report it.

    Unique long-option abbreviations accepted by :mod:`argparse` designate the
    same source here.

    The end-of-options marker leaves source-looking positional arguments to
    :mod:`argparse` instead of treating them as options.
    """
    if argv is None:
        argv = command_line_args()

    parser = setup_parser()
    option_strings = {
        option
        for action in parser._actions
        for option in action.option_strings
    }
    long_options = {option for option in option_strings
                    if option.startswith('--')}

    def is_source_option(option):
        if option == '--source':
            return True
        matches = {candidate for candidate in long_options
                   if candidate.startswith(option)}
        return matches == {'--source'}

    def is_recognized_option(argument):
        # This also recognizes unknown option-looking strings, as argparse
        # does when deciding that ``--source --unknown`` has no value, while
        # leaving negative numbers available as path values.
        return parser._parse_optional(argument) is not None

    source = None
    for i, arg in enumerate(argv):
        if arg == '--':
            break
        long_option, separator, attached = arg.partition('=')
        if arg == '-s' or (not separator and is_source_option(long_option)):
            if i + 1 == len(argv) or is_recognized_option(argv[i + 1]):
                raise SystemExit(
                    'error: argument -s/--source: expected one argument')
            source = Path(argv[i + 1])
        elif separator and (long_option == '-s'
                            or is_source_option(long_option)):
            source = Path(attached)
        elif arg.startswith('-s') and len(arg) > 2:
            source = Path(arg[2:])
    if source is not None:
        if not _path_is_directory(source):
            raise SystemExit(f"error: source directory {source} does not exist")
        return source.absolute()
    source = Path(os.environ.get('SAGE_DOC_SRC', 'src/doc'))
    if not _path_is_directory(source):
        # The default is relative to the root of a checkout, and the builder
        # runs from anywhere: fall back on the tree it was installed from.
        from sage.env import SAGE_DOC_SRC
        source = Path(SAGE_DOC_SRC)
    return source.absolute()


def sage_root_for(source_dir: Path):
    """
    Return a Sage source tree able to generate documentation sources, or
    ``None``.

    A distribution ships the generated sources and has no such tree.

    The tree holding ``source_dir`` comes first: the sources of a second
    checkout are generated from the package metadata of *that* checkout, not
    from the one that happens to be running.

    Every directory on the way up is tried, first of the path with its symbolic
    links resolved and then of the path as it was given.  The resolved tree has
    to win when a link kept below one checkout points into another one.
    ``<root>/src/doc`` is the usual place, but meson writes into
    ``<root>/build/<tag>/src/doc``.
    """
    from sage.env import SAGE_ROOT

    def ancestors(directory: Path):
        return [directory, *directory.parents]

    resolved = source_dir.resolve()
    candidates = ancestors(resolved)
    if resolved != source_dir:
        candidates += ancestors(source_dir)
    if SAGE_ROOT:
        candidates.append(Path(SAGE_ROOT))
    # This file sits in <root>/src/sage_docbuild.
    candidates.append(Path(__file__).parent.parent.parent)
    for candidate in candidates:
        if _path_is_regular_file(
                candidate / 'tools' / 'bootstrap-docs.py'):
            return candidate
    return None


def generate_doc_sources(source_dir: Path) -> None:
    """
    Write the documentation sources that the repository does not keep.

    The installation instructions, the page of every package and the options of
    the command line are generated from the package metadata and from the
    parser itself. A distribution ships them and meson writes them into its
    build directory, but a checkout has neither until something makes them,
    and Sphinx would only report an include that it cannot read.

    Which files the generator writes depends on the package metadata, so it is
    asked itself whether the tree is complete rather than probed for a file or
    two: a tree missing, or holding a truncated or an outdated copy of, any one
    of them builds into an error.

    Failing to write them is an error too. The build that follows would report
    an include that it cannot read, or document the packages of a Sage other
    than the one asked for, and neither says what went wrong.
    """
    import subprocess

    if not _path_is_directory(source_dir / 'en'):
        # Not a tree of documentation sources: writing hundreds of generated
        # files into it would be of no use to anyone.
        return
    root = sage_root_for(source_dir)
    if root is None:
        # A distribution ships the sources and carries no generator.
        return
    bootstrap = root / 'tools' / 'bootstrap-docs.py'
    options = source_dir / 'en' / 'reference' / 'repl' / 'options.txt'
    options_manifest = (source_dir / 'en' / 'installation' / '__pycache__'
                        / _OPTIONS_MANIFEST_NAME)

    # The generator reads the package metadata of the tree it writes into, and
    # finds that tree through the environment, which names the running one.
    env = {**os.environ, 'SAGE_ROOT': str(root), 'SAGE_SRC': str(root / 'src'),
           # argparse otherwise wraps its captured help according to whichever
           # terminal width happened to leak into this noninteractive build.
           'COLUMNS': _CLI_HELP_COLUMNS,
           # Bytecode caches are not CLI sources and must not change source
           # directory identities while the help subprocess is inspected.
           'PYTHONDONTWRITEBYTECODE': '1'}

    def run(command, capture=False):
        logger.warning('Generating documentation sources: %s',
                       ' '.join(str(part) for part in command))
        return subprocess.run(command, check=True, text=True, env=env,
                              capture_output=capture)

    try:
        # The directory comes first, so that a checkout whose generator is
        # older than --check writes the sources instead of taking the option
        # for the directory to write them into.
        check = subprocess.run([sys.executable, str(bootstrap),
                                str(source_dir), '--check'],
                               check=False, text=True, env=env,
                               capture_output=True)
        if check.stderr:
            logger.warning('%s', check.stderr.rstrip())
        if check.returncode:
            logger.debug('%s', check.stdout)
            run([sys.executable, str(bootstrap), str(source_dir)])
        # A complete installed source tree can be read-only.  Do not create or
        # open a writer lock unless the page actually needs to change; once a
        # lock is needed, check again in case another process won the race.
        if _needs_options_page(options, root, options_manifest):
            with _options_generation_lock(options_manifest):
                if not _needs_options_page(options, root, options_manifest):
                    return
                inputs, race = _options_input_state(root)
                usage = run([sys.executable, '-c', _CLI_HELP_PROGRAM,
                             str(root / 'src' / 'sage' / 'cli'), '--help'],
                            capture=True).stdout
                if not usage.strip():
                    raise OSError(
                        f'{root}: the command line described no options')
                after, after_race = _options_input_state(root)
                if (after, after_race) != (inputs, race):
                    raise OSError(
                        f'{root}: sage.cli changed while its help was read')
                _write_text_atomic(options, usage)
                manifest = {
                    'version': _OPTIONS_MANIFEST_VERSION,
                    'inputs': after,
                    'output': hashlib.sha256(
                        usage.encode('utf-8')).hexdigest(),
                }
                _write_text_atomic(
                    options_manifest,
                    json.dumps(manifest, indent=2, sort_keys=True) + '\n')
                final, final_race = _options_input_state(root)
                if (final, final_race) != (inputs, race):
                    raise OSError(
                        f'{root}: sage.cli changed while its help was written')
    except (OSError, subprocess.SubprocessError) as error:
        raise SystemExit('error: could not generate the documentation sources '
                         f'of {source_dir}: {error}')


# Print the ``--help`` of the command line whose package sits in the directory
# named by the first argument.  Putting the package on ``sys.path`` would not
# be enough: an installed Sage answers for ``sage.cli`` through a finder of its
# own, which comes before anything the path leads to, so this hands the import
# system the directory that was asked for and lets the rest of Sage be the one
# that is installed.
_CLI_HELP_PROGRAM = r'''
import importlib.util, os, sys

PACKAGE, DIRECTORY = 'sage.cli', sys.argv[1]
del sys.argv[1]


class Finder:
    def find_spec(self, fullname, path=None, target=None):
        if fullname != PACKAGE and not fullname.startswith(PACKAGE + '.'):
            return None
        rest = fullname[len(PACKAGE):].lstrip('.')
        base = os.path.join(DIRECTORY, *rest.split('.')) if rest else DIRECTORY
        initializer = os.path.join(base, '__init__.py')
        if os.path.isfile(initializer):
            return importlib.util.spec_from_file_location(
                fullname, initializer, submodule_search_locations=[base])
        if os.path.isfile(base + '.py'):
            return importlib.util.spec_from_file_location(fullname, base + '.py')
        # Returning None would let the editable-install finder answer from the
        # checkout that runs docbuild, silently mixing two sage.cli packages.
        raise ModuleNotFoundError(
            f'{fullname} is absent from the requested package at {DIRECTORY}')


for name in list(sys.modules):
    if name == PACKAGE or name.startswith(PACKAGE + '.'):
        del sys.modules[name]
sage_package = sys.modules.get('sage')
if sage_package is not None:
    vars(sage_package).pop('cli', None)
sys.meta_path.insert(0, Finder())
from sage.cli import main
sys.exit(main())
'''


def _options_input_state(root: Path) -> tuple[str, str]:
    """Return stable content and transient race digests for CLI help inputs."""

    cli = root / 'src' / 'sage' / 'cli'
    try:
        cli_info = cli.stat()
        initializer_info = (cli / '__init__.py').stat()
    except FileNotFoundError as error:
        raise FileNotFoundError(f'sage.cli does not exist in {root}') from error
    if (not stat.S_ISDIR(cli_info.st_mode)
            or not stat.S_ISREG(initializer_info.st_mode)):
        raise FileNotFoundError(f'sage.cli does not exist in {root}')
    content = hashlib.sha256(_CLI_HELP_PROGRAM.encode('utf-8'))
    race = hashlib.sha256()

    def add_file(label: str, path: Path) -> None:
        file_digest, identities = _input_file_digest(path)
        content.update(label.encode('utf-8', errors='surrogateescape'))
        content.update(b'\0')
        content.update(file_digest.encode('ascii'))
        content.update(b'\0')
        race.update(label.encode('utf-8', errors='surrogateescape'))
        race.update(b'\0')
        race.update(repr(identities).encode('ascii'))
        race.update(b'\0')

    def add_directory(label: str, path: Path) -> None:
        info = path.lstat()
        if not stat.S_ISDIR(info.st_mode):
            raise OSError(f'command-line source is not a directory: {path}')
        content.update(label.encode('utf-8', errors='surrogateescape'))
        content.update(b'\0directory\0')
        race.update(label.encode('utf-8', errors='surrogateescape'))
        race.update(b'\0')
        race.update(repr(_stat_identity(info)).encode('ascii'))
        race.update(b'\0')

    add_file('sage_docbuild/__main__.py', Path(__file__))
    add_directory('sage/cli', cli)
    content.update(_CLI_HELP_COLUMNS.encode('ascii'))
    content.update(b'\0')
    version = (f'{sys.implementation.name}:{sys.version_info.major}.'
               f'{sys.version_info.minor}')
    content.update(version.encode('ascii'))
    content.update(b'\0')
    for source in _walk_cli_tree(cli):
        label = f'sage/cli/{source.relative_to(cli).as_posix()}'
        info = source.lstat()
        if stat.S_ISDIR(info.st_mode):
            add_directory(label, source)
        elif source.suffix == '.py':
            add_file(label, source)
    return content.hexdigest(), race.hexdigest()


def _options_input_digest(root: Path) -> str:
    """Return the persistent content digest for CLI help inputs."""

    return _options_input_state(root)[0]


def _needs_options_page(options: Path, root: Path,
                        manifest_path: Path | None = None) -> bool:
    """
    Return whether the page listing the options of the command line has to be
    written again.

    It is the ``--help`` of the command line itself, so it goes out of date
    with the module implementing that command line, and an interrupted run
    leaves it empty.
    """
    if not _path_is_directory(options.parent):
        return False  # the document that includes it is not in this tree
    if manifest_path is None:
        # ``options`` normally is <source>/en/reference/repl/options.txt.
        try:
            source_dir = options.parents[3]
        except IndexError:
            return True
        manifest_path = (source_dir / 'en' / 'installation' / '__pycache__'
                         / _OPTIONS_MANIFEST_NAME)
    try:
        written = options.lstat()
        recorded = manifest_path.lstat()
    except OSError:
        return True
    if (not stat.S_ISREG(written.st_mode)
            or not stat.S_ISREG(recorded.st_mode)
            or not written.st_size):
        return True
    try:
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        inputs = _options_input_digest(root)
        return not (
            isinstance(manifest, dict)
            and manifest.get('version') == _OPTIONS_MANIFEST_VERSION
            and manifest.get('inputs') == inputs
            and manifest.get('output') == _file_digest(options)
        )
    except (OSError, ValueError, TypeError):
        return True


def help_documents():
    """
    Append and return a tabular list of documents available to the Sage
    documentation builder.
    """
    docs = get_documents(source_dir_for_help())
    s = "DOCUMENTs:\n"
    s += format_columns(docs)
    s += "\n"
    if 'reference' in docs:
        s += "Other valid document names take the form 'reference/DIR', where\n"
        s += "DIR is a subdirectory of src/doc/en/reference/.\n"
        s += "This builds just the specified part of the reference manual.\n"
    s += "DOCUMENT may also have the form 'file=/path/to/FILE', which builds\n"
    s += "the documentation for the specified file.\n"
    return s


def get_formats():
    """
    Return a list of output formats the Sage documentation builder
    will accept on the command-line.

    The formats are those of any builder, so this needs no document: building
    one here to ask it would need the options that are still being parsed.
    """
    formats = sorted(_output_formats() - {'html', 'pdf'})
    # The formats a reader asks for come first.
    return ['html', 'pdf'] + formats


def help_formats():
    """
    Append and return a tabular list of output formats available to
    the Sage documentation builder.
    """
    return "FORMATs:\n" + format_columns(get_formats())


def help_commands(name='all'):
    """
    Append and return a tabular list of commands, if any, the Sage
    documentation builder can run on the indicated document.  The
    default is to list all commands for all documents.
    """
    # To do: Generate the lists dynamically, using class attributes,
    # as with the Builders above.
    s = ""
    command_dict = {'reference': [
        'print_included_modules', 'print_modified_modules        (*)',
        'print_unincluded_modules', 'print_new_and_updated_modules (*)']}
    for doc in command_dict:
        if name == 'all' or doc == name:
            s += "COMMANDs for the DOCUMENT '" + doc + "':\n"
            s += format_columns(command_dict[doc])
            s += "(*) Since the last build.\n"
    return s


class help_message_long(argparse.Action):
    """
    Print an extended help message for the Sage documentation builder
    and exits.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        help_funcs = [help_usage, help_description, help_documents,
                      help_formats, help_commands]
        for f in help_funcs:
            print(f())
        parser.print_help()
        print(help_examples())
        sys.exit(0)


class help_message_short(argparse.Action):
    """
    Print a help message for the Sage documentation builder.

    The message includes command-line usage and a list of options.
    The message is printed only on the first call.  If error is True
    during this call, the message is printed only if the user hasn't
    requested a list (e.g., documents, formats, commands).
    """
    def __call__(self, parser, namespace, values, option_string=None):
        if not hasattr(namespace, 'printed_help'):
            parser.print_help()
            setattr(namespace, 'printed_help', 1)
        sys.exit(0)


class help_wrapper(argparse.Action):
    """
    A helper wrapper for command-line options to the Sage
    documentation builder that print lists, such as document names,
    formats, and document-specific commands.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        if option_string in ['-D', '--documents']:
            print(help_documents(), end="")
        if option_string in ['-F', '--formats']:
            print(help_formats(), end="")
        if self.dest == 'commands':
            print(help_commands(values), end="")
        setattr(namespace, 'printed_list', 1)
        sys.exit(0)


def setup_parser():
    """
    Set up and return a command-line ArgumentParser instance for the
    Sage documentation builder.
    """
    # Documentation: https://docs.python.org/library/argparse.html
    parser = argparse.ArgumentParser(prog='sage --docbuild',
                                     usage=help_usage(compact=True),
                                     description=help_description(compact=True),
                                     add_help=False)
    # Standard options. Note: We use explicit option.dest names
    # to avoid ambiguity.
    standard = parser.add_argument_group("Standard")
    standard.add_argument("-h", "--help", nargs=0, action=help_message_short,
                          help="show a help message and exit")
    standard.add_argument("-H", "--help-all", nargs=0, action=help_message_long,
                          help="show an extended help message and exit")
    standard.add_argument("-D", "--documents", nargs=0, action=help_wrapper,
                          help="list all available DOCUMENTs")
    standard.add_argument("-F", "--formats", nargs=0, action=help_wrapper,
                          help="list all output FORMATs")
    standard.add_argument("-C", "--commands", dest="commands",
                          type=str, metavar="DOC", action=help_wrapper,
                          help="list all COMMANDs for DOCUMENT DOC; use 'all' to list all")
    standard.add_argument("-i", "--inherited", dest="inherited",
                          action="store_true",
                          help="include inherited members in reference manual; may be slow, may fail for PDF output")
    standard.add_argument("-u", "--underscore", dest="underscore",
                          action="store_true",
                          help="include variables prefixed with '_' in reference manual; may be slow, may fail for PDF output")
    standard.add_argument("-j", "--mathjax", "--jsmath", dest="mathjax",
                          action="store_true",
                          help="ignored for backwards compatibility")
    standard.add_argument("--no-plot", dest="no_plot",
                          action="store_true",
                          help="do not include graphics auto-generated using the '.. plot' markup")
    standard.add_argument("--include-tests-blocks", dest="skip_tests", default=True,
                          action="store_false",
                          help="include TESTS blocks in the reference manual")
    standard.add_argument("--no-pdf-links", dest="no_pdf_links",
                          action="store_true",
                          help="do not include PDF links in DOCUMENT 'website'; FORMATs: html, json, pickle, web")
    standard.add_argument("--live-doc", dest="live_doc",
                          action="store_true",
                          help="make Sage code blocks live for html FORMAT")
    standard.add_argument("--warn-links", dest="warn_links",
                          action="store_true",
                          help="issue a warning whenever a link is not properly resolved; passes '-n' (nitpicky) to sphinx, except for first-pass (inventory) builds")
    standard.add_argument("--check-nested", dest="check_nested",
                          action="store_true",
                          help="check picklability of nested classes in DOCUMENT 'reference'")
    standard.add_argument("--no-prune-empty-dirs", dest="no_prune_empty_dirs",
                          action="store_true",
                          help="do not prune empty directories in the documentation source")
    standard.add_argument("--use-cdns", dest="use_cdns", default=False,
                          action="store_true",
                          help="assume internet connection and use CDNs; in particular, use MathJax CDN")
    standard.add_argument("-N", "--no-colors", dest="color",
                          action="store_false",
                          help="do not color output; does not affect children")
    standard.add_argument("-q", "--quiet", dest="verbose",
                          action="store_const", const=0,
                          help="work quietly; same as --verbose=0")
    standard.add_argument("-v", "--verbose", dest="verbose",
                          type=int, default=1, metavar="LEVEL",
                          action="store",
                          help="report progress at LEVEL=0 (quiet), 1 (normal), 2 (info), or 3 (debug); does not affect children")
    standard.add_argument("-s", "--source", dest="source_dir", type=Path,
                          default=None, metavar="DIR", action="store",
                          help="directory containing the documentation source files")
    standard.add_argument("-o", "--output", dest="output_dir", default=None,
                            type=Path,
                          metavar="DIR", action="store",
                          help="if DOCUMENT is a single file ('file=...'), write output to this directory")

    # Advanced options.
    advanced = parser.add_argument_group("Advanced",
                                         "Use these options with care.")
    advanced.add_argument("-S", "--sphinx-opts", dest="sphinx_opts",
                          type=str, metavar="OPTS",
                          action="store",
                          help="pass comma-separated OPTS to sphinx-build; must precede OPTS with '=', as in '-S=-q,-aE' or '-S=\"-q,-aE\"'")
    advanced.add_argument("-U", "--update-mtimes", dest="update_mtimes",
                          action="store_true",
                          help="before building reference manual, update modification times for auto-generated reST files")
    advanced.add_argument("-k", "--keep-going", dest="keep_going",
                          action="store_true",
                          help="Do not abort on errors but continue as much as possible after an error")
    advanced.add_argument("--all-documents", dest="all_documents",
                          type=str, metavar="ARG",
                          choices=['all', 'reference'],
                          help="if ARG is 'reference', list all subdocuments"
                          " of en/reference. If ARG is 'all', list all main"
                          " documents")
    parser.add_argument("document", nargs='?', type=str, metavar="DOCUMENT",
                        help="name of the document to build. It can be either one of the documents listed by -D or 'file=/path/to/FILE' to build documentation for this specific file.")
    parser.add_argument("format", nargs='?', type=str,
                        metavar="FORMAT or COMMAND", help='document output format (or command)')
    return parser


def setup_logger(verbose=1, color=True):
    """
    Set up a Python Logger instance for the Sage documentation builder.

    The optional argument sets logger's level and message format.
    """
    # Set up colors. Adapted from sphinx.cmdline.
    import sphinx.util.console as c
    if not color or not sys.stdout.isatty() or not c.color_terminal():
        c.nocolor()

    # Available colors: black, darkgray, (dark)red, dark(green),
    # brown, yellow, (dark)blue, purple, fuchsia, turquoise, teal,
    # lightgray, white.  Available styles: reset, bold, faint,
    # standout, underline, blink.

    # Set up log record formats.
    format_std = "%(message)s"
    formatter = logging.Formatter(format_std)

    # format_debug = "%(module)s #%(lineno)s %(funcName)s() %(message)s"
    fields = ['%(module)s', '#%(lineno)s', '%(funcName)s()', '%(message)s']
    colors = ['darkblue', 'darkred', 'brown', 'reset']
    styles = ['reset', 'reset', 'reset', 'reset']
    format_debug = ""
    for i in range(len(fields)):
        format_debug += c.colorize(styles[i], c.colorize(colors[i], fields[i]))
        if i != len(fields):
            format_debug += " "

    # Note: There's also Handler.setLevel().  The argument is the
    # lowest severity message that the respective logger or handler
    # will pass on.  The default levels are DEBUG, INFO, WARNING,
    # ERROR, and CRITICAL.  We use "WARNING" for normal verbosity and
    # "ERROR" for quiet operation.  It's possible to define custom
    # levels.  See the documentation for details.
    if verbose == 0:
        logger.setLevel(logging.ERROR)
    if verbose == 1:
        logger.setLevel(logging.WARNING)
    if verbose == 2:
        logger.setLevel(logging.INFO)
    if verbose == 3:
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter(format_debug)

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)


class IntersphinxCache:
    """
    Replace sphinx.ext.intersphinx.fetch_inventory by an in-memory
    cached version.
    """
    def __init__(self):
        self.inventories = {}
        self.real_fetch_inventory = sphinx.ext.intersphinx.fetch_inventory
        sphinx.ext.intersphinx.fetch_inventory = self.fetch_inventory

    def fetch_inventory(self, app, uri, inv):
        """
        Return the result of ``sphinx.ext.intersphinx.fetch_inventory()``
        from a cache if possible. Otherwise, call
        ``sphinx.ext.intersphinx.fetch_inventory()`` and cache the result.
        """
        t = (uri, inv)
        try:
            return self.inventories[t]
        except KeyError:
            i = self.real_fetch_inventory(app, uri, inv)
            self.inventories[t] = i
            return i


def main():
    # The options of SAGE_DOCBUILD_OPTS come before the ones written on the
    # command line, so that they are processed as options and overridden.  Note
    # that the args passed to parse_args() shouldn't include sys.argv[0].
    parser = setup_parser()
    args: BuildOptions = parser.parse_args(command_line_args()) # type: ignore

    # Check that the docs source directory exists.  Note that sage.env reads
    # the environment when it is imported, and that this function writes the
    # directories it settles on back into it, so nothing here may import it:
    # the extensions that do would see the values it is replacing.
    if args.source_dir is None:
        args.source_dir = Path(os.environ.get('SAGE_DOC_SRC', 'src/doc'))
    args.source_dir = args.source_dir.absolute()
    if not _path_is_directory(args.source_dir):
        parser.error(f"Source directory {args.source_dir} does not exist.")

    if args.all_documents:
        if args.all_documents == 'reference':
            docs = get_all_reference_documents(args.source_dir / 'en')
        elif args.all_documents == 'all':
            docs = get_all_documents(args.source_dir)
        else:
            parser.error(f"Unknown argument {args.all_documents} for --all-documents.")
        for d in docs:
            print(d.as_posix())
        sys.exit(0)

    # Check that the docs output directory exists.  This is also where a build
    # reads the inventories and the bibliography of the other documents from,
    # so a build redirected with -o resolves nothing against them.
    args.output_dir_given = args.output_dir is not None
    if args.output_dir is None:
        args.output_dir = Path(os.environ.get('SAGE_DOC', 'src/doc'))
    args.output_dir = args.output_dir.absolute()
    output_info = _path_status(args.output_dir)
    if output_info is not None and not stat.S_ISDIR(output_info.st_mode):
        parser.error(f"Output directory {args.output_dir} is not a directory.")

    # Get the name and type (target format) of the document we are
    # trying to build.
    name, typ = args.document, args.format
    if not name or not typ:
        parser.print_help()
        sys.exit(1)
    if name.startswith('file=') and typ not in _output_formats():
        parser.error(f"Unknown single-file output format {typ!r}.")

    # Set up module-wide logging.
    setup_logger(args.verbose, args.color)

    def excepthook(*exc_info):
        logger.error('Error building the documentation.', exc_info=exc_info)
        logger.info('''
Note: incremental documentation builds sometimes cause spurious
error messages. To be certain that these are real errors, run
"make doc-clean doc-uninstall" first and try again.''')

    sys.excepthook = excepthook

    # Set up the environment based on the command-line options
    if args.check_nested:
        os.environ['SAGE_CHECK_NESTED'] = 'True'
    if args.underscore:
        os.environ['SAGE_DOC_UNDERSCORE'] = "True"
    if args.sphinx_opts:
        # The option is documented as a comma-separated list, and the shell has
        # parsed the value already: a comma is the only separator left. Cutting
        # on whitespace too would break -S=-D,html_title=Sage Reference Manual,
        # and shlex.split() would eat the backslashes of
        # -S=-Dlatex_elements.preamble=\usepackage{microtype}.
        build_options.ALLSPHINXOPTS += [opt for opt in
                                        (piece.strip() for piece
                                         in args.sphinx_opts.split(','))
                                        if opt]
    if args.no_pdf_links:
        build_options.WEBSITESPHINXOPTS = ['-A', 'hide_pdf_links=1']
    if args.warn_links:
        build_options.WARN_LINKS = True
    if args.no_plot:
        os.environ['SAGE_SKIP_PLOT_DIRECTIVE'] = 'yes'
    if args.live_doc:
        os.environ['SAGE_LIVE_DOC'] = 'yes'
    if args.skip_tests:
        os.environ['SAGE_SKIP_TESTS_BLOCKS'] = 'True'
    if args.use_cdns:
        os.environ['SAGE_USE_CDNS'] = 'yes'
    os.environ['SAGE_DOC_SRC'] = str(args.source_dir)
    os.environ['SAGE_DOC'] = str(args.output_dir)

    build_options.ABORT_ON_ERROR = not args.keep_going

    # Set up Intersphinx cache
    _ = IntersphinxCache()

    # get_builder() reports a document that does not exist and stops, so it
    # comes first: writing hundreds of generated files into a tree only to
    # refuse to build anything from it would leave a mess behind.  A single
    # file is documented on its own and needs none of them.
    builder = get_builder(name, args)
    build = getattr(builder, typ)
    if not callable(build):
        raise AttributeError(
            f"{type(builder).__name__!r} object has no command {typ!r}")
    # A SingleFileBuilder prepares its output tree under the requested output
    # directory while it is constructed, so it may have created the directory
    # since the validation above.
    if _path_status(args.output_dir) is None:
        try:
            args.output_dir.mkdir(parents=True)
        except Exception as error:
            parser.error(
                f"Failed to create output directory {args.output_dir}: {error}")
    if not name.startswith('file='):
        generate_doc_sources(args.source_dir)

    if not args.no_prune_empty_dirs:
        # Delete empty directories. This is needed in particular for empty
        # directories due to "git checkout" which never deletes empty
        # directories it leaves behind. See Issue #20010.
        # Issue #31948: This is not parallelization-safe; use the option
        # --no-prune-empty-dirs to turn it off
        for dirpath, dirnames, filenames in os.walk(args.source_dir, topdown=False):
            if not dirnames + filenames:
                logger.warning('Deleting empty directory {0}'.format(dirpath))
                os.rmdir(dirpath)

    import sage.all  # TODO: Remove once all modules can be imported independently  # noqa: F401

    build()


if __name__ == '__main__':
    sys.exit(main())
