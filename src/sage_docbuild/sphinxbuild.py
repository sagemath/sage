# sage.doctest: needs sphinx
r"""
Sphinx build script

This is Sage's version of the ``sphinx-build`` script. We redirect ``stdout`` and
``stderr`` to our own logger, and remove some unwanted chatter.
"""
# ****************************************************************************
#       Copyright (C) 2013-2014 Volker Braun <vbraun.name@gmail.com>
#                     2013-2017 J. H. Palmieri <<palmieri@math.washington.edu>
#                     2013-2017 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#                     2014      Christopher Schwan <cschwan@students.uni-mainz.de>
#                     2014      Nicolas M. Thiéry <nthiery@users.sf.net>
#                     2015      Marc Mezzarobba <marc@mezzarobba.net>
#                     2015      André Apitzsch <andre.apitzsch@etit.tu-chemnitz.de>
#                     2018      Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import ast
import io
import os
import re
import stat
import sys
import tokenize
from pathlib import Path

import sphinx
import sphinx.cmd.build

from sage.env import SAGE_DOC_SRC, SAGE_SRC


# override the fancy multi-line formatting
def term_width_line(text):
    return text + '\n'


sphinx.util.console.term_width_line = term_width_line


# https://en.wikipedia.org/wiki/ANSI_escape_code
ANSI_ESCAPE_SEQUENCE = re.compile(r'''
    \x1b    # ESC
    \[      # CSI sequence starts
    [0-?]*  # parameter bytes
    [ -/]*  # intermediate bytes
    [@-~]   # final byte
    ''', re.VERBOSE)


# Warnings that do not indicate a problem with Sage's documentation.  Sphinx
# emits each of these without a warning type, so :confval:`suppress_warnings`
# cannot be used and we have to match the message itself.
IGNORED_WARNINGS = (
    re.compile("WARNING: favicon file 'favicon.ico' does not exist"),
    re.compile('WARNING: html_static_path entry .* does not exist'),
    re.compile('WARNING: while setting up extension'),
    re.compile('WARNING: .* is not referenced'),
    # A single-file build does not load the multidocs extension; see #29651.
    re.compile(r"WARNING: unknown config value 'multidoc_first_pass'"),
)

# A first pass runs before the inventories of the other documents exist, so
# references into them cannot be resolved yet.
FIRST_PASS_IGNORED_WARNINGS = (
    re.compile('WARNING: undefined label'),
)

# An inventory build writes no output: it resolves no citation and keeps no
# search index.
INVENTORY_IGNORED_WARNINGS = (
    re.compile('WARNING: citation not found:'),
    re.compile("WARNING: search index couldn't be loaded"),
)


# This is only a cheap prefilter.  Whether the matching line really defines a
# citation (rather than displaying one in a literal block, for example) is
# decided by the reStructuredText parser below.
CITATION_CANDIDATE = re.compile(rb'(?<!\.)\.\.[ \t]+\[')
CITATION_CANDIDATE_TEXT = re.compile(r'^[ \t]*\.\.[ \t]+\[', re.M)
INCLUDE_CANDIDATE = re.compile(
    rb'^[ \t]*\.\.[ \t]+include::', re.IGNORECASE | re.MULTILINE)


def _citation_names(text: str, source: str = '<citation scan>', *,
                    follow_includes: bool = False) -> set[str]:
    r"""
    Return the citation names that *text* actually defines.

    Parsing the text matters: a line that looks like a definition inside a
    literal block is an example, not a citation of the document containing it.

    EXAMPLES::

        sage: from sage_docbuild.sphinxbuild import _citation_names
        sage: text = ('Example::\n\n'
        ....:         '   .. [NotAReference] only code\n\n'
        ....:         '.. [Real] a reference\n')
        sage: _citation_names(text)
        {'Real'}
    """
    from docutils import nodes
    from docutils.frontend import get_default_settings
    from docutils.parsers.rst import Parser
    from docutils.utils import new_document

    settings = get_default_settings(Parser)
    settings.report_level = 5
    settings.halt_level = 6
    settings.warning_stream = io.StringIO()
    settings.file_insertion_enabled = follow_includes
    settings.raw_enabled = False
    document = new_document(source, settings)
    try:
        Parser().parse(text, document)
    except Exception:
        # A source that docutils cannot parse on its own cannot safely prove
        # that a warning is known.  Leaving its label out reports the warning;
        # guessing with a regex would hide it.
        return set()
    return {citation[0].astext()
            for citation in document.findall(nodes.citation)
            if citation.children}


def _python_docstrings(text: str, filename: str):
    """Yield the real docstrings of a Python source file."""
    try:
        tree = ast.parse(text, filename=filename)
    except (SyntaxError, ValueError):
        return
    documented = (ast.Module, ast.ClassDef, ast.FunctionDef,
                  ast.AsyncFunctionDef)
    for node in ast.walk(tree):
        if not isinstance(node, documented) or not node.body:
            continue
        first = node.body[0]
        if (isinstance(first, ast.Expr)
                and isinstance(first.value, ast.Constant)
                and isinstance(first.value.value, str)):
            yield first.value.value


def _cython_docstrings(text: str):
    """Yield module, class, function, and property docstrings from Cython."""
    # Python's tokenizer understands Cython string literals and indentation.
    # Keep track of which suites can own a docstring: accepting the first
    # string in every suite would mistake an ``if`` body for documentation.
    expects_docstring = [True]
    statement = []
    pending_suite = False
    docstring_parts = []

    def docstring_suite(tokens):
        words = [token.string for token in tokens]
        if ('def' in words or 'class' in words or 'cppclass' in words
                or words[:1] == ['property']):
            return True
        if words[:1] == ['cpdef']:
            return '(' in words
        return (words[:1] == ['cdef']
                and 'extern' not in words
                and '(' in words)

    try:
        tokens = tokenize.generate_tokens(io.StringIO(text).readline)
        for token in tokens:
            if token.type in (tokenize.ENCODING, tokenize.COMMENT, tokenize.NL):
                continue
            if token.type == tokenize.INDENT:
                expects_docstring.append(pending_suite)
                pending_suite = False
                statement.clear()
                continue
            if token.type == tokenize.DEDENT:
                if len(expects_docstring) > 1:
                    expects_docstring.pop()
                continue
            if token.type == tokenize.NEWLINE:
                if docstring_parts:
                    yield ''.join(docstring_parts)
                    docstring_parts.clear()
                pending_suite = bool(statement) and docstring_suite(statement)
                statement.clear()
                continue
            if token.type == tokenize.ENDMARKER:
                if docstring_parts:
                    yield ''.join(docstring_parts)
                break

            at_start = not statement
            statement.append(token)
            if at_start and token.type == tokenize.STRING and expects_docstring[-1]:
                try:
                    value = ast.literal_eval(token.string)
                except (SyntaxError, ValueError):
                    value = None
                if isinstance(value, str):
                    docstring_parts.append(value)
            elif (docstring_parts and token.type == tokenize.STRING
                  and all(part.type == tokenize.STRING
                          for part in statement[:-1])):
                try:
                    value = ast.literal_eval(token.string)
                except (SyntaxError, ValueError):
                    value = None
                if isinstance(value, str):
                    docstring_parts.append(value)
            expects_docstring[-1] = False
    except (IndentationError, tokenize.TokenError):
        return


def _walk_citation_tree(root: Path):
    """Yield a source tree without hiding traversal permission errors."""

    def raise_error(error):
        raise error

    for directory, subdirectories, filenames in os.walk(
            root, onerror=raise_error):
        subdirectories[:] = [name for name in subdirectories
                             if name != '__pycache__']
        base = Path(directory)
        yield from (base / name for name in filenames)


def _citation_path_status(path: Path):
    """Return a followed status, suppressing absence but no other error."""

    try:
        path.lstat()
    except FileNotFoundError:
        return None
    try:
        return path.stat()
    except FileNotFoundError as error:
        raise OSError(f'broken or moving citation source path: {path}') from error


def _module_citation_labels(source_roots, source_files=()) -> set[str]:
    """Return citation labels defined by docstrings in the given sources."""
    docstrings = []

    def collect(path):
        if path.suffix not in ('.py', '.pyx', '.pxd'):
            return
        try:
            encoded = path.read_bytes()
        except FileNotFoundError:
            return
        if CITATION_CANDIDATE.search(encoded) is None:
            return
        try:
            with tokenize.open(path) as source_file:
                text = source_file.read()
        except FileNotFoundError:
            return
        except (SyntaxError, UnicodeError):
            return
        strings = (_python_docstrings(text, str(path))
                   if path.suffix == '.py' else _cython_docstrings(text))
        docstrings.extend(docstring for docstring in strings
                          if CITATION_CANDIDATE_TEXT.search(docstring))

    for path in source_files:
        collect(path)
    for root in source_roots:
        for path in _walk_citation_tree(root):
            collect(path)
    # Parsing once is much faster than constructing one docutils document per
    # docstring.  A top-level transition separates the independent strings and
    # closes any literal block at the end of the previous one.
    return _citation_names('\n\n----\n\n'.join(docstrings))


def _documentation_citation_labels(doc_roots) -> tuple[set[str], bool]:
    """Return citation labels in the documentation and whether a bibliography exists."""
    labels: set[str] = set()
    have_bibliography = False
    for root in doc_roots:
        bibliography = root / 'en' / 'reference' / 'references' / 'index.rst'
        bibliography_info = _citation_path_status(bibliography)
        if bibliography_info is not None:
            if not stat.S_ISREG(bibliography_info.st_mode):
                raise OSError(
                    f'citation bibliography is not a regular file: '
                    f'{bibliography}')
            have_bibliography = True
        root_info = _citation_path_status(root)
        if root_info is None:
            continue
        if not stat.S_ISDIR(root_info.st_mode):
            raise OSError(f'documentation source is not a directory: {root}')
        for path in _walk_citation_tree(root):
            if path.suffix != '.rst':
                continue
            try:
                encoded = path.read_bytes()
            except FileNotFoundError:
                continue
            # An include may carry definitions even when this file has no
            # citation-looking line of its own.  literalinclude is deliberately
            # not special-cased: docutils keeps its contents literal.
            if (CITATION_CANDIDATE.search(encoded) is None
                    and INCLUDE_CANDIDATE.search(encoded) is None):
                continue
            try:
                text = encoded.decode('utf-8')
            except UnicodeDecodeError:
                text = encoded.decode('utf-8', 'replace')
            labels.update(_citation_names(
                text, str(path), follow_includes=True))
    return labels, have_bibliography


def citation_labels(single_file_path=None, single_file_source_root=None):
    r"""
    Return every citation label that the sources of Sage define, or ``None``
    when the bibliography of the reference manual is not among the sources.

    Most citations are collected in a reference manual, and the rest are
    written in the docstring defining them.  Both are parsed as reStructuredText
    so that examples and comments are not mistaken for definitions.  For a file
    of a second checkout, that checkout is searched as well as the active one.

    No answer is cached: a docbuild process can build more than one file, edit
    a source between builds, or switch the source root in its environment.

    EXAMPLES::

        sage: from sage_docbuild.sphinxbuild import citation_labels
        sage: labels = citation_labels()
        sage: labels is None or 'AB2007' in labels
        True
    """
    def distinct(paths):
        answer = []
        seen = set()
        for candidate in paths:
            if candidate is None:
                continue
            path = Path(candidate).absolute()
            try:
                identity = path.resolve()
            except OSError:
                identity = path
            info = _citation_path_status(path)
            if info is None:
                continue
            if identity not in seen and stat.S_ISDIR(info.st_mode):
                seen.add(identity)
                answer.append(path)
        return answer

    source_roots = [Path(os.environ.get('SAGE_SRC') or SAGE_SRC)]
    doc_roots = [Path(os.environ.get('SAGE_DOC_SRC') or SAGE_DOC_SRC)]
    source_files = []
    if single_file_path is not None and single_file_source_root is not None:
        file = Path(single_file_path).absolute()
        root = Path(single_file_source_root).absolute()
        try:
            relative = file.relative_to(root)
        except ValueError:
            relative = None
        # A one-component path is a standalone module.  Its parent is not a
        # Python source tree to walk (and is often a large temporary directory).
        if relative is not None and len(relative.parts) > 1:
            # The import root may be a home or temporary directory containing
            # unrelated projects.  Only the target's top-level package is a
            # source boundary that can be inferred from the file itself.
            try:
                target_in_known_root = any(
                    file.resolve().is_relative_to(known.resolve())
                    for known in source_roots)
            except OSError:
                target_in_known_root = any(
                    file.is_relative_to(known) for known in source_roots)
            if not target_in_known_root:
                source_roots.append(root / relative.parts[0])
            sage_info = _citation_path_status(root / 'sage')
            docbuild_info = _citation_path_status(root / 'sage_docbuild')
            configure_info = _citation_path_status(
                root.parent / 'configure.ac')
            if (sage_info is not None and stat.S_ISDIR(sage_info.st_mode)
                    and ((docbuild_info is not None
                          and stat.S_ISDIR(docbuild_info.st_mode))
                         or (configure_info is not None
                             and stat.S_ISREG(configure_info.st_mode)))):
                doc_roots.append(root / 'doc')
        else:
            source_files.append(file)

    source_roots = distinct(source_roots)
    # An explicit --source may be a small custom manual, but it must not hide
    # the bibliography of the Sage source tree that supplies the documented
    # modules.  A checkout conventionally keeps it beside ``SAGE_SRC``.
    doc_roots.extend(root / 'doc' for root in source_roots)
    doc_roots = distinct(doc_roots)
    labels, have_bibliography = _documentation_citation_labels(doc_roots)
    if not have_bibliography:
        # Without the bibliography, the citations that would be found are a
        # small part of the ones that Sage defines, and reporting the rest as
        # unknown would be worse than reporting none of them.
        return None
    labels.update(_module_citation_labels(source_roots, source_files))
    return frozenset(labels)


class KnownCitation:
    r"""
    Match the warning about a citation that Sage defines.

    A single file is documented without the multidocs extension, hence without
    the bibliography of the reference manual, so a citation of it - and most
    Sage files cite something - has nothing to resolve against; see #29651.
    Ignoring every such warning would hide the citations that are misspelled,
    which are the ones worth reporting, so the labels that the sources define
    are consulted instead; see :func:`citation_labels`.

    Instances answer to ``search()``, as the compiled patterns they are used
    among do.

    EXAMPLES::

        sage: from sage_docbuild.sphinxbuild import KnownCitation
        sage: known = KnownCitation({'Cohen1996'})
        sage: bool(known.search('WARNING: citation not found: Cohen1996'))
        True
        sage: bool(known.search('WARNING: citation not found: Cohen1997'))
        False
        sage: bool(known.search('WARNING: citation not found: cohen1996'))
        True

    Sphinx colors what it writes to a terminal, and the label is read past the
    escape sequences::

        sage: bool(known.search(
        ....:     '\x1b[91mWARNING: citation not found: Cohen1996\x1b[39;49;00m'))
        True
    """
    _pattern = re.compile(r'WARNING: citation not found: (\S+)')

    def __init__(self, labels):
        from docutils.nodes import fully_normalize_name

        self._labels = frozenset(map(fully_normalize_name, labels))

    def search(self, line):
        match = self._pattern.search(ANSI_ESCAPE_SEQUENCE.sub('', line))
        if match is None:
            return None
        from docutils.nodes import fully_normalize_name

        if fully_normalize_name(match.group(1)) in self._labels:
            return match
        return None


def single_file_ignored_warnings(single_file_path=None,
                                 single_file_source_root=None):
    r"""
    Return the warning patterns that a file documented on its own is built with.

    Its citations are the ones that :class:`KnownCitation` accounts for.  Where
    the sources carry no bibliography to tell one citation from another, all of
    them are ignored.

    EXAMPLES::

        sage: from sage_docbuild.sphinxbuild import single_file_ignored_warnings
        sage: all(hasattr(pattern, 'search')
        ....:     for pattern in single_file_ignored_warnings())
        True
    """
    labels = citation_labels(single_file_path, single_file_source_root)
    if labels is None:
        return (re.compile('WARNING: citation not found:'),)
    return (KnownCitation(labels),)


class SageSphinxLogger():
    r"""
    This implements the file object interface to serve as
    ``sys.stdout``/``sys.stderr`` replacement.
    """
    ansi_escape_sequence = ANSI_ESCAPE_SEQUENCE
    ansi_escape_sequence_color = re.compile(r'''
        \x1b    # ESC
        \[      # CSI sequence starts
        [0-9;]* # parameter bytes
                # intermediate bytes
        m       # final byte
        ''', re.VERBOSE)

    prefix_len = 9

    def __init__(self, stream, prefix, *, warnings_are_errors=False,
                 ignored_warnings=()):
        self._init_chatter(warnings_are_errors, ignored_warnings)
        self._stream = stream
        self._color = stream.isatty()
        prefix = prefix[0:self.prefix_len]
        prefix = ('[{0:' + str(self.prefix_len) + '}]').format(prefix)
        self._is_stdout = (stream.fileno() == 1)
        self._is_stderr = (stream.fileno() == 2)
        if self._is_stdout:
            color = 'darkgreen'
        elif self._is_stderr:
            color = 'red'
        else:
            color = 'lightgray'
        self._prefix = sphinx.util.console.colorize(color, prefix)
        # When we see an error in the log, we store it here and raise it at the
        # end of the file (sometimes the lines following the error still
        # contain valuable information.)
        self._error = None

    def _init_chatter(self, warnings_are_errors, ignored_warnings):
        # We drop any messages from the output that match these regular
        # expressions. These just bloat the output and do not contain any
        # information that we care about.
        self._useless_chatter = (
            re.compile(r'^$'),
            re.compile(r'^Running Sphinx'),
            re.compile(r'^updating environment: 0 added, 0 changed, 0 removed'),
            re.compile(r'^building \[.*\]: targets for 0 source files that are out of date'),
            re.compile(r'^building \[.*\]: targets for 0 po files that are out of date'),
            re.compile(r'^building \[.*\]: targets for 0 mo files that are out of date'),
            re.compile(r'^build succeeded'),  # We still have "Build finished."
            re.compile(r'^Saved pickle file: citations\.pickle'),
            re.compile(r'^Compiling|Copying|Merging|Writing'),
            re.compile(r'^compiling|copying|checking|dumping|executing|generating|linking|loading|looking|pickling|preparing|reading|writing'),
            re.compile(r'done'),
            re.compile(r'^WARNING:$'),
        )

        self._ignored_warnings = tuple(ignored_warnings)
        self._useless_chatter += self._ignored_warnings

        # replacements: pairs of regular expressions and their replacements,
        # to be applied to Sphinx output.
        self.replacements = [(re.compile('build succeeded, [0-9]+ warning[s]?.'),
                              'build succeeded.')]

        # Diagnostics that make the build fail. Sphinx reports a failed build
        # through its exit status, which :func:`runsphinx` checks; these
        # patterns catch the message levels that Sphinx and docutils emit
        # without failing the build by themselves.
        self._error_patterns = (re.compile('SEVERE'), re.compile('ERROR'))
        if warnings_are_errors:
            self._error_patterns += (re.compile('WARNING:'),)

    def _filter_out(self, line):
        if self._error is not None and self._is_stdout:
            # swallow non-errors after an error occurred
            return True
        line = re.sub(self.ansi_escape_sequence, '', line)
        line = line.strip()
        for regex in self._useless_chatter:
            if regex.search(line) is not None:
                return True
        return False

    def _check_errors(self, line):
        r"""
        Search for errors in line.

        EXAMPLES::

            sage: from sys import stdout
            sage: from sage_docbuild.sphinxbuild import SageSphinxLogger
            sage: logger = SageSphinxLogger(stdout, "doctesting")
            sage: logger._log_line("ERROR: something went wrong\n") # indirect doctest
            [doctestin] ERROR: something went wrong
            sage: logger.raise_errors()
            Traceback (most recent call last):
            ...
            OSError: ERROR: something went wrong

        A warning is an error only where the build type says so::

            sage: logger = SageSphinxLogger(stdout, "doctesting")
            sage: logger._log_line("WARNING: undefined label: 'foo'\n")
            [doctestin] WARNING: undefined label: 'foo'
            sage: logger.raise_errors()

            sage: logger = SageSphinxLogger(stdout, "doctesting",
            ....:                           warnings_are_errors=True)
            sage: logger._log_line("WARNING: undefined label: 'foo'\n")
            [doctestin] WARNING: undefined label: 'foo'
            sage: logger.raise_errors()
            Traceback (most recent call last):
            ...
            OSError: WARNING: undefined label: 'foo'

        Ignored warnings are neither errors nor printed::

            sage: from sage_docbuild.sphinxbuild import FIRST_PASS_IGNORED_WARNINGS
            sage: logger = SageSphinxLogger(stdout, "doctesting",
            ....:                           warnings_are_errors=True,
            ....:                           ignored_warnings=FIRST_PASS_IGNORED_WARNINGS)
            sage: logger._log_line("WARNING: undefined label: 'foo'\n")
            sage: logger.raise_errors()
        """
        if self._error is not None:
            return  # we already have found an error
        for error in self._error_patterns:
            if error.search(line) is not None:
                for ignored in self._ignored_warnings:
                    if ignored.search(line) is not None:
                        break
                else:
                    self._error = line
                    return

    def _log_line(self, line):
        r"""
        Write ``line`` to the output stream with some mangling.

        EXAMPLES::

            sage: from sys import stdout
            sage: from sage_docbuild.sphinxbuild import SageSphinxLogger
            sage: logger = SageSphinxLogger(stdout, "doctesting")
            sage: logger._log_line("building documentation…\n")
            [doctestin] building documentation…

        TESTS:

        Verify that :issue:`25160` has been resolved::

            sage: logger = SageSphinxLogger(stdout, "#25160")
            sage: import traceback
            sage: try:
            ....:     raise Exception("artificial exception")
            ....: except Exception:
            ....:     for line in traceback.format_exc().split('\n'):
            ....:         logger._log_line(line)
            [#25160   ] Traceback (most recent call last):
            [#25160   ]   File ...
            [#25160   ]     raise Exception("artificial exception")
            [#25160   ] Exception: artificial exception
        """
        skip_this_line = self._filter_out(line)
        self._check_errors(line)
        for (old, new) in self.replacements:
            line = old.sub(new, line)
        line = self._prefix + ' ' + line.rstrip() + '\n'
        if not self._color:
            line = self.ansi_escape_sequence_color.sub('', line)
        if not skip_this_line:
            # sphinx does produce messages in the current locals which
            # could be non-ascii
            # see https://github.com/sagemath/sage/issues/27706
            self._stream.write(line if isinstance(line, str) else line.encode('utf8'))
            self._stream.flush()

    def raise_errors(self):
        r"""
        Raise an exceptions if any errors have been found while parsing the
        Sphinx output.

        EXAMPLES::

            sage: from sys import stdout
            sage: from sage_docbuild.sphinxbuild import SageSphinxLogger
            sage: logger = SageSphinxLogger(stdout, "doctesting")
            sage: logger._log_line("This is a SEVERE error\n")
            [doctestin] This is a SEVERE error
            sage: logger.raise_errors()
            Traceback (most recent call last):
            ...
            OSError: This is a SEVERE error

        """
        if self._error is not None:
            raise OSError(self._error)

    _line_buffer = ''

    def _write(self, string):
        self._line_buffer += string
        lines = self._line_buffer.splitlines()
        for i, line in enumerate(lines):
            last = (i == len(lines) - 1)
            if last and not self._line_buffer.endswith('\n'):
                self._line_buffer = line
                return
            self._log_line(line)
        self._line_buffer = ''

    # file object interface follows

    closed = False
    encoding = None
    mode = 'w'
    name = '<log>'
    newlines = None
    softspace = 0

    def isatty(self):
        return True

    def close(self):
        if self._line_buffer != '':
            self._log_line(self._line_buffer)
            self._line_buffer = ''

    def flush(self):
        self._stream.flush()

    def write(self, str):
        try:
            self._write(str)
        except OSError:
            raise
        except Exception:
            import traceback
            traceback.print_exc(file=self._stream)

    def writelines(self, sequence):
        for line in sequence:
            self.write(line)


def runsphinx(argv, *, prefix=None, warnings_are_errors=True,
              first_pass=False, is_inventory=False, single_file=False,
              single_file_path=None, single_file_source_root=None):
    r"""
    Run ``sphinx-build`` with the arguments ``argv``, logging its output.

    INPUT:

    - ``argv`` -- list of arguments for ``sphinx-build``, without the program
      name; the last one is the output directory

    - ``prefix`` -- string used to tag the output lines; defaults to the name
      of the output directory

    - ``warnings_are_errors`` -- whether a warning makes the build fail. The
      LaTeX builder warns about markup that it cannot translate, which is not
      a reason to stop.

    - ``first_pass`` -- whether this build runs before the inventories of the
      other documents have been written

    - ``is_inventory`` -- whether this build only writes an inventory

    - ``single_file`` -- whether this build documents one file of its own,
      outside of any manual

    - ``single_file_path`` -- the file being documented, when ``single_file``
      is true

    - ``single_file_source_root`` -- the import root of ``single_file_path``

    A failed build raises :exc:`OSError`, either because Sphinx reported a
    nonzero exit status, or because it logged a diagnostic that Sage treats as
    an error.
    """
    ignored_warnings = IGNORED_WARNINGS
    if first_pass:
        ignored_warnings += FIRST_PASS_IGNORED_WARNINGS
    if single_file:
        ignored_warnings += single_file_ignored_warnings(
            single_file_path, single_file_source_root)
    if is_inventory:
        ignored_warnings += INVENTORY_IGNORED_WARNINGS

    if prefix is None:
        prefix = os.path.basename(argv[-1])

    saved_stdout = sys.stdout
    saved_stderr = sys.stderr

    original_filters = None
    if not sys.warnoptions:
        import warnings
        original_filters = warnings.filters[:]
        warnings.filterwarnings("ignore", category=DeprecationWarning, module='sphinx.util.inspect')

    try:
        sys.stdout = SageSphinxLogger(sys.stdout, prefix,
                                      warnings_are_errors=warnings_are_errors,
                                      ignored_warnings=ignored_warnings)
        sys.stderr = SageSphinxLogger(sys.stderr, prefix,
                                      warnings_are_errors=warnings_are_errors,
                                      ignored_warnings=ignored_warnings)
        # Note that this call as of early 2018 leaks memory. So make sure that
        # you don't call runsphinx() several times in a row. (i.e., you want to
        # fork() somewhere before this call.)
        # We don't use subprocess here, as we don't want to re-initialize Sage
        # for every docbuild as this takes a while.
        status = sphinx.cmd.build.build_main(argv)
        # Report a logged diagnostic in preference to the exit status: it says
        # what went wrong, while the status only says that something did.
        sys.stderr.raise_errors()
        sys.stdout.raise_errors()
        if status:
            raise OSError(f"sphinx-build exited with status {status}")
    finally:
        sys.stdout = saved_stdout
        sys.stderr = saved_stderr
        sys.stdout.flush()
        sys.stderr.flush()
        # A failed build must not leave the filter behind: the next build in
        # this process would silently drop the warnings it matches.
        if original_filters is not None:
            warnings.filters = original_filters[:]
