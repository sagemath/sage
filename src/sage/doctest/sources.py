"""
Classes for sources of doctests

This module defines various classes for sources from which doctests
originate, such as files, functions or database entries.

AUTHORS:

- David Roe (2012-03-27) -- initial version, based on Robert Bradshaw's code.
"""

# ****************************************************************************
#       Copyright (C) 2012 David Roe <roed.math@gmail.com>
#                          Robert Bradshaw <robertwb@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os
import sys
import re
import random
import doctest
from Cython.Utils import is_package_dir
from sage.cpython.string import bytes_to_str
from sage.repl.load import load
from sage.misc.lazy_attribute import lazy_attribute
from .parsing import SageDocTestParser
from .util import NestedName
from sage.structure.dynamic_class import dynamic_class
from sage.env import SAGE_SRC, SAGE_LIB

# Python file parsing
triple_quotes = re.compile(r"\s*[rRuU]*((''')|(\"\"\"))")
name_regex = re.compile(r".*\s(\w+)([(].*)?:")

# LaTeX file parsing
begin_verb = re.compile(r"\s*\\begin{verbatim}")
end_verb = re.compile(r"\s*\\end{verbatim}\s*(%link)?")
begin_lstli = re.compile(r"\s*\\begin{lstlisting}")
end_lstli = re.compile(r"\s*\\end{lstlisting}\s*(%link)?")
skip = re.compile(r".*%skip.*")

# ReST file parsing
link_all = re.compile(r"^\s*\.\.\s+linkall\s*$")
double_colon = re.compile(r"^(\s*).*::\s*$")
code_block = re.compile(r"^(\s*)[.][.]\s*code-block\s*::.*$")

whitespace = re.compile(r"\s*")
bitness_marker = re.compile('#.*(32|64)-bit')
bitness_value = '64' if sys.maxsize > (1 << 32) else '32'

# For neutralizing doctests
find_prompt = re.compile(r"^(\s*)(>>>|sage:)(.*)")

# For testing that enough doctests are created
sagestart = re.compile(r"^\s*(>>> |sage: )\s*[^#\s]")
untested = re.compile("(not implemented|not tested)")

# For parsing a PEP 0263 encoding declaration
pep_0263 = re.compile(br'^[ \t\v]*#.*?coding[:=]\s*([-\w.]+)')

# Source line number in warning output
doctest_line_number = re.compile(r"^\s*doctest:[0-9]")


def get_basename(path):
    """
    This function returns the basename of the given path, e.g. sage.doctest.sources or doc.ru.tutorial.tour_advanced

    EXAMPLES::

        sage: from sage.doctest.sources import get_basename
        sage: from sage.env import SAGE_SRC
        sage: import os
        sage: get_basename(os.path.join(SAGE_SRC,'sage','doctest','sources.py'))
        'sage.doctest.sources'
    """
    if path is None:
        return None
    if not os.path.exists(path):
        return path
    path = os.path.abspath(path)
    root = os.path.dirname(path)
    # If the file is in the sage library, we can use our knowledge of
    # the directory structure
    dev = SAGE_SRC
    sp = SAGE_LIB
    if path.startswith(dev):
        # there will be a branch name
        i = path.find(os.path.sep, len(dev))
        if i == -1:
            # this source is the whole library....
            return path
        root = path[:i]
    elif path.startswith(sp):
        root = path[:len(sp)]
    else:
        # If this file is in some python package we can see how deep
        # it goes by the presence of __init__.py files.
        while os.path.exists(os.path.join(root, '__init__.py')):
            root = os.path.dirname(root)
    fully_qualified_path = os.path.splitext(path[len(root) + 1:])[0]
    if os.path.split(path)[1] == '__init__.py':
        fully_qualified_path = fully_qualified_path[:-9]
    return fully_qualified_path.replace(os.path.sep, '.')


class DocTestSource(object):
    """
    This class provides a common base class for different sources of doctests.

    INPUT:

    - ``options`` -- a :class:`sage.doctest.control.DocTestDefaults`
      instance or equivalent.
    """
    def __init__(self, options):
        """
        Initialization.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','sources.py')
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: TestSuite(FDS).run()
        """
        self.options = options

    def __eq__(self, other):
        """
        Comparison is just by comparison of attributes.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','sources.py')
            sage: DD = DocTestDefaults()
            sage: FDS = FileDocTestSource(filename,DD)
            sage: FDS2 = FileDocTestSource(filename,DD)
            sage: FDS == FDS2
            True
        """
        if type(self) != type(other):
            return False
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        """
        Test for non-equality.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','sources.py')
            sage: DD = DocTestDefaults()
            sage: FDS = FileDocTestSource(filename,DD)
            sage: FDS2 = FileDocTestSource(filename,DD)
            sage: FDS != FDS2
            False
        """
        return not (self == other)

    def _process_doc(self, doctests, doc, namespace, start):
        """
        Appends doctests defined in ``doc`` to the list ``doctests``.

        This function is called when a docstring block is completed
        (either by ending a triple quoted string in a Python file,
        unindenting from a comment block in a ReST file, or ending a
        verbatim or lstlisting environment in a LaTeX file).

        INPUT:

        - ``doctests`` -- a running list of doctests to which the new
          test(s) will be appended.

        - ``doc`` -- a list of lines of a docstring, each including
          the trailing newline.

        - ``namespace`` -- a dictionary or
          :class:`sage.doctest.util.RecordingDict`, used in the
          creation of new :class:`doctest.DocTest`s.

        - ``start`` -- an integer, giving the line number of the start
          of this docstring in the larger file.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.parsing import SageDocTestParser
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','util.py')
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: doctests, _ = FDS.create_doctests({})
            sage: manual_doctests = []
            sage: for dt in doctests:
            ....:     FDS.qualified_name = dt.name
            ....:     FDS._process_doc(manual_doctests, dt.docstring, {}, dt.lineno-1)
            sage: doctests == manual_doctests
            True
        """
        docstring = "".join(doc)
        new_doctests = self.parse_docstring(docstring, namespace, start)
        sig_on_count_doc_doctest = "sig_on_count() # check sig_on/off pairings (virtual doctest)\n"
        for dt in new_doctests:
            if len(dt.examples) > 0 and not (hasattr(dt.examples[-1],'sage_source')
                                             and dt.examples[-1].sage_source == sig_on_count_doc_doctest):
                # Line number refers to the end of the docstring
                sigon = doctest.Example(sig_on_count_doc_doctest, "0\n", lineno=docstring.count("\n"))
                sigon.sage_source = sig_on_count_doc_doctest
                dt.examples.append(sigon)
            doctests.append(dt)

    def _create_doctests(self, namespace, tab_okay=None):
        """
        Creates a list of doctests defined in this source.

        This function collects functionality common to file and string
        sources, and is called by
        :meth:`FileDocTestSource.create_doctests`.

        INPUT:

        - ``namespace`` -- a dictionary or
          :class:`sage.doctest.util.RecordingDict`, used in the
          creation of new :class:`doctest.DocTest`s.

        - ``tab_okay`` -- whether tabs are allowed in this source.

        OUTPUT:

        - ``doctests`` -- a list of doctests defined by this source

        - ``extras`` -- a dictionary with ``extras['tab']`` either
          False or a list of linenumbers on which tabs appear.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.util import NestedName
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','sources.py')
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: FDS.qualified_name = NestedName('sage.doctest.sources')
            sage: doctests, extras = FDS._create_doctests({})
            sage: len(doctests)
            41
            sage: extras['tab']
            False
            sage: extras['line_number']
            False
        """
        if tab_okay is None:
            tab_okay = isinstance(self,TexSource)
        self._init()
        self.line_shift = 0
        self.parser = SageDocTestParser(self.options.optional,
                                        self.options.long)
        self.linking = False
        doctests = []
        in_docstring = False
        unparsed_doc = False
        doc = []
        start = None
        tab_locations = []
        contains_line_number = False
        for lineno, line in self:
            if doctest_line_number.search(line) is not None:
                contains_line_number = True
            if "\t" in line:
                tab_locations.append(str(lineno+1))
            if "SAGE_DOCTEST_ALLOW_TABS" in line:
                tab_okay = True
            just_finished = False
            if in_docstring:
                if self.ending_docstring(line):
                    in_docstring = False
                    just_finished = True
                    self._process_doc(doctests, doc, namespace, start)
                    unparsed_doc = False
                else:
                    bitness = bitness_marker.search(line)
                    if bitness:
                        if bitness.groups()[0] != bitness_value:
                            self.line_shift += 1
                            continue
                        else:
                            line = line[:bitness.start()] + "\n"
                    if self.line_shift and sagestart.match(line):
                        # We insert blank lines to make up for the removed lines
                        doc.extend(["\n"]*self.line_shift)
                        self.line_shift = 0
                    doc.append(line)
                    unparsed_doc = True
            if not in_docstring and (not just_finished or self.start_finish_can_overlap):
                # to get line numbers in linked docstrings correct we
                # append a blank line to the doc list.
                doc.append("\n")
                if not line.strip():
                    continue
                if self.starting_docstring(line):
                    in_docstring = True
                    if self.linking:
                        # If there's already a doctest, we overwrite it.
                        if len(doctests) > 0:
                            doctests.pop()
                        if start is None:
                            start = lineno
                            doc = []
                    else:
                        self.line_shift = 0
                        start = lineno
                        doc = []
        # In ReST files we can end the file without decreasing the indentation level.
        if unparsed_doc:
            self._process_doc(doctests, doc, namespace, start)

        extras = dict(tab=not tab_okay and tab_locations,
                      line_number=contains_line_number,
                      optionals=self.parser.optionals)
        if self.options.randorder is not None and self.options.randorder is not False:
            # we want to randomize even when self.randorder = 0
            random.seed(self.options.randorder)
            randomized = []
            while doctests:
                i = random.randint(0, len(doctests) - 1)
                randomized.append(doctests.pop(i))
            return randomized, extras
        else:
            return doctests, extras


class StringDocTestSource(DocTestSource):
    r"""
    This class creates doctests from a string.

    INPUT:

    - ``basename`` -- string such as 'sage.doctests.sources', going
      into the names of created doctests and examples.

    - ``source`` -- a string, giving the source code to be parsed for
      doctests.

    - ``options`` -- a :class:`sage.doctest.control.DocTestDefaults`
      or equivalent.

    - ``printpath`` -- a string, to be used in place of a filename
      when doctest failures are displayed.

    - ``lineno_shift`` -- an integer (default: 0) by which to shift
      the line numbers of all doctests defined in this string.

    EXAMPLES::

        sage: from sage.doctest.control import DocTestDefaults
        sage: from sage.doctest.sources import StringDocTestSource, PythonSource
        sage: from sage.structure.dynamic_class import dynamic_class
        sage: s = "'''\n    sage: 2 + 2\n    4\n'''"
        sage: PythonStringSource = dynamic_class('PythonStringSource',(StringDocTestSource, PythonSource))
        sage: PSS = PythonStringSource('<runtime>', s, DocTestDefaults(), 'runtime')
        sage: dt, extras = PSS.create_doctests({})
        sage: len(dt)
        1
        sage: extras['tab']
        []
        sage: extras['line_number']
        False

        sage: s = "'''\n\tsage: 2 + 2\n\t4\n'''"
        sage: PSS = PythonStringSource('<runtime>', s, DocTestDefaults(), 'runtime')
        sage: dt, extras = PSS.create_doctests({})
        sage: extras['tab']
        ['2', '3']

        sage: s = "'''\n    sage: import warnings; warnings.warn('foo')\n    doctest:1: UserWarning: foo \n'''"
        sage: PSS = PythonStringSource('<runtime>', s, DocTestDefaults(), 'runtime')
        sage: dt, extras = PSS.create_doctests({})
        sage: extras['line_number']
        True
    """
    def __init__(self, basename, source, options, printpath, lineno_shift=0):
        r"""
        Initialization

        TESTS::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import StringDocTestSource, PythonSource
            sage: from sage.structure.dynamic_class import dynamic_class
            sage: s = "'''\n    sage: 2 + 2\n    4\n'''"
            sage: PythonStringSource = dynamic_class('PythonStringSource',(StringDocTestSource, PythonSource))
            sage: PSS = PythonStringSource('<runtime>', s, DocTestDefaults(), 'runtime')
            sage: TestSuite(PSS).run()
        """
        self.qualified_name = NestedName(basename)
        self.printpath = printpath
        self.source = source
        self.lineno_shift = lineno_shift
        DocTestSource.__init__(self, options)

    def __iter__(self):
        """
        Iterating over this source yields pairs ``(lineno, line)``.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import StringDocTestSource, PythonSource
            sage: from sage.structure.dynamic_class import dynamic_class
            sage: s = "'''\n    sage: 2 + 2\n    4\n'''"
            sage: PythonStringSource = dynamic_class('PythonStringSource',(StringDocTestSource, PythonSource))
            sage: PSS = PythonStringSource('<runtime>', s, DocTestDefaults(), 'runtime')
            sage: for n, line in PSS:
            ....:     print("{} {}".format(n, line))
            0 '''
            1     sage: 2 + 2
            2     4
            3 '''
        """
        for lineno, line in enumerate(self.source.split('\n')):
            yield lineno + self.lineno_shift, line + '\n'

    def create_doctests(self, namespace):
        r"""
        Creates doctests from this string.

        INPUT:

        - ``namespace`` -- a dictionary or :class:`sage.doctest.util.RecordingDict`.

        OUTPUT:

        - ``doctests`` -- a list of doctests defined by this string

        - ``tab_locations`` -- either False or a list of linenumbers
          on which tabs appear.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import StringDocTestSource, PythonSource
            sage: from sage.structure.dynamic_class import dynamic_class
            sage: s = "'''\n    sage: 2 + 2\n    4\n'''"
            sage: PythonStringSource = dynamic_class('PythonStringSource',(StringDocTestSource, PythonSource))
            sage: PSS = PythonStringSource('<runtime>', s, DocTestDefaults(), 'runtime')
            sage: dt, tabs = PSS.create_doctests({})
            sage: for t in dt:
            ....:     print("{} {}".format(t.name, t.examples[0].sage_source))
            <runtime> 2 + 2
        """
        return self._create_doctests(namespace)


class FileDocTestSource(DocTestSource):
    """
    This class creates doctests from a file.

    INPUT:

    - ``path`` -- string, the filename

    - ``options`` -- a :class:`sage.doctest.control.DocTestDefaults`
      instance or equivalent.

    EXAMPLES::

        sage: from sage.doctest.control import DocTestDefaults
        sage: from sage.doctest.sources import FileDocTestSource
        sage: from sage.env import SAGE_SRC
        sage: import os
        sage: filename = os.path.join(SAGE_SRC,'sage','doctest','sources.py')
        sage: FDS = FileDocTestSource(filename,DocTestDefaults())
        sage: FDS.basename
        'sage.doctest.sources'

    TESTS::

        sage: TestSuite(FDS).run()

    ::

        sage: from sage.doctest.control import DocTestDefaults
        sage: from sage.doctest.sources import FileDocTestSource
        sage: filename = tmp_filename(ext=".txtt")
        sage: FDS = FileDocTestSource(filename,DocTestDefaults())
        Traceback (most recent call last):
        ...
        ValueError: unknown extension for the file to test (=...txtt),
        valid extensions are: .py, .pyx, .pxd, .pxi, .sage, .spyx, .tex, .rst

    """
    def __init__(self, path, options):
        """
        Initialization.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','sources.py')
            sage: FDS = FileDocTestSource(filename,DocTestDefaults(randorder=0))
            sage: FDS.options.randorder
            0
        """
        self.path = path
        DocTestSource.__init__(self, options)
        base, ext = os.path.splitext(path)
        valid_code_ext = ('.py', '.pyx', '.pxd', '.pxi', '.sage', '.spyx')
        if ext in valid_code_ext:
            self.__class__ = dynamic_class('PythonFileSource',(FileDocTestSource,PythonSource))
            self.encoding = "utf-8"
        elif ext == '.tex':
            self.__class__ = dynamic_class('TexFileSource',(FileDocTestSource,TexSource))
            self.encoding = "utf-8"
        elif ext == '.rst':
            self.__class__ = dynamic_class('RestFileSource',(FileDocTestSource,RestSource))
            self.encoding = "utf-8"
        else:
            valid_ext = ", ".join(valid_code_ext + ('.tex', '.rst'))
            raise ValueError("unknown extension for the file to test (={}),"
                    " valid extensions are: {}".format(path, valid_ext))

    def __iter__(self):
        r"""
        Iterating over this source yields pairs ``(lineno, line)``.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: filename = tmp_filename(ext=".py")
            sage: s = "'''\n    sage: 2 + 2\n    4\n'''"
            sage: with open(filename, 'w') as f:
            ....:     _ = f.write(s)
            sage: FDS = FileDocTestSource(filename, DocTestDefaults())
            sage: for n, line in FDS:
            ....:     print("{} {}".format(n, line))
            0 '''
            1     sage: 2 + 2
            2     4
            3 '''

        The encoding is "utf-8" by default::

            sage: FDS.encoding
            'utf-8'

        We create a file with a Latin-1 encoding without declaring it::

            sage: s = b"'''\nRegardons le polyn\xF4me...\n'''\n"
            sage: with open(filename, 'wb') as f:
            ....:     _ = f.write(s)
            sage: FDS = FileDocTestSource(filename, DocTestDefaults())
            sage: L = list(FDS)
            Traceback (most recent call last):
            ...
            UnicodeDecodeError: 'utf...8' codec can...t decode byte 0xf4 in position 18: invalid continuation byte

        This works if we add a PEP 0263 encoding declaration::

            sage: s = b"#!/usr/bin/env python\n# -*- coding: latin-1 -*-\n" + s
            sage: with open(filename, 'wb') as f:
            ....:     _ = f.write(s)
            sage: FDS = FileDocTestSource(filename, DocTestDefaults())
            sage: L = list(FDS)
            sage: FDS.encoding
            'latin-1'
        """
        with open(self.path, 'rb') as source:
            for lineno, line in enumerate(source):
                if lineno < 2:
                    match = pep_0263.search(line)
                    if match:
                        self.encoding = bytes_to_str(match.group(1), 'ascii')
                yield lineno, line.decode(self.encoding)

    @lazy_attribute
    def printpath(self):
        """
        Whether the path is printed absolutely or relatively depends on an option.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: root = os.path.realpath(os.path.join(SAGE_SRC,'sage'))
            sage: filename = os.path.join(root,'doctest','sources.py')
            sage: cwd = os.getcwd()
            sage: os.chdir(root)
            sage: FDS = FileDocTestSource(filename,DocTestDefaults(randorder=0,abspath=False))
            sage: FDS.printpath
            'doctest/sources.py'
            sage: FDS = FileDocTestSource(filename,DocTestDefaults(randorder=0,abspath=True))
            sage: FDS.printpath
            '.../sage/doctest/sources.py'
            sage: os.chdir(cwd)
        """
        if self.options.abspath:
            return os.path.abspath(self.path)
        else:
            relpath = os.path.relpath(self.path)
            if relpath.startswith(".." + os.path.sep):
                return self.path
            else:
                return relpath

    @lazy_attribute
    def basename(self):
        """
        The basename of this file source, e.g. sage.doctest.sources

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','rings','integer.pyx')
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: FDS.basename
            'sage.rings.integer'
        """
        return get_basename(self.path)

    @lazy_attribute
    def in_lib(self):
        """
        Whether this file is part of a package (i.e. is in a directory
        containing an ``__init__.py`` file).

        Such files aren't loaded before running tests.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC, 'sage', 'rings', 'integer.pyx')
            sage: FDS = FileDocTestSource(filename, DocTestDefaults())
            sage: FDS.in_lib
            True
            sage: filename = os.path.join(SAGE_SRC, 'sage', 'doctest', 'tests', 'abort.rst')
            sage: FDS = FileDocTestSource(filename, DocTestDefaults())
            sage: FDS.in_lib
            False

        You can override the default::

            sage: FDS = FileDocTestSource("hello_world.py",DocTestDefaults())
            sage: FDS.in_lib
            False
            sage: FDS = FileDocTestSource("hello_world.py",DocTestDefaults(force_lib=True))
            sage: FDS.in_lib
            True
        """
        # We need an explicit bool() because is_package_dir() returns
        # 1/None instead of True/False.
        return bool(self.options.force_lib or
                is_package_dir(os.path.dirname(self.path)))

    def create_doctests(self, namespace):
        r"""
        Return a list of doctests for this file.

        INPUT:

        - ``namespace`` -- a dictionary or :class:`sage.doctest.util.RecordingDict`.

        OUTPUT:

        - ``doctests`` -- a list of doctests defined in this file.

        - ``extras`` -- a dictionary

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','sources.py')
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: len(doctests)
            41
            sage: extras['tab']
            False

        We give a self referential example::

            sage: doctests[18].name
            'sage.doctest.sources.FileDocTestSource.create_doctests'
            sage: doctests[18].examples[10].source
            'doctests[Integer(18)].examples[Integer(10)].source\n'

        TESTS:

        We check that we correctly process results that depend on 32
        vs 64 bit architecture::

            sage: import sys
            sage: bitness = '64' if sys.maxsize > (1 << 32) else '32'
            sage: gp.get_precision() == 38
            False # 32-bit
            True  # 64-bit
            sage: ex = doctests[18].examples[13]
            sage: (bitness == '64' and ex.want == 'True  \n') or (bitness == '32' and ex.want == 'False \n')
            True

        We check that lines starting with a # aren't doctested::

            #sage: raise RuntimeError
        """
        if not os.path.exists(self.path):
            import errno
            raise IOError(errno.ENOENT, "File does not exist", self.path)
        base, filename = os.path.split(self.path)
        _, ext = os.path.splitext(filename)
        if not self.in_lib and ext in ('.py', '.pyx', '.sage', '.spyx'):
            cwd = os.getcwd()
            if base:
                os.chdir(base)
            try:
                load(filename, namespace) # errors raised here will be caught in DocTestTask
            finally:
                os.chdir(cwd)
        self.qualified_name = NestedName(self.basename)
        return self._create_doctests(namespace)

    def _test_enough_doctests(self, check_extras=True, verbose=True):
        """
        This function checks to see that the doctests are not getting
        unexpectedly skipped.  It uses a different (and simpler) code
        path than the doctest creation functions, so there are a few
        files in Sage that it counts incorrectly.

        INPUT:

        - ``check_extras`` -- bool (default True), whether to check if
          doctests are created that don't correspond to either a
          ``sage: `` or a ``>>> `` prompt.

        - ``verbose`` -- bool (default True), whether to print
          offending line numbers when there are missing or extra
          tests.

        TESTS::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.env import SAGE_SRC
            sage: cwd = os.getcwd()
            sage: os.chdir(SAGE_SRC)
            sage: import itertools
            sage: for path, dirs, files in itertools.chain(os.walk('sage'), os.walk('doc')): # long time
            ....:     path = os.path.relpath(path)
            ....:     dirs.sort(); files.sort()
            ....:     for F in files:
            ....:         _, ext = os.path.splitext(F)
            ....:         if ext in ('.py', '.pyx', '.pxd', '.pxi', '.sage', '.spyx', '.rst'):
            ....:             filename = os.path.join(path, F)
            ....:             FDS = FileDocTestSource(filename, DocTestDefaults(long=True, optional=True, force_lib=True))
            ....:             FDS._test_enough_doctests(verbose=False)
            There are 3 unexpected tests being run in sage/doctest/parsing.py
            There are 1 unexpected tests being run in sage/doctest/reporting.py
            sage: os.chdir(cwd)
        """
        expected = []
        rest = isinstance(self, RestSource)
        if rest:
            skipping = False
            in_block = False
            last_line = ''
        for lineno, line in self:
            if not line.strip():
                continue
            if rest:
                if line.strip().startswith(".. nodoctest"):
                    return
                # We need to track blocks in order to figure out whether we're skipping.
                if in_block:
                    indent = whitespace.match(line).end()
                    if indent <= starting_indent:
                        in_block = False
                        skipping = False
                if not in_block:
                    m1 = double_colon.match(line)
                    m2 = code_block.match(line.lower())
                    starting = (m1 and not line.strip().startswith(".. ")) or m2
                    if starting:
                        if ".. skip" in last_line:
                            skipping = True
                        in_block = True
                        starting_indent = whitespace.match(line).end()
                last_line = line
            if (not rest or in_block) and sagestart.match(line) and not ((rest and skipping) or untested.search(line.lower())):
                expected.append(lineno+1)
        actual = []
        tests, _ = self.create_doctests({})
        for dt in tests:
            if dt.examples:
                for ex in dt.examples[:-1]: # the last entry is a sig_on_count()
                    actual.append(dt.lineno + ex.lineno + 1)
        shortfall = sorted(set(expected).difference(set(actual)))
        extras = sorted(set(actual).difference(set(expected)))
        if len(actual) == len(expected):
            if not shortfall:
                return
            dif = extras[0] - shortfall[0]
            for e, s in zip(extras[1:],shortfall[1:]):
                if dif != e - s:
                    break
            else:
                print("There are %s tests in %s that are shifted by %s" % (len(shortfall), self.path, dif))
                if verbose:
                    print("    The correct line numbers are %s" % (", ".join(str(n) for n in shortfall)))
                return
        elif len(actual) < len(expected):
            print("There are %s tests in %s that are not being run" % (len(expected) - len(actual), self.path))
        elif check_extras:
            print("There are %s unexpected tests being run in %s" % (len(actual) - len(expected), self.path))
        if verbose:
            if shortfall:
                print("    Tests on lines %s are not run" % (", ".join(str(n) for n in shortfall)))
            if check_extras and extras:
                print("    Tests on lines %s seem extraneous" % (", ".join(str(n) for n in extras)))


class SourceLanguage:
    """
    An abstract class for functions that depend on the programming language of a doctest source.

    Currently supported languages include Python, ReST and LaTeX.
    """
    def parse_docstring(self, docstring, namespace, start):
        """
        Return a list of doctest defined in this docstring.

        This function is called by :meth:`DocTestSource._process_doc`.
        The default implementation, defined here, is to use the
        :class:`sage.doctest.parsing.SageDocTestParser` attached to
        this source to get doctests from the docstring.

        INPUT:

        - ``docstring`` -- a string containing documentation and tests.

        - ``namespace`` -- a dictionary or :class:`sage.doctest.util.RecordingDict`.

        - ``start`` -- an integer, one less than the starting line number

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.parsing import SageDocTestParser
            sage: from sage.doctest.util import NestedName
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','util.py')
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: doctests, _ = FDS.create_doctests({})
            sage: for dt in doctests:
            ....:     FDS.qualified_name = dt.name
            ....:     dt.examples = dt.examples[:-1] # strip off the sig_on() test
            ....:     assert(FDS.parse_docstring(dt.docstring,{},dt.lineno-1)[0] == dt)
        """
        return [self.parser.get_doctest(docstring, namespace, str(self.qualified_name),
                                        self.printpath, start + 1)]

class PythonSource(SourceLanguage):
    """
    This class defines the functions needed for the extraction of doctests from python sources.

    EXAMPLES::

        sage: from sage.doctest.control import DocTestDefaults
        sage: from sage.doctest.sources import FileDocTestSource
        sage: from sage.env import SAGE_SRC
        sage: import os
        sage: filename = os.path.join(SAGE_SRC,'sage','doctest','sources.py')
        sage: FDS = FileDocTestSource(filename,DocTestDefaults())
        sage: type(FDS)
        <class 'sage.doctest.sources.PythonFileSource'>
    """
    # The same line can't both start and end a docstring
    start_finish_can_overlap = False

    def _init(self):
        """
        This function is called before creating doctests from a Python source.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','sources.py')
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: FDS._init()
            sage: FDS.last_indent
            -1
        """
        self.last_indent = -1
        self.last_line = None
        self.quotetype = None
        self.paren_count = 0
        self.bracket_count = 0
        self.curly_count = 0
        self.code_wrapping = False

    def _update_quotetype(self, line):
        r"""
        Updates the track of what kind of quoted string we're in.

        We need to track whether we're inside a triple quoted
        string, since a triple quoted string that starts a line
        could be the end of a string and thus not the beginning of a
        doctest (see sage.misc.sageinspect for an example).

        To do this tracking we need to track whether we're inside a
        string at all, since ''' inside a string doesn't start a
        triple quote (see the top of this file for an example).

        We also need to track parentheses and brackets, since we only
        want to update our record of last line and indentation level
        when the line is actually over.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','sources.py')
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: FDS._init()
            sage: FDS._update_quotetype('\"\"\"'); print(" ".join(list(FDS.quotetype)))
            " " "
            sage: FDS._update_quotetype("'''"); print(" ".join(list(FDS.quotetype)))
            " " "
            sage: FDS._update_quotetype('\"\"\"'); print(FDS.quotetype)
            None
            sage: FDS._update_quotetype("triple_quotes = re.compile(\"\\s*[rRuU]*((''')|(\\\"\\\"\\\"))\")")
            sage: print(FDS.quotetype)
            None
            sage: FDS._update_quotetype("''' Single line triple quoted string \\''''")
            sage: print(FDS.quotetype)
            None
            sage: FDS._update_quotetype("' Lots of \\\\\\\\'")
            sage: print(FDS.quotetype)
            None
        """
        def _update_parens(start,end=None):
            self.paren_count += line.count("(",start,end) - line.count(")",start,end)
            self.bracket_count += line.count("[",start,end) - line.count("]",start,end)
            self.curly_count += line.count("{",start,end) - line.count("}",start,end)
        pos = 0
        while pos < len(line):
            if self.quotetype is None:
                next_single = line.find("'",pos)
                next_double = line.find('"',pos)
                if next_single == -1 and next_double == -1:
                    next_comment = line.find("#",pos)
                    if next_comment == -1:
                        _update_parens(pos)
                    else:
                        _update_parens(pos,next_comment)
                    break
                elif next_single == -1:
                    m = next_double
                elif next_double == -1:
                    m = next_single
                else:
                    m = min(next_single, next_double)
                next_comment = line.find('#',pos,m)
                if next_comment != -1:
                    _update_parens(pos,next_comment)
                    break
                _update_parens(pos,m)
                if m+2 < len(line) and line[m] == line[m+1] == line[m+2]:
                    self.quotetype = line[m:m+3]
                    pos = m+3
                else:
                    self.quotetype = line[m]
                    pos = m+1
            else:
                next = line.find(self.quotetype,pos)
                if next == -1:
                    break
                elif next == 0 or line[next-1] != '\\':
                    pos = next + len(self.quotetype)
                    self.quotetype = None
                else:
                    # We need to worry about the possibility that
                    # there are an even number of backslashes before
                    # the quote, in which case it is not escaped
                    count = 1
                    slashpos = next - 2
                    while slashpos >= pos and line[slashpos] == '\\':
                        count += 1
                        slashpos -= 1
                    if count % 2 == 0:
                        pos = next + len(self.quotetype)
                        self.quotetype = None
                    else:
                        # The possible ending quote was escaped.
                        pos = next + 1

    def starting_docstring(self, line):
        """
        Determines whether the input line starts a docstring.

        If the input line does start a docstring (a triple quote),
        then this function updates ``self.qualified_name``.

        INPUT:

        - ``line`` -- a string, one line of an input file

        OUTPUT:

        - either None or a Match object.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.util import NestedName
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','sources.py')
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: FDS._init()
            sage: FDS.starting_docstring("r'''")
            <...Match object...>
            sage: FDS.ending_docstring("'''")
            <...Match object...>
            sage: FDS.qualified_name = NestedName(FDS.basename)
            sage: FDS.starting_docstring("class MyClass(object):")
            sage: FDS.starting_docstring("    def hello_world(self):")
            sage: FDS.starting_docstring("        '''")
            <...Match object...>
            sage: FDS.qualified_name
            sage.doctest.sources.MyClass.hello_world
            sage: FDS.ending_docstring("    '''")
            <...Match object...>
            sage: FDS.starting_docstring("class NewClass(object):")
            sage: FDS.starting_docstring("    '''")
            <...Match object...>
            sage: FDS.ending_docstring("    '''")
            <...Match object...>
            sage: FDS.qualified_name
            sage.doctest.sources.NewClass
            sage: FDS.starting_docstring("print(")
            sage: FDS.starting_docstring("    '''Not a docstring")
            sage: FDS.starting_docstring("    ''')")
            sage: FDS.starting_docstring("def foo():")
            sage: FDS.starting_docstring("    '''This is a docstring'''")
            <...Match object...>
        """
        indent = whitespace.match(line).end()
        quotematch = None
        if self.quotetype is None and not self.code_wrapping:
            # We're not inside a triple quote and not inside code like
            # print(
            #     """Not a docstring
            #     """)

            if line[indent] != '#' and (indent == 0 or indent > self.last_indent):
                quotematch = triple_quotes.match(line)
                # It would be nice to only run the name_regex when
                # quotematch wasn't None, but then we mishandle classes
                # that don't have a docstring.
                if not self.code_wrapping and self.last_indent >= 0 and indent > self.last_indent:
                    name = name_regex.match(self.last_line)
                    if name:
                        name = name.groups()[0]
                        self.qualified_name[indent] = name
                    elif quotematch:
                        self.qualified_name[indent] = '?'
        self._update_quotetype(line)
        if line[indent] != '#' and not self.code_wrapping:
            self.last_line, self.last_indent = line, indent
        self.code_wrapping = not (self.paren_count == self.bracket_count == self.curly_count == 0)
        return quotematch

    def ending_docstring(self, line):
        r"""
        Determines whether the input line ends a docstring.

        INPUT:

        - ``line`` -- a string, one line of an input file.

        OUTPUT:

        - an object that, when evaluated in a boolean context, gives
          True or False depending on whether the input line marks the
          end of a docstring.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.util import NestedName
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','sources.py')
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: FDS._init()
            sage: FDS.quotetype = "'''"
            sage: FDS.ending_docstring("'''")
            <...Match object...>
            sage: FDS.ending_docstring('\"\"\"')
        """
        quotematch = triple_quotes.match(line)
        if quotematch is not None and quotematch.groups()[0] != self.quotetype:
            quotematch = None
        self._update_quotetype(line)
        return quotematch

    def _neutralize_doctests(self, reindent):
        r"""
        Return a string containing the source of ``self``, but with
        doctests modified so they are not tested.

        This function is used in creating doctests for ReST files,
        since docstrings of Python functions defined inside verbatim
        blocks screw up Python's doctest parsing.

        INPUT:

        - ``reindent`` -- an integer, the number of spaces to indent
          the result.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import StringDocTestSource, PythonSource
            sage: from sage.structure.dynamic_class import dynamic_class
            sage: s = "'''\n    sage: 2 + 2\n    4\n'''"
            sage: PythonStringSource = dynamic_class('PythonStringSource',(StringDocTestSource, PythonSource))
            sage: PSS = PythonStringSource('<runtime>', s, DocTestDefaults(), 'runtime')
            sage: print(PSS._neutralize_doctests(0))
            '''
                safe: 2 + 2
                4
            '''
        """
        neutralized = []
        in_docstring = False
        self._init()
        for lineno, line in self:
            if not line.strip():
                neutralized.append(line)
            elif in_docstring:
                if self.ending_docstring(line):
                    in_docstring = False
                neutralized.append(" "*reindent + find_prompt.sub(r"\1safe:\3",line))
            else:
                if self.starting_docstring(line):
                    in_docstring = True
                neutralized.append(" "*reindent + line)
        return "".join(neutralized)

class TexSource(SourceLanguage):
    """
    This class defines the functions needed for the extraction of
    doctests from a LaTeX source.

    EXAMPLES::

        sage: from sage.doctest.control import DocTestDefaults
        sage: from sage.doctest.sources import FileDocTestSource
        sage: filename = "sage_paper.tex"
        sage: FDS = FileDocTestSource(filename,DocTestDefaults())
        sage: type(FDS)
        <class 'sage.doctest.sources.TexFileSource'>
    """
    # The same line can't both start and end a docstring
    start_finish_can_overlap = False

    def _init(self):
        """
        This function is called before creating doctests from a Tex file.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: filename = "sage_paper.tex"
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: FDS._init()
            sage: FDS.skipping
            False
        """
        self.skipping = False

    def starting_docstring(self, line):
        r"""
        Determines whether the input line starts a docstring.

        Docstring blocks in tex files are defined by verbatim or
        lstlisting environments, and can be linked together by adding
        %link immediately after the \end{verbatim} or \end{lstlisting}.

        Within a verbatim (or lstlisting) block, you can tell Sage not to
        process the rest of the block by including a %skip line.

        INPUT:

        - ``line`` -- a string, one line of an input file

        OUTPUT:

        - a boolean giving whether the input line marks the
          start of a docstring (verbatim block).

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: filename = "sage_paper.tex"
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: FDS._init()

        We start docstrings with \begin{verbatim} or \begin{lstlisting}::

            sage: FDS.starting_docstring(r"\begin{verbatim}")
            True
            sage: FDS.starting_docstring(r"\begin{lstlisting}")
            True
            sage: FDS.skipping
            False
            sage: FDS.ending_docstring("sage: 2+2")
            False
            sage: FDS.ending_docstring("4")
            False

        To start ignoring the rest of the verbatim block, use %skip::

            sage: FDS.ending_docstring("%skip")
            True
            sage: FDS.skipping
            True
            sage: FDS.starting_docstring("sage: raise RuntimeError")
            False

        You can even pretend to start another verbatim block while skipping::

            sage: FDS.starting_docstring(r"\begin{verbatim}")
            False
            sage: FDS.skipping
            True

        To stop skipping end the verbatim block::

            sage: FDS.starting_docstring(r"\end{verbatim} %link")
            False
            sage: FDS.skipping
            False

        Linking works even when the block was ended while skipping::

            sage: FDS.linking
            True
            sage: FDS.starting_docstring(r"\begin{verbatim}")
            True
        """
        if self.skipping:
            if self.ending_docstring(line, check_skip=False):
                self.skipping = False
            return False
        return bool(begin_verb.match(line) or begin_lstli.match(line))

    def ending_docstring(self, line, check_skip=True):
        r"""
        Determines whether the input line ends a docstring.

        Docstring blocks in tex files are defined by verbatim or
        lstlisting environments, and can be linked together by adding
        %link immediately after the \end{verbatim} or \end{lstlisting}.

        Within a verbatim (or lstlisting) block, you can tell Sage not to
        process the rest of the block by including a %skip line.

        INPUT:

        - ``line`` -- a string, one line of an input file

        - ``check_skip`` -- boolean (default True), used internally in starting_docstring.

        OUTPUT:

        - a boolean giving whether the input line marks the
          end of a docstring (verbatim block).

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: filename = "sage_paper.tex"
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: FDS._init()
            sage: FDS.ending_docstring(r"\end{verbatim}")
            True
            sage: FDS.ending_docstring(r"\end{lstlisting}")
            True
            sage: FDS.linking
            False

        Use %link to link with the next verbatim block::

            sage: FDS.ending_docstring(r"\end{verbatim}%link")
            True
            sage: FDS.linking
            True

        %skip also ends a docstring block::

            sage: FDS.ending_docstring("%skip")
            True
        """
        m = end_verb.match(line)
        if m:
            if m.groups()[0]:
                self.linking = True
            else:
                self.linking = False
            return True
        m = end_lstli.match(line)
        if m:
            if m.groups()[0]:
                self.linking = True
            else:
                self.linking = False
            return True
        if check_skip and skip.match(line):
            self.skipping = True
            return True
        return False


class RestSource(SourceLanguage):
    """
    This class defines the functions needed for the extraction of
    doctests from ReST sources.

    EXAMPLES::

        sage: from sage.doctest.control import DocTestDefaults
        sage: from sage.doctest.sources import FileDocTestSource
        sage: filename = "sage_doc.rst"
        sage: FDS = FileDocTestSource(filename,DocTestDefaults())
        sage: type(FDS)
        <class 'sage.doctest.sources.RestFileSource'>
    """
    # The same line can both start and end a docstring
    start_finish_can_overlap = True

    def _init(self):
        """
        This function is called before creating doctests from a ReST file.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: filename = "sage_doc.rst"
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: FDS._init()
            sage: FDS.link_all
            False
        """
        self.link_all = False
        self.last_line = ""
        self.last_indent = -1
        self.first_line = False
        self.skipping = False

    def starting_docstring(self, line):
        """
        A line ending with a double colon starts a verbatim block in a ReST file,
        as does a line containing ``.. CODE-BLOCK:: language``.

        This function also determines whether the docstring block
        should be joined with the previous one, or should be skipped.

        INPUT:

        - ``line`` -- a string, one line of an input file

        OUTPUT:

        - either None or a Match object.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: filename = "sage_doc.rst"
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: FDS._init()
            sage: FDS.starting_docstring("Hello world::")
            True
            sage: FDS.ending_docstring("    sage: 2 + 2")
            False
            sage: FDS.ending_docstring("    4")
            False
            sage: FDS.ending_docstring("We are now done")
            True
            sage: FDS.starting_docstring(".. link")
            sage: FDS.starting_docstring("::")
            True
            sage: FDS.linking
            True
        """
        if link_all.match(line):
            self.link_all = True
        if self.skipping:
            end_block = self.ending_docstring(line)
            if end_block:
                self.skipping = False
            else:
                return False
        m1 = double_colon.match(line)
        m2 = code_block.match(line.lower())
        starting = (m1 and not line.strip().startswith(".. ")) or m2
        if starting:
            self.linking = self.link_all or '.. link' in self.last_line
            self.first_line = True
            m = m1 or m2
            indent = len(m.groups()[0])
            if '.. skip' in self.last_line:
                self.skipping = True
                starting = False
        else:
            indent = self.last_indent
        self.last_line, self.last_indent = line, indent
        return starting

    def ending_docstring(self, line):
        """
        When the indentation level drops below the initial level the
        block ends.

        INPUT:

        - ``line`` -- a string, one line of an input file

        OUTPUT:

        - a boolean, whether the verbatim block is ending.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: filename = "sage_doc.rst"
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: FDS._init()
            sage: FDS.starting_docstring("Hello world::")
            True
            sage: FDS.ending_docstring("    sage: 2 + 2")
            False
            sage: FDS.ending_docstring("    4")
            False
            sage: FDS.ending_docstring("We are now done")
            True
        """
        if not line.strip():
            return False
        indent = whitespace.match(line).end()
        if self.first_line:
            self.first_line = False
            if indent <= self.last_indent:
                # We didn't indent at all
                return True
            self.last_indent = indent
        return indent < self.last_indent

    def parse_docstring(self, docstring, namespace, start):
        r"""
        Return a list of doctest defined in this docstring.

        Code blocks in a REST file can contain python functions with
        their own docstrings in addition to in-line doctests.  We want
        to include the tests from these inner docstrings, but Python's
        doctesting module has a problem if we just pass on the whole
        block, since it expects to get just a docstring, not the
        Python code as well.

        Our solution is to create a new doctest source from this code
        block and append the doctests created from that source.  We
        then replace the occurrences of "sage:" and ">>>" occurring
        inside a triple quote with "safe:" so that the doctest module
        doesn't treat them as tests.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.parsing import SageDocTestParser
            sage: from sage.doctest.util import NestedName
            sage: filename = "sage_doc.rst"
            sage: FDS = FileDocTestSource(filename,DocTestDefaults())
            sage: FDS.parser = SageDocTestParser(set(['sage']))
            sage: FDS.qualified_name = NestedName('sage_doc')
            sage: s = "Some text::\n\n    def example_python_function(a, \
            ....:      b):\n        '''\n        Brief description \
            ....:      of function.\n\n        EXAMPLES::\n\n            \
            ....:      sage: test1()\n            sage: test2()\n        \
            ....:      '''\n        return a + b\n\n    sage: test3()\n\nMore \
            ....:      ReST documentation."
            sage: tests = FDS.parse_docstring(s, {}, 100)
            sage: len(tests)
            2
            sage: for ex in tests[0].examples:
            ....:     print(ex.sage_source)
            test3()
            sage: for ex in tests[1].examples:
            ....:     print(ex.sage_source)
            test1()
            test2()
            sig_on_count() # check sig_on/off pairings (virtual doctest)
        """
        PythonStringSource = dynamic_class("sage.doctest.sources.PythonStringSource",
                                           (StringDocTestSource, PythonSource))
        min_indent = self.parser._min_indent(docstring)
        pysource = '\n'.join([l[min_indent:] for l in docstring.split('\n')])
        inner_source = PythonStringSource(self.basename, pysource,
                                          self.options,
                                          self.printpath, lineno_shift=start+1)
        inner_doctests, _ = inner_source._create_doctests(namespace, True)
        safe_docstring = inner_source._neutralize_doctests(min_indent)
        outer_doctest = self.parser.get_doctest(safe_docstring, namespace,
                                                str(self.qualified_name),
                                                self.printpath, start + 1)
        return [outer_doctest] + inner_doctests

class DictAsObject(dict):
    """
    A simple subclass of dict that inserts the items from the initializing dictionary into attributes.

    EXAMPLES::

        sage: from sage.doctest.sources import DictAsObject
        sage: D = DictAsObject({'a':2})
        sage: D.a
        2
    """
    def __init__(self, attrs):
        """
        Initialization.

        INPUT:

        - ``attrs`` -- a dictionary.

        EXAMPLES::

            sage: from sage.doctest.sources import DictAsObject
            sage: D = DictAsObject({'a':2})
            sage: D.a == D['a']
            True
            sage: D.a
            2
        """
        super(DictAsObject, self).__init__(attrs)
        self.__dict__.update(attrs)

    def __setitem__(self, ky, val):
        """
        We preserve the ability to access entries through either the
        dictionary or attribute interfaces.

        EXAMPLES::

            sage: from sage.doctest.sources import DictAsObject
            sage: D = DictAsObject({})
            sage: D['a'] = 2
            sage: D.a
            2
        """
        super(DictAsObject, self).__setitem__(ky, val)
        try:
            super(DictAsObject, self).__setattr__(ky, val)
        except TypeError:
            pass

    def __setattr__(self, ky, val):
        """
        We preserve the ability to access entries through either the
        dictionary or attribute interfaces.

        EXAMPLES::

            sage: from sage.doctest.sources import DictAsObject
            sage: D = DictAsObject({})
            sage: D.a = 2
            sage: D['a']
            2
        """
        super(DictAsObject, self).__setitem__(ky, val)
        super(DictAsObject, self).__setattr__(ky, val)
