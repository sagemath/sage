# -*- coding: utf-8 -*-
r"""
Test that the Sage library uses Python 3 syntax

AUTHORS:

- Volker Braun (2017-02-11): initial version

- Frédéric Chapoton (2017-02-28): fixes

- Julian Rüth (2017-03-14): documentation changes

EXAMPLES::

    sage: from sage.tests.py3_syntax import Python3SyntaxTest
    sage: py3_syntax = Python3SyntaxTest('sage', 'sage_setup')
    sage: py3_syntax.run_tests('.py')   # long time
"""
# ****************************************************************************
#       Copyright (C) 2017 Volker Braun <vbraun.name@gmail.com>
#                     2017 Frédéric Chapoton <chapoton@math.univ-lyon1.fr>
#                     2017 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import print_function

import os
import itertools
import subprocess

from sage.env import SAGE_SRC
from sage.cpython.string import bytes_to_str


class SortedDirectoryWalkerABC(object):
    r"""
    Walks the directory tree in a reproducible manner.

    INPUT:

    - ``*directories`` -- roots of the directory trees to walk

    EXAMPLES::

        sage: from sage.tests.py3_syntax import SortedDirectoryWalkerABC
        sage: walker = SortedDirectoryWalkerABC('sage/tests')

    """
    def __init__(self, *directories):
        r"""
        TESTS::

            sage: from sage.tests.py3_syntax import Python3SyntaxTest, SortedDirectoryWalkerABC
            sage: isinstance(Python3SyntaxTest('sage/tests'), SortedDirectoryWalkerABC)
            True
        """
        self._directories = tuple(
            os.path.join(SAGE_SRC, d) for d in directories
        )

    def __iter__(self):
        r"""
        Return an iterator over the files in the directories.

        OUTPUT:

        Iterator over triples consisting of path, filename, and file
        extension. The iteration order only depends on the filenames
        and not on filesystem layout.

        EXAMPLES::

            sage: from sage.tests.py3_syntax import Python3SyntaxTest
            sage: name = 'sage/tests/books/computational-mathematics-with-sagemath'
            sage: test = Python3SyntaxTest(name)
            sage: next(iter(test))
            ('.../sage/tests/books/computational-mathematics-with-sagemath', 'README', '')
        """
        tree_walk = itertools.chain(*map(os.walk, self._directories))
        for path, _, files in tree_walk:
            path = os.path.relpath(path)
            files.sort()
            for filename in files:
                if filename.startswith('.'):
                    continue
                _, ext = os.path.splitext(filename)
                yield (path, filename, ext)

    def run_tests(self, *extensions):
        r"""
        Run :meth:`test` on each file with a matching extension.

        INPUT:

        - ``*extensions`` -- the file extensions (including the leading dot) to
          check

        EXAMPLES::

            sage: from sage.tests.py3_syntax import SortedDirectoryWalkerABC
            sage: walker = SortedDirectoryWalkerABC('sage/tests/french_book')
            sage: walker.run_tests('.py')
            Traceback (most recent call last):
            ...
            NotImplementedError: to be implemented in derived class
        """
        buf = []
        for (path, filename, ext) in self:
            if ext in extensions:
                fqn = os.path.join(path, filename)
                buf.append(fqn)
                if len(buf) > 20:
                    self.test(*buf)
                    buf = []
        self.test(*buf)

    def test(self, *filenames):
        r"""
        Test the given filenames.

        To be implemented in derived classes.

        INPUT:

        - ``*filenames`` -- the fully qualified filenames to check

        EXAMPLES::

            sage: from sage.tests.py3_syntax import SortedDirectoryWalkerABC
            sage: walker = SortedDirectoryWalkerABC()
            sage: from sage.env import SAGE_SRC
            sage: filename = os.path.join(SAGE_SRC, 'sage', 'tests', 'py3_syntax.py')
            sage: walker.test(filename)
            Traceback (most recent call last):
            ...
            NotImplementedError: to be implemented in derived class
        """
        raise NotImplementedError('to be implemented in derived class')


PYTHON3_ENV = dict(os.environ)
PYTHON3_ENV.pop('PYTHONPATH', None)

class Python3SyntaxTest(SortedDirectoryWalkerABC):
    r"""
    Walks the directory tree and tests that all Python files are valid Python 3
    syntax.

    INPUT:

    - ``*directories`` -- roots of the directory trees to walk

    EXAMPLES::

        sage: from sage.tests.py3_syntax import Python3SyntaxTest
        sage: walker = Python3SyntaxTest('sage/tests')
        sage: walker.run_tests('.py')   # long time

    """
    def test(self, *filenames):
        r"""
        Print an error message if the given files are not using valid Python 3
        syntax.

        INPUT:

        - ``*filenames`` -- the fully qualified filenames to check

        EXAMPLES::

            sage: import os, tempfile
            sage: src = tempfile.NamedTemporaryFile(suffix='.py', mode='w+', delete=False)
            sage: _ = src.write('print "invalid print statement"')
            sage: src.close()
            sage: from sage.tests.py3_syntax import Python3SyntaxTest
            sage: py3_syntax = Python3SyntaxTest()
            sage: py3_syntax.test(src.name)
            Invalid Python 3 syntax found:
            Missing parentheses in call to 'print'.
            Did you mean print("invalid print statement")? (...py, line 1)
            sage: os.unlink(src.name)
        """

        # compile all given files in memory, printing all errors
        # inspired by the py_compile module (but without writing to file)
        script = """
import sys
import importlib.machinery
rv = 0
for file in sys.argv[1:]:
    loader = importlib.machinery.SourceFileLoader('<sage_test>', file)
    source_bytes = loader.get_data(file)
    try:
        code = loader.source_to_code(source_bytes, file)
    except Exception as err:
        print(err)
        rv = 1
sys.exit(rv)
"""
        cmd = [
            'python3',
            '-c',
            script,
        ] + list(filenames)
        process = subprocess.Popen(
            cmd,
            env=PYTHON3_ENV,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        if process.returncode == 0:
            return
        print('Invalid Python 3 syntax found:')
        if stdout:
            print(bytes_to_str(stdout))
        if stderr:
            print(bytes_to_str(stderr))
