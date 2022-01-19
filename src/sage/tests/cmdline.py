r"""
This file contains some tests that Sage command line options actually
do something.

We test the following command line options:

test.py
/path/to/test.py
test.sage
/path/to/test.sage
test.spyx
/path/to/test.spyx
--advanced
-c
--cython
--dev
--ecl
--experimental
--fixdoctests
--gap
--gdb
--gp
-h
--help
--info
--ipython
--lisp
--maxima
--min
--mwrank
--optional
--preparse
--python
--python3
-q
--R
--root
--rst2ipynb
--ipynb2rst
--sh
--singular
--sqlite3
--standard
--startuptime
-t
-v
--zzfoobar (illegal option)


AUTHORS:

- Jeroen Demeyer (2010-11-20): initial version (:trac:`10300`)

"""
from subprocess import Popen, PIPE
import os
import sys
import select


def test_executable(args, input="", timeout=100.0, pydebug_ignore_warnings=False, **kwds):
    r"""
    Run the program defined by ``args`` using the string ``input`` on
    the standard input.

    INPUT:

    - ``args`` -- a list of program arguments, the first being the
      executable.

    - ``input`` -- a string serving as standard input.  Usually, this
      should end with a newline.

    - ``timeout`` -- if the program produces no output for ``timeout``
      seconds, a RuntimeError is raised.

    - ``pydebug_ignore_warnings`` -- boolean. Set the PYTHONWARNINGS environment variable to ignore
      Python warnings when on a Python debug build (`--with-pydebug`, e.g. from building with
      `SAGE_DEBUG=yes`). Debug builds do not install the default warning filters, which can break
      some doctests. Unfortunately the environment variable does not support regex message filters,
      so the filter will catch a bit more than the default filters. Hence we only enable it on debug
      builds.

    - ``**kwds`` -- Additional keyword arguments passed to the
      :class:`Popen` constructor.

    OUTPUT: a tuple ``(out, err, ret)`` with the standard output,
    standard error and exitcode of the program run.


    EXAMPLES::

        sage: from sage.tests.cmdline import test_executable
        sage: (out, err, ret) = test_executable(["cat"], "Hello World!")
        sage: out
        'Hello World!'
        sage: err
        ''
        sage: ret
        0

    We test the timeout option::

        sage: (out, err, ret) = test_executable(["sleep", "1"], timeout=0.1)
        Traceback (most recent call last):
        ...
        RuntimeError: timeout in test_executable()

    TESTS:

    Run Sage itself with various options::

        sage: (out, err, ret) = test_executable(["sage"], pydebug_ignore_warnings=True)
        sage: out.find(version()) >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage"], "3^33\n", pydebug_ignore_warnings=True)
        sage: out.find(version()) >= 0
        True
        sage: out.find("5559060566555523") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "-q"], "3^33\n", pydebug_ignore_warnings=True)
        sage: out.find(version()) >= 0
        False
        sage: out.find("5559060566555523") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "-c", "print(3^33)"])
        sage: print(out)
        5559060566555523
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--min", "-c", "print(3^33)"])
        sage: print(out)
        5559060566555523
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--startuptime"])
        sage: out.find("Slowest module import") >= 0
        True
        sage: err
        ''
        sage: ret
        0

    Test help::

        sage: (out, err, ret) = test_executable(["sage", "-h"])
        sage: out.find("Optional arguments:") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--help"])
        sage: out.find("Optional arguments:") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--advanced"])
        sage: out.find("run the Sage cleaner.") >= 0
        True
        sage: err
        ''
        sage: ret
        0
        sage: out.find("print the Sage root directory") >= 0 # optional - build
        True
        sage: out.find("regular expression search through the Sage") >= 0 # optional - build
        True

    Basic information about the Sage installation::

        sage: (out, err, ret) = test_executable(["sage", "-v"])
        sage: out.find(version()) >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--root"])  # optional - build
        sage: len(out) >= 2   # at least one character + newline; optional - build
        True
        sage: err  # optional - build
        ''
        sage: ret  # optional -build
        0

    Test ``sage --info [packages]`` and the equivalent
    ``sage -p --info --info [packages]`` (the doubling of ``--info``
    is intentional, that option should be idempotent)::

        sage: out, err, ret = test_executable(["sage", "--info", "sqlite"])  # optional - build
        sage: print(out)  # optional - build
        sqlite...
        SQLite is a software library that implements a self-contained,
        serverless, zero-configuration, transactional SQL database engine.
        ...
        sage: err  # optional - build
        ''
        sage: ret  # optional - build
        0

        sage: out, err, ret = test_executable(["sage", "-p", "--info", "--info", "sqlite"])  # optional - build
        sage: print(out)  # optional - build
        sqlite...
        SQLite is a software library that implements a self-contained,
        serverless, zero-configuration, transactional SQL database engine.
        ...
        sage: err  # optional - build
        ''
        sage: ret  # optional - build
        0

    Test ``sage-run`` on a Python file, both with an absolute and with a relative path::

        sage: dir = tmp_dir(); name = 'python_test_file.py'
        sage: fullname = os.path.join(dir, name)
        sage: F = open(fullname, 'w')
        sage: _ = F.write("print(3^33)\n")
        sage: F.close()
        sage: (out, err, ret) = test_executable(["sage", fullname])
        sage: print(out)
        34
        sage: err
        ''
        sage: ret
        0
        sage: (out, err, ret) = test_executable(["sage", name], cwd=dir)
        sage: print(out)
        34
        sage: err
        ''
        sage: ret
        0

    The same as above, but now with a ``.sage`` file.  This indirectly
    also tests the preparser::

        sage: dir = tmp_dir(); name = 'sage_test_file.sage'
        sage: fullname = os.path.join(dir, name)
        sage: F = open(fullname, 'w')
        sage: _ = F.write("k.<a> = GF(5^3); print(a^124)\n")
        sage: F.close()
        sage: (out, err, ret) = test_executable(["sage", fullname])
        sage: print(out)
        1
        sage: err
        ''
        sage: ret
        0
        sage: (out, err, ret) = test_executable(["sage", name], cwd=dir)
        sage: print(out)
        1
        sage: err
        ''
        sage: ret
        0

    Test running a ``.spyx`` file::

        sage: dir = tmp_dir(); name = 'sage_test_file.spyx'
        sage: fullname = os.path.join(dir, name)
        sage: F = open(fullname, 'w')
        sage: _ = F.write("from cysignals.signals cimport *\nfrom sage.rings.integer cimport Integer\ncdef long i, s = 0\nsig_on()\nfor i in range(1000): s += i\nsig_off()\nprint(Integer(s))")
        sage: F.close()
        sage: (out, err, ret) = test_executable(["sage", fullname], pydebug_ignore_warnings=True)
        sage: print(out)
        499500
        sage: import re
        sage: bool(re.match('Compiling.*spyx.*', err))
        True
        sage: ret
        0
        sage: (out, err, ret) = test_executable(["sage", name], cwd=dir, pydebug_ignore_warnings=True)
        sage: print(out)
        499500
        sage: bool(re.match('Compiling.*spyx.*', err))
        True
        sage: ret
        0

    Testing ``sage --preparse FILE`` and ``sage -t FILE``.  First create
    a file and preparse it::

        sage: s = "# -*- coding: utf-8 -*-\n'''This is a test file.\nAnd I am its doctest'''\ndef my_add(a):\n    '''\n    Add 2 to a.\n\n        EXAMPLES::\n\n            sage: my_add(2)\n            4\n        '''\n    return a + 2\n"
        sage: script = os.path.join(tmp_dir(), 'my_script.sage')
        sage: script_py = script + '.py'
        sage: F = open(script, 'w')
        sage: _ = F.write(s)
        sage: F.close()
        sage: (out, err, ret) = test_executable(["sage", "--preparse", script])
        sage: ret
        0
        sage: os.path.isfile(script_py)
        True

    Now test my_script.sage and the preparsed version my_script.sage.py::

        sage: (out, err, ret) = test_executable(["sage", "-t", script])
        sage: ret
        0
        sage: out.find("All tests passed!") >= 0
        True
        sage: (out, err, ret) = test_executable(["sage", "-t", script_py])
        sage: ret
        0
        sage: out.find("All tests passed!") >= 0
        True

    Test that the coding line and doctest are preserved::

        sage: Fpy = open(script_py, "r")
        sage: Fpy.readline()
        '# -*- coding: utf-8 -*-\n'
        sage: Fpy.readline()
        "'''This is a test file.\n"
        sage: Fpy.readline()
        "And I am its doctest'''\n"

    Now for a file which should fail tests::

        sage: s = s.replace('4', '5') # (2+2 != 5)
        sage: F = open(script, 'w')
        sage: _ = F.write(s)
        sage: F.close()
        sage: (out, err, ret) = test_executable(["sage", "-t", script])
        sage: ret
        1
        sage: out.find("1 item had failures:") >= 0
        True

    Test ``sage -t --debug -p 2`` on a ReST file, the ``-p 2`` should
    be ignored. In Pdb, we run the ``help`` command::

        sage: s = "::\n\n    sage: assert True is False\n    sage: 2 + 2\n    5"
        sage: script = tmp_filename(ext='.rst')
        sage: F = open(script, 'w')
        sage: _ = F.write(s)
        sage: F.close()
        sage: (out, err, ret) = test_executable([
        ....:     "sage", "-t", "--debug", "-p", "2", "--warn-long", "0", script], "help")
        sage: print(out)
        Debugging requires single-threaded operation, setting number of threads to 1.
        Running doctests with ID...
        Doctesting 1 file.
        sage -t ...
        **********************************************************************
        File "...", line 3, in ...
        Failed example:
            assert True is False
        Exception raised:
            Traceback (most recent call last):
            ...
            AssertionError
        > <doctest ...>(1)<module>()
        -> assert True is False
        (Pdb)
        Documented commands (type help <topic>):
        ========================================
        ...
        **********************************************************************
        File "...", line 4, in ...
        Failed example:
            2 + 2
        Expected:
            5
        Got:
            4
        **********************************************************************
        Previously executed commands:
            s...: assert True is False
        sage:
        <BLANKLINE>
        Returning to doctests...
        **********************************************************************
        1 item had failures:
           2 of   3 in ...
            [2 tests, 2 failures, ...]
        ...
        sage: ret
        1

    Now run a test for the fixdoctests script and, in particular, check that the
    issues raised in :trac:`10589` are fixed. We have to go to slightly silly
    lengths to doctest the output.::

        sage: test='r\"\"\"Add a doc-test for the fixdoctest command line option and, in particular, check that\n:trac:`10589` is fixed.\n\nEXAMPLES::\n\n    sage: 1+1              # incorrect output\n    3\n    sage: m=matrix(ZZ,3)   # output when none is expected\n    [0 0 0]\n    [0 0 0]\n    [1 0 0]\n    sage: (2/3)*m          # no output when it is expected\n    sage: mu=PartitionTuple([[4,4],[3,3,2,1],[1,1]])   # output when none is expected\n    [4, 4, 3, 3, 2, 1, 1]\n    sage: mu.pp()          # uneven indentation\n    ****\n    ****\n    sage: PartitionTuples.options(convention="French")\n    sage: mu.pp()         # fix doctest with uneven indentation\n    sage: PartitionTuples.options._reset()\n\"\"\"\n'
        sage: test_file = os.path.join(tmp_dir(), 'test_file.py')
        sage: F = open(test_file, 'w')
        sage: _ = F.write(test)
        sage: F.close()
        sage: (out, err, ret) = test_executable(["sage", "--fixdoctests", test_file])
        sage: with open(test_file, 'r') as f:
        ....:     fixed_test = f.read()
        sage: import difflib
        sage: list(difflib.unified_diff(test.splitlines(), fixed_test.splitlines()))[2:-1]
        ['@@ -4,18 +4,23 @@\n',
         ' EXAMPLES::',
         ' ',
         '     sage: 1+1              # incorrect output',
         '-    3',
         '+    2',
         '     sage: m=matrix(ZZ,3)   # output when none is expected',
         '+    sage: (2/3)*m          # no output when it is expected',
         '     [0 0 0]',
         '     [0 0 0]',
         '-    [1 0 0]',
         '-    sage: (2/3)*m          # no output when it is expected',
         '+    [0 0 0]',
         '     sage: mu=PartitionTuple([[4,4],[3,3,2,1],[1,1]])   # output when none is expected',
         '-    [4, 4, 3, 3, 2, 1, 1]',
         '     sage: mu.pp()          # uneven indentation',
         '-    ****',
         '-    ****',
         '+       ****   ***   *',
         '+       ****   ***   *',
         '+              **',
         '+              *',
         '     sage: PartitionTuples.options(convention="French")',
         '     sage: mu.pp()         # fix doctest with uneven indentation',
         '+    *',
         '+    **',
         '+    ****   ***   *',
         '+    ****   ***   *',
         '     sage: PartitionTuples.options._reset()']

    Test external programs being called by Sage::

        sage: (out, err, ret) = test_executable(["sage", "--sh"], "echo Hello World\nexit 42\n")
        sage: out.find("Hello World\n") >= 0
        True
        sage: ret
        42

        sage: (out, err, ret) = test_executable(["sage", "--sh", "-c", "echo Hello World; exit 42"])
        sage: out.find("Hello World\n") >= 0
        True
        sage: ret
        42

        sage: (out, err, ret) = test_executable(["sage", "--ipython"], "\n3**33\n", pydebug_ignore_warnings=True)
        sage: out.find("5559060566555523") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--python"], "print(3^33)\n")
        sage: out
        '34\n'
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--python3"], "print(3^33)\n")
        sage: out
        '34\n'
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--cython"])
        sage: print(err)
        Cython (http://cython.org) is a compiler for code written in the
        Cython language.  Cython is based on Pyrex by Greg Ewing.
        ...

        sage: def has_tty():
        ....:     try:
        ....:         os.open(os.ctermid(), os.O_RDONLY)
        ....:         return True
        ....:     except OSError:
        ....:         return False
        sage: (out, err, ret) = test_executable(["sage", "--ecl"], "(* 12345 54321)\n")
        sage: out.find("Embeddable Common-Lisp") >= 0
        True
        sage: out.find("670592745") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--lisp"], "(* 12345 54321)\n")
        sage: out.find("Embeddable Common-Lisp") >= 0
        True
        sage: out.find("670592745") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--gap", "-q"], "Size(SymmetricGroup(5));\n")
        sage: out
        '120\n'
        sage: err.replace('gap: halving pool size.', '').strip()
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--gdb"], 'quit\n')  # optional - gdb
        sage: out.find('(gdb) ') >= 0  # optional - gdb
        True
        sage: ret  # optional - gdb
        0

        sage: (out, err, ret) = test_executable(["sage", "--mwrank", "-v0", "-q", "-o"], "0 0 1 -7 6 0 0 0 0 0\n")
        sage: out
        'Curve [0,0,1,-7,6] :\tRank = 3\n[[3],[[1,-1],[-2,3],[-7/4,25/8]]]\n\n\n'
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--singular"], "12345*54321;\n")
        sage: out.find("A Computer Algebra System for Polynomial Computations") >= 0
        True
        sage: out.find("670592745") >= 0
        True
        sage: err
        ''
        sage: ret
        0

    Test GP using the ``-f`` option which prevents the reading of a ``.gprc``
    configuration file::

        sage: (out, err, ret) = test_executable(["sage", "--gp", "-f"], "3^33\nquit(42)\n")
        sage: out.find("PARI/GP") >= 0
        True
        sage: out.find("5559060566555523") >= 0
        True
        sage: err
        ''
        sage: ret
        42

    Some programs of which we check functionality using only ``--version``::

        sage: (out, err, ret) = test_executable(["sage", "--maxima", "--version"])
        sage: out.find("Maxima ") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--R", "--version"])  # optional - r
        sage: out.find("R version ") >= 0                                      # optional - r
        True
        sage: err                                                              # optional - r
        ''
        sage: ret                                                              # optional - r
        0

        sage: (out, err, ret) = test_executable(["sage", "--sqlite3", "--version"])
        sage: out.startswith("3.")
        True
        sage: err
        ''
        sage: ret
        0

    Check some things requiring an internet connection::

        sage: (out, err, ret) = test_executable(["sage", "--standard"])  # optional - internet
        sage: out.find("cython") >= 0  # optional - internet
        True
        sage: err  # optional - internet
        ''
        sage: ret  # optional - internet
        0

        sage: (out, err, ret) = test_executable(["sage", "--optional"])  # optional - internet
        sage: out.find("database_cremona_ellcurve") >= 0  # optional - internet
        True
        sage: err  # optional - internet
        ''
        sage: ret  # optional - internet
        0

        sage: (out, err, ret) = test_executable(["sage", "--experimental"])  # optional - internet
        sage: out.find("valgrind") >= 0  # optional - internet
        True
        sage: err  # optional - internet
        ''
        sage: ret  # optional - internet
        0

    Check an illegal command line option.  This outputs an error to stdout,
    but we allow stderr in case this changes in the future::

        sage: (out, err, ret) = test_executable(["sage", "--zzfoobar"])
        sage: (out+err).find("unknown option: --zzfoobar") >= 0
        True
        sage: ret > 0
        True

    Test ``sage --rst2ipynb file.rst`` on a ReST file::

        sage: s = "::\n\n    sage: 2^10\n    1024\n    sage: 2 + 2\n    4"
        sage: input = tmp_filename(ext='.rst')
        sage: with open(input, 'w') as F:
        ....:     _ = F.write(s)
        sage: L = ["sage", "--rst2ipynb", input]
        sage: (out, err, ret) = test_executable(L)           # optional - rst2ipynb
        sage: err                                            # optional - rst2ipynb
        ''
        sage: ret                                            # optional - rst2ipynb
        0
        sage: from json import loads                         # optional - rst2ipynb
        sage: d = loads(out)                                 # optional - rst2ipynb
        sage: sorted(d.keys())                               # optional - rst2ipynb
        ['cells', 'metadata', 'nbformat', 'nbformat_minor']
        sage: d['cells'][1]['source']                        # optional - rst2ipynb
        ['2^10']
        sage: d['cells'][2]['source']                        # optional - rst2ipynb
        ['2 + 2']

    Test ``sage --rst2ipynb file.rst file.ipynb`` on a ReST file::

        sage: s = "::\n\n    sage: 2^10\n    1024\n    sage: 2 + 2\n    4"
        sage: input = tmp_filename(ext='.rst')
        sage: output = tmp_filename(ext='.ipynb')
        sage: with open(input, 'w') as F:
        ....:     _ = F.write(s)
        sage: L = ["sage", "--rst2ipynb", input, output]
        sage: test_executable(L)                              # optional - rst2ipynb
        ('', '', 0)
        sage: import json                                     # optional - rst2ipynb
        sage: d = json.load(open(output,'r'))                 # optional - rst2ipynb
        sage: type(d)                                         # optional - rst2ipynb
        <class 'dict'>
        sage: sorted(d.keys())                                # optional - rst2ipynb
        ['cells', 'metadata', 'nbformat', 'nbformat_minor']
        sage: d['metadata']                                   # optional - rst2ipynb
        {'kernelspec': {'display_name': 'sagemath', 'name': 'sagemath'}}
        sage: d['cells'][1]['cell_type']                      # optional - rst2ipynb
        'code'

    Test ``sage --ipynb2rst file.ipynb file.rst`` on a ipynb file::

        sage: s = r'''{
        ....:  "cells": [
        ....:   {
        ....:    "cell_type": "code",
        ....:    "execution_count": 1,
        ....:    "metadata": {},
        ....:    "outputs": [
        ....:     {
        ....:      "data": {
        ....:       "text/plain": [
        ....:        "2"
        ....:       ]
        ....:      },
        ....:      "execution_count": 1,
        ....:      "metadata": {},
        ....:      "output_type": "execute_result"
        ....:     }
        ....:    ],
        ....:    "source": [
        ....:     "1+1"
        ....:    ]
        ....:   },
        ....:   {
        ....:    "cell_type": "code",
        ....:    "execution_count": null,
        ....:    "metadata": {},
        ....:    "outputs": [],
        ....:    "source": []
        ....:   }
        ....:  ],
        ....:  "metadata": {
        ....:   "kernelspec": {
        ....:    "display_name": "SageMath 8.3.beta4",
        ....:    "language": "",
        ....:    "name": "sagemath"
        ....:   },
        ....:   "language_info": {
        ....:    "codemirror_mode": {
        ....:     "name": "ipython",
        ....:     "version": 2
        ....:    },
        ....:    "file_extension": ".py",
        ....:    "mimetype": "text/x-python",
        ....:    "name": "python",
        ....:    "nbconvert_exporter": "python",
        ....:    "pygments_lexer": "ipython2",
        ....:    "version": "2.7.15"
        ....:   }
        ....:  },
        ....:  "nbformat": 4,
        ....:  "nbformat_minor": 2
        ....: }
        ....: '''
        sage: t = '.. escape-backslashes\n.. default-role:: math\n\n\n::\n\n    sage: 1+1\n    2\n\n\n\n\n'
        sage: input = tmp_filename(ext='.ipynb')
        sage: output = tmp_filename(ext='.rst')
        sage: with open(input, 'w') as F:
        ....:     _ = F.write(s)
        sage: L = ["sage", "--ipynb2rst", input, output]
        sage: _ = test_executable(L)                        # optional - pandoc
        sage: print(open(output, 'r').read() == t)          # optional - pandoc  # known bug #32697
        True
    """
    pexpect_env = dict(os.environ)
    try:
        del pexpect_env["TERM"]
    except KeyError:
        pass

    __with_pydebug = hasattr(sys, 'gettotalrefcount')   # This is a Python debug build (--with-pydebug) 
    if __with_pydebug and pydebug_ignore_warnings:
        pexpect_env['PYTHONWARNINGS'] = ','.join([
            'ignore::DeprecationWarning',
        ])

    kwds['encoding'] = kwds.pop('encoding', 'utf-8')

    p = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE, env=pexpect_env,
              **kwds)
    if input:
        p.stdin.write(input)

    p.stdin.close()
    fdout = p.stdout.fileno()
    fderr = p.stderr.fileno()
    out = []
    err = []

    while True:
        # Try reading from fdout and fderr
        rfd = []
        if fdout:
            rfd.append(fdout)
        if fderr:
            rfd.append(fderr)
        if len(rfd) == 0:
            break
        timeout = float(timeout)
        rlist = select.select(rfd, [], [], timeout)[0]

        if len(rlist) == 0:
            # Timeout!
            p.terminate()
            raise RuntimeError("timeout in test_executable()")
        if fdout in rlist:
            s = p.stdout.read(1024)
            if not s:
                fdout = None   # EOF
                p.stdout.close()
            out.append(s)
        if fderr in rlist:
            s = p.stderr.read(1024)
            if not s:
                fderr = None   # EOF
                p.stderr.close()
            err.append(s)

    return (''.join(out), ''.join(err), p.wait())
