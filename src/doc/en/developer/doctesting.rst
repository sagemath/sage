.. nodoctest

.. _chapter-doctesting:

===========================
Doctesting the Sage Library
===========================

Doctesting a function ensures that the function performs as claimed by
its documentation. Testing can be performed using one thread or
multiple threads. After compiling a source version of Sage, doctesting
can be run on the whole Sage library, on all modules under a given
directory, or on a specified module only. For the purposes of this
chapter, suppose we have compiled Sage 5.9 from source and the top
level Sage directory is::

    [jdemeyer@sage sage-6.0]$ pwd
    /scratch/jdemeyer/build/sage-6.0

See the section :ref:`chapter-testing` for information on Sage's
automated testing process. The general syntax for doctesting is as
follows. To doctest a module in the library of a version of Sage, use
this syntax::

    /path/to/sage-x.y.z/sage -t [--long] /path/to/sage-x.y.z/path/to/module.py[x]

where ``--long`` is an optional argument (see :ref:`section-options`
for more options). The version of ``sage`` used must match the version
of Sage containing the module we want to doctest. A Sage module can be
either a Python script (with the file extension ".py") or it can be a
Cython script, in which case it has the file extension ".pyx".


Testing a Module
================

Say we want to run all tests in the sudoku module
``sage/games/sudoku.py``. In a terminal window, first we ``cd`` to the
top level Sage directory of our local Sage installation. Now  we can
start doctesting as demonstrated in the following terminal session::

    [jdemeyer@sage sage-6.0]$ ./sage -t src/sage/games/sudoku.py
    Running doctests with ID 2012-07-03-03-36-49-d82849c6.
    Doctesting 1 file.
    sage -t src/sage/games/sudoku.py
        [103 tests, 3.6 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 4.8 seconds
        cpu time: 3.6 seconds
        cumulative wall time: 3.6 seconds

The numbers output by the test show that testing the sudoku module
takes about four seconds, while testing all specified modules took the
same amount of time; the total time required includes some startup
time for the code that runs the tests. In this case, we only tested
one module so it is not surprising that the total testing time is
approximately the same as the time required to test only that one
module. Notice that the syntax is::

    [jdemeyer@sage sage-6.0]$ ./sage -t src/sage/games/sudoku.py
    Running doctests with ID 2012-07-03-03-39-02-da6accbb.
    Doctesting 1 file.
    sage -t src/sage/games/sudoku.py
        [103 tests, 3.6 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 4.9 seconds
        cpu time: 3.6 seconds
        cumulative wall time: 3.6 seconds

but not::

    [jdemeyer@sage sage-6.0]$ ./sage -t sage/games/sudoku.py
    Running doctests with ID 2012-07-03-03-40-53-6cc4f29f.
    No files matching sage/games/sudoku.py
    No files to doctest

We can also first ``cd`` to the directory containing the module
``sudoku.py`` and doctest that module as follows::

    [jdemeyer@sage sage-6.0]$ cd src/sage/games/
    [jdemeyer@sage games]$ ls
    __init__.py  hexad.py       sudoku.py           sudoku_backtrack.pyx
    all.py       quantumino.py  sudoku_backtrack.c
    [jdemeyer@sage games]$ ../../../../sage -t sudoku.py
    Running doctests with ID 2012-07-03-03-41-39-95ebd2ff.
    Doctesting 1 file.
    sage -t sudoku.py
        [103 tests, 3.6 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 5.2 seconds
        cpu time: 3.6 seconds
        cumulative wall time: 3.6 seconds

In all of the above terminal sessions, we used a local installation of
Sage to test its own modules. Even if we have a system-wide Sage
installation, using that version to doctest the modules of a local
installation is a recipe for confusion.


Troubleshooting
===============

To doctest modules of a Sage installation, from a terminal window we
first ``cd`` to the top level directory of that Sage installation,
otherwise known as the ``SAGE_ROOT`` of that installation. When we
run tests, we use that particular Sage installation via the syntax
``./sage``; notice the "dot-forward-slash" at the front of
``sage``. This is a precaution against confusion that can arise when
our system has multiple Sage installations. For example, the following
syntax is acceptable because we explicitly specify the Sage
installation in the current ``SAGE_ROOT``::

    [jdemeyer@sage sage-6.0]$ ./sage -t src/sage/games/sudoku.py
    Running doctests with ID 2012-07-03-03-43-24-a3449f54.
    Doctesting 1 file.
    sage -t src/sage/games/sudoku.py
        [103 tests, 3.6 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 4.9 seconds
        cpu time: 3.6 seconds
        cumulative wall time: 3.6 seconds
    [jdemeyer@sage sage-6.0]$ ./sage -t "src/sage/games/sudoku.py"
    Running doctests with ID 2012-07-03-03-43-54-ac8ca007.
    Doctesting 1 file.
    sage -t src/sage/games/sudoku.py
        [103 tests, 3.6 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 4.9 seconds
        cpu time: 3.6 seconds
        cumulative wall time: 3.6 seconds

The following syntax is not recommended as we are using a system-wide
Sage installation (if it exists):

.. skip

::

    [jdemeyer@sage sage-6.0]$ sage -t src/sage/games/sudoku.py
    sage -t  "src/sage/games/sudoku.py"
    **********************************************************************
    File "/home/jdemeyer/sage/sage-6.0/src/sage/games/sudoku.py", line 515:
        sage: h.solve(algorithm='backtrack').next()
    Exception raised:
        Traceback (most recent call last):
          File "/usr/local/sage/local/bin/ncadoctest.py", line 1231, in run_one_test
            self.run_one_example(test, example, filename, compileflags)
          File "/usr/local/sage/local/bin/sagedoctest.py", line 38, in run_one_example
            OrigDocTestRunner.run_one_example(self, test, example, filename, compileflags)
          File "/usr/local/sage/local/bin/ncadoctest.py", line 1172, in run_one_example
            compileflags, 1) in test.globs
          File "<doctest __main__.example_13[4]>", line 1, in <module>
            h.solve(algorithm='backtrack').next()###line 515:
        sage: h.solve(algorithm='backtrack').next()
          File "/home/jdemeyer/.sage/tmp/sudoku.py", line 607, in solve
            for soln in gen:
          File "/home/jdemeyer/.sage/tmp/sudoku.py", line 719, in backtrack
            from sudoku_backtrack import backtrack_all
        ImportError: No module named sudoku_backtrack
    **********************************************************************
    [...more errors...]
    2 items had failures:
       4 of  15 in __main__.example_13
       2 of   8 in __main__.example_14
    ***Test Failed*** 6 failures.
    For whitespace errors, see the file /home/jdemeyer/.sage//tmp/.doctest_sudoku.py
             [21.1 s]

    ----------------------------------------------------------------------
    The following tests failed:


            sage -t  "src/sage/games/sudoku.py"
    Total time for all tests: 21.3 seconds

In this case, we received an error because the system-wide Sage
installation is a different (older) version than the one we are
using for Sage development.  Make sure you always test the files
with the correct version of Sage.

Parallel Testing Many Modules
=============================

So far we have used a single thread to doctest a module in the Sage
library. There are hundreds, even thousands of modules in the Sage
library. Testing them all using one thread would take a few
hours. Depending on our hardware, this could take up to six hours or
more. On a multi-core system, parallel doctesting can significantly
reduce the testing time. Unless we also want to use our computer
while doctesting in parallel, we can choose to devote all the cores
of our system for parallel testing.

Let us doctest all modules in a directory, first using a single thread
and then using four threads. For this example, suppose we want to test
all the modules under ``sage/crypto/``. We can use a syntax similar to
that shown above to achieve this::

    [jdemeyer@sage sage-6.0]$ ./sage -t src/sage/crypto
    Running doctests with ID 2012-07-03-03-45-40-7f837dcf.
    Doctesting 24 files.
    sage -t src/sage/crypto/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/boolean_function.pyx
        [252 tests, 4.4 s]
    sage -t src/sage/crypto/cipher.py
        [10 tests, 0.0 s]
    sage -t src/sage/crypto/classical.py
        [718 tests, 11.3 s]
    sage -t src/sage/crypto/classical_cipher.py
        [130 tests, 0.5 s]
    sage -t src/sage/crypto/cryptosystem.py
        [82 tests, 0.1 s]
    sage -t src/sage/crypto/lattice.py
        [1 tests, 0.0 s]
    sage -t src/sage/crypto/lfsr.py
        [31 tests, 0.1 s]
    sage -t src/sage/crypto/stream.py
        [17 tests, 0.1 s]
    sage -t src/sage/crypto/stream_cipher.py
        [114 tests, 0.2 s]
    sage -t src/sage/crypto/util.py
        [122 tests, 0.2 s]
    sage -t src/sage/crypto/block_cipher/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/block_cipher/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/block_cipher/miniaes.py
        [430 tests, 1.3 s]
    sage -t src/sage/crypto/block_cipher/sdes.py
        [290 tests, 0.9 s]
    sage -t src/sage/crypto/mq/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/mq/mpolynomialsystem.py
        [320 tests, 9.1 s]
    sage -t src/sage/crypto/mq/mpolynomialsystemgenerator.py
        [42 tests, 0.1 s]
    sage -t src/sage/crypto/mq/sbox.py
        [124 tests, 0.8 s]
    sage -t src/sage/crypto/mq/sr.py
        [435 tests, 5.5 s]
    sage -t src/sage/crypto/public_key/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/public_key/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/public_key/blum_goldwasser.py
        [135 tests, 0.2 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 38.1 seconds
        cpu time: 29.8 seconds
        cumulative wall time: 35.1 seconds

Now we do the same thing, but this time we also use the optional
argument ``--long``::

    [jdemeyer@sage sage-6.0]$ ./sage -t --long src/sage/crypto/
    Running doctests with ID 2012-07-03-03-48-11-c16721e6.
    Doctesting 24 files.
    sage -t --long src/sage/crypto/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/all.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/boolean_function.pyx
        [252 tests, 4.2 s]
    sage -t --long src/sage/crypto/cipher.py
        [10 tests, 0.0 s]
    sage -t --long src/sage/crypto/classical.py
        [718 tests, 10.3 s]
    sage -t --long src/sage/crypto/classical_cipher.py
        [130 tests, 0.5 s]
    sage -t --long src/sage/crypto/cryptosystem.py
        [82 tests, 0.1 s]
    sage -t --long src/sage/crypto/lattice.py
        [1 tests, 0.0 s]
    sage -t --long src/sage/crypto/lfsr.py
        [31 tests, 0.1 s]
    sage -t --long src/sage/crypto/stream.py
        [17 tests, 0.1 s]
    sage -t --long src/sage/crypto/stream_cipher.py
        [114 tests, 0.2 s]
    sage -t --long src/sage/crypto/util.py
        [122 tests, 0.2 s]
    sage -t --long src/sage/crypto/block_cipher/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/block_cipher/all.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/block_cipher/miniaes.py
        [430 tests, 1.1 s]
    sage -t --long src/sage/crypto/block_cipher/sdes.py
        [290 tests, 0.7 s]
    sage -t --long src/sage/crypto/mq/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/mq/mpolynomialsystem.py
        [320 tests, 7.5 s]
    sage -t --long src/sage/crypto/mq/mpolynomialsystemgenerator.py
        [42 tests, 0.1 s]
    sage -t --long src/sage/crypto/mq/sbox.py
        [124 tests, 0.7 s]
    sage -t --long src/sage/crypto/mq/sr.py
        [437 tests, 82.4 s]
    sage -t --long src/sage/crypto/public_key/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/public_key/all.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/public_key/blum_goldwasser.py
        [135 tests, 0.2 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 111.8 seconds
        cpu time: 106.1 seconds
        cumulative wall time: 108.5 seconds

Notice the time difference between the first set of tests and the
second set, which uses the optional argument ``--long``. Many tests in the
Sage library are flagged with ``# long time`` because these are known to
take a long time to run through. Without using the optional ``--long``
argument, the module ``sage/crypto/mq/sr.py`` took about five
seconds. With this optional argument, it required 82 seconds to run
through all tests in that module. Here is a snippet of a function in
the module ``sage/crypto/mq/sr.py`` with a doctest that has been flagged
as taking a long time::

    def test_consistency(max_n=2, **kwargs):
        r"""
        Test all combinations of ``r``, ``c``, ``e`` and ``n`` in ``(1,
        2)`` for consistency of random encryptions and their polynomial
        systems. `\GF{2}` and `\GF{2^e}` systems are tested. This test takes
        a while.

        INPUT:

        - ``max_n`` -- maximal number of rounds to consider (default: 2)
        - ``kwargs`` -- are passed to the SR constructor

        TESTS:

        The following test called with ``max_n`` = 2 requires a LOT of RAM
        (much more than 2GB).  Since this might cause the doctest to fail
        on machines with "only" 2GB of RAM, we test ``max_n`` = 1, which
        has a more reasonable memory usage. ::

            sage: from sage.crypto.mq.sr import test_consistency
            sage: test_consistency(1)  # long time (80s on sage.math, 2011)
            True
        """

Now we doctest the same directory in parallel using 4 threads::

    [jdemeyer@sage sage-6.0]$ ./sage -tp 4 src/sage/crypto/
    Running doctests with ID 2012-07-07-00-11-55-9b17765e.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 24 files using 4 threads.
    sage -t src/sage/crypto/boolean_function.pyx
        [252 tests, 3.8 s]
    sage -t src/sage/crypto/block_cipher/miniaes.py
        [429 tests, 1.1 s]
    sage -t src/sage/crypto/mq/sr.py
        [432 tests, 5.7 s]
    sage -t src/sage/crypto/mq/sbox.py
        [123 tests, 0.8 s]
    sage -t src/sage/crypto/block_cipher/sdes.py
        [289 tests, 0.6 s]
    sage -t src/sage/crypto/classical_cipher.py
        [123 tests, 0.4 s]
    sage -t src/sage/crypto/stream_cipher.py
        [113 tests, 0.1 s]
    sage -t src/sage/crypto/public_key/blum_goldwasser.py
        [134 tests, 0.1 s]
    sage -t src/sage/crypto/lfsr.py
        [30 tests, 0.1 s]
    sage -t src/sage/crypto/util.py
        [121 tests, 0.1 s]
    sage -t src/sage/crypto/cryptosystem.py
        [79 tests, 0.0 s]
    sage -t src/sage/crypto/stream.py
        [12 tests, 0.0 s]
    sage -t src/sage/crypto/mq/mpolynomialsystemgenerator.py
        [40 tests, 0.0 s]
    sage -t src/sage/crypto/cipher.py
        [3 tests, 0.0 s]
    sage -t src/sage/crypto/lattice.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/block_cipher/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/public_key/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/public_key/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/mq/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/block_cipher/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/mq/mpolynomialsystem.py
        [318 tests, 8.4 s]
    sage -t src/sage/crypto/classical.py
        [717 tests, 10.4 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 12.9 seconds
        cpu time: 30.5 seconds
        cumulative wall time: 31.7 seconds
    [jdemeyer@sage sage-6.0]$ ./sage -tp 4 --long src/sage/crypto/
    Running doctests with ID 2012-07-07-00-13-04-d71f3cd4.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 24 files using 4 threads.
    sage -t --long src/sage/crypto/boolean_function.pyx
        [252 tests, 3.7 s]
    sage -t --long src/sage/crypto/block_cipher/miniaes.py
        [429 tests, 1.0 s]
    sage -t --long src/sage/crypto/mq/sbox.py
        [123 tests, 0.8 s]
    sage -t --long src/sage/crypto/block_cipher/sdes.py
        [289 tests, 0.6 s]
    sage -t --long src/sage/crypto/classical_cipher.py
        [123 tests, 0.4 s]
    sage -t --long src/sage/crypto/util.py
        [121 tests, 0.1 s]
    sage -t --long src/sage/crypto/stream_cipher.py
        [113 tests, 0.1 s]
    sage -t --long src/sage/crypto/public_key/blum_goldwasser.py
        [134 tests, 0.1 s]
    sage -t --long src/sage/crypto/lfsr.py
        [30 tests, 0.0 s]
    sage -t --long src/sage/crypto/cryptosystem.py
        [79 tests, 0.0 s]
    sage -t --long src/sage/crypto/stream.py
        [12 tests, 0.0 s]
    sage -t --long src/sage/crypto/mq/mpolynomialsystemgenerator.py
        [40 tests, 0.0 s]
    sage -t --long src/sage/crypto/cipher.py
        [3 tests, 0.0 s]
    sage -t --long src/sage/crypto/lattice.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/block_cipher/all.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/public_key/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/mq/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/all.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/block_cipher/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/public_key/all.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/mq/mpolynomialsystem.py
        [318 tests, 9.0 s]
    sage -t --long src/sage/crypto/classical.py
        [717 tests, 10.5 s]
    sage -t --long src/sage/crypto/mq/sr.py
        [434 tests, 88.0 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 90.4 seconds
        cpu time: 113.4 seconds
        cumulative wall time: 114.5 seconds

As the number of threads increases, the total testing time
decreases.


.. _section-parallel-test-whole-library:

Parallel Testing the Whole Sage Library
=======================================

The main Sage library resides in the directory
``SAGE_ROOT/src/``. We can use the syntax described above
to doctest the main library using multiple threads. When doing release
management or patching the main Sage library, a release manager would
parallel test the library using 10 threads with the following command::

    [jdemeyer@sage sage-6.0]$ ./sage -tp 10 --long src/

Another way is run ``make ptestlong``, which builds Sage (if necessary),
builds the Sage documentation (if necessary), and then runs parallel
doctests.  This determines the number of threads by reading the
environment variable :envvar:`MAKE`: if it is set to ``make -j12``, then
use 12 threads.  If :envvar:`MAKE` is not set, then by default it uses
the number of CPU cores (as determined by the Python function
``multiprocessing.cpu_count()``) with a minimum of 2 and a maximum of 8.

In any case, this will test the Sage library with multiple threads::

    [jdemeyer@sage sage-6.0]$ make ptestlong

Any of the following commands would also doctest the Sage library or
one of its clones::

    make test
    make check
    make testlong
    make ptest
    make ptestlong

The differences are:

* ``make test`` and ``make check`` --- These two commands run the same
  set of tests. First the Sage standard documentation is tested,
  i.e. the documentation that resides in

  * ``SAGE_ROOT/src/doc/common``
  * ``SAGE_ROOT/src/doc/en``
  * ``SAGE_ROOT/src/doc/fr``

  Finally, the commands doctest the Sage library. For more details on
  these command, see the file ``SAGE_ROOT/Makefile``.

* ``make testlong`` --- This command doctests the standard
  documentation:

  * ``SAGE_ROOT/src/doc/common``
  * ``SAGE_ROOT/src/doc/en``
  * ``SAGE_ROOT/src/doc/fr``

  and then the Sage library. Doctesting is run with the optional
  argument ``--long``. See the file ``SAGE_ROOT/Makefile`` for further
  details.

* ``make ptest`` --- Similar to the commands ``make test`` and ``make
  check``. However, doctesting is run with the number of threads as
  described above for ``make ptestlong``.

* ``make ptestlong`` --- Similar to the command ``make ptest``, but
  using the optional argument ``--long`` for doctesting.


Beyond the Sage Library
=======================

Doctesting also works fine for files not in the Sage library.  For
example, suppose we have a Python script called
``my_python_script.py``::

    [mvngu@sage build]$ cat my_python_script.py
    from sage.all_cmdline import *   # import sage library

    def square(n):
        """
        Return the square of n.

        EXAMPLES::

            sage: square(2)
            4
        """
        return n**2

Then we can doctest it just as with Sage library files::

    [mvngu@sage sage-6.0]$ ./sage -t my_python_script.py
    Running doctests with ID 2012-07-07-00-17-56-d056f7c0.
    Doctesting 1 file.
    sage -t my_python_script.py
        [1 test, 0.0 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 2.2 seconds
        cpu time: 0.0 seconds
        cumulative wall time: 0.0 seconds

Doctesting can also be performed on Sage scripts. Say we have a Sage
script called ``my_sage_script.sage`` with the following content::

    [mvngu@sage sage-6.0]$ cat my_sage_script.sage
    def cube(n):
        r"""
        Return the cube of n.

        EXAMPLES::

            sage: cube(2)
            8
        """
        return n**3

Then we can doctest it just as for Python files::

    [mvngu@sage build]$ sage-6.0/sage -t my_sage_script.sage
    Running doctests with ID 2012-07-07-00-20-06-82ee728c.
    Doctesting 1 file.
    sage -t my_sage_script.sage
        [1 test, 0.0 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 2.5 seconds
        cpu time: 0.0 seconds
        cumulative wall time: 0.0 seconds

Alternatively, we can preparse it to convert it to a Python script,
and then doctest that::

    [mvngu@sage build]$ sage-6.0/sage --preparse my_sage_script.sage
    [mvngu@sage build]$ cat my_sage_script.py
    # This file was *autogenerated* from the file my_sage_script.sage.
    from sage.all_cmdline import *   # import sage library
    _sage_const_3 = Integer(3)
    def cube(n):
        r"""
        Return the cube of n.

        EXAMPLES::

            sage: cube(2)
            8
        """
        return n**_sage_const_3
    [mvngu@sage build]$ sage-6.0/sage -t my_sage_script.py
    Running doctests with ID 2012-07-07-00-26-46-2bb00911.
    Doctesting 1 file.
    sage -t my_sage_script.py
        [2 tests, 0.0 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 2.3 seconds
        cpu time: 0.0 seconds
        cumulative wall time: 0.0 seconds

Doctesting from Within Sage
===========================

You can run doctests from within Sage, which can be useful since you
don't have to wait for Sage to start.  Use the ``run_doctests``
function in the global namespace, passing it either a string or a module::

    sage: run_doctests(sage.coding.sd_codes)
    Doctesting /Users/roed/sage/sage-5.3/src/sage/coding/sd_codes.py
    Running doctests with ID 2012-07-07-04-32-36-81f3853b.
    Doctesting 1 file.
    sage -t /Users/roed/sage/sage-5.3/src/sage/coding/sd_codes.py
        [18 tests, 0.3 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 0.4 seconds
        cpu time: 0.2 seconds
        cumulative wall time: 0.3 seconds

.. _section-options:

Optional Arguments
==================

Run Long Tests
--------------

Use the ``--long`` flag to run doctests that have been marked with
the comment ``# long time``.

No doctest should take longer than a second or so, and longer doctests
(taking up to 30-60 seconds) should be marked as ``# long time``.
These tests are normally skipped in order to reduce the time spent
running tests::

    [roed@sage sage-6.0]$ sage -t src/sage/rings/tests.py
    Running doctests with ID 2012-06-21-16-00-13-40835825.
    Doctesting 1 file.
    sage -t tests.py
        [18 tests, 1.1 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 2.9 seconds
        cpu time: 0.9 seconds
        cumulative wall time: 1.1 seconds

In order to run the long tests as well, do the following::

    [roed@sage sage-6.0]$ sage -t --long src/sage/rings/tests.py
    Running doctests with ID 2012-06-21-16-02-05-d13a9a24.
    Doctesting 1 file.
    sage -t tests.py
        [20 tests, 34.7 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 46.5 seconds
        cpu time: 25.2 seconds
        cumulative wall time: 34.7 seconds

To find tests that take longer than the allowed time use the
``--warn-long`` flag.  Without any options it will cause tests to fail
if they take longer than 1.0 second::

    [roed@sage sage-6.0]$ sage -t --warn-long src/sage/rings/factorint.pyx
    Running doctests with ID 2012-07-14-03-27-03-2c952ac1.
    Doctesting 1 file.
    sage -t --warn-long src/sage/rings/factorint.pyx
    **********************************************************************
    File "src/sage/rings/factorint.pyx", line 125, in sage.rings.factorint.base_exponent
    Failed example:
        base_exponent(-4)
    Test ran for 4.09 s
    **********************************************************************
    File "src/sage/rings/factorint.pyx", line 153, in sage.rings.factorint.factor_aurifeuillian
    Failed example:
        fa(2^6+1)
    Test ran for 2.22 s
    **********************************************************************
    File "src/sage/rings/factorint.pyx", line 155, in sage.rings.factorint.factor_aurifeuillian
    Failed example:
        fa(2^58+1)
    Test ran for 2.22 s
    **********************************************************************
    File "src/sage/rings/factorint.pyx", line 163, in sage.rings.factorint.factor_aurifeuillian
    Failed example:
        fa(2^4+1)
    Test ran for 2.25 s
    **********************************************************************
    2 items had failures:
       1 of   6 in sage.rings.factorint.base_exponent
       3 of   8 in sage.rings.factorint.factor_aurifeuillian
        [25 tests, 4 failures, 10.9 s]
    ------------------------------------------------------------------------
    sage -t --warn-long src/sage/rings/factorint.pyx # 4 doctests failed
    ------------------------------------------------------------------------
    Total time for all tests: 16.1 seconds
        cpu time: 9.7 seconds
        cumulative wall time: 10.9 seconds

You can also pass in an explicit amount of time::

    [roed@sage sage-6.0]$ sage -t --long --warn-long 2.0 src/sage/rings/tests.py
    Running doctests with ID 2012-07-14-03-30-13-c9164c9d.
    Doctesting 1 file.
    sage -t --long --warn-long 2.0 tests.py
    **********************************************************************
    File "tests.py", line 240, in sage.rings.tests.test_random_elements
    Failed example:
        sage.rings.tests.test_random_elements(trials=1000)  # long time (5 seconds)
    Test ran for 13.36 s
    **********************************************************************
    File "tests.py", line 283, in sage.rings.tests.test_random_arith
    Failed example:
        sage.rings.tests.test_random_arith(trials=1000)   # long time (5 seconds?)
    Test ran for 12.42 s
    **********************************************************************
    2 items had failures:
       1 of   4 in sage.rings.tests.test_random_arith
       1 of   4 in sage.rings.tests.test_random_elements
        [20 tests, 2 failures, 26.3 s]
    ------------------------------------------------------------------------
    sage -t --long --warn-long 2.0 tests.py # 2 doctests failed
    ------------------------------------------------------------------------
    Total time for all tests: 27.6 seconds
        cpu time: 24.8 seconds
        cumulative wall time: 26.3 seconds

Run Optional Tests
------------------

You can run tests that require optional packages by using the
``--optional`` flag.  Obviously, you need to have installed the
necessary optional packages in order for these tests to succeed.  See
http://www.sagemath.org/packages/optional/ in order to download
optional packages.

By default, Sage only runs doctests that are not marked with the ``optional`` tag.  This is equivalent to running ::

    [roed@sage sage-6.0]$ sage -t --optional=sage src/sage/rings/real_mpfr.pyx
    Running doctests with ID 2012-06-21-16-18-30-a368a200.
    Doctesting 1 file.
    sage -t src/sage/rings/real_mpfr.pyx
        [819 tests, 7.0 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 8.4 seconds
        cpu time: 4.1 seconds
        cumulative wall time: 7.0 seconds

If you want to also run tests that require magma, you can do the following::

    [roed@sage sage-6.0]$ sage -t --optional=sage,magma src/sage/rings/real_mpfr.pyx
    Running doctests with ID 2012-06-21-16-18-30-a00a7319
    Doctesting 1 file.
    sage -t src/sage/rings/real_mpfr.pyx
        [823 tests, 8.4 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 9.6 seconds
        cpu time: 4.0 seconds
        cumulative wall time: 8.4 seconds

In order to just run the tests that are marked as requiring magma, omit ``sage``::

    [roed@sage sage-6.0]$ sage -t --optional=magma src/sage/rings/real_mpfr.pyx
    Running doctests with ID 2012-06-21-16-18-33-a2bc1fdf
    Doctesting 1 file.
    sage -t src/sage/rings/real_mpfr.pyx
        [4 tests, 2.0 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 3.2 seconds
        cpu time: 0.1 seconds
        cumulative wall time: 2.0 seconds

To run all tests, regardless of whether they are marked optional, pass ``all`` as the ``optional`` tag::

    [roed@sage sage-6.0]$ sage -t --optional=all src/sage/rings/real_mpfr.pyx
    Running doctests with ID 2012-06-21-16-31-18-8c097f55
    Doctesting 1 file.
    sage -t src/sage/rings/real_mpfr.pyx
        [865 tests, 11.2 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 12.8 seconds
        cpu time: 4.7 seconds
        cumulative wall time: 11.2 seconds

Running Tests in Parallel
-------------------------

If you're testing many files, you can get big speedups by using more
than one thread.  To run doctests in parallel use the ``--nthreads``
flag (``-p`` is a shortened version).  Pass in the number of threads
you would like to use (by default Sage just uses 1)::

    [roed@sage sage-6.0]$ sage -tp 2 src/sage/doctest/
    Running doctests with ID 2012-06-22-19-09-25-a3afdb8c.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 8 files using 2 threads.
    sage -t src/sage/doctest/control.py
        [114 tests, 4.6 s]
    sage -t src/sage/doctest/util.py
        [114 tests, 0.6 s]
    sage -t src/sage/doctest/parsing.py
        [187 tests, 0.5 s]
    sage -t src/sage/doctest/sources.py
        [128 tests, 0.1 s]
    sage -t src/sage/doctest/reporting.py
        [53 tests, 0.1 s]
    sage -t src/sage/doctest/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/doctest/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/doctest/forker.py
        [322 tests, 15.5 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 17.0 seconds
        cpu time: 4.2 seconds
        cumulative wall time: 21.5 seconds

Doctesting All of Sage
----------------------

To doctest the whole Sage library use the ``--all`` flag (``-a`` for
short).  In addition to testing the code in Sage's Python and Cython
files, this command will run the tests defined in Sage's documentation
as well as testing the Sage notebook::

    [roed@sage sage-6.0]$ sage -t -a
    Running doctests with ID 2012-06-22-19-10-27-e26fce6d.
    Doctesting entire Sage library.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 2020 files.
    sage -t /Users/roed/sage/sage-5.3/src/sage/plot/plot.py
        [304 tests, 69.0 s]
    ...

If you want to just run the notebook tests, use the ``--sagenb`` flag instead.


Debugging Tools
---------------

Sometimes doctests fail (that's why we run them after all).  There are
various flags to help when something goes wrong.  If a doctest
produces a Python error, then normally tests continue after reporting
that an error occurred.  If you use the flag ``--debug`` (``-d`` for
short) then you will drop into an interactive Python debugger whenever
a Python exception occurs.  As an example, I modified
:mod:`sage.schemes.elliptic_curves.constructor` to produce an error::

    [roed@sage sage-6.0]$ sage -t --debug src/sage/schemes/elliptic_curves/constructor.py
    Running doctests with ID 2012-06-23-12-09-04-b6352629.
    Doctesting 1 file.
    **********************************************************************
    File "sage.schemes.elliptic_curves.constructor", line 4, in sage.schemes.elliptic_curves.constructor
    Failed example:
        EllipticCurve([0,0])
    Exception raised:
        Traceback (most recent call last):
          File "/Users/roed/sage/sage-5.3/local/lib/python2.7/site-packages/sage/doctest/forker.py", line 573, in _run
            self.execute(example, compiled, test.globs)
          File "/Users/roed/sage/sage-5.3/local/lib/python2.7/site-packages/sage/doctest/forker.py", line 835, in execute
            exec compiled in globs
          File "<doctest sage.schemes.elliptic_curves.constructor[0]>", line 1, in <module>
            EllipticCurve([Integer(0),Integer(0)])
          File "/Users/roed/sage/sage-5.3/local/lib/python2.7/site-packages/sage/schemes/elliptic_curves/constructor.py", line 346, in EllipticCurve
            return ell_rational_field.EllipticCurve_rational_field(x, y)
          File "/Users/roed/sage/sage-5.3/local/lib/python2.7/site-packages/sage/schemes/elliptic_curves/ell_rational_field.py", line 216, in __init__
            EllipticCurve_number_field.__init__(self, Q, ainvs)
          File "/Users/roed/sage/sage-5.3/local/lib/python2.7/site-packages/sage/schemes/elliptic_curves/ell_number_field.py", line 159, in __init__
            EllipticCurve_field.__init__(self, [field(x) for x in ainvs])
          File "/Users/roed/sage/sage-5.3/local/lib/python2.7/site-packages/sage/schemes/elliptic_curves/ell_generic.py", line 156, in __init__
            "Invariants %s define a singular curve."%ainvs
        ArithmeticError: Invariants [0, 0, 0, 0, 0] define a singular curve.
    > /Users/roed/sage/sage-5.3/local/lib/python2.7/site-packages/sage/schemes/elliptic_curves/ell_generic.py(156)__init__()
    -> "Invariants %s define a singular curve."%ainvs
    (Pdb) l
    151                 if len(ainvs) == 2:
    152                     ainvs = [K(0),K(0),K(0)] + ainvs
    153                 self.__ainvs = tuple(ainvs)
    154                 if self.discriminant() == 0:
    155                     raise ArithmeticError, \
    156  ->                       "Invariants %s define a singular curve."%ainvs
    157                 PP = projective_space.ProjectiveSpace(2, K, names='xyz');
    158                 x, y, z = PP.coordinate_ring().gens()
    159                 a1, a2, a3, a4, a6 = ainvs
    160                 f = y**2*z + (a1*x + a3*z)*y*z \
    161                     - (x**3 + a2*x**2*z + a4*x*z**2 + a6*z**3)
    (Pdb) p ainvs
    [0, 0, 0, 0, 0]
    (Pdb) quit
    **********************************************************************
    1 items had failures:
       1 of   1 in sage.schemes.elliptic_curves.constructor
    ***Test Failed*** 1 failures.
    sage -t src/sage/schemes/elliptic_curves/constructor.py
        [64 tests, 89.2 s]
    ------------------------------------------------------------------------
    sage -t src/sage/schemes/elliptic_curves/constructor.py # 1 doctest failed
    ------------------------------------------------------------------------
    Total time for all tests: 90.4 seconds
        cpu time: 4.5 seconds
        cumulative wall time: 89.2 seconds

Sometimes an error might be so severe that it causes Sage to segfault
or hang.  In such a situation you have a number of options.  The
doctest framework will print out the output so far, so that at least
you know what test caused the problem (if you want this output to
appear in real time use the ``--verbose`` flag).  To have doctests run
under the control of gdb, use the ``--gdb`` flag::

    [roed@sage sage-6.0]$ sage -t --gdb src/sage/schemes/elliptic_curves/constructor.py
    gdb -x /home/roed/sage-6.0.b5/local/bin/sage-gdb-commands --args python /home/roed/sage-6.0.b5/local/bin/sage-runtests --serial --nthreads 1 --timeout 1048576 --optional sage --stats_path /home/roed/.sage/timings2.json src/sage/schemes/elliptic_curves/constructor.py
    GNU gdb 6.8-debian
    Copyright (C) 2008 Free Software Foundation, Inc.
    License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
    This is free software: you are free to change and redistribute it.
    There is NO WARRANTY, to the extent permitted by law.  Type "show copying"
    and "show warranty" for details.
    This GDB was configured as "x86_64-linux-gnu"...
    [Thread debugging using libthread_db enabled]
    [New Thread 0x7f10f85566e0 (LWP 6534)]
    Running doctests with ID 2012-07-07-00-43-36-b1b735e7.
    Doctesting 1 file.
    sage -t src/sage/schemes/elliptic_curves/constructor.py
        [67 tests, 5.8 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 15.7 seconds
        cpu time: 4.4 seconds
        cumulative wall time: 5.8 seconds

    Program exited normally.
    (gdb) quit


Sage also includes valgrind, and you can run doctests under various
valgrind tools to track down memory issues: the relevant flags are
``--valgrind`` (or ``--memcheck``), ``--massif``, ``--cachegrind`` and
``--omega``.  See http://wiki.sagemath.org/ValgrindingSage for more details.

Once you're done fixing whatever problems where revealed by the
doctests, you can rerun just those files that failed their most recent
test by using the ``--failed`` flag (``-f`` for short)::

    [roed@sage sage-6.0]$ sage -t -fa
    Running doctests with ID 2012-07-07-00-45-35-d8b5a408.
    Doctesting entire Sage library.
    Only doctesting files that failed last test.
    No files to doctest


Miscellaneous Options
---------------------

There are various other options that change the behavior of Sage's
doctesting code.

Show only first failure
^^^^^^^^^^^^^^^^^^^^^^^

The first failure in a file often causes a cascade of others, as
NameErrors arise from variables that weren't defined and tests fail
because old values of variables are used.  To only see the first
failure in each doctest block use the ``--initial`` flag (``-i`` for
short).

Show skipped optional tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To print a summary at the end of each file with the number of optional
tests skipped, use the ``--show-skipped`` flag::

   [roed@sage sage-6.0]$ sage -t --show-skipped src/sage/rings/finite_rings/integer_mod.pyx
   Running doctests with ID 2013-03-14-15-32-05-8136f5e3.
   Doctesting 1 file.
   sage -t sage/rings/finite_rings/integer_mod.pyx
       2 axiom tests not run
       1 cunningham test not run
       2 fricas tests not run
       1 long test not run
       3 magma tests not run
       [440 tests, 4.0 s]
   ----------------------------------------------------------------------
   All tests passed!
   ----------------------------------------------------------------------
   Total time for all tests: 4.3 seconds
       cpu time: 2.4 seconds
       cumulative wall time: 4.0 seconds

Running tests with iterations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes tests fail intermittently.  There are two options that allow
you to run tests repeatedly in an attempt to search for Heisenbugs.
The flag ``--global-iterations`` takes an integer and runs the whole
set of tests that many times serially::

    [roed@sage sage-6.0]$ sage -t --global-iterations 2 src/sage/sandpiles
    Running doctests with ID 2012-07-07-00-59-28-e7048ad9.
    Doctesting 3 files (2 global iterations).
    sage -t src/sage/sandpiles/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/sandpiles/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/sandpiles/sandpile.py
        [711 tests, 14.7 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 17.6 seconds
        cpu time: 13.2 seconds
        cumulative wall time: 14.7 seconds
    sage -t src/sage/sandpiles/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/sandpiles/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/sandpiles/sandpile.py
        [711 tests, 13.8 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 14.3 seconds
        cpu time: 26.4 seconds
        cumulative wall time: 28.5 seconds

You can also iterate in a different order: the ``--file-iterations``
flag runs the tests in each file ``N`` times before proceeding::

    [roed@sage sage-6.0]$ sage -t --file-iterations 2 src/sage/sandpiles
    Running doctests with ID 2012-07-07-01-01-43-8f954206.
    Doctesting 3 files (2 file iterations).
    sage -t src/sage/sandpiles/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/sandpiles/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/sandpiles/sandpile.py
        [1422 tests, 13.3 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 29.6 seconds
        cpu time: 12.7 seconds
        cumulative wall time: 13.3 seconds


Note that the reported results are the average time for all tests in
that file to finish.  If a failure in a file occurs, then the failure
is reported and testing proceeds with the next file.

Using a different timeout
^^^^^^^^^^^^^^^^^^^^^^^^^

On a slow machine the default timeout of 5 minutes may not be enough
for the slowest files.  Use the ``--timeout`` flag (``-T`` for short)
to set it to something else::

    [roed@sage sage-6.0]$ sage -tp 2 --all --timeout 1
    Running doctests with ID 2012-07-07-01-09-37-deb1ab83.
    Doctesting entire Sage library.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 2067 files using 2 threads.
    sage -t src/sage/schemes/elliptic_curves/ell_rational_field.py
        Timed out!
    ...

Using absolute paths
^^^^^^^^^^^^^^^^^^^^

By default filenames are printed using relative paths.  To use
absolute paths instead pass in the ``--abspath`` flag::

    [roed@sage sage-6.0]$ sage -t --abspath src/sage/doctest/control.py
    Running doctests with ID 2012-07-07-01-13-03-a023e212.
    Doctesting 1 file.
    sage -t /home/roed/sage-6.0/src/sage/doctest/control.py
        [133 tests, 4.7 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 7.1 seconds
        cpu time: 0.2 seconds
        cumulative wall time: 4.7 seconds


Testing changed files
^^^^^^^^^^^^^^^^^^^^^

If you are working on some files in the Sage library it can be
convenient to test only the files that have changed.  To do so use the
``--new`` flag, which tests files that have been modified or added
since the last commit::

    [roed@sage sage-6.0]$ sage -t --new
    Running doctests with ID 2012-07-07-01-15-52-645620ee.
    Doctesting files changed since last HG commit.
    Doctesting 1 file.
    sage -t src/sage/doctest/control.py
        [133 tests, 3.7 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 3.8 seconds
        cpu time: 0.1 seconds
        cumulative wall time: 3.7 seconds


Running tests in a random order
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, tests are run in the order in which they appear in the
file.  To run tests in a random order (which can reveal subtle bugs),
use the ``--randorder`` flag and pass in a random seed::

    [roed@sage sage-6.0]$ sage -t --new --randorder 127
    Running doctests with ID 2012-07-07-01-19-06-97c8484e.
    Doctesting files changed since last HG commit.
    Doctesting 1 file.
    sage -t src/sage/doctest/control.py
        [133 tests, 3.6 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 3.7 seconds
        cpu time: 0.2 seconds
        cumulative wall time: 3.6 seconds

Note that even with this option, the tests within a given doctest block are still run in order.

Testing external files
^^^^^^^^^^^^^^^^^^^^^^

When testing a file that's not part of the Sage library, the testing
code loads the globals from that file into the namespace before
running tests.  To model the behavior used on the Sage library instead
(where imports must be explicitly specified), use the ``--force-lib``
flag.

Auxilliary files
^^^^^^^^^^^^^^^^

To specify a logfile (rather than use the default which is created for
``sage -t --all``), use the ``--logfile`` flag::

    [roed@sage sage-6.0]$ sage -t --logfile test1.log src/sage/doctest/control.py
    Running doctests with ID 2012-07-07-01-25-49-e7c0e52d.
    Doctesting 1 file.
    sage -t src/sage/doctest/control.py
        [133 tests, 4.3 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 6.7 seconds
        cpu time: 0.1 seconds
        cumulative wall time: 4.3 seconds
    [roed@sage sage-6.0]$ cat test1.log
    Running doctests with ID 2012-07-07-01-25-49-e7c0e52d.
    Doctesting 1 file.
    sage -t src/sage/doctest/control.py
        [133 tests, 4.3 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 6.7 seconds
        cpu time: 0.1 seconds
        cumulative wall time: 4.3 seconds


To give a json file storing the timings for each file, use the
``--stats_path`` flag.  These statistics are used in sorting files so
that slower tests are run first (and thus multiple processes are
utilized most efficiently)::

    [roed@sage sage-6.0]$ sage -tp 2 --stats-path ~/.sage/timings2.json --all
    Running doctests with ID 2012-07-07-01-28-34-2df4251d.
    Doctesting entire Sage library.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 2067 files using 2 threads.
    ...
