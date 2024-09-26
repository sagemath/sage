.. nodoctest

.. highlight:: shell-session

.. _chapter-doctesting:

=======================
Running Sage's Doctests
=======================

Doctesting a function ensures that the function performs as claimed by
its documentation. Testing can be performed using one thread or
multiple threads. After compiling a source version of Sage, doctesting
can be run on the whole Sage library, on all modules under a given
directory, or on a specified module only. For the purposes of this
chapter, suppose we have compiled Sage from source and the top
level directory is::

    [jdemeyer@localhost sage]$ pwd
    /home/jdemeyer/sage

See the section :ref:`chapter-testing` for information on Sage's
automated testing process. The general syntax for doctesting is as
follows. To doctest a module in the library of a version of Sage, use
this syntax:

.. CODE-BLOCK:: text

    /path/to/sage_root/sage -t [--long] /path/to/sage_root/path/to/module.py[x]

where ``--long`` is an optional argument (see :ref:`section-options`
for more options). The version of ``sage`` used must match the version
of Sage containing the module we want to doctest. A Sage module can be
either a Python script (with the file extension ".py") or it can be a
Cython script, in which case it has the file extension ".pyx".


Testing a module
================

Say we want to run all tests in the sudoku module
``sage/games/sudoku.py``. In a terminal window, first we ``cd`` to the
top level Sage directory of our local Sage installation. Now  we can
start doctesting as demonstrated in the following terminal session::

    [jdemeyer@localhost sage]$ ./sage -t src/sage/games/sudoku.py
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

    [jdemeyer@localhost sage]$ ./sage -t src/sage/games/sudoku.py
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

    [jdemeyer@localhost sage]$ ./sage -t sage/games/sudoku.py
    Running doctests with ID 2012-07-03-03-40-53-6cc4f29f.
    No files matching sage/games/sudoku.py
    No files to doctest

We can also first ``cd`` to the directory containing the module
``sudoku.py`` and doctest that module as follows::

    [jdemeyer@localhost sage]$ cd src/sage/games/
    [jdemeyer@localhost games]$ ls
    __init__.py  hexad.py       sudoku.py           sudoku_backtrack.pyx
    all.py       quantumino.py  sudoku_backtrack.c
    [jdemeyer@localhost games]$ ../../../../sage -t sudoku.py
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

You can also run the Sage doctester as follows::

   [jdemeyer@localhost sage]$ ./sage -tox -e doctest -- src/sage/games/sudoku.py

See :ref:`chapter-tools` for more information about tox.


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

    [jdemeyer@localhost sage]$ ./sage -t src/sage/games/sudoku.py
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
    [jdemeyer@localhost sage]$ ./sage -t "src/sage/games/sudoku.py"
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

    [jdemeyer@localhost sage]$ sage -t src/sage/games/sudoku.py
    sage -t  "src/sage/games/sudoku.py"
    **********************************************************************
    File "/home/jdemeyer/sage/src/sage/games/sudoku.py", line 515:
        sage: next(h.solve(algorithm='backtrack'))
    Exception raised:
        Traceback (most recent call last):
          File "/usr/local/sage/local/bin/ncadoctest.py", line 1231, in run_one_test
            self.run_one_example(test, example, filename, compileflags)
          File "/usr/local/sage/local/bin/sagedoctest.py", line 38, in run_one_example
            OrigDocTestRunner.run_one_example(self, test, example, filename, compileflags)
          File "/usr/local/sage/local/bin/ncadoctest.py", line 1172, in run_one_example
            compileflags, 1) in test.globs
          File "<doctest __main__.example_13[4]>", line 1, in <module>
            next(h.solve(algorithm='backtrack'))###line 515:
        sage: next(h.solve(algorithm='backtrack'))
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

Parallel testing many modules
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

    [jdemeyer@localhost sage]$ ./sage -t src/sage/crypto
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
    sage -t src/sage/crypto/sbox.pyx
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

    [jdemeyer@localhost sage]$ ./sage -t --long src/sage/crypto/
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
    sage -t --long src/sage/crypto/sbox.pyx
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
as taking a long time:

.. CODE-BLOCK:: python

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

    [jdemeyer@localhost sage]$ ./sage -tp 4 src/sage/crypto/
    Running doctests with ID 2012-07-07-00-11-55-9b17765e.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 24 files using 4 threads.
    sage -t src/sage/crypto/boolean_function.pyx
        [252 tests, 3.8 s]
    sage -t src/sage/crypto/block_cipher/miniaes.py
        [429 tests, 1.1 s]
    sage -t src/sage/crypto/mq/sr.py
        [432 tests, 5.7 s]
    sage -t src/sage/crypto/sbox.pyx
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
    [jdemeyer@localhost sage]$ ./sage -tp 4 --long src/sage/crypto/
    Running doctests with ID 2012-07-07-00-13-04-d71f3cd4.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 24 files using 4 threads.
    sage -t --long src/sage/crypto/boolean_function.pyx
        [252 tests, 3.7 s]
    sage -t --long src/sage/crypto/block_cipher/miniaes.py
        [429 tests, 1.0 s]
    sage -t --long src/sage/crypto/sbox.pyx
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

Parallel testing the whole Sage library
=======================================

The main Sage library resides in the directory
:sage_root:`src/`. We can use the syntax described above
to doctest the main library using multiple threads. When doing release
management or patching the main Sage library, a release manager would
parallel test the library using 10 threads with the following command::

    [jdemeyer@localhost sage]$ ./sage -tp 10 --long src/

Another way is run ``make ptestlong``, which builds Sage (if necessary),
builds the Sage documentation (if necessary), and then runs parallel
doctests.  This determines the number of threads by reading the
environment variable :envvar:`MAKE`: if it is set to ``make -j12``, then
use 12 threads.  If :envvar:`MAKE` is not set, then by default it uses
the number of CPU cores (as determined by the Python function
:func:`multiprocessing.cpu_count`) with a minimum of 2 and a maximum of 8.
(When this runs under the control of the `GNU make jobserver
<https://www.gnu.org/software/make/manual/make.html#Parallel>`_, then Sage
will request as most this number of job slots.)

In any case, this will test the Sage library with multiple threads::

    [jdemeyer@localhost sage]$ make ptestlong

Any of the following commands would also doctest the Sage library or
one of its clones:

.. CODE-BLOCK:: text

    make test
    make check
    make testlong
    make ptest
    make ptestlong

The differences are:

* ``make test`` and ``make check`` --- These two commands run the same
  set of tests. First the Sage standard documentation is tested,
  i.e. the documentation that resides in

  * :sage_root:`src/doc/common`
  * :sage_root:`src/doc/en`
  * :sage_root:`src/doc/fr`

  Finally, the commands doctest the Sage library. For more details on
  these command, see the file :sage_root:`Makefile`.

* ``make testlong`` --- This command doctests the standard
  documentation:

  * :sage_root:`src/doc/common`
  * :sage_root:`src/doc/en`
  * :sage_root:`src/doc/fr`

  and then the Sage library. Doctesting is run with the optional
  argument ``--long``. See the file :sage_root:`Makefile` for further
  details.

* ``make ptest`` --- Similar to the commands ``make test`` and ``make
  check``. However, doctesting is run with the number of threads as
  described above for ``make ptestlong``.

* ``make ptestlong`` --- Similar to the command ``make ptest``, but
  using the optional argument ``--long`` for doctesting.

The underlying command for running these tests is ``sage -t --all``. For
example, ``make ptestlong`` executes the command
``sage -t -p --all --long --logfile=logs/ptestlong.log``. So if you want
to add extra flags when you run these tests, for example ``--verbose``,
you can execute
``sage -t -p --all --long --verbose --logfile=path/to/logfile``.
Some of the extra testing options are discussed here; run
``sage -t -h`` for a complete list.


Beyond the Sage library
=======================

Doctesting also works fine for files not in the Sage library.  For
example, suppose we have a Python script called
``my_python_script.py``::

    [mvngu@localhost sage]$ cat my_python_script.py
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

    [mvngu@localhost sage]$ ./sage -t my_python_script.py
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

    [mvngu@localhost sage]$ cat my_sage_script.sage
    def cube(n):
        r"""
        Return the cube of n.

        EXAMPLES::

            sage: cube(2)
            8
        """
        return n**3

Then we can doctest it just as for Python files::

    [mvngu@localhost sage]$ ./sage -t my_sage_script.sage
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

    [mvngu@localhost sage]$ ./sage --preparse my_sage_script.sage
    [mvngu@localhost sage]$ cat my_sage_script.sage.py
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
    [mvngu@localhost sage]$ ./sage -t my_sage_script.sage.py
    Running doctests with ID 2012-07-07-00-26-46-2bb00911.
    Doctesting 1 file.
    sage -t my_sage_script.sage.py
        [2 tests, 0.0 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 2.3 seconds
        cpu time: 0.0 seconds
        cumulative wall time: 0.0 seconds


Doctesting from within Sage
===========================

You can run doctests from within Sage, which can be useful since you
don't have to wait for Sage to start.  Use the ``run_doctests``
function in the global namespace, passing it either a string or a module:

.. CODE-BLOCK:: ipycon

    sage: run_doctests(sage.combinat.affine_permutation)
    Running doctests with ID 2018-02-07-13-23-13-89fe17b1.
    Git branch: develop
    Using --optional=sagemath_doc_html,sage
    Doctesting 1 file.
    sage -t /opt/sage/sage_stable/src/sage/combinat/affine_permutation.py
        [338 tests, 4.32 s]
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    Total time for all tests: 4.4 seconds
        cpu time: 3.6 seconds
        cumulative wall time: 4.3 seconds


.. _section-options:

Optional arguments
==================

Run long doctests
-----------------

Ideally, doctests should not take any noticeable amount of time. If
you really need longer-running doctests (anything beyond about one
second) then you should mark them as:

.. CODE-BLOCK:: text

    sage: my_long_test()  # long time

Even then, long doctests should ideally complete in 5 seconds or
less. We know that you (the author) want to show off the capabilities
of your code, but this is not the place to do so. Long-running tests
will sooner or later hurt our ability to run the testsuite. Really,
doctests should be as fast as possible while providing coverage for
the code.

Use the ``--long`` flag to run doctests that have been marked with the
comment ``# long time``. These tests are normally skipped in order to
reduce the time spent running tests::

    [roed@localhost sage]$ ./sage -t src/sage/rings/tests.py
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

    [roed@localhost sage]$ ./sage -t --long src/sage/rings/tests.py
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
``--warn-long`` flag.  Without any options it will cause tests to
print a warning if they take longer than 1.0 second. Note that this is
a warning, not an error::

    [roed@localhost sage]$ ./sage -t --warn-long src/sage/rings/factorint.pyx
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
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    Total time for all tests: 16.1 seconds
        cpu time: 9.7 seconds
        cumulative wall time: 10.9 seconds

You can also pass in an explicit amount of time::

    [roed@localhost sage]$ ./sage -t --long --warn-long 2.0 src/sage/rings/tests.py
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
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    Total time for all tests: 27.6 seconds
        cpu time: 24.8 seconds
        cumulative wall time: 26.3 seconds

Finally, you can disable any warnings about long tests with
``--warn-long 0``.

Doctests start from a random seed::

    [kliem@localhost sage]$ ./sage -t src/sage/doctest/tests/random_seed.rst
    Running doctests with ID 2020-06-23-23-22-59-49f37a55.
    ...
    Doctesting 1 file.
    sage -t --warn-long 89.5 --random-seed=112986622569797306072457879734474628454 src/sage/doctest/tests/random_seed.rst
    **********************************************************************
    File "src/sage/doctest/tests/random_seed.rst", line 3, in sage.doctest.tests.random_seed
    Failed example:
        randint(5, 10)
    Expected:
        9
    Got:
        8
    **********************************************************************
    1 item had failures:
       1 of   2 in sage.doctest.tests.random_seed
        [1 test, 1 failure, 0.00 s]
    ----------------------------------------------------------------------
    sage -t --warn-long 89.5 --random-seed=112986622569797306072457879734474628454 src/sage/doctest/tests/random_seed.rst  # 1 doctest failed
    ----------------------------------------------------------------------
    Total time for all tests: 0.0 seconds
        cpu time: 0.0 seconds
        cumulative wall time: 0.0 seconds

This seed can be set explicitly to reproduce possible failures::

    [kliem@localhost sage]$ ./sage -t --warn-long 89.5                              \
                              --random-seed=112986622569797306072457879734474628454 \
                              src/sage/doctest/tests/random_seed.rst
    Running doctests with ID 2020-06-23-23-24-28-14a52269.
    ...
    Doctesting 1 file.
    sage -t --warn-long 89.5 --random-seed=112986622569797306072457879734474628454 src/sage/doctest/tests/random_seed.rst
    **********************************************************************
    File "src/sage/doctest/tests/random_seed.rst", line 3, in sage.doctest.tests.random_seed
    Failed example:
        randint(5, 10)
    Expected:
        9
    Got:
        8
    **********************************************************************
    1 item had failures:
       1 of   2 in sage.doctest.tests.random_seed
        [1 test, 1 failure, 0.00 s]
    ----------------------------------------------------------------------
    sage -t --warn-long 89.5 --random-seed=112986622569797306072457879734474628454 src/sage/doctest/tests/random_seed.rst  # 1 doctest failed
    ----------------------------------------------------------------------
    Total time for all tests: 0.0 seconds
        cpu time: 0.0 seconds
        cumulative wall time: 0.0 seconds

It can also be set explicitly using the environment variable
``SAGE_DOCTEST_RANDOM_SEED``.


.. _section-optional-doctest-flag:

Run optional doctests
---------------------

You can run tests that require optional packages by using the
``--optional`` flag.  Obviously, you need to have installed the
necessary optional packages in order for these tests to succeed.

By default, Sage only runs doctests that are not marked with the ``optional`` tag.  This is equivalent to running ::

    [roed@localhost sage]$ ./sage -t --optional=sagemath_doc_html,sage \
                                  src/sage/rings/real_mpfr.pyx
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

    [roed@localhost sage]$ ./sage -t --optional=sagemath_doc_html,sage,magma \
                                  src/sage/rings/real_mpfr.pyx
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

In order to just run the tests that are marked as requiring magma, omit ``sage`` and ``sagemath_doc_html``::

    [roed@localhost sage]$ ./sage -t --optional=magma src/sage/rings/real_mpfr.pyx
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

If you want Sage to detect external software or other capabilities
(such as magma, latex, internet) automatically and run all of the
relevant tests, then add ``external``::

    [roed@localhost sage]$ ./sage -t --optional=external src/sage/rings/real_mpfr.pyx
    Running doctests with ID 2016-03-16-14-10-21-af2ebb67.
    Using --optional=external
    External software to be detected: cplex,gurobi,internet,latex,macaulay2,magma,maple,mathematica,matlab,octave,scilab
    Doctesting 1 file.
    sage -t --warn-long 28.0 src/sage/rings/real_mpfr.pyx
        [5 tests, 0.04 s]
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    Total time for all tests: 0.5 seconds
        cpu time: 0.0 seconds
        cumulative wall time: 0.0 seconds
    External software detected for doctesting: magma

To run all tests, regardless of whether they are marked optional, pass ``all`` as the ``optional`` tag::

    [roed@localhost sage]$ ./sage -t --optional=all src/sage/rings/real_mpfr.pyx
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


Running doctests in parallel
----------------------------

If you're testing many files, you can get big speedups by using more
than one thread.  To run doctests in parallel use the ``--nthreads``
flag (``-p`` is a shortened version).  Pass in the number of threads
you would like to use (by default Sage just uses 1)::

    [roed@localhost sage]$ ./sage -tp 2 src/sage/doctest/
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


Doctesting all of Sage
----------------------

To doctest the whole Sage library use the ``--all`` flag (``-a`` for
short).  In addition to testing the code in Sage's Python and Cython
files, this command will run the tests defined in Sage's documentation
as well as testing the Sage notebook::

    [roed@localhost sage]$ ./sage -t -a
    Running doctests with ID 2012-06-22-19-10-27-e26fce6d.
    Doctesting entire Sage library.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 2020 files.
    sage -t /Users/roed/sage/src/sage/plot/plot.py
        [304 tests, 69.0 s]
    ...


Debugging tools
---------------

Sometimes doctests fail (that's why we run them after all).  There are
various flags to help when something goes wrong.  If a doctest
produces a Python error, then normally tests continue after reporting
that an error occurred.  If you use the flag ``--debug`` (``-d`` for
short) then you will drop into an interactive Python debugger whenever
a Python exception occurs.  As an example, I modified
:mod:`sage.schemes.elliptic_curves.constructor` to produce an error::

    [roed@localhost sage]$ ./sage -t --debug \
                                 src/sage/schemes/elliptic_curves/constructor.py
    Running doctests with ID 2012-06-23-12-09-04-b6352629.
    Doctesting 1 file.
    **********************************************************************
    File "sage.schemes.elliptic_curves.constructor", line 4, in sage.schemes.elliptic_curves.constructor
    Failed example:
        EllipticCurve([0,0])
    Exception raised:
        Traceback (most recent call last):
          File ".../site-packages/sage/doctest/forker.py", line 573, in _run
            self.execute(example, compiled, test.globs)
          File ".../site-packages/sage/doctest/forker.py", line 835, in execute
            exec compiled in globs
          File "<doctest sage.schemes.elliptic_curves.constructor[0]>", line 1, in <module>
            EllipticCurve([Integer(0),Integer(0)])
          File ".../site-packages/sage/schemes/elliptic_curves/constructor.py", line 346, in EllipticCurve
            return ell_rational_field.EllipticCurve_rational_field(x, y)
          File ".../site-packages/sage/schemes/elliptic_curves/ell_rational_field.py", line 216, in __init__
            EllipticCurve_number_field.__init__(self, Q, ainvs)
          File ".../site-packages/sage/schemes/elliptic_curves/ell_number_field.py", line 159, in __init__
            EllipticCurve_field.__init__(self, [field(x) for x in ainvs])
          File ".../site-packages/sage/schemes/elliptic_curves/ell_generic.py", line 156, in __init__
            "Invariants %s define a singular curve."%ainvs
        ArithmeticError: Invariants [0, 0, 0, 0, 0] define a singular curve.
    > .../site-packages/sage/schemes/elliptic_curves/ell_generic.py(156)__init__()
    -> "Invariants %s define a singular curve."%ainvs
    (Pdb) l
    151                 if len(ainvs) == 2:
    152                     ainvs = [K(0),K(0),K(0)] + ainvs
    153                 self.__ainvs = tuple(ainvs)
    154                 if self.discriminant() == 0:
    155                     raise ArithmeticError(
    156  ->                     "Invariants %s define a singular curve."%ainvs)
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

    [roed@localhost sage]$ ./sage -t --gdb \
                                  src/sage/schemes/elliptic_curves/constructor.py
    exec gdb --eval-commands="run" --args /home/roed/sage/local/var/lib/sage/venv-python3.9/bin/python3 sage-runtests --serial --timeout=0 --stats-path=/home/roed/.sage/timings2.json --optional=pip,sage,sage_spkg src/sage/schemes/elliptic_curves/constructor.py
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

    [roed@localhost sage]$ ./sage -t -fa
    Running doctests with ID 2012-07-07-00-45-35-d8b5a408.
    Doctesting entire Sage library.
    Only doctesting files that failed last test.
    No files to doctest


Miscellaneous options
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

   [roed@localhost sage]$ ./sage -t --show-skipped \
                                 src/sage/rings/finite_rings/integer_mod.pyx
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

    [roed@localhost sage]$ ./sage -t --global-iterations 2 src/sage/sandpiles
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

    [roed@localhost sage]$ ./sage -t --file-iterations 2 src/sage/sandpiles
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

    [roed@localhost sage]$ ./sage -tp 2 --all --timeout 1
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

    [roed@localhost sage]$ ./sage -t --abspath src/sage/doctest/control.py
    Running doctests with ID 2012-07-07-01-13-03-a023e212.
    Doctesting 1 file.
    sage -t /home/roed/sage/src/sage/doctest/control.py
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

    [roed@localhost sage]$ ./sage -t --new
    Running doctests with ID 2012-07-07-01-15-52-645620ee.
    Doctesting files changed since last git commit.
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

    [roed@localhost sage]$ ./sage -t --new --randorder 127
    Running doctests with ID 2012-07-07-01-19-06-97c8484e.
    Doctesting files changed since last git commit.
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

When testing a file which is not part of a package (which is not in a
directory containing an ``__init__.py`` file), the testing
code loads the globals from that file into the namespace before
running tests.  To disable this behaviour (and require imports to be
explicitly specified), use the ``--force-lib`` option.

.. _section-doctest-auxiliary-files:

Auxiliary files
^^^^^^^^^^^^^^^

To specify a logfile (rather than use the default which is created for
``sage -t --all``), use the ``--logfile`` flag::

    [roed@localhost sage]$ ./sage -t --logfile test1.log src/sage/doctest/control.py
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
    [roed@localhost sage]$ cat test1.log
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


To give a json file storing the timings and pass/fail status for each file, use the
``--stats-path`` flag; the default location of this file is ``~/.sage/timings2.json``.
The doctester reads it if it exists, for the purpose of sorting the files
so that slower tests are run first (and thus multiple processes are utilized most
efficiently)::

    [roed@localhost sage]$ ./sage -tp 2 --stats-path ~/.sage/timings2.json --all
    Running doctests with ID 2012-07-07-01-28-34-2df4251d.
    Doctesting entire Sage library.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 2067 files using 2 threads.
    ...

At the end of the doctest run, Sage updates the json file if it exists or creates
a new one.

The recorded pass/fail status of the files can be used for running only those files
that failed their most recent test by using the ``--failed`` flag (``-f`` for short).

Using the option ``--baseline-stats-path known-test-failures.json``,
it is possible to distinguish files with known doctest failures
from new failures. The file ``known-test-failures.json`` should be
prepared in the same format as ``timings2.json``.

Source files marked as failed there will be marked as "[failed in baseline]"
failures in the doctest report; and if there are only baseline failures, no
new failures, then ``sage -t`` will exit with status code 0 (success).


.. _section-doctesting-venv:

Options for testing in virtual environments
-------------------------------------------

The distribution packages of the modularized Sage library can be tested in virtual environments.
Sage has infrastructure to create such virtual environments using ``tox``, which is explained
in detail in :ref:`section-modularized-doctesting`.  Our examples in this section
refer to this setting, but it applies the same to any user-created virtual environments.

The virtual environments, set up in directories such as
``pkgs/sagemath-standard/.tox/sagepython-sagewheels-nopypi-norequirements``
contain installations of built (non-editable) wheels.

To test all modules of Sage that are installed in a virtual environment,
use the option ``--installed`` (instead of ``--all``)::

    [mkoeppe@localhost sage]$ pkgs/sagemath-standard/.tox/sagepython-.../sage -t   \
                                -p4 --installed

This tests against the doctests as they appear in the installed copies of the files
(in ``site-packages/sage/...``).
Note that these installed copies should never be edited, as they can
be overwritten without warning.

When testing a modularized distribution package other than sagemath-standard,
the top-level module :mod:`sage.all` is not available.  Use the option ``--environment``
to select an appropriate top-level module::

    [mkoeppe@localhost sage]$ pkgs/sagemath-categories/.tox/sagepython-.../sage -t \
                                -p4 --environment sage.all__sagemath_categories    \
                                --installed

To test the installed modules against the doctests as they appear in the source
tree (``src/sage/...``)::

    [mkoeppe@localhost sage]$ pkgs/sagemath-categories/.tox/sagepython-.../sage -t \
                                -p4 --environment sage.all__sagemath_categories    \
                                src/sage/structure

Note that testing all doctests as they appear in the source tree does not make sense
because many of the source files may not be installed in the virtual environment.
Use the option ``--if-installed`` to skip the source files of all Python/Cython modules
that are not installed in the virtual environment::

    [mkoeppe@localhost sage]$ pkgs/sagemath-categories/.tox/sagepython-.../sage -t \
                                -p4 --environment sage.all__sagemath_categories    \
                                --if-installed src/sage/schemes

This option can also be combined with ``--all``::

    [mkoeppe@localhost sage]$ pkgs/sagemath-categories/.tox/sagepython-.../sage -t \
                                -p4 --environment sage.all__sagemath_categories    \
                                --if-installed --all


.. _section-fixdoctests:

The doctest fixer
=================

Sage provides a development tool that assists with updating doctests.


Updating doctest outputs
------------------------

By default, ``./sage --fixdoctests`` runs the doctester and replaces the expected outputs
of all examples by the actual outputs from the current version of Sage::

    [mkoeppe@localhost sage]$ ./sage --fixdoctests \
                                --overwrite src/sage/arith/weird.py

For example, when applied to this Python file::

  | r"""
  | ...
  |
  | EXAMPLES::
  |
  |     sage: 2 + 2
  |     5
  |     sage: factor("91")
  |     "7" * "13"
  | ...

the doctest fixer edits the file as follows::

  | r"""
  | ...
  |
  | EXAMPLES::
  |
  |     sage: 2 + 2
  |     4
  |     sage: factor("91")
  |     Traceback (most recent call last):
  |     ...
  |     TypeError: unable to factor '91'
  | ...

As this command edits the source file, it may be a good practice to first use ``git commit``
to save any changes made in the file.

After running the doctest fixer, it is a good idea to use ``git diff`` to check
all edits that the automated tool made.

An alternative to this workflow is to use the option ``--keep-both``. When expected and
actual output of an example differ, it duplicates the example, marking the two copies
``# optional - EXPECTED`` and ``# optional - GOT``. (Thus, when re-running the doctester,
neither of the two copies is run; this makes ``./sage --fixdoctests`` idempotent.)

When exceptions are expected by an example, it is standard practice to abbreviate
the tracebacks using ``...``.  The doctest fixer uses this abbreviation automatically
when formatting the actual output, as shown in the above example.
To disable it so that the details of the exception
can be inspected, use the option ``--full-tracebacks``. This is particularly useful
in combination with ``--keep-both``::

    [mkoeppe@localhost sage]$ ./sage --fixdoctests --keep-both --full-tracebacks \
                                --overwrite src/sage/arith/weird.py

This will give the following result on the above example::

  | r"""
  | ...
  |
  | EXAMPLES::
  |
  |     sage: 2 + 2                                 # optional - EXPECTED
  |     5
  |     sage: 2 + 2                                 # optional - GOT
  |     4
  |     sage: factor("91")                          # optional - EXPECTED
  |     "7" * "13"
  |     sage: factor("91")                          # optional - GOT
  |     Traceback (most recent call last):
  |     ...
  |     File "<doctest...>", line 1, in <module>
  |     factor("91")
  |     File ".../src/sage/arith/misc.py", line 2680, in factor
  |     raise TypeError("unable to factor {!r}".format(n))
  |     TypeError: unable to factor '91'
  | ...
  | """

To make sure that all doctests are updated, you may have to use the option ``--long``::

    [mkoeppe@localhost sage]$ ./sage --fixdoctests --long \
                                --overwrite src/sage/arith/weird.py

If you are not comfortable with allowing this tool to edit your source files, you can use
the option ``--no-overwrite``, which will create a new file with the extension ``.fixed``
instead of overwriting the source file::

    [mkoeppe@localhost sage]$ ./sage --fixdoctests \
                                --no-overwrite src/sage/arith/weird.py


.. _section-fixdoctests-optional-needs:

Managing ``# optional`` and ``# needs`` tags
--------------------------------------------

When a file uses a ``# sage.doctest: optional/needs FEATURE`` directive, the
doctest fixer automatically removes the redundant ``# optional/needs FEATURE``
tags from all ``sage:`` lines. Likewise, when a block-scoped tag
``sage: # optional/needs FEATURE`` is used, then the doctest fixer removes
redundant tags from all doctests in this scope. For example::

  | # sage.doctest: optional - sirocco, needs sage.rings.number_field
  | r"""
  | ...
  |
  | EXAMPLES::
  |
  |     sage: # needs sage.modules sage.rings.number_field
  |     sage: Q5 = QuadraticField(5)
  |     sage: V = Q5^42                                 # needs sage.modules
  |     sage: T = transmogrify(V)           # optional - bliss sirocco

is automatically transformed to::

  | # sage.doctest: optional - sirocco, needs sage.rings.number_field
  | r"""
  | ...
  |
  | EXAMPLES::
  |
  |     sage: # needs sage.modules
  |     sage: Q5 = QuadraticField(5)
  |     sage: V = Q5^42
  |     sage: T = transmogrify(V)               # optional - bliss

The doctest fixer also aligns the ``# optional/needs FEATURE`` tags on
individual doctests at a fixed set of tab stops.

The doctester may issue style warnings when ``# optional/needs`` tags are
repeated on a whole block of doctests, suggesting to use a block-scoped tag
instead. The doctest fixer makes these changes automatically.

There are situations in which the doctester and doctest fixer show too
much restraint and a manual intervention would improve the formatting
of the doctests. In the example below, the doctester does not issue a
style warning because the first doctest line does not carry the ``# needs``
tag::

  | EXAMPLES::
  |
  |     sage: set_verbose(-1)
  |     sage: P.<x,y,z> = ProjectiveSpace(QQbar, 2)     # needs sage.rings.number_field
  |     sage: C = Curve([x^3*y + 2*x^2*y^2 + x*y^3      # needs sage.rings.number_field
  |     ....:             + x^3*z + 7*x^2*y*z
  |     ....:             + 14*x*y^2*z + 9*y^3*z], P)
  |     sage: Q = P([0,0,1])                            # needs sage.rings.number_field
  |     sage: C.tangents(Q)                             # needs sage.rings.number_field
  |     [x + 4.147899035704788?*y,
  |      x + (1.426050482147607? + 0.3689894074818041?*I)*y,
  |      x + (1.426050482147607? - 0.3689894074818041?*I)*y]

To change this example, there are two approaches:

#. Just add the line ``sage: # needs sage.rings.number_field`` at
   the beginning and run the doctest fixer, which will remove the tags on the individual
   doctests that have now become redundant.

#. Insert a blank line after the first doctest line, splitting the block into two.
   Now the ``# needs`` tag is repeated on the whole second block, so running the doctest
   fixer will add a block-scoped tag and remove the individual tags::

     | EXAMPLES::
     |
     |     sage: set_verbose(-1)
     |
     |     sage: # needs sage.rings.number_field
     |     sage: P.<x,y,z> = ProjectiveSpace(QQbar, 2)
     |     sage: C = Curve([x^3*y + 2*x^2*y^2 + x*y^3
     |     ....:             + x^3*z + 7*x^2*y*z
     |     ....:             + 14*x*y^2*z + 9*y^3*z], P)
     |     sage: Q = P([0,0,1])
     |     sage: C.tangents(Q)
     |     [x + 4.147899035704788?*y,
     |      x + (1.426050482147607? + 0.3689894074818041?*I)*y,
     |      x + (1.426050482147607? - 0.3689894074818041?*I)*y]

In places where the doctester issues a doctest dataflow warning
(``Variable ... referenced here was set only in doctest marked '# optional - FEATURE'``),
the doctest fixer automatically adds the missing ``# optional/needs`` tags.

Sometimes code changes can make existing ``# optional/needs FEATURE`` tags unnecessary.
In an installation or virtual environment where ``FEATURE`` is not available,
you can invoke the doctest fixer with the option ``--probe FEATURE``.
Then it will run examples marked ``# optional/needs - FEATURE`` silently, and if the example
turns out to work anyway, the tag is automatically removed.

.. note::

   Probing works best when the doctests within a docstring do not reuse the same variable
   for different values.

To have the doctest fixer take care of the ``# optional/needs`` tags,
but not change the expected results of examples, use the option ``--only-tags``.
This mode is suitable for mostly unattended runs on many files.

With the option ``--verbose``, the doctest fixer shows the doctester's messages
one by one and reports the changes made.

.. warning::

   While the doctest fixer guarantees to preserve any comments that
   appear before ``# optional/needs`` and all parenthesized comments
   of the form ``# optional - FEATURE (EXPLANATION)``, any free-form comments
   that may be mixed with the doctest tags will be lost.

If you don't want to update any doctests, you can use the
option ``--no-test``. In this mode, the doctest fixer does not run
the doctester and only normalizes the style of the ``# optional`` tags.


Use in virtual environments
---------------------------

The doctest fixer can also run tests using the Sage doctester installed in
a virtual environment::

    [mkoeppe@localhost sage]$ ./sage --fixdoctests --overwrite                      \
                                --distribution sagemath-categories                  \
                                src/sage/geometry/schemes/generic/*.py

This command, using ``--distribution``, is equivalent to a command
that uses the more specific options ``--venv`` and ``--environment``::

    [mkoeppe@localhost sage]$ ./sage --fixdoctests --overwrite                      \
                                --venv pkgs/sagemath-categories/.tox/sagepython-... \
                                --environment sage.all__sagemath_categories
                                src/sage/geometry/schemes/generic/*.py

Either way, the options ``--keep-both``, ``--full-tracebacks``, and
``--if-installed`` are implied.

In this mode of operation, when the doctester encounters a global name
that is unknown in its virtual environment (:class:`NameError`),
the doctest fixer will look up the name in its own environment (typically
a full installation of the Sage library) and add a ``# needs ...`` tag
to the doctest.

Likewise, when the doctester runs into a :class:`ModuleNotFoundError`,
the doctest fixer will automatically add a ``# needs ...`` tag.

The switch ``--distribution`` can be repeated; the given distributions
will be tested in sequence.  Using ``--distribution all`` is equivalent
to a preset list of ``--distribution`` switches.  With the switch
``--fixed-point``, the doctest fixer runs the given distributions until
no more changes are made.


Updating baseline files
-----------------------

The modularized distribution packages ``pkgs/sagemath-categories`` and
``pkgs/sagemath-repl`` contain files ``known-test-failures*.json`` for use
with the option ``--baseline-stats-path``, see section
:ref:`section-doctest-auxiliary-files`.

After running the doctesters of the distributions, for example, via
``sage --fixdoctests``, you can use the test results stored in
``timings2.json`` files to update the ``known-test-failures*.json`` files.
This update can be done using the command::

    [mkoeppe@localhost sage]$ ./sage --fixdoctests --no-test                        \
                                --update-known-test-failures --distribution all
