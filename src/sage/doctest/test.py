# sage_setup: distribution = sagemath-repl
"""
Test the doctesting framework

Many tests (with expected failures or crashes) are run in a
subprocess, those tests can be found in the ``tests/`` subdirectory.

EXAMPLES::

    sage: import signal
    sage: import subprocess
    sage: import time
    sage: from sage.env import SAGE_SRC
    sage: tests_dir = os.path.join(SAGE_SRC, 'sage', 'doctest', 'tests')
    sage: tests_env = dict(os.environ)

Unset :envvar:`TERM` when running doctests, see :issue:`14370`::

    sage: try:
    ....:     del tests_env['TERM']
    ....: except KeyError:
    ....:     pass
    sage: kwds = {'cwd': tests_dir, 'env':tests_env}

Check that :issue:`2235` has been fixed::

    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:       "--random-seed=0", "--optional=sage", "longtime.rst"], **kwds)
    Running doctests...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=0 longtime.rst
    [0 tests, ...s]
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    ...
    0
    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "-l", "longtime.rst"], **kwds)
    Running doctests...
    Doctesting 1 file.
    sage -t --long --warn-long 0.0 --random-seed=0 longtime.rst
    [1 test, ...s]
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    ...
    0

Check handling of tolerances::

    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "tolerance.rst"], **kwds)
    Running doctests...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=0 tolerance.rst
    **********************************************************************
    File "tolerance.rst", line ..., in sage.doctest.tests.tolerance
    Failed example:
        print(":-(")    # abs tol 0.1
    Expected:
        :-)
    Got:
        :-(
    **********************************************************************
    File "tolerance.rst", line ..., in sage.doctest.tests.tolerance
    Failed example:
        print("1.0 2.0 3.0")  # abs tol 0.1
    Expected:
        4.0 5.0
    Got:
        1.0 2.0 3.0
    **********************************************************************
    File "tolerance.rst", line ..., in sage.doctest.tests.tolerance
    Failed example:
        print("Hello")  # abs tol 0.1
    Expected:
        1.0
    Got:
        Hello
    **********************************************************************
    File "tolerance.rst", line ..., in sage.doctest.tests.tolerance
    Failed example:
        print("1.0")  # abs tol 0.1
    Expected:
        Hello
    Got:
        1.0
    **********************************************************************
    File "tolerance.rst", line ..., in sage.doctest.tests.tolerance
    Failed example:
        print("Hello 1.1")  # abs tol 0.1
    Expected:
        Goodbye 1.0
    Got:
        Hello 1.1
    **********************************************************************
    File "tolerance.rst", line ..., in sage.doctest.tests.tolerance
    Failed example:
        print("Hello 1.0")  # rel tol 1e-6
    Expected:
        Goodbye 0.999999
    Got:
        Hello 1.0
    Tolerance exceeded:
        0.999999 vs 1.0, tolerance 2e-6 > 1e-6
    **********************************************************************
    File "tolerance.rst", line ..., in sage.doctest.tests.tolerance
    Failed example:
        print("Hello 1.0")  # rel tol 1e-6
    Expected:
        Hello ...
    Got:
        Hello 1.0
    Note: combining tolerance (# tol) with ellipsis (...) is not supported
    **********************************************************************
    ...
    1

Test the ``--initial`` option::

    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "-i", "initial.rst"], **kwds)
    Running doctests...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=0 initial.rst
    **********************************************************************
    File "initial.rst", line 4, in sage.doctest.tests.initial
    Failed example:
        a = binomiak(10,5)  # random to test that we still get the exception
    Exception raised:
        Traceback (most recent call last):
        ...
        NameError: name 'binomiak' is not defined
    **********************************************************************
    File "initial.rst", line 14, in sage.doctest.tests.initial
    Failed example:
        binomial(10,5)
    Expected:
        255
    Got:
        252
    **********************************************************************
    ...
    ----------------------------------------------------------------------
    sage -t  --warn-long 0.0 --random-seed=0 initial.rst  # 5 doctests failed
    ----------------------------------------------------------------------
    ...
    1

Test the ``--exitfirst`` option::

    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "--exitfirst", "initial.rst"], **kwds)
    Running doctests...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=0 initial.rst
    **********************************************************************
    File "initial.rst", line 4, in sage.doctest.tests.initial
    Failed example:
        a = binomiak(10,5)  # random to test that we still get the exception
    Exception raised:
        Traceback (most recent call last):
        ...
        NameError: name 'binomiak' is not defined
    **********************************************************************
    ...
    ----------------------------------------------------------------------
    sage -t  --warn-long 0.0 --random-seed=0 initial.rst  # 1 doctest failed
    ----------------------------------------------------------------------
    ...
    1

Test a timeout using the ``SAGE_TIMEOUT`` environment variable.  Also set
``CYSIGNALS_CRASH_NDEBUG`` to help ensure the test times out in a timely
manner (:issue:`26912`)::

    sage: from copy import deepcopy
    sage: kwds2 = deepcopy(kwds)
    sage: kwds2['env'].update({'SAGE_TIMEOUT': '1', 'CYSIGNALS_CRASH_NDEBUG': '1'})
    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "99seconds.rst"], **kwds2)
    Running doctests...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=0 99seconds.rst
        Timed out
    **********************************************************************
    Tests run before process (pid=...) timed out:
    ...
    ----------------------------------------------------------------------
    sage -t --warn-long 0.0 --random-seed=0 99seconds.rst  # Timed out
    ----------------------------------------------------------------------
    ...
    4

Test handling of ``KeyboardInterrupt`` in doctests::

    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "keyboardinterrupt.rst"], **kwds)
    Running doctests...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=0 keyboardinterrupt.rst
    **********************************************************************
    File "keyboardinterrupt.rst", line 11, in sage.doctest.tests.keyboardinterrupt
    Failed example:
        raise KeyboardInterrupt
    Exception raised:
        Traceback (most recent call last):
        ...
        KeyboardInterrupt
    **********************************************************************
    ...
    ----------------------------------------------------------------------
    sage -t --warn-long 0.0 --random-seed=0 keyboardinterrupt.rst  # 1 doctest failed
    ----------------------------------------------------------------------
    ...
    1

Interrupt the doctester::

    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "interrupt.rst"], **kwds)
    Running doctests...
    Doctesting 1 file.
    Killing test interrupt.rst
    ----------------------------------------------------------------------
    Doctests interrupted: 0/1 files tested
    ----------------------------------------------------------------------
    ...
    128

Interrupt the doctester (while parallel testing) when a doctest cannot
be interrupted. We also test that passing a ridiculous number of threads
doesn't hurt::

    sage: F = tmp_filename()
    sage: from copy import deepcopy
    sage: kwds2 = deepcopy(kwds)
    sage: kwds2['env']['DOCTEST_TEST_PID_FILE'] = F  # Doctester will write its PID in this file
    sage: subprocess.call(["sage", "-tp", "1000000", "--timeout=120",  # long time
    ....:      "--die_timeout=10", "--optional=sage",
    ....:      "--warn-long", "0", "99seconds.rst", "interrupt_diehard.rst"], **kwds2)
    Running doctests...
    Doctesting 2 files using 1000000 threads...
    Killing test 99seconds.rst
    Killing test interrupt_diehard.rst
    ----------------------------------------------------------------------
    Doctests interrupted: 0/2 files tested
    ----------------------------------------------------------------------
    ...
    128

Even though the doctester master process has exited, the child process
is still alive, but it should be killed automatically
after the ``die_timeout`` given above (10 seconds)::

    sage: # long time
    sage: pid = int(open(F).read())
    sage: time.sleep(2)
    sage: os.kill(pid, signal.SIGQUIT)  # 2 seconds passed => still alive
    sage: time.sleep(8)
    sage: os.kill(pid, signal.SIGQUIT)  # 10 seconds passed => dead  # random
    Traceback (most recent call last):
    ...
    ProcessLookupError: ...

If the child process is dead and removed, the last output should be as above.
However, the child process interrupted its parent process (see
``'interrupt_diehard.rst'``), and became an orphan process. Depending on the
system, an orphan process may eventually become a zombie process instead of
being removed, and then the last output would just be a blank. Hence the ``#
random`` tag.

Test a doctest failing with ``abort()``::

    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "abort.rst"], **kwds)
    Running doctests...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=0 abort.rst
        Killed due to abort
    **********************************************************************
    Tests run before process (pid=...) failed:
    ...
    ------------------------------------------------------------------------
    Unhandled SIGABRT: An abort() occurred.
    This probably occurred because a *compiled* module has a bug
    in it and is not properly wrapped with sig_on(), sig_off().
    Python will now terminate.
    ------------------------------------------------------------------------
    ...
    ----------------------------------------------------------------------
    sage -t --warn-long 0.0 --random-seed=0 abort.rst  # Killed due to abort
    ----------------------------------------------------------------------
    ...
    16

A different kind of crash::

    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "fail_and_die.rst"], **kwds)
    Running doctests...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=0 fail_and_die.rst
    **********************************************************************
    File "fail_and_die.rst", line 5, in sage.doctest.tests.fail_and_die
    Failed example:
        this_gives_a_NameError
    Exception raised:
        Traceback (most recent call last):
        ...
        NameError: name 'this_gives_a_NameError' is not defined
        Killed due to kill signal
    **********************************************************************
    Tests run before process (pid=...) failed:
    ...
    ----------------------------------------------------------------------
    sage -t --warn-long 0.0 --random-seed=0 fail_and_die.rst  # Killed due to kill signal
    ----------------------------------------------------------------------
    ...
    16

Test that ``sig_on_count`` is checked correctly::

    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "sig_on.rst"], **kwds)
    Running doctests...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=0 sig_on.rst
    **********************************************************************
    File "sig_on.rst", line 6, in sage.doctest.tests.sig_on
    Failed example:
        sig_on_count() # check sig_on/off pairings (virtual doctest)
    Expected:
        0
    Got:
        1
    **********************************************************************
    1 item had failures:
       1 of   5 in sage.doctest.tests.sig_on
        [3 tests, 1 failure, ...]
    ----------------------------------------------------------------------
    sage -t --warn-long 0.0 --random-seed=0 sig_on.rst  # 1 doctest failed
    ----------------------------------------------------------------------
    ...
    1

Test logfiles in serial and parallel mode (see :issue:`19271`)::

    sage: t = tmp_filename()
    sage: subprocess.call(["sage", "-t", "--serial", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "--logfile", t, "simple_failure.rst"],
    ....:      stdout=open(os.devnull, "w"), **kwds)
    1
    sage: print(open(t).read())  # long time
    Running doctests...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=0 simple_failure.rst
    **********************************************************************
    File "simple_failure.rst", line 7, in sage.doctest.tests.simple_failure
    Failed example:
        a * b
    Expected:
        20
    Got:
        15
    **********************************************************************
    1 item had failures:
       1 of   5 in sage.doctest.tests.simple_failure
        [4 tests, 1 failure, ...]
    ----------------------------------------------------------------------
    sage -t --warn-long 0.0 --random-seed=0 simple_failure.rst  # 1 doctest failed
    ----------------------------------------------------------------------
    ...

    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "--logfile", t, "simple_failure.rst"],
    ....:      stdout=open(os.devnull, "w"), **kwds)
    1
    sage: print(open(t).read())  # long time
    Running doctests...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=0 simple_failure.rst
    **********************************************************************
    File "simple_failure.rst", line 7, in sage.doctest.tests.simple_failure
    Failed example:
        a * b
    Expected:
        20
    Got:
        15
    **********************************************************************
    1 item had failures:
       1 of   5 in sage.doctest.tests.simple_failure
        [4 tests, 1 failure, ...]
    ----------------------------------------------------------------------
    sage -t --warn-long 0.0 --random-seed=0 simple_failure.rst  # 1 doctest failed
    ----------------------------------------------------------------------
    ...

Test the ``--debug`` option::

    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "--debug", "simple_failure.rst"],
    ....:      stdin=open(os.devnull), **kwds)
    Running doctests...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=0 simple_failure.rst
    **********************************************************************
    File "simple_failure.rst", line 7, in sage.doctest.tests.simple_failure
    Failed example:
        a * b
    Expected:
        20
    Got:
        15
    **********************************************************************
    Previously executed commands:
        s...: a = 3
        s...: b = 5
        s...: a + b
        8
    sage:
    <BLANKLINE>
    Returning to doctests...
    **********************************************************************
    1 item had failures:
       1 of   5 in sage.doctest.tests.simple_failure
        [4 tests, 1 failure, ...]
    ----------------------------------------------------------------------
    sage -t --warn-long 0.0 --random-seed=0 simple_failure.rst  # 1 doctest failed
    ----------------------------------------------------------------------
    ...
    1

Test running under gdb, without and with a timeout::

    sage: subprocess.call(["sage", "-t",  "--warn-long", "0",   # long time, optional: gdb
    ....:      "--random-seed=0", "--optional=sage", "--gdb", "1second.rst"],
    ....:      stdin=open(os.devnull), **kwds)
    exec gdb ...
    Running doctests...
    Doctesting 1 file...
    sage -t... 1second.rst...
        [2 tests, ... s]
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    ...
    0

gdb might need a long time to start up, so we allow 30 seconds::

    sage: subprocess.call(["sage", "-t",  "--warn-long", "0",   # long time, optional: gdb
    ....:      "--random-seed=0", "--optional=sage", "--gdb", "-T30", "99seconds.rst"],
    ....:      stdin=open(os.devnull), **kwds)
    exec gdb ...
    Running doctests...
        Timed out
    4

Test the ``--show-skipped`` option::

    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "--show-skipped", "show_skipped.rst"], **kwds)
    Running doctests ...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=0 show_skipped.rst
        2 tests not run due to known bugs
        1 gap test not run
        1 long test not run
        1 not tested test not run
        0 tests not run because we ran out of time
        [2 tests, ... s]
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    ...
    0

Optional tests are run correctly::

    sage: subprocess.call(["sage", "-t", "--warn-long", "0", "--long",  # long time
    ....:      "--random-seed=0", "--show-skipped", "--optional=sage,gap", "show_skipped.rst"], **kwds)
    Running doctests ...
    Doctesting 1 file.
    sage -t --long --warn-long 0.0 --random-seed=0 show_skipped.rst
        2 tests not run due to known bugs
        1 not tested test not run
        0 tests not run because we ran out of time
        [4 tests, ... s]
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    ...
    0

    sage: subprocess.call(["sage", "-t", "--warn-long", "0", "--long",  # long time
    ....:      "--random-seed=0", "--show-skipped", "--optional=gAp", "show_skipped.rst"], **kwds)
    Running doctests ...
    Doctesting 1 file.
    sage -t --long --warn-long 0.0 --random-seed=0 show_skipped.rst
        2 tests not run due to known bugs
        1 not tested test not run
        2 sage tests not run
        0 tests not run because we ran out of time
        [2 tests, ... s]
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    ...
    0

Test an invalid value for ``--optional``::

    sage: subprocess.call(["sage", "-t", "--warn-long", "0",
    ....:      "--random-seed=0", "--optional=bad-option", "show_skipped.rst"], **kwds)
    Traceback (most recent call last):
    ...
    ValueError: invalid optional tag 'bad-option'
    1

Test ``atexit`` support in the doctesting framework::

    sage: F = tmp_filename()
    sage: os.path.isfile(F)
    True
    sage: from copy import deepcopy
    sage: kwds2 = deepcopy(kwds)
    sage: kwds2['env']['DOCTEST_DELETE_FILE'] = F
    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "atexit.rst"], **kwds2)
    Running doctests...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=0 atexit.rst
        [3 tests, ... s]
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    ...
    0
    sage: os.path.isfile(F)  # long time
    False
    sage: try:
    ....:     os.unlink(F)
    ....: except OSError:
    ....:     pass

Test that random tests are reproducible::

    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=0", "--optional=sage", "random_seed.rst"], **kwds)
    Running doctests...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=0 random_seed.rst
    **********************************************************************
    File "random_seed.rst", line 3, in sage.doctest.tests.random_seed
    Failed example:
        randint(5, 10)
    Expected:
        9
    Got:
        5
    **********************************************************************
    1 item had failures:
       1 of   2 in sage.doctest.tests.random_seed
        [1 test, 1 failure, ...s]
    ----------------------------------------------------------------------
    sage -t --warn-long 0.0 --random-seed=0 random_seed.rst  # 1 doctest failed
    ----------------------------------------------------------------------
    ...
    1
    sage: subprocess.call(["sage", "-t", "--warn-long", "0",    # long time
    ....:      "--random-seed=1", "--optional=sage", "random_seed.rst"], **kwds)
    Running doctests...
    Doctesting 1 file.
    sage -t --warn-long 0.0 --random-seed=1 random_seed.rst
        [1 test, ...s]
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    ...
    0
"""
