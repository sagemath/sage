"""
Testing signal handling.

AUTHORS:

 - Jeroen Demeyer (2010-09-29): initial version (#10030)

 - Jeroen Demeyer (2013-11-04): wrap some tests within nogil (#15352)
"""
#*****************************************************************************
#       Copyright (C) 2010 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import signal

cdef extern from '../tests/c_lib.h':
    void ms_sleep(long ms) nogil
    void signal_after_delay(int signum, long ms) nogil
    void signals_after_delay(int signum, long ms, long interval, int n) nogil

cdef extern from *:
    ctypedef int volatile_int "volatile int"


include 'sage/ext/signals.pxi'
include 'sage/ext/interrupt.pxi'
include 'sage/ext/stdsage.pxi'
from cpython cimport PyErr_SetString


# Default delay in milliseconds before raising signals
cdef long DEFAULT_DELAY = 200


########################################################################
# C helper functions                                                   #
########################################################################
cdef void infinite_loop() nogil:
    while True:
        pass

cdef void infinite_malloc_loop() nogil:
    cdef size_t s = 1
    while True:
        sage_free(sage_malloc(s))
        s *= 2
        if (s > 1000000): s = 1

# Dereference a NULL pointer on purpose. This signals a SIGSEGV on most
# systems, but on older Mac OS X and possibly other systems, this
# signals a SIGBUS instead. In any case, this should give some signal.
cdef void dereference_null_pointer() nogil:
    cdef long* ptr = <long*>(0)
    ptr[0] += 1


########################################################################
# Python helper functions                                              #
########################################################################
class return_exception:
    """
    Decorator class which makes a function *return* an exception which
    is raised, to simplify doctests raising exceptions.

    EXAMPLES::

        sage: from sage.tests.interrupt import return_exception
        sage: @return_exception
        ....: def raise_interrupt():
        ....:     raise KeyboardInterrupt("just testing")
        sage: raise_interrupt()
        KeyboardInterrupt('just testing',)
    """
    def __init__ (self, func):
        self.func = func
    def __call__ (self, *args):
        try:
            return self.func(*args)
        except BaseException as e:
            return e

def interrupt_after_delay(ms_delay = 500):
    """
    Send an interrupt signal (``SIGINT``) to the Sage process
    after a delay of ``ms_delay`` milliseconds.

    INPUT:

    - ``ms_delay`` -- (default: 500) a nonnegative integer indicating
      how many milliseconds to wait before raising the interrupt signal.

    EXAMPLES:

    This function is meant to test interrupt functionality.  We
    demonstrate here how to test that the ``factor`` function can be
    interrupted::

        sage: import sage.tests.interrupt
        sage: try:
        ....:     sage.tests.interrupt.interrupt_after_delay()
        ....:     factor(10^1000 + 3)
        ....: except KeyboardInterrupt:
        ....:     print "Caught KeyboardInterrupt"
        Caught KeyboardInterrupt
    """
    signal_after_delay(SIGINT, ms_delay)


########################################################################
# Test basic macros from c_lib/headers/interrupt.h                     #
# Since these are supposed to work without the GIL, we do all tests    #
# (if possible) within a "with nogil" block                            #
########################################################################
def test_sig_off():
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_off()
    """
    with nogil:
        sig_on()
        sig_off()

@return_exception
def test_sig_on(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_on()
        KeyboardInterrupt()
    """
    with nogil:
        signal_after_delay(SIGINT, delay)
        sig_on()
        infinite_loop()

def test_sig_str(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_str()
        Traceback (most recent call last):
        ...
        RuntimeError: Everything ok!
    """
    with nogil:
        sig_str("Everything ok!")
        signal_after_delay(SIGABRT, delay)
        infinite_loop()

cdef c_test_sig_on_cython():
    sig_on()
    infinite_loop()

@return_exception
def test_sig_on_cython(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_on_cython()
        KeyboardInterrupt()
    """
    signal_after_delay(SIGINT, delay)
    c_test_sig_on_cython()

cdef int c_test_sig_on_cython_except() nogil except 42:
    sig_on()
    infinite_loop()

@return_exception
def test_sig_on_cython_except(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_on_cython_except()
        KeyboardInterrupt()
    """
    with nogil:
        signal_after_delay(SIGINT, delay)
        c_test_sig_on_cython_except()

cdef void c_test_sig_on_cython_except_all() nogil except *:
    sig_on()
    infinite_loop()

@return_exception
def test_sig_on_cython_except_all(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_on_cython_except_all()
        KeyboardInterrupt()
    """
    with nogil:
        signal_after_delay(SIGINT, delay)
        c_test_sig_on_cython_except_all()

@return_exception
def test_sig_check(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_check()
        KeyboardInterrupt()
    """
    signal_after_delay(SIGINT, delay)
    while True:
        with nogil:
            sig_check()

@return_exception
def test_sig_check_inside_sig_on(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_check_inside_sig_on()
        KeyboardInterrupt()
    """
    with nogil:
        signal_after_delay(SIGINT, delay)
        sig_on()
        while True:
            sig_check()


########################################################################
# Test sig_retry() and sig_error()                                     #
########################################################################
def test_sig_retry():
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_retry()
        10
    """
    cdef volatile_int v = 0

    with nogil:
        sig_on()
        if v < 10:
            v = v + 1
            sig_retry()
        sig_off()
    return v

@return_exception
def test_sig_retry_and_signal(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_retry_and_signal()
        KeyboardInterrupt()
    """
    cdef volatile_int v = 0

    with nogil:
        sig_on()
        if v < 10:
            v = v + 1
            sig_retry()
        signal_after_delay(SIGINT, delay)
        infinite_loop()

@return_exception
def test_sig_error():
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_error()
        ValueError('some error',)
    """
    sig_on()
    PyErr_SetString(ValueError, "some error")
    sig_error()


########################################################################
# Test no_except macros                                                #
########################################################################
def test_sig_on_no_except(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_on_no_except()
        42
    """
    if not sig_on_no_except():
        # We should never get here, because this sig_on_no_except()
        # will not catch a signal.
        print "Unexpected zero returned from sig_on_no_except()"
    sig_off()

    signal_after_delay(SIGINT, delay)
    if not sig_on_no_except():
        # We get here when we caught a signal.  An exception
        # has been raised, but Cython doesn't realize it yet.
        try:
            # Make Cython realize that there is an exception.
            # To Cython, it will look like the exception was raised on
            # the following line, so the try/except should work.
            cython_check_exception()
        except KeyboardInterrupt:
            return 42
        return 0 # fail
    infinite_loop()

def test_sig_str_no_except(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_str_no_except()
        Traceback (most recent call last):
        ...
        RuntimeError: Everything ok!
    """
    if not sig_on_no_except():
        # We should never get here, because this sig_on_no_except()
        # will not catch a signal.
        print "Unexpected zero returned from sig_on_no_except()"
    sig_off()

    if not sig_str_no_except("Everything ok!"):
        cython_check_exception()
        return 0 # fail
    signal_after_delay(SIGABRT, delay)
    infinite_loop()

@return_exception
def test_sig_check_no_except(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_check_no_except()
        KeyboardInterrupt()
    """
    with nogil:
        signal_after_delay(SIGINT, delay)
        while True:
            if not sig_check_no_except():
                cython_check_exception()
                break # fail


########################################################################
# Test different signals                                               #
########################################################################
def test_signal_segv(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_signal_segv()
        Traceback (most recent call last):
        ...
        SignalError: Segmentation fault
    """
    with nogil:
        sig_on()
        signal_after_delay(SIGSEGV, delay)
        infinite_loop()

def test_signal_fpe(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_signal_fpe()
        Traceback (most recent call last):
        ...
        FloatingPointError: Floating point exception
    """
    with nogil:
        sig_on()
        signal_after_delay(SIGFPE, delay)
        infinite_loop()

def test_signal_ill(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_signal_ill()
        Traceback (most recent call last):
        ...
        SignalError: Illegal instruction
    """
    with nogil:
        sig_on()
        signal_after_delay(SIGILL, delay)
        infinite_loop()

def test_signal_abrt(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_signal_abrt()
        Traceback (most recent call last):
        ...
        RuntimeError: Aborted
    """
    with nogil:
        sig_on()
        signal_after_delay(SIGABRT, delay)
        infinite_loop()

def test_signal_bus(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_signal_bus()
        Traceback (most recent call last):
        ...
        SignalError: Bus error
    """
    with nogil:
        sig_on()
        signal_after_delay(SIGBUS, delay)
        infinite_loop()

def test_signal_quit(long delay = DEFAULT_DELAY):
    """
    TESTS:

    We run Sage in a subprocess and make it raise a SIGQUIT under
    ``sig_on()``.  This should cause Sage to exit::

        sage: from subprocess import *
        sage: cmd = 'from sage.tests.interrupt import *; test_signal_quit()'
        sage: print Popen(['sage', '-c', cmd], stdout=PIPE, stderr=PIPE).communicate()[1]  # long time
        ---...---
    """
    # The sig_on() shouldn't make a difference for SIGQUIT
    with nogil:
        sig_on()
        signal_after_delay(SIGQUIT, delay)
        infinite_loop()


########################################################################
# Test with "true" errors (not signals raised by hand)                 #
########################################################################
def test_dereference_null_pointer():
    """
    TESTS:

    This test should result in either a Segmentation Fault or a Bus
    Error. ::

        sage: from sage.tests.interrupt import *
        sage: test_dereference_null_pointer()
        Traceback (most recent call last):
        ...
        SignalError: ...
    """
    with nogil:
        sig_on()
        dereference_null_pointer()

def unguarded_dereference_null_pointer():
    """
    TESTS:

    We run Sage in a subprocess and dereference a NULL pointer without
    using ``sig_on()``. This will crash Sage::

        sage: from subprocess import *
        sage: cmd = 'from sage.tests.interrupt import *; unguarded_dereference_null_pointer()'
        sage: print Popen(['sage', '-c', cmd], stdout=PIPE, stderr=PIPE).communicate()[1]  # long time
        ---...---
        Unhandled SIG...
        This probably occurred because a *compiled* component of Sage has a bug
        in it and is not properly wrapped with sig_on(), sig_off().
        Sage will now terminate.
        ------------------------------------------------------------------------
    """
    with nogil:
        dereference_null_pointer()

def test_abort():
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_abort()
        Traceback (most recent call last):
        ...
        RuntimeError: Aborted
    """
    with nogil:
        sig_on()
        abort()

def unguarded_abort():
    """
    TESTS:

    We run Sage in a subprocess and make it call abort()::

        sage: from subprocess import *
        sage: cmd = 'from sage.tests.interrupt import *; unguarded_abort()'
        sage: print Popen(['sage', '-c', cmd], stdout=PIPE, stderr=PIPE).communicate()[1]  # long time
        ---...---
        Unhandled SIGABRT: An abort() occurred in Sage.
        This probably occurred because a *compiled* component of Sage has a bug
        in it and is not properly wrapped with sig_on(), sig_off().
        Sage will now terminate.
        ------------------------------------------------------------------------
    """
    with nogil:
        abort()

def test_bad_str(long delay = DEFAULT_DELAY):
    """
    TESTS:

    We run Sage in a subprocess and induce an error during the signal handler::

        sage: from subprocess import *
        sage: cmd = 'from sage.tests.interrupt import *; test_bad_str()'
        sage: print Popen(['sage', '-c', cmd], stdout=PIPE, stderr=PIPE).communicate()[1]  # long time
        ---...---
        An error occured during signal handling.
        This probably occurred because a *compiled* component of Sage has a bug
        in it and is not properly wrapped with sig_on(), sig_off().
        Sage will now terminate.
        ------------------------------------------------------------------------
    """
    cdef char* s = <char*>(16)
    with nogil:
        sig_str(s)
        signal_after_delay(SIGILL, delay)
        infinite_loop()


########################################################################
# Test various usage scenarios for sig_on()/sig_off()                  #
########################################################################
@return_exception
def test_sig_on_cython_after_delay(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_on_cython_after_delay()
        KeyboardInterrupt()
    """
    with nogil:
        signal_after_delay(SIGINT, delay)
        ms_sleep(delay * 2)  # We get signaled during this sleep
        sig_on()             # The signal should be detected here
        abort()              # This should not be reached

def test_sig_on_inside_try(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_on_inside_try()
    """
    try:
        with nogil:
            sig_on()
            signal_after_delay(SIGABRT, delay)
            infinite_loop()
    except RuntimeError:
        pass

def test_interrupt_bomb(int n = 100, int p = 10):
    """
    Have `p` processes each sending `n` interrupts in very quick
    succession and see what happens :-)

    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_interrupt_bomb()  # long time (1200 + 5*p + 10*n milliseconds)
        Received ... interrupts
    """
    cdef int i

    # Spawn p processes, each sending n signals with an interval of 1 millisecond
    cdef long base_delay = DEFAULT_DELAY + 5*p
    for i in range(p):
        signals_after_delay(SIGINT, base_delay, 1, n)

    # Some time later (after the smoke clears up) send a SIGABRT,
    # which will raise RuntimeError.
    signal_after_delay(SIGABRT, base_delay + 10*n + 1000)
    i = 0
    while True:
        try:
            with nogil:
                sig_on()
                infinite_loop()
        except KeyboardInterrupt:
            i = i + 1
        except RuntimeError:
            break
    print "Received %i/%i interrupts"%(i,n*p)

# Special thanks to Robert Bradshaw for suggesting the try/finally
# construction. -- Jeroen Demeyer
def test_try_finally_signal(long delay = DEFAULT_DELAY):
    """
    Test a try/finally construct for sig_on() and sig_off(), raising
    a signal inside the ``try``.

    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_try_finally_signal()
        Traceback (most recent call last):
        ...
        RuntimeError: Aborted
    """
    sig_on()
    try:
        signal_after_delay(SIGABRT, delay)
        infinite_loop()
    finally:
        sig_off()

def test_try_finally_raise():
    """
    Test a try/finally construct for sig_on() and sig_off(), raising
    a Python exception inside the ``try``.

    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_try_finally_raise()
        Traceback (most recent call last):
        ...
        RuntimeError: Everything ok!
    """
    sig_on()
    try:
        raise RuntimeError, "Everything ok!"
    finally:
        sig_off()

def test_try_finally_return():
    """
    Test a try/finally construct for sig_on() and sig_off(), doing a
    normal ``return`` inside the ``try``.

    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_try_finally_return()
        'Everything ok!'
    """
    sig_on()
    try:
        return "Everything ok!"
    finally:
        sig_off()


########################################################################
# Test sig_block()/sig_unblock()                                       #
########################################################################
def test_sig_block(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_block()
        42
    """
    cdef volatile_int v = 0

    try:
        with nogil:
            sig_on()
            sig_block()
            signal_after_delay(SIGINT, delay)
            ms_sleep(delay * 2)  # We get signaled during this sleep
            v = 42
            sig_unblock()        # Here, the interrupt will be handled
            sig_off()
    except KeyboardInterrupt:
        return v

    # Never reached
    return 1

def test_sig_block_outside_sig_on(long delay = DEFAULT_DELAY):
    """
    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sig_block_outside_sig_on()
        'Success'
    """
    with nogil:
        signal_after_delay(SIGINT, delay)

        # sig_block()/sig_unblock() shouldn't do anything
        # since we're outside of sig_on()
        sig_block()
        ms_sleep(delay * 2)  # We get signaled during this sleep
        sig_unblock()

    try:
        sig_on()  # Interrupt caught here
    except KeyboardInterrupt:
        return "Success"
    abort()   # This should not be reached

def test_signal_during_malloc(long delay = DEFAULT_DELAY):
    """
    Test a signal arriving during a sage_malloc() or sage_free() call.
    Since these are wrapped with sig_block()/sig_unblock(), we should
    safely be able to interrupt them.

    TESTS::

        sage: from sage.tests.interrupt import *
        sage: for i in range(4):  # Several times to reduce chances of false positive
        ...       test_signal_during_malloc()
    """
    try:
        with nogil:
            signal_after_delay(SIGINT, delay)
            sig_on()
            infinite_malloc_loop()
    except KeyboardInterrupt:
        pass


########################################################################
# Benchmarking functions                                               #
########################################################################
def sig_on_bench():
    """
    Call ``sig_on()`` and ``sig_off()`` 1 million times.

    TESTS::

        sage: from sage.tests.interrupt import *
        sage: sig_on_bench()
    """
    cdef int i
    with nogil:
        for i in range(1000000):
            sig_on()
            sig_off()

def sig_check_bench():
    """
    Call ``sig_check()`` 1 million times.

    TESTS::

        sage: from sage.tests.interrupt import *
        sage: sig_check_bench()
    """
    cdef int i
    with nogil:
        for i in range(1000000):
            sig_check()


########################################################################
# Test SIGHUP                                                          #
########################################################################
@return_exception
def test_sighup(long delay = DEFAULT_DELAY):
    """
    Test a basic SIGHUP signal, which would normally exit Sage by
    raising ``SystemExit``.

    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sighup()
        SystemExit()
    """
    with nogil:
        signal_after_delay(SIGHUP, delay)
        while True:
            sig_check()

@return_exception
def test_sigterm_and_sigint(long delay = DEFAULT_DELAY):
    """
    Test a SIGHUP and a SIGINT arriving at essentially the same time.
    The SIGINT should be ignored and we should get a ``SystemExit``.

    TESTS::

        sage: from sage.tests.interrupt import *
        sage: test_sigterm_and_sigint()
        SystemExit()
    """
    with nogil:
        sig_on()
        sig_block()
        signal_after_delay(SIGHUP, delay)
        signal_after_delay(SIGINT, delay)
        # 3 sleeps to ensure both signals arrive
        ms_sleep(delay)
        ms_sleep(delay)
        ms_sleep(delay)
        sig_unblock()
        sig_off()

def test_graceful_exit():
    r"""
    TESTS:

    Start a Sage subprocess, spawn a child PARI/GP process and kill the
    Sage process.  The PARI/GP process should exit by itself. ::

        sage: from subprocess import *
        sage: from signal import *
        sage: P = Popen(['sage-ipython'], stdin=PIPE, stdout=PIPE, stderr=PIPE)  # long time
        sage: P.stdin.write('from sage.tests.interrupt import *\n')  # long time
        sage: P.stdin.write('test_graceful_exit()\n')  # long time

    Now read from the child until we read ``"GO"``.  This ensures that
    the child Sage process has properly started before we terminate it::

        sage: while "GO" not in P.stdout.readline(): pass  # long time
        sage: os.kill(P.pid, SIGHUP)  # long time
        sage: print 'stdout =', P.stdout.read()  # long time
        stdout = ...Exiting spawned PARI/GP interpreter process...
        sage: P.wait()  # long time
        0
    """
    # This code is executed in the subprocess
    import os, sys
    from sage.interfaces.gp import gp

    # Keep PARI/GP busy
    gp(0)  # Ensure PARI/GP is started
    gp._expect.sendline("factor(2^1000-3);")

    # Print something to synchronize with the parent
    print("GO")
    sys.stdout.flush()

    # Wait to be killed...
    sig_on()
    infinite_loop()
