r"""
Cython interface to the interrupt handling code

See ``src/sage/tests/interrupt.pyx`` for extensive tests.

AUTHORS:

- Jeroen Demeyer (2010-10-13): initial version

"""
#*****************************************************************************
#       Copyright (C) 2011 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include 'sage/ext/interrupt.pxi'
include 'sage/ext/signals.pxi'
from libc.stdio cimport freopen, stdin


class AlarmInterrupt(KeyboardInterrupt):
    """
    Exception class for :func:`alarm` timeouts.

    EXAMPLES::

        sage: raise AlarmInterrupt
        Traceback (most recent call last):
        ...
        AlarmInterrupt
        sage: from sage.ext.interrupt import do_raise_exception
        sage: import signal
        sage: do_raise_exception(signal.SIGALRM)
        Traceback (most recent call last):
        ...
        AlarmInterrupt
    """
    pass

class SignalError(BaseException):
    """
    Exception class for critical signals such as ``SIGSEGV``. Inherits
    from ``BaseException`` because these normally should not be handled.

    EXAMPLES::

        sage: from sage.ext.interrupt import do_raise_exception
        sage: import signal
        sage: do_raise_exception(signal.SIGSEGV)
        Traceback (most recent call last):
        ...
        SignalError: Segmentation fault
    """
    pass

cdef int sig_raise_exception(int sig, const char* msg) except 0:
    """
    Raise an exception for signal number ``sig`` with message ``msg``
    (or a default message if ``msg`` is ``NULL``).
    """
    if sig == SIGHUP or sig == SIGTERM:
        # Redirect stdin from /dev/null to close interactive sessions
        freopen("/dev/null", "r", stdin)
        # This causes Python to exit
        raise SystemExit
    if sig == SIGINT:
        raise KeyboardInterrupt
    if sig == SIGALRM:
        raise AlarmInterrupt
    if sig == SIGILL:
        if msg == NULL:
            msg = "Illegal instruction"
        raise SignalError(msg)
    if sig == SIGABRT:
        if msg == NULL:
            msg = "Aborted"
        raise RuntimeError(msg)
    if sig == SIGFPE:
        if msg == NULL:
            msg = "Floating point exception"
        raise FloatingPointError(msg)
    if sig == SIGBUS:
        if msg == NULL:
            msg = "Bus error"
        raise SignalError(msg)
    if sig == SIGSEGV:
        if msg == NULL:
            msg = "Segmentation fault"
        raise SignalError(msg)

    raise SystemError("unknown signal number %s"%sig)

def do_raise_exception(sig, msg=None):
    """
    Python version of :func:`sig_raise_exception`, just for doctesting.

    EXAMPLES::

        sage: from sage.ext.interrupt import do_raise_exception
        sage: import signal
        sage: do_raise_exception(signal.SIGFPE)
        Traceback (most recent call last):
        ...
        FloatingPointError: Floating point exception
        sage: do_raise_exception(signal.SIGBUS, "CUSTOM MESSAGE")
        Traceback (most recent call last):
        ...
        SignalError: CUSTOM MESSAGE
        sage: do_raise_exception(0)
        Traceback (most recent call last):
        ...
        SystemError: unknown signal number 0

    For interrupts, the message is ignored, see :trac:`17949`::

        sage: do_raise_exception(signal.SIGINT, "ignored")
        Traceback (most recent call last):
        ...
        KeyboardInterrupt
        sage: do_raise_exception(signal.SIGALRM, "ignored")
        Traceback (most recent call last):
        ...
        AlarmInterrupt
    """
    cdef const char* m
    if msg is None:
        m = NULL
    else:
        m = msg
    sig_raise_exception(sig, m)


def init_interrupts():
    """
    Initialize the Sage interrupt framework.

    This is normally done exactly once during Sage startup when
    importing this module.
    """
    # Set the Python-level interrupt handler. When a SIGINT occurs,
    # this will not be called directly. Instead, a SIGINT is caught by
    # the libcsage (c_lib) interrupt handler. If it happens during pure
    # Python code (not within sig_on()/sig_off()), the handler will set
    # Python's interrupt flag. Python regularly checks this and will
    # call its interrupt handler (which is the one we set now). This
    # handler issues a sig_check() which finally raises the
    # KeyboardInterrupt exception.
    import signal
    signal.signal(signal.SIGINT, sage_python_check_interrupt)

    setup_sage_signal_handler()
    _signals.raise_exception = sig_raise_exception


def sig_on_reset():
    """
    Return the current value of ``_signals.sig_on_count`` and set its
    value to zero. This is used by the doctesting framework.

    EXAMPLES::

        sage: from sage.ext.interrupt import sig_on_reset
        sage: cython('sig_on()'); sig_on_reset()
        1
        sage: sig_on_reset()
        0
    """
    cdef int s = _signals.sig_on_count
    _signals.sig_on_count = 0
    return s


def sage_python_check_interrupt(sig, frame):
    """
    Python-level interrupt handler for interrupts raised in Python
    code. This simply delegates to the interrupt handling code in
    libcsage (c_lib).
    """
    sig_check()


init_interrupts()
