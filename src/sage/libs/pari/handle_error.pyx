"""
Functions for handling PARI errors

AUTHORS:

- Peter Bruin (September 2013): initial version (#9640)

"""
include "sage/ext/stdsage.pxi"
include "sage/ext/interrupt.pxi"
include "decl.pxi"

from cpython cimport PyErr_Occurred


# Global variable to record error string
cdef str pari_error_string

cdef void _pari_init_error_handling():
    """
    Set up our code for handling PARI errors.

    TEST::

        sage: try:
        ....:     p = pari.polcyclo(-1)
        ....: except PariError as e:
        ....:     print e.errtext()
        ....:
          ***   argument must be positive in polcyclo.

    """
    global pari_error_string
    global cb_pari_err_recover
    global cb_pari_handle_exception
    pari_error_string = ""
    cb_pari_err_recover = _pari_err_recover
    cb_pari_handle_exception = _pari_handle_exception

cdef void _pari_check_warning():
    """
    Print pending PARI library warning messages to stderr.

    TEST::

        sage: pari('warning("test")')
          ***   user warning: test

    """
    global pari_error_string
    if pari_error_string != "":
        import sys
        sys.stderr.write(pari_error_string)
        pari_error_string = ""

cdef int _pari_handle_exception(long err) except 0:
    """
    Convert a PARI error into a Sage exception, unless the error was
    a stack overflow, in which case we enlarge the stack.

    This function is a callback from the PARI error handler.

    EXAMPLES::

        sage: pari('error("test")')
        Traceback (most recent call last):
        ...
        RuntimeError: PARI user exception
          ***   at top-level: error("test")
          ***                 ^-------------
          ***   user error: test

        sage: pari(1)/pari(0)
        Traceback (most recent call last):
        ...
        PariError: division by zero

    """
    if err == errpile:
        # PARI is out of memory.  We double the size of the PARI stack
        # and retry the computation.
        from sage.libs.pari.all import pari
        pari.allocatemem(silent=True)
        return 0

    if err == user:
        raise RuntimeError("PARI user exception\n%s" % pari_error_string)
    else:
        from sage.libs.pari.all import PariError
        raise PariError(err, pari_error_string)

cdef void _pari_err_recover(long err):
    """
    Reset the error string and jump back to ``sig_on()``, either to
    retry the code (in case of no error) or to make the already-raised
    exception known to Python.

    TEST:

    Perform a computation that requires doubling the default stack
    several times::

        sage: pari.allocatemem(2^12)
        PARI stack size set to 4096 bytes
        sage: x = pari('2^(2^26)')
        sage: x == 2^(2^26)
        True

    """
    global pari_error_string
    pari_error_string = ""
    if not PyErr_Occurred():
        # No exception was raised => retry the computation starting
        # from sig_on(). This can happen if we successfully enlarged the
        # PARI stack in _pari_handle_exception().
        sig_retry()
    else:
        # An exception was raised.  Jump to the signal-handling code
        # which will cause sig_on() to see the exception.
        sig_error()
