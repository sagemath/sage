"""
Functions for handling PARI errors

AUTHORS:

- Peter Bruin (September 2013): initial version (:trac:`9640`)

- Jeroen Demeyer (January 2015): use ``cb_pari_err_handle`` (:trac:`14894`)

"""

#*****************************************************************************
#       Copyright (C) 2013 Peter Bruin
#       Copyright (C) 2015 Jeroen Demeyer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from .paridecl cimport *
from .paripriv cimport *
include "cysignals/signals.pxi"

from cpython cimport PyErr_Occurred
from pari_instance cimport pari_instance


# We derive PariError from RuntimeError, for backward compatibility with
# code that catches the latter.
class PariError(RuntimeError):
    """
    Error raised by PARI
    """
    def errnum(self):
        r"""
        Return the PARI error number corresponding to this exception.

        EXAMPLES::

            sage: try:
            ....:     pari('1/0')
            ....: except PariError as err:
            ....:     print(err.errnum())
            31
        """
        return self.args[0]

    def errtext(self):
        """
        Return the message output by PARI when this error occurred.

        EXAMPLE::

            sage: try:
            ....:     pari('pi()')
            ....: except PariError as e:
            ....:     print(e.errtext())
            not a function in function call

        """
        return self.args[1]

    def errdata(self):
        """
        Return the error data (a ``t_ERROR`` gen) corresponding to this
        error.

        EXAMPLES::

            sage: try:
            ....:     pari(Mod(2,6))^-1
            ....: except PariError as e:
            ....:     E = e.errdata()
            sage: E
            error("impossible inverse in Fp_inv: Mod(2, 6).")
            sage: E.component(2)
            Mod(2, 6)
        """
        return self.args[2]

    def __repr__(self):
        r"""
        TESTS::

            sage: PariError(11)
            PariError(11)
        """
        return "PariError(%d)"%self.errnum()

    def __str__(self):
        r"""
        Return a suitable message for displaying this exception.

        This is simply the error text with certain trailing characters
        stripped.

        EXAMPLES::

            sage: try:
            ....:     pari('1/0')
            ....: except PariError as err:
            ....:     print(err)
            _/_: impossible inverse in gdiv: 0

        A syntax error::

            sage: pari('!@#$%^&*()')
            Traceback (most recent call last):
            ...
            PariError: syntax error, unexpected $undefined
        """
        return self.errtext().rstrip(" .:")


cdef void _pari_init_error_handling():
    """
    Set up our code for handling PARI errors.

    TESTS::

        sage: try:
        ....:     p = pari.polcyclo(-1)
        ....: except PariError as e:
        ....:     print(e.errtext())
        domain error in polcyclo: index <= 0

    Warnings still work just like in GP::

        sage: pari('warning("test")')
          ***   user warning: test
    """
    global cb_pari_err_handle
    global cb_pari_err_recover
    cb_pari_err_handle = _pari_err_handle
    cb_pari_err_recover = _pari_err_recover


cdef int _pari_err_handle(GEN E) except 0:
    """
    Convert a PARI error into a Sage exception, unless the error was
    a stack overflow, in which case we enlarge the stack.

    This function is a callback from the PARI error handler.

    EXAMPLES::

        sage: pari('error("test")')
        Traceback (most recent call last):
        ...
        PariError: error: user error: test
        sage: pari(1)/pari(0)
        Traceback (most recent call last):
        ...
        PariError: impossible inverse in gdiv: 0

    """
    cdef long errnum = E[1]

    sig_block()
    cdef char* errstr
    cdef char* s
    try:
        if errnum == e_STACK:
            # Custom error message for PARI stack overflow
            pari_error_string = "the PARI stack overflows (current size: {}; maximum size: {})\n"
            pari_error_string += "You can use pari.allocatemem() to change the stack size and try again"
            pari_error_string = pari_error_string.format(pari_mainstack.size, pari_mainstack.vsize)
        else:
            errstr = pari_err2str(E)
            pari_error_string = errstr.decode('ascii')
            pari_free(errstr)

        s = closure_func_err()
        if s is not NULL:
            pari_error_string = s.decode('ascii') + ": " + pari_error_string

        raise PariError(errnum, pari_error_string, pari_instance.new_gen_noclear(E))
    finally:
        sig_unblock()


cdef void _pari_err_recover(long errnum):
    """
    Reset the error string and jump back to ``sig_on()``, either to
    retry the code (in case of no error) or to make the already-raised
    exception known to Python.

    TEST:

    Perform a computation that requires doubling the default stack
    several times::

        sage: pari.allocatemem(2^12, 2^26)
        PARI stack size set to 4096 bytes, maximum size set to 67108864
        sage: x = pari('2^(2^26)')
        sage: x == 2^(2^26)
        True

    """
    # An exception was raised.  Jump to the signal-handling code
    # which will cause sig_on() to see the exception.
    sig_error()
