"""
NTL error handler

AUTHOR:

- Jeroen Demeyer (2015-02-15): initial version, see :trac:`17784`
"""

#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


include "sage/ext/interrupt.pxi"
from ntl_tools cimport SetErrorCallbackFunction
from cpython.exc cimport PyErr_SetString

class NTLError(RuntimeError):
    """
    Exceptions from the NTL library.

    EXAMPLES::

        sage: a = ntl.ZZX([0])
        sage: a.quo_rem(a)
        Traceback (most recent call last):
        ...
        NTLError: DivRem: division by zero
    """

cdef void NTL_error_callback(const char* s, void* context):
    PyErr_SetString(NTLError, s)
    sig_error()

def setup_NTL_error_callback():
    """
    Setup the NTL error handler callback.

    EXAMPLES::

        sage: from sage.libs.ntl.error import setup_NTL_error_callback
        sage: setup_NTL_error_callback()
    """
    SetErrorCallbackFunction(NTL_error_callback, NULL)
