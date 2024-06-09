# sage_setup: distribution = sagemath-ntl
# distutils: libraries = ntl gmp
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: libraries = NTL_LIBRARIES
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
"""
NTL error handler

AUTHOR:

- Jeroen Demeyer (2015-02-15): initial version, see :issue:`17784`

- Jeroen Demeyer (2015-07-09): use standard NTL ``ErrorMsgCallback``,
  see :issue:`18875`
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


from sage.libs.ntl.ntl_tools cimport ErrorMsgCallback
from sage.cpython.string cimport char_to_str


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


cdef void NTL_error_callback(const char* s) except *:
    raise NTLError(char_to_str(s))


def setup_NTL_error_callback():
    """
    Setup the NTL error handler callback.

    EXAMPLES::

        sage: from sage.libs.ntl.error import setup_NTL_error_callback
        sage: setup_NTL_error_callback()
    """
    global ErrorMsgCallback
    ErrorMsgCallback = NTL_error_callback
