# sage_setup: distribution = sagemath-flint
# ****************************************************************************
#       Copyright (C) 2021 Xavier Caruso <xavier.caruso@normalesup.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.rings.padics.precision_error import PrecisionError


cdef inline int ERROR_ABANDON = 1
cdef inline int ERROR_NOTDEFINED = 1 << 1
cdef inline int ERROR_PRECISION = 1 << 2
cdef inline int ERROR_OVERFLOW = 1 << 3
cdef inline int ERROR_NOTSQUARE = 1 << 4  # maybe we should have something more generic here
cdef inline int ERROR_INTEGRAL = 1 << 5
cdef inline int ERROR_DIVISION = 1 << 6
cdef inline int ERROR_CIRCULAR = 1 << 7

cdef inline int ERROR_UNEXPECTED = 1 << 30


def raise_error(error, permissive=False):
    r"""
    Raise an error according to the given error code.

    INPUT:

    - ``error`` -- integer; the error code

    - ``permissive`` -- boolean (default: ``False``); if ``True``,
      do not raise weak errors (precision, abandon)

    TESTS::

        sage: from sage.rings.padics.padic_relaxed_errors import *
        sage: raise_error(0)
        sage: raise_error(64)
        Traceback (most recent call last):
        ...
        ZeroDivisionError: denominator is not invertible

        sage: raise_error(1)
        Traceback (most recent call last):
        ...
        PrecisionError: computation has been abandoned; try to increase precision

        sage: raise_error(1, permissive=True)
    """
    if error & ERROR_UNEXPECTED:
        raise RuntimeError("error code = %s" % error)
    if error & ERROR_CIRCULAR:
        raise RecursionError("definition looks circular")
    if error & ERROR_DIVISION:
        raise ZeroDivisionError("denominator is not invertible")
    if error & ERROR_INTEGRAL:
        raise ValueError("not in the ring of integers")
    if error & ERROR_NOTSQUARE:
        raise ValueError("not a square")
    if not permissive:
        if error & ERROR_OVERFLOW:
            raise OverflowError
        if error & (ERROR_PRECISION | ERROR_NOTDEFINED):
            raise PrecisionError("not enough precision")
        if error & ERROR_ABANDON:
            raise PrecisionError("computation has been abandoned; try to increase precision")
