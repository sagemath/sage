r"""
Generic numerical approximation function
"""

# ****************************************************************************
#       Copyright (C) 2016 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.parent cimport Parent
from sage.structure.element cimport parent
cdef Parent CDF
from sage.rings.real_mpfr import RealField
from sage.rings.complex_mpfr import ComplexField
from sage.rings.complex_double import CDF


def numerical_approx_generic(x, prec):
    """
    Generic implementation of ``numerical_approx`` using coercion or
    conversion to a real or complex field.

    EXAMPLES::

        sage: from sage.arith.numerical_approx import numerical_approx_generic
        sage: numerical_approx_generic(pi, 20)
        3.1416
        sage: numerical_approx_generic(int(42), 20)
        42.000
        sage: numerical_approx_generic(float(4.2), 20)
        4.2000
    """
    P = parent(x)

    cdef Parent RR = RealField(prec)
    cmap = RR.coerce_map_from(P)
    if cmap is not None:
        return cmap(x)

    cdef Parent CC = ComplexField(prec)
    cmap = CC.coerce_map_from(P)
    if cmap is not None:
        return cmap(x)

    # Coercion didn't work: there are 3 possibilities:
    # (1) There is a coercion possible to a lower precision
    # (2) There is a conversion but no coercion
    # (3) The type doesn't convert at all

    # Figure out input precision to check for case (1)
    try:
        inprec = x.prec()
    except AttributeError:
        if prec > 53 and CDF.has_coerce_map_from(P):
            # If we can coerce to CDF, assume input precision was 53 bits
            inprec = 53
        else:
            # Otherwise, assume precision wasn't the issue
            inprec = prec

    if prec > inprec:
        raise TypeError("cannot approximate to a precision of %s bits, use at most %s bits" % (prec, inprec))

    # The issue is not precision, try conversion instead
    try:
        return RR(x)
    except (TypeError, ValueError):
        return CC(x)
