"""
Arithmetic functions using the arb library
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

from sage.libs.flint.types cimport ulong
from sage.libs.flint.fmpq cimport fmpq_t, fmpq_init, fmpq_clear, fmpq_get_mpq
from sage.libs.arb.bernoulli cimport bernoulli_fmpq_ui
from sage.libs.arb.acb_modular cimport acb_modular_hilbert_class_poly
from sage.rings.rational cimport Rational
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer_ring import ZZ


def bernoulli(n):
    """
    Return the ``n``-th Bernoulli number using ``arb``.

    INPUT:

    - ``n`` -- unsigned integer

    OUTPUT: a rational number

    EXAMPLES::

        sage: from sage.libs.arb.arith import bernoulli
        sage: bernoulli(24)
        -236364091/2730
        sage: bernoulli(0)
        1
        sage: bernoulli(1)
        -1/2

    TESTS::

        sage: bernoulli(-1)
        Traceback (most recent call last):
        ...
        OverflowError: can...t convert negative value to mp_limb_t
    """
    cdef ulong i = n
    cdef Rational q = <Rational>Rational.__new__(Rational)
    cdef fmpq_t x
    fmpq_init(x)
    bernoulli_fmpq_ui(x, i)
    fmpq_get_mpq(q.value, x)
    fmpq_clear(x)
    return q


def hilbert_class_polynomial(D):
    """
    Return the Hilbert class polynomial for discriminant ``D`` using ``arb``.

    If ``D`` is not a quadratic discriminant, the result is the zero polynomial.

    INPUT:

    - ``D`` -- negative integer

    OUTPUT: integer polynomial

    EXAMPLES::

        sage: from sage.libs.arb.arith import hilbert_class_polynomial
        sage: hilbert_class_polynomial(-163)
        x + 262537412640768000
        sage: hilbert_class_polynomial(-162)
        0
    """
    cdef long n = D
    cdef Polynomial_integer_dense_flint poly
    poly = PolynomialRing(ZZ, "x", implementation="FLINT")()
    acb_modular_hilbert_class_poly(poly._poly, n)
    return poly
