# sage_setup: distribution = sagemath-modules
"""
Helper code for ternary quadratic forms
"""

# ****************************************************************************
#       Copyright (C) 2012 Gustavo Rama
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import gcd, inverse_mod, xgcd
from sage.matrix.constructor import matrix
from sage.misc.prandom import randint
from sage.rings.finite_rings.integer_mod import mod
from sage.rings.integer_ring import ZZ


def red_mfact(a, b):
    """
    Auxiliary function for reduction that finds the reduction factor of integers `a`, `b`.

    INPUT:

    - ``a``, ``b`` -- integers

    OUTPUT: integer

    EXAMPLES::

        sage: from sage.quadratic_forms.ternary import red_mfact
        sage: red_mfact(0, 3)
        0
        sage: red_mfact(-5, 100)
        9
    """

    if a:
        return (-b + abs(a))//(2*a)
    else:
        return 0


def _reduced_ternary_form_eisenstein_with_matrix(a1, a2, a3, a23, a13, a12):
    """
    Find the coefficients of the equivalent unique reduced ternary form according to the conditions
    of Dickson's "Studies in the Theory of Numbers", pp164-171, and the transformation matrix.
    See TernaryQF.is_eisenstein_reduced for the conditions.

    EXAMPLES::

        sage: from sage.quadratic_forms.ternary import _reduced_ternary_form_eisenstein_with_matrix
        sage: Q = TernaryQF([293, 315, 756, 908, 929, 522])
        sage: qr, M = _reduced_ternary_form_eisenstein_with_matrix(293, 315, 756, 908, 929, 522)
        sage: qr
        (1, 2, 2, -1, 0, -1)
        sage: M
        [ -54  137  -38]
        [ -23   58  -16]
        [  47 -119   33]
        sage: Qr = TernaryQF(qr)
        sage: Qr.is_eisenstein_reduced()
        True
        sage: Q(M) == Qr
        True
    """
    # M = identity_matrix(3)
    # M = matrix(ZZ, 3, [m11, m12, m13, m21, m22, m23, m31, m32, m33])
    m11, m12, m13, m21, m22, m23, m31, m32, m33 = 1, 0, 0, 0, 1, 0, 0, 0, 1

    loop = True

    while loop:

        # adjust
        v = a1 + a2 + a23 + a13 + a12
        if (v < 0):
            # M *= matrix(ZZ, 3, [1, 0, 1, 0, 1, 1, 0, 0, 1])
            [m13] = [m11 + m12 + m13]
            [m23] = [m21 + m22 + m23]
            [m33] = [m31 + m32 + m33]
            a3 += v
            a23 += a12 + 2*a2
            a13 += a12 + 2*a1

        # cuadred 12
        m = red_mfact(a1, a12)
        # M *= matrix(ZZ, 3, [1, m, 0, 0, 1, 0, 0, 0, 1])
        [m12] = [m*m11 + m12]
        [m22] = [m*m21 + m22]
        [m32] = [m*m31 + m32]
        t = a1*m
        a12 += t
        a2 += a12*m
        a12 += t
        a23 += a13*m

        # cuadred 23
        m = red_mfact(a2, a23)
        # M *= matrix(ZZ, 3, [1, 0, 0, 0, 1, m, 0, 0, 1])
        [m13] = [m*m12 + m13]
        [m23] = [m*m22 + m23]
        [m33] = [m*m32 + m33]
        t = a2*m
        a23 += t
        a3 += a23*m
        a23 += t
        a13 += a12*m

        # cuadred 13
        m = red_mfact(a1, a13)
        # M *= matrix(ZZ, 3, [1, 0, m, 0, 1, 0, 0, 0, 1])
        [m13] = [m*m11 + m13]
        [m23] = [m*m21 + m23]
        [m33] = [m*m31 + m33]
        t = a1*m
        a13 += t
        a3 += a13*m
        a13 += t
        a23 += a12*m

        # order 12
        if a1 > a2 or (a1 == a2 and abs(a23) > abs(a13)):
            # M *= matrix(ZZ, 3, [0, -1, 0, -1, 0, 0, 0, 0, -1])
            [m11, m12, m13] = [-m12, -m11, -m13]
            [m21, m22, m23] = [-m22, -m21, -m23]
            [m31, m32, m33] = [-m32, -m31, -m33]
            [a1, a2] = [a2, a1]
            [a13, a23] = [a23, a13]

        # order 23
        if a2 > a3 or (a2 == a3 and abs(a13) > abs(a12)):
            # M *= matrix(ZZ, 3, [-1, 0, 0, 0, 0, -1, 0, -1, 0])
            [m11, m12, m13] = [-m11, -m13, -m12]
            [m21, m22, m23] = [-m21, -m23, -m22]
            [m31, m32, m33] = [-m31, -m33, -m32]
            [a2, a3] = [a3, a2]
            [a13, a12] = [a12, a13]

        # order 12
        if a1 > a2 or (a1 == a2 and abs(a23) > abs(a13)):
            # M *= matrix(ZZ, 3, [0, -1, 0, -1, 0, 0, 0, 0, -1])
            [m11, m12, m13] = [-m12, -m11, -m13]
            [m21, m22, m23] = [-m22, -m21, -m23]
            [m31, m32, m33] = [-m32, -m31, -m33]
            [a1, a2] = [a2, a1]
            [a13, a23] = [a23, a13]

        # signs
        if a23*a13*a12 > 0:
            # a23, a13, a12 positive

            if (a23 < 0):
                # M *= diagonal_matrix([-1, 1, 1])
                m11 = -m11
                m21 = -m21
                m31 = -m31
                a23 = -a23
            if (a13 < 0):
                # M *= diagonal_matrix([1, -1, 1])
                m12 = -m12
                m22 = -m22
                m32 = -m32
                a13 = -a13
            if (a12 < 0):
                # M *= diagonal_matrix([1, 1, -1])
                m13 = -m13
                m23 = -m23
                m33 = -m33
                a12 = -a12

        else:
            # a23, a13, a12 nonpositive

            [s1, s2, s3] = [a23 > 0, a13 > 0, a12 > 0]
            if (s1 + s2 + s3) % 2:
                if a23 == 0:
                    s1 = 1
                else:
                    if a13 == 0:
                        s2 = 1
                    else:
                        if (a12 == 0):
                            s3 = 1
            if s1:
                # M *= diagonal_matrix([-1, 1, 1])
                m11 = -m11
                m21 = -m21
                m31 = -m31
                a23 = -a23
            if s2:
                # M *= diagonal_matrix([1, -1, 1])
                m12 = -m12
                m22 = -m22
                m32 = -m32
                a13 = -a13
            if s3:
                # M *= diagonal_matrix([1, 1, -1])
                m13 = -m13
                m23 = -m23
                m33 = -m33
                a12 = -a12

        loop = not (abs(a23) <= a2 and abs(a13) <= a1 and abs(a12) <= a1 and a1 + a2 + a23 + a13 + a12 >= 0)

    # adj 3
    if a1 + a2 + a23 + a13 + a12 == 0 and 2*a1 + 2*a13 + a12 > 0:
        # M *= matrix(ZZ, 3, [-1, 0, 1, 0, -1, 1, 0, 0, 1])
        [m11, m12, m13] = [-m11, -m12, m11 + m12 + m13]
        [m21, m22, m23] = [-m21, -m22, m21 + m22 + m23]
        [m31, m32, m33] = [-m31, -m32, m31 + m32 + m33]
        # a3 += a1+a2+a23+a13+a12 = 0
        a23 = -2*a2 - a23 - a12
        a13 = -2*a1 - a13 - a12

    # adj 5.12
    if a1 == -a12 and a13 != 0:
        # M *= matrix(ZZ, 3, [-1, -1, 0, 0, -1, 0, 0, 0, 1])
        [m11, m12] = [-m11, -m11 - m12]
        [m21, m22] = [-m21, -m21 - m22]
        [m31, m32] = [-m31, -m31 - m32]
        # a2 += a1 + a12 = 0
        a23 = -a23 - a13
        a13 = -a13
        a12 = -a12  # = 2*a1 + a12

    # adj 5.13
    if a1 == -a13 and a12 != 0:
        # M *= matrix(ZZ, 3, [-1, 0, -1, 0, 1, 0, 0, 0, -1])
        [m11, m13] = [-m11, -m11 - m13]
        [m21, m23] = [-m21, -m21 - m23]
        [m31, m33] = [-m31, -m31 - m33]
        # a3 += a1 + a13 = 0
        a23 = -a23 - a12
        a13 = -a13  # = 2*a1 + a13
        a12 = -a12

    # adj 5.23
    if a2 == -a23 and a12 != 0:
        # M *= matrix(ZZ, 3, [1, 0, 0, 0, -1, -1, 0, 0, -1])
        [m12, m13] = [-m12, -m12 - m13]
        [m22, m23] = [-m22, -m22 - m23]
        [m32, m33] = [-m32, -m32 - m33]
        # a3 += a2 + a23 = 0
        a23 = -a23  # = 2*a2 + a23
        a13 = -a13 - a12
        a12 = -a12

    # adj 4.12
    if a1 == a12 and a13 > 2*a23:
        # M *= matrix(ZZ, 3, [-1, -1, 0, 0, 1, 0, 0, 0, -1])
        [m11, m12, m13] = [-m11, -m11 + m12, -m13]
        [m21, m22, m23] = [-m21, -m21 + m22, -m23]
        [m31, m32, m33] = [-m31, -m31 + m32, -m33]
        # a2 += a1 - a12 = 0
        a23 = -a23 + a13
        # a12 = 2*a1 - a12

    # adj 4.13
    if a1 == a13 and a12 > 2*a23:
        # M *= matrix(ZZ, 3, [-1, 0, -1, 0, -1, 0, 0, 0, 1])
        [m11, m12, m13] = [-m11, -m12, -m11 + m13]
        [m21, m22, m23] = [-m21, -m22, -m21 + m23]
        [m31, m32, m33] = [-m31, -m32, -m31 + m33]
        # a3 += a1 - a13 = 0
        a23 = -a23 + a12
        # a13 = 2*a1 - a13

    # adj 4.23
    if a2 == a23 and a12 > 2*a13:
        # M *= matrix(ZZ, 3, [-1, 0, 0, 0, -1, -1, 0, 0, 1])
        [m11, m12, m13] = [-m11, -m12, -m12 + m13]
        [m21, m22, m23] = [-m21, -m22, -m22 + m23]
        [m31, m32, m33] = [-m31, -m32, -m32 + m33]
        # a3 += a2 - a23 = 0
        # a23 = 2*a2 - a23
        a13 = -a13 + a12

    # order 12
    if a1 == a2 and abs(a23) > abs(a13):
        # M *= matrix(ZZ, 3, [0, -1, 0, -1, 0, 0, 0, 0, -1])
        [m11, m12, m13] = [-m12, -m11, -m13]
        [m21, m22, m23] = [-m22, -m21, -m23]
        [m31, m32, m33] = [-m32, -m31, -m33]
        [a1, a2] = [a2, a1]
        [a13, a23] = [a23, a13]

    # order 23
    if a2 == a3 and abs(a13) > abs(a12):
        # M *= matrix(ZZ, 3, [-1, 0, 0, 0, 0, -1, 0, -1, 0])
        [m11, m12, m13] = [-m11, -m13, -m12]
        [m21, m22, m23] = [-m21, -m23, -m22]
        [m31, m32, m33] = [-m31, -m33, -m32]
        [a13, a12] = [a12, a13]

    # order 12
    if a1 == a2 and abs(a23) > abs(a13):
        # M *= matrix(ZZ, 3, [0, -1, 0, -1, 0, 0, 0, 0, -1])
        [m11, m12, m13] = [-m12, -m11, -m13]
        [m21, m22, m23] = [-m22, -m21, -m23]
        [m31, m32, m33] = [-m32, -m31, -m33]
        [a13, a23] = [a23, a13]

    return (a1, a2, a3, a23, a13, a12), \
            matrix(ZZ, 3, (m11, m12, m13, m21, m22, m23, m31, m32, m33))


def _reduced_ternary_form_eisenstein_without_matrix(a1, a2, a3, a23, a13, a12):
    """
    Find the coefficients of the equivalent unique reduced ternary form according to the conditions
    of Dickson's "Studies in the Theory of Numbers", pp164-171.
    See TernaryQF.is_eisenstein_reduced for the conditions.

    EXAMPLES::

        sage: from sage.quadratic_forms.ternary import _reduced_ternary_form_eisenstein_without_matrix
        sage: Q = TernaryQF([293, 315, 756, 908, 929, 522])
        sage: qr = _reduced_ternary_form_eisenstein_without_matrix(293, 315, 756, 908, 929, 522)
        sage: qr
        (1, 2, 2, -1, 0, -1)
        sage: Qr = TernaryQF(qr)
        sage: Qr.is_eisenstein_reduced()
        True
    """
    loop = True

    while loop:

        # adjust
        v = a1+a2+a23+a13+a12
        if v < 0:
            a3 += v
            a23 += a12+2*a2
            a13 += a12+2*a1

        # cuadred 12
        m = red_mfact(a1, a12)
        t = a1*m
        a12 += t
        a2 += a12*m
        a12 += t
        a23 += a13*m

        # cuadred 23
        m = red_mfact(a2, a23)
        t = a2*m
        a23 += t
        a3 += a23*m
        a23 += t
        a13 += a12*m

        # cuadred 13
        m = red_mfact(a1, a13)
        t = a1*m
        a13 += t
        a3 += a13*m
        a13 += t
        a23 += a12*m

        # order 12
        if a1 > a2 or (a1 == a2 and abs(a23) > abs(a13)):
            [a1, a2] = [a2, a1]
            [a13, a23] = [a23, a13]

        # order 23
        if a2 > a3 or (a2 == a3 and abs(a13) > abs(a12)):
            [a2, a3] = [a3, a2]
            [a13, a12] = [a12, a13]

        # order 12
        if a1 > a2 or (a1 == a2 and abs(a23) > abs(a13)):
            [a1, a2] = [a2, a1]
            [a13, a23] = [a23, a13]

        # signs
        if a23*a13*a12 > 0:
            # a23, a13, a12 positive

            if a23 < 0:
                a23 = -a23
            if a13 < 0:
                a13 = -a13
            if a12 < 0:
                a12 = -a12

        else:
            # a23, a13, a12 nonpositive

            [s1, s2, s3] = [a23 > 0, a13 > 0, a12 > 0]
            if (s1+s2+s3) % 2:
                if a23 == 0:
                    s1 = 1
                elif a13 == 0:
                    s2 = 1
                elif a12 == 0:
                    s3 = 1
            if s1:
                a23 = -a23
            if s2:
                a13 = -a13
            if s3:
                a12 = -a12

        loop = not (abs(a23) <= a2 and abs(a13) <= a1 and abs(a12) <= a1 and a1+a2+a23+a13+a12 >= 0)

    # adj 3
    if a1+a2+a23+a13+a12 == 0 and 2*a1+2*a13+a12 > 0:
        # a3 += a1+a2+a23+a13+a12
        a23 = -2*a2-a23-a12
        a13 = -2*a1-a13-a12

    # adj 5.12
    if a1 == -a12 and a13 != 0:
        # a2 += a1+a12
        a23 = -a23-a13
        a13 = -a13
        a12 = -a12  # = 2*a1+a12

    # adj 5.13
    if a1 == -a13 and a12 != 0:
        # a3 += a1+a13
        a23 = -a23-a12
        a13 = -a13  # = 2*a1+a13
        a12 = -a12

    # adj 5.23
    if a2 == -a23 and a12 != 0:
        # a3 += a2+a23
        a23 = -a23  # = 2*a2+a23
        a13 = -a13-a12
        a12 = -a12

    # adj 4.12
    if a1 == a12 and a13 > 2*a23:
        # a 2 += a1-a12
        a23 = -a23 + a13
        # a12 = 2*a1 - a12

    # adj 4.13
    if a1 == a13 and a12 > 2*a23:
        # a3 += a1-a13
        a23 = -a23 + a12
        # a13 = 2*a1 - a13

    # adj 4.23
    if a2 == a23 and a12 > 2*a13:
        # a3 += a2-a23
        # a23 = 2*a2 - a23
        a13 = -a13 + a12

    # order 12
    if a1 == a2 and abs(a23) > abs(a13):
        [a1, a2] = [a2, a1]
        [a13, a23] = [a23, a13]

    # order 23
    if a2 == a3 and abs(a13) > abs(a12):
        [a13, a12] = [a12, a13]

    # order 12
    if a1 == a2 and abs(a23) > abs(a13):
        [a13, a23] = [a23, a13]

    return a1, a2, a3, a23, a13, a12


def primitivize(long long v0, long long v1, long long v2, p):
    """
    Given a 3-tuple `v` not singular mod `p`, return a primitive 3-tuple version of `v` mod `p`.

    EXAMPLES::

        sage: from sage.quadratic_forms.ternary import primitivize
        sage: primitivize(12, 13, 14, 5)
        (3, 2, 1)
        sage: primitivize(12, 13, 15, 5)
        (4, 1, 0)
    """
    if v2 % p:
        v2_inv = inverse_mod(v2, p)
        return v2_inv*v0 % p, v2_inv*v1 % p, 1
    elif v1 % p:
        return inverse_mod(v1, p)*v0 % p, 1, 0
    else:
        return 1, 0, 0


def evaluate(a, b, c, r, s, t, v):
    """
    Function to evaluate the ternary quadratic form `(a, b, c, r, s, t)` in a 3-tuple `v`.

    EXAMPLES::

        sage: from sage.quadratic_forms.ternary import evaluate
        sage: Q = TernaryQF([1, 2, 3, -1, 0, 0])
        sage: v = (1, -1, 19)
        sage: Q(v)
        1105
        sage: evaluate(1, 2, 3, -1, 0, 0, v)
        1105
    """
    return a*v[0]**2+b*v[1]**2+c*v[2]**2+r*v[2]*v[1]+s*v[2]*v[0]+t*v[1]*v[0]


def _find_zeros_mod_p_2(a, b, c, r, s, t):
    """
    Function to find the zeros mod 2 of a ternary quadratic form.

    EXAMPLES::

        sage: Q = TernaryQF([1, 2, 2, -1, 0, 0])
        sage: from sage.quadratic_forms.ternary import _find_zeros_mod_p_2
        sage: zeros = _find_zeros_mod_p_2(1, 2, 2, -1, 0, 0)
        sage: zeros
        [(0, 1, 0), (0, 0, 1), (1, 1, 1)]
        sage: Q((0, 1, 0))
        2
        sage: Q((0, 0, 1))
        2
        sage: Q((1, 1, 1))
        4
    """
    zeros = []
    v = (1, 0, 0)
    if evaluate(a, b, c, r, s, t, v) % 2 == 0:
        zeros.append(v)
    for i in range(2):
        v = (i, 1, 0)
        if evaluate(a, b, c, r, s, t, v) % 2 == 0:
            zeros.append(v)
    for i in range(2):
        for j in range(2):
            v = (i, j, 1)
            if evaluate(a, b, c, r, s, t, v) % 2 == 0:
                zeros.append(v)
    return zeros


def pseudorandom_primitive_zero_mod_p(a, b, c, r, s, t, p):
    """
    Find a zero of the form `(a, b, 1)` of the ternary quadratic form given by the coefficients `(a, b, c, r, s, t)`
    mod `p`, where `p` is a odd prime that doesn't divide the discriminant.

    EXAMPLES::

        sage: Q = TernaryQF([1, 2, 2, -1, 0, 0])
        sage: p = 1009
        sage: from sage.quadratic_forms.ternary import pseudorandom_primitive_zero_mod_p
        sage: v = pseudorandom_primitive_zero_mod_p(1, 2, 2, -1, 0, 0, p)               # needs sage.libs.pari
        sage: v[2]                                                                      # needs sage.libs.pari
        1
        sage: Q(v)%p                                                                    # needs sage.libs.pari
        0
    """
    # [a, b, c, r, s, t] = Q.coefficients()
    while True:

        r1 = randint(0, p-1)
        r2 = randint(0, p-1)
        alpha = (b*r1**2+t*r1+a) % p
        if alpha != 0:

            beta = (2*b*r1*r2+t*r2+r*r1+s) % p
            gamma = (b*r2**2+r*r2+c) % p
            disc = beta**2-4*alpha*gamma
            if mod(disc, p).is_square():

                z = (-beta+mod(disc, p).sqrt().lift())*(2*alpha).inverse_mod(p)
                return z % p, (r1*z+r2) % p, 1


def _find_zeros_mod_p_odd(long long a, long long b, long long c, long long r, long long s, long long t, long long p, v):
    """
    Find the zeros mod `p`, where `p` is an odd prime, of a ternary quadratic form given by its coefficients and a given zero of the form `v`.

    The prime `p` does not divide the discriminant of the form.

    EXAMPLES::

        sage: from sage.quadratic_forms.ternary import _find_zeros_mod_p_odd
        sage: Q = TernaryQF([1, 2, 2, -1, 0, 0])
        sage: p = 1009
        sage: v = (817, 974, 1)
        sage: Q(v)%1009
        0
        sage: zeros_1009 = _find_zeros_mod_p_odd(1, 2, 2, -1, 0, 0, 1009, v)
        sage: len(zeros_1009)
        1010
        sage: zeros_1009.sort()
        sage: zeros_1009[0]
        (0, 32, 1)
        sage: Q((0, 32, 1))
        2018
    """

    cdef long long a_i
    cdef long long c_i
    cdef long long a_inf
    cdef long long c_inf
    cdef long long i
    cdef long long l
    cdef long long x0
    cdef long long y0

    zeros = [v]
    x0 = v[0]
    y0 = v[1]
    more = False
    for i in range(p):
        a_i = (a*x0**2+b*i**2-2*b*i*y0+b*y0**2-t*i*x0+t*x0*y0-2*a*x0+t*i-t*y0+s*x0-r*i+r*y0+a+c-s) % p
        # b_i=(-2*b*i**2+2*b*i*y0+t*i*x0+2*a*x0-2*t*i+t*y0+r*i-2*a+s)%p
        if a_i == 0:
            w = ((x0-1) % p, (y0-i) % p, 1)
            if w == v:
                more = True
            else:
                zeros.append(w)
        else:
            c_i = (b*i**2+t*i+a) % p
            l = c_i*ZZ(a_i).inverse_mod(p)
            w = primitivize(l*(x0-1)+1, l*(y0-i)+i, l, p)
            if w[0] == v[0] and w[1] == v[1] and w[2] == v[2]:
                more = True
            else:
                zeros.append(w)
    if more:
        a_inf = (a*x0**2+b*y0**2+t*x0*y0-2*b*y0-t*x0+s*x0+r*y0+b+c-r) % p
        # b_inf=(2*b*y0+t*x0-2*b+r)%p
        if a_inf == 0:
            w = (x0 % p, (y0-1) % p, 1)
            zeros.append(w)
        else:
            c_inf = b % p
            l = c_inf*ZZ(a_inf).inverse_mod(p)
            w = primitivize(l*x0, l*(y0-1)+1, l, p)
            zeros.append(w)

    return zeros


def _find_zeros_mod_p(a, b, c, r, s, t, p):
    """
    Find the zeros mod `p` of the ternary quadratic form.

    The quadratic form is given by the coefficients `(a, b, c, r, s, t)`,
    and `p` is a prime that does not divide the discriminant of the form.

    EXAMPLES::

        sage: # needs sage.libs.pari
        sage: from sage.quadratic_forms.ternary import _find_zeros_mod_p
        sage: Q = TernaryQF([1, 2, 2, -1, 0, 0])
        sage: p = 1009
        sage: zeros_1009 = _find_zeros_mod_p(1, 2, 2, -1, 0, 0, p)
        sage: len(zeros_1009)
        1010
        sage: zeros_1009.sort()
        sage: zeros_1009[0]
        (0, 32, 1)
        sage: Q((0, 32, 1))
        2018
        sage: zeros_2 = _find_zeros_mod_p(1, 2, 2, -1, 0, 0, 2)
        sage: zeros_2
        [(0, 1, 0), (0, 0, 1), (1, 1, 1)]
    """
    if p == 2:
        return _find_zeros_mod_p_2(a, b, c, r, s, t)
    v = pseudorandom_primitive_zero_mod_p(a, b, c, r, s, t, p)
    return _find_zeros_mod_p_odd(a, b, c, r, s, t, p, v)


def _find_all_ternary_qf_by_level_disc(long long N, long long d):
    """
    Find the coefficients of all the reduced ternary quadratic forms given its discriminant d and level N.
    If N|4d and d|N^2, then it may be some forms with that discriminant and level.

    EXAMPLES::

        sage: from sage.quadratic_forms.ternary import _find_all_ternary_qf_by_level_disc
        sage: _find_all_ternary_qf_by_level_disc(44, 11)
        [(1, 1, 3, 0, -1, 0), (1, 1, 4, 1, 1, 1)]
        sage: _find_all_ternary_qf_by_level_disc(44, 11^2 * 16)
        [(3, 15, 15, -14, -2, -2), (4, 11, 12, 0, -4, 0)]
        sage: Q = TernaryQF([1, 1, 3, 0, -1, 0])
        sage: Q.is_eisenstein_reduced()
        True
        sage: Q.reciprocal_reduced()
        Ternary quadratic form with integer coefficients:
        [4 11 12]
        [0 -4 0]
        sage: _find_all_ternary_qf_by_level_disc(44, 22)
        []
        sage: _find_all_ternary_qf_by_level_disc(44, 33)
        Traceback (most recent call last):
        ...
        ValueError: There are no ternary forms of this level and discriminant
    """
    cdef long long a
    cdef long long b
    cdef long long c
    cdef long long r
    cdef long long s
    cdef long long t
    cdef long long m
    cdef long long g
    cdef long long u
    cdef long long v
    cdef long long g1
    cdef long long m_q
    cdef double a_max
    cdef double r_max
    cdef double b_max

    l = []

    if (4*d) % N:
        raise ValueError("There are no ternary forms of this level and discriminant")
    else:
        m = 4*d//N

    if (N**2) % d:
        raise ValueError("There are no ternary forms of this level and discriminant")

    a = 1
    a_max = (d/2.)**(1/3.)
    while a <= a_max:

        [g, u, v] = xgcd(4*a, m)
        g1 = (ZZ(g).squarefree_part()*g).sqrtrem()[0]
        t = 0
        while t <= a:

            alpha = max(a, m/4/a)
            b = alpha+((((u*t*t//g)-alpha)) % (m//g))
            b_max = (d/2.0/a)**(1/2.)
            while b <= b_max:

                s = 0
                beta = g//gcd(g, 2*t)
                while s <= a:

                    r = -b+((((2*s*t*u//g)+b)) % (m//g))
                    if s*t == 0:
                        r_max = 0
                    else:
                        r_max = b
                    while r <= r_max:

                        if (d-r*s*t+a*r*r+b*s*s) % (4*a*b-t**2) == 0:

                            c = (d-r*s*t+a*r*r+b*s*s)//(4*a*b-t**2)
                            if r <= 0:
                                is_reduced = True
                                if r < -b:
                                    is_reduced = False
                                elif not (b <= c and 0 <= a+b+r-s-t):
                                    is_reduced = False
                                elif a == b and abs(r) > abs(s):
                                    is_reduced = False
                                elif b == c and abs(s) > abs(t):
                                    is_reduced = False
                                elif a+b+r-s-t == 0 and 2*a-2*s-t > 0:
                                    is_reduced = False
                                elif a == t and s != 0:
                                    is_reduced = False
                                elif a == s and t != 0:
                                    is_reduced = False
                                elif b == -r and t != 0:
                                    is_reduced = False
                                if is_reduced:
                                    m_q = gcd((4*b*c-r**2, 4*a*c-s**2, 4*a*b-t**2, 2*s*t-4*a*r, -2*r*t+4*b*s, -2*r*s+4*c*t))
                                    if m_q == m:
                                        l.append((ZZ(a), ZZ(b), ZZ(c), ZZ(r), ZZ(-s), ZZ(-t)))
                            else:
                                is_reduced = True
                                if not (b <= c and 0 <= a+b+r+s+t):
                                    is_reduced = False
                                elif a == b and abs(r) > abs(s):
                                    is_reduced = False
                                elif b == c and abs(s) > abs(t):
                                    is_reduced = False
                                elif a+b+r+s+t == 0 and 2*a+2*s+t > 0:
                                    is_reduced = False
                                elif a == t and s > 2*r:
                                    is_reduced = False
                                elif a == s and t > 2*r:
                                    is_reduced = False
                                elif b == r and t > 2*s:
                                    is_reduced = False
                                if is_reduced:
                                    m_q = gcd((4*b*c-r**2, 4*a*c-s**2, 4*a*b-t**2, 2*s*t-4*a*r, 2*r*t-4*b*s, 2*r*s-4*c*t))
                                    if m_q == m:
                                        l.append((ZZ(a), ZZ(b), ZZ(c), ZZ(r), ZZ(s), ZZ(t)))
                        r += (m//g)
                    s += beta
                b += m//g
            t += g1

        a += 1
    return l


def _find_a_ternary_qf_by_level_disc(long long N, long long d):
    """
    Find the coefficients of a reduced ternary quadratic form given its discriminant d and level N.
    If N|4d and d|N^2, then it may be a form with that discriminant and level.

    EXAMPLES::

        sage: from sage.quadratic_forms.ternary import _find_a_ternary_qf_by_level_disc
        sage: _find_a_ternary_qf_by_level_disc(44, 11)
        (1, 1, 3, 0, -1, 0)
        sage: _find_a_ternary_qf_by_level_disc(44, 11^2 * 16)
        (3, 15, 15, -14, -2, -2)
        sage: Q = TernaryQF([1, 1, 3, 0, -1, 0])
        sage: Q.is_eisenstein_reduced()
        True
        sage: Q.level()
        44
        sage: Q.disc()
        11
        sage: _find_a_ternary_qf_by_level_disc(44, 22)
        sage: _find_a_ternary_qf_by_level_disc(44, 33)
        Traceback (most recent call last):
        ...
        ValueError: There are no ternary forms of this level and discriminant
    """
    cdef long long a
    cdef long long b
    cdef long long c
    cdef long long r
    cdef long long s
    cdef long long t
    cdef long long m
    cdef long long g
    cdef long long u
    cdef long long v
    cdef long long g1
    cdef long long m_q
    cdef double a_max
    cdef double r_max
    cdef double b_max

    if (4*d) % N:
        raise ValueError("There are no ternary forms of this level and discriminant")
    else:
        m = 4*d//N

    if (N**2) % d:
        raise ValueError("There are no ternary forms of this level and discriminant")

    a = 1
    a_max = (d/2.)**(1/3.)
    while a <= a_max:

        [g, u, v] = xgcd(4*a, m)
        g1 = (ZZ(g).squarefree_part()*g).sqrtrem()[0]
        t = 0
        while t <= a:

            alpha = max(a, m/4/a)
            b = alpha+((((u*t*t//g)-alpha)) % (m//g))
            b_max = (d/2.0/a)**(1/2.)
            while b <= b_max:

                s = 0
                beta = g//gcd(g, 2*t)
                while s <= a:

                    r = -b+((((2*s*t*u//g)+b)) % (m//g))
                    if s*t == 0:
                        r_max = 0
                    else:
                        r_max = b
                    while r <= r_max:

                        if (d-r*s*t+a*r*r+b*s*s) % (4*a*b-t**2) == 0:

                            c = (d-r*s*t+a*r*r+b*s*s)//(4*a*b-t**2)
                            if r <= 0:
                                is_reduced = True
                                if r < -b:
                                    is_reduced = False
                                elif not (b <= c and 0 <= a+b+r-s-t):
                                    is_reduced = False
                                elif a == b and abs(r) > abs(s):
                                    is_reduced = False
                                elif b == c and abs(s) > abs(t):
                                    is_reduced = False
                                elif a+b+r-s-t == 0 and 2*a-2*s-t > 0:
                                    is_reduced = False
                                elif a == t and s != 0:
                                    is_reduced = False
                                elif a == s and t != 0:
                                    is_reduced = False
                                elif b == -r and t != 0:
                                    is_reduced = False
                                if is_reduced:
                                    m_q = gcd((4*b*c-r**2, 4*a*c-s**2, 4*a*b-t**2, 2*s*t-4*a*r, -2*r*t+4*b*s, -2*r*s+4*c*t))
                                    if m_q == m:
                                        return ZZ(a), ZZ(b), ZZ(c), ZZ(r), ZZ(-s), ZZ(-t)
                            else:
                                is_reduced = True
                                if not (b <= c and 0 <= a+b+r+s+t):
                                    is_reduced = False
                                elif a == b and abs(r) > abs(s):
                                    is_reduced = False
                                elif b == c and abs(s) > abs(t):
                                    is_reduced = False
                                elif a+b+r+s+t == 0 and 2*a+2*s+t > 0:
                                    is_reduced = False
                                elif a == t and s > 2*r:
                                    is_reduced = False
                                elif a == s and t > 2*r:
                                    is_reduced = False
                                elif b == r and t > 2*s:
                                    is_reduced = False
                                if is_reduced:
                                    m_q = gcd((4*b*c-r**2, 4*a*c-s**2, 4*a*b-t**2, 2*s*t-4*a*r, 2*r*t-4*b*s, 2*r*s-4*c*t))
                                    if m_q == m:
                                        return ZZ(a), ZZ(b), ZZ(c), ZZ(r), ZZ(s), ZZ(t)
                        r += (m//g)
                    s += beta
                b += m//g
            t += g1

        a += 1


def extend(v):
    """
    Return the coefficients of a matrix `M` such that `M` has determinant gcd(`v`) and the first column is `v`.

    EXAMPLES::

        sage: from sage.quadratic_forms.ternary import extend
        sage: v = (6, 4, 12)
        sage: m = extend(v)
        sage: M = matrix(3, m); M
        [ 6  1  0]
        [ 4  1  0]
        [12  0  1]
        sage: M.det()
        2
        sage: v = (-12, 20, 30)
        sage: m = extend(v)
        sage: M = matrix(3, m)
        sage: M
        [-12   1   0]
        [ 20  -2   1]
        [ 30   0  -7]
        sage: M.det()
        2

    TESTS::

        sage: M = matrix(3, extend( (0, 0, 1) ))
        sage: M.det()
        1
        sage: M.column(0)
        (0, 0, 1)
        sage: M = matrix(3, extend( (0, 0, -2) ))
        sage: M.det()
        2
        sage: M.column(0)
        (0, 0, -2)
    """
    if v[0] == v[1] == 0:
        if v[2] < 0:
            return v[0], 0, 1, v[1], 1, 0, v[2], 0, 0
        else:
            return v[0], 1, 0, v[1], 0, 1, v[2], 0, 0

    b1 = xgcd(v[0], v[1])
    b2 = xgcd(b1[1], b1[2])
    b3 = xgcd(b1[0], v[2])

    return v[0], -b1[2], -b2[1]*b3[2], v[1], b1[1], -b2[2]*b3[2], v[2], 0, b3[1]


def _find_p_neighbor_from_vec(a, b, c, r, s, t, p, v, mat=False):
    """
    Find the coefficients of the reduced equivalent of the p-neighbor
    of the ternary quadratic given by Q = (a, b, c, r, s, t)  form associated
    to a given vector v satisfying:

    1. Q(v) = 0  mod p

    2. v is a non-singular point of the conic Q(v) = 0 mod p.

    Reference:  Gonzalo Tornaria's Thesis, Thrm 3.5, p34.

    EXAMPLES::

        sage: from sage.quadratic_forms.ternary import _find_p_neighbor_from_vec
        sage: Q = TernaryQF([1, 3, 3, -2, 0, -1])
        sage: Q
        Ternary quadratic form with integer coefficients:
        [1 3 3]
        [-2 0 -1]
        sage: Q.disc()
        29
        sage: v = (9, 7, 1)
        sage: v in Q.find_zeros_mod_p(11)                                               # needs sage.libs.pari
        True
        sage: q11, M = _find_p_neighbor_from_vec(1, 3, 3, -2, 0, -1, 11, v, mat=True)
        sage: Q11 = TernaryQF(q11)
        sage: Q11
        Ternary quadratic form with integer coefficients:
        [1 2 4]
        [-1 -1 0]
        sage: M = matrix(3, M)
        sage: M
        [    -1  -5/11   7/11]
        [     0 -10/11   3/11]
        [     0  -3/11  13/11]
        sage: Q(M) == Q11
        True
    """
    v0, w0, u0, v1, w1, u1, v2, w2, u2 = extend(v)

    m00 = 2*a*v0**2 + 2*b*v1**2 + 2*c*v2**2 + 2*r*v1*v2 + 2*s*v0*v2 + 2*t*v0*v1
    m11 = 2*a*w0**2 + 2*b*w1**2 + 2*t*w0*w1
    m22 = 2*a*u0**2 + 2*b*u1**2 + 2*c*u2**2 + 2*r*u1*u2 + 2*s*u0*u2 + 2*t*u0*u1
    m01 = 2*a*v0*w0 + 2*b*v1*w1 + r*v2*w1 + s*v2*w0 + t*v0*w1 + t*v1*w0
    m02 = 2*a*u0*v0 + 2*b*u1*v1 + 2*c*u2*v2 + r*u1*v2 + r*u2*v1 + s*u0*v2 + s*u2*v0 + t*u0*v1 + t*u1*v0
    m12 = 2*a*u0*w0 + 2*b*u1*w1 + r*u2*w1 + s*u2*w0 + t*u0*w1 + t*u1*w0

    if m02 % p:

        m0 = (-m00/m02/2) % p**2
        m1 = (-m01/m02) % p

        b00 = m0**2*m22/p**2 + 2*m0*m02/p**2 + m00/p**2
        b11 = m1**2*m22 + 2*m1*m12 + m11
        b22 = m22*p**2
        b01 = m0*m1*m22/p + m0*m12/p + m02*m1/p + m01/p
        b02 = m0*m22 + m02
        b12 = m1*m22*p + m12*p

        if mat:
            q, Mr = _reduced_ternary_form_eisenstein_with_matrix(ZZ(b00/2), ZZ(b11/2), ZZ(b22/2), ZZ(b12), ZZ(b02), ZZ(b01))
            r00, r01, r02, r10, r11, r12, r20, r21, r22 = Mr.list()
            t00 = p*r20*u0 + (m0*u0/p + v0/p)*r00 + (m1*u0 + w0)*r10
            t01 = p*r21*u0 + (m0*u0/p + v0/p)*r01 + (m1*u0 + w0)*r11
            t02 = p*r22*u0 + (m0*u0/p + v0/p)*r02 + (m1*u0 + w0)*r12
            t10 = p*r20*u1 + (m0*u1/p + v1/p)*r00 + (m1*u1 + w1)*r10
            t11 = p*r21*u1 + (m0*u1/p + v1/p)*r01 + (m1*u1 + w1)*r11
            t12 = p*r22*u1 + (m0*u1/p + v1/p)*r02 + (m1*u1 + w1)*r12
            t20 = p*r20*u2 + (m0*u2/p + v2/p)*r00 + (m1*u2 + w2)*r10
            t21 = p*r21*u2 + (m0*u2/p + v2/p)*r01 + (m1*u2 + w2)*r11
            t22 = p*r22*u2 + (m0*u2/p + v2/p)*r02 + (m1*u2 + w2)*r12
            return q, (t00, t01, t02, t10, t11, t12, t20, t21, t22)
        else:
            return _reduced_ternary_form_eisenstein_without_matrix(ZZ(b00/2), ZZ(b11/2), ZZ(b22/2), ZZ(b12), ZZ(b02), ZZ(b01))

    if m01 % p:

        m0 = (-m00/m01/2) % p**2
        m1 = (-m02/m01) % p

        b00 = m0**2*m11/p**2 + 2*m0*m01/p**2 + m00/p**2
        b11 = m1**2*m11 + 2*m1*m12 + m22
        b22 = m11*p**2
        b01 = m0*m1*m11/p + m0*m12/p + m01*m1/p + m02/p
        b02 = m0*m11 + m01
        b12 = m1*m11*p + m12*p

        if mat:
            q, Mr = _reduced_ternary_form_eisenstein_with_matrix(ZZ(b00/2), ZZ(b11/2), ZZ(b22/2), ZZ(b12), ZZ(b02), ZZ(b01))
            r00, r01, r02, r10, r11, r12, r20, r21, r22 = Mr.list()
            t00 = p*r20*w0 + (m0*w0/p + v0/p)*r00 + (m1*w0 + u0)*r10
            t01 = p*r21*w0 + (m0*w0/p + v0/p)*r01 + (m1*w0 + u0)*r11
            t02 = p*r22*w0 + (m0*w0/p + v0/p)*r02 + (m1*w0 + u0)*r12
            t10 = p*r20*w1 + (m0*w1/p + v1/p)*r00 + (m1*w1 + u1)*r10
            t11 = p*r21*w1 + (m0*w1/p + v1/p)*r01 + (m1*w1 + u1)*r11
            t12 = p*r22*w1 + (m0*w1/p + v1/p)*r02 + (m1*w1 + u1)*r12
            t20 = p*r20*w2 + (m0*w2/p + v2/p)*r00 + (m1*w2 + u2)*r10
            t21 = p*r21*w2 + (m0*w2/p + v2/p)*r01 + (m1*w2 + u2)*r11
            t22 = p*r22*w2 + (m0*w2/p + v2/p)*r02 + (m1*w2 + u2)*r12
            return q, (t00, t01, t02, t10, t11, t12, t20, t21, t22)
        else:
            return _reduced_ternary_form_eisenstein_without_matrix(ZZ(b00/2), ZZ(b11/2), ZZ(b22/2), ZZ(b12), ZZ(b02), ZZ(b01))


def _basic_lemma_vec(a, b, c, r, s, t, n):
    """
    Find a vector v such that the ternary quadratic form given by (a, b, c, r, s, t) evaluated at v is
    coprime with n a prime or 1.

    EXAMPLES::

        sage: from sage.quadratic_forms.ternary import _basic_lemma_vec
        sage: Q = TernaryQF([5, 2, 3, -1, 0, 0])
        sage: v = _basic_lemma_vec(5, 2, 3, -1, 0, 0, 5)
        sage: v
        (0, 1, 0)
        sage: Q(v)
        2
    """
    if n == 1:
        return 0, 0, 0

    if a % n:
        return 1, 0, 0
    elif b % n:
        return 0, 1, 0
    elif c % n:
        return 0, 0, 1

    if r % n:
        return 0, 1, 1
    elif s % n:
        return 1, 0, 1
    elif t % n:
        return 1, 1, 0

    raise ValueError("not primitive form")


def _basic_lemma(a, b, c, r, s, t, n):
    """
    Find a number represented by the ternary quadratic form given by the coefficients (a, b, c, r, s, t)
    and coprime to the prime n.

    EXAMPLES::

        sage: from sage.quadratic_forms.ternary import _basic_lemma
        sage: Q = TernaryQF([5, 2, 3, -1, 0, 0])
        sage: _basic_lemma(5, 2, 3, -1, 0, 0, 5)
        2
    """
    if n == 1:
        return 0

    if a % n:
        return a
    elif b % n:
        return b
    elif c % n:
        return c

    if r % n:
        return b + c + r
    elif s % n:
        return a + c + s
    elif t % n:
        return a + b + t

    raise ValueError("not primitive form")
