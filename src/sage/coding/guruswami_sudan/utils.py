r"""
Guruswami-Sudan utility methods

AUTHORS:

- Johan S. R. Nielsen, original implementation (see [Nie]_ for details)
- David Lucas, ported the original implementation in Sage
"""

#*****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#                     2015 Johan S. R. Nielsen <jsrn@jsrn.dk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.functions.other import binomial, floor, sqrt
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.arith.all import lcm
from sage.combinat.permutation import Permutation

def polynomial_to_list(p, len):
    r"""
    Returns ``p`` as a list of its coefficients of length ``len``.

    INPUT:

    - ``p`` -- a polynomial

    - ``len`` -- an integer. If ``len`` is smaller than the degree of ``p``, the
      returned list will be of size degree of ``p``, else it will be of size ``len``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import polynomial_to_list
        sage: F.<x> = GF(41)[]
        sage: p = 9*x^2 + 8*x + 37
        sage: polynomial_to_list(p, 4)
        [37, 8, 9, 0]
    """
    return list(p) + [0]*max(0, len-p.degree()-1)

def johnson_radius(n, d):
    r"""
    Returns the Johnson-radius for the code length `n` and the minimum distance `d`.

    The Johnson radius is defined as `n - \sqrt(n(n-d))`.

    INPUT:

    - ``n`` -- an integer, the length of the code
    - ``d`` -- an integer, the minimum distance of the code

    EXAMPLES::

        sage: sage.coding.guruswami_sudan.utils.johnson_radius(250, 181)
        -5*sqrt(690) + 250
    """
    return n - sqrt(n*(n-d))

def ligt(x):
    r"""
    Returns the least integer greater than ``x``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import ligt
        sage: ligt(41)
        42

    It works with any type of numbers (not only integers)::

        sage: ligt(41.041)
        42
    """
    return floor(x+1)

def gilt(x):
    r"""
    Returns the greatest integer smaller than ``x``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import gilt
        sage: gilt(43)
        42

    It works with any type of numbers (not only integers)::

        sage: gilt(43.041)
        43
    """
    if x in ZZ:
        return Integer(x-1)
    else:
        return floor(x)

def solve_degree2_to_integer_range(a,b,c):
    r"""
    Returns the greatest integer range `[i_1, i_2]` such that
    `i_1 > x_1` and `i_2 < x_2` where `x_1, x_2` are the two zeroes of the equation in `x`:
    `ax^2+bx+c=0`.

    If there is no real solution to the equation, it returns an empty range with negative coefficients.

    INPUT:

    - ``a``, ``b`` and ``c`` -- coefficients of a second degree equation, ``a`` being the coefficient of
      the higher degree term.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import solve_degree2_to_integer_range
        sage: solve_degree2_to_integer_range(1, -5, 1)
        (1, 4)

    If there is no real solution::

        sage: solve_degree2_to_integer_range(50, 5, 42)
        (-2, -1)
    """
    D = b**2 - 4*a*c
    if D < 0:
        return (-2,-1)
    sD = float(sqrt(D))
    minx, maxx = (-b-sD)/2.0/a , (-b+sD)/2.0/a
    mini, maxi = (ligt(minx), gilt(maxx))
    if mini > maxi:
        return (-2,-1)
    else:
        return (mini,maxi)

def apply_shifts(M, shifts):
    r"""
    Applies column shifts inplace to the polynomial matrix `M`.

    This is equivalent to multiplying the `n`th column of `M` with
    `x^{shifts[n]}`.

    INPUT:

    - ``M`` -- a polynomial matrix

    - ``shifts`` -- a list of non-negative integer shifts

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import apply_shifts
        sage: F.<x> = GF(7)[]
        sage: M = matrix(F, [[2*x^2 + x, 5*x^2 + 2*x + 1, 4*x^2 + x],\
                             [x^2 + 3*x + 3, 5*x^2 + 5*x + 1, 6*x^2 + 5*x + 4],\
                             [5*x^2 + 2*x + 4, 4*x^2 + 2*x, 5*x^2 + x + 2]])
        sage: shifts = [1, 2, 3]
        sage: apply_shifts(M, shifts)
        sage: M
        [          2*x^3 + x^2   5*x^4 + 2*x^3 + x^2           4*x^5 + x^4]
        [    x^3 + 3*x^2 + 3*x   5*x^4 + 5*x^3 + x^2 6*x^5 + 5*x^4 + 4*x^3]
        [  5*x^3 + 2*x^2 + 4*x         4*x^4 + 2*x^3   5*x^5 + x^4 + 2*x^3]
    """
    x = M.base_ring().gen()
    for j in range(M.ncols()):
        M.set_col_to_multiple_of_col(j,j, x**shifts[j])

def remove_shifts(M, shifts):
    r"""
    Removes the shifts inplace to the matrix `M` as they were introduced by
    :func:`apply_shifts`.

    If `M` was not earlier called with :func:`apply_shifts` using the same
    shifts, then the least significant coefficients of the entries of `M`,
    corresponding to how much we are shifting down, will be lost.

    INPUT:

    - ``M`` -- a polynomial matrix

    - ``shifts`` -- a list of non-negative integer shifts

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import remove_shifts
        sage: F.<x> = GF(7)[]
        sage: M = matrix(F, [[2*x^3 + x^2, 5*x^4 + 2*x^3 + x^2, 4*x^5 + x^4],\
                            [x^3 + 3*x^2 + 3*x, 5*x^4 + 5*x^3 + x^2, 6*x^5 + 5*x^4 + 4*x^3],\
                            [5*x^3 + 2*x^2 + 4*x, 4*x^4 + 2*x^3, 5*x^5 + x^4 + 2*x^3]])

        sage: shifts = [1, 2, 3]
        sage: remove_shifts(M, shifts)
        sage: M
        [      2*x^2 + x 5*x^2 + 2*x + 1       4*x^2 + x]
        [  x^2 + 3*x + 3 5*x^2 + 5*x + 1 6*x^2 + 5*x + 4]
        [5*x^2 + 2*x + 4     4*x^2 + 2*x   5*x^2 + x + 2]
    """
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            M[i,j] = M[i,j].shift(-shifts[j])

def _leading_position(v, shifts=None):
    r"""
    Returns the position of the highest-degree term of ``v``.

    This methods can manage shifted degree, by providing ``shift`` to it.

    In case of several positions having the same, highest degree, the right-most
    position is given.

    INPUT:

    - ``v`` -- a vector of polynomials

    - ``shifts`` -- (default: ``None``) a list of integer shifts to consider ``v`` under.
      If ``None``, all shifts are considered as ``0``.

    EXAMPLES::

    sage: from sage.coding.guruswami_sudan.utils import _leading_position
    sage: F.<x> = GF(7)[]
    sage: v = vector(F, [3*x^2 + 3*x + 4, 4*x + 3, 4*x^2 + 4*x + 5, x^2 + 2*x + 5, 3*x^2 + 5*x])
    sage: _leading_position(v)
    4
    """
    if not shifts:
        shifts=[0]*len(v)
    best=-1
    bestp=-1
    for p in range(0,len(v)):
        if not v[p].is_zero():
            vpdeg = v[p].degree() + shifts[p]
            if vpdeg >= best:
                best=vpdeg
                bestp = p
    if best==-1:
        return -1
    else:
        return bestp

def leading_term(v, shifts=None):
    r"""
    Returns the term of ``v`` with the highest degree.

    This methods can manage shifted degree, by providing ``shift`` to it.

    In case of several positions having the same, highest degree, the term with
    the right-most position is returned.

    INPUT:

    - ``v`` -- a vector of polynomials

    - ``shifts`` -- (default: ``None``) a list of integer shifts to consider ``v`` under.
      If ``None``, all shifts are considered as ``0``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import leading_term
        sage: F.<x> = GF(7)[]
        sage: v = vector(F, [3*x^2 + 3*x + 4, 4*x + 3, 4*x^2 + 4*x + 5, x^2 + 2*x + 5, 3*x^2 + 5*x])
        sage: leading_term(v)
        3*x^2 + 5*x
    """
    return v[_leading_position(v, shifts=shifts)]
