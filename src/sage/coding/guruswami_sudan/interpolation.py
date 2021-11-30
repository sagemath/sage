"""
Interpolation algorithms for the Guruswami-Sudan decoder

AUTHORS:

- Johan S. R. Nielsen, original implementation (see [Nie]_ for details)
- David Lucas, ported the original implementation in Sage
"""

# ****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#                     2015 Johan S. R. Nielsen <jsrn@jsrn.dk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.arith.misc import binomial
from sage.matrix.constructor import matrix
from sage.misc.misc_c import prod

####################### Linear algebra system solving ###############################
def _flatten_once(lstlst):
    r"""
    Flattens a list of list into a list, but only flattening one layer and
    returns a generator.

    This is similar to Python's ``flatten`` method, except that here, if you
    provide a list of lists of lists (and so on), it returns a list of lists
    (and so on) and not a list.

    INPUT:

    - ``lstlst`` -- a list of lists.

    EXAMPLES::

    sage: from sage.coding.guruswami_sudan.interpolation import _flatten_once
    sage: ll = [[1,2], [3,4], [5,6]]
    sage: list(_flatten_once(ll))
    [1, 2, 3, 4, 5, 6]
    """
    for lst in lstlst:
        for e in lst:
            yield e

#*************************************************************
#  Linear algebraic Interpolation algorithm, helper functions
#*************************************************************

def _monomial_list(maxdeg, l, wy):
    r"""
    Returns a list of all non-negative integer pairs `(i,j)` such that ``i + wy
    * j < maxdeg`` and ``j \geq l``.

    INPUT:

    - ``maxdeg``, ``l``, ``wy`` -- integers.

    OUTPUT:

    - a list of pairs of integers.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.interpolation import _monomial_list
        sage: list(_monomial_list(8, 1, 3))
        [(0, 0),
         (1, 0),
         (2, 0),
         (3, 0),
         (4, 0),
         (5, 0),
         (6, 0),
         (7, 0),
         (0, 1),
         (1, 1),
         (2, 1),
         (3, 1),
         (4, 1)]
    """
    for y in range(l + 1):
        for x in range(maxdeg - y * wy):
            yield (x, y)


def _interpolation_matrix_given_monomials(points, s, monomials):
    r"""
    Returns a matrix whose nullspace is a basis for all interpolation
    polynomials, each polynomial having its coefficients laid out according to
    the given list of monomials.

    The output is an `S \times T` matrix, where `T` is the length of
    ``monomials``, and `S = s(s+1)/2`. Its ``i``-th column will be the
    coefficients on the ``i``-th monomial in ``monomials``.

    INPUT:

    - ``points`` -- a list of pairs of field elements, the interpolation points.

    - ``s`` -- an integer, the multiplicity parameter from Guruswami-Sudan algorithm.

    - ``monomials`` -- a list of monomials, each represented by the powers as an integer pair `(i,j)`.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.interpolation import _interpolation_matrix_given_monomials
        sage: F = GF(11)
        sage: points = [ (F(0), F(1)), (F(1), F(5)) ]
        sage: s = 2
        sage: monomials = [(0, 0), (1, 0), (1, 1), (0, 2) ]
        sage: _interpolation_matrix_given_monomials(points, s, monomials)
        [ 1  0  0  1]
        [ 0  0  0  2]
        [ 0  1  1  0]
        [ 1  1  5  3]
        [ 0  0  1 10]
        [ 0  1  5  0]
    """
    def eqs_affine(x0, y0):
        r"""
        Make equation for the affine point x0, y0. Return a list of
        equations, each equation being a list of coefficients corresponding to
        the monomials in ``monomials``.
        """
        eqs = []
        for i in range(s):
            for j in range(s - i):
                eq = dict()
                for monomial in monomials:
                    ihat = monomial[0]
                    jhat = monomial[1]
                    if ihat >= i and jhat >= j:
                        icoeff = binomial(ihat, i) * x0**(ihat-i) \
                                    if ihat > i else 1
                        jcoeff = binomial(jhat, j) * y0**(jhat-j) \
                                    if jhat > j else 1
                        eq[monomial] = jcoeff * icoeff
                eqs.append([eq.get(monomial, 0) for monomial in monomials])
        return eqs
    return matrix(list(_flatten_once([eqs_affine(*point) for point in points])))


def _interpolation_max_weighted_deg(n, tau, s):
    """Return the maximal weighted degree allowed for an interpolation
    polynomial over `n` points, correcting `tau` errors and with multiplicity
    `s`

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.interpolation import _interpolation_max_weighted_deg
        sage: _interpolation_max_weighted_deg(10, 3, 5)
        35
    """
    return (n-tau) * s

def _interpolation_matrix_problem(points, tau, parameters, wy):
    r"""
    Returns the linear system of equations which ``Q`` should be a solution to.

    This linear system is returned as a matrix ``M`` and a list of monomials ``monomials``,
    where a vector in the right nullspace of ``M`` corresponds to an
    interpolation polynomial `Q`, by mapping the `t`'th element of such a vector
    to the coefficient to `x^iy^j`, where `(i,j)` is the `t`'th element of ``monomials``.

    INPUT:

    - ``points`` -- a list of interpolation points, as pairs of field elements.

    - ``tau`` -- an integer, the number of errors one wants to decode.

    - ``parameters`` -- (default: ``None``) a pair of integers, where:
        - the first integer is the multiplicity parameter of Guruswami-Sudan algorithm and
        - the second integer is the list size parameter.

    - ``wy`` -- an integer specifying the `y`-weighted degree that is to be
      minimised in the interpolation polynomial. In Guruswami-Sudan, this is
      `k-1`, where `k` is the dimension of the GRS code.

    EXAMPLES:

    The following parameters arise from Guruswami-Sudan decoding of an [6,2,5]
    GRS code over F(11) with multiplicity 2 and list size 4.

        sage: from sage.coding.guruswami_sudan.interpolation import _interpolation_matrix_problem
        sage: F = GF(11)
        sage: points = [(F(x), F(y)) for (x,y) in [(0, 5), (1, 1), (2, 4), (3, 6), (4, 3), (5, 3)]]
        sage: tau = 3
        sage: params = (2, 4)
        sage: wy = 1
        sage: _interpolation_matrix_problem(points, tau, params, wy)
        (
        [ 1  0  0  0  0  0  5  0  0  0  0  3  0  0  0  4  0  0  9  0]
        [ 0  0  0  0  0  0  1  0  0  0  0 10  0  0  0  9  0  0  5  0]
        [ 0  1  0  0  0  0  0  5  0  0  0  0  3  0  0  0  4  0  0  9]
        [ 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1]
        [ 0  0  0  0  0  0  1  1  1  1  1  2  2  2  2  3  3  3  4  4]
        [ 0  1  2  3  4  5  0  1  2  3  4  0  1  2  3  0  1  2  0  1]
        [ 1  2  4  8  5 10  4  8  5 10  9  5 10  9  7  9  7  3  3  6]
        [ 0  0  0  0  0  0  1  2  4  8  5  8  5 10  9  4  8  5  3  6]
        [ 0  1  4  1 10  3  0  4  5  4  7  0  5  9  5  0  9  3  0  3]
        [ 1  3  9  5  4  1  6  7 10  8  2  3  9  5  4  7 10  8  9  5]
        [ 0  0  0  0  0  0  1  3  9  5  4  1  3  9  5  9  5  4  6  7]
        [ 0  1  6  5  9  9  0  6  3  8 10  0  3  7  4  0  7  9  0  9]
        [ 1  4  5  9  3  1  3  1  4  5  9  9  3  1  4  5  9  3  4  5]
        [ 0  0  0  0  0  0  1  4  5  9  3  6  2  8 10  5  9  3  9  3]
        [ 0  1  8  4  3  4  0  3  2  1  9  0  9  6  3  0  5  7  0  4]
        [ 1  5  3  4  9  1  3  4  9  1  5  9  1  5  3  5  3  4  4  9]
        [ 0  0  0  0  0  0  1  5  3  4  9  6  8  7  2  5  3  4  9  1]
        [ 0  1 10  9  5  1  0  3  8  5  4  0  9  2  4  0  5  6  0  4], [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (0, 2), (1, 2), (2, 2), (3, 2), (0, 3), (1, 3), (2, 3), (0, 4), (1, 4)]
        )
    """
    s, l = parameters[0], parameters[1]
    monomials = list(_monomial_list(_interpolation_max_weighted_deg(len(points), tau, s), l, wy))
    M = _interpolation_matrix_given_monomials(points, s, monomials)
    return (M, monomials)


def gs_interpolation_linalg(points, tau, parameters, wy):
    r"""
    Compute an interpolation polynomial Q(x,y) for the Guruswami-Sudan algorithm
    by solving a linear system of equations.

    ``Q`` is a bivariate polynomial over the field of the points, such that the
    polynomial has a zero of multiplicity at least `s` at each of the points,
    where `s` is the multiplicity parameter. Furthermore, its ``(1,
    wy)``-weighted degree should be less than
    ``_interpolation_max_weighted_deg(n, tau, wy)``, where ``n`` is the number
    of points

    INPUT:

    - ``points`` -- a list of tuples ``(xi, yi)`` such that we seek ``Q`` with
      ``(xi,yi)`` being a root of ``Q`` with multiplicity ``s``.

    - ``tau`` -- an integer, the number of errors one wants to decode.

    - ``parameters`` -- (default: ``None``) a pair of integers, where:
        - the first integer is the multiplicity parameter of Guruswami-Sudan algorithm and
        - the second integer is the list size parameter.

    - ``wy`` -- an integer, the `y`-weight, where we seek ``Q`` of low
      ``(1,wy)`` weighted degree.

    EXAMPLES:

    The following parameters arise from Guruswami-Sudan decoding of an [6,2,5]
    GRS code over F(11) with multiplicity 2 and list size 4::

        sage: from sage.coding.guruswami_sudan.interpolation import gs_interpolation_linalg
        sage: F = GF(11)
        sage: points = [(F(x),F(y)) for (x,y) in [(0, 5), (1, 1), (2, 4), (3, 6), (4, 3), (5, 3)]]
        sage: tau = 3
        sage: params = (2, 4)
        sage: wy = 1
        sage: Q = gs_interpolation_linalg(points, tau, params, wy); Q
        4*x^5 - 4*x^4*y - 2*x^2*y^3 - x*y^4 + 3*x^4 - 4*x^2*y^2 + 5*y^4 - x^3 + x^2*y + 5*x*y^2 - 5*y^3 + 3*x*y - 2*y^2 + x - 4*y + 1

    We verify that the interpolation polynomial has a zero of multiplicity at least 2 in each point::

        sage: all( Q(x=a, y=b).is_zero() for (a,b) in points )
        True
        sage: x,y = Q.parent().gens()
        sage: dQdx = Q.derivative(x)
        sage: all( dQdx(x=a, y=b).is_zero() for (a,b) in points )
        True
        sage: dQdy = Q.derivative(y)
        sage: all( dQdy(x=a, y=b).is_zero() for (a,b) in points )
        True
    """
    M, monomials = _interpolation_matrix_problem(points, tau, parameters, wy)
    Ker = M.right_kernel()
    # Pick a non-zero element from the right kernel
    sol = Ker.basis()[0]
    # Construct the Q polynomial
    PF = M.base_ring()['x', 'y'] #make that ring a ring in <x>
    x, y = PF.gens()
    Q = sum([x**monomials[i][0] * y**monomials[i][1] * sol[i] for i in range(0, len(monomials))])
    return Q

####################### Lee-O'Sullivan's method ###############################

def lee_osullivan_module(points, parameters, wy):
    r"""
    Returns the analytically straight-forward basis for the `\GF q[x]` module
    containing all interpolation polynomials, as according to Lee and
    O'Sullivan.

    The module is constructed in the following way: Let `R(x)` be the Lagrange
    interpolation polynomial through the sought interpolation points `(x_i,
    y_i)`, i.e. `R(x_i) = y_i`. Let `G(x) = \prod_{i=1}^n (x-x_i)`. Then the
    `i`'th row of the basis matrix of the module is the coefficient-vector of
    the following polynomial in `\GF q[x][y]`:

        `P_i(x,y) = G(x)^{[i-s]} (y - R(x))^{i - [i-s]} y^{[i-s]}` ,

    where `[a]` for real `a` is `a` when `a > 0` and 0 otherwise. It is easily
    seen that `P_i(x,y)` is an interpolation polynomial, i.e. it is zero with
    multiplicity at least `s` on each of the points `(x_i, y_i)`.


    INPUT:

    - ``points`` -- a list of tuples ``(xi, yi)`` such that we seek ``Q`` with
      ``(xi,yi)`` being a root of ``Q`` with multiplicity ``s``.

    - ``parameters`` -- (default: ``None``) a pair of integers, where:
        - the first integer is the multiplicity parameter `s` of Guruswami-Sudan algorithm and
        - the second integer is the list size parameter.

    - ``wy`` -- an integer, the `y`-weight, where we seek ``Q`` of low
      ``(1,wy)`` weighted degree.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.interpolation import lee_osullivan_module
        sage: F = GF(11)
        sage: points = [(F(0), F(2)), (F(1), F(5)), (F(2), F(0)), (F(3), F(4)), (F(4), F(9))\
                , (F(5), F(1)), (F(6), F(9)), (F(7), F(10))]
        sage: params = (1, 1)
        sage: wy = 1
        sage: lee_osullivan_module(points, params, wy)
        [x^8 + 5*x^7 + 3*x^6 + 9*x^5 + 4*x^4 + 2*x^3 + 9*x   0]
        [ 10*x^7 + 4*x^6 + 9*x^4 + 7*x^3 + 2*x^2 + 9*x + 9   1]
    """
    s, l = parameters[0], parameters[1]
    F = points[0][0].parent()
    PF = F['x']
    x = PF.gens()[0]
    R = PF.lagrange_polynomial(points)
    G = prod(x - points[i][0] for i in range(0, len(points)))
    PFy = PF['y']
    y = PFy.gens()[0]
    ybasis = [(y-R)**i * G**(s-i) for i in range(0, s+1)] \
            + [y**(i-s) * (y-R)**s for i in range(s+1, l+1)]
    def pad(lst):
        return lst + [0]*(l+1-len(lst))
    modbasis = [pad(yb.coefficients(sparse=False)) for yb in ybasis]
    return matrix(PF, modbasis)

def gs_interpolation_lee_osullivan(points, tau, parameters, wy):
    r"""
    Returns an interpolation polynomial Q(x,y) for the given input using the
    module-based algorithm of Lee and O'Sullivan.

    This algorithm constructs an explicit `(\ell+1) \times (\ell+1)` polynomial
    matrix whose rows span the `\GF q[x]` module of all interpolation
    polynomials. It then runs a row reduction algorithm to find a low-shifted
    degree vector in this row space, corresponding to a low weighted-degree
    interpolation polynomial.

    INPUT:

    - ``points`` -- a list of tuples ``(xi, yi)`` such that we seek ``Q`` with
      ``(xi,yi)`` being a root of ``Q`` with multiplicity ``s``.

    - ``tau`` -- an integer, the number of errors one wants to decode.

    - ``parameters`` -- (default: ``None``) a pair of integers, where:
        - the first integer is the multiplicity parameter of Guruswami-Sudan algorithm and
        - the second integer is the list size parameter.

    - ``wy`` -- an integer, the `y`-weight, where we seek ``Q`` of low
      ``(1,wy)`` weighted degree.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.interpolation import gs_interpolation_lee_osullivan
        sage: F = GF(11)
        sage: points = [(F(0), F(2)), (F(1), F(5)), (F(2), F(0)), (F(3), F(4)), (F(4), F(9))\
                , (F(5), F(1)), (F(6), F(9)), (F(7), F(10))]
        sage: tau = 1
        sage: params = (1, 1)
        sage: wy = 1
        sage: Q = gs_interpolation_lee_osullivan(points, tau, params, wy)
        sage: Q / Q.lc() # make monic
        x^3*y + 2*x^3 - x^2*y + 5*x^2 + 5*x*y - 5*x + 2*y - 4
    """
    from .utils import _degree_of_vector
    s, l = parameters[0], parameters[1]
    F = points[0][0].parent()
    M = lee_osullivan_module(points, (s,l), wy)
    shifts = [i * wy for i in range(0,l+1)]
    Mnew = M.reduced_form(shifts=shifts)
    # Construct Q as the element of the row with the lowest weighted degree
    Qlist = min(Mnew.rows(), key=lambda r: _degree_of_vector(r, shifts))
    PFxy = F['x,y']
    xx, yy = PFxy.gens()
    Q = sum(yy**i * PFxy(Qlist[i]) for i in range(0,l+1))
    return Q
