# sage.doctest: needs scipy sage.graphs sage.groups
r"""
Riemann matrices and endomorphism rings of algebraic Riemann surfaces

This module provides a class, :class:`RiemannSurface`, to model the
Riemann surface determined by a plane algebraic curve over a subfield
of the complex numbers.

A homology basis is derived from the edges of a Voronoi cell decomposition based
on the branch locus. The pull-back of these edges to the Riemann surface
provides a graph on it that contains a homology basis.

The class provides methods for computing the Riemann period matrix of the
surface numerically, using a certified homotopy continuation method due to
[Kr2016]_.

The class also provides facilities for computing the endomorphism ring of the
period lattice numerically, by determining integer (near) solutions to the
relevant approximate linear equations.

One can also calculate the Abel-Jacobi map on the Riemann surface, and there
is basic functionality to interface with divisors of curves to facilitate this.

AUTHORS:

- Alexandre Zotine, Nils Bruin (2017-06-10): initial version
- Nils Bruin, Jeroen Sijsling (2018-01-05): algebraization, isomorphisms
- Linden Disney-Hogg, Nils Bruin (2021-06-23): efficient integration
- Linden Disney-Hogg, Nils Bruin (2022-09-07): Abel-Jacobi map

EXAMPLES:

We compute the Riemann matrix of a genus 3 curve::

    sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
    sage: R.<x,y> = QQ[]
    sage: f = x^4-x^3*y+2*x^3+2*x^2*y+2*x^2-2*x*y^2+4*x*y-y^3+3*y^2+2*y+1
    sage: S = RiemannSurface(f, prec=100)
    sage: M = S.riemann_matrix()

We test the usual properties, i.e., that the period matrix is symmetric and that
the imaginary part is positive definite::

    sage: all(abs(a) < 1e-20 for a in (M-M.T).list())
    True
    sage: iM = Matrix(RDF,3,3,[a.imag_part() for a in M.list()])
    sage: iM.is_positive_definite()
    True

We compute the endomorphism ring and check it has `\ZZ`-rank 6::

    sage: A = S.endomorphism_basis(80,8)
    sage: len(A) == 6
    True

In fact it is an order in a number field::

    sage: T.<t> = QQ[]
    sage: K.<a> = NumberField(t^6 - t^5 + 2*t^4 + 8*t^3 - t^2 - 5*t + 7)
    sage: all(len(a.minpoly().roots(K)) == a.minpoly().degree() for a in A)
    True

We can look at an extended example of the Abel-Jacobi functionality. We will
demonstrate a particular half-canonical divisor on Klein's Curve, known in
the literature::

    sage: f = x^3*y + y^3 + x
    sage: S = RiemannSurface(f, integration_method='rigorous')
    sage: BL = S.places_at_branch_locus(); BL
    [Place (x, y, y^2),
     Place (x^7 + 27/4, y + 4/9*x^5, y^2 + 4/3*x^3),
     Place (x^7 + 27/4, y - 2/9*x^5, y^2 + 1/3*x^3)]

We can read off out the output of ``places_at_branch_locus`` to choose our
divisor, and we can calculate the canonical divisor using curve functionality::

    sage: P0 = 1*BL[0]
    sage: from sage.schemes.curves.constructor import Curve
    sage: C = Curve(f)
    sage: F = C.function_field()
    sage: K = (F(x).differential()).divisor() - F(f.derivative(y)).divisor()
    sage: Pinf, Pinf_prime = C.places_at_infinity()
    sage: if K-3*Pinf-1*Pinf_prime: Pinf, Pinf_prime = (Pinf_prime, Pinf);
    sage: D = P0 + 2*Pinf - Pinf_prime

Note we could check using exact techniques that `2D = K`::

    sage: Z = K - 2*D
    sage: (Z.degree() == 0, len(Z.basis_differential_space()) == S.genus, len(Z.basis_function_space()) == 1)
    (True, True, True)

We can also check this using our Abel-Jacobi functions::

    sage: avoid = C.places_at_infinity()
    sage: Zeq, _ = S.strong_approximation(Z, avoid)
    sage: Zlist = S.divisor_to_divisor_list(Zeq)
    sage: AJ = S.abel_jacobi(Zlist)  # long time (1 second)
    sage: S.reduce_over_period_lattice(AJ).norm() < 1e-10  # long time
    True

REFERENCES:

The initial version of this code was developed alongside [BSZ2019]_.
"""
# ****************************************************************************
#       Copyright (C) 2017 Alexandre Zotine, Nils Bruin
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from scipy.spatial import Voronoi
from sage.arith.functions import lcm
from sage.arith.misc import GCD, algdep
from sage.ext.fast_callable import fast_callable
from sage.graphs.graph import Graph
from sage.groups.matrix_gps.finitely_generated import MatrixGroup
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.matrix.constructor import Matrix
from sage.matrix.special import block_matrix
from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten
from sage.misc.functional import numerical_approx
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_integer import IntegerLattice
from sage.numerical.gauss_legendre import integrate_vector, integrate_vector_N
from sage.rings.complex_mpfr import ComplexField, CDF
from sage.rings.function_field.constructor import FunctionField
from sage.rings.function_field.divisor import FunctionFieldDivisor
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RealField
from sage.schemes.curves.constructor import Curve
import sage.libs.mpmath.all as mpall

lazy_import('sage.rings.qqbar', 'number_field_elements_from_algebraics')


def voronoi_ghost(cpoints, n=6, CC=CDF):
    r"""
    Convert a set of complex points to a list of real tuples `(x,y)`, and
    appends n points in a big circle around them.

    The effect is that, with n >= 3, a Voronoi decomposition will have only
    finite cells around the original points. Furthermore, because the extra
    points are placed on a circle centered on the average of the given points,
    with a radius 3/2 times the largest distance between the center and the
    given points, these finite cells form a simply connected region.

    INPUT:

    - ``cpoints`` -- list of complex numbers

    OUTPUT:

    A list of real tuples `(x,y)` consisting of the original points and a set of
    points which surround them.

    EXAMPLES::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import voronoi_ghost
        sage: L = [1 + 1*I, 1 - 1*I, -1 + 1*I, -1 - 1*I]
        sage: voronoi_ghost(L)  # abs tol 1e-6
        [(1.0, 1.0),
         (1.0, -1.0),
         (-1.0, 1.0),
         (-1.0, -1.0),
         (2.121320343559643, 0.0),
         (1.0606601717798216, 1.8371173070873836),
         (-1.060660171779821, 1.8371173070873839),
         (-2.121320343559643, 2.59786816870648e-16),
         (-1.0606601717798223, -1.8371173070873832),
         (1.06066017177982, -1.8371173070873845)]
    """
    cpoints = [CC(c) for c in cpoints]
    average = sum(cpoints) / len(cpoints)
    if len(cpoints) == 1:
        radius = 1
    else:
        radius = 3 * max(abs(c - average) for c in cpoints) / 2
    z = CC.zeta(n)
    extra_points = [average + radius * z**i for i in range(n)]
    return [tuple(c) for c in cpoints + extra_points]


def bisect(L, t):
    r"""
    Find position in a sorted list using bisection.

    Given a list `L = [(t_0,...),(t_1,...),...(t_n,...)]` with increasing `t_i`,
    find the index i such that `t_i <= t < t_{i+1}` using bisection. The rest of
    the tuple is available for whatever use required.

    INPUT:

    - ``L`` -- list of tuples such that the first term of each tuple is a real
      number between 0 and 1. These real numbers must be increasing

    - ``t`` -- real number between `t_0` and `t_n`

    OUTPUT: integer i, giving the position in L where t would be in

    EXAMPLES:

    Form a list of the desired form, and pick a real number between 0 and 1::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import bisect
        sage: L = [(0.0, 'a'), (0.3, 'b'), (0.7, 'c'), (0.8, 'd'), (0.9, 'e'), (1.0, 'f')]
        sage: t = 0.5
        sage: bisect(L,t)
        1

    Another example which demonstrates that if t is equal to one of the t_i, it
    returns that index::

        sage: L = [(0.0, 'a'), (0.1, 'b'), (0.45, 'c'), (0.5, 'd'), (0.65, 'e'), (1.0, 'f')]
        sage: t = 0.5
        sage: bisect(L,t)
        3
    """
    # Defining starting indices for the loop.
    min = 0
    max = len(L) - 1
    # If the input t is not between 0 and 1, raise an error.
    if t < L[min][0] or t > L[max][0]:
        raise ValueError("value for t out of range")
    # Main loop.
    while min < max - 1:
        # Bisect.
        mid = (max + min) // 2
        # If it's equal, return the index we bisected to.
        if t == L[mid][0]:
            return mid
        # If it's smaller, then we're on the left side.
        elif t < L[mid][0]:
            max = mid
        # Otherwise we're on the right side.
        else:
            min = mid
    # Once the loop terminates, we return what the indices converged to.
    return min


def numerical_inverse(C):
    """
    Compute numerical inverse of a matrix via LU decomposition.

    INPUT:

    - ``C`` -- a real or complex invertible square matrix

    EXAMPLES::

        sage: C = matrix(CC, 3, 3, [-4.5606e-31 + 1.2326e-31*I,
        ....:                       -0.21313 + 0.24166*I,
        ....:                       -3.4513e-31 + 0.16111*I,
        ....:                       -1.0175 + 9.8608e-32*I,
        ....:                       0.30912 + 0.19962*I,
        ....:                       -4.9304e-32 + 0.39923*I,
        ....:                       0.96793 - 3.4513e-31*I,
        ....:                       -0.091587 + 0.19276*I,
        ....:                       3.9443e-31 + 0.38552*I])
        sage: from sage.schemes.riemann_surfaces.riemann_surface import numerical_inverse
        sage: 3e-16 < (C^-1*C-C^0).norm() < 1e-15
        True
        sage: (numerical_inverse(C)*C-C^0).norm() < 3e-16
        True
    """
    R = C.parent()
    prec = R.base_ring().prec()
    mpall.mp.prec = prec
    with mpall.workprec(prec):
        Cmp = mpall.matrix([mpall.sage_to_mpmath(list(c), prec) for c in C])
        PLU = mpall.lu(Cmp)
    P, L, U = (R([mpall.mpmath_to_sage(c, prec) for c in M]) for M in PLU)
    return U.inverse() * L.inverse() * P


class ConvergenceError(ValueError):
    r"""
    Error object suitable for raising and catching when Newton iteration fails.

    EXAMPLES::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import ConvergenceError
        sage: raise ConvergenceError("test")
        Traceback (most recent call last):
        ...
        ConvergenceError: test
        sage: isinstance(ConvergenceError(),ValueError)
        True
    """
    pass


def differential_basis_baker(f):
    r"""
    Compute a differential basis for a curve that is nonsingular outside (1:0:0),(0:1:0),(0:0:1).

    Baker's theorem tells us that if a curve has its singularities at the coordinate vertices and meets
    some further easily tested genericity criteria,
    then we can read off a basis for the regular differentials from the interior of the
    Newton polygon spanned by the monomials. While this theorem only applies to special plane curves
    it is worth implementing because the analysis is relatively cheap and it applies to a lot of
    commonly encountered curves (e.g., curves given by a hyperelliptic model). Other advantages include
    that we can do the computation over any exact base ring (the alternative Singular based method for
    computing the adjoint ideal requires the rationals), and that we can avoid being affected by subtle bugs
    in the Singular code.

    ``None`` is returned when ``f`` does not describe a curve of the relevant type. If ``f`` is of the relevant
    type, but is of genus `0` then ``[]`` is returned (which are both False values, but they are not equal).

    INPUT:

    - ``f`` -- a bivariate polynomial

    EXAMPLES::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import differential_basis_baker
        sage: R.<x,y> = QQ[]
        sage: f = x^3 + y^3 + x^5*y^5
        sage: differential_basis_baker(f)
        [y^2, x*y, x*y^2, x^2, x^2*y, x^2*y^2, x^2*y^3, x^3*y^2, x^3*y^3]
        sage: f = y^2 - (x-3)^2*x
        sage: differential_basis_baker(f) is None
        True
        sage: differential_basis_baker(x^2+y^2-1)
        []

    TESTS::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import differential_basis_baker
        sage: R.<x,y> = QQ[]
        sage: f = y^12 - x*(x - 1)^7
        sage: differential_basis_baker(f) is None
        True
    """
    k = f.base_ring()
    R = PolynomialRing(k, 3, "x,y,z")
    x, y, z = R.gens()
    F = f(x / z, y / z).numerator()
    W = [F] + [F.derivative(v) for v in R.gens()]
    # we check that the singularities lie at (1:0:0),(0:1:0),(0:0:1)
    # by checking that the eliminations of x, y, z result in
    # (principal) ideals generated by a monomial. This is a sufficient
    # condition, but not completely necessary.
    # It's cheap to check, though.
    for c in R.gens():
        B = GCD([W[i].resultant(W[j], c) for i in range(4) for j in range(i)])
        if len(B.monomials()) > 1:
            return None
    from sage.geometry.polyhedron.constructor import Polyhedron

    D = {(k[0], k[1]): v for k, v in f.monomial_coefficients().items()}
    P = Polyhedron(D)
    kT = k["t"]
    # here we check the additional genericity conditions: that the polynomials
    # along the edges of the Newton polygon are square-free.
    for e in P.bounded_edges():
        h = kT([D.get(tuple(c), 0) for c in Polyhedron(e).integral_points()])
        if not h.is_squarefree():
            return None
    x, y = f.parent().gens()
    return [
        x**(a[0] - 1) * y**(a[1] - 1)
        for a in P.integral_points()
        if P.interior_contains(a)
    ]


def find_closest_element(item, lst):
    r"""
    Return the index of the closest element of a list.

    Given ``List`` and ``item``, return the index of the element ``l`` of ``List``
    which minimises ``(item-l).abs()``. If there are multiple such elements, the
    first is returned.

    INPUT:

    - ``item`` -- value to minimize the distance to over the list

    - ``lst`` -- list to look for closest element in

    EXAMPLES::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import find_closest_element
        sage: i = 5
        sage: l = list(range(10))
        sage: i == find_closest_element(i, l)
        True

    Note that this method does no checks on the input, but will fail for inputs
    where the absolute value or subtraction do not make sense.
    """
    dists = [(item - l).abs() for l in lst]
    return dists.index(min(dists))


def reparameterize_differential_minpoly(minpoly, z0):
    r"""
    Rewrites a minimal polynomial to write is around `z_0`.

    Given a minimal polynomial `m(z,g)`, where `g` corresponds to a differential
    on the surface (that is, it is represented as a rational function, and
    implicitly carries a factor `dz`), we rewrite the minpoly in terms of
    variables `\bar{z}, \bar{g}` s.t now `\bar{z}=0 \Leftrightarrow z=z_0`.

    INPUT:

    - ``minpoly`` -- a polynomial in two variables, where the first variable
      corresponds to the base coordinate on the Riemann surface
    - ``z0`` -- complex number or infinity; the point about which to
      reparameterize

    OUTPUT: a polynomial in two variables giving the reparameterize minimal polynomial

    EXAMPLES:

    On the curve given by `w^2 - z^3 + 1 = 0`, we have differential
    `\frac{dz}{2w} = \frac{dz}{2\sqrt{z^3-1}}`
    with minimal polynomial `g^2(z^3-1) - 1/4=0`. We can make the substitution
    `\bar{z}=z^{-1}` to parameterise the differential about `z=\infty` as

    .. MATH::

        \frac{-\bar{z}^{-2} d\bar{z}}{2\sqrt{\bar{z}^{-3}-1}} = \frac{-d\bar{z}}{2\sqrt{\bar{z}(1-\bar{z}^3)}}.

    Hence the transformed differential should have minimal polynomial
    `\bar{g}^2 \bar{z} (1 - \bar{z}^3) - 1/4 = 0`, and we can check this::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, reparameterize_differential_minpoly
        sage: R.<z,w> = QQ[]
        sage: S = RiemannSurface(w^2-z^3+1)
        sage: minpoly = S._cohomology_basis_bounding_data[1][0][2]
        sage: z0 = Infinity
        sage: reparameterize_differential_minpoly(minpoly, z0)
        -zbar^4*gbar^2 + zbar*gbar^2 - 1/4

    We can further check that reparameterising about `0` is the identity
    operation::

        sage: reparameterize_differential_minpoly(minpoly, 0)(*minpoly.parent().gens()) == minpoly
        True

    .. NOTE::

        As part of the routine, when reparameterising about infinity, a
        rational function is reduced and then the numerator is taken. Over
        an inexact ring this is numerically unstable, and so it is advisable
        to only reparameterize about infinity over an exact ring.
    """
    P = minpoly.parent()
    F = PolynomialRing(P.base_ring(), [str(v) + "bar" for v in P.gens()])

    try:
        Inf = bool(z0 == z0.parent()(Infinity))
    except TypeError:
        Inf = False

    if Inf:
        F = F.fraction_field()
        mt = F(minpoly(F.gen(0)**(-1), -F.gen(0)**2 * F.gen(1)))
        mt.reduce()
        mt = mt.numerator()
    else:
        mt = minpoly(F.gen(0) + z0, F.gen(1))
    return mt


class RiemannSurface:
    r"""
    Construct a Riemann Surface. This is specified by the zeroes of a bivariate
    polynomial with rational coefficients `f(z,w) = 0`.

    INPUT:

    - ``f`` -- a bivariate polynomial with rational coefficients. The surface is
      interpreted as the covering space of the coordinate plane in the first
      variable.

    - ``prec`` -- the desired precision of computations on the surface in bits
      (default: 53)

    - ``certification`` -- boolean (default: ``True``); value indicating
      whether homotopy continuation is certified or not. Uncertified
      homotopy continuation can be faster.

    - ``differentials`` -- (default: ``None``) if specified, provides a list
      of polynomials `h` such that `h/(df/dw) dz` is a regular
      differential on the Riemann surface. This is taken as a basis of
      the regular differentials, so the genus is assumed to be equal
      to the length of this list. The results from the homology basis
      computation are checked against this value.  Providing this
      parameter makes the computation independent from Singular.  For
      a nonsingular plane curve of degree `d`, an appropriate set is
      given by the monomials of degree up to `d-3`.

    - ``integration_method`` -- (default: ``'rigorous'``). String specifying the
      integration method to use when calculating the integrals of differentials.
      The options are ``'heuristic'`` and ``'rigorous'``, the latter of
      which is often the most efficient.

    EXAMPLES::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
        sage: R.<z,w> = QQ[]
        sage: f = w^2 - z^3 + 1
        sage: RiemannSurface(f)
        Riemann surface defined by polynomial f = -z^3 + w^2 + 1 = 0, with 53 bits of precision

    Another Riemann surface with 100 bits of precision::

        sage: S = RiemannSurface(f, prec=100); S
        Riemann surface defined by polynomial f = -z^3 + w^2 + 1 = 0, with 100 bits of precision
        sage: S.riemann_matrix()^6 #abs tol 0.00000001
        [1.0000000000000000000000000000 - 1.1832913578315177081175928479e-30*I]

    We can also work with Riemann surfaces that are defined over fields with a
    complex embedding, but since the current interface for computing genus and
    regular differentials in Singular presently does not support extensions of
    `\QQ`, we need to specify a description of the differentials ourselves. We give
    an example of a CM elliptic curve::

        sage: Qt.<t> = QQ[]
        sage: K.<a> = NumberField(t^2-t+3,embedding=CC(0.5+1.6*I))
        sage: R.<x,y> = K[]
        sage: f = y^2 + y - (x^3 + (1-a)*x^2 - (2+a)*x - 2)
        sage: S = RiemannSurface(f, prec=100, differentials=[1])
        sage: A = S.endomorphism_basis()
        sage: len(A)
        2
        sage: all(len(T.minpoly().roots(K)) > 0 for T in A)
        True

    The ``'heuristic'`` integration method uses the method ``integrate_vector``
    defined in ``sage.numerical.gauss_legendre`` to compute integrals of differentials.
    As mentioned there, this works by iteratively doubling the number of nodes
    used in the quadrature, and uses a heuristic based on the rate at which the
    result is seemingly converging to estimate the error. The ``'rigorous'``
    method uses results from [Neu2018]_, and bounds the algebraic integrands on
    circular domains using Cauchy's form of the remainder in Taylor approximation
    coupled to Fujiwara's bound on polynomial roots (see Bruin-DisneyHogg-Gao,
    in preparation). Note this method of bounding on circular domains is also
    implemented in :meth:`_compute_delta`. The net result of this bounding is
    that one can know (an upper bound on) the number of nodes required to achieve
    a certain error. This means that for any given integral, assuming that the
    same number of nodes is required by both methods in order to achieve the
    desired error (not necessarily true in practice), approximately half
    the number of integrand evaluations are required. When the required number
    of nodes is high, e.g. when the precision required is high, this can make
    the ``'rigorous'`` method much faster. However, the ``'rigorous'`` method does
    not benefit as much from the caching of the ``nodes`` method over multiple
    integrals. The result of this is that, for calls of :meth:`matrix_of_integral_values`
    if the computation is 'fast', the heuristic method may outperform the
    rigorous method, but for slower computations the rigorous method can be much
    faster::

        sage: f = z*w^3 + z^3 + w
        sage: p = 53
        sage: Sh = RiemannSurface(f, prec=p, integration_method='heuristic')
        sage: Sr = RiemannSurface(f, prec=p, integration_method='rigorous')
        sage: from sage.numerical.gauss_legendre import nodes
        sage: import time
        sage: nodes.cache.clear()
        sage: ct = time.time()
        sage: Rh = Sh.riemann_matrix()
        sage: ct1 = time.time()-ct
        sage: nodes.cache.clear()
        sage: ct = time.time()
        sage: Rr = Sr.riemann_matrix()
        sage: ct2 = time.time()-ct
        sage: ct2/ct1  # random
        1.2429363969691192

    Note that for the above curve, the branch points are evenly distributed, and
    hence the implicit assumptions in the heuristic method are more sensible,
    meaning that a higher precision is required to see the heuristic method
    being significantly slower than the rigorous method. For a worse conditioned
    curve, this effect is more pronounced::

        sage: q = 1 / 10
        sage: f = y^2 - (x^2 - 2*x + 1 + q^2) * (x^2 + 2*x + 1 + q^2)
        sage: p = 500
        sage: Sh = RiemannSurface(f, prec=p, integration_method='heuristic')
        sage: Sr = RiemannSurface(f, prec=p, integration_method='rigorous')
        sage: nodes.cache.clear()
        sage: Rh = Sh.riemann_matrix()  # long time (8 seconds)
        sage: nodes.cache.clear()
        sage: Rr = Sr.riemann_matrix()  # long time (1 seconds)

    This disparity in timings can get increasingly worse, and testing has shown
    that even for random quadrics the heuristic method can be as bad as 30 times
    slower.

    TESTS:

    This elliptic curve has a relatively poorly conditioned set of branch
    points, so it challenges the path choice a bit. The code just verifies that
    the period is quadratic, because the curve has CM, but really the test is
    that the computation completes at all.::

        sage: prec = 50
        sage: Qx.<t> = QQ[]
        sage: CC = ComplexField(prec)
        sage: g = t^2-t-1
        sage: phiCC = g.roots(CC)[1][0]
        sage: K.<phi> = NumberField(g, embedding=phiCC)
        sage: R.<X,Y> = K[]
        sage: f = Y^2+X*Y+phi*Y-(X^3-X^2-2*phi*X+phi)
        sage: S = RiemannSurface(f,prec=prec, differentials=[1])
        sage: tau = S.riemann_matrix()[0, 0]
        sage: tau.algdep(6).degree() == 2
        True
    """

    def __init__(
        self,
        f,
        prec=53,
        certification=True,
        differentials=None,
        integration_method="rigorous"
    ):
        r"""
        TESTS::

            sage: R.<z,w> = QQ[]
            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: S = RiemannSurface(w^2 - z^3 + 1)
            sage: TestSuite(S).run() #not tested; Unclear what pickling strategy is best.
        """
        # Initializations.
        self._prec = prec
        self._certification = certification
        if not (integration_method == "heuristic" or integration_method == "rigorous"):
            raise ValueError("invalid integration method")
        self._integration_method = integration_method
        self._R = f.parent()
        if len(self._R.gens()) != 2:
            raise ValueError('only bivariate polynomials supported')
        if f.degree() <= 1:
            raise ValueError('equation must be of degree at least 2')
        z, w = self._R.gen(0), self._R.gen(1)
        self._CC = ComplexField(self._prec)
        self._RR = RealField(self._prec)
        self._CCz = PolynomialRing(self._CC, [self._R.gen(0)])
        self._CCw = PolynomialRing(self._CC, [self._R.gen(1)])
        self._RRz = PolynomialRing(self._RR, [self._R.gen(0)])
        self._curve = None
        self.f = f
        if differentials is not None:
            self._differentials = [self._R(a) for a in differentials]
            self.genus = len(self._differentials)
        else:
            B = differential_basis_baker(f)
            if B is not None:
                self._differentials = B
                self.genus = len(B)
            else:
                self._differentials = None
                self.genus = self._R.ideal(self.f).genus()
                if self.genus < 0:
                    raise ValueError(
                        "Singular reports negative genus. Specify differentials manually."
                    )
        self.degree = self.f.degree(w)
        self._dfdw = self.f.derivative(w)
        self._dfdz = self.f.derivative(z)
        self._discriminant = self.f.resultant(self._dfdw, w)
        # Coefficients of the polynomial for use in homotopy continuation.
        self._a0 = self._CCz(self.f.coefficient({w: self.degree})(self._CCz.gen(), 0))
        self._a0roots = self._a0.roots(multiplicities=False)
        self._aks = [
            self._CCz(self.f.coefficient({w: self.degree - k - 1})(self._CCz.gen(), 0))
            for k in range(self.degree)
        ]
        # Compute the branch locus. Takes the square-free part of the discriminant
        # because of numerical issues.
        self.branch_locus = []
        existing_factors = [x[0] for x in self._discriminant.factor()]
        for fac in existing_factors:
            self.branch_locus += self._CCz(fac(self._CCz.gen(), 0)).roots(
                multiplicities=False
            )
        self._f_branch_locus = self.branch_locus
        self._cohomology_basis_bounding_data = self._bounding_data(
            self.cohomology_basis(), exact=True
        )
        RBzg, bounding_data_list = self._cohomology_basis_bounding_data
        minpoly_list = [bd[2] for bd in bounding_data_list]
        # We now want to calculate the additional branchpoints associated to
        # the differentials.
        discriminants = [RBzg(self._discriminant(*RBzg.gens()))]
        for minpoly in minpoly_list:
            F = RBzg(minpoly)
            dF = F.derivative(RBzg.gen(1))
            discriminants += [F.resultant(dF, RBzg.gen(1))]
        combined_discriminant = lcm(discriminants)(*self._R.gens())
        self._differentials_branch_locus = []
        for x in combined_discriminant.factor():
            if x[0] not in existing_factors:
                self._differentials_branch_locus += self._CCz(
                    x[0](self._CCz.gen(), 0)
                ).roots(multiplicities=False)
        # We add these branchpoints to the existing.
        # self.branch_locus = self.branch_locus+self._differentials_branch_locus
        # We now want to also check whether Infinity is a branch point of any
        # of the differentials.
        # This will be useful when calculating the Abel-Jacobi map.
        minpoly_list = [
            reparameterize_differential_minpoly(mp, Infinity) for mp in minpoly_list
        ]
        discriminants = []
        for minpoly in minpoly_list:
            F = RBzg(minpoly)
            dF = F.derivative(RBzg.gen(1))
            discriminants += [F.resultant(dF, RBzg.gen(1))]
        discriminant_about_infinity = RBzg(lcm(discriminants))
        if discriminant_about_infinity(0, 0) == 0:
            self._differentials_branch_locus.append(self._CC(Infinity))

        # Voronoi diagram and the important points associated with it
        self.voronoi_diagram = Voronoi(voronoi_ghost(self.branch_locus, CC=self._CC))
        self._vertices = [self._CC(x0, y0) for x0, y0 in self.voronoi_diagram.vertices]
        self._wvalues = [self.w_values(z0) for z0 in self._vertices]
        # We arbitrarily, but sensibly, set the basepoint to be the rightmost vertex
        self._basepoint = (
            self._vertices.index(sorted(self._vertices, key=lambda z: z.real())[-1]),
            0,
        )
        self._Sn = SymmetricGroup(range(self.degree))
        self._L = {}
        self._integral_dict = {}
        self._fastcall_f = fast_callable(f, domain=self._CC)
        self._fastcall_dfdw = fast_callable(self._dfdw, domain=self._CC)
        self._fastcall_dfdz = fast_callable(self._dfdz, domain=self._CC)
        self._fastcall_cohomology_basis = [
            fast_callable(h, domain=self._CC) for h in self.cohomology_basis()
        ]

    def __repr__(self):
        r"""
        Return a string representation of the Riemann surface class.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: RiemannSurface(f)
            Riemann surface defined by polynomial f = -z^4 + w^2 + 1 = 0, with 53 bits of precision
        """
        s = (
            "Riemann surface defined by polynomial f = %s = 0, with %s bits of precision"
            % (self.f, self._prec)
        )
        return s

    def w_values(self, z0):
        r"""
        Return the points lying on the surface above ``z0``.

        INPUT:

        - ``z0`` -- complex number; a point in the complex z-plane

        OUTPUT:

        A set of complex numbers corresponding to solutions of `f(z_0,w) = 0`.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)

        Find the w-values above the origin, i.e. the solutions of `w^2 + 1 = 0`::

            sage: S.w_values(0)  # abs tol 1e-14
            [-1.00000000000000*I, 1.00000000000000*I]

        Note that typically the method returns a list of length ``self.degree``,
        but that at ramification points, this may no longer be true::

            sage: S.w_values(1)  # abs tol 1e-14
            [0.000000000000000]
        """
        return self.f(z0, self._CCw.gen(0)).roots(multiplicities=False)

    @cached_method
    def downstairs_edges(self):
        r"""
        Compute the edgeset of the Voronoi diagram.

        OUTPUT:

        A list of integer tuples corresponding to edges between vertices in the
        Voronoi diagram.

        EXAMPLES:

        Form a Riemann surface, one with a particularly simple branch locus::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 + z^3 - z^2
            sage: S = RiemannSurface(f)

        Compute the edges::

            sage: S.downstairs_edges()
            [(0, 1), (0, 5), (1, 4), (2, 3), (2, 4), (3, 5), (4, 5)]

        This now gives an edgeset which one could use to form a graph.

        .. NOTE::

            The numbering of the vertices is given by the Voronoi package.
        """
        # Because of how we constructed the Voronoi diagram, the first n points
        # correspond to the branch locus points.
        # The regions of these points are all of the edges which do not go off
        # to infinity, which are exactly the ones we want.
        n = len(self.branch_locus)
        desired_edges = [
            self.voronoi_diagram.regions[self.voronoi_diagram.point_region[i]]
            for i in range(n)
        ]
        # First construct the edges as a set because the regions will overlap
        # and we do not want to have two of the same edge.
        edges1 = set()
        for c in desired_edges:
            for j in range(len(c) - 1):
                edges1.add(frozenset((c[j], c[j + 1])))
            edges1.add(frozenset((c[0], c[-1])))
        # Then make it into a list and sort it.
        # The sorting is important - it will make computing the monodromy group
        # MUCH easier.
        # We orient all the edges so that we go from lower to higher
        # numbered vertex for the continuation.
        edges = [(i0, i1) if (i0 < i1) else (i1, i0) for (i0, i1) in edges1]
        edges.sort()
        return edges

    def downstairs_graph(self):
        r"""
        Return the Voronoi decomposition as a planar graph.

        The result of this routine can be useful to interpret the labelling of
        the vertices. See also :meth:`upstairs_graph`.

        OUTPUT: the Voronoi decomposition as a graph, with appropriate planar embedding

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: S.downstairs_graph()
            Graph on 11 vertices
        """
        G = Graph(self.downstairs_edges())
        G.set_pos(dict(enumerate(list(v) for v in self._vertices)))
        return G

    @cached_method
    def upstairs_graph(self):
        r"""
        Return the graph of the upstairs edges.

        This method can be useful for generating paths in the surface between points labelled
        by upstairs vertices, and verifying that a homology basis is likely computed correctly.
        See also :meth:`downstairs_graph`.

        OUTPUT:

        The homotopy-continued Voronoi decomposition as a graph, with appropriate 3D embedding.

        EXAMPLES::

            sage: R.<z,w> = QQ[]
            sage: S = Curve(w^2-z^4+1).riemann_surface()
            sage: G = S.upstairs_graph(); G
            Graph on 22 vertices
            sage: G.genus()
            1
            sage: G.is_connected()
            True
        """
        G = Graph(self.upstairs_edges(), immutable=True)
        G.set_pos(
            {
                (i, j): [
                    self._vertices[i].real(),
                    self._vertices[i].imag(),
                    self.w_values(self._vertices[i])[j].imag(),
                ]
                for i in range(len(self._vertices))
                for j in range(self.degree)
            },
            dim=3,
        )
        return G

    def _compute_delta(self, z1, epsilon, wvalues=None):
        r"""
        Compute a delta for homotopy continuation when moving along a path.

        INPUT:

        - ``z1`` -- complex number in the z-plane

        - ``epsilon`` -- real number which is the minimum distance between
          the w-values above ``z1``

        - ``wvalues`` -- list (default: ``None``); if specified, saves
          recomputation

        OUTPUT: a real number, which is a step size for moving along a path

        EXAMPLES:

        Form a Riemann Surface::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)

        Pick a point which lies on the Voronoi diagram, and compute an
        appropriate epsilon::

            sage: z1 = S._vertices[0]
            sage: currw = S.w_values(z1)
            sage: n = len(currw)
            sage: epsilon = min([abs(currw[i] - currw[n-j-1]) for i in range(n) for j in range(n-i-1)])/3
            sage: S._compute_delta(z1, epsilon)  # abs tol 1e-8
            0.152628501142363

        If the Riemann surface does not have certified homotopy continuation,
        then the delta will just be the minimum distance away from a branch
        point::

            sage: T = RiemannSurface(f, certification=False)
            sage: z1 = T._vertices[0]
            sage: currw = T.w_values(z1)
            sage: n = len(currw)
            sage: epsilon = min([abs(currw[i] - currw[n-j-1]) for i in range(n) for j in range(n-i-1)])/3
            sage: T._compute_delta(z1, epsilon)  # abs tol 1e-8
            0.381881307912987
        """
        if self._certification:
            if wvalues is None:
                wvalues = self.w_values(z1)
            # For computation of rho. Need the branch locus + roots of a0.
            badpoints = self._f_branch_locus + self._a0roots
            rho = min(abs(z1 - z) for z in badpoints) / 2
            Y = max(
                abs(self._fastcall_dfdz(z1, wi) / self._fastcall_dfdw(z1, wi))
                for wi in wvalues
            )

            # compute M
            upperbounds = [
                sum(ak[k] * (abs(z1) + rho)**k for k in range(ak.degree()))
                for ak in self._aks
            ]
            upperbounds.reverse()
            # If a0 is a constant polynomial, it is obviously bounded below.
            if not self._a0roots:
                lowerbound = self._CC(self._a0) / 2
            else:
                lowerbound = (
                    self._a0[self._a0.degree()]
                    * prod(abs((zk - z1) - rho) for zk in self._a0roots)
                    / 2
                )
            M = 2 * max(
                (upperbounds[k] / lowerbound).abs().nth_root(k + 1)
                for k in range(self.degree - 1)
            )
            return (
                rho
                * (
                    ((rho * Y - epsilon)**2 + 4 * epsilon * M).sqrt()
                    - (rho * Y + epsilon)
                )
                / (2 * M - 2 * rho * Y)
            )
        else:
            # Instead, we just compute the minimum distance between branch
            # points and the point in question.
            return min(abs(b - z1) for b in self._f_branch_locus) / 2

    def homotopy_continuation(self, edge):
        r"""
        Perform homotopy continuation along an edge of the Voronoi diagram using
        Newton iteration.

        INPUT:

        - ``edge`` -- tuple ``(z_start, z_end)`` indicating the straight line
          over which to perform the homotopy continutation

        OUTPUT:

        A list containing the initialised continuation data. Each entry in the
        list contains: the `t` values that entry corresponds to, a list of
        complex numbers corresponding to the points which are reached when
        continued along the edge when traversing along the direction of the
        edge, and a value ``epsilon`` giving the minimumdistance between the
        fibre values divided by 3. The ordering of these points indicates how
        they have been permuted due to the weaving of the curve.

        EXAMPLES:

        We check that continued values along an edge correspond (up to the
        appropriate permutation) to what is stored. Note that the permutation
        was originally computed from this data::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f)
            sage: edge1 = sorted(S.edge_permutations())[0]
            sage: sigma = S.edge_permutations()[edge1]
            sage: edge = [S._vertices[i] for i in edge1]
            sage: continued_values = S.homotopy_continuation(edge)[-1][1]
            sage: stored_values = S.w_values(S._vertices[edge1[1]])
            sage: all(abs(continued_values[i]-stored_values[sigma(i)]) < 1e-8 for i in range(3))
            True
        """
        z_start, z_end = edge
        z_start = self._CC(z_start)
        z_end = self._CC(z_end)
        ZERO = self._RR.zero()
        ONE = self._RR.one()
        datastorage = []
        path_length = abs(z_end - z_start)

        def path(t):
            return z_start * (1 - t) + z_end * t

        # Primary procedure.
        T = ZERO
        currw = self.w_values(path(T))
        n = len(currw)
        epsilon = (
            min([abs(currw[i] - currw[j]) for i in range(1, n) for j in range(i)]) / 3
        )
        datastorage += [(T, currw, epsilon)]
        while T < ONE:
            delta = self._compute_delta(path(T), epsilon, wvalues=currw) / path_length
            # Move along the path by delta.
            T += delta
            # If T exceeds 1, just set it to 1 and compute.
            if T > ONE:
                delta -= T - ONE
                T = ONE
            while True:
                try:
                    neww = self._determine_new_w(path(T), currw, epsilon)
                except ConvergenceError:
                    delta /= 2
                    T -= delta
                else:
                    break
            currw = neww
            epsilon = (
                min([abs(currw[i] - currw[j]) for i in range(1, n) for j in range(i)])
                / 3
            )
            datastorage += [(T, currw, epsilon)]
        return datastorage

    def _determine_new_w(self, z0, oldw, epsilon):
        r"""
        A procedure to Newton iterate a list of w-values simultaneously.

        Used primarily for moving along the surface for integration or
        homotopy continuation.

        INPUT:

        - ``z0`` -- complex number

        - ``oldw`` -- list of w-values which are presumed to be guesses of
          the w-values above ``z0``

        - ``epsilon`` -- the minimum distance between the points of ``oldw``
          divided by 3

        OUTPUT:

        A list of points the same length as ``oldw`` corresponding to the new
        Newton iterated points.

        However, if the Newton iteration exceeds the allotted attempts,
        or exits the ``epsilon`` ball, raises a convergence error.

        EXAMPLES:

        First, a trivial example where we guess exactly what the roots are::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: z0 = S._vertices[0]
            sage: epsilon = 0.1
            sage: oldw = S.w_values(z0)
            sage: neww = S._determine_new_w(z0, oldw, epsilon); neww #abs tol 0.00000001
                [-0.934613146929672 + 2.01088055918363*I,
                 0.934613146929672 - 2.01088055918363*I]

        Which should be exactly the same as the w-values we started with.::

            sage: abs(neww[0] - oldw[0]) #abs tol 0.00000001
            0.000000000000...
            sage: abs(neww[1] - oldw[1]) #abs tol 0.00000001
            0.000000000000...

        Here is an example where we exit the ``epsilon`` bound. This approach is
        based on the homotopy continuation procedure which traverses along a
        path and attempts Newton iteration::

            sage: g = z^3*w + w^3 + z
            sage: T = RiemannSurface(g)
            sage: z0 = T._vertices[2]*(0.9) + 0.3*I
            sage: epsilon = 0.5
            sage: oldw = T.w_values(T._vertices[2])
            sage: T._determine_new_w(z0, oldw, epsilon)
            Traceback (most recent call last):
            ...
            ConvergenceError: Newton iteration escaped neighbourhood

        .. NOTE::

            Algorithmically, this method is nearly identical to :meth:`_newton_iteration`,
            but this method takes a list of `w` values.
        """
        # Tools of Newton iteration.
        F = self._fastcall_f
        dF = self._fastcall_dfdw
        neww = []
        prec = self._CC.prec()
        # Iterate over all roots.
        for i in range(len(oldw)):
            delta = F(z0, oldw[i]) / dF(z0, oldw[i])
            Ndelta = delta.norm()
            wi = oldw[i] - delta
            # it is possible in theory that Newton iteration fails to
            # converge without escaping. We catch this by capping the
            # number of iterations by 100
            for j in range(100):
                # If we exceed the epsilon bound from homotopy continuation,
                # terminate.
                if abs(wi - oldw[i]) >= epsilon:
                    raise ConvergenceError("Newton iteration escaped neighbourhood")
                new_delta = F(z0, wi) / dF(z0, wi)
                Nnew_delta = new_delta.norm()
                # If we found the root exactly, or if delta only affects half the digits and
                # stops getting smaller, we decide that we have converged.
                if (new_delta == 0) or (
                    Nnew_delta >= Ndelta
                    and Ndelta.sign_mantissa_exponent()[2] + prec
                    < wi.norm().sign_mantissa_exponent()[2]
                ):
                    neww.append(wi)
                    break
                delta = new_delta
                Ndelta = Nnew_delta
                wi -= delta
            # If we run 100 iterations without a result, terminate.
            else:
                raise ConvergenceError(
                    "Newton iteration fails to converge after %s iterations" % j
                )
        return neww

    def _newton_iteration(self, z0, oldw, epsilon):
        r"""
        A non-vectorized Newton iteration procedure used for integration.

        INPUT:

        - ``z0`` -- complex number

        - ``oldw`` -- a w-value which is presumed to be a guess of one of
          the w-values above ``z0``

        - ``epsilon`` -- the minimum distance between the w-values divided by 3

        OUTPUT:

        A complex number, which should be a w-value above ``z0``.

        However, if the Newton iteration exceeds the allotted attempts,
        or exits the ``epsilon`` ball, raises a convergence error.

        EXAMPLES:

        First, a trivial example where we guess exactly what the root is::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: z0 = S._vertices[0]
            sage: epsilon = 0.1
            sage: oldw = S.w_values(z0)[0]
            sage: neww = S._newton_iteration(z0,oldw,epsilon); neww  # abs tol 0.00000001
            -0.934613146929672 + 2.01088055918363*I

        Which should be exactly the same as the w-value we started with::

            sage: oldw - neww  # abs tol 0.00000001
            0.000000000000000

        Here is an example where we exit the epsilon bound. This approach is
        based on the homotopy continuation procedure which traverses along a
        path and attempts Newton iteration::

            sage: g = z^3*w + w^3 + z
            sage: T = RiemannSurface(g)
            sage: z0 = T._vertices[2]*(0.9) + 0.3*I
            sage: epsilon = 0.5
            sage: oldw = T.w_values(T._vertices[2])[1]
            sage: T._newton_iteration(z0, oldw, epsilon)
            Traceback (most recent call last):
            ...
            ConvergenceError: Newton iteration escaped neighbourhood
        """
        F = self._fastcall_f
        dF = self._fastcall_dfdw
        prec = self._CC.prec()
        delta = F(z0, oldw) / dF(z0, oldw)
        Ndelta = delta.norm()
        neww = oldw - delta
        eps_squared = epsilon**2
        # it is possible in theory that Newton iteration fails to converge
        # without escaping. We catch this by capping the number of iterations
        # by 100
        for j in range(100):
            if (neww - oldw).norm() > eps_squared:
                raise ConvergenceError("Newton iteration escaped neighbourhood")
            new_delta = F(z0, neww) / dF(z0, neww)
            Nnew_delta = new_delta.norm()
            # If we found the root exactly, or if delta only affects half the digits and
            # stops getting smaller, we decide that we have converged.
            if (new_delta == 0) or (
                Nnew_delta >= Ndelta
                and Ndelta.sign_mantissa_exponent()[2] + prec
                < neww.norm().sign_mantissa_exponent()[2]
            ):
                return neww
            delta = new_delta
            Ndelta = Nnew_delta
            neww -= delta
        raise ConvergenceError("Newton iteration fails to converge")

    @cached_method
    def upstairs_edges(self):
        r"""
        Compute the edgeset of the lift of the downstairs graph onto the Riemann
        surface.

        OUTPUT:

        An edgeset between vertices (i, j), where i corresponds to the i-th
        point in the Voronoi diagram vertices, and j is the j-th w-value
        associated with that point.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 + z^3 - z^2
            sage: S = RiemannSurface(f)
            sage: edgeset = S.upstairs_edges()
            sage: len(edgeset) == S.degree*len(S.downstairs_edges())
            True
            sage: {(v[0],w[0]) for v,w in edgeset} == set(S.downstairs_edges())
            True
        """
        edgeset = []
        n = len(self._wvalues[0])
        # Lifts each edge individually.
        for e in self.downstairs_edges():
            i0, i1 = e
            d_edge = (self._vertices[i0], self._vertices[i1])
            # Epsilon for checking w-value later.
            val = self._wvalues[i1]
            epsilon = (
                min(
                    abs(val[i] - val[n - j - 1])
                    for i in range(n)
                    for j in range(n - i - 1)
                )
                / 3
            )
            # Homotopy continuation along e.
            self._L[e] = self.homotopy_continuation(d_edge)
            homotopycont = self._L[e][-1][1]
            for i in range(len(homotopycont)):
                # Checks over the w-values of the next point to check which it is.
                for j in range(len(self._wvalues[i1])):
                    if abs(homotopycont[i] - self._wvalues[i1][j]) < epsilon:
                        # Once it finds the appropriate w-value, adds the edge.
                        edgeset = edgeset + [[(i0, i), (i1, j)]]
                        continue
        return edgeset

    def _edge_permutation(self, edge):
        r"""
        Compute the permutation of the w-values above a point in the z-plane
        when moving along an edge in the Voronoi diagram.

        INPUT:

        - ``edge`` -- an edge on the Voronoi diagram

        OUTPUT:

        A permutation corresponding to how the roots interchange when moving
        along the edge.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f)

        Compute the edge permutation of (1, 2) on the Voronoi diagram::

            sage: S._edge_permutation((1, 2))
            (0,2,1)

        This indicates that while traversing along the direction of `(2, 9)`,
        the 2nd and 3rd layers of the Riemann surface are interchanging.
        """
        if edge in self.downstairs_edges():
            # find all upstairs edges that are lifts of the given
            # downstairs edge and store the corresponding indices at
            # start and end that label the branches upstairs.
            L = [
                (j0, j1)
                for ((i0, j0), (i1, j1)) in self.upstairs_edges()
                if edge == (i0, i1)
            ]
            # we should be finding exactly "degree" of these
            assert len(L) == self.degree
            # and as a corollary of how we construct them, the indices
            # at the start should be in order
            assert all(a == b[0] for a, b in enumerate(L))
            return self._Sn([j1 for j0, j1 in L])
        raise ValueError("edge not in Voronoi diagram")

    @cached_method
    def edge_permutations(self) -> dict:
        r"""
        Compute the permutations of branches associated to each edge.

        Over the vertices of the Voronoi decomposition around the branch locus,
        we label the fibres. By following along an edge, the lifts of the edge
        induce a permutation of that labelling.

        OUTPUT:

        A dictionary with as keys the edges of the Voronoi decomposition and as
        values the corresponding permutations.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 + z^2+1
            sage: S = RiemannSurface(f)
            sage: S.edge_permutations()
            {(0, 2): (),
             (0, 4): (),
             (1, 2): (),
             (1, 3): (0,1),
             (1, 6): (),
             (2, 0): (),
             (2, 1): (),
             (2, 5): (0,1),
             (3, 1): (0,1),
             (3, 4): (),
             (4, 0): (),
             (4, 3): (),
             (5, 2): (0,1),
             (5, 7): (),
             (6, 1): (),
             (6, 7): (),
             (7, 5): (),
             (7, 6): ()}
        """
        D = {e: self._edge_permutation(e) for e in self.downstairs_edges()}
        for (a, b), p in list(D.items()):
            D[(b, a)] = p**(-1)
        return D

    @cached_method
    def monodromy_group(self):
        r"""
        Compute local monodromy generators of the Riemann surface.

        For each branch point, the local monodromy is encoded by a permutation.
        The permutations returned correspond to positively oriented loops around
        each branch point, with a fixed base point. This means the generators
        are properly conjugated to ensure that together they generate the global
        monodromy. The list has an entry for every finite point stored in
        ``self.branch_locus``, plus an entry for the ramification above infinity.

        OUTPUT:

        A list of permutations, encoding the local monodromy at each branch
        point.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z, w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f)
            sage: G = S.monodromy_group(); G
            [(0,1,2), (0,1), (0,2), (1,2), (1,2), (1,2), (0,1), (0,2), (0,2)]

        The permutations give the local monodromy generators for the branch
        points::

            sage: list(zip(S.branch_locus + [unsigned_infinity], G)) #abs tol 0.0000001
            [(0.000000000000000, (0,1,2)),
             (-1.31362670141929, (0,1)),
             (-0.819032851784253 - 1.02703471138023*I, (0,2)),
             (-0.819032851784253 + 1.02703471138023*I, (1,2)),
             (0.292309440469772 - 1.28069133740100*I, (1,2)),
             (0.292309440469772 + 1.28069133740100*I, (1,2)),
             (1.18353676202412 - 0.569961265016465*I, (0,1)),
             (1.18353676202412 + 0.569961265016465*I, (0,2)),
             (Infinity, (0,2))]

        We can check the ramification by looking at the cycle lengths and verify
        it agrees with the Riemann-Hurwitz formula::

            sage: 2*S.genus-2 == -2*S.degree + sum(e-1 for g in G for e in g.cycle_type())
            True
        """
        n = len(self.branch_locus)
        G = Graph(self.downstairs_edges())
        # we get all the regions
        loops = [
            self.voronoi_diagram.regions[i][:]
            for i in self.voronoi_diagram.point_region
        ]
        # and construct their Voronoi centers as complex numbers
        centers = self.branch_locus + [
            self._CC(x, y) for x, y in self.voronoi_diagram.points[n:]
        ]
        for center, loop in zip(centers, loops):
            if -1 in loop:
                # for loops involving infinity we take the finite part of the path
                i = loop.index(-1)
                loop[:] = loop[i + 1 :] + loop[:i]
            else:
                # and for finite ones we close the paths
                loop.append(loop[0])
            # we make sure the loops are positively oriented wrt. their center
            v0 = self._vertices[loop[0]]
            v1 = self._vertices[loop[1]]
            M = Matrix([list(v0 - center), list(v1 - center)])
            if M.det() < 0:
                loop.reverse()

        # we stitch together the paths that are part of loops through
        # infinity. There should be a unique way of doing so.
        inf_loops = loops[n:]
        inf_path = inf_loops.pop()
        while inf_loops:
            inf_path += (inf_loops.pop())[1:]
        assert inf_path[0] == inf_path[-1]

        loops = loops[:n]
        loops.append(inf_path)

        P0 = loops[0][0]
        monodromy_gens = []
        edge_perms = self.edge_permutations()
        SG = self._Sn
        for c in loops:
            to_loop = G.shortest_path(P0, c[0])
            to_loop_perm = SG.prod(
                edge_perms[(to_loop[i], to_loop[i + 1])]
                for i in range(len(to_loop) - 1)
            )
            c_perm = SG.prod(edge_perms[(c[i], c[i + 1])] for i in range(len(c) - 1))
            monodromy_gens.append(to_loop_perm * c_perm * ~to_loop_perm)
        return monodromy_gens

    @cached_method
    def homology_basis(self):
        r"""
        Compute the homology basis of the Riemann surface.

        OUTPUT:

        A list of paths `L = [P_1, \dots, P_n]`. Each path `P_i` is of the form
        `(k, [p_1 ... p_m, p_1])`, where `k` is the number of times to traverse
        the path (if negative, to traverse it backwards), and the `p_i` are
        vertices of the upstairs graph.

        EXAMPLES:

        In this example, there are two paths that form the homology basis::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: g = w^2 - z^4 + 1
            sage: S = RiemannSurface(g)
            sage: S.homology_basis()  # random
            [[(1, [(3, 1), (5, 0), (9, 0), (10, 0), (2, 0), (4, 0),
                (7, 1), (10, 1), (3, 1)])],
             [(1, [(8, 0), (6, 0), (7, 0), (10, 0), (2, 0), (4, 0),
                (7, 1), (10, 1), (9, 1), (8, 0)])]]

        In order to check that the answer returned above is reasonable, we
        test some basic properties. We express the faces of the downstairs graph
        as ZZ-linear combinations of the edges and check that the projection
        of the homology basis upstairs projects down to independent linear
        combinations of an even number of faces::

            sage: dg = S.downstairs_graph()
            sage: edges = dg.edges(sort=True)
            sage: E = ZZ^len(edges)
            sage: edge_to_E = { e[:2]: E.gen(i) for i,e in enumerate(edges)}
            sage: edge_to_E.update({ (e[1],e[0]): -E.gen(i) for i,e in enumerate(edges)})
            sage: face_span = E.submodule([sum(edge_to_E[e] for e in f) for f in dg.faces()])
            sage: def path_to_E(path):
            ....:     k,P = path
            ....:     return k*sum(edge_to_E[(P[i][0],P[i+1][0])] for i in range(len(P)-1))
            sage: hom_basis = [sum(path_to_E(p) for p in loop) for loop in S.homology_basis()]
            sage: face_span.submodule(hom_basis).rank()
            2
            sage: [sum(face_span.coordinate_vector(b))%2 for b in hom_basis]
            [0, 0]
        """
        if self.genus == 0:
            return []

        cycles = self.upstairs_graph().cycle_basis()
        # Computing the Gram matrix.
        cn = len(cycles)
        # Forming a list of lists of zeroes.
        # Later this will be converted into a matrix.
        intersectionprod = [[0] * cn for _ in cycles]

        # as it turns out, in extreme examples argument computation
        # can be quite dominant so we cache this (since we may end up
        # using these values multiple times)
        direction_cache = {}

        def direction(center, neighbour):
            k = (center, neighbour)
            if k not in direction_cache:
                theta = (self._vertices[neighbour] - self._vertices[center]).argument()
                direction_cache[k] = theta
                return theta
            else:
                return direction_cache[k]

        # This loop will start at the entry (0,1), and proceed along the row up
        # til (0,cn-1).
        # Then it will go to entry (1,2), and proceed along the row, etc.
        for i in range(1, cn):
            for j in range(i):
                # Initializing the intersection product value.
                intsum = 0
                # Intersection of the edges
                intsec = set(cycles[i]).intersection(set(cycles[j]))
                for v in intsec:
                    # Get indices of the vertex in the cycles.
                    i0 = cycles[i].index(v)
                    i1 = cycles[j].index(v)
                    # Get the complex value of the vertex v.
                    center = cycles[i][i0][0]

                    # We are in the following situation:
                    # We have two paths a_in->v->a_out and
                    # b_in->v->b_out intersecting. We say they
                    # are "positively oriented" if the a-path
                    # and the b-path are oriented as the x and y axes, i.e.,
                    # if, when we walk around v in counter-clockwise direction,
                    # we encounter a_in,b_in,a_out,b_out.

                    # we can also have that b_in and/or b_out overlaps with
                    # a_in and/or a_out. If we just score the orientation of
                    # b_in and b_out individually, we can deal with this
                    # by just ignoring the overlapping vertex. The "half"
                    # score will be appropriately complemented at one of the
                    # next vertices.

                    a_in = cycles[i][i0 - 1][0]
                    a_out = cycles[i][(i0 + 1) % len(cycles[i])][0]
                    b_in = cycles[j][i1 - 1][0]
                    b_out = cycles[j][(i1 + 1) % len(cycles[j])][0]

                    # we can get the angles (and hence the rotation order)
                    # by taking the arguments of the differences.

                    a_in_arg = direction(center, a_in)
                    a_out_arg = direction(center, a_out)
                    b_in_arg = direction(center, b_in)
                    b_out_arg = direction(center, b_out)

                    # we make sure to test overlap on the indices, so no rounding
                    # problems occur with that.

                    if (b_in != a_in) and (b_in != a_out):
                        if (
                            (a_in_arg < b_in_arg < a_out_arg)
                            or (b_in_arg < a_out_arg < a_in_arg)
                            or (a_out_arg < a_in_arg < b_in_arg)
                        ):
                            intsum += 1
                        elif (
                            (a_out_arg < b_in_arg < a_in_arg)
                            or (b_in_arg < a_in_arg < a_out_arg)
                            or (a_in_arg < a_out_arg < b_in_arg)
                        ):
                            intsum -= 1
                        else:
                            raise RuntimeError("impossible edge orientation")
                    if (b_out != a_in) and (b_out != a_out):
                        if (
                            (a_in_arg < b_out_arg < a_out_arg)
                            or (b_out_arg < a_out_arg < a_in_arg)
                            or (a_out_arg < a_in_arg < b_out_arg)
                        ):
                            intsum -= 1
                        elif (
                            (a_out_arg < b_out_arg < a_in_arg)
                            or (b_out_arg < a_in_arg < a_out_arg)
                            or (a_in_arg < a_out_arg < b_out_arg)
                        ):
                            intsum += 1
                        else:
                            raise RuntimeError("impossible edge orientation")
                assert (intsum % 2) == 0
                intsum = intsum // 2
                intersectionprod[i][j] = intsum
                # Skew Symmetry
                intersectionprod[j][i] = -intsum
        Gmatrix = Matrix(intersectionprod)
        G_normalized, P = Gmatrix.symplectic_form()
        if G_normalized.rank() != 2 * self.genus:
            raise RuntimeError("rank of homology pairing mismatches twice stored genus")
        # Define the cycle sets.
        acycles = [[] for i in range(self.genus)]
        bcycles = [[] for i in range(self.genus)]
        # There are g a and b cycles.
        for i in range(self.genus):
            # Range over the size of the Gram matrix.
            for j in range(cn):
                # Forms the acycles and bcycles. If the entry in the
                # transformation matrix is nonzero, it adds the coefficient at
                # that entry, and the corresponding cycle. (also, forms it
                # into a loop)
                if P[i][j] != 0:
                    acycles[i] += [(P[i][j], list(cycles[j]) + [cycles[j][0]])]
                if P[self.genus + i][j] != 0:
                    bcycles[i] += [
                        (P[self.genus + i][j], list(cycles[j]) + [cycles[j][0]])
                    ]
        return acycles + bcycles

    def make_zw_interpolator(self, upstairs_edge, initial_continuation=None):
        r"""
        Given a downstairs edge for which continuation data has been initialised,
        return a function that computes `z(t), w(t)` , where `t` in `[0,1]` is a
        parametrization of the edge.

        INPUT:

        - ``upstairs_edge`` -- tuple ``((z_start, sb), (z_end,))`` giving the
          start and end values of the base coordinate along the straight-line
          path and the starting branch
        - ``initial_continuation`` -- list (optional);  output of
          ``homotopy_continuation`` initialising the continuation data

        OUTPUT:

        A tuple ``(g, d)``, where ``g`` is the function that computes the interpolation
        along the edge and ``d`` is the difference of the z-values of the end and
        start point.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: _ = S.homology_basis()
            sage: u_edge = [(0, 0), (1, 0)]
            sage: d_edge = tuple(u[0] for u in u_edge)
            sage: u_edge = [(S._vertices[i], j) for i, j in u_edge]
            sage: initial_continuation = S._L[d_edge]
            sage: g, d = S.make_zw_interpolator(u_edge, initial_continuation)
            sage: all(f(*g(i*0.1)).abs() < 1e-13 for i in range(10))
            True
            sage: abs((g(1)[0]-g(0)[0]) - d) < 1e-13
            True

        .. NOTE::

            The interpolator returned by this method can effectively hang if
            either ``z_start`` or ``z_end`` are branchpoints. In these situations
            it is better to take a different approach rather than continue to use
            the interpolator.
        """
        downstairs_edge = tuple(u[0] for u in upstairs_edge)
        z_start, z_end = downstairs_edge
        z_start = self._CC(z_start)
        z_end = self._CC(z_end)
        if initial_continuation is None:
            initial_continuation = self.homotopy_continuation(downstairs_edge)
        currL = initial_continuation
        windex = upstairs_edge[0][1]

        def w_interpolate(t):
            if t < 0 or t > 1:
                raise ValueError("t outside path range")
            if t == 0:
                return z_start, currL[0][1][windex]
            elif t == 1:
                return z_end, currL[-1][1][windex]
            while True:
                i = bisect(currL, t)
                t1, w1, epsilon = currL[i]
                w1 = w1[windex]
                t2, w2, _ = currL[i + 1]
                w2 = w2[windex]
                z0 = (1 - t) * z_start + t * z_end
                w0 = self._CC(((t2 - t) * w1 + (t - t1) * w2) / (t2 - t1))
                try:
                    desired_result = self._newton_iteration(z0, w0, epsilon)
                except ConvergenceError:
                    pass
                else:
                    return z0, desired_result
                # If we did not succeed, we insert a new point in our interpolation list
                tnew = t
                while True:
                    tnew = (t1 + tnew) / 2
                    znew = (1 - tnew) * z_start + tnew * z_end
                    try:
                        neww1 = self._determine_new_w(znew, currL[i][1], epsilon)
                    except ConvergenceError:
                        pass
                    else:
                        # When *no* ConvergenceError is raised, we have succeeded and we can exit
                        break
                # once the loop has succeeded we insert our new value
                t1 = tnew
                currL.insert(i + 1, (t1, neww1, epsilon))

        return w_interpolate, (z_end - z_start)

    def simple_vector_line_integral(self, upstairs_edge, differentials):
        r"""
        Perform vectorized integration along a straight path.

        INPUT:

        - ``upstairs_edge`` -- tuple; either a pair of integer tuples
          corresponding to an edge of the upstairs graph, or a tuple
          ``((z_start, sb), (z_end, ))`` as in the input of
          ``make_zw_interpolator``

        - ``differentials`` -- list of polynomials; a polynomial `g`
          represents the differential `g(z,w)/(df/dw) dz` where `f(z,w)=0` is
          the equation defining the Riemann surface

        OUTPUT: a complex number, the value of the line integral

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f); S
            Riemann surface defined by polynomial f = -z^4 + w^2 + 1 = 0, with 53 bits of precision

        Since we make use of data from homotopy continuation, we need to compute
        the necessary data::

            sage: M = S.riemann_matrix()
            sage: differentials = S.cohomology_basis()
            sage: S.simple_vector_line_integral([(0, 0), (1, 0)], differentials) #abs tol 0.00000001
            (1.14590610929717e-16 - 0.352971844594760*I)

        .. NOTE::

            Uses data that :meth:`homology_basis` initializes, and may give incorrect
            values if :meth:`homology_basis` has not initialized them. In practice
            it is more efficient to set ``differentials`` to a fast-callable version
            of differentials to speed up execution.
        """
        d_edge = tuple(u[0] for u in upstairs_edge)
        # Using a try-catch here allows us to retain a certain amount of back compatibility
        # for users.
        try:
            initial_continuation = self._L[d_edge]
            upstairs_edge = (
                (self._vertices[d_edge[0]], upstairs_edge[0][1]),
                (self._vertices[d_edge[1]],),
            )
        except KeyError:
            initial_continuation = self.homotopy_continuation(d_edge)
        w_of_t, Delta_z = self.make_zw_interpolator(upstairs_edge, initial_continuation)
        V = VectorSpace(self._CC, len(differentials))

        def integrand(t):
            zt, wt = w_of_t(t)
            dfdwt = self._fastcall_dfdw(zt, wt)
            return V([omega(zt, wt) / dfdwt for omega in differentials])

        return integrate_vector(integrand, self._prec) * Delta_z

    def cohomology_basis(self, option=1):
        r"""
        Compute the cohomology basis of this surface.

        INPUT:

        - ``option`` -- presently, this routine uses Singular's ``adjointIdeal``
          and passes the ``option`` parameter on. Legal values are 1, 2, 3 ,4,
          where 1 is the default. See the Singular documentation for the
          meaning. The backend for this function may change, and support for
          this parameter may disappear.

        OUTPUT:

        This returns a list of polynomials `g` representing the holomorphic
        differentials `g/(df/dw) dz`, where `f(z,w)=0` is the equation
        specifying the Riemann surface.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f)
            sage: S.cohomology_basis()
            [1, w, z]
        """
        if self.genus == 0:
            self._differentials = []
            return self._differentials  # [0]
        if self._differentials is None:
            # Computes differentials from the adjointIdeal using Singular
            # First we homogenize
            base = self.f.base_ring()
            # It's important we use a degree ordering; see below.
            R = self._R
            k = PolynomialRing(base, names='Z,W,U', order='degrevlex')
            dehom = k.Hom(R)([R.gen(0), R.gen(1), R.one()])
            fnew = self.f(k.gen(0) / k.gen(2), k.gen(1) / k.gen(2)).numerator()

            # We load the relevant functionality into singularlib
            import sage.libs.singular.function_factory

            sage.libs.singular.function_factory.lib("paraplanecurves.lib")
            adjointIdeal = sage.libs.singular.function.singular_function("adjointIdeal")
            libsing_options = sage.libs.singular.option.LibSingularVerboseOptions()

            # We compute the adjoint ideal (note we need to silence "redefine")
            redef_save = libsing_options["redefine"]
            try:
                libsing_options["redefine"] = False
                J = adjointIdeal(fnew, option)
            finally:
                libsing_options["redefine"] = redef_save

            # We are interested in the (degree-3) subspace of the adjoint ideal.
            # We compute this by intersecting with (Z,W,U)^(degree-3). Then the
            # lowest degree generators are a basis of the relevant subspace.
            d = fnew.total_degree()
            J2 = k.ideal(J).intersection(
                k.ideal([k.gen(0), k.gen(1), k.gen(2)])**(d - 3)
            )
            generators = [dehom(c) for c in J2.gens() if c.degree() == d - 3]
            if len(generators) != self.genus:
                raise ValueError(
                    "computed regular differentials do not match stored genus"
                )
            self._differentials = generators
        return self._differentials

    def _bounding_data(self, differentials, exact=False):
        r"""
        Compute the data required to bound a differential on a circle.

        Given a differential, one can bound it on a circular region using its
        derivative and its minimal polynomial (in the coordinate of the base).

        INPUT:

        - ``differentials`` -- list of polynomials in ``self._R`` giving
          the numerators of the differentials, as per the output of
          :meth:`cohomology_basis`

        - ``exact`` -- boolean (default: ``False``); whether to return the minimal
          polynomials over the exact base ring, or over ``self._CC``

        OUTPUT:

        A tuple ``(Rzg, [(g, dgdz, F, a0_info), ...])`` where each element of
        the list corresponds to an element of ``differentials``. Introducing the
        notation ``RBzg = PolynomialRing(self._R, ['z','g'])`` and
        ``CCzg = PolynomialRing(self._CC, ['z','g'])``, we have that:

        - ``Rzg`` is either ``RBzg`` or ``CCzg`` depending on the value of
          ``exact``,
        - ``g`` is the full rational function in ``self._R.fraction_field()``
          giving the differential,
        - ``dgdz`` is the derivative of ``g`` with respect to ``self._R.gen(0)``,
          written in terms of ``self._R.gen(0)`` and ``g``, hence laying in
          ``RBzg``,
        - ``F`` is the minimal polynomial of ``g`` over ``self._R.gen(0)``,
          laying in the polynomial ring ``Rzg``,
        - ``a0_info`` is a tuple ``(lc, roots)`` where ``lc`` and ``roots`` are
          the leading coefficient and roots of the polynomial in ``CCzg.gen(0)``
          that is the coefficient of the term of ``F`` of highest degree in
          ``CCzg.gen(1)``.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: f = y^2 - x^3 + 1
            sage: S = RiemannSurface(f)
            sage: differentials = S.cohomology_basis(); differentials
            [1]
            sage: S._dfdw
            2*y
            sage: S._bounding_data(differentials)
            (Multivariate Polynomial Ring in z, g over Complex Field with 53 bits of precision,
             [(1/(2*y),
               (-3*z^2*g)/(2*z^3 - 2),
               z^3*g^2 - g^2 - 0.250000000000000,
               (1.00000000000000,
                [1.00000000000000,
                 -0.500000000000000 - 0.866025403784439*I,
                 -0.500000000000000 + 0.866025403784439*I]))])

        Note that ``self._bounding_data(self._cohomology_basis(), exact=True)``
        is stored in ``self._cohomology_basis_bounding_data``::

            sage: S._cohomology_basis_bounding_data
            (Multivariate Polynomial Ring in z, g over Rational Field,
             [(1/(2*y),
               (-3*z^2*g)/(2*z^3 - 2),
               z^3*g^2 - g^2 - 1/4,
               (1.00000000000000,
                [1.00000000000000,
                 -0.500000000000000 - 0.866025403784439*I,
                 -0.500000000000000 + 0.866025403784439*I]))])
        """
        # This copies previous work by NB, outputting the zipped list required
        # for a certified line integral.
        RB = self._R.base_ring()
        P = PolynomialRing(RB, "Z")
        k = P.fraction_field()
        KP = PolynomialRing(k, "W")  # W->fraction field
        fZW = self.f(P.gen(0), KP.gen(0))
        L = k.extension(fZW, "Wb")
        dfdw_L = self._dfdw(P.gen(0), L.gen(0))
        integrand_list = [h / self._dfdw for h in differentials]
        # minpoly_univ gives the minimal polynomial for h, in variable x, with
        # coefficients given by polynomials with coefficients in P (i.e.
        # rational polynomials in Z).
        minpoly_univ = [
            (h(P.gen(0), L.gen(0)) / dfdw_L).minpoly().numerator()
            for h in differentials
        ]
        RBzg = PolynomialRing(RB, ["z", "g"])
        # The following line changes the variables in these minimal polynomials
        # as Z -> z, x -> G, then evaluates at G = QQzg.gens(1) ( = g )
        RBzgG = PolynomialRing(RBzg, "G")
        minpoly_list = [
            RBzgG([c(RBzg.gen(0)) for c in list(h)])(RBzg.gen(1)) for h in minpoly_univ
        ]
        # h(z,g)=0 --> dg/dz = - dhdz/dhdg
        dgdz_list = [
            -h.derivative(RBzg.gen(0)) / h.derivative(RBzg.gen(1)) for h in minpoly_list
        ]

        CCzg = PolynomialRing(self._CC, ["z", "g"])
        CCminpoly_list = [CCzg(h) for h in minpoly_list]

        a0_list = [P(h.leading_coefficient()) for h in minpoly_univ]
        # Note that because the field over which the Riemann surface is defined
        # is embedded into CC, it has characteristic 0, and so we know the
        # irreducible factors are all separable, i.e. the roots have multiplicity
        # one.
        a0_info = [
            (
                self._CC(a0.leading_coefficient()),
                flatten(
                    [
                        self._CCz(F).roots(multiplicities=False) * m
                        for F, m in a0.factor()
                    ]
                ),
            )
            for a0 in a0_list
        ]
        if exact:
            return RBzg, list(zip(integrand_list, dgdz_list, minpoly_list, a0_info))
        else:
            return CCzg, list(zip(integrand_list, dgdz_list, CCminpoly_list, a0_info))

    def rigorous_line_integral(self, upstairs_edge, differentials, bounding_data):
        r"""
        Perform vectorized integration along a straight path.

        Using the error bounds for Gauss-Legendre integration found in [Neu2018]_
        and a method for bounding an algebraic integrand on a circular domains
        using Cauchy's form of the remainder in Taylor approximation coupled to
        Fujiwara's bound on polynomial roots (see Bruin-DisneyHogg-Gao, in
        preparation), this method calculates (semi-)rigorously the integral of a
        list of differentials along an edge of the upstairs graph.

        INPUT:

        - ``upstairs_edge`` -- tuple; either a pair of integer tuples
          corresponding to an edge of the upstairs graph, or a tuple
          ``((z_start, sb), (z_end, ))`` as in the input of
          ``make_zw_interpolator``

        - ``differentials`` -- list of polynomials; a polynomial `g`
          represents the differential `g(z,w)/(df/dw) dz` where `f(z,w)=0` is
          the equation defining the Riemann surface

        - ``bounding_data`` -- tuple containing the data required for bounding
          the integrands. This should be in the form of the output from
          :meth:`_bounding_data`.

        OUTPUT: a complex number, the value of the line integral

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f); S
            Riemann surface defined by polynomial f = -z^4 + w^2 + 1 = 0, with 53 bits of precision

        Since we make use of data from homotopy continuation, we need to compute
        the necessary data::

            sage: _ = S.homology_basis()
            sage: differentials = S.cohomology_basis()
            sage: bounding_data = S._bounding_data(differentials)
            sage: S.rigorous_line_integral([(0,0), (1,0)], differentials, bounding_data)  # abs tol 1e-10
            (1.80277751848459e-16 - 0.352971844594760*I)

        .. NOTE::

            Uses data that ``homology_basis`` initializes, and may give incorrect
            values if :meth:`homology_basis` has not initialized them.

            Note also that the data of the differentials is contained within
            ``bounding_data``. It is, however, still advantageous to have this
            be a separate argument, as it lets the user supply a fast-callable
            version of the differentials, to significantly speed up execution
            of the integrand calls, and not have to re-calculate these
            fast-callables for every run of the function. This is also the benefit
            of representing the  differentials as a polynomial over a known
            common denominator.

        .. TODO::

            Note that bounding_data contains the information of the integrands,
            so one may want to check for consistency between ``bounding_data``
            and ``differentials``. If so one would not want to do so at the
            expense of speed.

            Moreover, the current implementation bounds along a line by
            splitting it up into segments, each of which can be covered entirely
            by a single circle, and then placing inside that the ellipse
            required to bound as per [Neu2018]_. This is reliably more efficient
            than the heuristic method, especially in poorly-conditioned cases
            where discriminant points are close together around the edges, but
            in the case where the branch locus is well separated, it can require
            slightly more nodes than necessary. One may want to include a method
            here to transition in this regime to an algorithm that covers the
            entire line with one ellipse, then bounds along that ellipse with
            multiple circles.
        """
        # Note that this, in its current formalism, makes no check that bounding
        # data at all corresponds to  the differentials given. The onus is then
        # on the design of other functions which use it.

        # CCzg is required to be known as we need to know the ring which the minpolys
        # lie in.
        CCzg, bounding_data_list = bounding_data
        CCz = CCzg.univariate_ring(CCzg.gen(1)).base_ring()

        d_edge = tuple(u[0] for u in upstairs_edge)
        # Using a try-catch here allows us to retain a certain amount of back
        # compatibility for users.
        try:
            initial_continuation = self._L[d_edge]
            upstairs_edge = (
                (self._vertices[d_edge[0]], upstairs_edge[0][1]),
                (self._vertices[d_edge[1]],),
            )
        except KeyError:
            initial_continuation = self.homotopy_continuation(d_edge)

        zwt, z1_minus_z0 = self.make_zw_interpolator(
            upstairs_edge, initial_continuation
        )
        z0 = zwt(0)[0]
        z1 = zwt(1)[0]

        # list of (centre, radius) pairs that still need to be processed
        # None is a sentinel value to indicate that the minimum number of
        # nodes required to integrate on the corresponding segment within
        # the required error tolerance is not yet known.
        ball_stack = [(self._RR(1 / 2), self._RR(1 / 2), None)]
        alpha = self._RR(912 / 1000)
        # alpha set manually for scaling purposes. Basic benchmarking shows
        # that ~0.9 is a sensible value.
        E_global = self._RR(2)**(-self._prec + 3)

        # Output will iteratively store the output of the integral.
        V = VectorSpace(self._CC, len(differentials))
        output = V(0)

        # The purpose of this loop is as follows: We know we will be using
        # Gauss-Legendre quadrature to do the integral, and results from [Neu2018]_
        # tell us an upper bound on the number of nodes required to achieve a
        # given error bound for this quadrature, provided we have a bound for
        # the integrand on a certain ellipse in the complex plane. The method
        # developed by Bruin and Gao that uses Cauchy and Fujiwara can bound an
        # algebraic integrand on a circular region. Hence we need a way to change
        # from bounding with an ellipse to bounding with a circle. The size of
        # these circles will be constrained by the distance to the nearest point
        # where the integrand blows up, i.e. the nearest branchpoint. Basic
        # benchmarking showed that it was in general a faster method to split
        # the original line segment into multiple smaller line segments, and
        # compute the contribution from each of the line segments bounding with
        # a single circle, the benefits mainly coming when the curve is poorly
        # conditioned s.t. the branch points are close together. The following
        # loop does exactly this, repeatedly bisecting a segment if it is not
        # possible to cover it entirely in a ball which encompasses an appropriate
        # ellipse.
        def local_N(ct, rt):
            cz = (1 - ct) * z0 + ct * z1  # This is the central z-value of our ball.
            distances = [(cz - b).abs() for b in self.branch_locus]
            rho_z = min(distances)
            rho_t = rho_z / (z1_minus_z0).abs()
            rho_t = alpha * rho_t + (1 - alpha) * rt  # sqrt(rho_t*rt) could also work
            rho_z = rho_t * (z1_minus_z0).abs()
            delta_z = (alpha * rho_t + (1 - alpha) * rt) * (z1_minus_z0).abs()
            # delta_z and delta_z^2 / (rho_z * (rho_z - delta_z)) are the two
            # prefactors that occur in the computation of the magnitude bound
            # M. delta_z should never be infinite, but the second factor could
            # be if rho_z - delta_z is 0. Mathematically it would never be 0
            # as we ensure rho_t > rt before running local_N, but the
            # floating point operations can ruin this.
            # The second prefactor is actually homogeneous in
            # z1_minus_z0.abs(), so we shall compute this factor without those
            # multiplications as a function of rho_t / rt which should thus be
            # more resistance to floating-point errors.
            pf2 = (alpha + (1 - alpha) * (rt / rho_t))**2 / (
                  (1 - alpha) * (1 - rt / rho_t)
            )
            expr = (
                rho_t / rt + ((rho_t / rt)**2 - 1).sqrt()
            )  # Note this is really exp(arcosh(rho_t/rt))
            Ni = 3
            cw = zwt(ct)[1]
            for g, dgdz, minpoly, (a0lc, a0roots) in bounding_data_list:
                z_1 = a0lc.abs() * prod((cz - r).abs() - rho_z for r in a0roots)
                n = minpoly.degree(CCzg.gen(1))
                ai_new = [
                    CCz(minpoly.coefficient({CCzg.gen(1): i}))(z=cz + self._CCz.gen(0))
                    for i in range(n)
                ]
                ai_pos = [self._RRz([c.abs() for c in h.list()]) for h in ai_new]
                m = [a(rho_z) / z_1 for a in ai_pos]
                l = len(m)
                M_tilde = 2 * max(
                    (m[i].abs())**(1 / self._RR(l - i)) for i in range(l)
                )
                cg = g(cz, cw)
                cdgdz = dgdz(cz, cg)
                M = delta_z * cdgdz.abs() + pf2 * M_tilde
                N_required = (
                    (M * (self._RR.pi() + 64 / (15 * (expr**2 - 1))) / E_global).log()
                    / (2 * expr.log())
                )
                if N_required.is_positive_infinity():
                    return 2**max(60, self._prec)
                Ni = max(Ni, N_required.ceil())
            return Ni

        while ball_stack:
            ct, rt, lN = ball_stack.pop()
            ncts = [ct - rt / 2, ct + rt / 2]
            nrt = rt / 2

            if lN is None:
                cz = (1 - ct) * z0 + ct * z1
                distances = [(cz - b).abs() for b in self.branch_locus]
                rho_z = min(distances)
                rho_t = rho_z / (z1_minus_z0).abs()

                if rho_t <= rt:
                    ball_stack.append((ncts[0], nrt, None))
                    ball_stack.append((ncts[1], nrt, None))
                    continue

                lN = local_N(ct, rt)

            nNs = [local_N(nct, nrt) for nct in ncts]

            if sum(nNs) < lN:
                ball_stack.append((ncts[0], nrt, nNs[0]))
                ball_stack.append((ncts[1], nrt, nNs[1]))
                continue

            if lN % 2 and not lN == 3:
                lN += 1

            ct_minus_rt = ct - rt
            two_rt = 2 * rt

            def integrand(t):
                zt, wt = zwt(ct_minus_rt + t * two_rt)
                dfdwt = self._fastcall_dfdw(zt, wt)
                return V([h(zt, wt) / dfdwt for h in differentials])

            output += two_rt * integrate_vector_N(integrand, self._prec, lN)

        return output * z1_minus_z0

    def matrix_of_integral_values(self, differentials, integration_method='heuristic'):
        r"""
        Compute the path integrals of the given differentials along the homology
        basis.

        The returned answer has a row for each differential. If the Riemann
        surface is given by the equation `f(z,w)=0`, then the differentials are
        encoded by polynomials g, signifying the differential `g(z,w)/(df/dw)
        dz`.

        INPUT:

        - ``differentials`` -- list of polynomials

        - ``integration_method`` -- (default: ``'heuristic'``) string specifying
          the integration method to use. The options are ``'heuristic'`` and
          ``'rigorous'``.

        OUTPUT:

        A matrix, one row per differential, containing the values of the path
        integrals along the homology basis of the Riemann surface.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(x^3 + y^3 + 1)
            sage: B = S.cohomology_basis()
            sage: m = S.matrix_of_integral_values(B)
            sage: parent(m)
            Full MatrixSpace of 1 by 2 dense matrices over Complex Field with 53 bits of precision
            sage: (m[0,0]/m[0,1]).algdep(3).degree() # curve is CM, so the period is quadratic
            2

        .. NOTE::

            If ``differentials is self.cohomology_basis()``, the calculations
            of the integrals along the edges are written to ``self._integral_dict``.
            This is as this data will be required when computing the Abel-Jacobi
            map, and so it is helpful to have is stored rather than recomputing.
        """
        cycles = self.homology_basis()

        def normalize_pairs(L):
            r"""
            Return a list of edges encoded by the path in L.
            The edges are normalized to be in the direction in which
            the homotopy continuation should have been computed along them.
            """
            R = []
            for i in range(len(L) - 1):
                if L[i][0] < L[i + 1][0]:
                    R.append((L[i], L[i + 1]))
                else:
                    R.append((L[i + 1], L[i]))
            return R

        occurring_edges = set()
        occurring_edges.update(*[normalize_pairs(p[1]) for h in cycles for p in h])

        if differentials is self.cohomology_basis():
            fcd = self._fastcall_cohomology_basis
            integral_dict = self._integral_dict
        else:
            fcd = [fast_callable(omega, domain=self._CC) for omega in differentials]
            integral_dict = {}

        if integration_method == "heuristic":
            line_int = lambda edge: self.simple_vector_line_integral(edge, fcd)
        elif integration_method == "rigorous":
            bd = self._bounding_data(differentials)
            line_int = lambda edge: self.rigorous_line_integral(edge, fcd, bd)
        else:
            raise ValueError("invalid integration method")

        integral_dict = {edge: line_int(edge) for edge in occurring_edges}

        rows = []
        for cycle in cycles:
            V = VectorSpace(self._CC, len(differentials)).zero()
            for multiplicity, loop in cycle:
                for i in range(len(loop) - 1):
                    if loop[i][0] < loop[i + 1][0]:
                        direction = 1
                        upstairs_edge = (loop[i], loop[i + 1])
                    else:
                        direction = -1
                        upstairs_edge = (loop[i + 1], loop[i])
                    V += (multiplicity * direction) * integral_dict[upstairs_edge]
            rows.append(V)
        return Matrix(rows).transpose()

    @cached_method
    def period_matrix(self):
        r"""
        Compute the period matrix of the surface.

        OUTPUT: a matrix of complex values

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f, prec=30)
            sage: M = S.period_matrix()

        The results are highly arbitrary, so it is hard to check if the result
        produced is correct. The closely related ``riemann_matrix`` is somewhat
        easier to test.::

            sage: parent(M)
            Full MatrixSpace of 3 by 6 dense matrices
             over Complex Field with 30 bits of precision
            sage: M.rank()
            3

        One can check that the two methods give similar answers::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: f = y^2 - x^3 + 1
            sage: S = RiemannSurface(f, integration_method='rigorous')
            sage: T = RiemannSurface(f, integration_method='heuristic')
            sage: RM_S = S.riemann_matrix()
            sage: RM_T = T.riemann_matrix()
            sage: (RM_S-RM_T).norm() < 1e-10
            True
        """
        differentials = self.cohomology_basis()
        return self.matrix_of_integral_values(differentials, self._integration_method)

    def riemann_matrix(self):
        r"""
        Compute the Riemann matrix.

        OUTPUT: a matrix of complex values

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f, prec=60)
            sage: M = S.riemann_matrix()

        The Klein quartic has a Riemann matrix with values in a quadratic
        field::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2 - x + 2)
            sage: all(len(m.algdep(6).roots(K)) > 0 for m in M.list())
            True
        """
        PeriodMatrix = self.period_matrix()
        Am = PeriodMatrix[0 : self.genus, 0 : self.genus]
        RM = (
            numerical_inverse(Am)
            * PeriodMatrix[0 : self.genus, self.genus : 2 * self.genus]
        )
        return RM

    def plot_paths(self):
        r"""
        Make a graphical representation of the integration paths.

        This returns a two dimensional plot containing the branch points (in red) and
        the integration paths (obtained from the Voronoi cells of the branch
        points). The integration paths are plotted by plotting the points that
        have been computed for homotopy continuation, so the density gives an
        indication of where numerically sensitive features occur.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2 - x^3 - x)
            sage: S.plot_paths()                                                        # needs sage.plot
            Graphics object consisting of 2 graphics primitives
        """
        from sage.plot.point import point2d

        P = []

        # trigger the computation of the homology basis, so that self._L is present
        self.homology_basis()

        for e in self._L.keys():
            z0 = self._vertices[e[0]]
            z1 = self._vertices[e[1]]

            def path(t):
                return (1 - t) * z0 + t * z1

            T = self._L[e]
            P += [path(t[0]) for t in T]
        return point2d(P, size=1) + point2d(self.branch_locus, color='red')

    def plot_paths3d(self, thickness=0.01):
        r"""
        Return the homology basis as a graph in 3-space.

        The homology basis of the surface is constructed by taking the Voronoi
        cells around the branch points and taking the inverse image of the edges
        on the Riemann surface. If the surface is given by the equation
        `f(z,w)`, the returned object gives the image of this graph in 3-space
        with coordinates `\left(\operatorname{Re}(z), \operatorname{Im}(z),
        \operatorname{Im}(w)\right)`.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2 - x^3 - x)
            sage: S.plot_paths3d()                                                      # needs sage.plot
            Graphics3d Object
        """
        from sage.plot.graphics import Graphics
        from sage.plot.plot3d.shapes2 import point3d, line3d

        P = Graphics()

        # trigger the computation of the homology basis, so that
        # self._L is present
        self.homology_basis()

        for e in self._L.keys():
            z0 = self._vertices[e[0]]
            z1 = self._vertices[e[1]]

            def path(t):
                z = (1 - t) * z0 + t * z1
                return (z.real_part(), z.imag_part())

            T = self._L[e]
            color = "blue"
            for i in range(self.degree):
                P += line3d(
                    [path(t[0]) + (t[1][i].imag_part(),) for t in T],
                    color=color,
                    thickness=thickness,
                )
            for z, ws in zip(self._vertices, self._wvalues):
                for w in ws:
                    P += point3d(
                        [z.real_part(), z.imag_part(), w.imag_part()],
                        color='purple',
                        size=20,
                    )
        return P

    def endomorphism_basis(self, b=None, r=None):
        r"""
        Numerically compute a `\ZZ`-basis for the endomorphism ring.

        Let `\left(I | M \right)` be the normalized period matrix (`M` is the
        `g\times g` :meth:`riemann_matrix`). We consider the system of matrix
        equations `MA + C = (MB + D)M` where `A, B, C, D` are `g\times g`
        integer matrices.  We determine small integer (near) solutions using LLL
        reductions.  These solutions are returned as `2g \times 2g` integer
        matrices obtained by stacking `\left(D | B\right)` on top of `\left(C |
        A\right)`.

        INPUT:

        - ``b`` -- integer (default provided); the equation coefficients are
          scaled by `2^b` before rounding to integers

        - ``r`` -- integer (default: ``b/4``); solutions that have all
          coefficients smaller than `2^r` in absolute value are reported as
          actual solutions

        OUTPUT:

        A list of `2g \times 2g` integer matrices that, for large enough ``r``
        and ``b-r``, generate the endomorphism ring.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(x^3 + y^3 + 1)
            sage: B = S.endomorphism_basis(); B #random
            [
            [1 0]  [ 0 -1]
            [0 1], [ 1  1]
            ]
            sage: sorted([b.minpoly().disc() for b in B])
            [-3, 1]
        """
        M = self.riemann_matrix()
        return integer_matrix_relations(M, M, b, r)

    def homomorphism_basis(self, other, b=None, r=None):
        r"""
        Numerically compute a `\ZZ`-basis for module of homomorphisms to a given
        complex torus.

        Given another complex torus (given as the analytic Jacobian of a Riemann
        surface), numerically compute a basis for the homomorphism module. The
        answer is returned as a list of `2g \times 2g` integer matrices `T=(D, B; C, A)`
        such that if the columns of `(I|M_1)` generate the lattice defining the
        Jacobian of the Riemann surface and the columns of `(I|M_2)` do this for
        the codomain, then approximately we have `(I|M_2)T=(D+M_2C)(I|M_1)`, i.e., up
        to a choice of basis for `\CC^g` as a complex vector space, we we
        realize `(I|M_1)` as a sublattice of `(I|M_2)`.

        INPUT:

        - ``b`` -- integer (default provided); the equation coefficients are
          scaled by `2^b` before rounding to integers

        - ``r`` -- integer (default: ``b/4``); solutions that have all
          coefficients smaller than `2^r` in absolute value are reported as
          actual solutions

        OUTPUT:

        A list of `2g \times 2g` integer matrices that, for large enough ``r``
        and ``b-r``, generate the homomorphism module.

        EXAMPLES::

            sage: S1 = EllipticCurve("11a1").riemann_surface()
            sage: S2 = EllipticCurve("11a3").riemann_surface()
            sage: [m.det() for m in S1.homomorphism_basis(S2)]
            [5]
        """
        M1 = self.riemann_matrix()
        M2 = other.riemann_matrix()
        return integer_matrix_relations(M2, M1, b, r)

    def tangent_representation_numerical(self, Rs, other=None):
        r"""
        Compute the numerical tangent representations corresponding to the
        homology representations in ``Rs``.

        The representations on homology ``Rs`` have to be given with respect to
        the symplectic homology basis of the Jacobian of ``self`` and ``other``.
        Such matrices can for example be obtained via
        :meth:`endomorphism_basis`.

        Let `P` and `Q` be the period matrices of ``self`` and ``other``. Then
        for a homology representation `R`, the corresponding tangential
        representation `T` satisfies `T P = Q R`.

        INPUT:

        - ``Rs`` -- set of matrices on homology to be converted to their
          tangent representations

        - ``other`` -- (default: ``self``) the codomain; another Riemann
          surface

        OUTPUT: the numerical tangent representations of the matrices in ``Rs``

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: A.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2 - (x^6 + 2*x^4 + 4*x^2 + 8), prec = 100)
            sage: P = S.period_matrix()
            sage: Rs = S.endomorphism_basis()
            sage: Ts = S.tangent_representation_numerical(Rs)
            sage: all(((T*P - P*R).norm() < 2^(-80)) for [T, R] in zip(Ts, Rs))
            True
        """
        if not other:
            other = self
        P = self.period_matrix()
        CCP = P.base_ring()
        g = self.genus
        Q = other.period_matrix()
        Ptsubinv = numerical_inverse((P.transpose())[list(range(g))])
        Ts = []
        for R in Rs:
            QRtsub = ((Q * R).transpose())[list(range(g))]
            Tt = Ptsubinv * QRtsub
            T = Tt.transpose().change_ring(CCP)
            Ts.append(T)
        return Ts

    def tangent_representation_algebraic(self, Rs, other=None, epscomp=None):
        r"""
        Compute the algebraic tangent representations corresponding to the
        homology representations in ``Rs``.

        The representations on homology ``Rs`` have to be given with respect to
        the symplectic homology basis of the Jacobian of ``self`` and ``other``.
        Such matrices can for example be obtained via
        :meth:`endomorphism_basis`.

        Let `P` and `Q` be the period matrices of ``self`` and ``other``. Then
        for a homology representation `R`, the corresponding tangential
        representation `T` satisfies `T P = Q R`.

        INPUT:

        - ``Rs`` -- set of matrices on homology to be converted to their
          tangent representations

        - ``other`` -- (default: ``self``) the codomain; another Riemann
          surface

        - ``epscomp`` -- real number (default: ``2^(-prec + 30)``). Used to
          determine whether a complex number is close enough to a root of a
          polynomial.

        OUTPUT: the algebraic tangent representations of the matrices in ``Rs``

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: A.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2 - (x^6 + 2*x^4 + 4*x^2 + 8), prec = 100)
            sage: Rs = S.endomorphism_basis()
            sage: Ts = S.tangent_representation_algebraic(Rs)
            sage: Ts[0].base_ring().maximal_order().discriminant() == 8
            True
        """
        if not epscomp:
            epscomp = 2**(-self._prec + 30)
        QQalg = QQ.algebraic_closure()

        def polynomialize_element(alpha):
            d = 1
            while True:
                d += 1
                dep = algdep(alpha, d, height_bound=10**d)
                if dep and dep(alpha) < epscomp:
                    return dep

        def algebraize_element(alpha):
            alphaPol = polynomialize_element(alpha)
            CC = alpha.parent()
            for tup in alphaPol.roots(QQalg):
                rt = tup[0]
                if (alpha - CC(rt)).abs() < epscomp:
                    return rt
            raise AssertionError("No close root found while algebraizing")

        def algebraize_matrices(Ts):
            nr = Ts[0].nrows()
            nc = Ts[0].ncols()
            TsAlg = [T.apply_map(algebraize_element) for T in Ts]
            elts = [x for TAl in TsAlg for x in TAl.list()]
            eltsAlg = number_field_elements_from_algebraics(elts)[1]
            L = eltsAlg[0].parent()
            TsAlgL = []
            for i in range(len(Ts)):
                TAlgL = [eltsAlg[j] for j in range(i * nr * nc, (i + 1) * nr * nc)]
                TsAlgL.append(Matrix(L, nr, nc, TAlgL))
            return TsAlgL

        Ts = self.tangent_representation_numerical(Rs, other=other)
        return algebraize_matrices(Ts)

    def rosati_involution(self, R):
        r"""
        Compute the Rosati involution of an endomorphism.

        The endomorphism in question should be given by its homology
        representation with respect to the symplectic basis of the Jacobian.

        INPUT:

        - ``R`` -- integral matrix

        OUTPUT: the result of applying the Rosati involution to ``R``

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: A.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2 - (x^6 + 2*x^4 + 4*x^2 + 8), prec = 100)
            sage: Rs = S.endomorphism_basis()
            sage: S.rosati_involution(S.rosati_involution(Rs[1])) == Rs[1]
            True
        """

        def standard_symplectic_matrix(n):
            one = Matrix.identity(n)
            zero = Matrix.zero(n)
            return Matrix.block([[zero, -one], [one, zero]])

        g = self.genus
        if not (R.nrows() == 2 * g == R.ncols()):
            raise AssertionError(
                "Matrix is not the homology representation of an endomorphism"
            )
        J = standard_symplectic_matrix(g)
        return -J * R.transpose() * J

    def symplectic_isomorphisms(self, other=None, hom_basis=None, b=None, r=None):
        r"""
        Numerically compute symplectic isomorphisms.

        INPUT:

        - ``other`` -- (default: ``self``) the codomain; another Riemann
          surface

        - ``hom_basis`` -- (default: ``None``) a `\ZZ`-basis of the
          homomorphisms from ``self`` to ``other``, as obtained from
          :meth:`homomorphism_basis`. If you have already calculated this
          basis, it saves time to pass it via this keyword argument. Otherwise
          the method will calculate it.

        - ``b`` -- integer (default provided); as for
          :meth:`homomorphism_basis`, and used in its invocation if
          (re)calculating said basis

        - ``r`` -- integer (default: ``b/4``);  as for
          :meth:`homomorphism_basis`, and used in its invocation if
          (re)calculating said basis

        OUTPUT:

        This returns the combinations of the elements of
        :meth:`homomorphism_basis` that correspond to symplectic
        isomorphisms between the Jacobians of ``self`` and ``other``.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: f = y^2 - (x^6 + 2*x^4 + 4*x^2 + 8)
            sage: X = RiemannSurface(f, prec=100)
            sage: P = X.period_matrix()
            sage: g = y^2 - (x^6 + x^4 + x^2 + 1)
            sage: Y = RiemannSurface(g, prec=100)
            sage: Q = Y.period_matrix()
            sage: Rs = X.symplectic_isomorphisms(Y)
            sage: Ts = X.tangent_representation_numerical(Rs, other = Y)
            sage: test1 = all(((T*P - Q*R).norm() < 2^(-80)) for [T, R] in zip(Ts, Rs))
            sage: test2 = all(det(R) == 1 for R in Rs)
            sage: test1 and test2
            True
        """
        if not other:
            other = self
        if hom_basis:
            Rs = hom_basis
        else:
            Rs = self.homomorphism_basis(other=other, b=b, r=r)
        r = len(Rs)
        g = self.genus
        A = PolynomialRing(QQ, r, "x")
        gensA = A.gens()
        # Use that the trace is positive definite; we could also put this as an
        # extra condition when determining the endomorphism basis to speed up
        # that calculation slightly
        R = sum(gensA[i] * Rs[i].change_ring(A) for i in range(r))
        tr = (R * self.rosati_involution(R)).trace()
        # Condition tr = 2 g creates ellipsoid
        M = Matrix(
            ZZ,
            r,
            r,
            [tr.derivative(gen1).derivative(gen2) for gen1 in gensA for gen2 in gensA],
        )
        vs = M.__pari__().qfminim(4 * g)[2].sage().transpose()
        vs = [v for v in vs if v * M * v == 4 * g]
        vs += [-v for v in vs]
        RsIso = []
        for v in vs:
            R = sum(v[i] * Rs[i] for i in range(r))
            if R * self.rosati_involution(R) == 1:
                RsIso.append(R)
        return RsIso

    def symplectic_automorphism_group(self, endo_basis=None, b=None, r=None):
        r"""
        Numerically compute the symplectic automorphism group as a permutation
        group.

        INPUT:

        - ``endo_basis`` -- (default: ``None``) a `\ZZ`-basis of the
          endomorphisms of ``self``, as obtained from
          :meth:`endomorphism_basis`. If you have already calculated this
          basis, it saves time to pass it via this keyword argument. Otherwise
          the method will calculate it.

        - ``b`` -- integer (default provided); as for
          :meth:`homomorphism_basis`, and used in its invocation if
          (re)calculating said basis

        - ``r`` -- integer (default: ``b/4``);  as for
          :meth:`homomorphism_basis`, and used in its invocation if
          (re)calculating said basis

        OUTPUT:

        The symplectic automorphism group of the Jacobian of the Riemann
        surface. The automorphism group of the Riemann surface itself can be
        recovered from this; if the curve is hyperelliptic, then it is
        identical, and if not, then one divides out by the central element
        corresponding to multiplication by -1.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: A.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2 - (x^6 + 2*x^4 + 4*x^2 + 8), prec = 100)
            sage: G = S.symplectic_automorphism_group()
            sage: G.as_permutation_group().is_isomorphic(DihedralGroup(4))
            True
        """
        RsAut = self.symplectic_isomorphisms(hom_basis=endo_basis, b=b, r=r)
        return MatrixGroup(RsAut)

    def __add__(self, other):
        r"""
        Return the disjoint union of the Riemann surface and the other argument.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y> = QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S1+S1
            Riemann surface sum with period lattice of rank 4
        """
        return RiemannSurfaceSum([self, other])

    def _integrate_differentials_iteratively(
        self, upstairs_edge, cutoff_individually=False, raise_errors=True, prec=None
    ):
        r"""
        Integrate the cohomology basis along a straight line edge.

        The cohomology basis is integrated along a straight line using a version
        of the double exponential quadrature. This method of integrating the
        cohomology basis is especially useful when integrating out to infinity,
        or near roots of ``self._dfdw``. In order to aid with convergence of the
        method, two main modification to a standard integrator are made, most
        importantly of which is the truncation of the integral near branch points,
        where the first term in the Puiseux series of the integrands are used to
        approximately bound the integral. The ``cutoff_individually`` parameter
        allows the user to set whether that truncation is uniform over all the
        integrands, which improves the complexity of the algorithm, but loses
        the ability to gain benefits where integrands vanish to a high order at
        the branchpoint.

        INPUT:

        - ``upstairs_edge`` -- tuple. A tuple of complex numbers of the form
          ``((z_start, w_start), z_end)`` specifying the path to integrate
          along, where ``z_start`` may be infinite, in which case ``w_start``
          must be an integer specifying the branch.

        - ``cutoff_individually`` -- boolean (default: ``False``); whether to
          truncate the integrand uniformly or not. If ``None``, then no
          truncation is applied

        - ``raise_errors`` -- boolean (default: ``True``); by default the code
          uses convergence errors to ensure any answers returned are accurate.
          This can be turned off to return answers faster that are not
          necessarily correct.

        - ``prec`` -- integer (default: ``self._prec``); the precision to try
          and achieve, defined as `2^{-\text{prec}+3}`

        OUTPUT:

        Tuple ``(I, gs)`` where ``I`` is the vector of integrals, and ``gs`` are
        the values of the differentials at ``z_end``.

        EXAMPLES:

        We know that for the surface given by `w^2-z^4-1` a cohomology basis is
        given by `\frac{dz}{2w}`. One can verify analytically that
        `\int_0^1 frac{dt}{\sqrt{1-t^4}}=\frac{\sqrt{\pi}\Gamma(5/4)}{\Gamma(3/4)}`,
        and we check this with the integrator, being careful with signs::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: S = RiemannSurface(w^2 + z^4 - 1, prec=100)
            sage: branch = 0
            sage: eps = S._RR(2)**(-S._prec)
            sage: z_start = 1 - eps
            sage: z_end = 0
            sage: w_start = S.w_values(z_start)[0]
            sage: s = sign(w_start)
            sage: u_edge = ((z_start, w_start), z_end)
            sage: J, _ = S._integrate_differentials_iteratively(u_edge)
            sage: bool(J[0] + s*S._RR(sqrt(pi)*gamma(5/4)/gamma(3/4)/2) < 1e-10)        # needs sage.symbolic
            True

        .. NOTE::

            The cutoff methodology is calculating the first term in the Puiseux
            series of the differentials about z_start. In future it may be
            desirable to extend this further and use the truncated Puiseux series
            entirely to integrate the differentials.
        """
        (z_start, w_start), z_end = upstairs_edge
        z_start = self._CC(z_start)
        z_end = self._CC(z_end)

        if z_end == self._CC(Infinity):
            raise NotImplementedError

        _, bounding_data_list = self._cohomology_basis_bounding_data
        mp_list = [bd[2] for bd in bounding_data_list]

        # Parameterise so zbar=0 corresponds to z=z_start
        mp_list = [reparameterize_differential_minpoly(mp, z_start) for mp in mp_list]

        # Depending on whether we have reparameterized about infinity or not,
        # we initialise some values we will need in the calculation, inclduing
        # the function `initalize', which at a given value of zbar, calculates
        # the starting value for the i-th differential so it can be iterated
        # from via homotopy continuation.
        if z_start == self._CC(Infinity):
            CCzg = PolynomialRing(self._CC, ["zbar", "gbar"])
            mp_list = [CCzg(mp) for mp in mp_list]
            J = 1 / z_end
            endscale = -(z_end**(-2))

            def initialise(z, i):
                DF = ComplexField(2 * self._prec)
                DFw = PolynomialRing(DF, "wbar")
                z = DF(z)
                R = DF(z**(-1))
                wR = DFw(self.f(R, DFw.gen(0))).roots(multiplicities=False)[w_start]
                newg = -(R**2) * self.cohomology_basis()[i](R, wR) / self._dfdw(R, wR)
                err = mp_list[i](z, newg).abs()
                if err > tau:
                    rs = mp_list[i](z, DFw.gen(0)).roots(multiplicities=False)
                    sb = find_closest_element(newg, rs)
                    newg = rs[sb]
                return newg

        else:
            CCzg = mp_list[0].parent()
            J = z_end - z_start
            endscale = 1

            def initialise(z, i):
                newg = self.cohomology_basis()[i](z_start, w_start) / self._dfdw(
                    z_start, w_start
                )
                err = mp_list[i](z, newg).abs()
                if err > tau:
                    rs = mp_list[i](z, self._CCw.gen(0)).roots(multiplicities=False)
                    sb = find_closest_element(newg, rs)
                    newg = rs[sb]
                return newg

        # As multiple calls of the minimal polynomial and it's derivative will
        # be required for the homotopy continuaiton, we create fast-callable
        # versions of these.
        fc_mp_list = [fast_callable(mp, domain=self._CC) for mp in mp_list]
        fc_dmp_list = [
            fast_callable(mp.derivative(CCzg.gen(1)), domain=self._CC) for mp in mp_list
        ]

        if prec is None:
            prec = self._prec
        # tau here is playing the role of the desired error.
        tau = self._RR(2)**(-prec + 3)
        one = self._RR(1)
        la = self._RR.pi() / 2

        # Cutoffs are used to allow us to not have to integrate as close into
        # a singularity as we might otherwise have to, by knowing that we can
        # truncate the integration interval and only introduce a finite error
        # that can be bounded by knowledge of the asymptotics of the integrands,
        # which we have from their minimal polynomials. This is really a
        # precursor to what would be ideal to implement eventually, namely
        # a method that uses Puiseux series to integrate into singularities.
        # We allow for cutoffs to be tailored to each integrand, or we take a
        # uniform value.
        if cutoff_individually is None:
            cutoffs = [0]
            cutoff_individually = False
        else:
            cutoffs = []
            A = PolynomialRing(self._CC, "xyz")
            aes = []
            for mp in mp_list:
                d = mp.monomial_coefficients()
                mp = sum(
                    [
                        d[k] * CCzg.gen(0)**k[0] * CCzg.gen(1)**k[1]
                        for k in d.keys()
                        if d[k].abs() > tau
                    ]
                )
                cst = min([iz for (iz, ig) in d.keys() if ig == 0])
                a = QQ(max([(cst - iz) / ig for (iz, ig) in d.keys() if ig > 0]))
                sum_coeffs = sum(
                    [
                        d[k] * A.gen(0)**k[1]
                        for k in d.keys()
                        if ((k[1] == 0 and k[0] == cst) or k[1] * a + k[0] - cst == 0)
                    ]
                )
                G = max([r.abs() for r in sum_coeffs.roots(multiplicities=False)])
                cutoffs.append(((a + 1) * tau / G)**(1 / self._CC(a + 1)) / J.abs())
                aes.append(a)
            cutoff_individually = bool(
                not all(ai <= 0 for ai in aes) and cutoff_individually
            )

        # The `raise_errors' variable toggles what we do in the event that
        # newton iteration hasn't converged to the desired precision in a
        # fixed number of steps, here set to 100.
        # If the default value of True is taken, then the failure to converge
        # raises an error. If the value of False is taken, this failure to
        # converge happens silently, thus allowing the user to get *an*
        # answer out of the integration, but numerical imprecision is to be
        # expected. As such, we set the maximum number of steps in the sequence
        # of DE integrations to be lower in the latter case.
        if raise_errors:
            n_steps = self._prec - 1
        else:
            n_steps = 15

        from sage.functions.log import lambert_w

        V = VectorSpace(self._CC, self.genus)
        h = one
        Nh = (-lambert_w(-1, -tau / 2) / la).log().ceil()
        h0 = Nh * h

        # Depending on how the cutoffs were defined, we now create the function
        # which calculates the integrand we want to integrate via double-
        # exponential methods. This will get the value at the next node by
        # homotopy-continuing from the last node value. There is also a slight
        # technical condition which implements the cutoffs.
        if cutoff_individually:
            z_fc_list = list(zip(fc_mp_list, fc_dmp_list))

            def fv(hj, previous_estimate_and_validity):
                u2 = la * hj.sinh()
                t = 1 / (2 * u2.exp() * u2.cosh())
                z0 = J * t
                outg = []
                valid = self.genus * [True]
                previous_estimate, validity = previous_estimate_and_validity
                for i in range(self.genus):
                    co = cutoffs[i]
                    pv = validity[i]
                    if t < co:
                        outg.append(0)
                        valid[i] = False
                    elif not pv:
                        outg.append(initialise(z0, i))
                    else:
                        F, dF = z_fc_list[i]
                        oldg = previous_estimate[i]
                        delta = F(z0, oldg) / dF(z0, oldg)
                        Ndelta = delta.norm()
                        newg = oldg - delta
                        for j in range(100):
                            new_delta = F(z0, newg) / dF(z0, newg)
                            Nnew_delta = new_delta.norm()
                            if (new_delta == 0) or (
                                Nnew_delta >= Ndelta
                                and (Ndelta.sign_mantissa_exponent()[2] + self._prec)
                                < newg.norm().sign_mantissa_exponent()[2]
                            ):
                                outg.append(newg)
                                break
                            delta = new_delta
                            Ndelta = Nnew_delta
                            newg -= delta
                        else:
                            if raise_errors:
                                raise ConvergenceError("Newton iteration fails to converge")
                            else:
                                outg.append(newg)
                fj = V(outg)
                u1 = la * hj.cosh()
                w = u1 / (2 * u2.cosh()**2)
                return (fj, valid), w * fj

            f0, v0 = fv(h0, (self.genus * [0], self.genus * [False]))
        else:
            cutoffs.append(1)
            cutoff = min(cutoffs)
            cutoff_z = J * cutoff
            J -= cutoff_z

            def fv(hj, previous_estimate):
                u2 = la * hj.sinh()
                t = 1 / (2 * u2.exp() * u2.cosh())
                z0 = cutoff_z + J * t
                outg = []
                for F, dF, oldg in zip(fc_mp_list, fc_dmp_list, previous_estimate):
                    delta = F(z0, oldg) / dF(z0, oldg)
                    Ndelta = delta.norm()
                    newg = oldg - delta
                    for j in range(100):
                        new_delta = F(z0, newg) / dF(z0, newg)
                        Nnew_delta = new_delta.norm()
                        if (new_delta == 0) or (
                            Nnew_delta >= Ndelta
                            and (Ndelta.sign_mantissa_exponent()[2] + self._prec)
                            < newg.norm().sign_mantissa_exponent()[2]
                        ):
                            outg.append(newg)
                            break
                        delta = new_delta
                        Ndelta = Nnew_delta
                        newg -= delta
                    else:
                        if raise_errors:
                            raise ConvergenceError("Newton iteration fails to converge")
                        else:
                            outg.append(newg)
                fj = V(outg)
                u1 = la * hj.cosh()
                w = u1 / (2 * u2.cosh()**2)
                return fj, w * fj

            u1, u2 = (la * h0.cosh(), la * h0.sinh())
            y, w = (1 / (2 * u2.exp() * u2.cosh()), u1 / (2 * u2.cosh() ** 2))
            z0 = cutoff_z + J * y
            f0 = [initialise(z0, i) for i in range(self.genus)]
            f0 = V(f0)
            v0 = w * f0

        D3_over_tau = v0.norm(Infinity)
        D4 = D3_over_tau
        results = []

        # we now calculate the integral via double-exponential methods
        # repeatedly halving the step size and then using a heuristic
        # convergence check. The maximum number of steps allowed is
        # currently set to make sure the step size does not fall below the
        # resolution set by the binary precision used.
        for k in range(n_steps):
            hj = h0
            val = v0
            fj = f0
            for j in range(2 * Nh):
                hj -= h
                try:
                    fj, v = fv(hj, fj)
                except ConvergenceError:
                    break
                D3_over_tau = max(v.norm(Infinity), D3_over_tau)
                val += v
            if j == 2 * Nh - 1:
                results.append(h * val)
                D4 = max(D4, v.norm(Infinity))
            if len(results) > 2:
                if results[-1] == results[-2] or results[2] == results[-3]:
                    D = tau
                else:
                    D1 = (results[-1] - results[-2]).norm(Infinity)
                    D2 = (results[-1] - results[-3]).norm(Infinity)
                    D = min(
                        one,
                        max(
                            D1**(D1.log() / D2.log()),
                            D2**2,
                            tau * D3_over_tau,
                            D4,
                            tau,
                        ),
                    )

                if D <= tau:
                    if cutoff_individually:
                        fj = fj[0]
                    return J * results[-1], endscale * fj
            h /= 2
            Nh *= 2
        # Note that throughout this loop there is a return statement, intended
        # to be activated when the sequence of integral approximations is
        # deemed to have converged by the heuristic error. If this has no
        # happened by the time we have gone through the process n_steps times,
        # we have one final error handle. Again, this will throw an error if
        # the raise_errors flag is true, but will just return the answer otherwise.
        if raise_errors:
            raise ConvergenceError("Newton iteration fails to converge")

        return (J * results[-1], endscale * fj)

    def _aj_based(self, P):
        r"""
        Return the Abel-Jacobi map to ``P`` from ``self._basepoint``.

        Computes a representative of the Abel-Jacobi map from ``self._basepoint``
        to ``P`` via a well-chosen vertex ``V``. The representative given will be
        dependent on the path chosen.

        INPUT:

        - ``P`` -- tuple. A pair giving the endpoint of the integral, either in
          the form ``(z, w)`` or ``(Infinity, branch)``, where in the latter case
          we are using the convention that the `w` value over `\infty` is given by
          the limit as ``z`` tends to `\infty` of ``self.w_values(z)[branch]``.

        OUTPUT: a vector of length ``self.genus``

        EXAMPLES:

        As the output of ``_aj_based`` is difficult to interpret due to its path
        dependency, we look at the output of :meth:`abel_jacobi`. We check for
        two hyperelliptic curves that the Abel-Jacobi map between two branch
        points is a 2-torsion point over the lattice. Note we must remember to
        reduce over the period lattice, as results are path dependent::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: p = 100
            sage: S = RiemannSurface(y^2-x^3+1, prec=p)
            sage: divisor = [(-1, (Infinity, 0)), (1, (1, 0))]
            sage: AJ = S.abel_jacobi(divisor)
            sage: AJx2 = [2*z for z in AJ]
            sage: vector(AJx2).norm()  # abs tol 1e-10
            2.4286506478875809114000865640
            sage: bool(S.reduce_over_period_lattice(AJx2).norm() < 1e-10)
            True
            sage: S = RiemannSurface(y^2-x^4+1, prec=p)
            sage: divisor = [(-1, (-1, 0)), (1, (1, 0))]
            sage: AJ = S.abel_jacobi(divisor)
            sage: AJx2 = [2*z for z in AJ]
            sage: bool(S.reduce_over_period_lattice(AJx2).norm() < 1e-10)
            True
        """
        #####
        fcd = self._fastcall_cohomology_basis

        if self._integration_method == "heuristic":
            line_int = lambda edge: self.simple_vector_line_integral(edge, fcd)
        else:
            bd = self._cohomology_basis_bounding_data
            line_int = lambda edge: self.rigorous_line_integral(edge, fcd, bd)
        #####
        B = self._basepoint
        zP, wP = P

        try:
            Inf = bool(zP == zP.parent()(Infinity))
        except TypeError:
            Inf = False

        if Inf:
            zV = self._vertices[B[0]]
            if zV == 0:
                zV += 1
            upstairs_edge = (P, zV)
            ci = bool(self._CC(Infinity) in self._differentials_branch_locus)
            AJ, endgs = self._integrate_differentials_iteratively(
                upstairs_edge, cutoff_individually=ci
            )
            AJ = -AJ
            g0e = endgs[0]
            ws = self.w_values(zV)
            g0s = [self.cohomology_basis()[0](zV, wi) / self._dfdw(zV, wi) for wi in ws]
            W_index = find_closest_element(g0e, g0s)
            if (
                g0e
                - self.cohomology_basis()[0](zV, ws[W_index])
                / self._dfdw(zV, ws[W_index])
            ).abs() > 1e-10:
                raise ConvergenceError(
                    "Integrand continuation failed to get representative values, higher precision required."
                )
            V_index = B[0]
        else:
            zP = self._CC(zP)
            wP = self._CC(wP)
            V_index = find_closest_element(zP, self._vertices)

            if zP == self._vertices[V_index]:
                W_index = find_closest_element(wP, self._wvalues[V_index])
                AJ = 0
            else:
                b_index = find_closest_element(zP, self.branch_locus)
                b = self.branch_locus[b_index]
                # bl = self.branch_locus+self._differentials_branch_locus
                # b_index = find_closest_element(zP, bl)
                # b = bl[b_index]

                scale = max(b.abs() for b in self.branch_locus)
                d1 = self._CC(1e-2) * scale

                # We choose the first vertex we want to go to.
                # If the closest vertex is closer than the nearest branch point, just take that vertex
                # otherwise we need something smarter.
                delta = self._RR(2)**(-self._prec + 1)
                if not (
                    (zP - self._vertices[V_index]).abs() < (zP - b).abs()
                    or (zP - b).abs() <= delta
                ):
                    region = self.voronoi_diagram.regions[
                        self.voronoi_diagram.point_region[b_index]
                    ]
                    args = [
                        (self._vertices[i] - zP).argument() - (b - zP).argument()
                        for i in region
                    ]
                    suitable_vertex_indices = [
                        region[i]
                        for i in range(len(region))
                        if args[i].abs() - self._RR.pi() / 2 >= -self._RR(1e-15)
                    ]
                    suitable_vertices = [
                        self._vertices[i] for i in suitable_vertex_indices
                    ]
                    if suitable_vertices == []:
                        raise ValueError(
                            "There is no satisfactory choice of V for zP={}".format(zP)
                        )
                    V_index = suitable_vertex_indices[
                        find_closest_element(zP, suitable_vertices)
                    ]
                #####
                zV = self._vertices[V_index]

                if (zP - b).abs() >= d1 or b in self._differentials_branch_locus:
                    wP_index = find_closest_element(wP, self.w_values(zP))
                    d_edge = (zP, zV)
                    u_edge = ((zP, wP_index), (zV,))
                    initial_continuation = self.homotopy_continuation(d_edge)
                    AJ = -line_int(u_edge)

                    w_end = initial_continuation[-1][1][wP_index]
                    W_index = find_closest_element(w_end, self._wvalues[V_index])
                else:
                    zs = zP
                    ws = wP

                    #####
                    # Here we need a block of code to change the vertex if the path
                    # from zP to zV would go through a ramification point of the integrands
                    fl = [
                        c
                        for c in self._differentials_branch_locus
                        if not c == self._CC(Infinity)
                    ]
                    ts = [
                        ((c - zP) * (zV - zP).conjugate()).real()
                        / (zP - zV).norm()**2
                        for c in fl
                    ]
                    ds = [
                        (fl[i] - zP - ts[i] * (zV - zP)).abs()
                        for i in range(len(ts))
                        if (ts[i] >= 0 and ts[i] <= 1)
                    ]
                    while len(ds) >= 1 and min(ds) < delta:
                        V_index = suitable_vertex_indices.pop()
                        zV = self._vertices[V_index]
                        ts = [
                            ((c - zP) * (zV - zP).conjugate()).real()
                            / (zP - zV).norm()**2
                            for c in fl
                        ]
                        ds = [
                            (fl[i] - zP - ts[i] * (zV - zP)).abs()
                            for i in range(len(ts))
                            if (ts[i] >= 0 and ts[i] <= 1)
                        ]
                    #####

                    while self._dfdw(zs, ws).abs() == 0:
                        zs = zs + delta * (zV - zs) / (zV - zs).abs() / 2
                        ws_list = self.w_values(zs)
                        wP_index = find_closest_element(ws, ws_list)
                        ws = ws_list[wP_index]
                    upstairs_edge = ((zs, ws), zV)
                    AJ, endgs = self._integrate_differentials_iteratively(
                        upstairs_edge, cutoff_individually=False
                    )
                    AJ = -AJ
                    g0e = endgs[0]

                    ws = self.w_values(zV)
                    g0s = [
                        self.cohomology_basis()[0](zV, wi) / self._dfdw(zV, wi)
                        for wi in ws
                    ]
                    W_index = find_closest_element(g0e, g0s)
                    if (
                        g0e
                        - self.cohomology_basis()[0](zV, ws[W_index])
                        / self._dfdw(zV, ws[W_index])
                    ).abs() > 1e-10:
                        raise ConvergenceError(
                            "Integrand continuation failed to get representative values, higher precision required."
                        )

        uV_index = (V_index, W_index)
        #####
        G = self.upstairs_graph()
        path = G.shortest_path(B, uV_index)
        edges = [(path[i], path[i + 1]) for i in range(len(path) - 1)]
        #####
        for e in edges:
            if e[1][0] > e[0][0]:
                s = 1
            else:
                s = -1
                e = tuple(reversed(e))
            try:
                AJ += s * self._integral_dict[e]
            except KeyError:
                Ie = line_int(e)
                self._integral_dict[e] = Ie
                AJ += s * Ie
        return AJ

    def abel_jacobi(self, divisor, verbose=False):
        r"""
        Return the Abel-Jacobi map of ``divisor``.

        Return a representative of the Abel-Jacobi map of a divisor with basepoint
        ``self._basepoint``.

        INPUT:

        - ``divisor`` -- list. A list with each entry a tuple of the form ``(v, P)``,
          where ``v`` is the valuation of the divisor at point ``P``, ``P`` as per
          the input to :meth:`_aj_based`.

        - ``verbose`` -- logical (default: ``False``); whether to report the progress
          of the computation, in terms of how many elements of the list ``divisor``
          have been completed

        OUTPUT: a vector of length ``self.genus``

        EXAMPLES:

        We can test that the Abel-Jacobi map between two branchpoints of a
        superelliptic curve of degree `p` is a `p`-torsion point in the Jacobian::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: p = 4
            sage: S = RiemannSurface(y^p-x^4+1, prec=100)
            sage: divisor = [(-1, (-1, 0)), (1, (1, 0))]
            sage: AJ = S.abel_jacobi(divisor)  # long time (15 seconds)
            sage: AJxp = [p*z for z in AJ]  # long time
            sage: bool(S.reduce_over_period_lattice(AJxp).norm()<1e-7)  # long time
            True
        """
        if isinstance(divisor, FunctionFieldDivisor):
            divisor = self.divisor_to_divisor_list(divisor)
        ans = 0
        n = len(divisor)
        for i in range(n):
            v, p = divisor[i]
            if verbose:
                print("starting computation for p = {}".format(p))
            ans += v * self._aj_based(p)
            if verbose:
                print(
                    "Done, {}% complete".format(numerical_approx(100 * (i + 1) / n, 11))
                )
        return ans

    def reduce_over_period_lattice(
        self, vector, method='ip', b=None, r=None, normalised=False
    ):
        r"""
        Reduce a vector over the period lattice.

        Given a vector of length ``self.genus``, this method returns a vector
        in the same orbit of the period lattice that is short. There are two
        possible methods, ``'svp'`` which returns a certified shortest vector,
        but can be much slower for higher genus curves, and ``'ip'``, which is
        faster but not guaranteed to return the shortest vector. In general the
        latter will perform well when the lattice basis vectors are of similar
        size.

        INPUT:

        - ``vector`` -- vector. A vector of length ``self.genus`` to reduce over
          the lattice

        - ``method`` -- string (default: ``'ip'``) specifying the method
          to use to reduce the vector; the options are ``'ip'`` and ``'svp'``

        - ``b`` -- integer (default provided); as for
          :meth:`homomorphism_basis`, and used in its invocation if
          (re)calculating said basis

        - ``r`` -- integer (default: ``b/4``);  as for
          :meth:`homomorphism_basis`, and used in its invocation if
          (re)calculating said basis

        - ``normalised`` -- logical (default: ``False``); whether to use the
          period matrix with the differentials normalised s.t. the `A`-matrix
          is the identity.

        OUTPUT:

        Complex vector of length ``self.genus`` in the same orbit as ``vector``
        in the lattice.

        EXAMPLES:

        We can check that the lattice basis vectors themselves are reduced to
        zero::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2 - x^5 + 1)
            sage: epsilon = S._RR(2)^(-S._prec+1)
            sage: for vector in S.period_matrix().columns():
            ....:     print(bool(S.reduce_over_period_lattice(vector).norm()<epsilon))
            True
            True
            True
            True

        We can also check that the method ``'svp'`` always gives a smaller norm
        than ``'ip'``::

            sage: for vector in S.period_matrix().columns():
            ....:     n1 = S.reduce_over_period_lattice(vector).norm()
            ....:     n2 = S.reduce_over_period_lattice(vector, method='svp').norm()
            ....:     print(bool(n2<=n1))
            True
            True
            True
            True
        """
        if not len(vector) == self.genus:
            raise ValueError("Input vector needs to be of length {}".format(self.genus))

        VR = VectorSpace(self._RR, 2 * self.genus)
        VC = VectorSpace(self._CC, self.genus)
        I = self._CC(0, 1)
        PM = self.period_matrix()

        if normalised:
            AM = PM[:, 0 : self.genus]
            AInv = numerical_inverse(AM)
            PM = AInv * PM

        if method == "svp":
            H = max(
                max(z.real_part().abs() for z in vector),
                max(z.imag_part().abs() for z in vector),
            )
            if b is None:
                b = self._prec - 5 - H.log2().floor()
            if r is None:
                r = b // 4
            S = 2**b
            if H * S > 2**(self._prec - 4):
                raise ValueError("insufficient precision for b=%s" % b)

            def C2Z(v):
                vR = [(S * z.real_part()).round() for z in v]
                vR += [(S * z.imag_part()).round() for z in v]
                return vR

            M = Matrix(
                ZZ, 2 * self.genus, 2 * self.genus, [C2Z(c) for c in PM.columns()]
            )
            u = C2Z(vector)
            L = IntegerLattice(M)
            u = VR(u) - VR(L.closest_vector(u))
            reduced = VC(
                [self._CC(u[i] + I * u[i + self.genus]) / S for i in range(self.genus)]
            )

        elif method == "ip":

            def C2R(v):
                return VR([z.real_part() for z in v] + [z.imag_part() for z in v])

            u = C2R(vector)
            basis_vecs = [C2R(c) for c in PM.columns()]
            M = Matrix([[ei.dot_product(ej) for ei in basis_vecs] for ej in basis_vecs])
            v_dot_e = VR([u.dot_product(e) for e in basis_vecs])
            coeffs = M.solve_right(v_dot_e)
            u -= sum([t.round() * e for t, e in zip(coeffs, basis_vecs)])
            reduced = VC(
                [self._CC(u[i] + I * u[i + self.genus]) for i in range(self.genus)]
            )
        else:
            raise ValueError("Must give a valid method.")

        return reduced

    def curve(self):
        r"""
        Return the curve from which this Riemann surface is obtained.

        Riemann surfaces explicitly obtained from a curve return that same object.
        For others, the curve is constructed and cached, so that an identical curve is
        returned upon subsequent calls.

        OUTPUT: curve from which Riemann surface is obtained

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: C = Curve( y^3+x^3-1)
            sage: S = C.riemann_surface()
            sage: S.curve() is C
            True
        """
        if self._curve is None:
            self._curve = Curve(self.f)
        return self._curve

    def places_at_branch_locus(self):
        r"""
        Return the places above the branch locus.

        Return a list of the of places above the branch locus. This must be
        done over the base ring, and so the places are given in terms of the
        factors of the discriminant. Currently, this method only works when
        ``self._R.base_ring() == QQ`` as for other rings, the function field
        for ``Curve(self.f)`` is not implemented. To go from these divisors to
        a divisor list, see :meth:`divisor_to_divisor_list`.

        OUTPUT:

        List of places of the functions field ``Curve(self.f).function_field()``.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(25*(x^4+y^4+1) - 34*(x^2*y^2+x^2+y^2))
            sage: S.places_at_branch_locus()
            [Place (x - 2, (x - 2)*y, y^2 - 17/5, y^3 - 17/5*y),
             Place (x + 2, (x + 2)*y, y^2 - 17/5, y^3 - 17/5*y),
             Place (x - 1/2, (x - 1/2)*y, y^2 - 17/20, y^3 - 17/20*y),
             Place (x + 1/2, (x + 1/2)*y, y^2 - 17/20, y^3 - 17/20*y),
             Place (x^4 - 34/25*x^2 + 1, y, y^2, y^3),
             Place (x^4 - 34/25*x^2 + 1, (x^4 - 34/25*x^2 + 1)*y, y^2 - 34/25*x^2 - 34/25, y^3 + (-34/25*x^2 - 34/25)*y)]
        """
        BP = []
        K = self._R.base_ring()
        if K is not QQ:
            raise NotImplementedError
        C = self.curve()
        KC = C.function_field()
        g0, g1 = self._R.gens()
        Kb = FunctionField(K, str(g0))
        MO = Kb.maximal_order()
        BP = []
        for x in self._discriminant.factor():
            fac = x[0](g0, 0)
            p0 = MO.ideal(fac).place()
            BP += KC.places_above(p0)
        return BP

    def strong_approximation(self, divisor, S):
        r"""
        Apply the method of strong approximation to a divisor.

        As described in [Neu2018]_, apply the method of strong approximation to
        ``divisor`` with list of places to avoid ``S``. Currently, this method
        only works when ``self._R.base_ring() == QQ`` as for other rings, the function
        field for ``Curve(self.f)`` is not implemented.

        INPUT:

        - ``divisor`` -- an element of ``Curve(self.f).function_field().divisor_group()``

        - ``S`` -- list of places to avoid

        OUTPUT:

        A tuple ``(D, B)``, where ``D`` is a new divisor, linearly equivalent
        to ``divisor``, but not intersecting ``S``, and ``B`` is a list of tuples
        ``(v, b)`` where ``b`` are the functions giving the linear equivalence,
        added with multiplicity ``v``.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2-x^3+1)
            sage: avoid = Curve(S.f).places_at_infinity()
            sage: D = 1*avoid[0]
            sage: S.strong_approximation(D, avoid)
            (- Place (x - 2, (x - 2)*y)
              + Place (x - 1, y)
              + Place (x^2 + x + 1, y),
             [(1, (1/(x - 2))*y)])
        """
        # One would standardly expect to run this with
        #     S = Curve(self.f).places_at_infinity()
        # or
        #     S = Curve(self.f).places_at_infinity()+self.places_at_branch_locus()
        #
        # To avoid current implementation issues with going between divisors
        # and divisor lists, we implement a method that handles only divisors
        K = self._R.base_ring()
        if not K == QQ:
            raise NotImplementedError
        C = self.curve()
        KC = C.function_field()
        g0, g1 = self._R.gens()
        Kb = FunctionField(K, str(g0))
        MO = Kb.maximal_order()

        D_base = -sum(S)

        rr = self._vertices[self._basepoint[0]].real()
        rr = rr.ceil()
        Fac = g0 - K(rr)
        p0 = MO.ideal(Fac).place()
        q0 = KC.places_above(p0)[0]

        new_divisor = divisor
        B = []
        for p in divisor.support():
            if p in S:
                v = divisor.valuation(p)
                i = S.index(p)
                Q = S[i]
                D = D_base + Q
                if not D:
                    ios = self.genus
                else:
                    ios = len(D.basis_differential_space())
                while ios > 0:
                    D += q0
                    ios = len(D.basis_differential_space())
                LD = D.function_space()
                V = LD[0]
                a = LD[1]
                b = 0
                for s in S:
                    LDps = (D + s).function_space()
                    Vps = LDps[0]
                    ebd = [LDps[2](a(g)) for g in V.gens()]
                    U = Vps.span(ebd)
                    Quot = Vps.quotient(U)
                    bs = LDps[1](Quot.lift(Quot.basis()[0]))
                    b += bs
                B.append((v, b))
                new_divisor += v * b.divisor()
        return new_divisor, B

    def divisor_to_divisor_list(self, divisor, eps=None):
        r"""
        Turn a divisor into a list for :meth:`abel_jacobi`.

        Given ``divisor`` in ``Curve(self.f).function_field().divisor_group()``,
        consisting of places above finite points in the base, return an equivalent
        divisor list suitable for input into :meth:`abel_jacboi`.

        INPUT:

        - ``divisor`` -- an element of ``Curve(self.f).function_field().divisor_group()``
        - ``eps`` -- real number (optional); tolerance used to determine whether a complex
          number is close enough to a root of a polynomial

        OUTPUT:

        A list with elements of the form ``(v, (z, w))`` representing the finite places.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2-x^3+1)
            sage: D = sum(S.places_at_branch_locus())
            sage: S.divisor_to_divisor_list(D)
            [(1, (1.00000000000000, 0.000000000000000)),
             (1, (-0.500000000000000 - 0.866025403784439*I, 0.000000000000000)),
             (1, (-0.500000000000000 + 0.866025403784439*I, 0.000000000000000))]

        .. TODO::

            Currently this method can only handle places above finite points in
            the base. It would be useful to extend this to allow for places at
            infinity.
        """
        # If this error bound is too restrictive, this method might fail and
        # not return. One might want to change the way this error is handled.
        if not eps:
            eps = self._RR(2)**(-self._prec + 3)
        dl = []

        PZ = PolynomialRing(self._R.base(), "z").fraction_field()
        RF = PolynomialRing(PZ, "w")

        for d in divisor.support():
            if d.is_infinite_place():
                raise NotImplementedError(
                    "Conversion of infinite places not implemented yet."
                )
            v = divisor.valuation(d)
            gs = d._prime.gens()

            g0 = self._R(gs[0])
            gis = [
                sum([PZ(gi.list()[i]) * RF.gen()**i for i in range(len(gi.list()))])
                for gi in gs[1:]
            ]

            rs = self._CCz(g0).roots()
            rys = []

            for r, m in rs:
                ys = []
                for gi in gis:
                    # This test is a bit clunky, it surely can be made more efficient.
                    if len(ys):
                        ers = min([gi(y, r).abs() for y in ys])
                    else:
                        ers = 1

                    if not ers <= eps:
                        poly = self._CCw(gi(self._CCw.gen(0), r))
                        if poly == 0:
                            nys = []
                        else:
                            nys = poly.roots()
                        ys.extend(ny[0] for ny in nys)
                        rys.extend((v * m * n, (r, y)) for y, n in nys)

            if rys:
                dl.extend(rys)
            else:
                for r, m in rs:
                    ys = self._CCw(self.f(r, self._CCw.gen(0))).roots()
                    dl.extend([(v * m * n, (r, y)) for y, n in ys])
        if not sum([v[0] for v in dl]) == divisor.degree():
            raise ValueError(
                "numerical instability, list of wrong degree, returning list {}".format(
                    dl
                )
            )
        return dl


def integer_matrix_relations(M1, M2, b=None, r=None):
    r"""
    Determine integer relations between complex matrices.

    Given two square matrices with complex entries of size `g`, `h` respectively,
    numerically determine an (approximate) `\ZZ`-basis for the `2g \times 2h` matrices
    with integer entries of the shape `(D, B; C, A)` such that `B+M_1*A=(D+M_1*C)*M2`.
    By considering real and imaginary parts separately we obtain `2gh` equations
    with real coefficients in `4gh` variables. We scale the coefficients by a
    constant `2^b` and round them to integers, in order to obtain an integer
    system of equations. Standard application of LLL allows us to determine near
    solutions.

    The user can specify the parameter `b`, but by default the system will
    choose a `b` based on the size of the coefficients and the precision with
    which they are given.

    INPUT:

    - ``M1`` -- square complex valued matrix

    - ``M2`` -- square complex valued matrix of same size as ``M1``

    - ``b`` -- integer (default provided); the equation coefficients are scaled
      by `2^b` before rounding to integers

    - ``r`` -- integer (default: ``b/4``); the vectors found by LLL that satisfy
      the scaled equations to within `2^r` are reported as solutions

    OUTPUT:

    A list of `2g \times 2h` integer matrices that, for large enough `r`, `b-r`,
    generate the `\ZZ`-module of relevant transformations.

    EXAMPLES::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import integer_matrix_relations
        sage: M1 = M2 = matrix(CC, 2, 2, [CC(d).sqrt() for d in [2,-3,-3,-6]])
        sage: T = integer_matrix_relations(M1,M2)
        sage: id = parent(M1)(1)
        sage: M1t = [id.augment(M1) * t for t in T]
        sage: [((m[:,:2]^(-1)*m)[:,2:]-M2).norm() < 1e-13 for m in M1t]
        [True, True]
    """
    if not (M1.is_square() and M2.is_square()):
        raise ValueError("matrices need to be square")
    prec = min(M1.base_ring().precision(), M2.base_ring().precision())
    H = max(
        max(abs(m.real_part()) for m in M1.list() + M2.list()),
        max(abs(m.imag_part()) for m in M1.list() + M2.list()),
    )
    if b is None:
        b = prec - 5 - H.log2().floor()
    if r is None:
        r = b // 4
    S = 2**b
    if H * S > 2**(prec - 4):
        raise ValueError("insufficient precision for b=%s" % b)
    g1 = M1.ncols()
    g2 = M2.ncols()
    CC = (
        M1.base_ring()
        if (M1.base_ring().precision() <= M2.base_ring().precision())
        else M2.base_ring()
    )
    V = ["%s%s" % (n, i) for n in ["a", "b", "c", "d"] for i in range(1, 1 + g1 * g2)]
    R = PolynomialRing(CC, V)
    vars = R.gens()
    A = Matrix(R, g1, g2, vars[: g1 * g2])
    B = Matrix(R, g1, g2, vars[g1 * g2 : 2 * g1 * g2])
    C = Matrix(R, g1, g2, vars[2 * g1 * g2 : 3 * g1 * g2])
    D = Matrix(R, g1, g2, vars[3 * g1 * g2 : 4 * g1 * g2])
    W = ((M1 * A + B) - (M1 * C + D) * M2).list()
    vars = R.gens()
    mt = Matrix(ZZ, [[1 if i == j else 0 for j in range(4 * g1 * g2)] +
      [(S * w.monomial_coefficient(vi).real_part()).round() for w in W] +
      [(S * w.monomial_coefficient(vi).imag_part()).round() for w in W]
                     for i, vi in enumerate(vars)])
    # we compute an LLL-reduced basis of this lattice:
    mtL = mt.LLL()

    def vectomat(v):
        A = Matrix(g1, g2, v[: g1 * g2].list())
        B = Matrix(g1, g2, v[g1 * g2 : 2 * g1 * g2].list())
        C = Matrix(g1, g2, v[2 * g1 * g2 : 3 * g1 * g2].list())
        D = Matrix(g1, g2, v[3 * g1 * g2 : 4 * g1 * g2].list())
        return D.augment(B).stack(C.augment(A))

    c = 2**r
    return [vectomat(v) for v in mtL if all(a.abs() <= c for a in v[g1 * g2 :])]


class RiemannSurfaceSum(RiemannSurface):
    r"""
    Represent the disjoint union of finitely many Riemann surfaces.

    Rudimentary class to represent disjoint unions of Riemann surfaces. Exists
    mainly (and this is the only functionality actually implemented) to
    represents direct products of the complex tori that arise as analytic
    Jacobians of Riemann surfaces.

    INPUT:

    - ``L`` -- list of RiemannSurface objects

    EXAMPLES::

        sage: _.<x> = QQ[]
        sage: SC = HyperellipticCurve(x^6-2*x^4+3*x^2-7).riemann_surface(prec=60)
        sage: S1 = HyperellipticCurve(x^3-2*x^2+3*x-7).riemann_surface(prec=60)
        sage: S2 = HyperellipticCurve(1-2*x+3*x^2-7*x^3).riemann_surface(prec=60)
        sage: len(SC.homomorphism_basis(S1+S2))
        2
    """

    def __init__(self, L):
        r"""
        TESTS::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y> = QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S2 = RiemannSurface(y^2-x^3-x-5)
            sage: S = RiemannSurfaceSum([S1,S2])
            sage: S.riemann_matrix() == S1.riemann_matrix().block_sum(S2.riemann_matrix())
            True
        """
        if not all(isinstance(l, RiemannSurface) for l in L):
            raise ValueError("summands must be RiemannSurface objects")
        prec = min(l._prec for l in L)
        self._prec = prec
        self.genus = sum(s.genus for s in L)
        it = iter(L)
        s = next(it)
        g = s.genus
        PM = s.period_matrix()
        PM1 = PM[:g, :g]
        PM2 = PM[:g, g : 2 * g]
        tau = s.riemann_matrix()
        for s in it:
            g = s.genus
            PM = s.period_matrix()
            PM1 = PM1.block_sum(PM[:g, :g])
            PM2 = PM2.block_sum(PM[:g, g : 2 * g])
            tau = tau.block_sum(s.riemann_matrix())
        self.PM = block_matrix([[PM1, PM2]], subdivide=False)
        self.tau = tau

    def period_matrix(self):
        r"""
        Return the period matrix of the surface.

        This is just the diagonal block matrix constructed from the period
        matrices of the constituents.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y> = QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S2 = RiemannSurface(y^2-x^3-x-5)
            sage: S = RiemannSurfaceSum([S1,S2])
            sage: S1S2 = S1.period_matrix().block_sum(S2.period_matrix())
            sage: S.period_matrix() == S1S2[[0,1],[0,2,1,3]]
            True
        """
        return self.PM

    def riemann_matrix(self):
        r"""
        Return the normalized period matrix of the surface.

        This is just the diagonal block matrix constructed from the Riemann
        matrices of the constituents.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y> = QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S2 = RiemannSurface(y^2-x^3-x-5)
            sage: S = RiemannSurfaceSum([S1,S2])
            sage: S.riemann_matrix() == S1.riemann_matrix().block_sum(S2.riemann_matrix())
            True
        """
        return self.tau

    def __repr__(self) -> str:
        r"""
        Return string describing Riemann surface sum.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y> = QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S2 = RiemannSurface(y^2-x^3-x-5)
            sage: RiemannSurfaceSum([S1,S2])
            Riemann surface sum with period lattice of rank 4
        """
        return "Riemann surface sum with period lattice of rank " + str(2 * self.genus)

    def __add__(self, other):
        r"""
        Return the disjoint union of the Riemann surface and the other argument.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y> = QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S1+S1+S1
            Riemann surface sum with period lattice of rank 6
        """
        return RiemannSurfaceSum([self, other])
