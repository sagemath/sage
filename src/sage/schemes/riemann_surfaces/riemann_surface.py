r"""
Computation of Riemann matrices and endomorphism rings of algebraic Riemann surfaces.

This module provides a class, RiemannSurface, to model the Riemann surface
determined by a plane algebraic curve over a subfield of the complex numbers.

A homology basis is derived from the edges of a Voronoi cell decomposition based on
the branch locus. The pull-back of these edges to the Riemann surface provides
a graph on it that contains a homology basis.

The class provides methods for computing the Riemann period matrix of the
surface numerically, using a certified homotopy continuation method due to
[Kr2016].

The class also provides facilities for computing the endomorphism ring of the
period lattice numerically, by determining integer (near) solutions to the relevant
approximate linear equations.

AUTHORS:

- Alexandre Zotine, Nils Bruin (2017-06-10): initial version

EXAMPLES:

We compute the Riemann matrix of a genus 3 curve::

    sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
    sage: R.<x,y>=QQ[]
    sage: f=x^4-x^3*y+2*x^3+2*x^2*y+2*x^2-2*x*y^2+4*x*y-y^3+3*y^2+2*y+1
    sage: S=RiemannSurface(f,prec=100)
    sage: M=S.riemann_matrix()

We test the usual properties, i.e., that the period matrix is symmetric and that
the imaginary part is positive definite::

    sage: all(abs(a) < 1e-20 for a in (M-M.T).list())
    True
    sage: iM=Matrix(RDF,3,3,[a.imag_part() for a in M.list()])
    sage: iM.is_positive_definite()
    True

We compute the endomorphism ring and check it has `\ZZ`-rank 6::

    sage: A=S.endomorphism_basis(80,8)
    sage: len(A) == 6
    True

In fact it is an order in a number field::

    sage: T.<t>=QQ[]
    sage: K.<a>=NumberField(t^6 - t^5 + 2*t^4 + 8*t^3 - t^2 - 5*t + 7)
    sage: all(len(a.minpoly().roots(K)) == a.minpoly().degree() for a in A)
    True

"""

#*****************************************************************************
#       Copyright (C) 2017 Alexandre Zotine, Nils Bruin
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from scipy.spatial import Voronoi, voronoi_plot_2d
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.complex_field import ComplexField, CDF
from sage.rings.real_mpfr import RealField
from sage.rings.real_double import RDF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.arith.srange import srange
from sage.ext.fast_callable import fast_callable
from sage.graphs.graph import Graph
from sage.matrix.constructor import Matrix
from sage.modules.free_module import VectorSpace
from sage.numerical.gauss_legendre import integrate_vector
from sage.misc.misc_c import prod
import operator

def voronoi_ghost(cpoints, n=6, CC=CDF):
    r"""
    Convert a set of complex points to a list of real tuples `(x,y)`,
    and appends n points in a big circle around them.

    The effect is that, with n >= 3, a Voronoi decomposition will have only
    finite cells around the original points. Furthermore, because
    the extra points are placed on a circle centered on the average of the given
    points, with a radius 3/2 times the largest distance between the center and
    the given points, these finite cells form a simply connected region.

    INPUT:

    - ``cpoints`` -- a list of complex numbers

    OUTPUT:

    A list of real tuples `(x,y)` consisting of the original points and a set
    of points which surround them.

    EXAMPLES::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import voronoi_ghost
        sage: L = [1 + 1*I, 1 - 1*I, -1 + 1*I, -1 - 1*I]
        sage: voronoi_ghost(L) # abs tol 1e-6
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
    average = sum(cpoints)/len(cpoints)
    if len(cpoints) == 1:
        radius = 1
    else:
        radius = 3*max(abs(c-average) for c in cpoints)/2
    z = CC.zeta(n)
    extra_points = [average+radius*z**i for i in range(n)]
    return [(c.real_part(),c.imag_part()) for c in cpoints+extra_points]

def bisect(L,t):
    r"""
    Find position in a sorted list using bisection.

    Given a list `L = [(t_0,...),(t_1,...),...(t_n,...)]` with
    increasing t_i, find the index i such that `t_i <= t < t_{i+1}` using bisection.
    The rest of the tuple is available for whatever use required.

    INPUT:

    - ``L`` -- A list of tuples such that the first term of each tuple is a real
      number between 0 and 1. These real numbers must be increasing.

    - ``t`` -- A real number between `t_0` and `t_n`.

    OUTPUT:

    An integer i, giving the position in L where t would be in

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
    while (min < max-1):
        # Bisect.
        mid = (max+min)//2
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

class RiemannSurface(object):
    r"""
    Construct a Riemann Surface. This is specified by the zeroes of a bivariate
    polynomial with rational coefficients `f(z,w) = 0`.

    INPUT:

    - ``f`` -- a bivariate polynomial with rational coefficients.
      The surface is interpreted as the covering space of the
      coordinate plane in the first variable.

    - ``prec`` -- the desired precision of computations on the surface
      in bits (default: 53)

    - ``certification`` -- a boolean (default: True) value indicating whether
      homotopy continuation is certified or not. Uncertified homotopy continuation
      can be faster.

    - ``differentials`` -- (default: None). If specified, provides a list of
      polynomials `h` such that `h/(df/dw) dz` is a regular differential on the
      Riemann surface. This is taken as a basis of the regular differentials, so
      the genus is assumed to be equal to the length of this list. The results from
      the homology basis computation are checked against this value. Providing this
      parameter makes the computation independent from Singular. For a nonsingular
      plane curve of degree `d`, an appropriate set is given by the monomials of degree
      up to `d-3`.

    EXAMPLES::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
        sage: R.<z,w> = QQ[]
        sage: f = w^2 - z^3 + 1
        sage: RiemannSurface(f)
        Riemann surface defined by polynomial f = -z^3 + w^2 + 1 = 0, with 53 bits of precision

    Another Riemann surface with 100 bits of precision::

        sage: S = RiemannSurface(f, prec=100); S
        Riemann surface defined by polynomial f = -z^3 + w^2 + 1 = 0, with 100 bits of precision
        sage: S.riemann_matrix() #abs tol 0.00000001
        [0.500000000000000000000000... + 0.866025403784438646763723...*I]

    We can also work with Riemann surfaces that are defined over fields with a
    complex embedding, but since the current interface for computing genus and
    regular differentials in Singular presently does not support extensions of QQ,
    we need to specify a description of the differentials ourselves. We give an
    example of a CM elliptic curve::

        sage: Qt.<t> = QQ[]
        sage: K.<a> = NumberField(t^2-t+3,embedding=CC(0.5+1.6*I))
        sage: R.<x,y> = K[]
        sage: f = y^2+y-(x^3+(1-a)*x^2-(2+a)*x-2)
        sage: S = RiemannSurface(f,prec=100,differentials=[1])
        sage: A = S.endomorphism_basis()
        sage: len(A)
        2
        sage: all( len(T.minpoly().roots(K)) > 0 for T in A)
        True

    TESTS:

    This elliptic curve has a relatively poorly conditioned set of branch points,
    so it challenges the path choice a bit. The code just verifies that the period is quadratic,
    because the curve has CM, but really the test is that the computation completes at all.::

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
    def __init__(self, f, prec=53, certification=True, differentials=None):
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
        self._R = f.parent()
        if len(self._R.gens()) != 2:
            raise ValueError('only bivariate polynomials supported.')
        if f.degree() <= 1:
            raise ValueError('equation must be of degree at least 2.')
        z, w = self._R.gen(0), self._R.gen(1)
        self._CC = ComplexField(self._prec)
        self._RR = RealField(self._prec)
        self._CCz = PolynomialRing(self._CC, [self._R.gen(0)])
        self._CCw = PolynomialRing(self._CC, [self._R.gen(1)])
        self.f = f
        if differentials is not None:
            self._differentials = [self._R(a) for a in differentials]
            self.genus = len(self._differentials)
        else:
            self._differentials = None
            self.genus = self._R.ideal(self.f).genus()
            if self.genus < 0:
                raise ValueError("Singular reports negative genus. Specify differentials manually.")
        self.degree = self.f.degree(w)
        self._dfdw = self.f.derivative(w)
        self._dfdz = self.f.derivative(z)
        self._discriminant = self.f.resultant(self._dfdw,w)
        # Coefficients of the polynomial for use in homotopy continuation.
        self._a0 = self._CCz(self.f.coefficient({w:self.degree})(self._CCz.gen(),0))
        self._a0roots = self._a0.roots(multiplicities=False)
        self._aks = [self._CCz(self.f.coefficient({w:self.degree - k - 1})
                               (self._CCz.gen(),0)) for k in range(self.degree)]
        # Compute the branch locus. Takes the square-free part of the discriminant
        # because of numerical issues.
        self.branch_locus = []
        for x in self._discriminant.factor():
            self.branch_locus += self._CCz(x[0](self._CCz.gen(),0)).roots(multiplicities=False)
        # Voronoi diagram and the important points associated with it
        self.voronoi_diagram = Voronoi(voronoi_ghost(self.branch_locus,CC=self._CC))
        self._vertices = [self._CC(x0,y0) for x0,y0 in self.voronoi_diagram.vertices]
        self._wvalues = [self.w_values(z0) for z0 in self._vertices]
        self._Sn = SymmetricGroup(srange(self.degree))
        self._L = dict()
        self._PM = None
        self._fastcall_f = fast_callable(f,domain=self._CC)
        self._fastcall_dfdw = fast_callable(self._dfdw,domain=self._CC)
        self._fastcall_dfdz = fast_callable(self._dfdz,domain=self._CC)

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
        s = 'Riemann surface defined by polynomial f = %s = 0, with %s bits of precision'%(self.f, self._prec)
        return s

    def w_values(self, z0):
        r"""
        Returns the points lying on the surface above ``z0``.

        INPUT:

        - ``z0`` -- (complex) a point in the complex z-plane.

        OUTPUT:

        A set of complex numbers corresponding to solutions of `f(z_0,w) = 0`.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)

        Find the w-values above the origin, i.e. the solutions of `w^2 + 1 = 0`::

            sage: S.w_values(0) # abs tol 1e-14
            [-1.00000000000000*I, 1.00000000000000*I]

        """
        return self.f(z0,self._CCw.gen(0)).roots(multiplicities=False)

    @cached_method
    def downstairs_edges(self):
        r"""
        Compute the edgeset of the Voronoi diagram.

        OUTPUT:

        A list of integer tuples corresponding to edges between vertices
        in the Voronoi diagram.

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
        # The regions of these points are all of the edges which don't go off
        # to infinity, which are exactly the ones we want.
        n = len(self.branch_locus)
        desired_edges = [self.voronoi_diagram.regions[self.voronoi_diagram.point_region[i]] for i in range(n)]
        # First construct the edges as a set because the regions will overlap
        # and we don't want to have two of the same edge.
        edges1 = set()
        for c in desired_edges:
            for j in range(len(c)-1):
                edges1.add(frozenset((c[j],c[j+1])))
            edges1.add(frozenset((c[0],c[-1])))
        # Then make it into a list and sort it.
        # The sorting is important - it will make computing the monodromy group
        # MUCH easier.
        # We orient all the edges so that we go from lower to higher
        # numbered vertex for the continuation.
        edges = [(i0,i1) if (i0 < i1) else (i1,i0) for (i0,i1) in edges1]
        edges.sort()
        return edges

    def downstairs_graph(self):
        r"""
        Retun the Voronoi decomposition as a planar graph.

        The result of this routine can be useful to interpret the labelling
        of the vertices.

        OUTPUT:

        The Voronoi decomposition as a graph, with appropriate planar embedding.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: S.downstairs_graph()
            Graph on 11 vertices

        Similarly one can form the graph of the upstairs edges, which is visually
        rather less attractive but can be instructive to verify that a homology
        basis is likely correctly computed.::

            sage: G=Graph(S.upstairs_edges()); G
            Graph on 22 vertices
            sage: G.is_planar()
            False
            sage: G.genus()
            1
            sage: G.is_connected()
            True

        """
        G=Graph(self.downstairs_edges())
        G.set_pos(dict(enumerate([list(v) for v in self._vertices])))
        return G

    def _compute_delta(self, z1, epsilon, wvalues=None):
        r"""
        Compute a delta for homotopy continuation when moving along a path.

        INPUT:

        - ``z1`` -- a complex number in the z-plane

        - ``epsilon`` -- a real number, which is the minimum distance between
          the w-values above ``z1``

        - ``wvalues`` -- a list (default: None). If specified, saves recomputation.

        OUTPUT:

        A real number, which is a step size for moving along a path.

        EXAMPLES:

        Form a Riemann Surface::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)

        Pick a point which lies on the voronoi diagram, and compute an
        appropriate epsilon::

            sage: z1 = S._vertices[0]
            sage: currw = S.w_values(z1)
            sage: n = len(currw)
            sage: epsilon = min([abs(currw[i] - currw[n-j-1]) for i in range(n) for j in range(n-i-1)])/3
            sage: S._compute_delta(z1, epsilon) # abs tol 1e-8
            0.152628501142363

        If the Riemann surface doesn't have certified homotopy continuation,
        then the delta will just be the minimum distance away from a branch
        point::

            sage: T = RiemannSurface(f, certification=False)
            sage: z1 = T._vertices[0]
            sage: currw = T.w_values(z1)
            sage: n = len(currw)
            sage: epsilon = min([abs(currw[i] - currw[n-j-1]) for i in range(n) for j in range(n-i-1)])/3
            sage: T._compute_delta(z1, epsilon) # abs tol 1e-8
            0.381881307912987

        """
        if self._certification:
            if wvalues is None:
                wvalues = self.w_values(z1)
            # For computation of rho. Need the branch locus + roots of a0.
            badpoints = self.branch_locus + self._a0roots
            rho = min(abs(z1 - z) for z in badpoints)/2
            Y = max(abs(self._fastcall_dfdz(z1,wi)/self._fastcall_dfdw(z1,wi)) for wi in wvalues)

            # compute M
            upperbounds = [sum(ak[k]*(abs(z1) + rho)**k for k in range(ak.degree())) for ak in self._aks]
            upperbounds.reverse()
            # If a0 is a constant polynomial, it is obviously bounded below.
            if self._a0roots == []:
                lowerbound = self._CC(self._a0)/2
            else:
                lowerbound = self._a0[self._a0.degree()]*prod(abs((zk - z1) - rho) for zk in self._a0roots)/2
            M = 2*max(abs((upperbounds[k]/lowerbound))**(1/(k+1)) for k in range(self.degree-1))
            return rho*( ((rho*Y - epsilon)**2 + 4*epsilon*M).sqrt() - (rho*Y + epsilon))/(2*M - 2*rho*Y)
        else:
            # Instead, we just compute the minimum distance between branch
            # points and the point in question.
            return min([abs(b-z1) for b in self.branch_locus])/2

    def homotopy_continuation(self, edge):
        r"""
        Perform homotopy continuation along an edge of the Voronoi diagram
        using Newton iteration.

        INPUT:

        - ``edge`` -- a tuple of integers indicating an edge of the Voronoi
          diagram

        OUTPUT:

        A list of complex numbers corresponding to the points which are reached
        when traversing along the direction of the edge. The ordering of these
        points indicates how they have been permuted due to the weaving of the
        curve.

        EXAMPLES:

        We check that continued values along an edge correspond (up to the appropriate
        permutation) to what is stored. Note that the permutation was originally
        computed from this data.

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f)
            sage: edge1 = S.edge_permutations().keys()[0]
            sage: sigma = S.edge_permutations()[edge1]
            sage: continued_values = S.homotopy_continuation(edge1)
            sage: stored_values = S.w_values(S._vertices[edge1[1]])
            sage: all( abs(continued_values[i]-stored_values[sigma(i)]) < 1e-8 for i in range(3))
            True

        """
        i0, i1 = edge
        ZERO = self._RR.zero()
        ONE = self._RR.one()
        datastorage = []
        z_start = self._CC(self._vertices[i0])
        z_end = self._CC(self._vertices[i1])
        path_length = abs(z_end - z_start)
        def path(t):
            return z_start*(1-t) + z_end*t
        # Primary procedure.
        T = ZERO
        currw = self.w_values(path(T))
        n = len(currw)
        epsilon = min([abs(currw[i] - currw[j]) for i in range(1,n) for j in range(i)])/3
        datastorage += [(T,currw,epsilon)]
        while T < ONE:
            delta = self._compute_delta(path(T), epsilon, wvalues=currw)/path_length
            # Move along the path by delta.
            T += delta
            # If T exceeds 1, just set it to 1 and compute.
            if T > ONE:
                delta -= (T-ONE)
                T = ONE
            while True:
                try:
                    neww = self._determine_new_w(path(T),currw,epsilon)
                except ConvergenceError:
                    delta /= 2
                    T -= delta
                else:
                    break
            currw = neww
            epsilon = min([abs(currw[i] - currw[j]) for i in range(1,n) for j in range(i)])/3
            datastorage += [(T,currw,epsilon)]
        self._L[edge] = datastorage
        return currw

    def _determine_new_w(self, z0, oldw, epsilon):
        r"""
        A procedure to Newton iterate a list of w-values simultaneously.

        Used primarily for moving along the surface for integration or
        homotopy continuation.

        INPUT:

        - ``z0`` -- a complex number

        - ``oldw`` -- a list of w-values which are presumed to be guesses of
          the w-values above ``z0``.

        - ``epsilon`` -- the minimum distance between the points of ``oldw``
          divided by 3

        OUTPUT:

        A list of points the same length as ``oldw`` corresponding to the new
        newton iterated points.

        However, if the newton iteration exceedes the alloted attempts, or
        exits the ``epsilon`` ball, raises a convergence error.

        EXAMPLES:

        First, a trivial example where we guess exactly what the roots are::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: z0 = S._vertices[0]
            sage: epsilon = 0.1
            sage: oldw = S.w_values(z0)
            sage: neww = S._determine_new_w(z0,oldw,epsilon); neww #abs tol 0.00000001
                [-0.934613146929672 + 2.01088055918363*I,
                 0.934613146929672 - 2.01088055918363*I]

        Which should be exactly the same as the w-values we started with.::

            sage: abs(neww[0] - oldw[0]) #abs tol 0.00000001
            0.000000000000...
            sage: abs(neww[1] - oldw[1]) #abs tol 0.00000001
            0.000000000000...

        Here is an example where we exit the ``epsilon`` bound. This approach
        is based on the homotopy continuation procedure which traverses along
        a path and attempts newton iteration::

            sage: g = z^3*w + w^3 + z
            sage: T = RiemannSurface(g)
            sage: z0 = T._vertices[2]*(0.9) - T._vertices[15]*(0.1)
            sage: epsilon = 0.5
            sage: oldw = T.w_values(T._vertices[2])
            sage: T._determine_new_w(z0,oldw,epsilon)
            [-0.562337685361648 + 0.151166007149998*I,
             0.640201585779414 - 1.48567225836436*I,
             -0.0778639004177661 + 1.33450625121437*I]

        """
        # Tools of newton iteration.
        F = self._fastcall_f
        dF = self._fastcall_dfdw
        neww = []
        prec = self._CC.prec()
        # Iterate over all roots.
        for i in range(len(oldw)):
            delta = F(z0,oldw[i])/dF(z0,oldw[i])
            Ndelta = delta.norm()
            wi = oldw[i]-delta
            #it is possible in theory that Newton iteration fails to converge
            #without escaping. We catch this by capping the number of iterations
            #by 100
            for j in range(100):
                # If we exceed the epsilon bound from homotopy continuation,
                # terminate.
                if abs(wi - oldw[i]) >= epsilon:
                    raise ConvergenceError("Newton iteration escaped neighbourhood")
                new_delta = F(z0,wi)/dF(z0,wi)
                Nnew_delta = new_delta.norm()
                # If we found the root exactly, or if delta only affects half the digits and
                # stops getting smaller, we decide that we have converged.
                if (new_delta == 0) or (Nnew_delta>=Ndelta and
                            Ndelta.sign_mantissa_exponent()[2]+prec < wi.norm().sign_mantissa_exponent()[2]):
                    neww.append(wi)
                    break
                delta=new_delta
                Ndelta=Nnew_delta
                wi-=delta
            # If we run 100 iterations without a result, terminate.
            else:
                raise ConvergenceError("Newton iteration fails to converge after %s iterations" % j)
        return neww

    def _newton_iteration(self, z0, oldw, epsilon):
        r"""
        A non-vectorized Newton iteration procedure used for integration.

        INPUT:

        - ``z0`` -- a complex number.

        - ``oldw`` -- a w-value which is presumed to be a guess of one of
          the w-values above ``z0``.

        - ``epsilon`` -- the minimum distance between the w-values divided by 3.

        OUTPUT:

        A complex number, which should be a w-value above ``z0``.

        However, if the Newton iteration exceedes the alloted attempts, or
        exits the ``epsilon`` ball, raises a convergence error.

        EXAMPLES:

        First, a trivial example where we guess exactly what the root is::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: z0 = S._vertices[0]
            sage: epsilon = 0.1
            sage: oldw = S.w_values(z0)[0]
            sage: neww = S._newton_iteration(z0,oldw,epsilon); neww #abs tol 0.00000001
            -0.934613146929672 + 2.01088055918363*I

        Which should be exactly the same as the w-value we started with::

            sage: oldw - neww #abs tol 0.00000001
            0.000000000000000

        Here is an example where we exit the epsilon bound. This approach
        is based on the homotopy continuation procedure which traverses along
        a path and attempts newton iteration::

            sage: g = z^3*w + w^3 + z
            sage: T = RiemannSurface(g)
            sage: z0 = T._vertices[2]*(0.9) - T._vertices[15]*(0.1)
            sage: epsilon = 0.5
            sage: oldw = T.w_values(T._vertices[2])[0]
            sage: T._newton_iteration(z0, oldw, epsilon)
            -0.562337685361648 + 0.151166007149998*I

        """
        F = self._fastcall_f
        dF = self._fastcall_dfdw
        prec = self._CC.prec()
        delta = F(z0,oldw)/dF(z0,oldw)
        Ndelta = delta.norm()
        neww = oldw-delta
        eps_squared = epsilon**2
        #it is possible in theory that Newton iteration fails to converge
        #without escaping. We catch this by capping the number of iterations
        #by 100
        for j in range(100):
            if (neww-oldw).norm() > eps_squared:
                raise ConvergenceError("Newton iteration escaped neighbourhood")
            new_delta = F(z0,neww)/dF(z0,neww)
            Nnew_delta = new_delta.norm()
            # If we found the root exactly, or if delta only affects half the digits and
            # stops getting smaller, we decide that we have converged.
            if (new_delta == 0) or (Nnew_delta>=Ndelta and
                    Ndelta.sign_mantissa_exponent()[2]+prec < neww.norm().sign_mantissa_exponent()[2]):
                return neww
            delta = new_delta
            Ndelta = Nnew_delta
            neww-=delta
        raise ConvergenceError("Newton iteration fails to converge")

    @cached_method
    def upstairs_edges(self):
        r"""
        Compute the edgeset of the lift of the downstairs graph onto the
        Riemann surface.

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
            # Epsilon for checking w-value later.
            epsilon = min([abs(self._wvalues[i1][i] - self._wvalues[i1][n-j-1]) for i in range(n) for j in range(n-i-1)])/3
            # Homotopy continuation along e.
            homotopycont = self.homotopy_continuation(e)
            for i in range(len(homotopycont)):
                # Checks over the w-values of the next point to check which it is.
                for j in range(len(self._wvalues[i1])):
                    if abs(homotopycont[i] - self._wvalues[i1][j]) < epsilon:
                        # Once it finds the appropriate w-value, adds the edge.
                        edgeset = edgeset + [[(i0,i),(i1,j)]]
                        continue
        return edgeset

    def _edge_permutation(self, edge):
        r"""
        Compute the permutation of the w-values above a point in the z-plane
        when moving along an edge in the Voronoi diagram.

        INPUT:

        - ``edge`` -- an edge on the Voronoi diagram

        OUTPUT:

        A permutation corresponding to how the roots interchange when
        moving along the edge.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f)

        Compute the edge permutation of (1,2) on the Voronoi diagram::

            sage: S._edge_permutation((1,2))
            (0,2,1)

        This indicates that while traversing along the direction of `(5,16)`,
        the 2nd and 3rd layers of the Riemann surface are interchanging.
        """
        if edge in self.downstairs_edges():
            #find all upstairs edges that are lifts of the given downstairs edge
            #and store the corresponding indices at start and end that label the
            #branches upstairs.
            L = [(j0,j1) for ((i0,j0),(i1,j1)) in self.upstairs_edges() if edge==(i0,i1)]
            #we should be finding exactly "degree" of these
            assert len(L) == self.degree
            #and as a corollary of how we construct them, the indices at the start
            #should be in order
            assert all(a==b[0] for a,b in enumerate(L))
            return self._Sn([j1 for j0,j1 in L])
        else:
            raise ValueError('edge not in Voronoi diagram')

    @cached_method
    def edge_permutations(self):
        r"""
        Compute the permutations of branches associated to each edge

        Over the vertices of the Voronoi decomposition around the branch
        locus, we label the fibres. By following along an edge, the lifts
        of the edge induce a permutation of that labelling.

        OUTPUT:

        A dictionary with as keys the edges of the Voronoi decomposition
        and as values the corresponding permutations.

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
        D=dict( (e,self._edge_permutation(e)) for e in self.downstairs_edges())
        for e in list(D.keys()):
            D[(e[1],e[0])]=D[e]**(-1)
        return D

    @cached_method
    def monodromy_group(self):
        r"""
        Compute local monodromy generators of the riemann surface.

        For each branch point, the local monodromy is encoded by a permutation.
        The permutations returned correspond to positively oriented loops around
        each branch point, with a fixed base point. This means the generators are
        properly conjugated to ensure that together they generate the global monodromy.
        The list has an entry for every finite point stored in `self.branch_locus`, plus an entry
        for the ramification above infinity.

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

        The permutations give the local monodromy generators for the branch points::

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
        #we get all the regions
        loops = [self.voronoi_diagram.regions[i][:] for i in self.voronoi_diagram.point_region]
        #and construct their Voronoi centers as complex numbers
        centers = self.branch_locus + [self._CC(x,y) for x,y in self.voronoi_diagram.points[n:]]
        for center, loop in zip(centers,loops):
            if -1 in loop:
                #for loops involving infinity we take the finite part of the path
                i = loop.index(-1)
                loop[:] = loop[i+1:]+loop[:i]
            else:
                #and for finite ones we close the paths
                loop.append(loop[0])
            #we make sure the loops are positively oriented wrt. their center
            v0 = self._vertices[loop[0]]
            v1 = self._vertices[loop[1]]
            M = Matrix([list(v0-center),list(v1-center)])
            if M.det() < 0:
                loop.reverse()

        #we stitch together the paths that are part of loops through
        #infinity. There should be a unique way of doing so.
        inf_loops = loops[n:]
        inf_path = inf_loops.pop()
        while (inf_loops):
            inf_path += (inf_loops.pop())[1:]
        assert inf_path[0] == inf_path[-1]

        loops = loops[:n]
        loops.append(inf_path)

        P0 = loops[0][0]
        monodromy_gens = []
        edge_perms = self.edge_permutations()
        for c in loops:
            to_loop = G.shortest_path(P0,c[0])
            to_loop_perm = reduce(
                operator.mul,
                (edge_perms[(to_loop[i],to_loop[i+1])]
                    for i in range(len(to_loop)-1)),
                self._Sn(()))
            c_perm = reduce(
                operator.mul,
                (edge_perms[(c[i],c[i+1])]
                    for i in range(len(c)-1)),
                self._Sn(()))
            monodromy_gens.append(to_loop_perm*c_perm*to_loop_perm**(-1))
        return monodromy_gens

    @cached_method
    def homology_basis(self):
        r"""
        Compute the homology basis of the Riemann surface.

        OUTPUT:

        A list of paths `L = [P_1, \dots, P_n]`.

        Each path `P_i` is of the form `(k, [p_1 ... p_m, p_1])`, where
        `k` is the number of times to traverse the path (if negative, to
        traverse it backwards), and the `p_i` are vertices of the
        upstairs graph.

        EXAMPLES:

        In this example, there are two paths that form the homology basis::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: g = w^2 - z^4 + 1
            sage: S = RiemannSurface(g)
            sage: S.homology_basis()
            [[(1,
               [(3, 1),
                (5, 0),
                (9, 0),
                (10, 0),
                (2, 0),
                (4, 0),
                (7, 1),
                (10, 1),
                (3, 1)])],
             [(1,
               [(8, 0),
                (6, 0),
                (7, 0),
                (10, 0),
                (2, 0),
                (4, 0),
                (7, 1),
                (10, 1),
                (9, 1),
                (8, 0)])]]

        """
        if self.genus == 0:
            return []

        edgesu = self.upstairs_edges()
        cycles = Graph(edgesu).cycle_basis()
        # Computing the Gram matrix.
        cn = len(cycles)
        # Forming a list of lists of zeroes. Later this will be converted into a
        # matrix.
        intersectionprod = [[0 for c in cycles] for c in cycles]
        # This loop will start at the entry (0,1), and proceed along the row up
        # til (0,cn-1).
        # Then it will go to entry (1,2), and proceed along the row, etc.
        for i in range(1,cn):
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
                    vd = self._vertices[cycles[i][i0][0]]

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

                    a_in=cycles[i][i0-1][0]
                    a_out=cycles[i][(i0+1)%len(cycles[i])][0]
                    b_in=cycles[j][i1-1][0]
                    b_out=cycles[j][(i1+1)%len(cycles[j])][0]

                    # we can get the angles (and hence the rotation order)
                    # by taking the arguments of the differences.

                    a_in_arg=(self._vertices[a_in]-vd).argument()
                    a_out_arg=(self._vertices[a_out]-vd).argument()
                    b_in_arg=(self._vertices[b_in]-vd).argument()
                    b_out_arg=(self._vertices[b_out]-vd).argument()

                    # we make sure to test overlap on the indices, so no rounding
                    # problems occur with that.

                    if (b_in != a_in) and (b_in != a_out):
                        if ((a_in_arg<b_in_arg<a_out_arg) or
                                (b_in_arg<a_out_arg<a_in_arg) or
                                (a_out_arg<a_in_arg<b_in_arg)):
                            intsum += 1
                        elif ((a_out_arg<b_in_arg<a_in_arg) or
                                (b_in_arg<a_in_arg<a_out_arg) or
                                (a_in_arg<a_out_arg<b_in_arg)):
                            intsum -= 1
                        else:
                            raise RuntimeError("impossible edge orientation")
                    if (b_out != a_in) and (b_out != a_out):
                        if ((a_in_arg<b_out_arg<a_out_arg) or
                                (b_out_arg<a_out_arg<a_in_arg) or
                                (a_out_arg<a_in_arg<b_out_arg)):
                            intsum -= 1
                        elif ((a_out_arg<b_out_arg<a_in_arg) or
                                (b_out_arg<a_in_arg<a_out_arg) or
                                (a_in_arg<a_out_arg<b_out_arg)):
                            intsum += 1
                        else:
                            raise RuntimeError("impossible edge orientation")
                assert (intsum%2) == 0
                intsum = intsum//2
                intersectionprod[i][j] = intsum
                # Skew Symmetry
                intersectionprod[j][i] = -intsum
        Gmatrix = Matrix(intersectionprod)
        G_normalized,P = Gmatrix.symplectic_form()
        if G_normalized.rank() != 2*self.genus:
            raise RuntimeError("rank of homology pairing mismatches twice stored genus")
        # Define the cycle sets.
        acycles = [[] for i in range(self.genus)]
        bcycles = [[] for i in range(self.genus)]
        # There are g a and b cycles.
        for i in range(self.genus):
            # Range over the size of the Gram matrix.
            for j in range(cn):
                # Forms the acycles and bcycles. If the entry in the
                # transformation matrix is non-zero, it adds the coefficient at
                # that entry, and the corresponding cycle. (also, forms it
                # into a loop)
                if P[i][j] != 0:
                    acycles[i] += [(P[i][j],[x for x in cycles[j]]+[cycles[j][0]])]
                if P[self.genus + i][j] != 0:
                    bcycles[i] += [(P[self.genus + i][j],[x for x in cycles[j]]+[cycles[j][0]])]
        return acycles + bcycles

    def make_zw_interpolator(self, upstairs_edge):
        r"""
        Given an upstairs edge for which continuation data has been stored,
        return a function that computes `z(t),w(t)` , where t in `[0,1]` is a
        parametrization of the edge.

        INPUT:

        - ``upstairs_edge`` -- a pair of integer tuples indicating an edge on
          the upstairs graph of the surface

        OUTPUT:

        A tuple (g, d), where g is the function that computes the interpolation
        along the edge and d is the difference of the z-values of the end and start point.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: _ = S.homology_basis()
            sage: g,d = S.make_zw_interpolator([(0,0),(1,0)]);
            sage: all(f(*g(i*0.1)).abs() < 1e-13for i in range(10))
            True
            sage: abs((g(1)[0]-g(0)[0]) - d) < 1e-13
            True

        """
        eindex = tuple(u[0] for u in upstairs_edge)
        i0, i1 = eindex
        z_start = self._vertices[i0]
        z_end = self._vertices[i1]
        currL = self._L[eindex]
        windex = upstairs_edge[0][1]
        def w_interpolate(t):
            if t < 0 or t > 1:
                raise ValueError("t outside path range")
            if t == 0:
                return z_start,currL[0][1][windex]
            elif t == 1:
                return z_end,currL[-1][1][windex]
            while True:
                i = bisect(currL,t)
                t1, w1 ,epsilon = currL[i]
                w1 = w1[windex]
                t2, w2, _ = currL[i+1]
                w2 = w2[windex]
                z0 = (1-t)*z_start+t*z_end
                w0 = self._CC(((t2-t)*w1+(t-t1)*w2)/(t2-t1))
                try:
                    desired_result = self._newton_iteration(z0,w0,epsilon)
                except ConvergenceError:
                    pass
                else:
                    return z0,desired_result
                #If we did not succeed, we insert a new point in our interpolation list
                tnew=t
                while True:
                    tnew = (t1 + tnew)/2
                    znew = (1-tnew)*self._vertices[i0]+tnew*self._vertices[i1]
                    try:
                        neww1 = self._determine_new_w(znew,currL[i][1],epsilon)
                    except ConvergenceError:
                        pass
                    else:
                        #When *no* ConvergenceError is raised, we have succeeded and we can exit
                        break
                #once the loop has succeeded we insert our new value
                t1 = tnew
                self._L[eindex].insert(i+1,(t1,neww1,epsilon))
        return w_interpolate,(z_end-z_start)

    def simple_vector_line_integral(self, upstairs_edge, differentials):
        r"""
        Perfom vectorized integration along a straight path.

        INPUT:

        - ``upstairs_edge`` -- a pair of integer tuples corresponding to an
          edge of the upstairs graph.

        - ``differentials`` -- a list of polynomials; a polynomial `g`
          represents the differential `g(z,w)/(df/dw) dz` where `f(z,w)=0` is
          the equation defining the Riemann surface.

        OUTPUT:

        A complex number, the value of the line integral.

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
            sage: S.simple_vector_line_integral([(0,0),(1,0)], differentials) #abs tol 0.00000001
            (1.14590610929717e-16 - 0.352971844594760*I)

        ..NOTE::

            Uses data that "homology_basis" initializes.
        """
        w_of_t,Delta_z = self.make_zw_interpolator(upstairs_edge)
        V = VectorSpace(self._CC,self.genus)
        def integrand(t):
            zt,wt = w_of_t(t)
            dfdwt = self._fastcall_dfdw(zt,wt)
            return V([omega(zt,wt)/dfdwt for omega in differentials])

        I=integrate_vector(integrand,self._prec)*Delta_z
        return I

    def cohomology_basis(self, option=1):
        r"""
        Compute the cohomology basis of this surface.

        INPUT:

        - ``option`` -- Presently, this routine uses Singular's ``adjointIdeal``
            and passes the ``option`` parameter on. Legal values are 1, 2, 3 ,4,
            where 1 is the default. See the Singular documentation for the meaning.
            The backend for this function may change, and support for this parameter may
            disappear.

        OUTPUT:

        Returns a list of polynomials `g` representing the holomorphic
        differentials `g/(df/dw) dz`, where `f(z,w)=0` is the equation
        specifying the Riemann surface.

        EXAMPLES:

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f)
            sage: S.cohomology_basis()
            [1, w, z]
        """
        if self.genus == 0:
            self._differentials = []
            return self._differentials[0]
        if self._differentials is None:
            # Computes differentials from the adjointIdeal using Singular
            # First we homogenize
            base = self.f.base_ring()
            # It's important we use a degree ordering; see below.
            R = self._R
            k = PolynomialRing(base,names="Z,W,U",order="degrevlex")
            dehom = k.Hom(R)([R.gen(0),R.gen(1),R.one()])
            fnew = self.f(k.gen(0)/k.gen(2),k.gen(1)/k.gen(2)).numerator()

            # We load the relevant functionality into singularlib
            import sage.libs.singular.function_factory
            sage.libs.singular.function_factory.lib("paraplanecurves.lib")
            adjointIdeal = sage.libs.singular.function.singular_function("adjointIdeal")
            libsing_options=sage.libs.singular.option.LibSingularVerboseOptions()

            # We compute the adjoint ideal (note we need to silence "redefine")
            redef_save = libsing_options['redefine']
            try:
                libsing_options['redefine'] = False
                J = adjointIdeal(fnew,option)
            finally:
                libsing_options['redefine'] = redef_save

            # We are interested in the (degree-3) subspace of the adjoint ideal.
            # We compute this by intersecting with (Z,W,U)^(degree-3). Then the
            # lowest degree generators are a basis of the relevant subspace.
            d=fnew.total_degree()
            J2 = k.ideal(J).intersection(k.ideal([k.gen(0),k.gen(1),k.gen(2)])**(d-3))
            generators = [dehom(c) for c in J2.gens() if c.degree() == d-3]
            if len(generators) != self.genus:
                raise ValueError("computed regular differentials do not match stored genus")
            self._differentials = generators
        return self._differentials

    def matrix_of_integral_values(self, differentials):
        r"""
        Compute the path integrals of the given differentials along the homology basis.

        The returned answer has a row for each differential. If the Riemann surface is
        given by the equation `f(z,w)=0`, then the differentials are encoded by polynomials
        g, signifying the differential `g(z,w)/(df/dw) dz`.

        INPUT:

        - ``differentials`` -- a list of polynomials.

        OUTPUT:

        A matrix, one row per differential, containing the values of the path integrals along
        the homology basis of the Riemann surface.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(x^3 + y^3 + 1)
            sage: B = S.cohomology_basis()
            sage: S.matrix_of_integral_values(B) #abs tol 1e-12
            [   0.883319375142725 - 1.52995403705719*I 1.76663875028545 + 5.55111512312578e-17*I]

        """
        cycles = self.homology_basis()
        def normalize_pairs(L):
            r"""
            Returns a list of edges encoded by the path in L.
            The edges are normalized to be in the direction in which
            the homotopy continuation should have been computed along them.
            """
            R=[]
            for i in range(len(L)-1):
                if L[i][0]<L[i+1][0]:
                    R.append((L[i],L[i+1]))
                else:
                    R.append((L[i+1],L[i]))
            return R
        occurring_edges = set()
        occurring_edges.update(*[normalize_pairs(p[1]) for h in cycles for p in h])
        integral_dict=dict()
        for upstairs_edge in occurring_edges:
            integral_dict[upstairs_edge]=self.simple_vector_line_integral(upstairs_edge,differentials)
        rows=[]
        for cycle in cycles:
            V = VectorSpace(self._CC,self.genus).zero()
            for multiplicity,loop in cycle:
                for i in range(len(loop)-1):
                    if loop[i][0]<loop[i+1][0]:
                        direction=1
                        upstairs_edge=(loop[i],loop[i+1])
                    else:
                        direction=-1
                        upstairs_edge=(loop[i+1],loop[i])
                    V+=(multiplicity*direction)*integral_dict[upstairs_edge]
            rows.append(V)
        return Matrix(rows).transpose()

    @cached_method
    def period_matrix(self):
        r"""
        Compute the period matrix of the surface.

        OUTPUT:

        A matrix of complex values.

        EXAMPLES:

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f, prec=30)
            sage: M = S.period_matrix()

        The results are highly arbitrary, so it is hard to check if the result produced is
        correct. The closely related `riemann matrix` is somewhat easier to test.

            sage: parent(M)
            Full MatrixSpace of 3 by 6 dense matrices over Complex Field with 30 bits of precision
            sage: M.rank()
            3

        """
        differentials = self.cohomology_basis()
        differentials = [fast_callable(omega,domain=self._CC)
                for omega in self.cohomology_basis()]
        PM=self.matrix_of_integral_values(differentials)
        return PM

    def riemann_matrix(self):
        r"""
        Compute the Riemann matrix.

        OUTPUT:

        A matrix of complex values

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f, prec=60)
            sage: M = S.riemann_matrix()

        The Klein quartic has a Riemann matrix with values is a quadratic
        field::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2-x+2)
            sage: all(len(m.algdep(6).roots(K)) > 0 for m in M.list())
            True

        """
        PeriodMatrix = self.period_matrix()
        Am = PeriodMatrix[0:self.genus,0:self.genus]
        RM = (Am.inverse())*PeriodMatrix[0:self.genus,self.genus:2*self.genus]
        return RM

    def plot_paths(self):
        r"""
        Make a graphical representation of the integration paths.

        Returns a two dimensional plot containing the branch points (in red)
        and the integration paths (obtained from the Voronoi cells of the
        branch points). The integration paths are plotted by plotting the points
        that have been computed for homotopy continuation, so the density
        gives an indication of where numerically sensitive features occur.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2 - x^3 - x)
            sage: S.plot_paths()
            Graphics object consisting of 2 graphics primitives

        """
        from sage.plot.point import point2d
        P=[]

        #trigger the computation of the homology basis, so that self._L is present
        self.homology_basis()

        for e in self._L.keys():
            z0=self._vertices[e[0]]
            z1=self._vertices[e[1]]
            def path(t):
                z=(1-t)*z0+t*z1
                return z
            T=self._L[e]
            P+=[path(t[0]) for t in T]
        plt=point2d(P,size=1)+point2d(self.branch_locus,color="red")
        return plt

    def plot_paths3d(self,thickness=0.01):
        r"""
        Return the homology basis as a graph in 3-space.

        The homology basis of the surface is constructed by taking the Voronoi
        cells around the branch points and taking the inverse image of the
        edges on the Riemann surface. If the surface is given by the equation
        `f(z,w)`, the returned object gives the image of this graph in 3-space
        with coordinates
        `\left(\operatorname{Re}(z), \operatorname{Im}(z), \operatorname{Im}(w)\right)`.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2-x^3-x)
            sage: S.plot_paths3d()
            Graphics3d Object

        """
        from sage.plot.graphics import Graphics
        from sage.plot.plot3d.shapes2 import point3d, line3d
        P = Graphics()

        #trigger the computation of the homology basis, so that self._L is present
        self.homology_basis()

        for e in self._L.keys():
            z0 = self._vertices[e[0]]
            z1 = self._vertices[e[1]]
            def path(t):
                z = (1-t)*z0+t*z1
                return (z.real_part(),z.imag_part())
            T = self._L[e]
            color = "blue"
            for i in range(self.degree):
                P += line3d([path(t[0])+(t[1][i].imag_part(),) for t in T],color=color,thickness=thickness)
            for z,ws in zip(self._vertices,self._wvalues):
                for w in ws:
                    P += point3d([z.real_part(),z.imag_part(),w.imag_part()],color="purple", size=20)
        return P

    def endomorphism_basis(self, b=None, r=None):
        r"""
        Numerically compute a `\ZZ`-basis for the endomorphism ring.

        Let `\left(I | M \right)` be the normalized period matrix (`M` is the
        `g\times g` :meth:`riemann_matrix`).
        We consider the system of matrix equations `MA + C = (MB + D)M` where
        `A, B, C, D` are `g\times g` integer matrices.  We determine small
        integer (near) solutions using LLL reductions.  These solutions are
        returned as `2g \times 2g` integer matrices obtained by stacking
        `\left(D | B\right)` on top of `\left(C | A\right)`.

        INPUT:

        - ``b`` -- integer (default provided). The equation coefficients
          are scaled by `2^b` before rounding to integers.

        - ``r`` -- integer (default: ``b/4``). Solutions that have all
          coefficients smaller than `2^r` in absolute value are reported as
          actual solutions.

        OUTPUT:

        A list of `2g \times 2g` integer matrices that, for large enough ``r``
        and ``b-r``, generate the endomorphism ring.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(x^3 + y^3 + 1)
            sage: S.endomorphism_basis()
            [
            [1 0]  [ 0 -1]
            [0 1], [ 1  1]
            ]

        """
        M = self.riemann_matrix()
        return integer_matrix_relations(M,M,b,r)

    def homomorphism_basis(self, other, b=None, r=None):
        r"""
        Numerically compute a `\ZZ`-basis for module of homomorphisms to
        a given complex torus.

        Given another complex torus (given as the analytic Jacobian of a
        Riemann surface), numerically compute a basis for the homomorphism
        module. The answer is returned as a list of 2g x 2g integer matrices
        T=(D, B; C, A)
        such that if the columns of (I|M1) generate the lattice defining
        the Jacobian of the Riemann surface and the columns of (I|M2) do this
        for the codomain, then approximately we have
        (I|M2)T=(D+M2C)(I|M1),
        i.e., up to a choice of basis for `\CC^g` as a complex vector space, we
        we realize (I|M1) as a sublattice of (I|M2).

        INPUT:

        - ``b`` -- integer (default provided). The equation coefficients
          are scaled by `2^b` before rounding to integers.

        - ``r`` -- integer (default: ``b/4``). Solutions that have all
          coefficients smaller than `2^r` in absolute value are reported as
          actual solutions.

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
        return integer_matrix_relations(M2,M1,b,r)

    def __add__(self,other):
        r"""
        Return the disjoint union of the Riemann surface and the other argument.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y>=QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S1+S1
            Riemann surface sum with period lattice of rank 4

        """
        return RiemannSurfaceSum([self,other])

def integer_matrix_relations(M1,M2,b=None,r=None):
    r"""
    Determine integer relations between complex matrices

    Given two g x g matrices with complex entries, numerically determine
    an (approximate) ZZ-basis for the 2g x 2g matrices with integer entries
    of the shape (D, B; C, A) such that
    B+M1*A=(D+M1*C)*M2
    By considering real and imaginary parts separately we obtain `2g^2`
    equations with real coefficients in `4g^2` variables. We scale the
    coefficients by a constant `2^b` and round them to integers, in order
    to obtain an integer system of equations. Standard application of LLL
    allows us to determine near solutions.

    The user can specify the parameter `b`, but by default the system will
    choose a `b` based on the size of the coefficients and the precision
    with which they are given.

    INPUT:

    - `M1` -- square complex valued matrix

    - `M2` -- square complex valued matrix of same size as M1

    - ``b`` -- integer (default provided). The equation coefficients
      are scaled by `2^b` before rounding to integers.

    - ``r`` -- integer (default: ``b/4``). The vectors found by LLL that satisfy
      the scaled equations to withing `2^r` are reported as solutions.

    OUTPUT:

    A list of 2g x 2g integer matrices that, for large enough `r`, `b-r`,
    generate the ZZ-module of relevant transformations.

    EXAMPLES::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import integer_matrix_relations
        sage: M1=M2=matrix(CC,2,2,[sqrt(d) for d in [2,-3,-3,-6]])
        sage: T=integer_matrix_relations(M1,M2)
        sage: id=parent(M1)(1)
        sage: M1t=[id.augment(M1) * t for t in T]
        sage: [((m[:,:2]^(-1)*m)[:,2:]-M2).norm() < 1e-13 for m in M1t]
        [True, True]

    """
    if not(M1.ncols()==M2.ncols() and M1.nrows()==M1.nrows() and M2.nrows()==M2.nrows()):
        raise ValueError("matrices need to be square of same dimensions")
    prec = min(M1.base_ring().precision(),M2.base_ring().precision())
    H = max(max( abs(m.real_part()) for m in M1.list()+M2.list()), max( abs(m.imag_part()) for m in M1.list()+M2.list()))
    if b is None:
        b = prec-5-H.log2().floor()
    if r is None:
        r = b//4
    S = 2**b
    if H*S > 2**(prec-4):
        raise ValueError("insufficient precision for b=%s"%b)
    g = M1.ncols()
    CC = M1.base_ring()
    V = ["%s%s"%(n,i) for n in ["a","b","c","d"] for i in srange(1,1+g**2)]
    R = PolynomialRing(CC,V)
    A = Matrix(R,g,g,V[:g**2])
    B = Matrix(R,g,g,V[g**2:2*g**2])
    C = Matrix(R,g,g,V[2*g**2:3*g**2])
    D = Matrix(R,g,g,V[3*g**2:4*g**2])
    W = ((M1*A+B) - (M1*C+D)*M2).list()
    vars = R.gens()
    mt = Matrix(ZZ,[[1 if i==j else 0 for j in range(4*g**2)] +
      [(S*w.monomial_coefficient(vars[i]).real_part()).round() for w in W] +
      [(S*w.monomial_coefficient(vars[i]).imag_part()).round() for w in W] for i in range(len(vars))])
    # we compute an LLL-reduced basis of this lattice:
    mtL = mt.LLL()
    def vectomat(v,g):
        A = Matrix(g,g,v[:g**2].list())
        B = Matrix(g,g,v[g**2:2*g**2].list())
        C = Matrix(g,g,v[2*g**2:3*g**2].list())
        D = Matrix(g,g,v[3*g**2:4*g**2].list())
        return D.augment(B).stack(C.augment(A))
    c = 2**r
    return [vectomat(v,g) for v in mtL if all(a.abs() <= c for a in v[4*g**2:])]

class RiemannSurfaceSum(RiemannSurface):
    r"""
    Represent the disjoint union of finitely many Riemann surfaces.

    Rudimentary class to represent disjoint unions of Riemann surfaces.
    Exists mainly (and this is the only functionality actually
    implemented) to represents direct products of the complex tori that
    arise as analytic Jacobians of Riemann surfaces.

    INPUT:

    - L -- list of RiemannSurface objects

    EXAMPLES::

        sage: _.<x> = QQ[]
        sage: SC = HyperellipticCurve(x^6-2*x^4+3*x^2-7).riemann_surface(prec=60)
        sage: S1 = HyperellipticCurve(x^3-2*x^2+3*x-7).riemann_surface(prec=60)
        sage: S2 = HyperellipticCurve(1-2*x+3*x^2-7*x^3).riemann_surface(prec=60)
        sage: len(SC.homomorphism_basis(S1+S2))
        2

    """
    def __init__(self,L):
        r"""
        TESTS::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y>=QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S2 = RiemannSurface(y^2-x^3-x-5)
            sage: S = RiemannSurfaceSum([S1,S2])
            sage: S.riemann_matrix() == S1.riemann_matrix().block_sum(S2.riemann_matrix())
            True

        """
        if not all(isinstance(l,RiemannSurface) for l in L):
            raise ValueError("summands must be RiemannSurface objects")
        prec=min(l._prec for l in L)
        self._prec=prec
        it=iter(L)
        M=next(it).riemann_matrix()
        for s in it:
            M=M.block_sum(s.riemann_matrix())
        self.M=M

    def riemann_matrix(self):
        r"""
        Return the normalized period matrix of the surface

        This is just the diagonal block matrix constructed from the Riemann matrices
        of the constituents.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y>=QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S2 = RiemannSurface(y^2-x^3-x-5)
            sage: S = RiemannSurfaceSum([S1,S2])
            sage: S.riemann_matrix() == S1.riemann_matrix().block_sum(S2.riemann_matrix())
            True

        """
        return self.M

    def __repr__(self):
        r"""
        Return string describing Riemann surface sum.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y>=QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S2 = RiemannSurface(y^2-x^3-x-5)
            sage: RiemannSurfaceSum([S1,S2])
            Riemann surface sum with period lattice of rank 4

        """
        return "Riemann surface sum with period lattice of rank "+repr(2*self.M.ncols())

    def __add__(self,other):
        r"""
        Return the disjoint union of the Riemann surface and the other argument.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y>=QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S1+S1+S1
            Riemann surface sum with period lattice of rank 6

        """
        return RiemannSurfaceSum([self,other])
