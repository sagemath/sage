# sage_setup: distribution = sagemath-schemes
r"""
Dynamical systems for products of projective spaces

This class builds on the product projective space class.
The main constructor functions are given by ``DynamicalSystem`` and
``DynamicalSystem_projective``. The constructors function can take either
polynomials or a morphism from which to construct a dynamical system.

The must be specified.

EXAMPLES::

    sage: P1xP1.<x,y,u,v> = ProductProjectiveSpaces(QQ, [1, 1])
    sage: DynamicalSystem_projective([x^2*u, y^2*v, x*v^2, y*u^2], domain=P1xP1)
    Dynamical System of Product of projective spaces P^1 x P^1 over Rational Field
      Defn: Defined by sending (x : y , u : v) to
            (x^2*u : y^2*v , x*v^2 : y*u^2).
"""
# ****************************************************************************
#     Copyright (C) 2014 Ben Hutz <bn4941@gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
# as published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
#     https://www.gnu.org/licenses/
# ****************************************************************************
from copy import copy
from sage.dynamics.arithmetic_dynamics.generic_ds import DynamicalSystem
from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective
from sage.rings.integer_ring import ZZ
from sage.rings.quotient_ring import QuotientRing_generic
from sage.schemes.product_projective.morphism import ProductProjectiveSpaces_morphism_ring


class DynamicalSystem_product_projective(DynamicalSystem,
                                         ProductProjectiveSpaces_morphism_ring):
    r"""
    The class of dynamical systems on products of projective spaces.

    .. WARNING::

        You should not create objects of this class directly because
        no type or consistency checking is performed. The preferred
        method to construct such dynamical systems is to use
        :func:`~sage.dynamics.arithmetic_dynamics.generic_ds.DynamicalSystem_projective`
        function.

    INPUT:

    - ``polys`` -- list of `n_1 + \cdots + n_r` multi-homogeneous polynomials, all
      of which should have the same parent

    - ``domain`` -- a projective scheme embedded in
      `P^{n_1-1} \times \cdots \times P^{n_r-1}`

    EXAMPLES::

        sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
        sage: DynamicalSystem_projective([x^2, y^2, z^2, w^2, u^2], domain=T)
        Dynamical System of Product of projective spaces P^2 x P^1 over Rational Field
          Defn: Defined by sending (x : y : z , w : u) to (x^2 : y^2 : z^2 , w^2 : u^2).
    """

    def __init__(self, polys, domain):
        r"""
        The Python constructor.

        See :class:`DynamicalSystem` for details.

        EXAMPLES::

            sage: T.<x,y,w,u> = ProductProjectiveSpaces([1, 1], QQ)
            sage: DynamicalSystem_projective([x^2, y^2, w^2, u^2], domain=T)
            Dynamical System of Product of projective spaces P^1 x P^1 over Rational Field
              Defn: Defined by sending (x : y , w : u) to
                    (x^2 : y^2 , w^2 : u^2).
        """
        DynamicalSystem.__init__(self, polys, domain)

    def _call_with_args(self, P, check=True):
        r"""
        Make dynamical systems of products of projective spaces callable.

        INPUT:

        - ``P`` -- a point in the domain

        - ``check`` -- boolean (default: ``True``); whether or not to perform
          the input checks on the image point

        OUTPUT: the image point in the codomain

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: F = DynamicalSystem_projective([x^2*u, y^2*w, z^2*u, w^2, u^2], domain=T)
            sage: F(T([2, 1, 3, 0, 1]))
            (4/9 : 0 : 1 , 0 : 1)
        """
        if check:
            from sage.schemes.product_projective.point import ProductProjectiveSpaces_point_ring
            if not isinstance(P, ProductProjectiveSpaces_point_ring):
                try:
                    P = self.domain()(P)
                except (TypeError, NotImplementedError):
                    raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented" % (P, self.domain()))
            elif self.domain() != P.codomain():
                raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented" % (P, self.domain()))

        A = self.domain()
        Q = list(P)
        newP = [f(Q) for f in self.defining_polynomials()]
        return A.point(newP, check)

    def nth_iterate(self, P, n, normalize=False):
        r"""
        Return the ``n``-th iterate of ``P`` by this dynamical system.

        If ``normalize`` is ``True``, then the coordinates are
        automatically normalized.

        .. TODO:: Is there a more efficient way to do this?

        INPUT:

        - ``P`` -- a point in ``self.domain()``

        - ``n`` -- positive integer

        - ``normalize`` -- boolean (default: ``False``)

        OUTPUT: a point in ``self.codomain()``

        EXAMPLES::

            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1, 2], QQ)
            sage: f = DynamicalSystem_projective([a^3, b^3 + a*b^2, x^2, y^2 - z^2, z*y],
            ....:                                domain=Z)
            sage: P = Z([1, 1, 1, 1, 1])
            sage: f.nth_iterate(P, 3)
            (1/1872 : 1 , 1 : 1 : 0)

        ::

            sage: Z.<a,b,x,y> = ProductProjectiveSpaces([1, 1], ZZ)
            sage: f = DynamicalSystem_projective([a*b, b^2, x^3 - y^3, y^2*x], domain=Z)
            sage: P = Z([2, 6, 2, 4])
            sage: f.nth_iterate(P, 2, normalize=True)
            (1 : 3 , 407 : 112)
        """
        if P.codomain() != self.domain():
            raise TypeError("point is not defined over domain of function")
        n = ZZ(n)
        if n < 0:
            raise TypeError("must be a forward orbit")
        if n == 0:
            return self
        else:
            Q = self(P)
            if normalize:
                Q.normalize_coordinates()
            for i in range(2,n+1):
                Q = self(Q)
                if normalize:
                    Q.normalize_coordinates()
            return Q

    def orbit(self, P, N, **kwds):
        r"""
        Return the orbit of `P` by this dynamical system.

        Let `F` be this dynamical system. If `N` is an integer return `[P,F(P),\ldots,F^N(P)]`.

        If `N` is a list or tuple `N = [m, k]` return
        `[F^m(P),\ldots,F^k(P)]`.  Automatically normalize the points
        if ``normalize`` is ``True``. Perform the checks on point
        initialize if ``check`` is ``True``.

        INPUT:

        - ``P`` -- a point in ``self.domain()``

        - ``N`` -- nonnegative integer or list or tuple of two nonnegative integers

        kwds:

        - ``check`` -- boolean (default: ``True``)

        - ``normalize`` -- boolean (default: ``False``)

        OUTPUT: list of points in ``self.codomain()``

        EXAMPLES::

            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1, 2], QQ)
            sage: f = DynamicalSystem_projective([a^3, b^3 + a*b^2, x^2, y^2 - z^2, z*y],
            ....:                                domain=Z)
            sage: P = Z([1, 1, 1, 1, 1])
            sage: f.orbit(P, 3)
            [(1 : 1 , 1 : 1 : 1), (1/2 : 1 , 1 : 0 : 1),
             (1/12 : 1 , -1 : 1 : 0), (1/1872 : 1 , 1 : 1 : 0)]

        ::

            sage: Z.<a,b,x,y> = ProductProjectiveSpaces([1, 1], ZZ)
            sage: f = DynamicalSystem_projective([a*b, b^2, x^3 - y^3, y^2*x], domain=Z)
            sage: P = Z([2, 6, 2, 4])
            sage: f.orbit(P, 3, normalize=True)
            [(1 : 3 , 1 : 2), (1 : 3 , -7 : 4),
             (1 : 3 , 407 : 112), (1 : 3 , 66014215 : 5105408)]
        """
        if P.codomain() != self.domain():
            raise TypeError("point is not defined over domain of function")
        if not isinstance(N, (list,tuple)):
            N = [0, N]
        try:
            N[0] = ZZ(N[0])
            N[1] = ZZ(N[1])
        except TypeError:
            raise TypeError("orbit bounds must be integers")
        if N[0] < 0 or N[1] < 0:
            raise TypeError("orbit bounds must be nonnegative")
        if N[0] > N[1]:
            return []

        Q = copy(P)
        check = kwds.pop("check", True)
        normalize = kwds.pop("normalize", False)

        if normalize:
            Q.normalize_coordinates()
        for i in range(1, N[0]+1):
            Q = self(Q, check)
            if normalize:
                Q.normalize_coordinates()
        orb = [Q]
        for i in range(N[0]+1, N[1]+1):
            Q = self(Q, check)
            if normalize:
                Q.normalize_coordinates()
            orb.append(Q)
        return orb

    def nth_iterate_map(self, n):
        r"""
        Return the ``n``-th iterate of this dynamical system.

        ALGORITHM:

        Uses a form of successive squaring to reduce computations.


        .. TODO:: This could be improved.

        INPUT:

        - ``n`` -- positive integer

        OUTPUT: a dynamical system of products of projective spaces

        EXAMPLES::

            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1 , 2], QQ)
            sage: f = DynamicalSystem_projective([a^3, b^3, x^2, y^2, z^2], domain=Z)
            sage: f.nth_iterate_map(3)
            Dynamical System of Product of projective spaces P^1 x P^2 over
            Rational Field
              Defn: Defined by sending (a : b , x : y : z) to
                    (a^27 : b^27 , x^8 : y^8 : z^8).
        """
        E = self.domain()
        D = int(n)
        if D < 0:
            raise TypeError("iterate number must be a nonnegative integer")
        F = list(self._polys)
        Coord_ring = E.coordinate_ring()
        if isinstance(Coord_ring, QuotientRing_generic):
            PHI = [gen.lift() for gen in Coord_ring.gens()]
        else:
            PHI = list(Coord_ring.gens())

        while D:
            if D & 1:
                PHI = [poly(*F) for poly in PHI]
            if D > 1:  # avoid extra iterate
                F = [poly(*F) for poly in F]  # 'square'
            D >>= 1
        return DynamicalSystem_projective(PHI, domain=self.domain())


class DynamicalSystem_product_projective_field(DynamicalSystem_product_projective):
    pass


class DynamicalSystem_product_projective_finite_field(DynamicalSystem_product_projective_field):

    def cyclegraph(self):
        r"""
        Return the digraph of all orbits of this morphism mod `p`.

        OUTPUT: a digraph

        EXAMPLES::

            sage: P.<a,b,c,d> = ProductProjectiveSpaces(GF(3), [1,1])
            sage: f = DynamicalSystem_projective([a^2, b^2, c^2, d^2], domain=P)
            sage: f.cyclegraph()                                                        # needs sage.graphs
            Looped digraph on 16 vertices

        ::

            sage: P.<a,b,c,d> = ProductProjectiveSpaces(GF(5), [1,1])
            sage: f = DynamicalSystem_projective([a^2, b^2, c, d], domain=P)
            sage: f.cyclegraph()                                                        # needs sage.graphs
            Looped digraph on 36 vertices

        ::

            sage: P.<a,b,c,d,e> = ProductProjectiveSpaces(GF(2), [1,2])
            sage: f = DynamicalSystem_projective([a^2, b^2, c, d, e], domain=P)
            sage: f.cyclegraph()                                                        # needs sage.graphs
            Looped digraph on 21 vertices

        .. TODO:: Dynamical systems for subschemes of product projective spaces needs work.
                  Thus this is not implemented for subschemes.
        """
        V = []
        E = []
        from sage.schemes.product_projective.space import ProductProjectiveSpaces_ring
        if isinstance(self.domain(), ProductProjectiveSpaces_ring):
            for P in self.domain():
                V.append(str(P))
                Q = self(P)
                E.append([str(Q)])
        else:
            raise NotImplementedError("cyclegraph for product projective spaces not implemented for subschemes")
        from sage.graphs.digraph import DiGraph
        return DiGraph(dict(zip(V, E)), loops=True)
