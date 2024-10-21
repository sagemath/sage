# sage_setup: distribution = sagemath-categories
r"""
Points for products of projective spaces

This class builds on the projective space class and its point and morphism classes.

EXAMPLES:

We construct products projective spaces of various dimensions over the same ring.::

    sage: P1xP1.<x,y, u,v> = ProductProjectiveSpaces(QQ, [1, 1])
    sage: P1xP1([2, 1, 3, 1])
    (2 : 1 , 3 : 1)
"""
# ****************************************************************************
# Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#                    Ben Hutz <bn4941@gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
# as published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
# https://www.gnu.org/licenses/
# ****************************************************************************
from copy import copy

import sage.rings.abc

from sage.categories.integral_domains import IntegralDomains
from sage.categories.number_fields import NumberFields
from sage.rings.fraction_field import FractionField
from sage.rings.integer_ring import ZZ
from sage.schemes.generic.morphism import SchemeMorphism
from sage.schemes.generic.morphism import SchemeMorphism_point
from sage.structure.sequence import Sequence
from sage.structure.richcmp import richcmp


class ProductProjectiveSpaces_point_ring(SchemeMorphism_point):
    r"""
    The class of points on products of projective spaces.

    The components are projective space points.

    EXAMPLES::

        sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
        sage: T.point([1, 2, 3, 4, 5])
        (1/3 : 2/3 : 1 , 4/5 : 1)
    """
    def __init__(self, parent, polys, check=True):
        r"""
        The Python constructor.

        INPUT:

        - ``parent`` -- Hom-set

        - ``polys`` -- anything that defines a point in the class

        - ``check`` -- boolean (default: ``True``); whether or not to perform
          input checks

        EXAMPLES::

            sage: P1.<x0,x1,x2> = ProjectiveSpace(QQ, 2)
            sage: P2 = ProjectiveSpace(QQ, 3, 'y')
            sage: T = ProductProjectiveSpaces([P1, P2])
            sage: Q1 = P1(1, 1, 1)
            sage: Q2 = P2(1, 2, 3, 4)
            sage: T([Q1, Q2])
            (1 : 1 : 1 , 1/4 : 1/2 : 3/4 : 1)

        ::

            sage: T = ProductProjectiveSpaces([2, 2, 2], GF(5), 'x')
            sage: T.point([1, 2, 3, 4, 5, 6, 7, 8, 9])
            (2 : 4 : 1 , 4 : 0 : 1 , 3 : 2 : 1)

        ::

            sage: T.<x,y,z,w> = ProductProjectiveSpaces([1, 1], GF(5))
            sage: X = T.subscheme([x - y, z - 2*w])
            sage: X([1, 1, 2, 1])
            (1 : 1 , 2 : 1)
        """
        polys = copy(polys)
        SchemeMorphism.__init__(self, parent)
        if all(isinstance(P, SchemeMorphism_point) for P in polys):
            if check:
                Q = []
                self._points = []
                for i in range(len(polys)):
                    if polys[i].codomain() != parent.codomain().ambient_space()[i]:
                        raise ValueError("points must be in correct projective spaces")
                    Q += list(polys[i])
                    self._points.append(polys[i])
                parent.codomain()._check_satisfies_equations(Q)
            self._points = polys
        else:
            R = parent.codomain().ambient_space().base_ring()
            polys = Sequence(polys, R)
            N = parent.codomain().ambient_space().dimension_relative_components()
            if check:
                parent.codomain()._check_satisfies_equations(polys)
            splitpolys = self.codomain().ambient_space()._factors(polys)
            self._points = [parent.codomain().ambient_space()[i].point(splitpolys[i], check) for i in range(len(N))]

    def __getitem__(self, i):
        r"""
        Return the ``i``-th coordinate point.

        INPUT:

        - ``i`` -- integer

        OUTPUT: the projective space point that is the ``i``-th coordinate

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([2, 2, 2], GF(5), 'x')
            sage: P = T([1, 0, 1, 1, 0, 0, 0, 0, 1])
            sage: P[1]
            (1 : 0 : 0)
            sage: P[1].codomain()
            Projective Space of dimension 2 over Finite Field of size 5
            sage: P[1][0]
            1
        """
        return self._points[i]

    def _repr_(self):
        r"""
        Return a string representation of this point.

        OUTPUT: string

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([2, 2], ZZ, 'x')
            sage: P = T([1, 2, 3, 4, 5, 6])
            sage: P._repr_()
            '(1 : 2 : 3 , 4 : 5 : 6)'
        """
        return '(%s)' % (" , ".join((" : ".join(repr(f) for f in Q))
                                    for Q in self._points))

    def _richcmp_(self, right, op):
        r"""
        Compare two points in products of projective spaces.

        INPUT:

        - ``other`` -- another point

        OUTPUT: boolean

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([2, 2], ZZ, 'x')
            sage: P = T([1, 2, 3, 4, 5, 6])
            sage: Q = T([2, 4, 6, 4, 5, 6])
            sage: P == Q
            True

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1, 1, 1], ZZ, 'x')
            sage: P = T([1, 2, 3, 4, 5, 6])
            sage: Q = T([2, 4, 6, 4, 1, 0])
            sage: P != Q
            True

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1, 1, 1], GF(5), 'x')
            sage: P = T([3, 2, 3, 4, 1, 0])
            sage: Q = T([1, 2, 3, 4, 3, 1])
            sage: P > Q
            True

        ::

            sage: T = ProductProjectiveSpaces([1, 1, 1], GF(5), 'x')
            sage: P = T([1, 2, 3, 4, 1, 0])
            sage: Q = T([1, 2, 3, 4, 3, 0])
            sage: P == Q
            True

        ::

            sage: T = ProductProjectiveSpaces([1, 1, 1], GF(5), 'x')
            sage: P = T([1, 2, 3, 4, 1, 0])
            sage: Q = T([1, 2, 3, 4, 3, 1])
            sage: P < Q
            True
        """
        #needed for Digraph
        if not isinstance(right, (ProductProjectiveSpaces_point_ring)):
            return NotImplemented
        else:
            return richcmp(self._points, right._points, op)

    def __copy__(self):
        r"""
        Return a copy of this point.

        OUTPUT: a point in the same space as third point

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1, 1], QQ, 'x')
            sage: P = T([2, 1, 0, 1])
            sage: Q = P.__copy__()
            sage: P is Q
            False
            sage: P == Q
            True
        """
        P = [copy(self[i]) for i in range(self.codomain().ambient_space().num_components())]
        return (self.codomain().point(P, False))

    def __iter__(self):
        r"""
        Iterate over the coordinates of the point.

        OUTPUT: an iterator

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1, 1], QQ, 'x')
            sage: P = T([2, 1, 0, 1])
            sage: iter = P.__iter__()
            sage: next(iter)
            2
            sage: next(iter)
            1
            sage: list(P)
            [2, 1, 0, 1]
        """
        return (x for P in self._points for x in P._coords)

    def __len__(self):
        """
        Return the total number of coordinates in ``self``.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1, 1], QQ, 'x')
            sage: P = T([2, 1, 0, 1])
            sage: len(P)
            4
        """
        image = self.codomain().ambient_space()
        return image.dimension() + image.num_components()

    def __hash__(self):
        """
        Compute the hash value of this point.

        OUTPUT: integer

        EXAMPLES::

            sage: PP = ProductProjectiveSpaces(Zmod(6), [1, 1])
            sage: H = hash(PP([5, 1, 3, 4]))

        ::

            sage: PP = ProductProjectiveSpaces(ZZ, [1, 2])
            sage: hash(PP([1, 1, 2, 2, 2])) == hash(PP([1, 1, 1, 1, 1]))
            True

        ::

            sage: PP = ProductProjectiveSpaces(QQ, [1, 1])
            sage: hash(PP([1/7, 1, 2, 1])) == hash((1/7, 1, 2, 1))
            True

        ::

            sage: PP = ProductProjectiveSpaces(GF(7), [1, 1, 1])
            sage: hash(PP([4, 1, 5, 4, 6, 1])) == hash((4, 1, 5, 4, 6, 1))
            False
            sage: hash(PP([4, 1, 5, 4, 6, 1])) == hash((4, 1, 3, 1, 6, 1))
            True
        """
        R = self.codomain().base_ring()
        # if there is a fraction field normalize the point so that
        # equal points have equal hash values
        if R in IntegralDomains():
            P = self.change_ring(FractionField(R))
            P.normalize_coordinates()
            return hash(tuple(P))
        # if there is no good way to normalize return
        # a constant value
        return hash(self.codomain())

    def normalize_coordinates(self):
        r"""
        Remove common factors (componentwise) from the coordinates of this point (including `-1`).

        OUTPUT: none

        EXAMPLES::

            sage: T.<x,y,z,u,v,w> = ProductProjectiveSpaces([2, 2], ZZ)
            sage: P = T.point([5, 10, 15, 4, 2, 6]);
            sage: P.normalize_coordinates()
            sage: P
            (1 : 2 : 3 , 2 : 1 : 3)
        """
        for i in range(self.codomain().ambient_space().num_components()):
            self[i].normalize_coordinates()

    def dehomogenize(self, L):
        r"""
        Dehomogenize `k`-th point at `L[k]`-th coordinate.

        This function computes the appropriate affine patch using ``L``
        and then returns the dehomogenized point on of this affine space.

        INPUT:

        - ``L`` -- list of nonnegative integers

        OUTPUT: :class:`SchemeMorphism_point_affine`

        EXAMPLES::

            sage: PP = ProductProjectiveSpaces([2, 2, 2], QQ, 'x')
            sage: A = PP([2, 4, 6, 23, 46, 23, 9, 3, 1])
            sage: A.dehomogenize([0, 1, 2])
            (2, 3, 1/2, 1/2, 9, 3)

        ::

            sage: # needs sage.rings.real_mpfr sage.symbolic
            sage: PP.<a,b,x,y,z> = ProductProjectiveSpaces([1, 2], CC)
            sage: X = PP.subscheme([a^2 + b^2])
            sage: P = X([2, 2*i, -3, 6*i, 3 - 6*i])
            sage: P.dehomogenize([1,0])
            (-1.00000000000000*I, -2.00000000000000*I, -1.00000000000000 + 2.00000000000000*I)

        ::

            sage: PP = ProductProjectiveSpaces([1, 1], ZZ)
            sage: A = PP([0,1,2,4])
            sage: A.dehomogenize([0,0])
            Traceback (most recent call last):
            ...
            ValueError: can...t dehomogenize at 0 coordinate
        """
        PP = self.codomain()
        A = PP.affine_patch(L)
        pt = []
        for i in range(PP.ambient_space().num_components()):
            pt.extend(self[i].dehomogenize(L[i]))
        return A(pt)

    def scale_by(self, t):
        r"""
        Scale the coordinates of the point by ``t``, done componentwise.

        A :exc:`TypeError` occurs if the point is not in the base ring
        of the codomain after scaling.

        INPUT:

        - ``t`` -- a ring element

        EXAMPLES::

            sage: T.<x, y, z, u, v, w> = ProductProjectiveSpaces([1, 1, 1], ZZ)
            sage: P = T.point([5, 10, 15, 4, 2, 6]);
            sage: P.scale_by([2, 1, 1])
            sage: P
            (10 : 20 , 15 : 4 , 2 : 6)
        """
        if not isinstance(t, (tuple, list)):
            raise TypeError("%s must be a list or tuple" % t)
        if len(t) != self.codomain().ambient_space().num_components():
            raise TypeError("%s must have same number of components as %r" % (t, self))
        for i in range(self.codomain().ambient_space().num_components()):
            self[i].scale_by(t[i])

    def change_ring(self, R, **kwds):
        r"""
        Return a new :class:`ProductProjectiveSpaces_point` which is this point coerced to ``R``.

        If the keyword ``check`` is ``True``, then the initialization checks are performed.
        The user may specify the embedding into ``R`` with a keyword.

        INPUT:

        - ``R`` -- ring

        kwds:

        - ``check`` -- boolean

        - ``embedding`` -- field embedding from the base ring of this point to ``R``

        OUTPUT: :class:`ProductProjectiveSpaces_point`

        EXAMPLES::

            sage: T.<x,y,z,u,v,w> = ProductProjectiveSpaces([1, 1, 1], ZZ)
            sage: P = T.point([5, 3, 15, 4, 2, 6])
            sage: P.change_ring(GF(3))
            (1 : 0 , 0 : 1 , 1 : 0)
        """
        check = kwds.get('check', True)
        S = self.codomain().change_ring(R)
        Q = [P.change_ring(R, **kwds) for P in self._points]
        return S.point(Q, check)

    def global_height(self, prec=None):
        r"""
        Return the absolute logarithmic height of the point.

        This function computes the maximum of global height of each
        component point in the product. Global height of component
        point is computed using function for projective point.

        INPUT:

        - ``prec`` -- desired floating point precision (default:
          default RealField precision)

        OUTPUT: a real number

        EXAMPLES::

            sage: PP = ProductProjectiveSpaces(QQ, [2, 2], 'x')
            sage: Q = PP([1, 7, 5, 18, 2, 3])
            sage: Q.global_height()                                                     # needs sage.symbolic
            2.89037175789616

        ::

            sage: PP = ProductProjectiveSpaces(ZZ, [1, 1], 'x')
            sage: A = PP([-30, 2, 1, 6])
            sage: A.global_height()                                                     # needs sage.symbolic
            2.70805020110221

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: k.<w> = NumberField(x^2 + 5)
            sage: PP = ProductProjectiveSpaces(k, [1, 2], 'y')
            sage: Q = PP([3, 5*w + 1, 1, 7*w, 10])
            sage: Q.global_height()
            2.75062910527236

        ::

            sage: PP = ProductProjectiveSpaces(QQbar, [1, 1], 'x')                      # needs sage.rings.number_field
            sage: Q = PP([1, QQbar(sqrt(2)), QQbar(5^(1/3)), QQbar(3^(1/3))])           # needs sage.rings.number_field sage.symbolic
            sage: Q.global_height()                                                     # needs sage.rings.number_field sage.symbolic
            0.536479304144700
        """
        K = self.codomain().base_ring()
        if K not in NumberFields() and K != ZZ and not isinstance(K, (sage.rings.abc.Order, sage.rings.abc.AlgebraicField)):
            raise TypeError("must be over a number field or a number field order or QQbar")

        n = self.codomain().ambient_space().num_components()
        return max(self[i].global_height(prec=prec) for i in range(n))

    def local_height(self, v, prec=None):
        r"""
        Return the maximum of the local height of the coordinates of this point.

        This function computes the maximum of local height of each
        component point in the product. Local height of component
        point is computed using function for projective point.

        INPUT:

        - ``v`` -- a prime or prime ideal of the base ring

        - ``prec`` -- desired floating point precision (default:
          default RealField precision)

        OUTPUT: a real number

        EXAMPLES::

            sage: PP = ProductProjectiveSpaces(QQ, [1, 1], 'x')
            sage: A = PP([11, 5, 10, 2])
            sage: A.local_height(5)                                                     # needs sage.rings.real_mpfr
            1.60943791243410

        ::

            sage: P = ProductProjectiveSpaces(QQ, [1, 2], 'x')
            sage: Q = P([1, 4, 1/2, 2, 32])
            sage: Q.local_height(2)                                                     # needs sage.rings.real_mpfr
            4.15888308335967
        """
        K = FractionField(self.domain().base_ring())
        if K not in NumberFields():
            raise TypeError("must be over a number field or a number field order")

        n = self.codomain().ambient_space().num_components()
        return max(self[i].local_height(v, prec=prec) for i in range(n))


class ProductProjectiveSpaces_point_field(ProductProjectiveSpaces_point_ring):

    def intersection_multiplicity(self, X):
        r"""
        Return the intersection multiplicity of the codomain of this point and subscheme ``X`` at this point.

        This uses the subscheme implementation of intersection_multiplicity. This point must be a point
        on a subscheme of a product of projective spaces.

        INPUT:

        - ``X`` -- a subscheme in the same ambient space as the codomain of this point

        OUTPUT: integer

        EXAMPLES::

            sage: PP.<x,y,z,u,v> = ProductProjectiveSpaces(QQ, [2, 1])
            sage: X = PP.subscheme([y^2*z^3*u - x^5*v])
            sage: Y = PP.subscheme([u^3 - v^3, x - y])
            sage: Q = X([0,0,1,1,1])
            sage: Q.intersection_multiplicity(Y)                                        # needs sage.libs.singular
            2
        """
        from sage.schemes.product_projective.space import ProductProjectiveSpaces_ring
        if isinstance(self.codomain(), ProductProjectiveSpaces_ring):
            raise TypeError("this point must be a point on a subscheme of a product of projective spaces")
        return self.codomain().intersection_multiplicity(X, self)

    def multiplicity(self):
        r"""
        Return the multiplicity of this point on its codomain.

        This uses the subscheme implementation of multiplicity. This point must be a point
        on a subscheme of a product of projective spaces.

        OUTPUT: integer

        EXAMPLES::

            sage: PP.<x,y,z,w,u,v,t> = ProductProjectiveSpaces(QQ, [3, 2])
            sage: X = PP.subscheme([x^8*t - y^8*t + z^5*w^3*v])
            sage: Q1 = X([1,1,0,0,-1,-1,1])
            sage: Q1.multiplicity()                                                     # needs sage.libs.singular
            1
            sage: Q2 = X([0,0,0,1,0,1,1])
            sage: Q2.multiplicity()                                                     # needs sage.libs.singular
            5
            sage: Q3 = X([0,0,0,1,1,0,0])
            sage: Q3.multiplicity()                                                     # needs sage.libs.singular
            6
        """
        from sage.schemes.product_projective.space import ProductProjectiveSpaces_ring
        if isinstance(self.codomain(), ProductProjectiveSpaces_ring):
            raise TypeError("this point must be a point on a subscheme of a product of projective spaces")
        return self.codomain().multiplicity(self)


class ProductProjectiveSpaces_point_finite_field(ProductProjectiveSpaces_point_field):
    pass
