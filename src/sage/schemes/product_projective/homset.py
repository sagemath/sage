r"""
Set of homomorphisms
"""

#*****************************************************************************
# Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#                    Ben Hutz <bn4941@gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
# as published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
# http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.fields import Fields
from sage.categories.number_fields import NumberFields
from sage.misc.mrange import xmrange
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField
from sage.schemes.generic.homset import SchemeHomset_points

class SchemeHomset_points_product_projective_spaces_ring(SchemeHomset_points):
    r"""
    Set of rational points of a product of projective spaces.

    INPUT:

    See :class:`~sage.schemes.generic.homset.SchemeHomset_generic`.

    EXAMPLES::

        sage: from sage.schemes.product_projective.homset import SchemeHomset_points_product_projective_spaces_ring
        sage: SchemeHomset_points_product_projective_spaces_ring(Spec(QQ), \
        ProductProjectiveSpaces([1, 1], QQ, 'z'))
        Set of rational points of Product of projective spaces P^1 x P^1 over Rational Field
        """

    def _element_constructor_(self, v, **kwds):
        r"""
        The element constructor.

        INPUT:

        - ``v`` -- anything that determines a scheme morphism in the Hom-set.

        OUTPUT:

        The scheme morphism determined by ``v``.

        EXAMPLES::

            sage: P = ProductProjectiveSpaces([1, 1], ZZ, 'z')
            sage: Q = P([4, 6, 6, 2]); Q
            (4 : 6 , 6 : 2)
            sage: type(Q)
            <class 'sage.schemes.product_projective.point.ProductProjectiveSpaces_point_ring'>
            sage: P(QQ)._element_constructor_([4, 2, 2, 0])
            (4 : 2 , 2 : 0)
        """
        return self.codomain()._point(self, v, **kwds)

class SchemeHomset_points_product_projective_spaces_field(SchemeHomset_points_product_projective_spaces_ring):
    def points(self, B=0, prec=53):
        r"""
        Return some or all rational points of a projective scheme.

        Over a finite field, all points are returned. Over an infinite field, all points satisfying the bound
        are returned. For a zero-dimensional subscheme, all points are returned regardless of whether the base
        ring is a field or not.

        INPUT:

        - `B` -- integer (optional, default=0). The bound for the
          coordinates.
        - ``prec`` - the precision to use to compute the elements of bounded height for number fields.

        OUTPUT:

        - a list of rational points of a projective scheme.

        .. WARNING::

           In the current implementation, the output of the [Doyle-Krumm]_ algorithm
           cannot be guaranteed to be correct due to the necessity of floating point
           computations. In some cases, the default 53-bit precision is
           considerably lower than would be required for the algorithm to
           generate correct output.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProductProjectiveSpaces([1, 1], QQ)
            sage: X = P.subscheme([x - y, z^2 - 2*w^2])
            sage: X(P.base_ring()).points()
            []

        ::

            sage: u = QQ['u'].0
            sage: P.<x,y,z,w> = ProductProjectiveSpaces([1,1], NumberField(u^2 - 2, 'v'))
            sage: X = P.subscheme([x^2 - y^2, z^2 - 2*w^2])
            sage: X(P.base_ring()).points()
            [(-1 : 1 , -v : 1), (1 : 1 , v : 1), (1 : 1 , -v : 1), (-1 : 1 , v : 1)]

        ::

            sage: u = QQ['u'].0
            sage: K = NumberField(u^2 + 1, 'v')
            sage: P.<x,y,z,w> = ProductProjectiveSpaces([1, 1], K)
            sage: P(K).points(1)
            [(0 : 1 , 0 : 1), (0 : 1 , v : 1), (0 : 1 , -1 : 1), (0 : 1 , -v : 1), (0 : 1 , 1 : 1),
            (0 : 1 , 1 : 0), (v : 1 , 0 : 1), (v : 1 , v : 1), (v : 1 , -1 : 1), (v : 1 , -v : 1),
            (v : 1 , 1 : 1), (v : 1 , 1 : 0), (-1 : 1 , 0 : 1), (-1 : 1 , v : 1), (-1 : 1 , -1 : 1),
            (-1 : 1 , -v : 1), (-1 : 1 , 1 : 1), (-1 : 1 , 1 : 0), (-v : 1 , 0 : 1), (-v : 1 , v : 1),
            (-v : 1 , -1 : 1), (-v : 1 , -v : 1), (-v : 1 , 1 : 1), (-v : 1 , 1 : 0), (1 : 1 , 0 : 1),
            (1 : 1 , v : 1), (1 : 1 , -1 : 1), (1 : 1 , -v : 1), (1 : 1 , 1 : 1), (1 : 1 , 1 : 0),
            (1 : 0 , 0 : 1), (1 : 0 , v : 1), (1 : 0 , -1 : 1), (1 : 0 , -v : 1), (1 : 0 , 1 : 1),
            (1 : 0 , 1 : 0)]

        ::

            sage: P.<x,y,z,u,v> = ProductProjectiveSpaces([2, 1], GF(3))
            sage: P(P.base_ring()).points()
            [(0 : 0 : 1 , 0 : 1), (1 : 0 : 1 , 0 : 1), (2 : 0 : 1 , 0 : 1), (0 : 1 : 1 , 0 : 1), (1 : 1 : 1 , 0 : 1),
            (2 : 1 : 1 , 0 : 1), (0 : 2 : 1 , 0 : 1), (1 : 2 : 1 , 0 : 1), (2 : 2 : 1 , 0 : 1), (0 : 1 : 0 , 0 : 1),
            (1 : 1 : 0 , 0 : 1), (2 : 1 : 0 , 0 : 1), (1 : 0 : 0 , 0 : 1), (0 : 0 : 1 , 1 : 1), (1 : 0 : 1 , 1 : 1),
            (2 : 0 : 1 , 1 : 1), (0 : 1 : 1 , 1 : 1), (1 : 1 : 1 , 1 : 1), (2 : 1 : 1 , 1 : 1), (0 : 2 : 1 , 1 : 1),
            (1 : 2 : 1 , 1 : 1), (2 : 2 : 1 , 1 : 1), (0 : 1 : 0 , 1 : 1), (1 : 1 : 0 , 1 : 1), (2 : 1 : 0 , 1 : 1),
            (1 : 0 : 0 , 1 : 1), (0 : 0 : 1 , 2 : 1), (1 : 0 : 1 , 2 : 1), (2 : 0 : 1 , 2 : 1), (0 : 1 : 1 , 2 : 1),
            (1 : 1 : 1 , 2 : 1), (2 : 1 : 1 , 2 : 1), (0 : 2 : 1 , 2 : 1), (1 : 2 : 1 , 2 : 1), (2 : 2 : 1 , 2 : 1),
            (0 : 1 : 0 , 2 : 1), (1 : 1 : 0 , 2 : 1), (2 : 1 : 0 , 2 : 1), (1 : 0 : 0 , 2 : 1), (0 : 0 : 1 , 1 : 0),
            (1 : 0 : 1 , 1 : 0), (2 : 0 : 1 , 1 : 0), (0 : 1 : 1 , 1 : 0), (1 : 1 : 1 , 1 : 0), (2 : 1 : 1 , 1 : 0),
            (0 : 2 : 1 , 1 : 0), (1 : 2 : 1 , 1 : 0), (2 : 2 : 1 , 1 : 0), (0 : 1 : 0 , 1 : 0), (1 : 1 : 0 , 1 : 0),
            (2 : 1 : 0 , 1 : 0), (1 : 0 : 0 , 1 : 0)]
        """
        X = self.codomain()

        from sage.schemes.product_projective.space import is_ProductProjectiveSpaces
        if not is_ProductProjectiveSpaces(X) and X.base_ring() in Fields():
            # no points
            if X.dimension() == -1:
                return []
            # if X is zero-dimensional
            if X.dimension() == 0:
                points = set()
                # find points from all possible affine patches
                for I in xmrange([n + 1 for n in X.ambient_space().dimension_relative_components()]):
                    [Y,phi] = X.affine_patch(I, True)
                    aff_points = Y.rational_points()
                    for PP in aff_points:
                        points.add(phi(PP))
                return list(points)
        R = self.value_ring()
        points = []
        if R in NumberFields():
            if not B > 0:
                raise TypeError("a positive bound B (= %s) must be specified"%B)
            for P in X.ambient_space().points_of_bounded_height(B, prec):
                try:
                    points.append(X(P))
                except TypeError:
                    pass
            return points
        if is_FiniteField(R):
            for P in X.ambient_space().rational_points():
                try:
                    points.append(X(P))
                except TypeError:
                    pass
            return points
        else:
            raise TypeError("unable to enumerate points over %s"%R)
