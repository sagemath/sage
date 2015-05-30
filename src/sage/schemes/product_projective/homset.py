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

from sage.schemes.generic.homset import SchemeHomset_points

class SchemeHomset_points_product_projective_spaces_ring(SchemeHomset_points):
    r"""
    Set of rational points of a product of projective spaces.

    INPUT:

    See :class:`~sage.schemes.generic.homset.SchemeHomset_generic`.

    EXAMPLES::

        sage: from sage.schemes.product_projective.homset import SchemeHomset_points_product_projective_spaces_ring
        sage: SchemeHomset_points_product_projective_spaces_ring(Spec(QQ), ProductProjectiveSpaces([1,1],QQ,'z'))
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

            sage: P = ProductProjectiveSpaces([1,1],QQ,'z')
            sage: Q = P([1,2,2,3]); Q
            (1/2 : 1 , 2/3 : 1)
            sage: type(Q)
            <class 'sage.schemes.product_projective.point.ProductProjectiveSpaces_point_ring'>
            sage: P(QQ)._element_constructor_([1,2,2,0])
            (1/2 : 1 , 1 : 0)
        """
        return self.codomain()._point(self, v, **kwds)
