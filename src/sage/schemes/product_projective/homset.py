r"""
Set of homomorphisms

AUTHORS:

- Volker Braun and Ben Hutz (2014): initial version

- Raghukul Raman (2018): code cleanup and added support for rational field
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
from sage.misc.misc_c import prod
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField
from sage.rings.rational_field import is_RationalField
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme
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
    def points(self, **kwds):
        r"""
        Return some or all rational points of a projective scheme.

        Over a finite field, all points are returned. Over an infinite field, all points satisfying the bound
        are returned. For a zero-dimensional subscheme, all points are returned regardless of whether the base
        ring is a field or not.

        For number fields, this uses the
        Doyle-Krumm algorithm 4 (algorithm 5 for imaginary quadratic) for
        computing algebraic numbers up to a given height [DK2013]_ or
        uses the chinese remainder theorem and points modulo primes
        for larger bounds.

        The algorithm requires floating point arithmetic, so the user is
        allowed to specify the precision for such calculations.
        Additionally, due to floating point issues, points
        slightly larger than the bound may be returned. This can be controlled
        by lowering the tolerance.


        INPUT:

        - ``bound`` - a real number

        - ``tolerance`` - a rational number in (0,1] used in doyle-krumm algorithm-4

        - ``precision`` - the precision to use for computing the elements of bounded height of number fields.

        - ``algorithm`` - either 'sieve' or 'enumerate' algorithms can be used over `QQ`. If
          not specified, enumerate is used only for small height bounds.

        OUTPUT:

        - a list of rational points of a projective scheme

        EXAMPLES::

            sage: P.<x,y,z,w> = ProductProjectiveSpaces([1, 1], QQ)
            sage: X = P.subscheme([x - y, z^2 - 2*w^2])
            sage: X(P.base_ring()).points()
            []

        ::

            sage: u = QQ['u'].0
            sage: P.<x,y,z,w> = ProductProjectiveSpaces([1,1], NumberField(u^2 - 2, 'v'))
            sage: X = P.subscheme([x^2 - y^2, z^2 - 2*w^2])
            sage: sorted(X(P.base_ring()).points())
            [(-1 : 1 , -v : 1), (-1 : 1 , v : 1), (1 : 1 , -v : 1), (1 : 1 , v : 1)]

        ::

            sage: u = QQ['u'].0
            sage: K = NumberField(u^2 + 1, 'v')
            sage: P.<x,y,z,w> = ProductProjectiveSpaces([1, 1], K)
            sage: P(K).points(bound=1)
            [(-1 : 1 , -1 : 1), (-1 : 1 , -v : 1), (-1 : 1 , 0 : 1), (-1 : 1 , v : 1),
            (-1 : 1 , 1 : 0), (-1 : 1 , 1 : 1), (-v : 1 , -1 : 1), (-v : 1 , -v : 1),
            (-v : 1 , 0 : 1), (-v : 1 , v : 1), (-v : 1 , 1 : 0), (-v : 1 , 1 : 1),
            (0 : 1 , -1 : 1), (0 : 1 , -v : 1), (0 : 1 , 0 : 1), (0 : 1 , v : 1),
            (0 : 1 , 1 : 0), (0 : 1 , 1 : 1), (v : 1 , -1 : 1), (v : 1 , -v : 1),
            (v : 1 , 0 : 1), (v : 1 , v : 1), (v : 1 , 1 : 0), (v : 1 , 1 : 1),
            (1 : 0 , -1 : 1), (1 : 0 , -v : 1), (1 : 0 , 0 : 1), (1 : 0 , v : 1),
            (1 : 0 , 1 : 0), (1 : 0 , 1 : 1), (1 : 1 , -1 : 1), (1 : 1 , -v : 1),
            (1 : 1 , 0 : 1), (1 : 1 , v : 1), (1 : 1 , 1 : 0), (1 : 1 , 1 : 1)]

        ::

            sage: P.<x,y,z,u,v> = ProductProjectiveSpaces([2, 1], GF(3))
            sage: P(P.base_ring()).points()
            [(0 : 0 : 1 , 0 : 1), (0 : 0 : 1 , 1 : 0), (0 : 0 : 1 , 1 : 1), (0 : 0 : 1 , 2 : 1),
            (0 : 1 : 0 , 0 : 1), (0 : 1 : 0 , 1 : 0), (0 : 1 : 0 , 1 : 1), (0 : 1 : 0 , 2 : 1),
            (0 : 1 : 1 , 0 : 1), (0 : 1 : 1 , 1 : 0), (0 : 1 : 1 , 1 : 1), (0 : 1 : 1 , 2 : 1),
            (0 : 2 : 1 , 0 : 1), (0 : 2 : 1 , 1 : 0), (0 : 2 : 1 , 1 : 1), (0 : 2 : 1 , 2 : 1),
            (1 : 0 : 0 , 0 : 1), (1 : 0 : 0 , 1 : 0), (1 : 0 : 0 , 1 : 1), (1 : 0 : 0 , 2 : 1),
            (1 : 0 : 1 , 0 : 1), (1 : 0 : 1 , 1 : 0), (1 : 0 : 1 , 1 : 1), (1 : 0 : 1 , 2 : 1),
            (1 : 1 : 0 , 0 : 1), (1 : 1 : 0 , 1 : 0), (1 : 1 : 0 , 1 : 1), (1 : 1 : 0 , 2 : 1),
            (1 : 1 : 1 , 0 : 1), (1 : 1 : 1 , 1 : 0), (1 : 1 : 1 , 1 : 1), (1 : 1 : 1 , 2 : 1), 
            (1 : 2 : 1 , 0 : 1), (1 : 2 : 1 , 1 : 0), (1 : 2 : 1 , 1 : 1), (1 : 2 : 1 , 2 : 1),
            (2 : 0 : 1 , 0 : 1), (2 : 0 : 1 , 1 : 0), (2 : 0 : 1 , 1 : 1), (2 : 0 : 1 , 2 : 1),
            (2 : 1 : 0 , 0 : 1), (2 : 1 : 0 , 1 : 0), (2 : 1 : 0 , 1 : 1), (2 : 1 : 0 , 2 : 1),
            (2 : 1 : 1 , 0 : 1), (2 : 1 : 1 , 1 : 0), (2 : 1 : 1 , 1 : 1), (2 : 1 : 1 , 2 : 1),
            (2 : 2 : 1 , 0 : 1), (2 : 2 : 1 , 1 : 0), (2 : 2 : 1 , 1 : 1), (2 : 2 : 1 , 2 : 1)]

        ::

            sage: PP.<x,y,z,u,v> = ProductProjectiveSpaces([2,1], QQ)
            sage: X = PP.subscheme([x + y, u*u-v*u])
            sage: X.rational_points(bound=2)
            [(-2 : 2 : 1 , 0 : 1),
             (-2 : 2 : 1 , 1 : 1),
             (-1 : 1 : 0 , 0 : 1),
             (-1 : 1 : 0 , 1 : 1),
             (-1 : 1 : 1 , 0 : 1),
             (-1 : 1 : 1 , 1 : 1),
             (-1/2 : 1/2 : 1 , 0 : 1),
             (-1/2 : 1/2 : 1 , 1 : 1),
             (0 : 0 : 1 , 0 : 1),
             (0 : 0 : 1 , 1 : 1),
             (1/2 : -1/2 : 1 , 0 : 1),
             (1/2 : -1/2 : 1 , 1 : 1),
             (1 : -1 : 1 , 0 : 1),
             (1 : -1 : 1 , 1 : 1),
             (2 : -2 : 1 , 0 : 1),
             (2 : -2 : 1 , 1 : 1)]

        better to enumerate with low codimension::

            sage: PP.<x,y,z,u,v,a,b,c> = ProductProjectiveSpaces([2,1,2], QQ)
            sage: X = PP.subscheme([x*u^2*a, b*z*u*v,z*v^2*c ])
            sage: len(X.rational_points(bound=1, algorithm='enumerate'))
            232
         """
        B = kwds.pop('bound', 0)
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
        if is_RationalField(R):
            if not B > 0:
                raise TypeError("a positive bound B (= %s) must be specified"%B)
            alg = kwds.pop('algorithm', None)
            if alg is None:
                # sieve should only be called for subschemes and if the bound is not very small
                N = prod([k+1 for k in X.ambient_space().dimension_relative_components()])
                if isinstance(X, AlgebraicScheme_subscheme) and B**N > 5000:
                    from sage.schemes.product_projective.rational_point import sieve
                    return sieve(X, B)
                else:
                    from sage.schemes.product_projective.rational_point import enum_product_projective_rational_field
                    return enum_product_projective_rational_field(self, B)
            elif alg == 'sieve':
                from sage.schemes.product_projective.rational_point import sieve
                return sieve(X, B)
            elif alg == 'enumerate':
                from sage.schemes.product_projective.rational_point import enum_product_projective_rational_field
                return enum_product_projective_rational_field(self, B)
            else:
                raise ValueError("algorithm must be 'sieve' or 'enumerate'")
        elif R in NumberFields():
            if not B > 0:
                raise TypeError("a positive bound B (= %s) must be specified"%B)
            from sage.schemes.product_projective.rational_point import enum_product_projective_number_field
            return enum_product_projective_number_field(self, bound=B)
        elif is_FiniteField(R):
            from sage.schemes.product_projective.rational_point import enum_product_projective_finite_field
            return enum_product_projective_finite_field(self)
        else:
            raise TypeError("unable to enumerate points over %s" % R)
