"""
Hom-sets of weighted projective schemes

AUTHORS:

- Gareth Ma (2024): initial version, based on unweighted version.
"""

# *****************************************************************************
#        Copyright (C) 2024 Gareth Ma <grhkm21@gmail.com>
#        Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#        Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from sage.schemes.generic.homset import SchemeHomset_points


class SchemeHomset_points_weighted_projective_ring(SchemeHomset_points):
    """
    Set of rational points of a weighted projective variety over a ring.

    INPUT:

    See :class:`SchemeHomset_generic`.

    EXAMPLES::

        sage: WP = WeightedProjectiveSpace(ZZ, [1, 3, 1]); WP
        Weighted Projective Space of dimension 2 with weights (1, 3, 1) over Integer Ring
        sage: WP.point_homset()
        Set of rational points of Weighted Projective Space of dimension 2 with weights (1, 3, 1) over Integer Ring
        sage: type(WP.point_homset())
        <class 'sage.schemes.weighted_projective.weighted_projective_homset.SchemeHomset_points_weighted_projective_ring_with_category'>
        sage: WP(2, 24, 2)
        (2 : 24 : 2)
        sage: WP(2, 24, 2) == WP(1, 3, 1)
        True
    """

    def points(self, **__):
        """
        Return some or all rational points of this weighted projective scheme.

        For dimension 0 subschemes points are determined through a groebner
        basis calculation. For schemes or subschemes with dimension greater than 1
        points are determined through enumeration up to the specified bound.

        TODO: modify implementation from projective space
        """
        raise NotImplementedError("enumerating points on weighted projective scheme is not implemented")


class SchemeHomset_points_weighted_projective_field(SchemeHomset_points):
    """
    Set of rational points of a weighted projective variety over a field.

    Placeholder class.

    EXAMPLES::

        sage: WP = WeightedProjectiveSpace(QQ, [1, 3, 1]); WP
        Weighted Projective Space of dimension 2 with weights (1, 3, 1) over Rational Field
        sage: WP.point_homset()
        Set of rational points of Weighted Projective Space of dimension 2 with weights (1, 3, 1) over Rational Field
        sage: type(WP.point_homset())
        <class 'sage.schemes.weighted_projective.weighted_projective_homset.SchemeHomset_points_weighted_projective_field_with_category'>
        sage: WP(2, 24, 2)
        (1 : 3 : 1)
        sage: WP(2, 24, 2) == WP(1, 3, 1)
        True
    """
    pass
