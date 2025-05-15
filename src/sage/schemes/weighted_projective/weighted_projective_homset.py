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

    EXAMPLES::

        sage: W = WeightedProjectiveSpace([1, 3, 1], QQ)
        sage: W.an_element().parent() is W.point_homset()
        True
        sage: W.point_homset()
        Set of rational points of Weighted Projective Space of dimension 2 with weights (1, 3, 1) over Rational Field
        sage: type(W.point_homset())
        <class 'sage.schemes.weighted_projective.weighted_projective_homset.SchemeHomset_points_weighted_projective_ring_with_category'>
        sage: W(2, 24, 2)
        (1 : 3 : 1)
        sage: W(2, 24, 2) == W(1, 3, 1)
        True

    ::

        sage: R.<x> = QQ[]
        sage: H = HyperellipticCurveSmoothModel(x^6 + x - 17); H
        Hyperelliptic Curve over Rational Field defined by y^2 = x^6 + x - 17
        Page: P = H(2, -7); P
        (2 : -7 : 1)
        sage: type(P)
        <class 'sage.schemes.weighted_projective.weighted_projective_point.SchemeMorphism_point_weighted_projective_ring'>
    """
