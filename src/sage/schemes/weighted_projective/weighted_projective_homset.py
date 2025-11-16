"""
Hom-sets of weighted projective schemes

AUTHORS:

- Gareth Ma (2024): initial version, based on unweighted version.
"""

# *****************************************************************************
#        Copyright (C) 2006 William Stein <wstein@gmail.com>
#        Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#        Copyright (C) 2024 Gareth Ma <grhkm21@gmail.com>
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

    See :class:`SchemeHomset_points`.

    EXAMPLES::

        sage: W = WeightedProjectiveSpace([3, 4, 5], QQ)
        sage: W.point_homset()
        Set of rational points of Weighted Projective Space of dimension 2 with weights (3, 4, 5) over Rational Field
        sage: W.an_element().parent() is W.point_homset()
        True
    """
