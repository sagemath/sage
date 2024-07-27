# sage_setup: distribution = sagemath-schemes
"""
Plane curves
"""

# *****************************************************************************
#
#   Sage: Open Source Mathematical Software
#
#       Copyright (C) 2005 William Stein <was@math.harvard.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.schemes.curves.constructor import Curve
from sage.schemes.curves.projective_curve import Hasse_bounds

from sage.misc.lazy_import import lazy_import

lazy_import('sage.schemes.curves.plane_curve_arrangement', 'PlaneCurveArrangements')

lazy_import('sage.schemes.curves.plane_curve_arrangement', 'AffinePlaneCurveArrangements')

lazy_import('sage.schemes.curves.plane_curve_arrangement', 'ProjectivePlaneCurveArrangements')
