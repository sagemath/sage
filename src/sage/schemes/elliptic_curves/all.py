"""
Exported elliptic curves functionality
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .constructor import (EllipticCurve,
                         EllipticCurve_from_c4c6,
                         EllipticCurve_from_j,
                         EllipticCurve_from_cubic,
                         EllipticCurves_with_good_reduction_outside_S)

from sage.misc.lazy_import import lazy_import
lazy_import('sage.schemes.elliptic_curves.jacobian', 'Jacobian')

lazy_import('sage.schemes.elliptic_curves.ell_rational_field',
            ['cremona_curves', 'cremona_optimal_curves'])

from .cm import ( cm_orders,
                 cm_j_invariants,
                 cm_j_invariants_and_orders,
                 hilbert_class_polynomial )

lazy_import('sage.schemes.elliptic_curves.ec_database', 'elliptic_curves')

from .kodaira_symbol import KodairaSymbol

from .ell_curve_isogeny import EllipticCurveIsogeny, isogeny_codomain_from_kernel

from .heegner import heegner_points, heegner_point
