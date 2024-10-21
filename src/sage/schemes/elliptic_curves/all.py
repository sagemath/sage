# sage_setup: distribution = sagemath-schemes
"""
Exported elliptic curves functionality
"""

# *****************************************************************************
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
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.schemes.elliptic_curves.constructor import (EllipticCurve,
                                                      EllipticCurve_from_c4c6,
                                                      EllipticCurve_from_j,
                                                      EllipticCurve_from_cubic,
                                                      EllipticCurves_with_good_reduction_outside_S)

from sage.misc.lazy_import import lazy_import
lazy_import('sage.schemes.elliptic_curves.jacobian', 'Jacobian')

lazy_import('sage.schemes.elliptic_curves.ell_finite_field', 'special_supersingular_curve')
lazy_import('sage.schemes.elliptic_curves.ell_rational_field',
            ['cremona_curves', 'cremona_optimal_curves'])

from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_with_prime_order
from sage.schemes.elliptic_curves.cm import (cm_orders,
                                             cm_j_invariants,
                                             cm_j_invariants_and_orders,
                                             hilbert_class_polynomial)

lazy_import('sage.schemes.elliptic_curves.ec_database', 'elliptic_curves')

from sage.schemes.elliptic_curves.kodaira_symbol import KodairaSymbol

from sage.schemes.elliptic_curves.ell_curve_isogeny import EllipticCurveIsogeny, isogeny_codomain_from_kernel

lazy_import('sage.schemes.elliptic_curves.mod_poly', 'classical_modular_polynomial')

from sage.schemes.elliptic_curves.heegner import heegner_points, heegner_point
del lazy_import
