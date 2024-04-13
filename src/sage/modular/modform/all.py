# sage_setup: distribution = sagemath-schemes
#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
#########################################################################

from sage.modular.modform.constructor import ModularForms, CuspForms, EisensteinForms, Newforms, Newform

from sage.modular.modform.eis_series import eisenstein_series_qexp, eisenstein_series_lseries

from sage.modular.modform.half_integral import half_integral_weight_modform_basis

from sage.modular.modform.theta import theta_qexp, theta2_qexp

from sage.misc.lazy_import import lazy_import

lazy_import('sage.modular.modform.j_invariant', 'j_invariant_qexp')
lazy_import('sage.modular.modform.vm_basis', ['victor_miller_basis', 'delta_qexp'])

from sage.modular.modform.hecke_operator_on_qexp import hecke_operator_on_qexp, hecke_operator_on_basis

from sage.modular.modform.numerical import NumericalEigenforms as numerical_eigenforms

from sage.modular.modform.element import delta_lseries

from sage.modular.modform.ring import ModularFormsRing
del lazy_import
