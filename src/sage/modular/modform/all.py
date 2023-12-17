#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

from .constructor import ModularForms, CuspForms, EisensteinForms, Newforms, Newform

from .eis_series import eisenstein_series_qexp, eisenstein_series_lseries

from .half_integral import half_integral_weight_modform_basis

from .theta import theta_qexp, theta2_qexp

from sage.misc.lazy_import import lazy_import

lazy_import('sage.modular.modform.j_invariant', 'j_invariant_qexp')
lazy_import('sage.modular.modform.vm_basis', ['victor_miller_basis', 'delta_qexp'])

from .hecke_operator_on_qexp import hecke_operator_on_qexp, hecke_operator_on_basis

from .numerical import NumericalEigenforms as numerical_eigenforms

from .element import delta_lseries

from .ring import ModularFormsRing
