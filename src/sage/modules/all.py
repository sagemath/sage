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

from sage.modules.free_module import FreeModule, VectorSpace, span

from sage.modules.free_quadratic_module import (
    FreeQuadraticModule,
    QuadraticSpace,
    InnerProductSpace,
)

from sage.modules.free_module_element import (
    vector,
    free_module_element,
    zero_vector,
    random_vector,
)

from sage.modules.vector_space_morphism import linear_transformation

from sage.modules.with_basis.all import *

from sage.misc.lazy_import import lazy_import

lazy_import("sage.modules.filtered_vector_space", "FilteredVectorSpace")
lazy_import("sage.modules.multi_filtered_vector_space", "MultiFilteredVectorSpace")
lazy_import("sage.modules.free_quadratic_module_integer_symmetric", "IntegralLattice")
lazy_import("sage.modules.torsion_quadratic_module", "TorsionQuadraticForm")
del lazy_import
