"""
Algebras
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
from sage.misc.lazy_import import lazy_import

import sage.algebras.catalog as algebras

from .quantum_groups.all import *
from .quatalg.all import *

# Algebra base classes
from .algebra import Algebra
from .free_algebra import FreeAlgebra
from .free_algebra_quotient import FreeAlgebraQuotient

from .steenrod.all import *
from .lie_algebras.all import *
from .quantum_groups.all import *
from .lie_conformal_algebras.all import *

from .finite_dimensional_algebras.all import FiniteDimensionalAlgebra

lazy_import('sage.algebras.group_algebra', 'GroupAlgebra')

lazy_import('sage.algebras.iwahori_hecke_algebra', 'IwahoriHeckeAlgebra')
lazy_import('sage.algebras.affine_nil_temperley_lieb',
            'AffineNilTemperleyLiebTypeA')
lazy_import('sage.algebras.nil_coxeter_algebra', 'NilCoxeterAlgebra')
lazy_import('sage.algebras.schur_algebra', ['SchurAlgebra',
                                            'SchurTensorModule'])

lazy_import('sage.algebras.hall_algebra', 'HallAlgebra')

lazy_import('sage.algebras.jordan_algebra', 'JordanAlgebra')

lazy_import('sage.algebras.shuffle_algebra', 'ShuffleAlgebra')

from .clifford_algebra import CliffordAlgebra, ExteriorAlgebra
from .weyl_algebra import DifferentialWeylAlgebra

lazy_import('sage.algebras.commutative_dga', 'GradedCommutativeAlgebra')

lazy_import('sage.algebras.rational_cherednik_algebra',
            'RationalCherednikAlgebra')

lazy_import('sage.algebras.tensor_algebra', 'TensorAlgebra')

lazy_import('sage.algebras.q_system', 'QSystem')

lazy_import('sage.algebras.cluster_algebra', 'ClusterAlgebra')

lazy_import('sage.algebras.yangian', 'Yangian')
