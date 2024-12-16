"""
Algebras
"""
# ****************************************************************************
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
# ****************************************************************************
from sage.misc.lazy_import import lazy_import

# old-style class for associative algebras, use Parent instead
from sage.rings.ring import Algebra

import sage.algebras.catalog as algebras

from sage.algebras.quatalg.all import *
from sage.algebras.steenrod.all import *
from sage.algebras.fusion_rings.all import *
from sage.algebras.lie_algebras.all import *
from sage.algebras.quantum_groups.all import *
from sage.algebras.lie_conformal_algebras.all import *

# Algebra base classes
from sage.algebras.free_algebra import FreeAlgebra
from sage.algebras.free_algebra_quotient import FreeAlgebraQuotient

from sage.algebras.finite_dimensional_algebras.all import FiniteDimensionalAlgebra

lazy_import('sage.algebras.group_algebra', 'GroupAlgebra')

lazy_import('sage.algebras.iwahori_hecke_algebra', 'IwahoriHeckeAlgebra')
lazy_import('sage.algebras.affine_nil_temperley_lieb',
            'AffineNilTemperleyLiebTypeA')
lazy_import('sage.algebras.nil_coxeter_algebra', 'NilCoxeterAlgebra')
lazy_import('sage.algebras.schur_algebra', ['SchurAlgebra',
                                            'SchurTensorModule'])

lazy_import('sage.algebras.hall_algebra', 'HallAlgebra')

lazy_import('sage.algebras.jordan_algebra', 'JordanAlgebra')

lazy_import('sage.algebras.octonion_algebra', 'OctonionAlgebra')

lazy_import('sage.algebras.shuffle_algebra', 'ShuffleAlgebra')

from sage.algebras.clifford_algebra import CliffordAlgebra, ExteriorAlgebra
from sage.algebras.weyl_algebra import DifferentialWeylAlgebra

lazy_import('sage.algebras.commutative_dga', 'GradedCommutativeAlgebra')

lazy_import('sage.algebras.rational_cherednik_algebra',
            'RationalCherednikAlgebra')

lazy_import('sage.algebras.tensor_algebra', 'TensorAlgebra')

lazy_import('sage.algebras.q_system', 'QSystem')

lazy_import('sage.algebras.cluster_algebra', 'ClusterAlgebra')

lazy_import('sage.algebras.yangian', 'Yangian')

lazy_import('sage.algebras.flag_algebras', ['CombinatorialTheory', 
                                            'GraphTheory', 
                                            'DiGraphTheory', 
                                            'ThreeGraphTheory', 
                                            'TournamentTheory', 
                                            'PermutationTheory', 
                                            'OEGraphTheory', 
                                            'OVGraphTheory', 
                                            'RamseyGraphTheory'])
