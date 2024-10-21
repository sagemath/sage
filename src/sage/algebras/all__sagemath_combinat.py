# sage_setup: distribution = sagemath-combinat
from sage.misc.lazy_import import lazy_import

# Algebra base classes
lazy_import('sage.algebras.free_algebra', 'FreeAlgebra')
lazy_import('sage.algebras.free_algebra_quotient', 'FreeAlgebraQuotient')

from sage.algebras.steenrod.all import *
from sage.algebras.quantum_groups.all import *

lazy_import('sage.algebras.iwahori_hecke_algebra', 'IwahoriHeckeAlgebra')
lazy_import('sage.algebras.affine_nil_temperley_lieb',
            'AffineNilTemperleyLiebTypeA')
lazy_import('sage.algebras.nil_coxeter_algebra', 'NilCoxeterAlgebra')
lazy_import('sage.algebras.schur_algebra', ['SchurAlgebra',
                                            'SchurTensorModule'])

lazy_import('sage.algebras.hall_algebra', 'HallAlgebra')

lazy_import('sage.algebras.jordan_algebra', 'JordanAlgebra')

lazy_import('sage.algebras.shuffle_algebra', 'ShuffleAlgebra')

lazy_import('sage.algebras.commutative_dga', 'GradedCommutativeAlgebra')

lazy_import('sage.algebras.rational_cherednik_algebra',
            'RationalCherednikAlgebra')

lazy_import('sage.algebras.tensor_algebra', 'TensorAlgebra')

lazy_import('sage.algebras.q_system', 'QSystem')

lazy_import('sage.algebras.cluster_algebra', 'ClusterAlgebra')

lazy_import('sage.algebras.yangian', 'Yangian')

del lazy_import
