from sage.misc.lazy_import import lazy_import

import sage.algebras.catalog as algebras

from .quatalg.all import *
from .steenrod.all import *
from .fusion_rings.all import *
from .lie_algebras.all import *
from .quantum_groups.all import *
from .lie_conformal_algebras.all import *

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

from .clifford_algebra import CliffordAlgebra, ExteriorAlgebra
from .weyl_algebra import DifferentialWeylAlgebra

lazy_import('sage.algebras.commutative_dga', 'GradedCommutativeAlgebra')

lazy_import('sage.algebras.rational_cherednik_algebra',
            'RationalCherednikAlgebra')

lazy_import('sage.algebras.tensor_algebra', 'TensorAlgebra')

lazy_import('sage.algebras.q_system', 'QSystem')

lazy_import('sage.algebras.cluster_algebra', 'ClusterAlgebra')

lazy_import('sage.algebras.yangian', 'Yangian')
