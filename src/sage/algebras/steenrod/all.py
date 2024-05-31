# sage_setup: distribution = sagemath-combinat
"""
The Steenrod algebra
"""
from sage.misc.lazy_import import lazy_import
lazy_import('sage.algebras.steenrod.steenrod_algebra', ['SteenrodAlgebra', 'Sq'])
del lazy_import
