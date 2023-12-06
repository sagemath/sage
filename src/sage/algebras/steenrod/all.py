# sage_setup: distribution = sagemath-combinat
"""
The Steenrod algebra
"""
from sage.misc.lazy_import import lazy_import
lazy_import('sage.algebras.steenrod.steenrod_algebra', ['SteenrodAlgebra', 'Sq'])
lazy_import('sage.algebras.steenrod.steenrod_algebra_bases',
            'steenrod_algebra_basis',
            deprecation=(32647, 'removed from namespace'))
del lazy_import
