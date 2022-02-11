"""
The Steenrod algebra
"""
from .steenrod_algebra import SteenrodAlgebra, Sq
from sage.misc.lazy_import import lazy_import
lazy_import('sage.algebras.steenrod.steenrod_algebra_bases',
            'steenrod_algebra_basis',
            deprecation=(32647, 'removed from namespace'))
