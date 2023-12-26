# sage_setup: distribution = sagemath-combinat
r"""
Top level of the distribution package sagemath-combinat

This distribution makes the following features available::

    sage: from sage.features.sagemath import *
    sage: sage__combinat().is_present()
    FeatureTestResult('sage.combinat', True)
    sage: sage__sat().is_present()
    FeatureTestResult('sage.sat', True)
"""

try:  # extra
    from sage.all__sagemath_modules import *
except ImportError:
    pass

from .all__sagemath_categories import *

try:  # extra
    from sage.all__sagemath_graphs import *
except ImportError:
    pass

from sage.algebras.all__sagemath_combinat import *
from sage.combinat.all__sagemath_combinat import *
from sage.libs.all__sagemath_combinat import *
from sage.rings.all__sagemath_combinat import *
from sage.monoids.all import *
from sage.games.all import *
from sage.sat.all import *
