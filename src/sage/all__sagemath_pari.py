# sage_setup: distribution = sagemath-pari
r"""
Top level of the distribution package sagemath-pari

This distribution makes the following features available::

    sage: from sage.features.sagemath import *
    sage: sage__libs__pari().is_present()
    FeatureTestResult('sage.libs.pari', True)
    sage: sage__rings__finite_rings().is_present()
    FeatureTestResult('sage.rings.finite_rings', True)
"""

from .all__sagemath_categories import *

from sage.libs.all__sagemath_pari import *
