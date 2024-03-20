# sage_setup: distribution = sagemath-flint
r"""
Top level of the distribution package sagemath-flint

This distribution makes the following features available::

    sage: from sage.features.sagemath import *
    sage: sage__libs__flint().is_present()
    FeatureTestResult('sage.libs.flint', True)
    sage: sage__rings__complex_interval_field().is_present()
    FeatureTestResult('sage.rings.complex_interval_field', True)
    sage: sage__rings__number_field().is_present()
    FeatureTestResult('sage.rings.number_field', True)
    sage: sage__rings__real_interval_field().is_present()
    FeatureTestResult('sage.rings.real_interval_field', True)
"""

from .all__sagemath_ntl import *

from .libs.all__sagemath_flint import *

from .rings.all__sagemath_flint import *

from sage.rings.qqbar import _init_qqbar
_init_qqbar()
