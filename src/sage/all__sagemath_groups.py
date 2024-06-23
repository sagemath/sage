# sage_setup: distribution = sagemath-groups
r"""
Top level of the distribution package sagemath-groups

This distribution makes the following feature available::

    sage: from sage.features.sagemath import *
    sage: sage__groups().is_present()
    FeatureTestResult('sage.groups', True)
"""

from .all__sagemath_modules import *
from .all__sagemath_gap import *

try:  # extra
    from sage.all__sagemath_combinat import *
except ImportError:
    pass

from sage.groups.all__sagemath_groups import *
