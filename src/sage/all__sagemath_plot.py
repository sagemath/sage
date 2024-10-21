# sage_setup: distribution = sagemath-plot
r"""
Top level of the distribution package sagemath-plot

This distribution makes the following feature available::

    sage: from sage.features.sagemath import *
    sage: sage__plot().is_present()
    FeatureTestResult('sage.plot', True)
"""
from .all__sagemath_modules import *

from sage.plot.all import *
from sage.plot.plot3d.all     import *
