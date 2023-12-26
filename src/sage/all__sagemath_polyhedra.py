# sage_setup: distribution = sagemath-polyhedra
r"""
Top level of the distribution package sagemath-polyhedra

This distribution makes the following features available::

    sage: from sage.features.sagemath import *
    sage: sage__geometry__polyhedron().is_present()
    FeatureTestResult('sage.geometry.polyhedron', True)
    sage: sage__numerical__mip().is_present()
    FeatureTestResult('sage.numerical.mip', True)
"""

from .all__sagemath_modules import *

try:  # extra
    from sage.all__sagemath_combinat import *
except ImportError:
    pass

try:  # extra
    from sage.all__sagemath_graphs import *
except ImportError:
    pass

try:  # extra
    from sage.all__sagemath_groups import *
except ImportError:
    pass

try:  # extra
    from sage.all__sagemath_flint import *
except ImportError:
    pass

try:  # extra
    from sage.all__sagemath_plot import *
except ImportError:
    pass

from sage.geometry.all__sagemath_polyhedra import *
from sage.geometry.triangulation.all import *
from sage.numerical.all import *
from sage.game_theory.all import *
from sage.schemes.all__sagemath_polyhedra import *
