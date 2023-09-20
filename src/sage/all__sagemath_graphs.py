# sage_setup: distribution = sagemath-graphs
r"""
Top level of the distribution package sagemath-graphs

This distribution makes the following feature available::

    sage: from sage.features.sagemath import *
    sage: sage__graphs().is_present()
    FeatureTestResult('sage.graphs', True)
"""

from .all__sagemath_categories import *

try:  # extra
    from sage.all__sagemath_modules import *
except ImportError:
    pass

try:  # extra
    from sage.all__sagemath_polyhedra import *
except ImportError:
    pass

from sage.graphs.all     import *

from sage.topology.all   import *

from sage.combinat.all__sagemath_graphs import *

from sage.sandpiles.all import *

from sage.knots.all import *
