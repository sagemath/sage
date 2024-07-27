# sage_setup: distribution = sagemath-modules
r"""
Top level of the distribution package sagemath-modules

This distribution makes the following features available::

    sage: from sage.features.sagemath import *
    sage: sage__modules().is_present()
    FeatureTestResult('sage.modules', True)
    sage: sage__rings__real_mpfr().is_present()
    FeatureTestResult('sage.rings.real_mpfr', True)
    sage: sage__rings__complex_double().is_present()
    FeatureTestResult('sage.rings.complex_double', True)
"""

from .all__sagemath_categories import *

try:  # extra
    from sage.all__sagemath_flint import *
except ImportError:
    pass

try:  # extra
    from sage.all__sagemath_pari import *
except ImportError:
    pass

from sage.misc.all__sagemath_modules import *
from sage.rings.all__sagemath_modules import *
from sage.combinat.all__sagemath_modules import *
from sage.algebras.all__sagemath_modules import *
from sage.modules.all import *
from sage.matrix.all import *
from sage.groups.all__sagemath_modules import *
from sage.geometry.all__sagemath_modules import *
from sage.homology.all__sagemath_modules import *
from sage.tensor.all import *
from sage.matroids.all import *
from sage.quadratic_forms.all__sagemath_modules import *
from sage.coding.all import *
from sage.crypto.all import *
from sage.stats.all import *
from sage.probability.all import *
from sage.calculus.all__sagemath_modules import *
from sage.numerical.all__sagemath_modules import *

import sage.crypto.mq as mq
import sage.stats.all as stats

true = True
false = False
