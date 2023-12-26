# sage_setup: distribution = sagemath-combinat

from sage.rings.all__sagemath_categories import *

from sage.misc.lazy_import import lazy_import

# Lazy Laurent series ring
lazy_import('sage.rings.lazy_series_ring', ['LazyLaurentSeriesRing', 'LazyPowerSeriesRing',
                                            'LazySymmetricFunctions', 'LazyDirichletSeriesRing'])

del lazy_import
