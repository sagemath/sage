# sage_setup: distribution = sagemath-gap

from sage.misc.lazy_import import lazy_import

lazy_import('sage.geometry.ribbon_graph', 'RibbonGraph')
del lazy_import
