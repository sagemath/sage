# sage_setup: distribution = sagemath-pari
from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.valuation.gauss_valuation', 'GaussValuation')
lazy_import('sage.rings.valuation', 'valuations_catalog', 'valuations')
lazy_import('sage.rings.valuation.value_group', 'DiscreteValueGroup')
del lazy_import
