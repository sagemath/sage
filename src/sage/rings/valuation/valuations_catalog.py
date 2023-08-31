from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.padics.padic_valuation', 'pAdicValuation')
lazy_import('sage.rings.function_field.valuation', 'FunctionFieldValuation')
lazy_import('sage.rings.valuation.gauss_valuation', 'GaussValuation')
lazy_import('sage.rings.valuation.trivial_valuation', ['TrivialDiscretePseudoValuation', 'TrivialPseudoValuation', 'TrivialValuation'])
lazy_import('sage.rings.valuation.limit_valuation', 'LimitValuation')

del lazy_import
