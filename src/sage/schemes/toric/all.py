from sage.misc.lazy_import import lazy_import

lazy_import('sage.schemes.toric.weierstrass', 'WeierstrassForm')
lazy_import('sage.schemes.toric.variety', ['AffineToricVariety', 'ToricVariety'])
lazy_import('sage.schemes.toric.library', 'toric_varieties')
lazy_import('sage.schemes.toric.fano_variety', 'CPRFanoToricVariety')
lazy_import('sage.schemes.toric.ideal', 'ToricIdeal')
