# sage_setup: distribution = sagemath-eclib

from sage.misc.lazy_import import lazy_import

lazy_import('sage.libs.eclib.constructor', 'CremonaModularSymbols')
lazy_import('sage.libs.eclib.interface', ['mwrank_EllipticCurve', 'mwrank_MordellWeil'])
lazy_import('sage.libs.eclib.mwrank', 'get_precision', 'mwrank_get_precision')
lazy_import('sage.libs.eclib.mwrank', 'set_precision', 'mwrank_set_precision')
lazy_import('sage.libs.eclib.mwrank', 'initprimes', 'mwrank_initprimes')

del lazy_import
