from .all__sagemath_combinat import *
from .all__sagemath_gap import *
from .all__sagemath_ntl import *
from .all__sagemath_pari import *

try:
    from .all__sagemath_symbolics import *
except ImportError:
    pass

from sage.misc.lazy_import import lazy_import

lazy_import('sage.libs.eclib.constructor', 'CremonaModularSymbols')
lazy_import('sage.libs.eclib.interface', ['mwrank_EllipticCurve', 'mwrank_MordellWeil'])
lazy_import('sage.libs.eclib.mwrank', 'get_precision', 'mwrank_get_precision')
lazy_import('sage.libs.eclib.mwrank', 'set_precision', 'mwrank_set_precision')
lazy_import('sage.libs.eclib.mwrank', 'initprimes', 'mwrank_initprimes')

lazy_import('sage.libs.flint.qsieve', 'qsieve')

del lazy_import
