# code exports

from .fano_variety import CPRFanoToricVariety
from .ideal import ToricIdeal
from .library import toric_varieties
from .variety import AffineToricVariety, ToricVariety


from sage.misc.lazy_import import lazy_import
lazy_import('sage.schemes.toric.weierstrass', 'WeierstrassForm')
