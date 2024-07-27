# sage_setup: distribution = sagemath-flint
from sage.misc.lazy_import import lazy_import
lazy_import('sage.rings.polynomial.weil.weil_polynomials', 'WeilPolynomials')
del lazy_import
