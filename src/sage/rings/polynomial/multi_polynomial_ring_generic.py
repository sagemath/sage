"""
TESTS::

    sage: R.<x,y> = QQbar[]
    sage: from sage.rings.polynomial.multi_polynomial_ring_generic import MPolynomialRing_generic
    doctest:...: DeprecationWarning: the module sage.rings.polynomial.multi_polynomial_ring_generic is deprecated, use sage.rings.polynomial.multi_polynomial_ring_base instead.
    See https://trac.sagemath.org/25563 for details.
    sage: isinstance(R, MPolynomialRing_generic)
    True
"""

from sage.misc.superseded import deprecation
deprecation(25563, "the module sage.rings.polynomial.multi_polynomial_ring_generic is deprecated, "
                   "use sage.rings.polynomial.multi_polynomial_ring_base instead.")

from .multi_polynomial_ring_base import MPolynomialRing_base
MPolynomialRing_generic = MPolynomialRing_base
