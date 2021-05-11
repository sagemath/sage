r"""
Deprecated in favor of :mod:`sage.rings.complex_mpfr`

TESTS::

    sage: from sage.rings.complex_field import ComplexField
    doctest:warning
    ...
    DeprecationWarning: the complex_field module is deprecated, please use sage.rings.complex_mpfr
    See http://trac.sagemath.org/24483 for details.
    sage: ComplexField()
    Complex Field with 53 bits of precision
"""

from sage.misc.superseded import deprecation
from sage.rings.complex_mpfr import *
deprecation(24483, "the complex_field module is deprecated, please use sage.rings.complex_mpfr")
