r"""
Deprecated in favor of :mod:`sage.rings.complex_mpfr`

TESTS::

    sage: from sage.rings.complex_number import ComplexNumber
    doctest:warning
    ...
    DeprecationWarning: the complex_number module is deprecated, please use sage.rings.complex_mpfr
    See http://trac.sagemath.org/24483 for details.
"""

from sage.misc.superseded import deprecation
from sage.rings.complex_mpfr import *
deprecation(24483, "the complex_number module is deprecated, please use sage.rings.complex_mpfr")

