r"""
Deprecated import from :trac:`24483`.

TESTS::

    sage: from sage.rings.complex_field import ComplexField
    sage: ComplexField()
    doctest:warning
    ...
    Importing ComplexField from here is deprecated. If you need to use it, please import it directly from sage.rings.complex_mpfr
    See http://trac.sagemath.org/24483 for details.
    Complex Field with 53 bits of precision
"""

from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.complex_mpfr', 'ComplexField', deprecation=24483)

