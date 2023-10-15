r"""
Deprecated module

TESTS::

    sage: from sage.libs.flint.fmpz_poly import Fmpz_poly
    sage: Fmpz_poly([1, 1])
    doctest:warning
    ...
    DeprecationWarning:
    Importing Fmpz_poly from here is deprecated; please use "from sage.libs.flint.fmpz_poly_sage import Fmpz_poly" instead.
    See https://github.com/sagemath/sage/issues/36449 for details.
    2  1 1
"""

from sage.misc.lazy_import import lazy_import

lazy_import('sage.libs.flint.fmpz_poly_sage', 'Fmpz_poly',  deprecation=36449)

del lazy_import
