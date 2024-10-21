# sage_setup: distribution = sagemath-flint

r"""
Deprecated modules.

Functions were moved in ulong_extras_sage.pyx

TESTS::

    sage: from sage.libs.flint.ulong_extras import n_factor_to_list
    sage: n_factor_to_list(60, 20)
    doctest:warning
    ...
    DeprecationWarning:
    Importing n_factor_to_list from here is deprecated; please use "from sage.libs.flint.ulong_extras_sage import n_factor_to_list" instead.
    See https://github.com/sagemath/sage/issues/36449 for details.
    [(2, 2), (3, 1), (5, 1)]
"""

from sage.misc.lazy_import import LazyImport

n_factor_to_list = LazyImport('sage.libs.flint.ulong_extras_sage', 'n_factor_to_list', deprecation=36449)
