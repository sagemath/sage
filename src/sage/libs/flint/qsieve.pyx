r"""
Deprecated module.

Functions were moved in qsieve_sage.pyx.

TESTS::

    sage: from sage.libs.flint.qsieve import qsieve
    sage: qsieve(1000)
    doctest:warning
    ...
    DeprecationWarning:
    Importing qsieve from here is deprecated; please use "from sage.libs.flint.qsieve_sage import qsieve" instead.
    See https://github.com/sagemath/sage/issues/36449 for details.
    [(2, 3), (5, 3)]
"""

from sage.misc.lazy_import import LazyImport

qsieve = LazyImport('sage.libs.flint.qsieve_sage', 'qsieve', deprecation=36449)
