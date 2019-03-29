"""
TESTS:

Test the deprecation warnings::

    sage: tests.CompleteMatchings
    doctest:warning
    ...
    DeprecationWarning:
    Importing CompleteMatchings from here is deprecated. If you need to use it, please import it directly from sage.tests.arxiv_0812_2725
    See https://trac.sagemath.org/27337 for details.
    <function CompleteMatchings at ...>
    sage: tests.modsym
    doctest:warning
    ...
    DeprecationWarning:
    Importing modsym from here is deprecated. If you need to use it, please import it directly from sage.modular.modsym.tests
    See https://trac.sagemath.org/27337 for details.
    <class ...sage.modular.modsym.tests.Test...>
"""

from sage.misc.lazy_import import lazy_import
lazy_import('sage.modular.modsym.tests', 'Test', as_='modsym', deprecation=27337)
lazy_import('sage.tests.arxiv_0812_2725', '*', deprecation=27337)
