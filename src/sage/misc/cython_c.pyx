"""
TESTS::

    sage: from sage.misc.cython_c import cython_compile
    doctest:...: DeprecationWarning: the module sage.misc.cython_c has been moved to sage.misc.cython
    See http://trac.sagemath.org/24105 for details.
    sage: cython_compile('''print("Hello World!")''')
    Hello World!
"""
from sage.misc.superseded import deprecation
deprecation(24105, "the module sage.misc.cython_c has been moved to sage.misc.cython")

from .cython import cython_compile
