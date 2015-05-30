"""
Deprecated preparser module

TESTS::

    sage: import sage.misc.preparser
    sage: sage.misc.preparser.is_loadable_filename('foo.sage')
    doctest:...: DeprecationWarning: 
    Importing is_loadable_filename from here is deprecated. If you need to use it, please import it directly from sage.repl.load
    See http://trac.sagemath.org/17396 for details.
    True
"""

from sage.misc.lazy_import import lazy_import

lazy_import('sage.repl.preparse', '*', deprecation=17396)
lazy_import('sage.repl.load', '*', deprecation=17396)
