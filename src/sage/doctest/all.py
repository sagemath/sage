# sage_setup: distribution = sagemath-repl
from sage.misc.lazy_import import lazy_import
lazy_import('sage.doctest.control', 'run_doctests')
del lazy_import
