# sage_setup: distribution = sagemath-categories

from sage.parallel.decorate import parallel, fork
from sage.misc.lazy_import import lazy_import
lazy_import('sage.parallel.parallelism', 'Parallelism')
del lazy_import
