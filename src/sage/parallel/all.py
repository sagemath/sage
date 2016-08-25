from __future__ import absolute_import

from .decorate import parallel, fork
from sage.misc.lazy_import import lazy_import
lazy_import('sage.parallel.parallelism', 'Parallelism')
