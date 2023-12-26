# sage_setup: distribution = sagemath-flint
"""
Interface to FLINT's ``qsieve_factor()``. This used to interact
with an external "QuadraticSieve" program, but its functionality has
been absorbed into flint.
"""

from sage.misc.lazy_import import LazyImport

qsieve = LazyImport('sage.libs.flint.qsieve_sage', 'qsieve', deprecation=36449)
