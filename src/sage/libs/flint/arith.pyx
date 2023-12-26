# sage_setup: distribution = sagemath-flint
# distutils: extra_compile_args = -D_XPG6
"""
FLINT Arithmetic Functions
"""

from sage.misc.lazy_import import LazyImport

bell_number = LazyImport('sage.libs.flint.arith_sage', 'bell_number', deprecation=36449)
bernoulli_number = LazyImport('sage.libs.flint.arith_sage', 'bernoulli_number', deprecation=36449)
euler_number = LazyImport('sage.libs.flint.arith_sage', 'euler_number', deprecation=36449)
stirling_number_1 = LazyImport('sage.libs.flint.arith_sage', 'stirling_number_1', deprecation=36449)
stirling_number_2 = LazyImport('sage.libs.flint.arith_sage', 'stirling_number_2', deprecation=36449)
number_of_partitions = LazyImport('sage.libs.flint.arith_sage', 'number_of_partitions', deprecation=36449)
dedekind_sum = LazyImport('sage.libs.flint.arith_sage', 'dedekind_sum', deprecation=36449)
harmonic_number = LazyImport('sage.libs.flint.arith_sage', 'harmonic_number', deprecation=36449)
