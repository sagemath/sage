# sage_setup: distribution = sagemath-flint

r"""
Deprecated module.

Functions were moved in arith_sage.pyx

TESTS::

    sage: from sage.libs.flint.arith import bell_number, bernoulli_number, euler_number, stirling_number_1, stirling_number_2, number_of_partitions, dedekind_sum, harmonic_number
    sage: bell_number(4)
    doctest:warning
    ...
    DeprecationWarning:
    Importing bell_number from here is deprecated; please use "from sage.libs.flint.arith_sage import bell_number" instead.
    See https://github.com/sagemath/sage/issues/36449 for details.
    15
    sage: bernoulli_number(4)
    doctest:warning
    ...
    DeprecationWarning:
    Importing bernoulli_number from here is deprecated; please use "from sage.libs.flint.arith_sage import bernoulli_number" instead.
    See https://github.com/sagemath/sage/issues/36449 for details.
    -1/30
    sage: euler_number(4)
    doctest:warning
    ...
    DeprecationWarning:
    Importing euler_number from here is deprecated; please use "from sage.libs.flint.arith_sage import euler_number" instead.
    See https://github.com/sagemath/sage/issues/36449 for details.
    5
    sage: stirling_number_1(2, 4)
    doctest:warning
    ...
    DeprecationWarning:
    Importing stirling_number_1 from here is deprecated; please use "from sage.libs.flint.arith_sage import stirling_number_1" instead.
    See https://github.com/sagemath/sage/issues/36449 for details.
    0
    sage: stirling_number_2(2, 4)
    doctest:warning
    ...
    DeprecationWarning:
    Importing stirling_number_2 from here is deprecated; please use "from sage.libs.flint.arith_sage import stirling_number_2" instead.
    See https://github.com/sagemath/sage/issues/36449 for details.
    0
    sage: number_of_partitions(4)
    doctest:warning
    ...
    DeprecationWarning:
    Importing number_of_partitions from here is deprecated; please use "from sage.libs.flint.arith_sage import number_of_partitions" instead.
    See https://github.com/sagemath/sage/issues/36449 for details.
    5
    sage: dedekind_sum(4, 5)
    doctest:warning
    ...
    DeprecationWarning:
    Importing dedekind_sum from here is deprecated; please use "from sage.libs.flint.arith_sage import dedekind_sum" instead.
    See https://github.com/sagemath/sage/issues/36449 for details.
    -1/5
    sage: harmonic_number(4)
    doctest:warning
    ...
    DeprecationWarning:
    Importing harmonic_number from here is deprecated; please use "from sage.libs.flint.arith_sage import harmonic_number" instead.
    See https://github.com/sagemath/sage/issues/36449 for details.
    25/12
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
