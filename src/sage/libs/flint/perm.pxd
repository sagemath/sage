# distutils: libraries = flint
# distutils: depends = flint/perm.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    long * _perm_init(long n)
    # Initialises the permutation for use.

    void _perm_clear(long *vec)
    # Clears the permutation.

    void _perm_set(long *res, const long *vec, long n)
    # Sets the permutation ``res`` to the same as the permutation ``vec``.

    void _perm_set_one(long *vec, long n)
    # Sets the permutation to the identity permutation.

    void _perm_inv(long *res, const long *vec, long n)
    # Sets ``res`` to the inverse permutation of ``vec``.
    # Allows aliasing of ``res`` and ``vec``.

    void _perm_compose(long *res, const long *vec1, const long *vec2, long n)
    # Forms the composition `\pi_1 \circ \pi_2` of two permutations
    # `\pi_1` and `\pi_2`.  Here, `\pi_2` is applied first, that is,
    # `(\pi_1 \circ \pi_2)(i) = \pi_1(\pi_2(i))`.
    # Allows aliasing of ``res``, ``vec1`` and ``vec2``.

    int _perm_parity(const long *vec, long n)
    # Returns the parity of ``vec``, 0 if the permutation is even and 1 if
    # the permutation is odd.

    int _perm_randtest(long *vec, long n, flint_rand_t state)
    # Generates a random permutation vector of length `n` and returns
    # its parity, 0 or 1.
    # This function uses the Knuth shuffle algorithm to generate a uniformly
    # random permutation without retries.

    int _perm_print(const long * vec, long n)
    # Prints the permutation vector of length `n` to ``stdout``.
