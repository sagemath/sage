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

    slong * _perm_init(slong n)
    # Initialises the permutation for use.

    void _perm_clear(slong *vec)
    # Clears the permutation.

    void _perm_set(slong *res, const slong *vec, slong n)
    # Sets the permutation ``res`` to the same as the permutation ``vec``.

    void _perm_set_one(slong *vec, slong n)
    # Sets the permutation to the identity permutation.

    void _perm_inv(slong *res, const slong *vec, slong n)
    # Sets ``res`` to the inverse permutation of ``vec``.
    # Allows aliasing of ``res`` and ``vec``.

    void _perm_compose(slong *res, const slong *vec1, const slong *vec2, slong n)
    # Forms the composition `\pi_1 \circ \pi_2` of two permutations
    # `\pi_1` and `\pi_2`.  Here, `\pi_2` is applied first, that is,
    # `(\pi_1 \circ \pi_2)(i) = \pi_1(\pi_2(i))`.
    # Allows aliasing of ``res``, ``vec1`` and ``vec2``.

    int _perm_parity(const slong *vec, slong n)
    # Returns the parity of ``vec``, 0 if the permutation is even and 1 if
    # the permutation is odd.

    int _perm_randtest(slong *vec, slong n, flint_rand_t state)
    # Generates a random permutation vector of length `n` and returns
    # its parity, 0 or 1.
    # This function uses the Knuth shuffle algorithm to generate a uniformly
    # random permutation without retries.

    int _perm_print(const slong * vec, slong n)
    # Prints the permutation vector of length `n` to ``stdout``.
