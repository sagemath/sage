# distutils: libraries = flint
# distutils: depends = flint/fmpzi.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fmpzi_init(fmpzi_t x)

    void fmpzi_clear(fmpzi_t x)

    void fmpzi_swap(fmpzi_t x, fmpzi_t y)

    void fmpzi_zero(fmpzi_t x)

    void fmpzi_one(fmpzi_t x)

    void fmpzi_set(fmpzi_t res, const fmpzi_t x)

    void fmpzi_set_si_si(fmpzi_t res, slong a, slong b)

    void fmpzi_print(const fmpzi_t x)

    void fmpzi_randtest(fmpzi_t res, flint_rand_t state, mp_bitcnt_t bits)

    int fmpzi_equal(const fmpzi_t x, const fmpzi_t y)

    int fmpzi_is_zero(const fmpzi_t x)

    int fmpzi_is_one(const fmpzi_t x)

    int fmpzi_is_unit(const fmpzi_t x)

    slong fmpzi_canonical_unit_i_pow(const fmpzi_t x)

    void fmpzi_canonicalise_unit(fmpzi_t res, const fmpzi_t x)

    slong fmpzi_bits(const fmpzi_t x)

    void fmpzi_norm(fmpz_t res, const fmpzi_t x)

    void fmpzi_conj(fmpzi_t res, const fmpzi_t x)

    void fmpzi_neg(fmpzi_t res, const fmpzi_t x)

    void fmpzi_add(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)

    void fmpzi_sub(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)

    void fmpzi_sqr(fmpzi_t res, const fmpzi_t x)

    void fmpzi_mul(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)

    void fmpzi_pow_ui(fmpzi_t res, const fmpzi_t x, ulong exp)

    void fmpzi_divexact(fmpzi_t q, const fmpzi_t x, const fmpzi_t y)
    # Sets *q* to the quotient of *x* and *y*, assuming that the
    # division is exact.

    void fmpzi_divrem(fmpzi_t q, fmpzi_t r, const fmpzi_t x, const fmpzi_t y)
    # Computes a quotient and remainder satisfying
    # `x = q y + r` with `N(r) \le N(y)/2`, with a canonical
    # choice of remainder when breaking ties.

    void fmpzi_divrem_approx(fmpzi_t q, fmpzi_t r, const fmpzi_t x, const fmpzi_t y)
    # Computes a quotient and remainder satisfying
    # `x = q y + r` with `N(r) < N(y)`, with an implementation-defined,
    # non-canonical choice of remainder.

    slong fmpzi_remove_one_plus_i(fmpzi_t res, const fmpzi_t x)
    # Divide *x* exactly by the largest possible power `(1+i)^k`
    # and return the exponent *k*.

    void fmpzi_gcd_euclidean(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
    void fmpzi_gcd_euclidean_improved(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
    void fmpzi_gcd_binary(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
    void fmpzi_gcd_shortest(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
    void fmpzi_gcd(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
    # Computes the GCD of *x* and *y*. The result is in canonical
    # unit form.
    # The *euclidean* version is a straightforward implementation
    # of Euclid's algorithm. The *euclidean_improved* version is
    # optimized by performing approximate divisions.
    # The *binary* version uses a (1+i)-ary analog of the binary
    # GCD algorithm for integers [Wei2000]_.
    # The *shortest* version finds the GCD as the shortest vector in a lattice.
    # The default version chooses an algorithm automatically.
