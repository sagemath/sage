# distutils: libraries = flint
# distutils: depends = flint/nmod_poly_factor.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void nmod_poly_factor_init(nmod_poly_factor_t fac)
    # Initialises ``fac`` for use. An ``nmod_poly_factor_t``
    # represents a polynomial in factorised form as a product of
    # polynomials with associated exponents.

    void nmod_poly_factor_clear(nmod_poly_factor_t fac)
    # Frees all memory associated with ``fac``.

    void nmod_poly_factor_realloc(nmod_poly_factor_t fac, slong alloc)
    # Reallocates the factor structure to provide space for
    # precisely ``alloc`` factors.

    void nmod_poly_factor_fit_length(nmod_poly_factor_t fac, slong len)
    # Ensures that the factor structure has space for at
    # least ``len`` factors.  This function takes care
    # of the case of repeated calls by always at least
    # doubling the number of factors the structure can hold.

    void nmod_poly_factor_set(nmod_poly_factor_t res, const nmod_poly_factor_t fac)
    # Sets ``res`` to the same factorisation as ``fac``.

    void nmod_poly_factor_print(const nmod_poly_factor_t fac)
    # Prints the entries of ``fac`` to standard output.

    void nmod_poly_factor_insert(nmod_poly_factor_t fac, const nmod_poly_t poly, slong exp)
    # Inserts the factor ``poly`` with multiplicity ``exp`` into
    # the factorisation ``fac``.
    # If ``fac`` already contains ``poly``, then ``exp`` simply
    # gets added to the exponent of the existing entry.

    void nmod_poly_factor_concat(nmod_poly_factor_t res, const nmod_poly_factor_t fac)
    # Concatenates two factorisations.
    # This is equivalent to calling :func:`nmod_poly_factor_insert`
    # repeatedly with the individual factors of ``fac``.
    # Does not support aliasing between ``res`` and ``fac``.

    void nmod_poly_factor_pow(nmod_poly_factor_t fac, slong exp)
    # Raises ``fac`` to the power ``exp``.

    ulong nmod_poly_remove(nmod_poly_t f, const nmod_poly_t p)
    # Removes the highest possible power of ``p`` from ``f`` and
    # returns the exponent.

    bint nmod_poly_is_irreducible(const nmod_poly_t f)
    # Returns 1 if the polynomial ``f`` is irreducible, otherwise returns 0.

    bint nmod_poly_is_irreducible_ddf(const nmod_poly_t f)
    # Returns 1 if the polynomial ``f`` is irreducible, otherwise returns 0.
    # Uses fast distinct-degree factorisation.

    bint nmod_poly_is_irreducible_rabin(const nmod_poly_t f)
    # Returns 1 if the polynomial ``f`` is irreducible, otherwise returns 0.
    # Uses Rabin irreducibility test.

    bint _nmod_poly_is_squarefree(mp_srcptr f, slong len, nmod_t mod)
    # Returns 1 if ``(f, len)`` is squarefree, and 0 otherwise. As a
    # special case, the zero polynomial is not considered squarefree.
    # There are no restrictions on the length.

    bint nmod_poly_is_squarefree(const nmod_poly_t f)
    # Returns 1 if ``f`` is squarefree, and 0 otherwise. As a special
    # case, the zero polynomial is not considered squarefree.

    void nmod_poly_factor_squarefree(nmod_poly_factor_t res, const nmod_poly_t f)
    # Sets ``res`` to a square-free factorization of ``f``.

    bint nmod_poly_factor_equal_deg_prob(nmod_poly_t factor, flint_rand_t state, const nmod_poly_t pol, slong d)
    # Probabilistic equal degree factorisation of ``pol`` into
    # irreducible factors of degree ``d``. If it passes, a factor is
    # placed in factor and 1 is returned, otherwise 0 is returned and
    # the value of factor is undetermined.
    # Requires that ``pol`` be monic, non-constant and squarefree.

    void nmod_poly_factor_equal_deg(nmod_poly_factor_t factors, const nmod_poly_t pol, slong d)
    # Assuming ``pol`` is a product of irreducible factors all of
    # degree ``d``, finds all those factors and places them in factors.
    # Requires that ``pol`` be monic, non-constant and squarefree.

    void nmod_poly_factor_distinct_deg(nmod_poly_factor_t res, const nmod_poly_t poly, slong * const *degs)
    # Factorises a monic non-constant squarefree polynomial ``poly``
    # of degree n into factors `f[d]` such that for `1 \leq d \leq n`
    # `f[d]` is the product of the monic irreducible factors of ``poly``
    # of degree `d`. Factors `f[d]` are stored in ``res``, and the degree `d`
    # of the irreducible factors is stored in ``degs`` in the same order
    # as the factors.
    # Requires that ``degs`` has enough space for ``(n/2)+1 * sizeof(slong)``.

    void nmod_poly_factor_distinct_deg_threaded(nmod_poly_factor_t res, const nmod_poly_t poly, slong * const *degs)
    # Multithreaded version of :func:`nmod_poly_factor_distinct_deg`.

    void nmod_poly_factor_cantor_zassenhaus(nmod_poly_factor_t res, const nmod_poly_t f)
    # Factorises a non-constant polynomial ``f`` into monic irreducible
    # factors using the Cantor-Zassenhaus algorithm.

    void nmod_poly_factor_berlekamp(nmod_poly_factor_t res, const nmod_poly_t f)
    # Factorises a non-constant, squarefree polynomial ``f`` into monic
    # irreducible factors using the Berlekamp algorithm.

    void nmod_poly_factor_kaltofen_shoup(nmod_poly_factor_t res, const nmod_poly_t poly)
    # Factorises a non-constant polynomial ``f`` into monic irreducible
    # factors using the fast version of Cantor-Zassenhaus algorithm proposed by
    # Kaltofen and Shoup (1998). More precisely this algorithm uses a
    # “baby step/giant step” strategy for the distinct-degree factorization
    # step. If :func:`flint_get_num_threads` is greater than one
    # :func:`nmod_poly_factor_distinct_deg_threaded` is used.

    mp_limb_t nmod_poly_factor_with_berlekamp(nmod_poly_factor_t res, const nmod_poly_t f)
    # Factorises a general polynomial ``f`` into monic irreducible factors
    # and returns the leading coefficient of ``f``, or 0 if ``f``
    # is the zero polynomial.
    # This function first checks for small special cases, deflates ``f``
    # if it is of the form `p(x^m)` for some `m > 1`, then performs a
    # square-free factorisation, and finally runs Berlekamp on all the
    # individual square-free factors.

    mp_limb_t nmod_poly_factor_with_cantor_zassenhaus(nmod_poly_factor_t res, const nmod_poly_t f)
    # Factorises a general polynomial ``f`` into monic irreducible factors
    # and returns the leading coefficient of ``f``, or 0 if ``f``
    # is the zero polynomial.
    # This function first checks for small special cases, deflates ``f``
    # if it is of the form `p(x^m)` for some `m > 1`, then performs a
    # square-free factorisation, and finally runs Cantor-Zassenhaus on all the
    # individual square-free factors.

    mp_limb_t nmod_poly_factor_with_kaltofen_shoup(nmod_poly_factor_t res, const nmod_poly_t f)
    # Factorises a general polynomial ``f`` into monic irreducible factors
    # and returns the leading coefficient of ``f``, or 0 if ``f``
    # is the zero polynomial.
    # This function first checks for small special cases, deflates ``f``
    # if it is of the form `p(x^m)` for some `m > 1`, then performs a
    # square-free factorisation, and finally runs Kaltofen-Shoup on all the
    # individual square-free factors.

    mp_limb_t nmod_poly_factor(nmod_poly_factor_t res, const nmod_poly_t f)
    # Factorises a general polynomial ``f`` into monic irreducible factors
    # and returns the leading coefficient of ``f``, or 0 if ``f``
    # is the zero polynomial.
    # This function first checks for small special cases, deflates ``f``
    # if it is of the form `p(x^m)` for some `m > 1`, then performs a
    # square-free factorisation, and finally runs either Cantor-Zassenhaus
    # or Berlekamp on all the individual square-free factors.
    # Currently Cantor-Zassenhaus is used by default unless the modulus is 2, in
    # which case Berlekamp is used.

    void _nmod_poly_interval_poly_worker(void* arg_ptr)
    # Worker function to compute interval polynomials in distinct degree
    # factorisation. Input/output is stored in
    # ``nmod_poly_interval_poly_arg_t``.
