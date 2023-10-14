# distutils: libraries = flint
# distutils: depends = flint/mpoly.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void mpoly_ctx_init(mpoly_ctx_t ctx, slong nvars, const ordering_t ord)
    # Initialize a context for specified number of variables and ordering.

    void mpoly_ctx_clear(mpoly_ctx_t mctx)
    # Clean up any space used by a context object.

    ordering_t mpoly_ordering_randtest(flint_rand_t state)
    # Return a random ordering. The possibilities are ``ORD_LEX``,
    # ``ORD_DEGLEX`` and ``ORD_DEGREVLEX``.

    void mpoly_ctx_init_rand(mpoly_ctx_t mctx, flint_rand_t state, slong max_nvars)
    # Initialize a context with a random choice for the ordering.

    int mpoly_ordering_isdeg(const mpoly_ctx_t ctx)
    # Return 1 if the ordering of the given context is a degree ordering (deglex or degrevlex).

    int mpoly_ordering_isrev(const mpoly_ctx_t cth)
    # Return 1 if the ordering of the given context is a reverse ordering (currently only
    # degrevlex).

    void mpoly_ordering_print(ordering_t ord)
    # Print a string (either "lex", "deglex" or "degrevlex") to standard
    # output, corresponding to the given ordering.

    void mpoly_monomial_add(ulong * exp_ptr, const ulong * exp2, const ulong * exp3, slong N)
    # Set ``(exp_ptr, N)`` to the sum of the monomials ``(exp2, N)`` and
    # ``(exp3, N)``, assuming ``bits <= FLINT_BITS``

    void mpoly_monomial_add_mp(ulong * exp_ptr, const ulong * exp2, const ulong * exp3, slong N)
    # Set ``(exp_ptr, N)`` to the sum of the monomials ``(exp2, N)`` and
    # ``(exp3, N)``.

    void mpoly_monomial_sub(ulong * exp_ptr, const ulong * exp2, const ulong * exp3, slong N)
    # Set ``(exp_ptr, N)`` to the difference of the monomials ``(exp2, N)`` and ``(exp3, N)``, assuming ``bits <= FLINT_BITS``

    void mpoly_monomial_sub_mp(ulong * exp_ptr, const ulong * exp2, const ulong * exp3, slong N)
    # Set ``(exp_ptr, N)`` to the difference of the monomials ``(exp2, N)`` and ``(exp3, N)``.

    int mpoly_monomial_overflows(ulong * exp2, slong N, ulong mask)
    # Return true if any of the fields of the given monomial ``(exp2, N)`` has
    # overflowed (or is negative). The ``mask`` is a word with the high bit of
    # each field set to 1. In other words, the function returns 1 if any word of
    # ``exp2`` has any of the nonzero bits in ``mask`` set. Assumes that
    # ``bits <= FLINT_BITS``.

    int mpoly_monomial_overflows_mp(ulong * exp_ptr, slong N, flint_bitcnt_t bits)
    # Return true if any of the fields of the given monomial ``(exp_ptr, N)``
    # has overflowed. Assumes that ``bits >= FLINT_BITS``.

    int mpoly_monomial_overflows1(ulong exp, ulong mask)
    # As per ``mpoly_monomial_overflows`` with ``N = 1``.

    void mpoly_monomial_set(ulong * exp2, const ulong * exp3, slong N)
    # Set the monomial ``(exp2, N)`` to ``(exp3, N)``.

    void mpoly_monomial_swap(ulong * exp2, ulong * exp3, slong N)
    # Swap the words in ``(exp2, N)`` and ``(exp3, N)``.

    void mpoly_monomial_mul_ui(ulong * exp2, const ulong * exp3, slong N, ulong c)
    # Set the words of ``(exp2, N)`` to the words of ``(exp3, N)``
    # multiplied by ``c``.

    bint mpoly_monomial_is_zero(const ulong * exp, slong N)
    # Return 1 if ``(exp, N)`` is zero.

    bint mpoly_monomial_equal(const ulong * exp2, const ulong * exp3, slong N)
    # Return 1 if the monomials ``(exp2, N)`` and ``(exp3, N)`` are equal.

    void mpoly_get_cmpmask(ulong * cmpmask, slong N, ulong bits, const mpoly_ctx_t mctx)
    # Get the mask ``(cmpmask, N)`` for comparisons.
    # ``bits`` should be set to the number of bits in the exponents
    # to be compared. Any function that compares monomials should use this
    # comparison mask.

    bint mpoly_monomial_lt(const ulong * exp2, const ulong * exp3, slong N, const ulong * cmpmask)
    # Return 1 if ``(exp2, N)`` is less than ``(exp3, N)``.

    bint mpoly_monomial_gt(const ulong * exp2, const ulong * exp3, slong N, const ulong * cmpmask)
    # Return 1 if ``(exp2, N)`` is greater than ``(exp3, N)``.

    int mpoly_monomial_cmp(const ulong * exp2, const ulong * exp3, slong N, const ulong * cmpmask)
    # Return `1` if ``(exp2, N)`` is greater than, `0` if it is equal to and
    # `-1` if it is less than ``(exp3, N)``.

    int mpoly_monomial_divides(ulong * exp_ptr, const ulong * exp2, const ulong * exp3, slong N, ulong mask)
    # Return 1 if the monomial ``(exp3, N)`` divides ``(exp2, N)``. If so
    # set ``(exp_ptr, N)`` to the quotient monomial. The ``mask`` is a word
    # with the high bit of each bit field set to 1. Assumes that
    # ``bits <= FLINT_BITS``.

    int mpoly_monomial_divides_mp(ulong * exp_ptr, const ulong * exp2, const ulong * exp3, slong N, flint_bitcnt_t bits)
    # Return 1 if the monomial ``(exp3, N)`` divides ``(exp2, N)``. If so
    # set ``(exp_ptr, N)`` to the quotient monomial. Assumes that
    # ``bits >= FLINT_BITS``.

    int mpoly_monomial_divides1(ulong * exp_ptr, const ulong exp2, const ulong exp3, ulong mask)
    # As per ``mpoly_monomial_divides`` with ``N = 1``.

    int mpoly_monomial_divides_tight(slong e1, slong e2, slong * prods, slong num)
    # Return 1 if the monomial ``e2`` divides the monomial ``e1``, where
    # the monomials are stored using factorial representation. The array
    # ``(prods, num)`` should consist of `1`, `b_1, b_1\times b_2, \ldots`,
    # where the `b_i` are the bases of the factorial number representation.

    flint_bitcnt_t mpoly_exp_bits_required_ui(const ulong * user_exp, const mpoly_ctx_t mctx)
    # Returns the number of bits required to store ``user_exp`` in packed
    # format. The returned number of bits includes space for a zeroed signed bit.

    flint_bitcnt_t mpoly_exp_bits_required_ffmpz(const fmpz * user_exp, const mpoly_ctx_t mctx)
    # Returns the number of bits required to store ``user_exp`` in packed
    # format. The returned number of bits includes space for a zeroed signed bit.

    flint_bitcnt_t mpoly_exp_bits_required_pfmpz(fmpz * const * user_exp, const mpoly_ctx_t mctx)
    # Returns the number of bits required to store ``user_exp`` in packed
    # format. The returned number of bits includes space for a zeroed signed bit.

    void mpoly_max_fields_ui_sp(ulong * max_fields, const ulong * poly_exps, slong len, ulong bits, const mpoly_ctx_t mctx)
    # Compute the field-wise maximum of packed exponents from ``poly_exps``
    # of length ``len`` and unpack the result into ``max_fields``.
    # The maximums are assumed to fit a ulong.

    void mpoly_max_fields_fmpz(fmpz * max_fields, const ulong * poly_exps, slong len, ulong bits, const mpoly_ctx_t mctx)
    # Compute the field-wise maximum of packed exponents from ``poly_exps``
    # of length ``len`` and unpack the result into ``max_fields``.

    void mpoly_max_degrees_tight(slong * max_exp, ulong * exps, slong len, slong * prods, slong num)
    # Return an array of ``num`` integers corresponding to the maximum degrees
    # of the exponents in the array of exponent vectors ``(exps, len)``,
    # assuming that the exponent are packed in a factorial representation. The
    # array ``(prods, num)`` should consist of `1`, `b_1`,
    # `b_1\times b_2, \ldots`, where the `b_i` are the bases of the factorial
    # number representation. The results are stored in the array ``max_exp``,
    # with the entry corresponding to the most significant base of the factorial
    # representation first in the array.

    int mpoly_monomial_exists(slong * index, const ulong * poly_exps, const ulong * exp, slong len, slong N, const ulong * cmpmask)
    # Returns true if the given exponent vector ``exp`` exists in the array of
    # exponent vectors ``(poly_exps, len)``, otherwise, returns false. If the
    # exponent vector is found, its index into the array of exponent vectors is
    # returned. Otherwise, ``index`` is set to the index where this exponent
    # could be inserted to preserve the ordering. The index can be in the range
    # ``[0, len]``.

    void mpoly_search_monomials(slong ** e_ind, ulong * e, slong * e_score, slong * t1, slong * t2, slong *t3, slong lower, slong upper, const ulong * a, slong a_len, const ulong * b, slong b_len, slong N, const ulong * cmpmask)
    # Given packed exponent vectors ``a`` and ``b``, compute a packed
    # exponent ``e`` such that the number of monomials in the cross product
    # ``a`` X ``b`` that are less than or equal to ``e`` is between
    # ``lower`` and ``upper``. This number is stored in ``e_store``. If
    # no such monomial exists, one is chosen so that the number of monomials is as
    # close as possible. This function assumes that ``1`` is the smallest
    # monomial and needs three arrays ``t1``, ``t2``, and ``t3`` of the
    # size as ``a`` for workspace. The parameter ``e_ind`` is set to one
    # of ``t1``, ``t2``, and ``t3`` and gives the locations of the
    # monomials in ``a`` X ``b``.

    int mpoly_term_exp_fits_ui(ulong * exps, ulong bits, slong n, const mpoly_ctx_t mctx)
    # Return whether every entry of the exponent vector of index `n` in
    # ``exps`` fits into a ``ulong``.

    int mpoly_term_exp_fits_si(ulong * exps, ulong bits, slong n, const mpoly_ctx_t mctx)
    # Return whether every entry of the exponent vector of index `n` in
    # ``exps`` fits into a ``slong``.

    void mpoly_get_monomial_ui(ulong * exps, const ulong * poly_exps, ulong bits, const mpoly_ctx_t mctx)
    # Convert the packed exponent ``poly_exps`` of bit count ``bits`` to a
    # monomial from the user's perspective. The exponents are assumed to fit
    # a ulong.

    void mpoly_get_monomial_ffmpz(fmpz * exps, const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
    # Convert the packed exponent ``poly_exps`` of bit count ``bits`` to a
    # monomial from the user's perspective.

    void mpoly_get_monomial_pfmpz(fmpz ** exps, const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
    # Convert the packed exponent ``poly_exps`` of bit count ``bits`` to a
    # monomial from the user's perspective.

    void mpoly_set_monomial_ui(ulong * exp1, const ulong * exp2, ulong bits, const mpoly_ctx_t mctx)
    # Convert the user monomial ``exp2`` to packed format using ``bits``.

    void mpoly_set_monomial_ffmpz(ulong * exp1, const fmpz * exp2, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
    # Convert the user monomial ``exp2`` to packed format using ``bits``.

    void mpoly_set_monomial_pfmpz(ulong * exp1, fmpz * const * exp2, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
    # Convert the user monomial ``exp2`` to packed format using ``bits``.

    void mpoly_pack_vec_ui(ulong * exp1, const ulong * exp2, ulong bits, slong nfields, slong len)
    # Packs a vector ``exp2`` into \{exp1} using a bit count of ``bits``.
    # No checking is done to ensure that the vector actually fits
    # into ``bits`` bits. The number of fields in each vector is
    # ``nfields`` and the total number of vectors to unpack is ``len``.

    void mpoly_pack_vec_fmpz(ulong * exp1, const fmpz * exp2, flint_bitcnt_t bits, slong nfields, slong len)
    # Packs a vector ``exp2`` into \{exp1} using a bit count of ``bits``.
    # No checking is done to ensure that the vector actually fits
    # into ``bits`` bits. The number of fields in each vector is
    # ``nfields`` and the total number of vectors to unpack is ``len``.

    void mpoly_unpack_vec_ui(ulong * exp1, const ulong * exp2, ulong bits, slong nfields, slong len)
    # Unpacks vector ``exp2`` of bit count ``bits`` into ``exp1``.
    # The number of fields in each vector is
    # ``nfields`` and the total number of vectors to unpack is ``len``.

    void mpoly_unpack_vec_fmpz(fmpz * exp1, const ulong * exp2, flint_bitcnt_t bits, slong nfields, slong len)
    # Unpacks vector ``exp2`` of bit count ``bits`` into ``exp1``.
    # The number of fields in each vector is
    # ``nfields`` and the total number of vectors to unpack is ``len``.

    int mpoly_repack_monomials(ulong * exps1, ulong bits1, const ulong * exps2, ulong bits2, slong len, const mpoly_ctx_t mctx)
    # Convert an array of length ``len`` of exponents ``exps2`` packed
    # using bits ``bits2`` into an array ``exps1`` using bits ``bits1``.
    # No checking is done to ensure that the result fits into bits ``bits1``.

    void mpoly_pack_monomials_tight(ulong * exp1, const ulong * exp2, slong len, const slong * mults, slong num, slong bits)
    # Given an array of possibly packed exponent vectors ``exp2`` of length
    # ``len``, where each field of each exponent vector is packed into the
    # given number of bits, return the corresponding array of monomial vectors
    # packed using a factorial numbering scheme. The "bases" for the factorial
    # numbering scheme are given as an array of integers ``mults``, the first
    # entry of which corresponds to the field of least significance in each
    # input exponent vector. Obviously the maximum exponent to be packed must be
    # less than the corresponding base in ``mults``.
    # The number of multipliers is given by ``num``. The code only considers
    # least significant ``num`` fields of each exponent vectors and ignores
    # the rest. The number of ignored fields should be passed in ``extras``.

    void mpoly_unpack_monomials_tight(ulong * e1, ulong * e2, slong len, slong * mults, slong num, slong bits)
    # Given an array of exponent vectors ``e2`` of length ``len`` packed
    # using a factorial numbering scheme, unpack the monomials into an array
    # ``e1`` of exponent vectors in standard packed format, where each field
    # has the given number of bits. The "bases" for the factorial
    # numbering scheme are given as an array of integers ``mults``, the first
    # entry of which corresponds to the field of least significance in each
    # exponent vector.

    void mpoly_main_variable_terms1(slong * i1, slong * n1, const ulong * exp1, slong l1, slong len1, slong k, slong num, slong bits)
    # Given an array of exponent vectors ``(exp1, len1)``, each exponent
    # vector taking one word of space, with each exponent being packed into the
    # given number of bits, compute ``l1`` starting offsets ``i1`` and
    # lengths ``n1`` (which may be zero) to break the exponents into chunks.
    # Each chunk consists of exponents have the same degree in the main variable.
    # The index of the main variable is given by `k`. The variables are indexed
    # from the variable of least significance, starting from `0`. The value
    # ``l1`` should be the degree in the main variable, plus one.

    int _mpoly_heap_insert(mpoly_heap_s * heap, ulong * exp, void * x, slong * next_loc, slong * heap_len, slong N, const ulong * cmpmask)
    # Given a heap, insert a new node `x` corresponding to the given exponent
    # into the heap. Heap elements are ordered by the exponent ``(exp, N)``,
    # with the largest element at the head of the heap. A pointer to the current
    # heap length must be passed in via ``heap_len``. This will be updated by
    # the function. Note that the index 0 position in the heap is not used, so
    # the length is always one greater than the number of elements.

    void _mpoly_heap_insert1(mpoly_heap1_s * heap, ulong exp, void * x, slong * next_loc, slong * heap_len, ulong maskhi)
    # As per ``_mpoly_heap_insert`` except that ``N = 1``, and
    # ``maskhi = cmpmask[0]``.

    void * _mpoly_heap_pop(mpoly_heap_s * heap, slong * heap_len, slong N, const ulong * cmpmask)
    # Pop the head of the heap. It is cast to a ``void *``. A pointer to the
    # current heap length must be passed in via ``heap_len``. This will be
    # updated by the function. Note that the index 0 position in the heap is not
    # used, so the length is always one greater than the number of elements. The
    # ``maskhi`` and ``masklo`` values are zero except for degrevlex
    # ordering, where they are as per the monomial comparison operations above.

    void * _mpoly_heap_pop1(mpoly_heap1_s * heap, slong * heap_len, ulong maskhi)
    # As per ``_mpoly_heap_pop1`` except that ``N = 1``, and
    # ``maskhi = cmpmask[0]``.
