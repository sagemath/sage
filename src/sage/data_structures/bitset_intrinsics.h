/* sage_setup: distribution = sagemath-categories
 */
#ifndef SAGE_BITSET_INTRINSICS
#define SAGE_BITSET_INTRINSICS

#include "gmp.h"
#include <stdio.h>

#if __POPCNT__ || __BMI__
    #include <immintrin.h>
#endif

#if __AVX2__
    #include <immintrin.h>
#elif __AVX__
    #include <immintrin.h>
    #include <emmintrin.h>
#elif __SSE4_1__
    #include <emmintrin.h>
    #include <smmintrin.h>
#endif


// This file contains functions of ``bitset_base.pxd``
// that can be optimized using intrinsics.

const mp_bitcnt_t LIMB_SIZE = sizeof(mp_limb_t);

#if __AVX__
    const mp_bitcnt_t ALIGNMENT = 32;
#elif __SSE4_1__
    const mp_bitcnt_t ALIGNMENT = 16;
#else
    const mp_bitcnt_t ALIGNMENT = sizeof(mp_limb_t);
#endif

// We assume that ``GMP_LIMB_BITS`` is a divisor of 64
#define LIMBS_PER_64 64/GMP_LIMB_BITS

static inline mp_bitcnt_t _set_non_zero(mp_limb_t* bits, mp_bitcnt_t* non_zero_chunks, mp_bitcnt_t limbs){
    /*
    Set the non zero positions of ``bits``.

    Return the number of non zero chunks.

    .. WARNING::

        It is assumed that ``bits`` is overaligned
        according to ``ALIGNMENT`` and its length is a multiple
        of ``ALIGNMENT``.
    */
    mp_bitcnt_t i;
    mp_bitcnt_t pos = 0;
#if __AVX__
    for(i = 0; i < limbs; i +=4*LIMBS_PER_64){
        __m256i A = _mm256_load_si256((const __m256i*)&bits[i]);
        if (!_mm256_testz_si256(A, A)){
            non_zero_chunks[pos] = i;
            pos++;
        }
    }
#elif __SSE4_1__
    for(i = 0; i < limbs; i +=2*LIMBS_PER_64){
        __m128i A = _mm_load_si128((const __m128i*)&bits[i]);
        if (!_mm_testz_si128(A, A)){
            non_zero_chunks[pos] = i;
            pos++;
        }
    }
#else
    for(i = 0; i < limbs; i++){
        if (bits[i]){
            non_zero_chunks[pos] = i;
            pos++;
        }
    }
#endif
    return pos;
}

/*
#############################################################################
# Bitset Comparison
#############################################################################
*/

enum cmpop_t {EQUAL, SUBSET, DISJOINT};

static inline int _bitset_isempty(mp_limb_t* bits, mp_bitcnt_t limbs){
    /*
    Test whether ``bits`` is empty.  Return ``True`` (i.e., ``1``) if the set is
    empty, ``False`` (i.e., ``0``) otherwise.
    */
    // First check lowest limb
    if (bits[0])
        return 0;
    if (limbs == 1)
        return 1;
    // Compare bits to itself shifted by 1 limb. If these compare equal,
    // all limbs must be 0.
    return mpn_cmp(bits+1, bits, limbs-1) == 0;
}

#if __AVX__
static inline int cmp_on_chunk(__m256i A, __m256i B, enum cmpop_t cmpop ){
    /*
    Test ``A cmpop B``,

    where ``cmpop`` is one of ``EQUAL``, ``SUBSET``, ``DISJOINT``.
    */
    if (cmpop == EQUAL)
        return _mm256_testc_si256(A, B) && _mm256_testc_si256(B, A);
    if (cmpop == SUBSET)
        return _mm256_testc_si256(B, A); // Need to be opposite order!
    // ``cmpop == DISJOINT``
    return _mm256_testz_si256(A, B);
}
#elif __SSE4_1__
static inline int cmp_on_chunk(__m128i A, __m128i B, enum cmpop_t cmpop){
    /*
    Test ``A cmpop B``,

    where ``cmpop`` is one of ``EQUAL``, ``SUBSET``, ``DISJOINT``.
    */
    if (cmpop == EQUAL)
        return _mm_testc_si128(A, B) && _mm_testc_si128(B, A);
    if (cmpop == SUBSET)
        return _mm_testc_si128(B, A); // Need to be opposite order!
    // ``cmpop == DISJOINT``
    return _mm_testz_si128(A, B);
}
#endif

static inline int cmp_on_limb(mp_limb_t A, mp_limb_t B, enum cmpop_t cmpop){
    /*
    Test ``A cmpop B``,

    where ``cmpop`` is one of ``EQUAL``, ``SUBSET``, ``DISJOINT``.
    */
    if (cmpop == EQUAL)
        return A == B;
    if (cmpop == SUBSET)
        return !(A & ~B);
    // ``cmpop == DISJOINT``
    return !(A & B);
}

static inline int cmp_with_mpn(mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs, enum cmpop_t cmpop){
    /*
    Test ``a cmpop b``,

    where ``cmpop`` is one of ``EQUAL``, ``SUBSET``, ``DISJOINT``.

    Use gmp if possible.
    */
    if (cmpop == EQUAL)
        return mpn_cmp(a, b, limbs) == 0;
    mp_bitcnt_t i;
    for(i = 0; i < limbs; i++){
        if (!cmp_on_limb(a[i], b[i], cmpop))
            return 0;
    }
    return 1;
}

static inline int _bitset_cmp(mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs, enum cmpop_t cmpop){
    /*
    Test ``a cmpop b``,

    where ``cmpop`` is one of ``EQUAL``, ``SUBSET``, ``DISJOINT``.
    */
#if __AVX__
    mp_bitcnt_t i;
    for(i = 0; i + (4*LIMBS_PER_64 - 1) < limbs; i += 4*LIMBS_PER_64){
        __m256i A = _mm256_loadu_si256((const __m256i*)&a[i]);
        __m256i B = _mm256_loadu_si256((const __m256i*)&b[i]);
        if (!cmp_on_chunk(A, B, cmpop))
            return 0;
    }
    for(; i < limbs; i++){
        if (!cmp_on_limb(a[i], b[i], cmpop))
                return 0;
    }
    return 1;
#elif __SSE4_1__
    mp_bitcnt_t i;
    for(i = 0; i + (2*LIMBS_PER_64 - 1) < limbs; i += 2*LIMBS_PER_64){
        __m128i A = _mm_loadu_si128((const __m128i*)&a[i]);
        __m128i B = _mm_loadu_si128((const __m128i*)&b[i]);
        if (!cmp_on_chunk(A, B, cmpop))
            return 0;
    }
    for(; i < limbs; i++){
        if (!cmp_on_limb(a[i], b[i], cmpop))
                return 0;
    }
    return 1;
#else
    return cmp_with_mpn(a, b, limbs, cmpop);
#endif
}

static inline int _sparse_bitset_cmp(mp_limb_t* a, mp_bitcnt_t* a_non_zero_chunks, mp_bitcnt_t a_n_non_zero_chunks, mp_limb_t* b, enum cmpop_t cmpop){
    /*
    Test ``a cmpop b``,

    where ``cmpop`` is one of ``SUBSET``, ``DISJOINT``.

    Only the non zero chunks of ``a`` are checked.

    .. WARNING::

        It is assumed that ``a`` and ``b`` are overaligned
        according to ``ALIGNMENT`` and their length are a multiple
        of ``ALIGNMENT``.
    */

    // For checking equality, it does not suffice to check only
    // the nonzero chunks of ``a``.
    assert(cmpop != EQUAL);
    mp_bitcnt_t i,j;
    for(j = 0; j < a_n_non_zero_chunks; j++){
        i = a_non_zero_chunks[j];
#if __AVX__
        __m256i A = _mm256_load_si256((const __m256i*)&a[i]);
        __m256i B = _mm256_load_si256((const __m256i*)&b[i]);
        if (!cmp_on_chunk(A, B, cmpop))
            return 0;

#elif __SSE4_1__
        __m128i A = _mm_load_si128((const __m128i*)&a[i]);
        __m128i B = _mm_load_si128((const __m128i*)&b[i]);
        if (!cmp_on_chunk(A, B, cmpop))
            return 0;
#else
        if (!cmp_on_limb(a[i], b[i], cmpop))
            return 0;
#endif
    }
    return 1;
}

/*
#############################################################################
# Bitset Searching
#############################################################################
*/

static inline long _bitset_first_in_limb_nonzero(mp_limb_t limb){
    /*
    Given a nonzero limb of a bitset, return the index of the first
    nonzero bit.
    */
#if (__BMI__) && (GMP_LIMB_BITS == 64) && (INTPTR_MAX == INT64_MAX)
    // If available we use intrinsics trailing zero count.
    // Also we check that ``mp_bitcnt_t`` is in fact 8 bytes long
    // and the architecture is 64-bit.
    return _tzcnt_u64(limb);
#else
    return mpn_scan1(&limb, 0);
#endif
}

static inline long _bitset_first_in_limb(mp_limb_t limb){
    /*
    Given a limb of a bitset, return the index of the first nonzero
    bit. If there are no bits set in the limb, return -1.
    */
    if (limb == 0)
        return -1;
    return _bitset_first_in_limb_nonzero(limb);
}

static inline long _bitset_len(mp_limb_t* bits, mp_bitcnt_t limbs){
    /*
    Calculate the number of items in the set (i.e., the number of nonzero bits).
    */
#if (__POPCNT__) && (GMP_LIMB_BITS == 64) && (INTPTR_MAX == INT64_MAX)
    // If available we use intrinsics popcount.
    // Also we check that ``mp_bitcnt_t`` is in fact 8 bytes long
    // and the architecture is 64-bit.
    mp_bitcnt_t i;
    uint64_t count = 0;
    for (i=0; i<limbs; i++){
        count += _mm_popcnt_u64(bits[i]);
    }
    return count;
#else
    return mpn_popcount(bits, limbs);
#endif
}

/*
#############################################################################
# Bitset Arithmetic
#############################################################################
*/

enum operation_t {AND, OR, ANDNOT, XOR};

#if __AVX2__
static inline __m256i operation_on_chunk(__m256i A, __m256i B, enum operation_t operation){
    /*
    Return ``A operation B``,

    where ``operation`` is one of ``AND``, ``OR``,
    ``ANDNOT``, ``XOR``.
    */
    if (operation == AND)
        return _mm256_and_si256(A, B);
    if (operation == OR)
        return _mm256_or_si256(A, B);
    if (operation == ANDNOT)
        return _mm256_andnot_si256(B, A);
    // ``operation == XOR``
    return _mm256_xor_si256(A, B);
}
#elif __SSE4_1__
static inline __m128i operation_on_chunk(__m128i A, __m128i B, enum operation_t operation){
    /*
    Return ``A operation B``,

    where ``operation`` is one of ``AND``, ``OR``,
    ``ANDNOT``, ``XOR``.
    */
    if (operation == AND)
        return _mm_and_si128(A, B);
    if (operation == OR)
        return _mm_or_si128(A, B);
    if (operation == ANDNOT)
        return _mm_andnot_si128(B, A);
    // ``operation == XOR``
    return _mm_xor_si128(A, B);
}
#endif

static inline mp_limb_t operation_on_limb(mp_limb_t A, mp_limb_t B, enum operation_t operation){
    /*
    Return ``A operation B``,

    where ``operation`` is one of ``AND``, ``OR``,
    ``ANDNOT``, ``XOR``.
    */
    if (operation == AND)
        return A & B;
    if (operation == OR)
        return A | B;
    if (operation == ANDNOT)
        return A & ~B;
    // ``operation == XOR``
    return A ^ B;
}

static inline void mpn_operation(mp_limb_t* dst, mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs, enum operation_t operation){
    /*
    Set ``dst`` to ``a operation b``,

    where ``operation`` is one of ``AND``, ``OR``,
    ``ANDNOT``, ``XOR``.

    Use gmp.
    */
    if (operation == AND)
        mpn_and_n(dst, a, b, limbs);
    if (operation == OR)
        mpn_ior_n(dst, a, b, limbs);
    if (operation == ANDNOT)
        mpn_andn_n(dst, a, b, limbs);
    if (operation == XOR)
        mpn_xor_n(dst, a, b, limbs);
}

static inline void _bitset_operation(mp_limb_t* dst, mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs, enum operation_t operation){
    /*
    Set ``dst`` to ``a operation b``,

    where ``operation`` is one of ``AND``, ``OR``,
    ``ANDNOT``, ``XOR``.
    */
#if __AVX2__
    mp_bitcnt_t i;
    for(i = 0; i + (4*LIMBS_PER_64 - 1) < limbs; i +=4*LIMBS_PER_64){
        __m256i A = _mm256_loadu_si256((const __m256i*)&a[i]);
        __m256i B = _mm256_loadu_si256((const __m256i*)&b[i]);
        __m256i D = operation_on_chunk(A, B, operation);
        _mm256_storeu_si256((__m256i*)&dst[i], D);
    }
    for(; i < limbs; i++){
        dst[i] = operation_on_limb(a[i], b[i], operation);
    }
#elif __SSE4_1__
    mp_bitcnt_t i;
    for(i = 0; i + (2*LIMBS_PER_64 - 1) < limbs; i +=2*LIMBS_PER_64){
        __m128i A = _mm_loadu_si128((const __m128i*)&a[i]);
        __m128i B = _mm_loadu_si128((const __m128i*)&b[i]);
        __m128i D = operation_on_chunk(A, B, operation);
        _mm_storeu_si128((__m128i*)&dst[i], D);
    }
    for(; i < limbs; i++){
        dst[i] = operation_on_limb(a[i], b[i], operation);
    }
#else
    mpn_operation(dst, a, b, limbs, operation);
#endif
}

/*
#############################################################################
# Bitset Arithmetic for sparse bitsets.
#############################################################################
*/

static inline mp_bitcnt_t _sparse_bitset_operation(mp_limb_t* dst, mp_bitcnt_t* dst_non_zero_chunks, mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs, enum operation_t operation){
    /*
    Set ``dst`` to ``a operation b``,

    where ``operation`` is one of ``AND``, ``OR``,
    ``ANDNOT``, ``XOR``.

    Set dst_non_zero_chunks according to the non zero chunks of dst.

    Return the number of non zero chunks.

    .. WARNING::

        It is assumed that ``a``, ``b``, ``dst`` are overaligned
        according to ``ALIGNMENT`` and their length are a multiple
        of ``ALIGNMENT``.
    */
    mp_bitcnt_t i;
    mp_bitcnt_t pos = 0;
#if __AVX2__
    for(i = 0; i < limbs; i +=4*LIMBS_PER_64){
        __m256i A = _mm256_load_si256((const __m256i*)&a[i]);
        __m256i B = _mm256_load_si256((const __m256i*)&b[i]);
        __m256i D = operation_on_chunk(A, B, operation);
        _mm256_store_si256((__m256i*)&dst[i], D);
        if (!cmp_on_chunk(D, D, DISJOINT)){
            dst_non_zero_chunks[pos] = i;
            pos++;
        }
    }
#elif __AVX__
    // In this case we can only perform the
    // operation on 128 bits each.
    // However, chunks are of size 256 bits,
    // because comparison can be performed on this
    // size.
    for(i = 0; i < limbs; i +=4*LIMBS_PER_64){
        __m128i A1 = _mm_load_si128((const __m128i*)&a[i]);
        __m128i B1 = _mm_load_si128((const __m128i*)&b[i]);
        __m128i D1 = operation_on_chunk(A1, B1, operation);
        _mm_store_si128((__m128i*)&dst[i], D1);
        __m128i A2 = _mm_load_si128((const __m128i*)&a[i+2*LIMBS_PER_64]);
        __m128i B2 = _mm_load_si128((const __m128i*)&b[i+2*LIMBS_PER_64]);
        __m128i D2 = operation_on_chunk(A2, B2, operation);
        _mm_store_si128((__m128i*)&dst[i+2*LIMBS_PER_64], D2);

        if (!_mm_testz_si128(D1, D1) || !_mm_testz_si128(D2, D2)){
            dst_non_zero_chunks[pos] = i;
            pos++;
        }
    }
#elif __SSE4_1__
    for(i = 0; i < limbs; i +=2*LIMBS_PER_64){
        __m128i A = _mm_load_si128((const __m128i*)&a[i]);
        __m128i B = _mm_load_si128((const __m128i*)&b[i]);
        __m128i D = operation_on_chunk(A, B, operation);
        _mm_store_si128((__m128i*)&dst[i], D);
        if (!cmp_on_chunk(D, D, DISJOINT)){
            dst_non_zero_chunks[pos] = i;
            pos++;
        }
    }
#else
    for(i = 0; i < limbs; i++){
        dst[i] = operation_on_limb(a[i], b[i], operation);
        if (dst[i]){
            dst_non_zero_chunks[pos] = i;
            pos++;
        }
    }
#endif
    return pos;
}

#endif
