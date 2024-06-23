# sage_setup: distribution = sagemath-flint
# sage.doctest: needs sage.graphs
"""
Boilerplate functions for a cython implementation of elements of path algebras

AUTHORS:

- Simon King (2015-08)
"""

# ****************************************************************************
#     Copyright (C) 2015 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport check_malloc, check_allocarray, sig_free
from cysignals.signals cimport sig_check, sig_on, sig_off

from cpython.ref cimport *
from cython.operator cimport predecrement as predec, postincrement as postinc

from sage.data_structures.bitset_base cimport *
from sage.structure.richcmp cimport richcmp_not_equal, rich_to_bool
from sage.libs.gmp.mpn cimport mpn_cmp
from libc.stdlib cimport free

cdef extern from *:  # Defined by Cython
    int unlikely(int) nogil
    int likely(int) nogil

########################################
#
#  Allocation and Deallocation of monomials
#
# Monomials are expensive, hence, copying will just be done by increasing a
# reference counter.

# Create a monomial by copying the given bounded integer sequence
cdef bint mon_create(path_mon_t out, biseq_t Mon, long Pos, mp_size_t L_len, mp_size_t S_len) except -1:
    biseq_init_copy(out.path, Mon)
    out.pos = Pos
    out.l_len = L_len
    out.s_len = S_len

# The following is only used in the free-list for terms.
# It changes an existing monomial in-place (which should NEVER
# be done on a monomial that is in use), re-allocating memory
# and filling it with a copy of the given bounded integer sequence.
cdef bint mon_realloc(path_mon_t out, biseq_t Mon, long Pos, mp_size_t L_len, mp_size_t S_len) except -1:
    biseq_dealloc(out.path)
    sig_check()
    biseq_init_copy(out.path, Mon)
    out.pos = Pos
    out.l_len = L_len
    out.s_len = S_len

# Create a monomial without copying the given bounded integer sequence
cdef bint mon_create_keep(path_mon_t out, biseq_t Mon, long Pos, mp_size_t L_len, mp_size_t S_len) except -1:
    out.path[0] = Mon[0]
    out.pos = Pos
    out.l_len = L_len
    out.s_len = S_len

# The following is only used in the free-list for terms.
# It changes an existing monomial in-place (which should NEVER
# be done on a monomial that is in use), re-allocating memory
# and filling it with the given bounded integer sequence (not a copy).
cdef bint mon_realloc_keep(path_mon_t out, biseq_t Mon, long Pos, mp_size_t L_len, mp_size_t S_len) noexcept:
    biseq_dealloc(out.path)
    out.path[0] = Mon[0]
    out.pos = Pos
    out.l_len = L_len
    out.s_len = S_len
    return True

cdef inline bint mon_copy(path_mon_t out, path_mon_t M) except -1:
    out.pos = M.pos
    out.l_len = M.l_len
    out.s_len = M.s_len
    biseq_init_copy(out.path, M.path)

# Deallocate the monomial, which means to decrease the reference count,
# or to actually deallocate the data if there is no reference left.
cdef inline void mon_free(path_mon_t M) noexcept:
    biseq_dealloc(M.path)

# Linearisation
cdef inline tuple mon_pickle(path_mon_t M):
    return (bitset_pickle(M.path.data) if M.path.length>0 else (),
            M.path.itembitsize, M.path.length, M.pos, M.l_len, M.s_len)

# De-linearisation
cdef bint mon_unpickle(path_mon_t out, tuple data) except -1:
    cdef tuple bitset_data
    cdef mp_bitcnt_t itembitsize
    cdef mp_size_t length
    cdef long Pos
    cdef mp_size_t L_len
    cdef mp_size_t S_len
    bitset_data, itembitsize, length, Pos, L_len, S_len = data
    out.path.itembitsize = itembitsize
    out.path.mask_item = limb_lower_bits_up(itembitsize)
    out.path.length = length

    # bitset_unpickle assumes that out.path.data is initialised.
    bitset_init(out.path.data, GMP_LIMB_BITS)
    if bitset_data:
        sig_on()
        bitset_unpickle(out.path.data, bitset_data)
        sig_off()
    out.pos = Pos
    out.l_len = L_len
    out.s_len = S_len


########################################
#
#  Monomial orders---we only use degree orders

# Negative degree reverse lexicographic ordering
cdef int negdegrevlex(path_mon_t M1, path_mon_t M2) except -2:
    # a*s_i*b<c*s_j*d <=>
    # 1. deg(a*b) > deg(c*d), otherwise
    # 2. deg(a) > deg(c) (note that one of them may be -1), otherwise
    # 3. deg(s_i) < deg(s_j), otherwise
    # 4. a*s_i*b <_revlex c*s_j*d, otherwise
    # 5. i<j
    cdef mp_size_t l1 = M1.path.length + M2.s_len  # sic!
    cdef mp_size_t l2 = M2.path.length + M1.s_len
    if l1 != l2:
        if l2 < l1:
            return -1
        return 1
    if M2.l_len != M1.l_len:
        if M2.l_len < M1.l_len:
            return -1
        return 1
    if M1.s_len != M2.s_len:
        if M1.s_len < M2.s_len:
            return -1
        return 1
    # mpn_cmp does comparison of long integers. If the two long integers have
    # the same number of digits (this is the case her), it is the same as
    # lexicographic comparison of the numbers. The highest digit corresponds
    # to the right-most item in the path. Hence, it becomes
    # reverse-lexicographic order.
    sig_on()
    cdef int c = mpn_cmp(M1.path.data.bits, M2.path.data.bits, M1.path.data.limbs)
    sig_off()
    if c!=0:
        return c
    if M1.pos != M2.pos:
        if M1.pos < M2.pos:
            return -1
        return 1
    return 0

# Degree reverse lexicographic ordering
cdef int degrevlex(path_mon_t M1, path_mon_t M2) except -2:
    # a*s_i*b<c*s_j*d <=>
    # 1. deg(a*b) < deg(c*d), otherwise
    # 2. deg(a) < deg(c) (note that one of them may be -1), otherwise
    # 3. deg(s_i) > deg(s_j), otherwise
    # 4. a*s_i*b <_revlex c*s_j*d, otherwise
    # 5. i<j
    cdef mp_size_t l1 = M1.path.length + M2.s_len  # sic!
    cdef mp_size_t l2 = M2.path.length + M1.s_len
    if l2 != l1:
        if l2 < l1:
            return 1
        return -1
    if M2.l_len != M1.l_len:
        if M2.l_len < M1.l_len:
            return 1
        return -1
    if M1.s_len != M2.s_len:
        if M1.s_len < M2.s_len:
            return 1
        return -1
    # mpn_cmp does comparison of long integers. If the two long integers have
    # the same number of digits (this is the case her), it is the same as
    # lexicographic comparison of the numbers. The highest digit corresponds
    # to the right-most item in the path. Hence, it becomes
    # reverse-lexicographic order.
    sig_on()
    cdef int c = mpn_cmp(M1.path.data.bits, M2.path.data.bits, M1.path.data.limbs)
    sig_off()
    if c!=0:
        return c
    if M1.pos != M2.pos:
        if M1.pos < M2.pos:
            return -1
        return 1
    return 0

# Negative degree lexicographic ordering
cdef int negdeglex(path_mon_t M1, path_mon_t M2) except -2:
    # a*s_i*b<c*s_j*d <=>
    # 1. deg(a*b) > deg(c*d), otherwise
    # 2. deg(a) > deg(c) (note that one of them may be -1), otherwise
    # 3. deg(s_i) < deg(s_j), otherwise
    # 4. a*s_i*b <_lex c*s_j*d, otherwise
    # 5. i<j
    cdef mp_size_t l1 = M1.path.length + M2.s_len  # sic!
    cdef mp_size_t l2 = M2.path.length + M1.s_len
    cdef size_t item1, item2
    if l2 != l1:
        if l2 < l1:
            return -1
        return 1
    if M2.l_len != M1.l_len:
        if M2.l_len < M1.l_len:
            return -1
        return 1
    if M1.s_len != M2.s_len:
        if M1.s_len < M2.s_len:
            return -1
        return 1
    for index in range(M1.path.length):
        item1 = biseq_getitem(M1.path, index)
        item2 = biseq_getitem(M2.path, index)
        sig_check()
        if item1 != item2:
            if item1 < item2:
                return -1
            return 1
    if M1.pos != M2.pos:
        if M1.pos < M2.pos:
            return -1
        return 1
    return 0

# Degree lexicographic ordering
cdef int deglex(path_mon_t M1, path_mon_t M2) except -2:
    # a*s_i*b<c*s_j*d <=>
    # 1. deg(a*b) < deg(c*d), otherwise
    # 2. deg(a) < deg(c) (note that one of them may be -1), otherwise
    # 3. deg(s_i) > deg(s_j), otherwise
    # 4. a*s_i*b <_lex c*s_j*d, otherwise
    # 5. i<j
    cdef mp_size_t l1 = M1.path.length + M2.s_len  # sic!
    cdef mp_size_t l2 = M2.path.length + M1.s_len
    cdef size_t item1, item2
    if l2 != l1:
        if l2 < l1:
            return 1
        return -1
    if M2.l_len != M1.l_len:
        if M2.l_len < M1.l_len:
            return 1
        return -1
    if M1.s_len != M2.s_len:
        if M1.s_len < M2.s_len:
            return 1
        return -1
    for index in range(M1.path.length):
        item1 = biseq_getitem(M1.path, index)
        item2 = biseq_getitem(M2.path, index)
        sig_check()
        if item1 != item2:
            if item1 < item2:
                return -1
            return 1
    if M1.pos != M2.pos:
        if M1.pos < M2.pos:
            return -1
        return 1
    return 0

########################################
#
#  Allocation and Deallocation of terms
###########################
# We use a freelist for terms

cdef struct freelist_t:
    path_term_t **pool
    size_t used
cdef size_t poolsize = 5000   # The freelist contains at most that many terms.

cdef freelist_t *freelist = <freelist_t*>check_malloc(sizeof(freelist_t))
freelist.used = 0
freelist.pool = <path_term_t**>check_allocarray(poolsize, sizeof(path_term_t*))

# Deallocate the term, and return the pointer .nxt, without using kill list
cdef inline path_term_t *term_free_force(path_term_t *T) noexcept:
    mon_free(T.mon)
    cdef path_term_t *out = T.nxt
    sig_free(T)
    return out

cdef class _FreeListProtector:
    """
    The purpose of this class is to deallocate our freelist
    of path algebra terms. When its only instance is deleted (which
    should only happen when a SageMath session ends), then the
    freelist is cleared.
    """
    def __dealloc__(self):
        """
        TESTS::

            sage: s = Sage()
            sage: s.eval("P = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'))")
            ''
            sage: s.eval("P.inject_variables()")
            'Defining e_1, x, y, z'
            sage: s.eval("x*y+y*z*x")
            'x*y + y*z*x'
            sage: s.quit()   # indirect doctest
        """
        cdef size_t i
        for i in range(freelist.used):
            term_free_force(freelist.pool[i])
            sig_check()
        sig_free(freelist.pool)
        sig_free(freelist)

_freelist_protector = _FreeListProtector()

# Put the term on the freelist (unless the list is full),
# and return the pointer .nxt
cdef inline path_term_t *term_free(path_term_t *T) noexcept:
    if T.coef!=NULL:
        Py_XDECREF(T.coef)
    if likely(freelist.used < poolsize):
        freelist.pool[postinc(freelist.used)] = T
        return T.nxt
    return term_free_force(T)

# Create a term by copying the given bounded integer sequence,
# with the given coefficient
cdef path_term_t *term_create(object coef, biseq_t Mon, long Pos, mp_size_t L_len, mp_size_t S_len) except NULL:
    cdef path_term_t *out
    if likely(freelist.used > 0):
        out = freelist.pool[predec(freelist.used)]
        mon_realloc(out.mon, Mon, Pos, L_len, S_len)
    else:
        out = <path_term_t*>check_malloc(sizeof(path_term_t))
        mon_create(out.mon, Mon, Pos, L_len, S_len)
    Py_INCREF(coef)
    out.coef = <PyObject*>coef
    out.nxt = NULL
    return out

# Create a term without copying the given bounded integer sequence
cdef path_term_t *term_create_keep(object coef, biseq_t Mon, long Pos, mp_size_t L_len, mp_size_t S_len) except NULL:
    cdef path_term_t *out
    if likely(freelist.used) > 0:
        out = freelist.pool[predec(freelist.used)]
        mon_realloc_keep(out.mon, Mon, Pos, L_len, S_len)
    else:
        out = <path_term_t*>check_malloc(sizeof(path_term_t))
        mon_create_keep(out.mon, Mon, Pos, L_len, S_len)
    Py_INCREF(coef)
    out.coef = <PyObject*>coef
    # out.nxt = NULL  # to be taken care of externally
    return out

# Create a term with a given coefficient, but empty monomial
cdef path_term_t *term_create_blank(object coef) except NULL:
    cdef path_term_t *out
    if likely(freelist.used > 0):
        out = freelist.pool[predec(freelist.used)]
        mon_free(out.mon)
    else:
        out = <path_term_t*>check_malloc(sizeof(path_term_t))
    Py_INCREF(coef)
    out.coef = <PyObject*>coef
    # out.nxt = NULL  # to be taken care of externally
    return out

######################################################################
######################################################################

# Copy a term; recall that copying the underlying monomial
# just means to increase its reference count. However,
# the copied TERM is new.
# The .nxt attribute is NOT defined on the copy of the term.
cdef path_term_t *term_copy(path_term_t *T) except NULL:
    cdef path_term_t *out
    if likely(freelist.used > 0):
        out = freelist.pool[predec(freelist.used)]
        mon_free(out.mon)
    else:
        out = <path_term_t*>check_malloc(sizeof(path_term_t))
    sig_on()
    mon_copy(out.mon, T.mon)
    sig_off()
    Py_XINCREF(T.coef)
    out.coef = T.coef
    # out.nxt is supposed to be taken care of externally
    return out

# Create a copy of T and recursively of T.nxt
cdef path_term_t *term_copy_recursive(path_term_t *T) except NULL:
    cdef path_term_t *out = term_copy(T)
    cdef path_term_t *first = out
    T = T.nxt
    while T!=NULL:
        out.nxt = term_copy(T)
        out = out.nxt
        T = T.nxt
    out.nxt = NULL
    return first

# Hash of a term; probably not a good one.
cdef inline long term_hash(path_term_t *T) noexcept:
    return (<long>hash(<object>T.coef)+(T.mon.l_len<<5)+(T.mon.pos<<10))^bitset_hash(T.mon.path.data)

# Recall that a monomial a*I*b (with I a generator of a free module)
# is encoded by a path a*s*b for some monomial s that refers to a
# so-called Schreyer ordering. The total degree of a*I*b is the length
# of a plus the length of b.
cdef inline mp_size_t term_total_degree(path_term_t *T) noexcept:
    return T.mon.path.length-T.mon.s_len

# Linearisation
cdef inline tuple term_pickle(path_term_t *T):
    return (<object>T.coef, mon_pickle(T.mon))

# De-linearisation
cdef inline path_term_t *term_unpickle(object coef, tuple mon_data) except NULL:
    cdef path_term_t *out = term_create_blank(coef)
    mon_unpickle(out.mon, mon_data)
    return out

########################################
#
#  Multiplication of monomials

# Return T*p, for a path p and a monomial T.
cdef bint mon_mul_path(path_mon_t out, path_mon_t T, biseq_t p) except -1:
    if unlikely(p.length == 0):
        return mon_copy(out, T)
    biseq_init_concat(out.path, T.path, p)
    out.pos = T.pos
    out.l_len = T.l_len
    out.s_len = T.s_len

# Return p*T, for a path p and a monomial T.
cdef bint path_mul_mon(path_mon_t out, biseq_t p, path_mon_t T) except -1:
    if unlikely(p.length == 0):
        return mon_copy(out, T)
    biseq_init_concat(out.path, p, T.path)
    out.pos = T.pos
    out.l_len = 0 if T.pos==-1 else T.l_len+p.length
    out.s_len = T.s_len

# Return p*T*q, for paths p,q and a monomial T.
cdef bint path_mul_mon_mul_path(path_mon_t out, biseq_t p, path_mon_t T, biseq_t q) except -1:
    # .l_len and .s_len are taken care of externally!
    if unlikely(p.length==0 and q.length==0):
        return mon_copy(out, T)
    if unlikely(p.length == 0):
        return mon_mul_path(out, T, q)
    if unlikely(q.length == 0):
        return path_mul_mon(out, p, T)
    out.pos = T.pos
    if unlikely(T.path.length == 0):
        biseq_init_concat(out.path, p, q)
        return True
    cdef mp_size_t pTlength = p.length + T.path.length
    cdef mp_size_t res_length = pTlength + q.length
    biseq_init(out.path, res_length, p.itembitsize)
    if res_length == 0:
        return False
    cdef mp_bitcnt_t pTsize = p.data.size+T.path.data.size
    sig_on()
    bitset_lshift(out.path.data, q.data, pTsize)
    cdef mp_bitcnt_t p_offset = p.data.size % GMP_LIMB_BITS
    # p_limbs gives the index of the limb that will store the first bit of the
    # shifted version of T.
    cdef mp_bitcnt_t p_limbs = (p.data.limbs - 1) if p_offset>0 else p.data.limbs

    # pT_limbs gives the index of the last limb used to store p+T
    cdef mp_bitcnt_t pT_limbs = (pTsize-1)//GMP_LIMB_BITS
    if ((T.path.data.size-1) % GMP_LIMB_BITS) + p_offset >= GMP_LIMB_BITS:
        # We shift all limbs of T. The highest bits of the highest limbs are
        # pushed out and returned by mpn_lshift. We need to assign them to the
        # beginning of the last limb that is (partially) occupied by p+T
        out.path.data.bits[pT_limbs] |= mpn_lshift(out.path.data.bits+p_limbs,
                                                   T.path.data.bits,
                                                   T.path.data.limbs, p_offset)
    else:
        if T.path.data.limbs>1:
            # If we would move all limbs of T, then the result would override
            # the lowest limb of the shifted copy of q. We thus only move all
            # but the last limb of T, assigning to the beginning of the last
            # limb of p+T the bits that have been pushed out.
            out.path.data.bits[pT_limbs] |= mpn_lshift(out.path.data.bits+p_limbs,
                                                       T.path.data.bits,
                                                       T.path.data.limbs-1,
                                                       p_offset)
            # Last, we need to move the last limb of T (which is only
            # partially occupied), namely into the spot between the previously
            # moved parts of T and the beginning of the shifted copy of q.
            out.path.data.bits[pT_limbs] |= (T.path.data.bits[T.path.data.limbs-1]<<p_offset)
        else:
            out.path.data.bits[p_limbs] |= (T.path.data.bits[T.path.data.limbs-1]<<p_offset)
    bitset_or(out.path.data, out.path.data, p.data)
    sig_off()

########################################
#  Addition and scaling of terms

# Return -T
cdef path_term_t *term_neg(path_term_t *T) except NULL:
    cdef path_term_t *out
    if likely(freelist.used > 0):
        out = freelist.pool[predec(freelist.used)]
        mon_free(out.mon)
    else:
        out = <path_term_t*>check_malloc(sizeof(path_term_t))
    cdef object coef = -<object>T.coef
    out.coef = <PyObject*>coef
    Py_INCREF(coef)
    mon_copy(out.mon, T.mon)
    # out.nxt is supposed to be taken care of externally
    return out

# Return -T, and recurse over T.nxt
cdef path_term_t *term_neg_recursive(path_term_t *T) except NULL:
    cdef path_term_t *out = term_neg(T)
    cdef path_term_t *first = out
    T = T.nxt
    while T!=NULL:
        sig_check()
        out.nxt = term_neg(T)
        out = out.nxt
        T = T.nxt
    out.nxt = NULL
    return first

# Return coef*T
cdef path_term_t *term_scale(path_term_t *T, object coef) except NULL:
    cdef path_term_t *out
    if likely(freelist.used > 0):
        out = freelist.pool[predec(freelist.used)]
        mon_free(out.mon)
    else:
        out = <path_term_t*>check_malloc(sizeof(path_term_t))
    cdef object new_coef = coef*<object>T.coef
    if new_coef:
        out.coef = <PyObject*>new_coef
        Py_INCREF(new_coef)
        mon_copy(out.mon, T.mon)
    else:
        out.coef = NULL
    # out.nxt is supposed to be taken care of externally
    return out

# Return coef*T and recurse over T.nxt
cdef path_term_t *term_scale_recursive(path_term_t *T, object coef) except NULL:
    cdef path_term_t *out = term_scale(T, coef)
    cdef path_term_t *first = out
    T = T.nxt
    while T!=NULL:
        sig_check()
        out.nxt = term_scale(T, coef)
        if out.nxt.coef == NULL:
            term_free(out.nxt)
            out.nxt = NULL
        else:
            out = out.nxt
        T = T.nxt
    out.nxt = NULL
    return first

# Return T1*T2.
# An error is raised if both T1 and T2 belong to a free module over a
# path algebra (but not to the path algebra itself). Hence, this function
# implements multiplication of a term in a path algebra, and the action
# of a term of a path algebra on a term in a free module.
cdef path_term_t *term_mul_term(path_term_t *T1, path_term_t *T2) except NULL:
    cdef mp_size_t new_l_len
    cdef long new_pos
    cdef mp_size_t new_s_len
    if T1.mon.pos!=-1:
        if T2.mon.pos!=-1:
            raise ValueError("we cannot multiply two module elements")
        new_l_len = T1.mon.l_len
        new_pos = T1.mon.pos
        new_s_len = T1.mon.s_len
    elif T2.mon.pos!=-1:
        new_l_len = T2.mon.l_len+T1.mon.path.length
        new_pos = T2.mon.pos
        new_s_len = T2.mon.s_len
    else:
        new_l_len = 0
        new_pos = -1
        new_s_len = 0
    cdef object new_coef = (<object>T1.coef)*(<object>T2.coef)

    cdef path_term_t *out
    if likely(freelist.used > 0):
        out = freelist.pool[predec(freelist.used)]
        if new_coef:
            out.coef = <PyObject*>(new_coef)
            Py_INCREF(new_coef)
            biseq_dealloc(out.mon.path)
            biseq_init_concat(out.mon.path, T1.mon.path, T2.mon.path)
        else:
            out.coef = NULL
    else:
        out = <path_term_t*>check_malloc(sizeof(path_term_t))
        if new_coef:
            out.coef = <PyObject*>(new_coef)
            Py_INCREF(new_coef)
            biseq_init_concat(out.mon.path, T1.mon.path, T2.mon.path)
        else:
            out.coef = NULL
    out.mon.pos = new_pos
    out.mon.l_len = new_l_len
    out.mon.s_len = new_s_len
    out.nxt = NULL
    return out

########################################
#
#  Basics for polynomials

# Create an empty polynomial
cdef inline path_poly_t *poly_create() except NULL:
    cdef path_poly_t *out = <path_poly_t*>check_malloc(sizeof(path_poly_t))
    out.lead = NULL
    out.nterms = 0
    return out

# Deallocate all terms of the polynomial, but NOT the polynomial itself
cdef inline void poly_dealloc(path_poly_t *P) noexcept:
    cdef path_term_t *T = P.lead
    while T!=NULL:
        T = term_free(T)

# Deallocate all terms of the polynomial, and free the chunk of memory
# used by the polynomial.
cdef inline void poly_free(path_poly_t *P) noexcept:
    poly_dealloc(P)
    sig_free(P)

# Fill "out" with a copy of the terms of P. Note that previous contents
# of "out" will NOT be freed---this function should thus only be called
# when "out" is empty.
cdef inline bint poly_icopy(path_poly_t *out, path_poly_t *P) except -1:
    cdef path_term_t *T = P.lead
    out.nterms = P.nterms
    out.lead = term_copy_recursive(T)
    return True

# Fill "out" with a copy of the terms of -P. Note that previous contents
# of "out" will NOT be freed---this function should thus only be called
# when "out" is empty.
cdef inline bint poly_icopy_neg(path_poly_t *out, path_poly_t *P) except -1:
    cdef path_term_t *T = P.lead
    out.nterms = P.nterms
    out.lead = term_neg_recursive(T)
    return True

# Fill "out" with a copy of the terms of coef*P. Note that previous contents
# of "out" will NOT be freed---this function should thus only be called
# when "out" is empty.
cdef bint poly_icopy_scale(path_poly_t *out, path_poly_t *P, object coef) except -1:
    cdef path_term_t *T = P.lead
    cdef path_term_t *res = term_scale(T, coef)
    out.nterms = 0
    out.lead = NULL
    while res.coef == NULL:
        sig_check()
        sig_free(res)
        T = T.nxt
        if T == NULL:
            return True
        res = term_scale(T, coef)
    out.lead = res
    out.nterms += 1
    T = T.nxt
    while T != NULL:
        sig_check()
        res.nxt = term_scale(T, coef)
        if res.nxt.coef == NULL:
            sig_free(res.nxt)
        else:
            res = res.nxt
            out.nterms += 1
        T = T.nxt
    if res != NULL:
        res.nxt = NULL
    return True

# Linearisation of a path polynomials
cdef list poly_pickle(path_poly_t *P):
    cdef list L = []
    cdef path_term_t *T = P.lead
    while T != NULL:
        L.append(term_pickle(T))
        T = T.nxt
    return L

# De-linearisation
cdef bint poly_inplace_unpickle(path_poly_t *P, list data) except -1:
    cdef tuple term_data
    cdef object coef
    cdef path_term_t *T
    if not data:
        P.nterms = 0
        P.lead = NULL
        return True
    P.nterms = len(data)
    coef, term_data = data.pop(0)
    P.lead = term_unpickle(coef, term_data)
    T = P.lead
    for coef, term_data in data:
        T.nxt = term_unpickle(coef, term_data)
        T = T.nxt
    T.nxt = NULL
    return True

############################################
#
#  Polynomial arithmetics

# Rich comparison of P1 and P2, using the given monomial ordering cmp_terms.
# Return a boolean.
cdef bint poly_richcmp(path_poly_t *P1, path_poly_t *P2, path_order_t cmp_terms, int op) noexcept:
    cdef path_term_t *T1 = P1.lead
    cdef path_term_t *T2 = P2.lead
    cdef int c
    cdef object t1
    cdef object t2
    while T1 != NULL and T2 != NULL:
        sig_check()
        c = cmp_terms(T1.mon, T2.mon)
        if c:
            return rich_to_bool(op, c)

        t1 = <object>T1.coef
        t2 = <object>T2.coef
        if t1 != t2:
            return richcmp_not_equal(t1, t2, op)
        T1 = T1.nxt
        T2 = T2.nxt
    if T1 == NULL:
        if T2 == NULL:
            return rich_to_bool(op, 0)
        return rich_to_bool(op, -1)
    return rich_to_bool(op, 1)

# Hash of a polynomial. Probably not a very strong hash.
cdef inline long poly_hash(path_poly_t *P) noexcept:
    cdef path_term_t *T = P.lead
    cdef long out = 0
    while T != NULL:
        out = out<<7 | (out>>(sizeof(long)-7))
        out += term_hash(T)
        T = T.nxt
    return out

# Change T1 inplace to T1+T2.coeff*T1. If the new coefficient is zero,
# then T1.coef becomes NULL
cdef inline void term_iadd(path_term_t *T1, path_term_t *T2) noexcept:
    cdef object coef = <object>(T1.coef) + <object>(T2.coef)
    Py_XDECREF(T1.coef)
    if coef:
        Py_INCREF(coef)
        T1.coef = <PyObject*>coef
    else:
        T1.coef = NULL

# Change P inplace to P+T. It is assumed that initially the terms of P are
# decreasingly sorted wrt. cmp_terms, and then it is guaranteed that they
# are decreasingly sorted wrt. cmp_terms after adding T.
# The addition is "destructive" for T, which means that one MUST NOT
# call term_free(T) after the addition!
cdef bint poly_iadd_term_d(path_poly_t *P, path_term_t *T, path_order_t cmp_terms) except -1:
    if P.lead == NULL:
        P.nterms += 1
        T.nxt = NULL
        P.lead = T
        return True
    cdef path_term_t *tmp = P.lead
    cdef int c
    cdef object coef
    c = cmp_terms(tmp.mon, T.mon)
    if c==-1:
        # The poly's lead term is smaller than T. Hence, we need to prepend
        # it.
        P.nterms += 1
        T.nxt = tmp
        P.lead = T
        return True
    elif c==0:
        sig_on()
        term_iadd(tmp, T)
        term_free(T)
        if tmp.coef==NULL:
            P.nterms -= 1
            P.lead = term_free(tmp)
        elif <object>(tmp.coef)==0:
            sig_off()
            raise RuntimeError("this should never happen")
        sig_off()
        return True
    while True:
        # At this point, we have tmp>T.
        #
        # We need to append the term, or continue until we can
        # insert/append
        sig_check()
        if tmp.nxt == NULL:
            P.nterms += 1
            T.nxt = NULL
            tmp.nxt = T
            return True
        c = cmp_terms(tmp.nxt.mon, T.mon)
        if c==-1:
            P.nterms += 1
            T.nxt = tmp.nxt
            tmp.nxt = T
            return True
        elif c==0:
            term_iadd(tmp.nxt, T)
            term_free(T)
            if tmp.nxt.coef==NULL:
                P.nterms -= 1
                tmp.nxt = term_free(tmp.nxt)
            elif <object>(tmp.coef)==0:
                raise RuntimeError("this should never happen")
            return True
        # otherwise, tmp is still larger than T. Hence, move to the next term
        # of P.
        tmp = tmp.nxt

# Change P1 inplace to P1+P2. It is assumed that both P1's and P2's terms
# are decreasingly sorted wrt. cmp_terms, and after the operation P1's terms
# will still be decreasingly sorted. After the operation, one MUST NOT
# call poly_free(P2)!
cdef bint poly_iadd_d(path_poly_t *P1, path_poly_t *P2, path_order_t cmp_terms) except -1:
    # Terms of P2 will be moved to P1, so that deallocation of P2 will not be
    # needed.
    if P1.lead == NULL:
        P1.lead = P2.lead
        P1.nterms = P2.nterms
        P2.nterms = 0
        P2.lead = NULL
        return 1
    if P2.lead == NULL:
        return 1
    cdef path_term_t *T1 = P1.lead
    cdef path_term_t *T2 = P2.lead
    cdef path_term_t *prev = NULL
    cdef int c
    cdef object new_coef
    while True:
        # Is one of the summands consumed already? Then we can use an easier
        # method.
        sig_check()
        if T1 == NULL:
            if prev==NULL:
                P1.lead = T2
            else:
                prev.nxt = T2
            P1.nterms += P2.nterms
            P2.nterms = 0
            P2.lead = NULL
            return 1
        elif T2 == NULL:
            if P2.nterms != 0:
                print("term counting of second summand was wrong! " +
                      str(P2.nterms))
            P2.lead = NULL
            return 1
        c = cmp_terms(T1.mon, T2.mon)
        if c == 0:
            # T1==T2 --- We need to add
            new_coef = <object>(T1.coef)+<object>(T2.coef)
            if new_coef:
                Py_INCREF(new_coef)
                Py_XDECREF(T1.coef)
                T1.coef = <PyObject*>new_coef
                prev = T1
                T1 = T1.nxt
            else:
                T1 = term_free(T1)
                if prev==NULL:
                    P1.lead = T1
                else:
                    prev.nxt = T1
                P1.nterms -= 1
            P2.nterms -= 1
            T2 = term_free(T2)
        elif c == 1:
            # We move the prev/T1 through P1, until prev>T2>=T1. But now,
            # T1>T2. So, we move prev/T1 further down.
            prev = T1
            T1 = T1.nxt
        else:
            # prev > T2 > T1. Hence, we insert T2, without copying
            if prev==NULL:
                P1.lead = T2
                T2 = T2.nxt
                prev = P1.lead
            else:
                prev.nxt = T2
                T2 = T2.nxt
                prev = prev.nxt
            prev.nxt = T1
            P1.nterms += 1
            P2.nterms -= 1

# Return P1+P2 (a new polynomial, and after the operation it is still safe
# to call poly_free(P2)). Both P1's and P2's terms are supposed to be
# decreasingly sorted wrt. cmp_terms, and so will be the terms of P1+P2.
cdef path_poly_t *poly_add(path_poly_t *P1, path_poly_t *P2, path_order_t cmp_terms) except NULL:
    cdef path_poly_t *out = poly_create()
    cdef path_term_t *T1 = P1.lead
    cdef path_term_t *T2 = P2.lead
    cdef path_term_t *T = NULL
    cdef path_term_t *res
    cdef size_t count1, count2  # How many terms of P1/P2 have been considered?
    count1 = 0
    count2 = 0
    cdef object coef
    cdef int c
    while True:
        sig_check()
        if T1 == NULL:
            out.nterms += (P2.nterms-count2)
            if T == NULL:
                if T2 == NULL:
                    out.lead = NULL
                else:
                    out.lead = term_copy_recursive(T2)
            else:
                if T2 == NULL:
                    T.nxt = NULL
                else:
                    T.nxt = term_copy_recursive(T2)
            return out
        if T2 == NULL:
            out.nterms += (P1.nterms-count1)
            if T == NULL:
                out.lead = term_copy_recursive(T1)
            else:
                T.nxt = term_copy_recursive(T1)
            return out

        c = cmp_terms(T1.mon, T2.mon)
        if c == 1:
            if T == NULL:
                out.lead = term_copy(T1)
                T = out.lead
            else:
                T.nxt = term_copy(T1)
                T = T.nxt
            T1 = T1.nxt
            count1 += 1
            out.nterms += 1
        elif c == -1:
            if T == NULL:
                out.lead = term_copy(T2)
                T = out.lead
            else:
                T.nxt = term_copy(T2)
                T = T.nxt
            T2 = T2.nxt
            count2 += 1
            out.nterms += 1
        else:
            coef = (<object>T1.coef)+(<object>T2.coef)
            if coef:
                out.nterms += 1
                if T == NULL:
                    out.lead = term_create(coef, T1.mon.path, T1.mon.pos, T1.mon.l_len, T1.mon.s_len)
                    T = out.lead
                else:
                    T.nxt = term_create(coef, T1.mon.path, T1.mon.pos, T1.mon.l_len, T1.mon.s_len)
                    T = T.nxt
            count1 += 1
            count2 += 1
            T1 = T1.nxt
            T2 = T2.nxt

# Return P1-P2 (a new polynomial, and after the operation it is still safe
# to call poly_free(P2)). Both P1's and P2's terms are supposed to be
# decreasingly sorted wrt. cmp_terms, and so will be the terms of P1-P2.
cdef path_poly_t *poly_sub(path_poly_t *P1, path_poly_t *P2, path_order_t cmp_terms) except NULL:
    cdef path_poly_t *out = poly_create()
    cdef path_term_t *T1 = P1.lead
    cdef path_term_t *T2 = P2.lead
    cdef path_term_t *T = NULL
    cdef path_term_t *res
    cdef size_t count1, count2  # How many terms of P1/P2 have been considered?
    count1 = 0
    count2 = 0
    cdef object coef
    cdef int c
    while True:
        sig_check()
        if T1 == NULL:
            out.nterms += (P2.nterms-count2)
            if T == NULL:
                if T2 == NULL:
                    out.lead = NULL
                else:
                    out.lead = term_neg_recursive(T2)
            else:
                if T2 == NULL:
                    T.nxt = NULL
                else:
                    T.nxt = term_neg_recursive(T2)
            return out
        if T2 == NULL:
            out.nterms += (P1.nterms-count1)
            if T == NULL:
                out.lead = term_copy_recursive(T1)
            else:
                T.nxt = term_copy_recursive(T1)
            return out

        c = cmp_terms(T1.mon, T2.mon)
        if c == 1:
            if T == NULL:
                out.lead = term_copy(T1)
                T = out.lead
            else:
                T.nxt = term_copy(T1)
                T = T.nxt
            T1 = T1.nxt
            count1 += 1
            out.nterms += 1
        elif c == -1:
            if T == NULL:
                out.lead = term_neg(T2)
                T = out.lead
            else:
                T.nxt = term_neg(T2)
                T = T.nxt
            T2 = T2.nxt
            count2 += 1
            out.nterms += 1
        else:
            coef = (<object>T1.coef)-(<object>T2.coef)
            if coef:
                out.nterms += 1
                if T == NULL:
                    out.lead = term_create(coef, T1.mon.path, T1.mon.pos, T1.mon.l_len, T1.mon.s_len)
                    T = out.lead
                else:
                    T.nxt = term_create(coef, T1.mon.path, T1.mon.pos, T1.mon.l_len, T1.mon.s_len)
                    T = T.nxt
            count1 += 1
            count2 += 1
            T1 = T1.nxt
            T2 = T2.nxt

#
#  In-place addition of a multiple of a polynomial
# Replace P1 by P1+coef*P2*R. Return a pointer to the first term of P1
# that may be involved in a change when calling the function again with
# P1, P2 and a cofactor that is smaller than R wrt. cmp_terms.
# The return value should then be provided as argument "P1start" of the
# next function call.
#
# We return P1start if P2.lead is NULL. Otherwise, if P1.lead becomes NULL
# during addition, then we return P2.lead.
#
# Let m be a monomial of P2. If it is of module-type (i.e., m.pos!=-1),
# then it is clear what m2=m*R is (m2.l_len and m2.s_len are obtained from
# m.l_len and m.s_len). Otherwise, however, it could be that we
# want m2=m*R to denote a module-type monomial. In that case, "m2.l_len"
# and "m2.s_len" are given by the respective arguments.

cdef path_term_t *poly_iadd_lmul(path_poly_t *P1, object coef, path_poly_t *P2, biseq_t R, path_order_t cmp_terms, long pos, mp_size_t l_len, mp_size_t s_len, path_term_t *P1start) except NULL:
    if not coef or P2.lead==NULL:
        return P1start
    cdef path_mon_t new_mon
    cdef object new_coef
    cdef path_term_t *prev = NULL
    cdef path_term_t *T1
    if P1start == NULL:
        T1 = P1.lead
    else:
        T1 = P1start
    cdef path_term_t *T2 = P2.lead
    cdef int c
    cdef path_term_t *out = P1start
    while T2!=NULL:
        sig_check()
        new_coef = coef*<object>(T2.coef)
        if not new_coef:
            T2 = T2.nxt
            continue

        if T2.mon.pos!=-1:
            mon_mul_path(new_mon, T2.mon, R)
            new_mon.pos = T2.mon.pos
            new_mon.l_len = T2.mon.l_len
            new_mon.s_len = T2.mon.s_len
        else:
            mon_mul_path(new_mon, T2.mon, R)
            new_mon.pos = pos
            new_mon.l_len = l_len
            new_mon.s_len = s_len
        # Now new_term is T2*R
        # We go down in P1 until we may append, insert or add
        while T1!=NULL:
            sig_check()
            c = cmp_terms(T1.mon, new_mon)
            if c!=1:
                break
            prev = T1
            T1 = prev.nxt
        if T1==NULL:
            # We need to append to P1
            if prev==NULL:
                P1.lead = term_create_blank(new_coef)
                P1.lead.mon[0] = new_mon[0]
                prev = P1.lead
            else:
                prev.nxt = term_create_blank(new_coef)
                prev.nxt.mon[0] = new_mon[0]
                prev = prev.nxt
            if T2 == P2.lead:
                out = prev
            prev.nxt = NULL
            P1.nterms += 1
        elif c==-1:
            # We need to insert between prev and T1
            if prev==NULL:
                P1.lead = term_create_blank(new_coef)
                P1.lead.mon[0] = new_mon[0]
                prev = P1.lead
            else:
                prev.nxt = term_create_blank(new_coef)
                prev.nxt.mon[0] = new_mon[0]
                prev = prev.nxt
            if T2 == P2.lead:
                out = prev
            prev.nxt = T1
            P1.nterms += 1
        else:
            # we add the coefficients and see what happens
            new_coef+=<object>(T1.coef)
            mon_free(new_mon)
            if new_coef:
                Py_INCREF(new_coef)
                Py_XDECREF(T1.coef)
                T1.coef = <PyObject*>new_coef
                if T2 == P2.lead:
                    out = T1
                prev = T1
                T1 = T1.nxt
            else:
                P1.nterms -= 1
                T1 = term_free(T1)
                if prev==NULL:
                    P1.lead = T1
                else:
                    prev.nxt = T1
                if T2 == P2.lead:
                    out = prev
        T2 = T2.nxt
    if out == NULL:
        return P2.lead
    return out

########################################
#
#  Basics for homogeneous polynomials

# Create an empty polynomial whose to-be-inserted terms
# have start- and end-points of the given integer labels
# start and end.
cdef path_homog_poly_t *homog_poly_create(int start, int end) except NULL:
    cdef path_homog_poly_t *out = <path_homog_poly_t*>check_malloc(sizeof(path_homog_poly_t))
    out.poly = poly_create()
    out.start = start
    out.end = end
    out.nxt = NULL
    return out

# Create a new homogeneous polynomial from a given polynomial P.
# It is assumed that all terms of P have start- and end-points
# with integer labels start and end. This assumption is NOT checked!
# P is inserted, not copied.
cdef path_homog_poly_t *homog_poly_init_poly(int start, int end, path_poly_t *P) except NULL:
    cdef path_homog_poly_t *out = <path_homog_poly_t*>check_malloc(sizeof(path_homog_poly_t))
    out.poly = P
    out.start = start
    out.end = end
    out.nxt = NULL
    return out

# L provides a list of pairs (P,coef), where P is a path
# and coef the coefficient of the corresponding term.
# With this function, one can create elements of a path algebra (pos==-1),
# or elements of free modules over a path algebra in summand pos.
# In either case, the length of left cofactors of each monomial will
# be zero, and also the length of Schreyer monomials will be zero.
cdef path_homog_poly_t *homog_poly_init_list(int start, int end, list L, path_order_t cmp_terms, long pos) except NULL:
    cdef path_homog_poly_t * out = homog_poly_create(start, end)
    cdef QuiverPath P
    for P, coef in L:
        poly_iadd_term_d(out.poly, term_create(coef, P._path, pos, 0, 0), cmp_terms)
    return out

cdef void homog_poly_free(path_homog_poly_t *P) noexcept:
    cdef path_homog_poly_t *nxt
    while P!=NULL:
        nxt = P.nxt
        poly_free(P.poly)
        sig_free(P)
        P = nxt

# Return a copy of H
cdef path_homog_poly_t *homog_poly_copy(path_homog_poly_t *H) except NULL:
    cdef path_homog_poly_t *out
    cdef path_homog_poly_t *tmp
    if H == NULL:
        raise ValueError("the polynomial to be copied is the NULL pointer")
    out = homog_poly_create(H.start, H.end)
    poly_icopy(out.poly, H.poly)
    tmp = out
    H = H.nxt
    while H != NULL:
        sig_check()
        tmp.nxt = homog_poly_create(H.start, H.end)
        tmp = tmp.nxt
        poly_icopy(tmp.poly, H.poly)
        H = H.nxt
    return out

# Linearisation
cdef list homog_poly_pickle(path_homog_poly_t *H):
    cdef list L = []
    while H != NULL:
        L.append((H.start, H.end, poly_pickle(H.poly)))
        H = H.nxt
    return L

# De-linearisation
cdef path_homog_poly_t *homog_poly_unpickle(list data) except NULL:
    # ASSUMPTION: data is not empty
    cdef int start, end
    cdef list poly_data
    cdef path_homog_poly_t *out
    start, end, poly_data = data.pop(0)
    out = homog_poly_create(start, end)
    poly_inplace_unpickle(out.poly, poly_data)
    cdef path_homog_poly_t *tmp = out
    for start, end, poly_data in data:
        sig_check()
        tmp.nxt = homog_poly_create(start, end)
        tmp = tmp.nxt
        poly_inplace_unpickle(tmp.poly, poly_data)
    return out

# Return -H
cdef path_homog_poly_t *homog_poly_neg(path_homog_poly_t *H) except NULL:
    cdef path_homog_poly_t *out
    cdef path_homog_poly_t *tmp
    if H == NULL:
        raise ValueError("the polynomial to be copied is the NULL pointer")
    out = homog_poly_create(H.start, H.end)
    poly_icopy_neg(out.poly, H.poly)
    tmp = out
    H = H.nxt
    while H != NULL:
        sig_check()
        tmp.nxt = homog_poly_create(H.start, H.end)
        tmp = tmp.nxt
        poly_icopy_neg(tmp.poly, H.poly)
        H = H.nxt
    return out

# Return coef*H
cdef path_homog_poly_t *homog_poly_scale(path_homog_poly_t *H, object coef) except NULL:
    # The first component may be zero, all other zero components are removed.
    cdef path_homog_poly_t *out
    cdef path_homog_poly_t *tmp
    if H == NULL:
        raise ValueError("the polynomial to be copied is the NULL pointer")
    out = homog_poly_create(H.start, H.end)
    poly_icopy_scale(out.poly, H.poly, coef)
    tmp = out
    H = H.nxt
    while H != NULL:
        sig_check()
        tmp.nxt = homog_poly_create(H.start, H.end)
        poly_icopy_scale(tmp.nxt.poly, H.poly, coef)
        if tmp.nxt.poly.nterms == 0:
            homog_poly_free(tmp.nxt)
            tmp.nxt = NULL
        else:
            tmp = tmp.nxt
        H = H.nxt
    return out

cdef path_homog_poly_t *homog_poly_get_predecessor_of_component(path_homog_poly_t *H, int s, int e) noexcept:
    # Search through H.nxt.nxt... and return the pointer C to a component of H
    # such that either C.nxt.start==s and C.nxt.end==e, or the component for
    # (s,e) should be inserted between C and C.nxt. Return NULL if H==NULL or
    # (s,e) should be inserted in front of H.
    if H == NULL:
        return NULL
    if H.start > s:
        return NULL
    elif H.start == s and H.end >= e:
        return NULL
    while True:
        sig_check()
        if H.nxt == NULL:
            return H
        if H.nxt.start == s:
            if H.nxt.end >= e:
                return H
        elif H.nxt.start > s:
            return H
        H = H.nxt
