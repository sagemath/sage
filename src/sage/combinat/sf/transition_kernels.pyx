r"""
Kostka and Murnaghan-Nakayama kernels for the classical bases

This module computes transitions between the Schur, complete homogeneous,
elementary, monomial and power sum bases one basis element at a time.
Partitions are plain tuples of Python integers and the results are
dictionaries mapping partitions to coefficients.

The combinatorial kernels are:

- Kostka numbers `K_{\lambda\mu}`, computed by adding (Pieri rule) or
  removing horizontal strips;

- inverse Kostka numbers, computed by removing (or adding) special border
  strips, i.e., by expanding the Jacobi-Trudi determinant;

- irreducible character values `\chi^\lambda(\mu)`, computed with the
  Murnaghan-Nakayama rule by adding or removing border strips, which are
  located using beta-numbers.

The module also provides direct transitions between the complete
homogeneous, elementary and power sum bases and between the power sum and
monomial bases, products of monomial symmetric functions, the exponents of
monomial symmetric polynomials, the Schur expansion of the modified
Hall-Littlewood functions and the semistandard tableaux of given shape and
content.

Intermediate results are memoized; call :func:`clear_caches` to free them.
The dictionaries returned by the public functions are fresh copies.

.. WARNING::

    The arguments are not validated; partitions must be given as weakly
    decreasing sequences of nonnegative integers.

AUTHORS:

- Mike Hansen (2026): initial version
"""

from cpython.ref cimport Py_INCREF
from cpython.tuple cimport PyTuple_New, PyTuple_SET_ITEM
from cysignals.memory cimport check_allocarray, check_calloc, sig_free

from math import factorial

from sage.rings.integer import Integer

cdef dict _kostka_cache = {}
cdef dict _char_cache = {}
cdef dict _h_to_s_cache = {}
cdef dict _p_to_s_cache = {}
cdef dict _s_to_m_cache = {}
cdef dict _m_to_s_cache = {}
cdef dict _s_to_h_cache = {}
cdef dict _partitions_cache = {}
cdef dict _gen_cache = {}
cdef dict _mult_cache = {}
cdef dict _p_to_m_cache = {}
cdef dict _m_to_p_cache = {}
cdef dict _mm_cache = {}
cdef dict _jing_cache = {}
cdef dict _qp_cache = {}
cdef object _ZZt = None


def clear_caches():
    r"""
    Clear all memoized intermediate results.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import clear_caches, kostka_number
        sage: kostka_number((2, 1), (1, 1, 1))
        2
        sage: clear_caches()
    """
    for c in (_kostka_cache, _char_cache, _h_to_s_cache, _p_to_s_cache,
              _s_to_m_cache, _m_to_s_cache, _s_to_h_cache, _partitions_cache,
              _gen_cache, _mult_cache, _p_to_m_cache, _m_to_p_cache,
              _mm_cache, _jing_cache, _qp_cache):
        c.clear()


##############################################################################
# Low level helpers
##############################################################################

cdef inline tuple _to_tuple(int* a, Py_ssize_t k):
    """
    Return the tuple ``(a[0], ..., a[k-1])`` with trailing zeros removed.
    """
    while k > 0 and a[k - 1] == 0:
        k -= 1
    cdef tuple t = PyTuple_New(k)
    cdef Py_ssize_t i
    cdef object x
    for i in range(k):
        x = a[i]
        Py_INCREF(x)
        PyTuple_SET_ITEM(t, i, x)
    return t


cdef tuple _normalize(mu):
    """
    Return ``mu`` as a weakly decreasing tuple of positive integers.
    """
    return tuple(sorted([int(p) for p in mu if p], reverse=True))


cdef bint _dominates(tuple lam, tuple mu) noexcept:
    """
    Return whether ``lam`` dominates ``mu``; both have the same size.
    """
    cdef Py_ssize_t i, L = len(lam)
    cdef long s = 0, t = 0
    for i in range(len(mu)):
        if i < L:
            s += <long> lam[i]
        t += <long> mu[i]
        if s < t:
            return False
    return True


cdef tuple _conjugate(tuple lam):
    cdef Py_ssize_t L = len(lam), i, j
    if not L:
        return ()
    cdef int m = lam[0]
    cdef int* c = <int*> check_calloc(m, sizeof(int))
    try:
        for i in range(L):
            for j in range(<int> lam[i]):
                c[j] += 1
        return _to_tuple(c, m)
    finally:
        sig_free(c)


cdef inline void _add_term(dict D, key, c) noexcept:
    """
    Add ``c`` to ``D[key]``, removing the entry if it cancels.
    """
    v = D.get(key, 0) + c
    if v:
        D[key] = v
    else:
        D.pop(key, None)


##############################################################################
# Horizontal strips
##############################################################################

cdef int _add_hstrips(int* lam, int* nu, Py_ssize_t L, Py_ssize_t i, int rem,
                      list out) except -1:
    # lam[0..L-1] is a partition with sentinel lam[L] = 0; rows i..L of nu
    # are still to be chosen and rem cells remain to be added. The rows
    # after row i can absorb at most lam[i] cells.
    cdef int a, cap, lo
    if i == L:
        if L == 0 or rem <= lam[L - 1]:
            nu[L] = rem
            out.append(_to_tuple(nu, L + 1))
        return 0
    cap = rem if i == 0 else lam[i - 1] - lam[i]
    if cap > rem:
        cap = rem
    lo = rem - lam[i]
    if lo < 0:
        lo = 0
    for a in range(cap, lo - 1, -1):
        nu[i] = lam[i] + a
        _add_hstrips(lam, nu, L, i + 1, rem - a, out)
    return 0


cdef int _remove_hstrips(int* lam, int* nu, Py_ssize_t L, Py_ssize_t i,
                         int rem, list out) except -1:
    # The rows after row i can lose at most lam[i + 1] cells.
    cdef int a, cap, lo
    if i == L:
        if rem == 0:
            out.append(_to_tuple(nu, L))
        return 0
    cap = lam[i] - lam[i + 1]
    if cap > rem:
        cap = rem
    lo = rem - lam[i + 1]
    if lo < 0:
        lo = 0
    for a in range(lo, cap + 1):
        nu[i] = lam[i] - a
        _remove_hstrips(lam, nu, L, i + 1, rem - a, out)
    return 0


cdef list _hstrips(tuple lam, int r, bint add):
    """
    Return the partitions obtained from ``lam`` by adding (or removing)
    a horizontal strip of size ``r``.
    """
    cdef Py_ssize_t L = len(lam), i
    cdef int* a = <int*> check_allocarray(L + 1, sizeof(int))
    cdef int* b
    try:
        b = <int*> check_allocarray(L + 1, sizeof(int))
    except MemoryError:
        sig_free(a)
        raise
    cdef list out = []
    try:
        for i in range(L):
            a[i] = lam[i]
        a[L] = 0
        if add:
            _add_hstrips(a, b, L, 0, r, out)
        else:
            _remove_hstrips(a, b, L, 0, r, out)
    finally:
        sig_free(a)
        sig_free(b)
    return out


##############################################################################
# Border strips via beta-numbers
##############################################################################

cdef list _border_strips(tuple lam, int r, bint add):
    r"""
    Return the pairs ``(nu, sign)`` where ``nu`` is obtained from ``lam``
    by adding (or removing) a border strip of size ``r`` and ``sign`` is
    `(-1)^{height}`.

    With `k` beads at positions `\beta_i = \lambda_i + k - 1 - i`, adding
    a border strip moves a bead from `b` to the empty position `b + r`; the
    height is the number of beads strictly in between.
    """
    cdef Py_ssize_t L = len(lam)
    cdef Py_ssize_t k = L + r if add else L
    cdef Py_ssize_t i, j, q, h
    cdef int b, t
    cdef list out = []
    if k == 0:
        return out
    cdef int* part = <int*> check_calloc(3 * k, sizeof(int))
    cdef int* beta = part + k
    cdef int* nu = part + 2 * k
    try:
        for i in range(L):
            part[i] = lam[i]
        for i in range(k):
            beta[i] = part[i] + <int> (k - 1 - i)
        for i in range(k):
            b = beta[i]
            if add:
                t = b + r
                j = i - 1
                while j >= 0 and beta[j] < t:
                    j -= 1
                if j >= 0 and beta[j] == t:
                    continue
                h = i - 1 - j
                for q in range(k):
                    nu[q] = part[q]
                nu[j + 1] = part[i] + r - <int> h
                for q in range(j + 2, i + 1):
                    nu[q] = part[q - 1] + 1
            else:
                t = b - r
                if t < 0:
                    continue
                j = i + 1
                while j < k and beta[j] > t:
                    j += 1
                if j < k and beta[j] == t:
                    continue
                h = j - i - 1
                for q in range(k):
                    nu[q] = part[q]
                for q in range(i, j - 1):
                    nu[q] = part[q + 1] - 1
                nu[j - 1] = part[i] - r + <int> h
            out.append((_to_tuple(nu, k), -1 if h & 1 else 1))
    finally:
        sig_free(part)
    return out


##############################################################################
# Partitions
##############################################################################

cdef int _dominated_rec(long* Lam, int n, int* mu, Py_ssize_t i, int maxpart,
                        int s, list out) except -1:
    cdef int p, top
    if s == n:
        out.append(_to_tuple(mu, i))
        return 0
    top = maxpart
    if n - s < top:
        top = n - s
    if Lam[i] - s < top:
        top = <int> (Lam[i] - s)
    for p in range(top, 0, -1):
        mu[i] = p
        _dominated_rec(Lam, n, mu, i + 1, p, s + p, out)
    return 0


cdef list _dominated_partitions(tuple lam):
    """
    Return the partitions dominated by ``lam`` in reverse lexicographic
    order.
    """
    cdef int n = sum(lam)
    cdef Py_ssize_t L = len(lam), i
    cdef list out = []
    if n == 0:
        return [()]
    cdef long* Lam = <long*> check_allocarray(n, sizeof(long))
    cdef int* mu
    try:
        mu = <int*> check_allocarray(n, sizeof(int))
    except MemoryError:
        sig_free(Lam)
        raise
    try:
        Lam[0] = lam[0]
        for i in range(1, n):
            Lam[i] = Lam[i - 1] + (<long> lam[i] if i < L else 0)
        _dominated_rec(Lam, n, mu, 0, n, 0, out)
    finally:
        sig_free(Lam)
        sig_free(mu)
    return out


cdef list _partitions(int n):
    cdef list P = _partitions_cache.get(n)
    if P is None:
        P = _dominated_partitions((n,) if n else ())
        _partitions_cache[n] = P
    return P


def z(mu):
    r"""
    Return `z_\mu = \prod_i i^{m_i} m_i!`, the size of the centralizer of
    a permutation of cycle type ``mu``.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import z
        sage: z((2, 1, 1)), z(()), z((3,))
        (4, 1, 3)
    """
    mu = _normalize(mu)
    cdef Py_ssize_t i = 0, j, k = len(mu)
    res = 1
    while i < k:
        j = i
        while j < k and mu[j] == mu[i]:
            j += 1
        res *= mu[i] ** (j - i) * factorial(j - i)
        i = j
    return Integer(res)


cdef list _partitions_z(int n):
    """
    Return the pairs ``(mu, z(mu))`` for the partitions ``mu`` of ``n``.
    """
    key = ('z', n)
    cdef list P = _partitions_cache.get(key)
    if P is None:
        P = [(mu, z(mu)) for mu in _partitions(n)]
        _partitions_cache[key] = P
    return P


##############################################################################
# Kostka numbers
##############################################################################

cdef object _kostka(tuple lam, tuple mu):
    # mu weakly decreasing, |lam| == |mu|; remove the smallest part of mu
    # as a horizontal strip.
    cdef Py_ssize_t k = len(mu)
    if not _dominates(lam, mu):
        return 0
    if k <= 1 or len(lam) <= 1:
        return 1
    key = (lam, mu)
    c = _kostka_cache.get(key)
    if c is not None:
        return c
    cdef tuple rest = mu[:k - 1]
    total = 0
    for nu in _hstrips(lam, mu[k - 1], False):
        total += _kostka(<tuple> nu, rest)
    _kostka_cache[key] = total
    return total


def kostka_number(lam, mu):
    r"""
    Return the Kostka number `K_{\lambda\mu}`, the number of semistandard
    tableaux of shape ``lam`` and content ``mu``.

    The content ``mu`` may be any composition.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import kostka_number
        sage: kostka_number((2, 1), (1, 1, 1))
        2
        sage: kostka_number((3, 2, 1), (1, 2, 0, 3))
        1
        sage: kostka_number((3, 2, 1), (1, 1, 1, 1, 1, 1))
        16
        sage: kostka_number((2, 2), (3, 1))
        0

    TESTS::

        sage: all(kostka_number(la, mu)
        ....:     == SemistandardTableaux(la, mu).cardinality()
        ....:     for n in range(7) for la in Partitions(n) for mu in Partitions(n))
        True
    """
    lam = tuple(int(p) for p in lam if p)
    mu = _normalize(mu)
    if sum(lam) != sum(mu):
        return Integer(0)
    return Integer(_kostka(lam, mu))


cdef dict _h_to_s(tuple mu):
    # Iterated Pieri rule; the largest parts are added first so that the
    # cached prefixes are shared by many partitions.
    D = _h_to_s_cache.get(mu)
    if D is not None:
        return D
    cdef Py_ssize_t k = len(mu)
    if k == 0:
        D = {(): 1}
    else:
        D = {}
        r = mu[k - 1]
        for lam, c in _h_to_s(mu[:k - 1]).items():
            for nu in _hstrips(lam, r, True):
                D[nu] = D.get(nu, 0) + c
    _h_to_s_cache[mu] = D
    return D


cdef dict _s_to_m(tuple lam):
    D = _s_to_m_cache.get(lam)
    if D is not None:
        return D
    D = {}
    for mu in _dominated_partitions(lam):
        c = _kostka(lam, <tuple> mu)
        if c:
            D[mu] = c
    _s_to_m_cache[lam] = D
    return D


##############################################################################
# Inverse Kostka numbers via special border strips
##############################################################################

# A special border strip of `\lambda` (of length `\ell`) is the border strip
# containing the lowest cell of the first column and ending at the end of
# row `i`. It has size `\lambda_i + \ell - i` and height `\ell - i`, and
# removing it leaves
# `\nu = (\lambda_1, \ldots, \lambda_{i-1}, \lambda_{i+1} - 1, \ldots, \lambda_\ell - 1)`.
# Expanding the Jacobi-Trudi determinant along its last row gives
#
#     s_lam = sum_i (-1)^{ell - i} h_{lam_i + ell - i} s_nu,
#
# so the inverse Kostka numbers count signed special border strip tabloids
# (Egecioglu-Remmel).

cdef list _remove_special_strips(tuple lam):
    """
    Return the triples ``(nu, a, sign)`` obtained by removing a special
    border strip of size ``a`` from ``lam``.
    """
    cdef Py_ssize_t L = len(lam), i, r
    cdef list out = []
    if not L:
        return out
    cdef int* part = <int*> check_allocarray(2 * L, sizeof(int))
    cdef int* nu = part + L
    try:
        for i in range(L):
            part[i] = lam[i]
        for i in range(L):
            for r in range(i):
                nu[r] = part[r]
            for r in range(i, L - 1):
                nu[r] = part[r + 1] - 1
            out.append((_to_tuple(nu, L - 1), part[i] + <int> (L - 1 - i),
                        -1 if (L - 1 - i) & 1 else 1))
    finally:
        sig_free(part)
    return out


cdef list _add_special_strips(tuple nu, int a):
    """
    Return the pairs ``(lam, sign)`` such that removing a special border
    strip of size ``a`` from ``lam`` gives ``nu``.
    """
    cdef Py_ssize_t m = len(nu), L, i, r
    cdef int li, below
    cdef list out = []
    cdef int* part = <int*> check_calloc(2 * (m + a) + 1, sizeof(int))
    cdef int* lam = part + m + a
    try:
        for i in range(m):
            part[i] = nu[i]
        # lam has length L and the strip ends in row i (0-based), so that
        # lam_i = a - (L - 1 - i) and lam_r = nu_{r-1} + 1 for r > i.
        for L in range(m + 1, m + a + 1):
            for i in range(m + 1 if m + 1 < L else L):
                li = a - <int> (L - 1 - i)
                below = part[i] + 1
                if li < below or (i > 0 and li > part[i - 1]):
                    continue
                for r in range(i):
                    lam[r] = part[r]
                lam[i] = li
                for r in range(i + 1, L):
                    lam[r] = part[r - 1] + 1
                out.append((_to_tuple(lam, L), -1 if (L - 1 - i) & 1 else 1))
    finally:
        sig_free(part)
    return out


cdef tuple _insert_part(tuple mu, a):
    cdef Py_ssize_t j = 0, k = len(mu)
    while j < k and mu[j] >= a:
        j += 1
    return mu[:j] + (a,) + mu[j:]


cdef dict _s_to_h(tuple lam):
    D = _s_to_h_cache.get(lam)
    if D is not None:
        return D
    if not lam:
        D = {(): 1}
    else:
        D = {}
        for nu, a, sign in _remove_special_strips(lam):
            for mu, c in _s_to_h(<tuple> nu).items():
                _add_term(D, _insert_part(<tuple> mu, a), c if sign > 0 else -c)
    _s_to_h_cache[lam] = D
    return D


cdef dict _m_to_s(tuple mu):
    # The coefficient of s_lam in m_mu is the coefficient of h_mu in s_lam,
    # so m_mu = sum_a A_a(m_{mu - a}) over the distinct parts a of mu, where
    # A_a adds a special border strip of size a.
    D = _m_to_s_cache.get(mu)
    if D is not None:
        return D
    cdef Py_ssize_t j, k = len(mu)
    if not k:
        D = {(): 1}
    else:
        D = {}
        for j in range(k):
            if j and mu[j] == mu[j - 1]:
                continue
            for nu, c in _m_to_s(mu[:j] + mu[j + 1:]).items():
                for lam, sign in _add_special_strips(<tuple> nu, mu[j]):
                    _add_term(D, lam, c if sign > 0 else -c)
    _m_to_s_cache[mu] = D
    return D


def h_to_s(mu):
    r"""
    Return the Schur expansion of `h_\mu`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import h_to_s
        sage: sorted(h_to_s((2, 1)).items())
        [((2, 1), 1), ((3,), 1)]
    """
    return dict(_h_to_s(_normalize(mu)))


def e_to_s(mu):
    r"""
    Return the Schur expansion of `e_\mu`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import e_to_s
        sage: sorted(e_to_s((2, 1)).items())
        [((1, 1, 1), 1), ((2, 1), 1)]
    """
    return {_conjugate(<tuple> la): c
            for la, c in _h_to_s(_normalize(mu)).items()}


def s_to_m(lam):
    r"""
    Return the monomial expansion of `s_\lambda`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import s_to_m
        sage: sorted(s_to_m((2, 1)).items())
        [((1, 1, 1), 2), ((2, 1), 1)]
    """
    return dict(_s_to_m(_normalize(lam)))


def m_to_s(mu):
    r"""
    Return the Schur expansion of `m_\mu`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import m_to_s
        sage: sorted(m_to_s((2, 1)).items())
        [((1, 1, 1), -2), ((2, 1), 1)]
    """
    return dict(_m_to_s(_normalize(mu)))


def s_to_h(lam):
    r"""
    Return the complete homogeneous expansion of `s_\lambda`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import s_to_h
        sage: sorted(s_to_h((2, 1)).items())
        [((2, 1), 1), ((3,), -1)]
    """
    return dict(_s_to_h(_normalize(lam)))


def s_to_e(lam):
    r"""
    Return the elementary expansion of `s_\lambda`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import s_to_e
        sage: sorted(s_to_e((2, 1, 1)).items())
        [((3, 1), 1), ((4,), -1)]
    """
    return dict(_s_to_h(_conjugate(_normalize(lam))))


##############################################################################
# Murnaghan-Nakayama
##############################################################################

cdef object _char(tuple lam, tuple mu):
    # |lam| == |mu|, mu weakly decreasing; remove the largest part first.
    cdef Py_ssize_t k = len(mu)
    if k == 0 or len(lam) == 1:
        return 1
    if lam[0] == 1:
        return -1 if (len(lam) - k) & 1 else 1
    key = (lam, mu)
    c = _char_cache.get(key)
    if c is not None:
        return c
    cdef tuple rest = mu[1:]
    total = 0
    for nu, sign in _border_strips(lam, mu[0], False):
        if sign > 0:
            total += _char(<tuple> nu, rest)
        else:
            total -= _char(<tuple> nu, rest)
    _char_cache[key] = total
    return total


def character_value(lam, mu):
    r"""
    Return `\chi^\lambda(\mu)`, the value of the irreducible character of
    the symmetric group indexed by ``lam`` on the class of cycle type ``mu``.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import character_value
        sage: [character_value((2, 1), mu) for mu in [(1, 1, 1), (2, 1), (3,)]]
        [2, 0, -1]

    TESTS:

    The first orthogonality relation and the degrees::

        sage: from sage.combinat.sf.transition_kernels import z
        sage: P = [tuple(la) for la in Partitions(6)]
        sage: all(sum(character_value(la, mu) * character_value(nu, mu) / z(mu)
        ....:         for mu in P) == (la == nu)
        ....:     for la in P for nu in P)
        True
        sage: all(character_value(la, (1,) * 7) == StandardTableaux(la).cardinality()
        ....:     for la in Partitions(7))
        True
    """
    lam = tuple(int(p) for p in lam if p)
    mu = _normalize(mu)
    if sum(lam) != sum(mu):
        raise ValueError("lam and mu must have the same size")
    return Integer(_char(lam, mu))


cdef dict _p_to_s(tuple mu):
    D = _p_to_s_cache.get(mu)
    if D is not None:
        return D
    cdef Py_ssize_t k = len(mu)
    if k == 0:
        D = {(): 1}
    else:
        D = {}
        r = mu[k - 1]
        for lam, c in _p_to_s(mu[:k - 1]).items():
            for nu, sign in _border_strips(lam, r, True):
                _add_term(D, nu, c if sign > 0 else -c)
    _p_to_s_cache[mu] = D
    return D


def p_to_s(mu):
    r"""
    Return the Schur expansion of `p_\mu`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import p_to_s
        sage: sorted(p_to_s((2, 1)).items())
        [((1, 1, 1), -1), ((3,), 1)]
    """
    return dict(_p_to_s(_normalize(mu)))


def s_to_p(lam):
    r"""
    Return the power sum expansion of `s_\lambda`, with rational
    coefficients `\chi^\lambda(\mu) / z_\mu`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import s_to_p
        sage: sorted(s_to_p((2, 1)).items())
        [((1, 1, 1), 1/3), ((3,), -1/3)]
    """
    lam = _normalize(lam)
    D = {}
    for mu, zmu in _partitions_z(sum(lam)):
        c = _char(lam, <tuple> mu)
        if c:
            D[mu] = Integer(c) / zmu
    return D


##############################################################################
# Multiplicative bases
##############################################################################

cdef enum:
    H_TO_P = 0
    E_TO_P = 1
    P_TO_H = 2
    P_TO_E = 3
    H_TO_E = 4  # also e_k in terms of h, by symmetry


cdef object _mult_factorials(tuple mu):
    r"""
    Return `\prod_i m_i(\mu)!`.
    """
    cdef Py_ssize_t i = 0, j, k = len(mu)
    res = 1
    while i < k:
        j = i
        while j < k and mu[j] == mu[i]:
            j += 1
        res *= factorial(j - i)
        i = j
    return res


cdef dict _generator(int kind, int k):
    r"""
    Return the expansion of the generator of degree ``k``:

    - `h_k = \sum_\mu p_\mu / z_\mu`
    - `e_k = \sum_\mu \epsilon_\mu p_\mu / z_\mu`
    - `p_k = \sum_\mu (-1)^{\ell(\mu)-1} k (\ell(\mu)-1)! / \prod_i m_i! \, h_\mu`
    - `p_k = \sum_\mu \epsilon_\mu k (\ell(\mu)-1)! / \prod_i m_i! \, e_\mu`
    - `h_k = \sum_\mu \epsilon_\mu \ell(\mu)! / \prod_i m_i! \, e_\mu`

    where `\epsilon_\mu = (-1)^{k - \ell(\mu)}`.
    """
    key = (kind, k)
    cdef dict D = _gen_cache.get(key)
    if D is not None:
        return D
    D = {}
    cdef Py_ssize_t l
    cdef int sign
    for mu in _partitions(k):
        l = len(<tuple> mu)
        sign = -1 if (k - l) & 1 else 1
        # the power sum generators are scaled by k! to stay integral
        if kind == H_TO_P:
            c = factorial(k) // z(mu)
        elif kind == E_TO_P:
            c = sign * factorial(k) // z(mu)
        elif kind == P_TO_H:
            c = ((-1 if (l - 1) & 1 else 1) * k * factorial(l - 1)
                 // _mult_factorials(<tuple> mu))
        elif kind == P_TO_E:
            c = sign * k * factorial(l - 1) // _mult_factorials(<tuple> mu)
        else:
            c = sign * factorial(l) // _mult_factorials(<tuple> mu)
        D[mu] = c
    _gen_cache[key] = D
    return D


cdef tuple _union(tuple a, tuple b):
    return tuple(sorted(a + b, reverse=True))


cdef dict _multiplicative(int kind, tuple lam):
    key = (kind, lam)
    cdef dict D = _mult_cache.get(key)
    if D is not None:
        return D
    cdef Py_ssize_t k = len(lam)
    if not k:
        D = {(): 1}
    else:
        D = {}
        g = _generator(kind, lam[k - 1])
        for a, c in _multiplicative(kind, lam[:k - 1]).items():
            for b, d in g.items():
                _add_term(D, _union(<tuple> a, <tuple> b), c * d)
    _mult_cache[key] = D
    return D


def h_to_p(mu):
    r"""
    Return the power sum expansion of `h_\mu`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import h_to_p
        sage: sorted(h_to_p((2,)).items())
        [((1, 1), 1/2), ((2,), 1/2)]
    """
    return _unscale(_multiplicative(H_TO_P, _normalize(mu)), _normalize(mu))


def e_to_p(mu):
    r"""
    Return the power sum expansion of `e_\mu`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import e_to_p
        sage: sorted(e_to_p((2,)).items())
        [((1, 1), 1/2), ((2,), -1/2)]
    """
    return _unscale(_multiplicative(E_TO_P, _normalize(mu)), _normalize(mu))


cdef dict _unscale(dict D, tuple lam):
    r"""
    Divide the coefficients of ``D`` by `\prod_i \lambda_i!`.
    """
    den = Integer(1)
    for p in lam:
        den *= factorial(p)
    return {nu: Integer(c) / den for nu, c in D.items()}


def p_to_h(mu):
    r"""
    Return the complete homogeneous expansion of `p_\mu`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import p_to_h
        sage: sorted(p_to_h((2,)).items())
        [((1, 1), -1), ((2,), 2)]
    """
    return dict(_multiplicative(P_TO_H, _normalize(mu)))


def p_to_e(mu):
    r"""
    Return the elementary expansion of `p_\mu`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import p_to_e
        sage: sorted(p_to_e((2,)).items())
        [((1, 1), 1), ((2,), -2)]
    """
    return dict(_multiplicative(P_TO_E, _normalize(mu)))


def h_to_e(mu):
    r"""
    Return the elementary expansion of `h_\mu`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import h_to_e
        sage: sorted(h_to_e((2,)).items())
        [((1, 1), 1), ((2,), -1)]
    """
    return dict(_multiplicative(H_TO_E, _normalize(mu)))


def e_to_h(mu):
    r"""
    Return the complete homogeneous expansion of `e_\mu`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import e_to_h
        sage: sorted(e_to_h((2, 1)).items())
        [((1, 1, 1), 1), ((2, 1), -1)]
    """
    return dict(_multiplicative(H_TO_E, _normalize(mu)))


##############################################################################
# Power sums and monomials
##############################################################################

cdef list _p_times_m(int k, tuple nu):
    r"""
    Return the pairs ``(rho, c)`` with `p_k m_\nu = \sum c\, m_\rho`; in the
    first pair ``k`` is a new part.

    The coefficient of `m_\rho` is the multiplicity in `\rho` of the part
    that was created or increased.
    """
    cdef Py_ssize_t j, n = len(nu)
    cdef tuple rho = _insert_part(nu, k)
    cdef list out = [(rho, rho.count(k))]
    prev = None
    for j in range(n):
        a = nu[j]
        if a == prev:
            continue
        prev = a
        rho = _insert_part(nu[:j] + nu[j + 1:], a + k)
        out.append((rho, rho.count(a + k)))
    return out


cdef dict _p_to_m(tuple mu):
    cdef dict D = _p_to_m_cache.get(mu)
    if D is not None:
        return D
    cdef Py_ssize_t k = len(mu)
    if not k:
        D = {(): 1}
    else:
        D = {}
        r = mu[k - 1]
        for nu, c in _p_to_m(mu[:k - 1]).items():
            for rho, d in _p_times_m(r, <tuple> nu):
                _add_term(D, rho, c * d)
    _p_to_m_cache[mu] = D
    return D


cdef dict _m_to_p(tuple lam):
    # With k = lam[0] and lam = rest + (k,):
    # mult_k(lam) m_lam = p_k m_rest - (the other terms of p_k m_rest),
    # and the other terms have fewer parts.
    cdef dict D = _m_to_p_cache.get(lam)
    if D is not None:
        return D
    cdef Py_ssize_t idx
    if not lam:
        D = {(): 1}
    else:
        k = lam[0]
        rest = lam[1:]
        D = {}
        for mu, c in _m_to_p(rest).items():
            _add_term(D, _insert_part(<tuple> mu, k), c)
        terms = _p_times_m(k, rest)
        for idx in range(1, len(terms)):
            rho, c = terms[idx]
            for mu, d in _m_to_p(<tuple> rho).items():
                _add_term(D, mu, -c * d)
        mk = Integer(lam.count(k))
        if mk != 1:
            D = {mu: c / mk for mu, c in D.items()}
    _m_to_p_cache[lam] = D
    return D


def p_to_m(mu):
    r"""
    Return the monomial expansion of `p_\mu`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import p_to_m
        sage: sorted(p_to_m((1, 1)).items())
        [((1, 1), 2), ((2,), 1)]
    """
    return dict(_p_to_m(_normalize(mu)))


def m_to_p(lam):
    r"""
    Return the power sum expansion of `m_\lambda`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import m_to_p
        sage: sorted(m_to_p((1, 1)).items())
        [((1, 1), 1/2), ((2,), -1/2)]
    """
    return dict(_m_to_p(_normalize(lam)))


##############################################################################
# Monomial products
##############################################################################

# The coefficient of m_nu in m_lam m_mu counts the pairs (alpha, beta) of
# rearrangements of lam and mu, padded with zeros to length
# L = len(lam) + len(mu), with alpha + beta = nu. Fixing alpha to lam gives
# a count c' with c_nu = c' * |orbit(lam)| / |orbit(nu)|. To compute c', the
# positions with equal values of lam form blocks, and each block receives a
# multiset of values of mu.

cdef int _mm_blocks(list blocks, Py_ssize_t bi, list vals, int* counts,
                    Py_ssize_t nvals, list parts, object weight,
                    dict acc) except -1:
    if bi == len(blocks):
        nu = tuple(sorted([p for p in parts if p], reverse=True))
        acc[nu] = acc.get(nu, 0) + weight
        return 0
    size = (<tuple> blocks[bi])[1]
    return _mm_values(blocks, bi, 0, size, vals, counts, nvals, parts,
                      weight * factorial(size), acc)


cdef int _mm_values(list blocks, Py_ssize_t bi, Py_ssize_t j, int rem,
                    list vals, int* counts, Py_ssize_t nvals, list parts,
                    object weight, dict acc) except -1:
    # Choose how many of the remaining positions of block bi receive the
    # value vals[j]; weight accumulates the number of arrangements.
    cdef int x, top
    cdef Py_ssize_t q
    if rem == 0:
        return _mm_blocks(blocks, bi + 1, vals, counts, nvals, parts,
                          weight, acc)
    if j == nvals:
        return 0
    w = (<tuple> blocks[bi])[0] + vals[j]
    top = counts[j] if counts[j] < rem else rem
    for x in range(top + 1):
        counts[j] -= x
        for q in range(x):
            parts.append(w)
        _mm_values(blocks, bi, j + 1, rem - x, vals, counts, nvals, parts,
                   weight // factorial(x), acc)
        counts[j] += x
        if x:
            del parts[len(parts) - x:]
    return 0


cdef list _blocks(tuple lam):
    """
    Return the pairs ``(value, multiplicity)`` of ``lam``.
    """
    cdef Py_ssize_t i = 0, j, k = len(lam)
    cdef list out = []
    while i < k:
        j = i
        while j < k and lam[j] == lam[i]:
            j += 1
        out.append((lam[i], j - i))
        i = j
    return out


cdef dict _monomial_product(tuple lam, tuple mu):
    if not lam:
        return {mu: 1}
    if not mu:
        return {lam: 1}
    if lam < mu:
        lam, mu = mu, lam
    key = (lam, mu)
    cdef dict D = _mm_cache.get(key)
    if D is not None:
        return D
    cdef Py_ssize_t ll = len(lam), lm = len(mu), L = ll + lm, i
    cdef list blocks = _blocks(lam) + [(0, lm)]
    cdef list mblocks = _blocks(mu) + [(0, ll)]
    cdef Py_ssize_t nvals = len(mblocks)
    cdef list vals = [b[0] for b in mblocks]
    cdef dict acc = {}
    cdef int* counts = <int*> check_allocarray(nvals, sizeof(int))
    try:
        for i in range(nvals):
            counts[i] = mblocks[i][1]
        _mm_blocks(blocks, 0, vals, counts, nvals, [], 1, acc)
    finally:
        sig_free(counts)
    denominator = _mult_factorials(lam) * factorial(lm)
    D = {nu: (c * _mult_factorials(<tuple> nu) * factorial(L - len(nu))
              // denominator)
         for nu, c in acc.items()}
    _mm_cache[key] = D
    return D


def monomial_product(lam, mu):
    r"""
    Return the monomial expansion of `m_\lambda m_\mu`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import monomial_product
        sage: sorted(monomial_product((1,), (1,)).items())
        [((1, 1), 2), ((2,), 1)]
        sage: sorted(monomial_product((2, 1), (2, 1)).items())
        [((2, 2, 1, 1), 4), ((2, 2, 2), 6), ((3, 2, 1), 2), ((3, 3), 2),
         ((4, 1, 1), 2), ((4, 2), 1)]
        sage: monomial_product((), (2, 1))
        {(2, 1): 1}
    """
    return dict(_monomial_product(_normalize(lam), _normalize(mu)))


cdef inline tuple _exact_tuple(int* a, Py_ssize_t k):
    cdef tuple t = PyTuple_New(k)
    cdef Py_ssize_t i
    cdef object x
    for i in range(k):
        x = a[i]
        Py_INCREF(x)
        PyTuple_SET_ITEM(t, i, x)
    return t


def monomial_exponents(mu, int n):
    r"""
    Return the exponent vectors of the monomials of `m_\mu(x_1, \ldots, x_n)`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import monomial_exponents
        sage: sorted(monomial_exponents((2, 1), 3))
        [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
        sage: monomial_exponents((1, 1, 1), 2)
        []
        sage: monomial_exponents((), 2)
        [(0, 0)]
    """
    cdef tuple m = _normalize(mu)
    cdef Py_ssize_t k = len(m), i, j
    cdef int tmp
    cdef list out = []
    if k > n:
        return out
    if n == 0:
        return [()]
    cdef int* a = <int*> check_calloc(n, sizeof(int))
    try:
        # start from the smallest arrangement and step through the others
        # in lexicographic order
        for i in range(k):
            a[n - 1 - i] = m[i]
        while True:
            out.append(_exact_tuple(a, n))
            i = n - 2
            while i >= 0 and a[i] >= a[i + 1]:
                i -= 1
            if i < 0:
                break
            j = n - 1
            while a[j] <= a[i]:
                j -= 1
            tmp = a[i]
            a[i] = a[j]
            a[j] = tmp
            i += 1
            j = n - 1
            while i < j:
                tmp = a[i]
                a[i] = a[j]
                a[j] = tmp
                i += 1
                j -= 1
    finally:
        sig_free(a)
    return out


##############################################################################
# Hall-Littlewood
##############################################################################

cdef object _polynomial_ring():
    global _ZZt
    if _ZZt is None:
        from sage.rings.integer_ring import ZZ
        _ZZt = ZZ['t']
    return _ZZt


cdef dict _jing(int m, tuple nu):
    # B_m(s_nu) = sum_{i,j} (-1)^i t^j h_{m+i+j} e_i^perp h_j^perp s_nu
    key = (m, nu)
    cdef dict D = _jing_cache.get(key)
    if D is not None:
        return D
    cdef int d = sum(nu), i, j, sign
    cdef dict terms = {}
    cdef dict poly
    for j in range(d + 1):
        for kappa in _hstrips(nu, j, False):
            kc = _conjugate(<tuple> kappa)
            for i in range(d - j + 1):
                sign = -1 if i & 1 else 1
                for rhoc in _hstrips(kc, i, False):
                    rho = _conjugate(<tuple> rhoc)
                    for lam in _hstrips(rho, m + i + j, True):
                        poly = terms.get(lam)
                        if poly is None:
                            poly = terms[lam] = {}
                        poly[j] = poly.get(j, 0) + sign
    R = _polynomial_ring()
    D = {}
    for lam, poly in terms.items():
        p = R(poly)
        if p:
            D[lam] = p
    _jing_cache[key] = D
    return D


cdef dict _qp_to_s(tuple mu):
    cdef dict D = _qp_cache.get(mu)
    if D is not None:
        return D
    if not mu:
        D = {(): _polynomial_ring().one()}
    else:
        D = {}
        for nu, c in _qp_to_s(mu[1:]).items():
            for lam, p in _jing(mu[0], <tuple> nu).items():
                _add_term(D, lam, c * p)
    _qp_cache[mu] = D
    return D


def hall_littlewood_qp_to_s(mu):
    r"""
    Return the Schur expansion of the modified Hall-Littlewood function
    `Q'_\mu`, whose coefficients are the Kostka-Foulkes polynomials
    `K_{\lambda\mu}(t)`.

    This uses Jing's vertex operators:
    `Q'_\mu = B_{\mu_1} \cdots B_{\mu_\ell}(1)` where
    `B_m = \sum_{i,j \geq 0} (-1)^i t^j h_{m+i+j} e_i^\perp h_j^\perp`.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import hall_littlewood_qp_to_s
        sage: sorted(hall_littlewood_qp_to_s((2, 1)).items())
        [((2, 1), 1), ((3,), t)]

    TESTS::

        sage: from sage.combinat.sf.kfpoly import kfpoly
        sage: all(hall_littlewood_qp_to_s(mu).get(tuple(la), 0) == kfpoly(la, mu)
        ....:     for n in range(1, 7) for mu in Partitions(n) for la in Partitions(n))
        True
    """
    return dict(_qp_to_s(_normalize(mu)))


##############################################################################
# Semistandard tableaux
##############################################################################

cdef int _ssyt_chains(tuple lam, list content, list rests, Py_ssize_t k,
                      list chain, list out) except -1:
    # Remove the cells containing the letter k, a horizontal strip. The
    # remaining shape must dominate the remaining content, so that it can
    # still be filled.
    if k == 0:
        out.append(tuple(chain))
        return 0
    cdef tuple rest = rests[k - 1]
    for nu in _hstrips(lam, content[k - 1], False):
        if _dominates(<tuple> nu, rest):
            chain[k - 1] = nu
            _ssyt_chains(<tuple> nu, content, rests, k - 1, chain, out)
    return 0


def semistandard_tableaux(shape, content):
    r"""
    Return the semistandard tableaux of shape ``shape`` and content
    ``content``, as lists of rows.

    The content may contain zeros. The tableaux are found by removing the
    horizontal strips formed by the largest entries, and are sorted
    lexicographically by their rows, which is the order used by Symmetrica.

    EXAMPLES::

        sage: from sage.combinat.sf.transition_kernels import kostka_number, semistandard_tableaux
        sage: semistandard_tableaux((3, 2, 1), (2, 2, 2))
        [[[1, 1, 2], [2, 3], [3]], [[1, 1, 3], [2, 2], [3]]]
        sage: semistandard_tableaux((3, 1), (1, 0, 2, 1))
        [[[1, 3, 3], [4]], [[1, 3, 4], [3]]]
        sage: semistandard_tableaux((2, 2), (3, 1))
        []
        sage: semistandard_tableaux((), ())
        [[]]

    TESTS::

        sage: all(len(semistandard_tableaux(la, mu)) == kostka_number(la, mu)
        ....:     for n in range(7) for la in Partitions(n) for mu in Compositions(n))
        True
    """
    cdef tuple lam = tuple(int(p) for p in shape if p)
    cdef list cont = [int(c) for c in content]
    cdef Py_ssize_t K = len(cont), i, r, c
    if sum(lam) != sum(cont):
        return []
    cdef list rests = [_normalize(cont[:i]) for i in range(K)]
    cdef list out = []
    _ssyt_chains(lam, cont, rests, K, [None] * K, out)

    cdef list res = []
    cdef list rows
    cdef tuple prev, cur
    for chain in out:
        shapes = chain + (lam,)
        rows = [[] for _ in lam]
        for i in range(1, K + 1):
            prev = shapes[i - 1]
            cur = shapes[i]
            letter = Integer(i)
            for r in range(len(cur)):
                for c in range(prev[r] if r < len(prev) else 0, cur[r]):
                    rows[r].append(letter)
        res.append(rows)
    res.sort()
    return res
