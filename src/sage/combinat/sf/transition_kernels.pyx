r"""
Kostka and Murnaghan-Nakayama kernels for the classical bases

This module computes transitions between the Schur, complete homogeneous,
elementary, monomial and power sum bases one basis element at a time.
Partitions are plain tuples of Python integers and the results are
dictionaries mapping partitions to coefficients.

The two combinatorial kernels are:

- Kostka numbers `K_{\lambda\mu}`, computed by adding (Pieri rule) or
  removing horizontal strips;

- irreducible character values `\chi^\lambda(\mu)`, computed with the
  Murnaghan-Nakayama rule by adding or removing border strips, which are
  located using beta-numbers.

The inverse Kostka transitions are obtained by unitriangularity in
dominance order.

Intermediate results are memoized; call :func:`clear_caches` to free them.
The dictionaries returned by the public functions are fresh copies.

.. WARNING::

    This is a prototype; the arguments are not validated.

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
              _s_to_m_cache, _m_to_s_cache, _s_to_h_cache, _partitions_cache):
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
    """
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
