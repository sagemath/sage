# distutils: libraries = gmp
# distutils: extra_compile_args = -D_XPG6

r"""
Lists of Manin symbols over `\QQ`, elements of `\mathbb{P}^1(\ZZ/N\ZZ)`
"""

from cysignals.memory cimport check_allocarray, sig_free
from cysignals.signals cimport sig_check

from sage.misc.search cimport search
from sage.structure.richcmp cimport rich_to_bool

cimport sage.rings.fast_arith
import sage.rings.fast_arith
cdef sage.rings.fast_arith.arith_int arith_int
cdef sage.rings.fast_arith.arith_llong arith_llong
arith_int = sage.rings.fast_arith.arith_int()
arith_llong = sage.rings.fast_arith.arith_llong()

ctypedef long long llong

###############################################################
#
# Int p1_normalize and p1list; a non-int version is below,
# which should be used for N > 46340
#
################################################################

cdef int c_p1_normalize_int(int N, int u, int v,
                            int* uu, int* vv, int* ss,
                            int compute_s) except -1:
    r"""
    Compute the canonical representative of
    `\mathbb{P}^1(\ZZ/N\ZZ)` equivalent to
    `(u,v)` along with a transforming scalar.

    INPUT:

    - ``N`` -- integer

    - ``u`` -- integer

    - ``v`` -- integer

    OUTPUT:

    If `\gcd(u,v,N) = 1`, then returns

    - ``uu`` -- integer

    - ``vv`` -- integer

    - ``ss`` -- integer such that `(ss*uu, ss*vv)` is congruent to `(u,v)` (mod `N`);

    If `\gcd(u,v,N) \not= 1`, returns ``0, 0, 0``.

    If ``compute_s`` is 0, ``s`` is not computed.
    """
    cdef int d, k, g, s, t, min_v, min_t, Ng, vNg
    if N == 1:
        uu[0] = 0
        vv[0] = 0
        ss[0] = 1
        return 0

    if N <= 0 or 46340 < N:
        raise OverflowError("Modulus is too large (must be <= 46340)")

    u = u % N
    v = v % N
    if u < 0:
        u += N
    if v < 0:
        v += N
    if u == 0:
        uu[0] = 0
        if arith_int.c_gcd_int(v,N) == 1:
            vv[0] = 1
        else:
            vv[0] = 0
        ss[0] = v
        return 0

    g = arith_int.c_xgcd_int(u, N, &s, &t)
    s = s % N
    if s < 0:
        s += N
    if arith_int.c_gcd_int(g, v) != 1:
        uu[0] = 0
        vv[0] = 0
        ss[0] = 0
        return 0

    # Now g = s*u + t*N, so s is a "pseudo-inverse" of u mod N
    # Adjust s modulo N/g so it is coprime to N.
    if g!=1:
        d = N/g
        while arith_int.c_gcd_int(s,N) != 1:
            s = (s+d) % N

    # Multiply [u,v] by s; then [s*u,s*v] = [g,s*v] (mod N)
    u = g
    v = (s*v) % N

    min_v = v
    min_t = 1
    if g!=1:
        Ng = N/g
        vNg = (v*Ng) % N
        t = 1
        for k in range(2, g + 1):
            v = (v + vNg) % N
            t = (t + Ng) % N
            if v<min_v and arith_int.c_gcd_int(t,N)==1:
                min_v = v
                min_t = t
    v = min_v
    if u < 0:
        u += N
    if v < 0:
        v += N
    uu[0] = u
    vv[0] = v
    if compute_s:
        ss[0] = arith_int.c_inverse_mod_int(s*min_t, N)
    return 0


def p1_normalize_int(N, u, v):
    r"""
    Compute the canonical representative of
    `\mathbb{P}^1(\ZZ/N\ZZ)` equivalent to
    `(u,v)` along with a transforming scalar.

    INPUT:

    - ``N`` -- integer

    - ``u`` -- integer

    - ``v`` -- integer

    OUTPUT:

    If `\gcd(u,v,N) = 1`, then returns

    - ``uu`` -- integer

    - ``vv`` -- integer

    - ``ss`` -- integer such that `(ss*uu, ss*vv)` is congruent to `(u,v)` (mod `N`);

    If `\gcd(u,v,N) \not= 1`, returns ``0, 0, 0``.

    EXAMPLES::

        sage: from sage.modular.modsym.p1list import p1_normalize_int
        sage: p1_normalize_int(90,7,77)
        (1, 11, 7)
        sage: p1_normalize_int(90,7,78)
        (1, 24, 7)
        sage: (7*24-78*1) % 90
        0
        sage: (7*24) % 90
        78
    """
    cdef int uu, vv, ss
    c_p1_normalize_int(N, u % N, v % N, &uu, &vv, &ss, 1)
    return (uu, vv, ss)


def p1list_int(int N):
    r"""
    Return a list of the normalized elements of
    `\mathbb{P}^1(\ZZ/N\ZZ)`.

    INPUT:

    - ``N`` -- integer (the level or modulus)

    EXAMPLES::

        sage: from sage.modular.modsym.p1list import p1list_int
        sage: p1list_int(6)
        [(0, 1),
        (1, 0),
        (1, 1),
        (1, 2),
        (1, 3),
        (1, 4),
        (1, 5),
        (2, 1),
        (2, 3),
        (2, 5),
        (3, 1),
        (3, 2)]

    ::

        sage: p1list_int(120)
        [(0, 1),
        (1, 0),
        (1, 1),
        (1, 2),
        (1, 3),
        ...
        (30, 7),
        (40, 1),
        (40, 3),
        (40, 11),
        (60, 1)]
    """
    cdef int g, u, v, s, c, d, h, d1, cmax
    cdef object lst

    if N == 1:
        return [(0, 0)]

    lst = [(0,1)]
    c = 1
    for d in range(N):
        lst.append((c, d))

    cmax = N // 2
    if N % 2:   # N odd, max divisor is <= N/3
        if N % 3:  # N not a multiple of 3 either, max is N/5
            cmax = N // 5
        else:
            cmax = N // 3

    for c in range(2, cmax + 1):
        if N % c == 0:  # c is a proper divisor
            h = N // c
            g = arith_int.c_gcd_int(c, h)
            for d in range(1, h + 1):
                sig_check()
                if arith_int.c_gcd_int(d, g) == 1:
                    d1 = d
                    while arith_int.c_gcd_int(d1, c) != 1:
                        d1 += h
                    c_p1_normalize_int(N, c, d1, &u, &v, &s, 0)
                    lst.append((u,v))
    lst.sort()
    return lst


###############################################################
#
# The following is a version of the three functions
#
#      c_p1_normalize_int, p1_normalize_int, and p1list_int
#
# but the gcd's are done using GMP, so there are no overflow
# worries.   Also, the one multiplication is done using double
# precision.
#
################################################################

cdef int c_p1_normalize_llong(int N, int u, int v,
                              int* uu, int* vv, int* ss,
                              int compute_s) except -1:
    r"""
    c_p1_normalize_llong(N, u, v):

    Computes the canonical representative of
    `\mathbb{P}^1(\ZZ/N\ZZ)` equivalent to `(u,v)` along
    with a transforming scalar 's' (if compute_s is 1).

    INPUT:

    - ``N`` -- integer (the modulus or level)

    - ``u`` -- integer (the first coordinate of (u:v))

    - ``v`` -- integer (the second coordinate of (u:v))

    - ``compute_s`` -- boolean (int)

    OUTPUT:

    If `\gcd(u,v,N) = 1`, then returns

    - ``uu`` -- integer

    - ``vv`` -- integer

    - ``ss`` -- integer such that `(ss*uu, ss*vv)` is equivalent to `(u,v)` mod `N`

    If `\gcd(u,v,N) \not= 1`, returns ``0, 0, 0``.

    EXAMPLES::

        sage: from sage.modular.modsym.p1list import p1_normalize_int
        sage: p1_normalize_int(90,7,77)
        (1, 11, 7)
        sage: p1_normalize_int(90,7,78)
        (1, 24, 7)
        sage: (7*24-78*1) % 90
        0
        sage: (7*24) % 90
        78

    TESTS:

    This test reflects :issue:`20932`::

        sage: N = 3*61379
        sage: import sage.modular.modsym.p1list as p1list
        sage: p1 = p1list.P1List(N) # not tested -- too long
        sage: p1.normalize_with_scalar(21, -1) # not tested -- too long
        (3, 105221, 7)
    """
    cdef int d, k, g, s, t, min_v, min_t, Ng, vNg
    cdef llong ll_s, ll_t, ll_N
    if N == 1:
        uu[0] = 0
        vv[0] = 0
        ss[0] = 1
        return 0

    ll_N = <long> N

    #if N<=0 or N >= 2**31:
    #    raise OverflowError("Modulus is too large (must be < 46340)")
    #    return -1

    u = u % N
    v = v % N
    if u < 0:
        u += N
    if v < 0:
        v += N
    if u == 0:
        uu[0] = 0
        if arith_int.c_gcd_int(v,N) == 1:
            vv[0] = 1
        else:
            vv[0] = 0
        ss[0] = v
        return 0

    #g = xgcd_int_llong(u, N, &s, &t)
    g = <int> arith_llong.c_xgcd_longlong(u, N, &ll_s, &ll_t)
    s = <int>(ll_s % ll_N)
    t = <int>(ll_t % ll_N)
    s = s % N
    if s < 0:
        s += N
    if arith_int.c_gcd_int(g, v) != 1:
        uu[0] = 0
        vv[0] = 0
        ss[0] = 0
        return 0

    # Now g = s*u + t*N, so s is a "pseudo-inverse" of u mod N
    # Adjust s modulo N/g so it is coprime to N.
    if g!=1:
        d = N/g
        while arith_int.c_gcd_int(s,N) != 1:
            s = (s+d) % N

    # Multiply [u,v] by s; then [s*u,s*v] = [g,s*v] (mod N)
    u = g
    # v = (s*v) % N
    v = <int> (((<llong>s) * (<llong>v)) % ll_N)

    min_v = v
    min_t = 1
    if g!=1:
        Ng = N/g
        vNg = <int> ((<llong>v * <llong> Ng) % ll_N)
        t = 1
        for k in range(2, g + 1):
            v = (v + vNg) % N
            t = (t + Ng) % N
            if v < min_v and arith_int.c_gcd_int(t, N) == 1:
                min_v = v
                min_t = t
    v = min_v
    if u < 0:
        u += N
    if v < 0:
        v += N
    uu[0] = u
    vv[0] = v
    if compute_s:
        ss[0] = <int> (arith_llong.c_inverse_mod_longlong((<llong> s)*(<llong> min_t), N) % ll_N)
    return 0


def p1_normalize_llong(N, u, v):
    r"""
    Compute the canonical representative of
    `\mathbb{P}^1(\ZZ/N\ZZ)` equivalent to
    `(u,v)` along with a transforming scalar.

    INPUT:

    - ``N`` -- integer

    - ``u`` -- integer

    - ``v`` -- integer

    OUTPUT:

    If `\gcd(u,v,N) = 1`, then returns

    - ``uu`` -- integer

    - ``vv`` -- integer

    - ``ss`` -- integer such that `(ss*uu, ss*vv)` is equivalent to `(u,v)` mod `N`

    If `\gcd(u,v,N) \not= 1`, returns ``0, 0, 0``.

    EXAMPLES::

        sage: from sage.modular.modsym.p1list import p1_normalize_llong
        sage: p1_normalize_llong(90000,7,77)
        (1, 11, 7)
        sage: p1_normalize_llong(90000,7,78)
        (1, 77154, 7)
        sage: (7*77154-78*1) % 90000
        0
        sage: (7*77154) % 90000
        78
    """
    cdef int uu, vv, ss
    c_p1_normalize_llong(N, u % N, v % N, &uu, &vv, &ss, 1)
    return (uu, vv, ss)


def p1list_llong(int N):
    r"""
    Return a list of the normalized elements of
    `\mathbb{P}^1(\ZZ/N\ZZ)`, as a plain list of
    2-tuples.

    INPUT:

    - ``N`` -- integer (the level or modulus)

    EXAMPLES::

        sage: from sage.modular.modsym.p1list import p1list_llong
        sage: N = 50000
        sage: L = p1list_llong(50000)
        sage: len(L) == N*prod([1+1/p for p,e in N.factor()])
        True
        sage: L[0]
        (0, 1)
        sage: L[len(L)-1]
        (25000, 1)

    TESTS:

    This test shows that :issue:`20932` has been resolved::

        sage: import sage.modular.modsym.p1list as p1list
        sage: [(i,j) for (i,j) in p1list.P1List(103809) if i != 1 and i != 3] # not tested -- too long
        [(0, 1), (34603, 1), (34603, 2), (34603, 3)]
    """
    cdef int g, u, v, s, c, d, h, d1, cmax
    if N == 1:
        return [(0, 0)]

    lst = [(0,1)]
    c = 1
    for d in range(N):
        lst.append((c, d))

    cmax = N // 2
    if N % 2:   # N odd, max divisor is <= N/3
        if N % 3:  # N not a multiple of 3 either, max is N/5
            cmax = N // 5
        else:
            cmax = N // 3

    for c in range(2, cmax + 1):
        if N % c == 0:  # c is a proper divisor
            h = N // c
            g = arith_int.c_gcd_int(c, h)
            for d in range(1, h + 1):
                if arith_int.c_gcd_int(d, g) == 1:
                    sig_check()
                    d1 = d
                    while arith_int.c_gcd_int(d1, c) != 1:
                        d1 += h
                    c_p1_normalize_llong(N, c, d1, &u, &v, &s, 0)
                    lst.append((u, v))
    lst.sort()
    return lst


def p1list(N):
    r"""
    Return the elements of the projective line modulo `N`,
    `\mathbb{P}^1(\ZZ/N\ZZ)`, as a plain list of 2-tuples.

    INPUT:

    - ``N`` -- integer; a positive integer (less than 2^31)

    OUTPUT:

    A list of the elements of the projective line
    `\mathbb{P}^1(\ZZ/N\ZZ)`, as plain 2-tuples.

    EXAMPLES::

        sage: from sage.modular.modsym.p1list import p1list
        sage: list(p1list(7))
        [(0, 1), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6)]
        sage: N = 23456
        sage: len(p1list(N)) == N*prod([1+1/p for p,e in N.factor()])
        True
    """
    if N <= 0:
        raise ValueError("N must be a positive integer")
    if N <= 46340:
        return p1list_int(N)
    if N <= 2147483647:
        return p1list_llong(N)
    raise OverflowError("p1list not defined for such large N")


def p1_normalize(int N, int u, int v):
    r"""
    Compute the canonical representative of
    `\mathbb{P}^1(\ZZ/N\ZZ)` equivalent to
    `(u,v)` along with a transforming scalar.

    INPUT:

    - ``N`` -- integer

    - ``u`` -- integer

    - ``v`` -- integer

    OUTPUT:

    If `\gcd(u,v,N) = 1`, then returns

    - ``uu`` -- integer

    - ``vv`` -- integer

    - ``ss`` -- integer such that `(ss*uu, ss*vv)` is equivalent to `(u,v)` mod `N`

    If `\gcd(u,v,N) \not= 1`, returns ``0, 0, 0``.

    EXAMPLES::

        sage: from sage.modular.modsym.p1list import p1_normalize
        sage: p1_normalize(90,7,77)
        (1, 11, 7)
        sage: p1_normalize(90,7,78)
        (1, 24, 7)
        sage: (7*24-78*1) % 90
        0
        sage: (7*24) % 90
        78

    ::

        sage: from sage.modular.modsym.p1list import p1_normalize
        sage: p1_normalize(50001,12345,54322)
        (3, 4667, 4115)
        sage: (12345*4667-54321*3) % 50001
        3
        sage: 4115*3 % 50001
        12345
        sage: 4115*4667 % 50001 == 54322 % 50001
        True
    """
    if N <= 46340:
        return p1_normalize_int(N, u, v)
    if N <= 2147483647:
        return p1_normalize_llong(N, u, v)
    else:
        raise OverflowError


cdef int p1_normalize_xgcdtable(int N, int u, int v,
                                int compute_s,
                                int *t_g, int *t_a, int *t_b,
                                int* uu, int* vv, int* ss) except -1:
    """
    INPUT:

    - ``N``, ``u``, ``v`` -- integers

    - ``compute_s`` -- do not compute s if ``compute_s == 0``

    - ``t_g``, ``t_a``, ``t_b`` -- integer arrays

    OUTPUT: ``uu, vv, ss`` -- reduced representative and normalizing
    scalar
    """
    cdef int d, k, g, s, t, min_v, min_t, Ng, vNg
    if N == 1:
        uu[0] = 0
        vv[0] = 0
        ss[0] = 1
        return 0

    if N <= 0 or 46340 < N:
        raise OverflowError("Modulus is too large (must be <= 46340)")

    u = u % N
    v = v % N
    if u < 0:
        u += N
    if v < 0:
        v += N
    if u == 0:
        uu[0] = 0
        if t_g[v] == 1:    # "if arith_int.c_gcd_int(v,N) == 1"
            vv[0] = 1
        else:
            vv[0] = 0
        ss[0] = v
        return 0

    #  WAS: "g = arith_int.c_xgcd_int(u, N, &s, &t)"
    g = t_g[u]
    s = t_a[u]
    t = t_b[u]
    s = s % N
    if s < 0:
        s += N
    if g != 1 and arith_int.c_gcd_int(g, v) != 1:
        uu[0] = 0
        vv[0] = 0
        ss[0] = 0
        return 0

    # Now g = s*u + t*N, so s is a "pseudo-inverse" of u mod N
    # Adjust s modulo N/g so it is coprime to N.
    if g!=1:
        d = N/g
        while t_g[s] != 1:  # while arith_int.c_gcd_int(s,N) != 1:
            s = (s+d) % N

    # Multiply [u,v] by s; then [s*u,s*v] = [g,s*v] (mod N)
    u = g
    v = (s*v) % N

    min_v = v
    min_t = 1
    if g!=1:
        Ng = N/g
        vNg = (v*Ng) % N
        t = 1
        for k in range(2, g + 1):
            v = (v + vNg) % N
            t = (t + Ng) % N
            if v < min_v and t_g[t] == 1:     # arith_int.c_gcd_int(t,N)==1:
                min_v = v
                min_t = t
    v = min_v
    if u < 0:
        u += N
    if v < 0:
        v += N
    uu[0] = u
    vv[0] = v
    if compute_s:
        # ss[0] = arith_int.c_inverse_mod_int(s*min_t, N);
        ss[0] = t_a[(s*min_t) % N]
    return 0


cdef class P1List():
    r"""
    The class for `\mathbb{P}^1(\ZZ/N\ZZ)`, the projective line modulo `N`.

    EXAMPLES::

        sage: P = P1List(12); P
        The projective line over the integers modulo 12
        sage: list(P)
        [(0, 1), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (1, 10), (1, 11), (2, 1), (2, 3), (2, 5), (3, 1), (3, 2), (3, 4), (3, 7), (4, 1), (4, 3), (4, 5), (6, 1)]

    Saving and loading works.

    ::

        sage: loads(dumps(P)) == P
        True
    """
    def __init__(self, int N):
        r"""
        The constructor for the class P1List.

        INPUT:

        - ``N`` -- positive integer (the modulus or level)

        OUTPUT: a P1List object representing `\mathbb{P}^1(\ZZ/N\ZZ)`

        EXAMPLES::

            sage: L = P1List(120) # indirect doctest
            sage: L
            The projective line over the integers modulo 120
        """
        cdef int i

        self.__N = N
        if N <= 46340:
            self.__list = p1list_int(N)
            self._normalize = c_p1_normalize_int
        elif N <= 2147483647:
            self.__list = p1list_llong(N)
            self._normalize = c_p1_normalize_llong
        else:
            raise OverflowError("p1list not defined for such large N.")
        self.__list.sort()
        self.__end_hash = dict([(x,i) for i, x in enumerate(self.__list[N+1:])])

        # Allocate memory for xgcd table.
        self.g = NULL
        self.s = NULL
        self.t = NULL
        self.g = <int*>check_allocarray(N, sizeof(int))
        self.s = <int*>check_allocarray(N, sizeof(int))
        self.t = <int*>check_allocarray(N, sizeof(int))

        # Initialize xgcd table
        cdef llong ll_s, ll_t, ll_N = N

        if N <= 46340:
            for i in range(N):
                self.g[i] = arith_int.c_xgcd_int(i, N, &self.s[i], &self.t[i])
        else:
            for i in range(N):
                self.g[i] = arith_llong.c_xgcd_longlong(i, N, &ll_s, &ll_t)
                self.s[i] = <int>(ll_s % ll_N)
                self.t[i] = <int>(ll_t % ll_N)

    def __dealloc__(self):
        """
        Deallocate memory for an object of the class P1List.
        """
        sig_free(self.g)
        sig_free(self.s)
        sig_free(self.t)

    def __richcmp__(self, other, int op):
        """
        Comparison function for objects of the class P1List.

        The order is the same as for the underlying modulus.

        EXAMPLES::

            sage: P1List(10) < P1List(11)
            True
            sage: P1List(12) > P1List(11)
            True
            sage: P1List(11) == P1List(11)
            True

            sage: t = [P1List(N) for N in [100,23,45]]
            sage: [L.N() for L in t]
            [100, 23, 45]
            sage: t.sort()
            sage: [L.N() for L in t]
            [23, 45, 100]
        """
        if not isinstance(other, P1List):
            return NotImplemented
        cdef P1List S = <P1List> self
        cdef P1List O = <P1List> other
        if S.__N < O.__N:
            return rich_to_bool(op, -1)
        elif S.__N > O.__N:
            return rich_to_bool(op, 1)
        return rich_to_bool(op, 0)

    def __reduce__(self):
        """
        Helper function used in pickling.

        EXAMPLES::

            sage: L = P1List(8)
            sage: L.__reduce__()
            (<... 'sage.modular.modsym.p1list.P1List'>, (8,))
        """
        return type(self), (self.__N, )

    def __getitem__(self, n):
        r"""
        Standard indexing/slicing function for the class ``P1List``.

        EXAMPLES::

            sage: L = P1List(8)
            sage: list(L) # indirect doctest
            [(0, 1),
            (1, 0),
            ...
            (2, 3),
            (4, 1)]
            sage: L[4:8] # indirect doctest
            [(1, 3), (1, 4), (1, 5), (1, 6)]
        """
        if isinstance(n, slice):
            start, stop, step = n.indices(len(self))
            return self.__list[start:stop:step]
        else:
            return self.__list[n]

    def __len__(self):
        """
        Return the length of this P1List.

        EXAMPLES::

            sage: L = P1List(8)
            sage: len(L) # indirect doctest
            12
        """
        return len(self.__list)

    def __repr__(self):
        """
        Return the string representation of this P1List.

        EXAMPLES::

            sage: L = P1List(8)
            sage: str(L)            # indirect doctest
            'The projective line over the integers modulo 8'
        """
        return "The projective line over the integers modulo %s" % self.__N

    def lift_to_sl2z(self, int i):
        r"""
        Lift the `i`-th element of this P1list to an element of
        `SL(2,\ZZ)`.

        If the `i`-th element is `(c,d)`, this function computes and
        returns a list `[a,b, c',d']` that defines a 2x2 matrix
        with determinant 1 and integer entries, such that `c=c'` (mod
        `N`) and `d=d'` (mod `N`).

        INPUT:

        - ``i`` -- integer (the index of the element to lift)

        EXAMPLES::

            sage: p = P1List(11)
            sage: p.list()[3]
            (1, 2)

            sage: p.lift_to_sl2z(3)
            [0, -1, 1, 2]

        AUTHORS:

        - Justin Walker
        """
        cdef int c, d, N

        if self.__N == 1:
            return [1,0,0,1]

        c, d = self.__list[i]
        N = self.__N
        # No overflow: This was adjudicated during init...
        if N <= 46340:
            return lift_to_sl2z_int(c, d, self.__N)
        elif N <= 2147483647:
            return lift_to_sl2z_llong(c, d, self.__N)
        else:
            raise OverflowError("N too large")

    def apply_I(self, int i):
        r"""
        Return the index of the result of applying the matrix
        `I=[-1,0;0,1]` to the `i`-th element of this P1List.

        INPUT:

        - ``i`` -- integer (the index of the element to act on)

        EXAMPLES::

            sage: L = P1List(120)
            sage: L[10]
            (1, 9)
            sage: L.apply_I(10)
            112
            sage: L[112]
            (1, 111)
            sage: L.normalize(-1,9)
            (1, 111)

        This operation is an involution::

            sage: all(L.apply_I(L.apply_I(i)) == i for i in range(len(L)))
            True
        """
        cdef int u, v, uu, vv, ss
        u,v = self.__list[i]
        self._normalize(self.__N, -u, v, &uu, &vv, &ss, 0)
        _, j = search(self.__list, (uu,vv))
        return j

    def apply_S(self, int i):
        r"""
        Return the index of the result of applying the matrix
        `S=[0,-1;1,0]` to the `i`-th element of this P1List.

        INPUT:

        - ``i`` -- integer (the index of the element to act on)

        EXAMPLES::

            sage: L = P1List(120)
            sage: L[10]
            (1, 9)
            sage: L.apply_S(10)
            159
            sage: L[159]
            (3, 13)
            sage: L.normalize(-9,1)
            (3, 13)

        This operation is an involution::

            sage: all(L.apply_S(L.apply_S(i)) == i for i in range(len(L)))
            True
        """
        cdef int u, v, uu, vv, ss
        u,v = self.__list[i]
        self._normalize(self.__N, -v, u, &uu, &vv, &ss, 0)
        _, j = search(self.__list, (uu,vv))
        return j

    def apply_T(self, int i):
        r"""
        Return the index of the result of applying the matrix
        `T=[0,1;-1,-1]` to the `i`-th element of this P1List.

        INPUT:

        - ``i`` -- integer (the index of the element to act on)

        EXAMPLES::

            sage: L = P1List(120)
            sage: L[10]
            (1, 9)
            sage: L.apply_T(10)
            157
            sage: L[157]
            (3, 10)
            sage: L.normalize(9,-10)
            (3, 10)

        This operation has order three::

            sage: all(L.apply_T(L.apply_T(L.apply_T(i))) == i for i in range(len(L)))
            True
        """
        cdef int u, v, uu, vv, ss
        u,v = self.__list[i]
        self._normalize(self.__N, v, -u-v, &uu, &vv, &ss, 0)
        _, j = search(self.__list, (uu,vv))
        return j

    cpdef index(self, int u, int v):
        r"""
        Return the index of the class of `(u,v)` in the fixed list
        of representatives of
        `\mathbb{P}^1(\ZZ/N\ZZ)`.

        INPUT:

        - ``u``, ``v`` -- integers with `\gcd(u,v,N)=1`

        OUTPUT:

        - ``i`` -- the index of `u`, `v`, in the P1list

        EXAMPLES::

            sage: L = P1List(120)
            sage: L[100]
            (1, 99)
            sage: L.index(1,99)
            100
            sage: all(L.index(L[i][0],L[i][1])==i for i in range(len(L)))
            True
        """
        if self.__N == 1:
            # there is exactly 1 class [(0,0)].
            return 0
        cdef int uu, vv, ss
        p1_normalize_xgcdtable(self.__N, u, v, 0, self.g, self.s, self.t, &uu, &vv, &ss)
        if uu == 1:
            return vv + 1
        elif uu == 0:
            if vv == 0:
                return -1
            return 0
        try:
            return self.__end_hash[(uu,vv)] + self.__N + 1
        except KeyError:
            return -1

    cdef index_and_scalar(self, int u, int v, int* i, int* s):
        r"""
        Compute the index of the class of `(u,v)` in the fixed list
        of representatives of `\mathbb{P}^1(\ZZ/N\ZZ)` and scalar s
        such that (s*a,s*b) = (u,v), with (a,b) the normalized
        representative.

        INPUT:

        - ``u``, ``v`` -- integers with `\gcd(u,v,N)=1`

        OUTPUT:

        - ``i`` -- the index of `u`, `v`, in the P1list

        - ``s`` -- normalizing scalar.
        """
        if self.__N == 1:
            # there is exactly 1 class [(0,0)].
            i[0] = 0
            s[0] = 1
            return

        cdef int uu, vv
        p1_normalize_xgcdtable(self.__N, u, v, 1, self.g, self.s, self.t, &uu, &vv, s)
        if uu == 1:
            i[0] = vv + 1
            return
        elif uu == 0:
            if vv == 0:
                i[0] = -1
                return
            i[0] = 0
            return
        try:
            i[0] = self.__end_hash[(uu,vv)] + self.__N + 1
            return
        except KeyError:
            i[0] = -1
            return

    def index_of_normalized_pair(self, int u, int v):
        r"""
        Return the index of the class of `(u,v)` in the fixed list
        of representatives of
        `\mathbb{P}^1(\ZZ/N\ZZ)`.

        INPUT:

        - ``u``, ``v`` -- integers with `\gcd(u,v,N)=1`, normalized so they lie
          in the list

        OUTPUT:

        - ``i`` -- the index of `(u:v)`, in the P1list

        EXAMPLES::

            sage: L = P1List(120)
            sage: L[100]
            (1, 99)
            sage: L.index_of_normalized_pair(1,99)
            100
            sage: all(L.index_of_normalized_pair(L[i][0],L[i][1])==i for i in range(len(L)))
            True
        """
        t, i = search(self.__list, (u,v))
        if t:
            return i
        return -1

    def list(self):
        r"""
        Return the underlying list of this :class:`P1List` object.

        EXAMPLES::

            sage: L = P1List(8)
            sage: type(L)
            <... 'sage.modular.modsym.p1list.P1List'>
            sage: type(L.list())
            <... 'list'>
        """
        return self.__list

    def normalize(self, int u, int v):
        r"""
        Return a normalised element of `\mathbb{P}^1(\ZZ/N\ZZ)`.

        INPUT:

        - ``u``, ``v`` -- integers with `\gcd(u,v,N)=1`

        OUTPUT: a 2-tuple ``(uu,vv)`` where `(uu:vv)` is a *normalized*
        representative of `(u:v)`

        NOTE: See also normalize_with_scalar() which also returns the
        normalizing scalar.

        EXAMPLES::

            sage: L = P1List(120)
            sage: (u,v) = (555555555,7777)
            sage: uu,vv = L.normalize(555555555,7777)
            sage: (uu,vv)
            (15, 13)
            sage: (uu*v-vv*u) % L.N() == 0
            True
        """
        cdef int uu, vv, ss
        self._normalize(self.__N, u, v, &uu, &vv, &ss, 0)
        return (uu,vv)

    def normalize_with_scalar(self, int u, int v):
        r"""
        Return a normalised element of `\mathbb{P}^1(\ZZ/N\ZZ)`, together with
        the normalizing scalar.

        INPUT:

        - ``u``, ``v`` -- integers with `\gcd(u,v,N)=1`

        OUTPUT:

        - a 3-tuple ``(uu,vv,ss)`` where `(uu:vv)` is a *normalized*
          representative of `(u:v)`, and `ss` is a scalar such that
          `(ss*uu, ss*vv) = (u,v)` (mod `N`).

        EXAMPLES::

            sage: L = P1List(120)
            sage: (u,v) = (555555555,7777)
            sage: uu,vv,ss = L.normalize_with_scalar(555555555,7777)
            sage: (uu,vv)
            (15, 13)
            sage: ((ss*uu-u)%L.N(), (ss*vv-v)%L.N())
            (0, 0)
            sage: (uu*v-vv*u) % L.N() == 0
            True
        """
        cdef int uu, vv, ss
        self._normalize(self.__N, u, v, &uu, &vv, &ss, 1)
        return (uu, vv, ss)

    def N(self):
        """
        Return the level or modulus of this P1List.

        EXAMPLES::

            sage: L = P1List(120)
            sage: L.N()
            120
        """
        return self.__N


cdef class export:
    cdef int c_p1_normalize_int(self, int N, int u, int v,
                                int* uu, int* vv, int* ss,
                                int compute_s) except -1:
        return c_p1_normalize_int(N, u, v, uu, vv, ss, compute_s)

    cdef int c_p1_normalize_llong(self, int N, int u, int v,
                                  int* uu, int* vv, int* ss,
                                  int compute_s) except -1:
        return c_p1_normalize_llong(N, u, v, uu, vv, ss, compute_s)


def lift_to_sl2z_int(int c, int d, int N):
    r"""
    Lift a pair `(c, d)` to an element of `SL(2, \ZZ)`.

    `(c,d)` is assumed to be an element of
    `\mathbb{P}^1(\ZZ/N\ZZ)`. This function computes
    and returns a list `[a, b, c', d']` that defines a 2x2
    matrix, with determinant 1 and integer entries, such that
    `c=c'` (mod `N`) and `d=d'` (mod `N`).

    INPUT:

    - ``c``, ``d``, ``N`` -- integers such that `\gcd(c,d,N)=1`

    EXAMPLES::

        sage: from sage.modular.modsym.p1list import lift_to_sl2z_int
        sage: lift_to_sl2z_int(2,6,11)
        [1, 8, 2, 17]
        sage: m = Matrix(Integers(),2,2,lift_to_sl2z_int(2,6,11)); m
        [ 1  8]
        [ 2 17]

    AUTHOR:

    - Justin Walker
    """
    cdef int z1, z2, g, m

    if N == 1:
        return [1, 0, 0, 1]

    if c == 0:
        if d == 1:
            return [1, 0, 0, 1]
        if d == N - 1:
            return [-1, 0, 0, -1]
        c = N

    g = arith_int.c_xgcd_int(c, d, &z1, &z2)

    # We're lucky: z1*c + z2*d = 1.
    if g==1:
        return [z2, -z1, c, d]

    # Have to try harder.
    m = c

    # compute prime-to-d part of m.
    while True:
        g = arith_int.c_gcd_int(m,d)
        if g == 1:
            break
        m = m / g

    # compute prime-to-N part of m.
    while True:
        g = arith_int.c_gcd_int(m,N)
        if g == 1:
            break
        m = m / g
    d += N * m
    g = arith_int.c_xgcd_int(c, d, &z1, &z2)

    if g != 1:
        raise ValueError("input must have gcd 1")

    return [z2, -z1, c, d]


def lift_to_sl2z_llong(llong c, llong d, int N):
    r"""
    Lift a pair `(c, d)` (modulo `N`) to an element of `SL(2, \ZZ)`.

    `(c,d)` is assumed to be an element of
    `\mathbb{P}^1(\ZZ/N\ZZ)`. This function computes and
    returns a list `[a, b, c', d']` that defines a 2x2 matrix,
    with determinant 1 and integer entries, such that `c=c'` (mod `N`)
    and `d=d'` (mod `N`).

    INPUT:

    - ``c``, ``d``, ``N`` -- integers such that `\gcd(c,d,N)=1`

    EXAMPLES::

        sage: from sage.modular.modsym.p1list import lift_to_sl2z_llong
        sage: lift_to_sl2z_llong(2,6,11)
        [1, 8, 2, 17]
        sage: m = Matrix(Integers(),2,2,lift_to_sl2z_llong(2,6,11)); m
        [ 1  8]
        [ 2 17]

    AUTHOR:

    - Justin Walker
    """
    cdef llong z1, z2, g, m

    if N == 1:
        return [1, 0, 0, 1]

    if c == 0:
        if d == 1:
            return [1, 0, 0, 1]
        if d == N - 1:
            return [-1, 0, 0, -1]
        c = N

    g = arith_llong.c_xgcd_longlong(c, d, &z1, &z2)

    # We're lucky: z1*c + z2*d = 1.
    if g==1:
        return [z2, -z1, c, d]

    # Have to try harder.
    m = c

    # compute prime-to-d part of m.
    while True:
        g = arith_llong.c_gcd_longlong(m,d)
        if g == 1:
            break
        m = m / g

    # compute prime-to-N part of m.
    while True:
        g = arith_llong.c_gcd_longlong(m,N)
        if g == 1:
            break
        m = m / g
    d += N * m
    g = arith_llong.c_xgcd_longlong(c, d, &z1, &z2)

    if g != 1:
        raise ValueError("input must have gcd 1")

    return [z2, -z1, c, d]


def lift_to_sl2z(c, d, N):
    r"""
    Return a list of Python ints `[a,b,c',d']` that are the entries of a
    2x2 matrix with determinant 1 and lower two entries congruent to
    `c,d` modulo `N`.

    INPUT:

    - ``c``, ``d``, ``N`` -- python ints or longs such that `\gcd(c,d,N)=1`

    EXAMPLES::

        sage: lift_to_sl2z(2,3,6)
        [1, 1, 2, 3]
        sage: lift_to_sl2z(2,3,6000000)
        [1, 1, 2, 3]

    You will get a :exc:`ValueError` exception if the input is invalid.
    Note that here gcd(15,6,24)=3::

        sage: lift_to_sl2z(15,6,24)
        Traceback (most recent call last):
        ...
        ValueError: input must have gcd 1

    This function is not implemented except for N at most 2**31::

        sage: lift_to_sl2z(1,1,2^32)
        Traceback (most recent call last):
        ...
        NotImplementedError: N too large

    TESTS::

        sage: lift_to_sl2z(0, 0, 1)
        [1, 0, 0, 1]
    """
    if N <= 46340:
        return lift_to_sl2z_int(c,d,N)
    elif N <= 2147483647:
        return lift_to_sl2z_llong(c,d,N)
    else:
        raise NotImplementedError("N too large")


def _make_p1list(n):
    """
    Helper function used in pickling.

    Not intended for end-users.

    EXAMPLES::

        sage: from sage.modular.modsym.p1list import _make_p1list
        sage: _make_p1list(3)
        doctest:...: DeprecationWarning: _make_p1list() is deprecated
        See https://github.com/sagemath/sage/issues/25848 for details.
        The projective line over the integers modulo 3
    """
    from sage.misc.superseded import deprecation_cython as deprecation
    deprecation(25848, '_make_p1list() is deprecated')
    return P1List(n)
