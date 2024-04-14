# sage_setup: distribution = sagemath-categories
# cython: binding=True
"""
Fast computation of combinatorial functions (Cython + mpz)

Currently implemented:

- Stirling numbers of the second kind
- iterators for set partitions
- iterator for Lyndon words
- iterator for perfect matchings
- conjugate of partitions

AUTHORS:

- Fredrik Johansson (2010-10): Stirling numbers of second kind
- Martin Rubey and Travis Scrimshaw (2018): iterators for set partitions,
  Lyndon words, and perfect matchings
"""

from cysignals.memory cimport check_allocarray, sig_free

from sage.libs.gmp.all cimport *
from sage.rings.integer cimport Integer
from sage.misc.lazy_import import LazyImport

set_partition_iterator = LazyImport('sage.combinat.set_partition_iterator', 'set_partition_iterator', deprecation=35741)
set_partition_iterator_blocks = LazyImport('sage.combinat.set_partition_iterator', 'set_partition_iterator_blocks', deprecation=35741)
linear_extension_iterator = LazyImport('sage.combinat.posets.linear_extension_iterator', 'linear_extension_iterator', deprecation=35741)


cdef void mpz_addmul_alt(mpz_t s, mpz_t t, mpz_t u, unsigned long parity) noexcept:
    """
    Set s = s + t*u * (-1)^parity
    """
    if parity & 1:
        mpz_submul(s, t, u)
    else:
        mpz_addmul(s, t, u)


cdef mpz_stirling_s2(mpz_t s, unsigned long n, unsigned long k):
    """
    Set s = S(n,k) where S(n,k) denotes a Stirling number of the
    second kind.

    Algorithm: S(n,k) = (sum_{j=0}^k (-1)^(k-j) C(k,j) j^n) / k!

    TODO: compute S(n,k) efficiently for large n when n-k is small
    (e.g. when k > 20 and n-k < 20)
    """
    cdef mpz_t t, u
    cdef mpz_t *bc
    cdef unsigned long j, max_bc
    # Some important special cases
    if k+1 >= n:
        # Upper triangle of n\k table
        if k > n:
            mpz_set_ui(s, 0)
        elif n == k:
            mpz_set_ui(s, 1)
        elif k+1 == n:
            # S(n,n-1) = C(n,2)
            mpz_set_ui(s, n)
            mpz_mul_ui(s, s, n-1)
            mpz_tdiv_q_2exp(s, s, 1)
    elif k <= 2:
        # Leftmost three columns of n\k table
        if k == 0:
            mpz_set_ui(s, 0)
        elif k == 1:
            mpz_set_ui(s, 1)
        elif k == 2:
            # 2^(n-1)-1
            mpz_set_ui(s, 1)
            mpz_mul_2exp(s, s, n-1)
            mpz_sub_ui(s, s, 1)
    # Direct sequential evaluation of the sum
    elif n < 200:
        mpz_init(t)
        mpz_init(u)
        mpz_set_ui(t, 1)
        mpz_set_ui(s, 0)
        for j in range(1, k//2+1):
            mpz_mul_ui(t, t, k+1-j)
            mpz_tdiv_q_ui(t, t, j)
            mpz_set_ui(u, j)
            mpz_pow_ui(u, u, n)
            mpz_addmul_alt(s, t, u, k+j)
            if 2*j != k:
                # Use the fact that C(k,j) = C(k,k-j)
                mpz_set_ui(u, k-j)
                mpz_pow_ui(u, u, n)
                mpz_addmul_alt(s, t, u, j)
        # Last term not included because loop starts from 1
        mpz_set_ui(u, k)
        mpz_pow_ui(u, u, n)
        mpz_add(s, s, u)
        mpz_fac_ui(t, k)
        mpz_tdiv_q(s, s, t)
        mpz_clear(t)
        mpz_clear(u)
    # Only compute odd powers, saving about half of the time for large n.
    # We need to precompute binomial coefficients since they will be accessed
    # out of order, adding overhead that makes this slower for small n.
    else:
        mpz_init(t)
        mpz_init(u)
        max_bc = (k+1)//2
        bc = <mpz_t*> check_allocarray(max_bc+1, sizeof(mpz_t))
        mpz_init_set_ui(bc[0], 1)
        for j in range(1, max_bc+1):
            mpz_init_set(bc[j], bc[j-1])
            mpz_mul_ui(bc[j], bc[j], k+1-j)
            mpz_tdiv_q_ui(bc[j], bc[j], j)
        mpz_set_ui(s, 0)
        for j in range(1, k+1, 2):
            mpz_set_ui(u, j)
            mpz_pow_ui(u, u, n)
            # Process each 2^p * j, where j is odd
            while True:
                if j > max_bc:
                    mpz_addmul_alt(s, bc[k-j], u, k+j)
                else:
                    mpz_addmul_alt(s, bc[j], u, k+j)
                j *= 2
                if j > k:
                    break
                mpz_mul_2exp(u, u, n)
        for j in range(max_bc+1):   # careful: 0 ... max_bc
            mpz_clear(bc[j])
        sig_free(bc)
        mpz_fac_ui(t, k)
        mpz_tdiv_q(s, s, t)
        mpz_clear(t)
        mpz_clear(u)


def _stirling_number2(n, k):
    """
    Python wrapper of mpz_stirling_s2.

        sage: from sage.combinat.combinat_cython import _stirling_number2
        sage: _stirling_number2(3, 2)
        3

    This is wrapped again by stirling_number2 in combinat.py.
    """
    cdef Integer s = Integer.__new__(Integer)
    mpz_stirling_s2(s.value, n, k)
    return s


#####################################################################
#  Lyndon word iterator

def lyndon_word_iterator(Py_ssize_t n, Py_ssize_t k):
    r"""
    Generate the Lyndon words of fixed length ``k`` with ``n`` letters.

    The resulting Lyndon words will be words represented as lists
    whose alphabet is ``range(n)`` (`= \{0, 1, \ldots, n-1\}`).

    ALGORITHM:

    The iterative FKM Algorithm 7.2 from [Rus2003]_.

    EXAMPLES::

        sage: from sage.combinat.combinat_cython import lyndon_word_iterator
        sage: list(lyndon_word_iterator(4, 2))
        [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
        sage: list(lyndon_word_iterator(2, 4))
        [[0, 0, 0, 1], [0, 0, 1, 1], [0, 1, 1, 1]]

    TESTS::

        sage: from sage.combinat.combinat_cython import lyndon_word_iterator
        sage: list(lyndon_word_iterator(6, 1))
        [[0], [1], [2], [3], [4], [5]]
        sage: list(lyndon_word_iterator(5, 0))
        []
        sage: list(lyndon_word_iterator(1, 1000))
        []
        sage: list(lyndon_word_iterator(1, 1))
        [[0]]
    """
    cdef Py_ssize_t i, j
    if k == 0:
        return
    if k == 1:
        for i in range(n):
            yield [i]
        return
    if n == 1:
        return

    cdef list a = [0] * (k+1)
    i = k
    while i != 0:
        a[i] += 1
        for j in range(1, k-i+1):
            a[j + i] = a[j]
        if k == i:
            yield a[1:]
        i = k
        while a[i] == n - 1:
            i -= 1


#  Perfect matchings iterator

def perfect_matchings_iterator(Py_ssize_t n):
    r"""
    Iterate over all perfect matchings with ``n`` parts.

    This iterates over all perfect matchings of `\{0, 1, \ldots, 2n-1\}`
    using a Gray code for fixed-point-free involutions due to Walsh [Wal2001]_.

    EXAMPLES::

        sage: from sage.combinat.combinat_cython import perfect_matchings_iterator
        sage: list(perfect_matchings_iterator(1))
        [[(0, 1)]]
        sage: list(perfect_matchings_iterator(2))
        [[(0, 1), (2, 3)], [(0, 2), (1, 3)], [(0, 3), (1, 2)]]

        sage: list(perfect_matchings_iterator(0))
        [[]]

    REFERENCES:

    - [Wal2001]_
    """
    if n == 0:
        yield []
        return

    cdef Py_ssize_t i, x, y, g, j, J
    cdef Py_ssize_t* e = <Py_ssize_t*> check_allocarray(2*n, sizeof(Py_ssize_t))
    for i in range(2*n):
        e[i] = i
    cdef Py_ssize_t* f = <Py_ssize_t*> check_allocarray(2*n, sizeof(Py_ssize_t))
    for i in range(2*n):
        if i % 2 == 0:
            f[i] = i + 1
        else:
            f[i] = i - 1
    cdef bint odd = False

    yield convert(f, n)
    while e[0] != n - 1:
        i = e[0]
        if odd:
            x = 2 * i
        else:
            x = i

        y = f[x]
        g = y - x - 1
        if g % 2 == odd:
            g += 1
            j = y + 1
        else:
            g -= 1
            j = y-1
        J = f[j]
        f[y] = J
        f[J] = y
        f[x] = j
        f[j] = x
        odd = not odd
        e[0] = 0
        if g == 0 or g == 2 * (n-i-1):
            e[i] = e[i+1]
            e[i+1] = i + 1

        yield convert(f, n)

    sig_free(e)
    sig_free(f)


cdef list convert(Py_ssize_t* f, Py_ssize_t n):
    """
    Convert a list ``f`` representing a fixed-point free involution
    to a set partition.
    """
    cdef list ret = []
    cdef Py_ssize_t i
    for i in range(2*n):
        if i < f[i]:
            ret.append((i, f[i]))
    return ret


#####################################################################
#  Set partition composition

def set_partition_composition(tuple sp1, tuple sp2):
    r"""
    Return a tuple consisting of the composition of the set partitions
    ``sp1`` and ``sp2`` and the number of components removed from the middle
    rows of the graph.

    EXAMPLES::

        sage: from sage.combinat.combinat_cython import set_partition_composition
        sage: sp1 = ((1,-2),(2,-1))
        sage: sp2 = ((1,-2),(2,-1))
        sage: p, c = set_partition_composition(sp1, sp2)
        sage: (SetPartition(p), c) == (SetPartition([[1,-1],[2,-2]]), 0)                # needs sage.combinat
        True
    """
    cdef int num_loops = 0  # The number of loops removed
    cdef list diagram = []  # The resulting composite diagram
    # sp1 is on top of sp2
    # positive values on top and negative on bottom
    cdef set remaining_top = set(sp1)
    cdef list remaining_bot = list(sp2)
    cdef list cur = []
    cdef tuple temp, top, to_remove
    cdef list block = []
    cdef Py_ssize_t i
    while remaining_bot or cur:
        if not cur:
            cur = list(remaining_bot.pop())
            block = []
        while cur:
            val = cur.pop()
            if val > 0:
                # Find what it is connected to in sp1
                to_remove = ()
                for top in remaining_top:
                    if -val in top:
                        to_remove = top
                        for entry in top:
                            if entry < 0:
                                # Check to see if that makes a new connection with
                                #   something in sp2.
                                # We go through this in reverse order so that when we
                                #   pop an element off, we do not need to update i.
                                for i in reversed(range(len(remaining_bot))):
                                    temp = <tuple> remaining_bot[i]
                                    if -entry in temp:
                                        remaining_bot.pop(i)
                                        cur.extend(temp)
                                        continue
                            else:
                                block.append(entry)
                        break
                if to_remove:
                    remaining_top.remove(to_remove)
            else:
                block.append(val)
        if not cur:
            if not block:
                num_loops += 1
            else:
                diagram.append(tuple(block))

    # Everything else should be completely contained in the top block
    assert all(val > 0 for top in remaining_top for val in top)
    diagram.extend(remaining_top)

    return (tuple(diagram), num_loops)


def conjugate(p):
    """
    Return the conjugate partition associated to the partition ``p``
    as a list.

    EXAMPLES::

        sage: from sage.combinat.combinat_cython import conjugate
        sage: conjugate([2,2])
        [2, 2]
        sage: conjugate([6,3,1])
        [3, 2, 2, 1, 1, 1]
    """
    cdef Py_ssize_t j, l
    cdef list conj
    l = len(p)
    if l == 0:
        return []
    conj = [l] * p[-1]
    for j in range(l - 1, 0, -1):
        conj.extend([j] * (p[j - 1] - p[j]))
    return conj
