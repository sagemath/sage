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

cimport cython

from cysignals.memory cimport check_allocarray, sig_free

from sage.libs.gmp.all cimport *
from sage.rings.integer cimport Integer

cdef void mpz_addmul_alt(mpz_t s, mpz_t t, mpz_t u, unsigned long parity):
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
## Lyndon word iterator

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

#####################################################################
## Set partition iterators

@cython.wraparound(False)
@cython.boundscheck(False)
cdef list from_word(list w, list base_set):
    cdef list sp = []
    cdef Py_ssize_t i
    cdef Py_ssize_t b
    for i in range(len(w)):
        b = <Py_ssize_t> (w[i])
        x = base_set[i]
        if len(sp) <= b:
            sp.append([x])
        else:
            sp[b].append(x)
    return sp

@cython.wraparound(False)
@cython.boundscheck(False)
def set_partition_iterator(base_set):
    """
    A fast iterator for the set partitions of the base set, which
    returns lists of lists instead of set partitions types.

    EXAMPLES::

        sage: from sage.combinat.combinat_cython import set_partition_iterator
        sage: list(set_partition_iterator([1,-1,x]))
        [[[1, -1, x]],
         [[1, -1], [x]],
         [[1, x], [-1]],
         [[1], [-1, x]],
         [[1], [-1], [x]]]
    """
    cdef list base = list(base_set)

    # Knuth, TAOCP 4A 7.2.1.5, Algorithm H
    cdef Py_ssize_t N = len(base)
    # H1: initialize
    cdef list a = [0] * N
    if N <= 1:
        yield from_word(a, base)
        return

    cdef list b = [1] * N
    cdef Py_ssize_t j
    cdef Py_ssize_t last = N - 1
    while True:
        # H2: visit
        yield from_word(a, base)
        if a[last] == b[last]:
            # H4: find j
            j = N - 2
            while a[j] == b[j]:
                j -= 1
            # H5: increase a_j
            if j == 0:
                break
            a[j] += 1
            # H6: zero out a_{j+1},...,a_{n-1}
            b[last] = b[j] + int(a[j] == b[j])
            j += 1
            while j < N - 1:
                a[j] = 0
                b[j] = b[last]
                j += 1
            a[last] = 0
        else:
            # H3: increase a_{n-1}
            a[last] += 1

@cython.wraparound(False)
@cython.boundscheck(False)
def _set_partition_block_gen(Py_ssize_t n, Py_ssize_t k, list a):
    r"""
    Recursively generate set partitions of ``n`` with fixed block
    size ``k`` using Algorithm 4.23 from [Rus2003]_.
    ``a`` is a list of size ``n``.

    EXAMPLES::

        sage: from sage.combinat.combinat_cython import _set_partition_block_gen
        sage: a = list(range(3))
        sage: for p in _set_partition_block_gen(3, 2, a):
        ....:     print(p)
        [0, 1, 0]
        [0, 1, 1]
        [0, 0, 1]
    """
    cdef Py_ssize_t i
    if n == k:
        yield a
        return

    for i in range(k):
        a[n-1] = i
        for P in _set_partition_block_gen(n-1, k, a):
            yield P
        a[n-1] = n-1
    if k > 1:
        a[n-1] = k-1
        for P in _set_partition_block_gen(n-1, k-1, a):
            yield P
        a[n-1] = n-1

@cython.wraparound(False)
@cython.boundscheck(False)
def set_partition_iterator_blocks(base_set, Py_ssize_t k):
    """
    A fast iterator for the set partitions of the base set into the
    specified number of blocks, which returns lists of lists
    instead of set partitions types.

    EXAMPLES::

        sage: from sage.combinat.combinat_cython import set_partition_iterator_blocks
        sage: list(set_partition_iterator_blocks([1,-1,x], 2))
        [[[1, x], [-1]], [[1], [-1, x]], [[1, -1], [x]]]
    """
    cdef list base = list(base_set)
    cdef Py_ssize_t n = len(base)
    cdef list a = list(range(n))
    # TODO: implement _set_partition_block_gen as an iterative algorithm
    for P in _set_partition_block_gen(n, k, a):
        yield from_word(<list> P, base)

## Perfect matchings iterator

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
## Linear extension iterator

from copy import copy
def _linear_extension_prepare(D):
    r"""
    The preprocessing routine in Figure 7 of "Generating Linear
    Extensions Fast" by Preusse and Ruskey.

    INPUT:

    - ``D``, the Hasse diagram of a poset

    OUTPUT:

    - a triple ``(le, a, b)``, where ``le`` is the first linear
      extension, and ``a`` and ``b`` are lists such that ``a[i]`` and
      ``b[i]`` are minimal elements of ``D`` after removing ``a[:i]``
      and ``b[:i]``.

    TESTS::

        sage: from sage.combinat.combinat_cython import _linear_extension_prepare
        sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })._hasse_diagram
        sage: _linear_extension_prepare(D)
        ([0, 1, 2, 3, 4], [1, 3], [2, 4])

    """
    dag_copy = copy(D) # this copy is destroyed during preparation
    le = []
    a = []
    b = []

    # the preprocessing routine found in Figure 7 of
    # "Generating Linear Extensions Fast" by
    # Pruesse and Ruskey
    while dag_copy.num_verts() != 0:
        # find all the minimal elements of dag_copy
        minimal_elements = dag_copy.sources()
        if not minimal_elements:
            raise ValueError("the digraph must be acyclic to have linear extensions")
        elif len(minimal_elements) == 1:
            le.append(minimal_elements[0])
            dag_copy.delete_vertex(minimal_elements[0])
        else:
            ap = minimal_elements[0]
            bp = minimal_elements[1]
            a.append(ap)
            b.append(bp)
            le.append(ap)
            le.append(bp)
            dag_copy.delete_vertex(ap)
            dag_copy.delete_vertex(bp)

    return (le, a, b)

@cython.wraparound(False)
@cython.boundscheck(False)
cdef void _linear_extension_switch(list _le, list _a, list _b, list _is_plus, Py_ssize_t i):
    """
    This implements the ``Switch`` procedure described on page 7
    of "Generating Linear Extensions Fast" by Pruesse and Ruskey.

    If ``i == -1``, then the sign is changed.  Otherwise, then
    ``_a[i]`` and ``_b[i]`` are transposed.

    """
    cdef Py_ssize_t a_index, b_index
    if i == -1:
        _is_plus[0] = not _is_plus[0]
    else:
        a = _a[i]; b = _b[i]
        a_index = _le.index(a)
        b_index = _le.index(b)
        _le[a_index] = b
        _le[b_index] = a
        _b[i] = a; _a[i] = b

@cython.wraparound(False)
@cython.boundscheck(False)
cdef bint _linear_extension_right_a(_D, list _le, list _a, list _b, Py_ssize_t i):
    """
    Return ``True`` if and only if ``_a[i]`` is incomparable with the
    element to its right in ``_le`` and the element to the right is
    not ``_b[i]``.

    This is the ``Right`` function described on page 8 of
    "Generating Linear Extensions Fast" by Pruesse and Ruskey.

    ::

        sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })._hasse_diagram             # not tested
        sage: _linear_extension_right_a(D, [0, 1, 2, 4, 3], [1, 4], [2, 3], 0)  # not tested
        False
        sage: _linear_extension_right_a(D, [0, 1, 2, 4, 3], [1, 4], [2, 3], 1)  # not tested
        False

    """
    cdef Py_ssize_t yindex
    x = _a[i]
    yindex = _le.index(x) + 1
    if yindex >= len(_le):
        return False
    y = _le[yindex]
    return y != _b[i] and _D.are_incomparable(x, y)

@cython.wraparound(False)
@cython.boundscheck(False)
cdef bint _linear_extension_right_b(_D, list _le, list _a, list _b, Py_ssize_t i):
    """
    Return True if and only if ``_b[i]`` is incomparable with the
    elements to its right in ``_le``.

    This is the ``Right`` function described on page 8 of
    "Generating Linear Extensions Fast" by Pruesse and Ruskey.

    ::

        sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })._hasse_diagram             # not tested
        sage: _linear_extension_right_b(D, [0, 1, 2, 4, 3], [1, 4], [2, 3], 0)  # not tested
        False
        sage: _linear_extension_right_b(D, [0, 1, 2, 4, 3], [1, 4], [2, 3], 1)  # not tested
        False

    """
    cdef Py_ssize_t yindex
    x = _b[i]
    yindex = _le.index(x) + 1
    if yindex >= len(_le):
        return False
    y = _le[yindex]
    return _D.are_incomparable(x, y)

@cython.wraparound(False)
@cython.boundscheck(False)
def _linear_extension_gen(_D, list _le, list _a, list _b, list _is_plus, Py_ssize_t i):
    """
    This a Python version of the GenLE routine found in Figure 8
    of "Generating Linear Extensions Fast" by Pruesse and Ruskey.

    TESTS::

        sage: from sage.combinat.combinat_cython import _linear_extension_prepare, _linear_extension_gen
        sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })._hasse_diagram
        sage: le, a, b = _linear_extension_prepare(D)
        sage: [e for e in _linear_extension_gen(D, le, a, b, [True], len(a)-1)]
        [[0, 2, 1, 3, 4]]

    """
    cdef int mra, mrb, mla
    cdef Py_ssize_t index, index1
    cdef bint typical
    if i == -1:
        return

    for e in _linear_extension_gen(_D, _le, _a, _b, _is_plus, i-1):
        yield e
    mrb = 0
    typical = False
    while _linear_extension_right_b(_D, _le, _a, _b, i):
        mrb += 1
        # move_right
        index = _le.index(_b[i]); index1 = index + 1
        _le[index] = _le[index1]
        _le[index1] = _b[i]
        if _is_plus[0]:
            yield _le[:]

        for e in _linear_extension_gen(_D, _le, _a, _b, _is_plus, i-1):
            yield e
        mra = 0
        while _linear_extension_right_a(_D, _le, _a, _b, i):
            typical = True
            mra += 1
            # move_right
            index = _le.index(_a[i]); index1 = index+1
            _le[index] = _le[index1]
            _le[index1] = _a[i]
            if _is_plus[0]:
                yield _le[:]

            for e in _linear_extension_gen(_D, _le, _a, _b, _is_plus, i-1):
                yield e

        if typical:
            _linear_extension_switch(_le, _a, _b, _is_plus, i-1)
            if _is_plus[0]:
                yield _le[:]

            for e in _linear_extension_gen(_D, _le, _a, _b, _is_plus, i-1):
                yield e
            if mrb % 2 == 1:
                mla = mra - 1
            else:
                mla = mra + 1
            for _ in range(mla):
                # move_left
                index = _le.index(_a[i]); index1 = index-1
                _le[index] = _le[index1]
                _le[index1] = _a[i]
                if _is_plus[0]:
                    yield _le[:]

                for e in _linear_extension_gen(_D, _le, _a, _b, _is_plus, i-1):
                    yield e

    if typical and (mrb % 2 == 1):
        # move_left
        index = _le.index(_a[i]); index1 = index-1
        _le[index] = _le[index1]
        _le[index1] = _a[i]
        if _is_plus[0]:
            yield _le[:]
    else:
        _linear_extension_switch(_le, _a, _b, _is_plus, i-1)
        if _is_plus[0]:
            yield _le[:]
    for e in _linear_extension_gen(_D, _le, _a, _b, _is_plus, i-1):
        yield e
    for _ in range(mrb):
        # move_left
        index = _le.index(_b[i]); index1 = index-1
        _le[index] = _le[index1]
        _le[index1] = _b[i]
        if _is_plus[0]:
            yield _le[:]

        for e in _linear_extension_gen(_D, _le, _a, _b, _is_plus, i-1):
            yield e


def linear_extension_iterator(D):
    """
    Iterate over the linear extensions of the poset.

    The list ``_le`` keeps track of the current linear extensions.  The
    boolean variable ``is_plus`` keeps track of the "sign".

    INPUT:

    - ``D``, the Hasse diagram of a poset.

    .. WARNING::

        It is assumed that ``D`` is not modified while the linear
        extensions are generated.

    EXAMPLES::

        sage: from sage.combinat.combinat_cython import linear_extension_iterator
        sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })._hasse_diagram
        sage: list(linear_extension_iterator(D))
        [[0, 1, 2, 3, 4],
         [0, 2, 1, 3, 4],
         [0, 2, 1, 4, 3],
         [0, 2, 4, 1, 3],
         [0, 1, 2, 4, 3]]

        sage: D = posets.BooleanLattice(3)._hasse_diagram
        sage: len(list(linear_extension_iterator(D)))
        48

        sage: D = posets.AntichainPoset(9)._hasse_diagram
        sage: len(list(linear_extension_iterator(D))) == factorial(9)           # long time
        True
    """
    _le, _a, _b = _linear_extension_prepare(D)
    _max_pair = len(_a) - 1
    _is_plus = [True] # this is modified by _linear_extension_switch

    yield _le[:]
    for e in _linear_extension_gen(D, _le, _a, _b, _is_plus, _max_pair):
        yield e
    _linear_extension_switch(_le, _a, _b, _is_plus, _max_pair)
    if _is_plus[0]:
        yield _le[:]
    for e in _linear_extension_gen(D, _le, _a, _b, _is_plus, _max_pair):
        yield e

#####################################################################
## Set partition composition

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
        sage: (SetPartition(p), c) == (SetPartition([[1,-1],[2,-2]]), 0)
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
    assert all(all(val > 0 for val in top) for top in remaining_top)
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
