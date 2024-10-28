# cython: binding=True
r"""
Fast set partition iterators
"""
cimport cython


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

        sage: from sage.combinat.set_partition_iterator import set_partition_iterator
        sage: list(set_partition_iterator([1,-1,x]))                                    # needs sage.symbolic
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

        sage: from sage.combinat.set_partition_iterator import _set_partition_block_gen
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

        sage: from sage.combinat.set_partition_iterator import set_partition_iterator_blocks
        sage: list(set_partition_iterator_blocks([1,-1,x], 2))                          # needs sage.symbolic
        [[[1, x], [-1]], [[1], [-1, x]], [[1, -1], [x]]]
    """
    cdef list base = list(base_set)
    cdef Py_ssize_t n = len(base)
    cdef list a = list(range(n))
    # TODO: implement _set_partition_block_gen as an iterative algorithm
    if n < k:
        return
    for P in _set_partition_block_gen(n, k, a):
        yield from_word(<list> P, base)
