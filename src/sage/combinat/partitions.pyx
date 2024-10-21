# sage_setup: distribution = sagemath-combinat
r"""
Iterators over the partitions of an integer

The iterators generate partitions in either increasing or decreasing
lexicographic orders and a partition `P` is represented either in ascending
(`P_i \leq P_{i+1}`) or descending (`P_i \geq P_{i+1}`) orders::

    sage: from sage.combinat.partitions import ZS1_iterator
    sage: for p in ZS1_iterator(4):
    ....:     print(p)
    [4]
    [3, 1]
    [2, 2]
    [2, 1, 1]
    [1, 1, 1, 1]
    sage: from sage.combinat.partitions import AccelDesc_iterator
    sage: for p in AccelDesc_iterator(4):
    ....:     print(p)
    [4]
    [3, 1]
    [2, 2]
    [2, 1, 1]
    [1, 1, 1, 1]
    sage: from sage.combinat.partitions import ZS2_iterator
    sage: for p in ZS2_iterator(4):
    ....:     print(p)
    [1, 1, 1, 1]
    [2, 1, 1]
    [2, 2]
    [3, 1]
    [4]
    sage: from sage.combinat.partitions import AccelAsc_iterator
    sage: for p in AccelAsc_iterator(4):
    ....:     print(p)
    [1, 1, 1, 1]
    [1, 1, 2]
    [1, 3]
    [2, 2]
    [4]

For each of these iterators, this module also provides a ``next`` method that
takes a partition as input and return the next partition in the corresponding
ordering::

    sage: from sage.combinat.partitions import ZS1_next
    sage: ZS1_next([2, 2])
    [2, 1, 1]
    sage: from sage.combinat.partitions import AccelDesc_next
    sage: AccelDesc_next([2, 2])
    [2, 1, 1]
    sage: from sage.combinat.partitions import ZS2_next
    sage: ZS2_next([2, 2])
    [3, 1]
    sage: from sage.combinat.partitions import AccelAsc_next
    sage: AccelAsc_next([2, 2])
    [4]

It is also possible to iterate over the partitions of bounded length::

    sage: from sage.combinat.partitions import ZS1_iterator_nk
    sage: for p in ZS1_iterator_nk(6, 3):
    ....:     print(p)
    [6]
    [5, 1]
    [4, 2]
    [4, 1, 1]
    [3, 3]
    [3, 2, 1]
    [2, 2, 2]

AUTHOR:

- William Stein (2007-07-28): initial version
- Jonathan Bober (2007-07-28): wrote the program ``partitions_c.cc``
  that does all the actual heavy lifting.
- David Coudert (2024-06-01): reshape method :meth:`ZS1_iterator` to ease the
  implementation of :meth:`ZS1_next` and add iterators and next methods based on
  ``ZS2``, ``AccelAsc`` and ``AccelDesc`` from [ZS1998]_ and [KS2012]_.
"""

# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#                     2007 Jonathan Bober <jwbober@gmail.com>
#                     2024 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


cdef inline ZS1_step(list P, int n, int *m, int *h):
    r"""
    Compute the partition following ``P`` in the ordering of the ZS1 algorithm.

    This is a helper method for methods :meth:`ZS1_iterator` and :meth:`ZS1_next`.
    Partition ``P`` is modified in place.

    INPUT:

    - ``P`` -- list of size `n` storing a partition of `n`

    - ``n`` -- integer; the sum of the elements of the partition stored in ``P``

    - ``m`` -- pointer to the index of the last element of the partition

    - ``h`` -- pointer to the last value larger than 1 in ``P``

    .. WARNING::

        The method assumes that input parameters are valid and consistent with
        the ZS1 algorithm as proposed in [ZS1998]_. It modifies the content of
        ``P`` and the values pointed by ``m`` and ``h``.

    EXAMPLES::

        sage: from sage.combinat.partitions import ZS1_iterator
        sage: it = ZS1_iterator(4)
        sage: next(it)
        [4]
        sage: type(_)
        <class 'list'>
    """
    # Invariants
    # (A) P[:m+1] is a partition of n.
    # (B) P[h+1:] is an array of n-(h+1) ones.
    # (C) P[i] > 1 for each i <= h.
    # (D) 0 <= h <= m.
    cdef int r, t
    cdef int mm = m[0]
    cdef int hh = h[0]
    if P[hh] == 2:
        mm += 1
        P[hh] = 1
        hh -= 1
    else:
        t = mm - hh + 1
        r = P[hh] - 1
        P[hh] = r
        while t >= r:
            hh += 1
            P[hh] = r
            t -= r
        if t == 0:
            mm = hh
        else:
            mm = hh + 1
            if t > 1:
                hh += 1
                P[hh] = t
    m[0] = mm
    h[0] = hh


def ZS1_iterator(int n):
    r"""
    Return an iterator over the partitions of ``n``.

    The partitions are generated in the decreasing lexicographic order and each
    partition is represented as a list ``P`` in descending order (i.e.,
    `P_i \geq P_{i+1}`). The method yields lists and not objects of type
    :class:`~sage.combinat.partition.Partition`.

    This is an implementation of the ZS1 algorithm found in [ZS1998]_.

    .. SEEALSO::

        - :meth:`sage.combinat.partitions.ZS2_iterator`
        - :meth:`sage.combinat.partitions.AccelDesc_iterator`

    EXAMPLES::

        sage: from sage.combinat.partitions import ZS1_iterator
        sage: for p in ZS1_iterator(4):
        ....:     print(p)
        [4]
        [3, 1]
        [2, 2]
        [2, 1, 1]
        [1, 1, 1, 1]
        sage: next(ZS1_iterator(4))
        [4]
        sage: type(_)
        <class 'list'>
    """
    # Easy cases.
    if n < 0:
        return
    if not n:
        yield []
        return
    cdef list P = [1] * n
    P[0] = n
    cdef int m = 0
    cdef int h = 0
    yield [n]
    while P[0] != 1:
        ZS1_step(P, n, &m, &h)
        yield P[:m+1]


def ZS1_next(list P):
    r"""
    Return the partition after ``P`` in the ordering of the ZS1 algorithm.

    INPUT:

    - ``P`` -- list encoding a partition of an integer `n` in descending order
      (i.e., `P_i \geq P_{i+1}`)

    EXAMPLES::

        sage: from sage.combinat.partitions import ZS1_iterator, ZS1_next
        sage: P = [4]
        sage: while P:
        ....:     print(P)
        ....:     P = ZS1_next(P)
        [4]
        [3, 1]
        [2, 2]
        [2, 1, 1]
        [1, 1, 1, 1]
        sage: A = [list(p) for p in Partitions(7)]
        sage: all(ZS1_next(p) == q for p, q in zip(A, A[1:]))
        True
    """
    if P[0] == 1:
        return
    if P[-1] != 1:
        return P[:-1] + [P[-1] - 1, 1]

    cdef int n = sum(P)
    cdef int m = len(P) - 1  # index of the last element in P
    # Search for the index h of the last value > 1 in P
    cdef int h = m
    while P[h] == 1:
        h -= 1
    # Let Q be such that Q[:m+1] == P and Q[h+1:] is an array of n-(h+1) ones
    cdef list Q = P + [1] * (n - m - 1)
    ZS1_step(Q, n, &m, &h)
    return Q[:m+1]


def ZS1_iterator_nk(int n, int k):
    r"""
    An iterator for the partitions of `n` of length at most `k` (in the
    decreasing lexicographic order) which returns lists and not objects of type
    :class:`~sage.combinat.partition.Partition`.

    The algorithm is a mild variation on :func:`ZS1_iterator`;
    I would not vow for its speed.

    EXAMPLES::

        sage: from sage.combinat.partitions import ZS1_iterator_nk
        sage: it = ZS1_iterator_nk(4, 3)
        sage: next(it)
        [4]
        sage: type(_)
        <class 'list'>
    """
    # Easy cases.
    if n <= 0:
        if n == 0 and k >= 0:
            yield []
        return
    if k <= 1:
        if k == 1:
            yield [n]
        return
    x = [1]*k
    x[0] = n

    cdef int m = 0
    cdef int h = 0
    cdef int r, t
    yield [n]
    while x[0] != 1:
        # Loop invariants at this point:
        # (A) x[:m+1] is a partition of n.
        # (B) x[h+1:m+1] is an array of m-h ones.
        # (C) x[i] > 1 for each i <= h.
        # (D) 0 <= h <= m < k.
        # Note that x[m+1:] might contain leftover from
        # previous steps; we don't clean up after ourselves.
        if x[h] == 2 and m + 1 < k:
            # We have a 2 in the partition, and the space to
            # spread it into two 1s.
            m += 1
            x[h] = 1
            x[m] = 1
            h -= 1
            yield x[:m+1]
        else:
            t = m - h + 1  # 1 + "the number of 1s to the right of x[h] that belong to the partition"
            r = x[h] - 1

            # This loop finds the largest h such that x[:h] can be completed
            # with integers smaller-or-equal to r=x[h]-1 into a partition of n.
            #
            # We decrement h until it becomes possible.
            while t > (k-h-1) * r:
                # Loop invariants:
                # t = n - sum(x[:h+1]) + 1;
                # r = x[h] - 1; x[h] > 1.
                if h == 0:
                    # No way to make the current partition
                    # lexicographically smaller.
                    return
                h -= 1
                t += r + 1
                r = x[h] - 1
            # Decrement x[h] from r + 1 to r, and replace
            # x[h+1:] by the lexicographically highest array
            # it could possibly be. This means replacing
            # x[h+1:] by the array [r, r, r, ..., r, s],
            # where s is the residue of t modulo r (or
            # nothing if that residue is 0).
            x[h] = r
            while t >= r:
                # Loop invariants: t = n - sum(x[:h+1]) + 1;
                # r = x[h] > 1.
                h += 1
                x[h] = r
                t -= r
            if t == 0:
                m = h
            else:
                m = h + 1
                if t > 1:
                    h += 1
                x[m] = t
            yield x[:m+1]
    # free(x)


cdef inline ZS2_step(list P, int n, int *m, int *h):
    r"""
    Compute the partition following ``P`` in the ordering of the ZS2 algorithm.

    This is a helper method for methods :meth:`ZS2_iterator` and :meth:`ZS2_next`.
    Partition ``P`` is modified in place.

    INPUT:

    - ``P`` -- list of size `n` storing a partition of `n`

    - ``n`` -- integer; the sum of the elements of the partition stored in ``P``

    - ``m`` -- pointer to the index of the last element of the partition

    - ``h`` -- pointer to the last value larger than 1 in ``P``

    .. WARNING::

        The method assumes that input parameters are valid and consistent with
        the ZS2 algorithm as proposed in [ZS1998]_. It modifies the content of
        ``P`` and the values pointed by ``m`` and ``h``.

    EXAMPLES::

        sage: from sage.combinat.partitions import ZS2_iterator
        sage: it = ZS2_iterator(4)
        sage: next(it)
        [1, 1, 1, 1]
        sage: type(_)
        <class 'list'>
    """
    # Invariants
    # (A) P[:m+1] is a partition of n.
    # (B) P[h+1:] is an array of n-(h+1) ones.
    # (C) P[i] > 1 for each i <= h.
    # (D) 0 <= h <= m.
    cdef int mm = m[0]
    cdef int hh = h[0]
    cdef int t, r
    if mm - hh > 1:
        hh += 1
        P[hh] = 2
        mm -= 1
    else:
        t = mm - 2
        while P[t] == P[mm - 1]:
            P[t] = 1
            t -= 1
        hh = t + 1
        P[hh] = P[mm - 1] + 1
        r = P[mm] + P[mm - 1] * (mm - hh - 1)
        P[mm] = 1
        if mm - hh > 1:
            P[mm - 1] = 1
        mm = hh + r - 1
    m[0] = mm
    h[0] = hh


def ZS2_iterator(int n):
    r"""
    Return an iterator over the partitions of ``n``.

    The partitions are generated in the increasing lexicographic order and each
    partition is represented as a list in descending order (i.e., `p_i \geq
    p_{i+1}`).

    This is an implementation of the ZS2 algorithm found in [ZS1998]_.

    .. SEEALSO::

        :meth:`sage.combinat.partitions.ZS1_iterator`

    EXAMPLES::

        sage: from sage.combinat.partitions import ZS2_iterator
        sage: for p in ZS2_iterator(4):
        ....:     print(p)
        [1, 1, 1, 1]
        [2, 1, 1]
        [2, 2]
        [3, 1]
        [4]
        sage: next(ZS2_iterator(4))
        [1, 1, 1, 1]
        sage: type(_)
        <class 'list'>
    """
    # Easy cases.
    if n < 0:
        return
    yield [1] * n
    if n <= 1:
        return
    cdef list P = [1]*n
    P[0] = 2
    cdef int h = 0
    cdef int m = n - 2  # index of the last element of the partition
    yield P[:m+1]
    while P[0] != n:
        ZS2_step(P, n, &m, &h)
        yield P[:m+1]


def ZS2_next(list P):
    r"""
    Return the partition after ``P`` in the ordering of the ZS2 algorithm.

    INPUT:

    - ``P`` -- list encoding a partition of an integer `n` in descending order
      (i.e., `P_i \geq P_{i+1}`)

    EXAMPLES::

        sage: from sage.combinat.partitions import ZS2_iterator, ZS2_next
        sage: P = [1, 1, 1, 1]
        sage: while P:
        ....:     print(P)
        ....:     P = ZS2_next(P)
        [1, 1, 1, 1]
        [2, 1, 1]
        [2, 2]
        [3, 1]
        [4]
    """
    if len(P) <= 1:
        return
    if P[0] == 1:
        return [2] + P[:-2]
    cdef int n = sum(P)
    cdef int m = len(P) - 1  # index of the last element in P
    # Search for the index h of the last value > 1 in P
    cdef int h = m
    while P[h] == 1:
        h -= 1
    # Let Q be such that Q[:m+1] == P and Q[h+1:] is an array of n-(h+1) ones
    cdef list Q = P + [1] * (n - m - 1)
    ZS2_step(Q, n, &m, &h)
    return Q[:m+1]


cdef inline AccelDesc_step(list P, int n, int* m, int* h):
    r"""
    Compute the partition following ``P`` in the ordering of the AccelDesc algorithm.

    This is a helper method for methods :meth:`AccelDesc_iterator` and
    :meth:`AccelDesc_next`. Partition ``P`` is modified in place.

    INPUT:

    - ``P`` -- list of size `n` storing a partition of `n`

    - ``n`` -- integer; the sum of the elements of the partition stored in ``P``

    - ``m`` -- pointer to the index of the last element of the partition

    - ``h`` -- pointer to the last value larger than 1 in ``P``

    .. WARNING::

        The method assumes that input parameters are valid and consistent with
        the AccelDesc algorithm as proposed in [ZS2012]_. It modifies the
        content of ``P`` and the values pointed by ``m`` and ``h``.

    EXAMPLES::

        sage: from sage.combinat.partitions import AccelDesc_iterator
        sage: it = AccelDesc_iterator(4)
        sage: next(it)
        [4]
        sage: type(_)
        <class 'list'>
    """
    # Invariants
    # (A) P[:m+1] is a partition of n.
    # (B) P[h+1:] is an array of n-(h+1) ones.
    # (C) P[i] > 1 for each i <= h.
    # (D) 0 <= h <= m.
    cdef int mm = m[0]
    cdef int hh = h[0]
    cdef int r, t
    if P[hh] == 2:
        mm += 1
        P[hh] = 1
        hh -= 1
    else:
        r = P[hh] - 1
        t = mm - hh + 1
        P[hh] = r
        while t >= r:
            hh += 1
            P[hh] = r
            t -= r
        if not t:
            mm = hh
        else:
            mm = hh + 1
            if t > 1:
                hh += 1
                P[hh] = t
    m[0] = mm
    h[0] = hh


def AccelDesc_iterator(int n):
    r"""
    Return an iterator over the partitions of ``n``.

    The partitions are generated in the increasing lexicographic order and each
    partition is represented as a list in descending order (i.e., `p_i \geq
    p_{i+1}`).

    This is an implementation of the AccelDesc algorithm found in [KS2012]_.

    .. SEEALSO::

        :meth:`sage.combinat.partitions.ZS1_iterator`

    EXAMPLES::

        sage: from sage.combinat.partitions import AccelDesc_iterator
        sage: for p in AccelDesc_iterator(4):
        ....:     print(p)
        [4]
        [3, 1]
        [2, 2]
        [2, 1, 1]
        [1, 1, 1, 1]
        sage: next(AccelDesc_iterator(4))
        [4]
        sage: type(_)
        <class 'list'>

    Check that :meth:`ZS1_iterator` and :meth:`AccelDesc_iterator` generate
    partitions in the same order::

        sage: from sage.combinat.partitions import ZS1_iterator
        sage: from sage.misc.prandom import randint
        sage: n = randint(1, 50)
        sage: all(p == q for p, q in zip(ZS1_iterator(n), AccelDesc_iterator(n)))  # long time
        True
    """
    # Easy cases.
    if n < 0:
        return
    if n <= 1:
        yield [1] * n
        return
    cdef list P = [1] * n
    P[0] = n
    yield [n]
    cdef int m = 0
    cdef int h = 0
    while h >= 0:
        AccelDesc_step(P, n, &m, &h)
        yield P[:m+1]


def AccelDesc_next(list P):
    r"""
    Return the partition after ``P`` in the ordering of the AccelDesc algorithm.

    INPUT:

    - ``P`` -- list encoding a partition of an integer `n` in descending order
      (i.e., `P_i \geq P_{i+1}`)

    EXAMPLES::

        sage: from sage.combinat.partitions import AccelDesc_iterator, AccelDesc_next
        sage: P = [4]
        sage: while P:
        ....:     print(P)
        ....:     P = AccelDesc_next(P)
        [4]
        [3, 1]
        [2, 2]
        [2, 1, 1]
        [1, 1, 1, 1]
        sage: A = [list(p) for p in Partitions(7)]
        sage: all(AccelDesc_next(p) == q for p, q in zip(A, A[1:]))
        True
    """
    if P[0] == 1:
        return
    if P[-1] != 1:
        return P[:-1] + [P[-1] - 1, 1]

    cdef int n = sum(P)
    cdef int m = len(P) - 1  # index of the last element in P
    # Search for the index h of the last value > 1 in P
    cdef int h = m
    while P[h] == 1:
        h -= 1
    # Let Q be such that Q[:m+1] == P and Q[h+1:] is an array of n-(h+1) ones
    cdef list Q = P + [1]*(n - m - 1)
    AccelDesc_step(Q, n, &m, &h)
    return Q[:m+1]


def AccelAsc_iterator(int n):
    r"""
    Return an iterator over the partitions of ``n``.

    The partitions are generated in the increasing lexicographic order and each
    partition is represented as a list in ascending order (i.e., `p_i \leq
    p_{i+1}`).

    This is an implementation of the ``AccelAsc`` algorithm found in [KS2012]_.

    .. SEEALSO::

        :meth:`sage.combinat.partitions.ZS1_iterator`

    EXAMPLES::

        sage: from sage.combinat.partitions import AccelAsc_iterator
        sage: for p in AccelAsc_iterator(4):
        ....:     print(p)
        [1, 1, 1, 1]
        [1, 1, 2]
        [1, 3]
        [2, 2]
        [4]
        sage: next(AccelAsc_iterator(4))
        [1, 1, 1, 1]
        sage: type(_)
        <class 'list'>
    """
    # Easy cases.
    if n < 0:
        return
    if n <= 1:
        yield [1] * n
        return
    cdef list P = [0] * n
    cdef int k = 1
    cdef int y = n - 1
    cdef int x, ell
    while k:
        k -= 1
        x = P[k] + 1
        while 2 * x <= y:
            P[k] = x
            y -= x
            k += 1
        ell = k + 1
        while x <= y:
            P[k] = x
            P[ell] = y
            yield P[:ell + 1]
            x += 1
            y -= 1
        y += x - 1
        P[k] = y + 1
        yield P[:k + 1]


def AccelAsc_next(list P):
    r"""
    Return the partition after ``P`` in the ordering of the ``AccelAsc`` algorithm.

    INPUT:

    - ``P`` -- list encoding a partition of an integer `n` in ascending order
      (i.e., `P_i \leq P_{i+1}`)

    EXAMPLES::

        sage: from sage.combinat.partitions import AccelAsc_next
        sage: P = [1, 1, 1, 1]
        sage: while P:
        ....:     print(P)
        ....:     P = AccelAsc_next(P)
        [1, 1, 1, 1]
        [1, 1, 2]
        [1, 3]
        [2, 2]
        [4]
    """
    if len(P) <= 1:
        return

    cdef int n = sum(P)
    cdef int m = len(P) - 1  # index of the last element in P
    # Let Q be such that Q[:m+1] == P and Q[m+1:] is an array of n-(m+1) ones
    cdef list Q = P + [0]*(n - m - 1)
    cdef int x, y
    cdef int k = m - 1
    x = Q[k] + 1
    y = P[m] - 1
    while 2 * x <= y:
        Q[k] = x
        y -= x
        k += 1
    if x <= y:
        Q[k] = x
        Q[k + 1] = y
        return Q[:k + 2]
    y += x - 1
    Q[k] = y + 1
    return Q[:k + 1]
