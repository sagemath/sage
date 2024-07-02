r"""
Brent Yorgey's fast algorithm for integer vector (multiset) partitions.

ALGORITHM:

Brent Yorgey, Generating Multiset Partitions, The Monad Reader, Issue 8,
September 2007, p. 5.

https://wiki.haskell.org/The_Monad.Reader/Previous_issues

AUTHORS:

- D\. K\. Sunko (2020-02-19): initial version
- F\. Chapoton (2020-02-22): conversion to iterators and shorter doctests and
  doc tweaks
- T\. Scrimshaw (2020-03-06): Cython optimizations and doc tweaks
"""
################################################################################
#          Copyright (C) 2020 Denis Sunko <dks@phy.hr>                         #
#          Copyright (C) 2020 Frédéric Chapoton <chapoton@math.unistra.fr>     #
#          Copyright (C) 2020 Travis Scrimshaw <tscrims@gmail.com>             #
#                                                                              #
# This program is free software: you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation, either version 2 of the License, or            #
# (at your option) any later version. The text of the license is found at      #
#                  https://www.gnu.org/licenses/                               #
################################################################################
#
# To understand the code below, consult the ALGORITHM.

cdef list vector_halve(list v):
    r"""
    Return the vector halfway (lexicographically) between ``v`` and zero.

    Internal part of :func:`fast_vector_partitions`.

    INPUT:

    - ``v`` -- list of nonnegative integers, understood as a vector

    OUTPUT:

    A list, understood as the integer vector halfway down the list of
    lexicographically ordered vectors between ``v`` and zero.

    EXAMPLES::

        sage: from sage.combinat.fast_vector_partitions import vector_halve  # not tested
        sage: vector_halve([1, 2, 3, 4, 5, 6, 7, 8, 9])  # not tested
        [0, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: vector_halve([2, 4, 6, 8, 5, 6, 7, 8, 9])  # not tested
        [1, 2, 3, 4, 2, 6, 7, 8, 9]

    .. NOTE::

        For vectors, ``v = a + b`` implies ``v = b + a``, which means that a
        downward search for such splittings, starting with ``v = v + 0``, need
        only look as far as some "v / 2", given precise meaning here.

        A similar logic is to stop the search for divisors of ``N`` at
        ``sqrt(N)``, halving the exponents in the prime decomposition.
        However, here "v / 2" does not mean halving each coordinate.
    """
    cdef list result = list(v)  # make a copy
    cdef Py_ssize_t i, vv
    for i in range(len(v)):
        vv = <Py_ssize_t> v[i]
        result[i] = vv // 2
        if vv % 2:
            # the less significant part is just copied
            return result
    return result


def recursive_within_from_to(list m, list s, list e, bint useS, bint useE):
    r"""
    Iterate over a lexicographically ordered list of lists ``v`` satisfying
    ``e <= v <= s`` and ``v <|= m`` as vectors.

    Internal part of :func:`fast_vector_partitions`.

    INPUT:

    - ``m`` -- list of nonnegative integers, understood as a vector
    - ``s`` -- list of nonnegative integers, understood as a vector
    - ``e`` -- list of nonnegative integers, understood as a vector
    - ``useS`` -- boolean
    - ``useE`` -- boolean

    EXAMPLES::

        sage: from sage.combinat.fast_vector_partitions import recursive_within_from_to
        sage: list(recursive_within_from_to([1, 2, 3],[1, 2, 2],[1, 1, 1],True,True))
        [[1, 2, 2], [1, 2, 1], [1, 2, 0], [1, 1, 3], [1, 1, 2], [1, 1, 1]]

    .. NOTE::

        The flags ``useS`` and ``useE`` are used to implement the condition
        efficiently. Because testing it loops over the vector, re-testing
        at each step as the vector is parsed is inefficient: all but the last
        comparison have been done cumulatively already. This code tests
        only for the last one, using the flags to accumulate information
        from previous calls.

    .. WARNING::

        Expects to be called with ``s <|= m``.

        Expects to be called first with ``useS == useE == True``.
    """
    cdef Py_ssize_t start, end, x
    cdef bint useSS, useEE

    if useS:
        start = <Py_ssize_t> s[0]
    else:
        start = <Py_ssize_t> m[0]

    if useE:
        end = <Py_ssize_t> e[0]
    else:
        end = 0

    if len(m) == 1:
        # We use this style of Cython code for now in order to get this to convert
        #   to an optimized pure C for loop. See Cython github issue #532.
        # for x in range(start, end - 1, -1):
        for x from start >= x >= end by 1:
            yield [x]  # we know the answer for singletons
    else:
        # We use this style of Cython code for now in order to get this to convert
        #   to an optimized pure C for loop. See Cython github issue #532.
        # for x in range(start, end - 1, -1):
        for x from start >= x >= end by 1:
            useSS = useS and x == <Py_ssize_t> s[0]
            useEE = useE and x == <Py_ssize_t> e[0]
            for o in recursive_within_from_to(m[1:], s[1:], e[1:], useSS, useEE):
                yield [x] + o


def within_from_to(list m, list s, list e):
    r"""
    Iterate over a lexicographically ordered list of lists ``v`` satisfying
    ``e <= v <= s`` and ``v <|= m`` as vectors.

    Internal part of :func:`fast_vector_partitions`.

    INPUT:

    - ``m`` -- list of nonnegative integers, understood as a vector
    - ``s`` -- list of nonnegative integers, understood as a vector
    - ``e`` -- list of nonnegative integers, understood as a vector

    EXAMPLES::

        sage: from sage.combinat.fast_vector_partitions import within_from_to
        sage: list(within_from_to([1, 2, 3], [1, 2, 2], [1, 1, 1]))
        [[1, 2, 2], [1, 2, 1], [1, 2, 0], [1, 1, 3], [1, 1, 2], [1, 1, 1]]

    .. NOTE::

        The input ``s`` will be "clipped" internally if it does not satisfy
        the condition ``s <|= m``.

        To understand the input check, some line art is helpful. Assume
        that ``(a,b)`` are the two least significant coordinates of some
        vector. Say::

            e = (2,3), s = (7,6), m = (9,8).

        In the figure, these values are denoted by ``E``, ``S``, and ``M``,
        while the letter ``X`` stands for all other allowed values of
        ``v = (a,b)``::

            b ^
              |
            8 --------X---X---X---X---X-----------M
              |                                   |
            7 -       X   X   X   X   X           |
              |                                   |
            6 -       X   X   X   X   X   S       |
              |                                   |
            5 -       X   X   X   X   X   X       |
              |                                   |
            4 -       X   X   X   X   X   X       |
              |                                   |
            3 -       E   X   X   X   X   X       |
              |                                   |
            2 -           X   X   X   X   X       |
              |                                   |
            1 -           X   X   X   X   X       |
              |                                   |
            0 ----|---|---X---X---X---X---X---|---|--->
              0   1   2   3   4   5   6   7   8   9   a

        If ``S`` moves horizontally, the full-height columns fill the box in
        until ``S`` reaches ``M``, at which point it remains the limit in the
        b-direction as it moves out of the box, while M takes over as the
        limit in the a-direction, so the ``M``-column remains filled only up to
        ``S``, no matter how much ``S`` moves further to the right.

        If ``S`` moves vertically, its column will be filled to the top of the
        box, but it remains the relevant limit in the a-direction, while ``M``
        takes over in the b-direction as ``S`` goes out of the box upwards.

        Both behaviors are captured by using the smaller coordinate of ``S``
        and ``M``, whenever ``S`` is outside the box defined by M. The input
        will be "clipped" accordingly in that case.

    .. WARNING::

        The "clipping" behavior is transparent to the user, but may be
        puzzling when comparing outputs with the function
        :func:`recursive_within_from_to` which has no input protection.
    """
    cdef list ss = s
    # if s is not in the box defined by m, we must clip:
    cdef Py_ssize_t i, j
    for i in range(len(m)):  # should have the same length as s
        if s[i] > m[i]:
            ss = list(ss)  # make a copy
            ss[i] = m[i]
            for j in range(i+1, len(m)):
                if ss[j] > m[j]:
                    ss[j] = m[j]
            break
    if e > ss:
        return
    yield from recursive_within_from_to(m, ss, e, True, True)


cdef inline list vector_sub(list a, list b):
    """
    Return ``a - b`` considered as vectors.

    This assumes ``len(b) >= len(a)``.
    """
    cdef Py_ssize_t i
    cdef list ret = []
    for i in range(len(a)):
        ret.append((<Py_ssize_t> a[i]) - (<Py_ssize_t> b[i]))
    return ret


def recursive_vector_partitions(list v, list vL):
    r"""
    Iterate over a lexicographically ordered list of lists, each list
    representing a vector partition of ``v``, such that no part of any
    partition is lexicographically smaller than ``vL``.

    Internal part of :func:`fast_vector_partitions`.

    INPUT:

    - ``v`` -- list of nonnegative integers, understood as a vector
    - ``vL`` -- list of nonnegative integers, understood as a vector

    EXAMPLES::

        sage: from sage.combinat.fast_vector_partitions import recursive_vector_partitions
        sage: list(recursive_vector_partitions([2, 2, 2],[1, 1, 1]))
        [[[2, 2, 2]], [[1, 1, 1], [1, 1, 1]]]
        sage: list(recursive_vector_partitions([2, 2, 2],[1, 1, 0]))
        [[[2, 2, 2]], [[1, 1, 1], [1, 1, 1]], [[1, 1, 0], [1, 1, 2]]]
        sage: list(recursive_vector_partitions([2, 2, 2],[1, 0, 1]))
        [[[2, 2, 2]],
         [[1, 1, 1], [1, 1, 1]],
         [[1, 1, 0], [1, 1, 2]],
         [[1, 0, 2], [1, 2, 0]],
         [[1, 0, 1], [1, 2, 1]]]
    """
    cdef list v_minus_vv, pp, vv
    yield [v]
    for vv in within_from_to(v, vector_halve(v), vL):
        v_minus_vv = vector_sub(v, vv)
        for pp in recursive_vector_partitions(v_minus_vv, vv):
            yield [vv] + pp


def fast_vector_partitions(v, min_vals=None):
    r"""
    Brent Yorgey's fast algorithm for integer vector (multiset) partitions.

    INPUT:

    - ``v`` -- list of nonnegative integers, understood as the vector
      to be partitioned

    - ``min_vals`` -- (optional) list of nonnegative integers, of same
      length as ``v``

    OUTPUT: list of lists, each representing a vector partition of ``v``

    If ``min_vals`` is given, only partitions with parts ``p >= min_vals`` in
    the lexicographic ordering will appear.

    If ``min_vals`` is given and ``len(min_vals) != len(v)``, an error
    is raised.

    EXAMPLES:

    The older the computer, the more impressive the comparison::

        sage: from sage.combinat.fast_vector_partitions import fast_vector_partitions
        sage: fastvparts = list(fast_vector_partitions([3, 3, 3]))
        sage: vparts = list(VectorPartitions([3, 3, 3]))
        sage: vparts == fastvparts[::-1]
        True
        sage: len(fastvparts)
        686
        sage: list(fast_vector_partitions([1, 2, 3], min_vals=[0, 1, 1]))
        [[[1, 2, 3]],
         [[0, 2, 3], [1, 0, 0]],
         [[0, 2, 2], [1, 0, 1]],
         [[0, 2, 1], [1, 0, 2]],
         [[0, 2, 0], [1, 0, 3]],
         [[0, 1, 3], [1, 1, 0]],
         [[0, 1, 2], [1, 1, 1]],
         [[0, 1, 1], [1, 1, 2]],
         [[0, 1, 1], [0, 1, 2], [1, 0, 0]],
         [[0, 1, 1], [0, 1, 1], [1, 0, 1]]]
        sage: L1 = list(fast_vector_partitions([5, 7, 6], min_vals=[1, 3, 2]))
        sage: L1 == list(VectorPartitions([5, 7, 6], min=[1, 3, 2]))[::-1]
        True

    .. NOTE::

        The partitions are returned as an iterator.

        In this documentation, ``a <|= b`` means ``a[i] <= b[i]`` for all ``i``
        (notation following B. Yorgey's paper). It is the monomial partial
        ordering in Dickson's lemma: ``a <|= b`` iff ``x^a`` divides ``x^b`` as
        monomials.

    .. WARNING::

        The ordering of the partitions is reversed with respect to the output of
        Sage class :class:`~sage.combinat.vector_partition.VectorPartitions`.
    """
    if min_vals is None:
        min_vals = (len(v) - 1) * [0] + [1]  # lexicographically smallest vector > 0
    if len(v) != len(min_vals):
        raise ValueError("the length of v and min_vals must be equal")
    return recursive_vector_partitions(list(v), list(min_vals))
