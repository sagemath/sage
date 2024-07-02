# sage.doctest: needs sage.combinat
r"""
Sage code for computing k-distant crossing numbers.

This code accompanies the article :arxiv:`0812.2725`. It is being
submitted because of a suggestion from
http://groups.google.com/group/sage-support/msg/3ea7ed2eeab0824a.

Right now, this code only computes k-dcrossings. If you are only
interested in the distribution, this is good enough because the extended
Kasraoui-Zeng involution tells us the distribution of k-dcrossings and
k-dnestings is symmetric. It would be nice, though, to have a function
which actually performed that involution.

AUTHORS:

- Dan Drake (2008-12-15): initial version.

EXAMPLES:

The example given in the paper. Note that in this format, we omit fixed
points since they cannot create any sort of crossing. ::

    sage: from sage.tests.arxiv_0812_2725 import *
    sage: dcrossing([(1,5), (2,4), (4,9), (6,12), (7,10), (10,11)])
    3
"""
#*****************************************************************************
# Copyright (C) 2008 Dan Drake <ddrake@member.ams.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or (at
# your option) any later version.
#
# See https://www.gnu.org/licenses/.
#*****************************************************************************
from sage.combinat.set_partition import SetPartitions as SetPartitions


def CompleteMatchings(n):
    """
    Return a generator for the complete matchings of the set [1..n].

    INPUT:

    - ``n`` -- nonnegative integer

    OUTPUT:

    A generator for the complete matchings of the set [1..n], or,
    what is basically the same thing, complete matchings of the
    graph K_n. Each complete matching is represented by a list of
    2-element tuples.

    EXAMPLES:

    There are 3 complete matchings on 4 vertices::

        sage: from sage.tests.arxiv_0812_2725 import *
        sage: [m for m in CompleteMatchings(4)]
        [[(3, 4), (1, 2)], [(2, 4), (1, 3)], [(2, 3), (1, 4)]]

    There are no complete matchings on an odd number of vertices; the
    number of complete matchings on an even number of vertices is a
    double factorial::

        sage: from sage.tests.arxiv_0812_2725 import *
        sage: [len([m for m in CompleteMatchings(n)]) for n in [0..8]]
        [1, 0, 1, 0, 3, 0, 15, 0, 105]

    The exact behavior of CompleteMatchings(n) if n is not a nonnegative
    integer depends on what [1..n] returns, and also on what range(1,
    len([1..n])) is.
    """
    yield from matchingsset(list(range(1, n + 1)))


def matchingsset(L):
    """
    Return a generator for complete matchings of the sequence L.

    This is not really meant to be called directly, but rather by
    CompleteMatchings().

    INPUT:

    - ``L`` -- a sequence. Lists, tuples, et cetera; anything that
      supports ``len()`` and slicing should work.

    OUTPUT:

    A generator for complete matchings on K_n, where n is the length
    of L and vertices are labeled by elements of L. Each matching is
    represented by a list of 2-element tuples.

    EXAMPLES::

        sage: from sage.tests.arxiv_0812_2725 import *
        sage: [m for m in matchingsset(('a', 'b', 'c', 'd'))]
        [[('c', 'd'), ('a', 'b')], [('b', 'd'), ('a', 'c')], [('b', 'c'), ('a', 'd')]]

    There's only one matching of the empty set/list/tuple: the empty
    matching. ::

        sage: [m for m in matchingsset(())]
        [[]]
    """
    if not L:
        yield []
    else:
        for k in range(1, len(L)):
            for m in matchingsset(L[1: k] + L[k + 1:]):
                yield m + [(L[0], L[k])]


def dcrossing(m_):
    """
    Return the largest k for which the given matching or set
    partition has a k-distant crossing.

    INPUT:

    - ``m`` -- a matching or set partition, as a list of 2-element tuples
      representing the edges. You'll need to call ``setp_to_edges()`` on
      the objects returned by ``SetPartitions()`` to put them into the
      proper format.

    OUTPUT:

    The largest k for which the object has a k-distant crossing.
    Matchings and set partitions with no crossings at all yield -1.

    EXAMPLES:

    The main example from the paper::

        sage: from sage.tests.arxiv_0812_2725 import *
        sage: dcrossing(setp_to_edges(Set(map(Set, [[1,5],[2,4,9],[3],[6,12],[7,10,11],[8]]))))
        3

    A matching example::

        sage: from sage.tests.arxiv_0812_2725 import *
        sage: dcrossing([(4, 7), (3, 6), (2, 5), (1, 8)])
        2

    TESTS:

    The empty matching and set partition are noncrossing::

        sage: dcrossing([])
        -1
        sage: dcrossing(Set([]))
        -1

    One edge::

        sage: dcrossing([Set((1,2))])
        -1
        sage: dcrossing(Set([Set((1,2))]))
        -1

    Set partition with block of size >= 3 is always at least
    0-dcrossing::

        sage: dcrossing(setp_to_edges(Set([Set((1,2,3))])))
        0
    """
    d = -1
    m = list(m_)
    while m:
        e1_ = m.pop()
        for e2_ in m:
            e1, e2 = sorted(e1_), sorted(e2_)
            if (e1[0] < e2[0] and e2[0] <= e1[1] and e1[1] < e2[1] and
                    e1[1] - e2[0] > d):
                d = e1[1] - e2[0]
            if (e2[0] < e1[0] and e1[0] <= e2[1] and e2[1] < e1[1] and
                    e2[1] - e1[0] > d):
                d = e2[1] - e1[0]
    return d


def setp_to_edges(p):
    """
    Transform a set partition into a list of edges.

    INPUT:

    - ``p`` -- a Sage set partition

    OUTPUT:

    A list of non-loop edges of the set partition. As this code just
    works with crossings, we can ignore the loops.

    EXAMPLES:

    The main example from the paper::

        sage: from sage.tests.arxiv_0812_2725 import *
        sage: sorted(setp_to_edges(Set(map(Set, [[1,5],[2,4,9],[3],[6,12],[7,10,11],[8]]))))
        [[1, 5], [2, 4], [4, 9], [6, 12], [7, 10], [10, 11]]
    """
    q = (sorted(b) for b in p)
    return [b[n: n + 2] for b in q for n in range(len(b) - 1)]


def dcrossvec_setp(n):
    """
    Return a list with the distribution of k-dcrossings on set partitions of [1..n].

    INPUT:

    - ``n`` -- nonnegative integer

    OUTPUT:

    A list whose k'th entry is the number of set partitions p for
    which dcrossing(p) = k. For example, let L = dcrossvec_setp(3).
    We have L = [1, 0, 4]. L[0] is 1 because there's 1 partition of
    [1..3] that has 0-dcrossing: [(1, 2, 3)].

    One tricky bit is that noncrossing matchings get put at the end,
    because L[-1] is the last element of the list. Above, we have
    L[-1] = 4 because the other four set partitions are all
    d-noncrossing. Because of this, you should not think of the last
    element of the list as having index n-1, but rather -1.

    EXAMPLES::

        sage: from sage.tests.arxiv_0812_2725 import *
        sage: dcrossvec_setp(3)
        [1, 0, 4]

        sage: dcrossvec_setp(4)
        [5, 1, 0, 9]

    The one set partition of 1 element is noncrossing, so the last
    element of the list is 1::

        sage: dcrossvec_setp(1)
        [1]
    """
    vec = [0] * n
    for p in SetPartitions(n):
        vec[dcrossing(setp_to_edges(p))] += 1
    return vec


def dcrossvec_cm(n):
    """
    Return a list with the distribution of k-dcrossings on complete matchings on n vertices.

    INPUT:

    - ``n`` -- nonnegative integer

    OUTPUT:

    A list whose k'th entry is the number of complete matchings m
    for which dcrossing(m) = k. For example, let L =
    dcrossvec_cm(4). We have L = [0, 1, 0, 2]. L[1] is 1 because
    there's one matching on 4 vertices that is 1-dcrossing: [(2, 4),
    (1, 3)]. L[0] is zero because dcrossing() returns the *largest*
    k for which the matching has a dcrossing, and 0-dcrossing is
    equivalent to 1-dcrossing for complete matchings.

    One tricky bit is that noncrossing matchings get put at the end,
    because L[-1] is the last element of the list. Because of this, you
    should not think of the last element of the list as having index
    n-1, but rather -1.

    If n is negative, you get silly results. Don't use them in your
    next paper. :)

    EXAMPLES:

    The single complete matching on 2 vertices has no crossings, so the
    only nonzero entry of the list (the last entry) is 1::

        sage: from sage.tests.arxiv_0812_2725 import *
        sage: dcrossvec_cm(2)
        [0, 1]

    Similarly, the empty matching has no crossings::

        sage: dcrossvec_cm(0)
        [1]

    For odd n, there are no complete matchings, so the list has all
    zeros::

        sage: dcrossvec_cm(5)
        [0, 0, 0, 0, 0]

        sage: dcrossvec_cm(4)
        [0, 1, 0, 2]
    """
    vec = [0] * max(n, 1)
    for m in CompleteMatchings(n):
        vec[dcrossing(m)] += 1
    return vec


def tablecolumn(n, k):
    """
    Return column n of Table 1 or 2 from the paper :arxiv:`0812.2725`.

    INPUT:

    - ``n`` -- positive integer

    - ``k`` -- integer for which table you want: Table 1 is complete
      matchings, Table 2 is set partitions

    OUTPUT:

    The n'th column of the table as a list. This is basically just the
    partial sums of dcrossvec_{cm,setp}(n).

    table2column(1, 2) incorrectly returns [], instead of [1], but you
    probably don't need this function to work through n = 1.

    EXAMPLES:

    Complete matchings::

        sage: from sage.tests.arxiv_0812_2725 import *
        sage: tablecolumn(2, 1)
        [1]

        sage: tablecolumn(6, 1)
        [5, 5, 11, 14, 15]

    Set partitions::

        sage: tablecolumn(5, 2)
        [21, 42, 51, 52]

        sage: tablecolumn(2, 2)
        [2]
    """
    if k == 1:
        v = dcrossvec_cm(n)
    else:
        v = dcrossvec_setp(n)
    i = v[-1]
    return [i + sum(v[:k]) for k in range(len(v) - 1)]
