"""
Iterators over the partitions of an integer

AUTHOR:

- William Stein (2007-07-28): initial version
- Jonathan Bober (2007-07-28): wrote the program ``partitions_c.cc``
  that does all the actual heavy lifting.
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#       Copyright (C) 2007 Jonathan Bober <jwbober@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

def ZS1_iterator(int n):
    """
    A fast iterator for the partitions of ``n`` (in the decreasing
    lexicographic order) which returns lists and not objects of type
    :class:`~sage.combinat.partition.Partition`.

    This is an implementation of the ZS1 algorithm found in
    [ZS98]_.

    REFERENCES:

    .. [ZS98] Antoine Zoghbi, Ivan Stojmenovic,
       *Fast Algorithms for Generating Integer Partitions*,
       Intern. J. Computer Math., Vol. 70., pp. 319--332.
       http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.42.1287

    EXAMPLES::

        sage: from sage.combinat.partitions import ZS1_iterator
        sage: it = ZS1_iterator(4)
        sage: next(it)
        [4]
        sage: type(_)
        <class 'list'>
    """
    # Easy cases.
    if n < 0:
        return
    if n == 0:
        yield []
        return
    x = [1]*n
    x[0] = n

    cdef int m = 0
    cdef int h = 0
    cdef int r, t
    yield [n]
    while x[0] != 1:
        # Loop invariants at this point:
        # (A) x[:m+1] is a partition of n.
        # (B) x[h+1:] is an array of n-(h+1) ones.
        # (C) x[i] > 1 for each i <= h.
        # (D) 0 <= h <= m.
        if x[h] == 2:
            m += 1
            x[h] = 1
            h -= 1
        else:
            t = m - h + 1
            r = x[h] - 1
            x[h] = r
            while t >= r:
                h += 1
                x[h] = r
                t -= r
            if t == 0:
                m = h
            else:
                m = h + 1
                if t > 1:
                    h += 1
                    x[h] = t
        #yield [x[i] for i in range(m+1)]
        yield x[:m+1]
    #free(x)

def ZS1_iterator_nk(int n, int k):
    """
    An iterator for the partitions of ``n`` of length at most ``k`` (in the
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
    #free(x)
