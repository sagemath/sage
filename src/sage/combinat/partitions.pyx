"""
Number of partitions of integer

AUTHOR:

- William Stein (2007-07-28): initial version
- Jonathan Bober (2007-07-28): wrote the program ``partitions_c.cc``
  that does all the actual heavy lifting.
"""

import sys

from sage.libs.gmp.types cimport mpz_t

cdef extern from "partitions_c.h":
    int part(mpz_t answer, unsigned int n)
    int test(bint longtest, bint forever)

include "sage/ext/interrupt.pxi"

from sage.rings.integer cimport Integer

def number_of_partitions(n):
    """
    Returns the number of partitions of the integer `n`.

    EXAMPLES::

        sage: from sage.combinat.partitions import number_of_partitions
        sage: number_of_partitions(3)
        3
        sage: number_of_partitions(10)
        42
        sage: number_of_partitions(40)
        37338
        sage: number_of_partitions(100)
        190569292
        sage: number_of_partitions(100000)
        27493510569775696512677516320986352688173429315980054758203125984302147328114964173055050741660736621590157844774296248940493063070200461792764493033510116079342457190155718943509725312466108452006369558934464248716828789832182345009262853831404597021307130674510624419227311238999702284408609370935531629697851569569892196108480158600569421098519

    TESTS::

        sage: n = 500 + randint(0,500)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1500 + randint(0,1500)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 100000000 + randint(0,100000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0  # long time (4s on sage.math, 2011)
        True

    Another consistency test for `n` up to 500::

        sage: len([n for n in [1..500] if number_of_partitions(n) != Partitions(n).cardinality(algorithm='pari')])
        0

    """
    n = Integer(n)
    if n < 0:
        raise ValueError("n (=%s) must be a nonnegative integer"%n)
    elif n <= 1:
        return Integer(1)  # part hangs on n=1 as input.
    if n >= Integer('4294967296'):
        raise ValueError("input must be a nonnegative integer less than 4294967296.")
    cdef unsigned int nn = n

    cdef Integer ans = Integer(0)

    sig_on()
    part(ans.value, nn)
    sig_off()

    return ans

def run_tests(bint longtest=False, bint forever=False):
    """
    Test Bober's algorithm.

    EXAMPLES::

        sage: from sage.combinat.partitions import run_tests
        sage: run_tests(False, False)
        Computing p(1)... OK.
        ...
        Done.
    """
    sig_on()
    error = test(longtest, forever)
    sig_off()
    print "Done."
    if error:
        return error

def ZS1_iterator(int n):
    """
    A fast iterator for the partitions of ``n`` (in the decreasing
    lexicographic order) which returns lists and not objects of type
    :class:`~sage.combinat.partition.Partition`.

    This is an implementation of the ZS1 algorithm found in
    [ZS98]_.

    REFERENCES:

    .. [ZS98] Antoine Zoghbi, Ivan Stojmenovic,
       *Fast Algorithms for Generating Integer Partitons*,
       Intern. J. Computer Math., Vol. 70., pp. 319--332.
       http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.42.1287

    EXAMPLES::

        sage: from sage.combinat.partitions import ZS1_iterator
        sage: it = ZS1_iterator(4)
        sage: it.next()
        [4]
        sage: type(_)
        <type 'list'>
    """
    # Easy cases.
    if n < 0:
        return
    if n == 0:
        yield []
        return
    #cdef int *x = <int*>malloc(sizeof(int) *n)
    #x[0] = n
    #cdef int i
    #for i in range(1, n):
    #    x[i] = 1
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
        #yield [x[i] for i in xrange(m+1)]
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
        sage: it.next()
        [4]
        sage: type(_)
        <type 'list'>
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
    #cdef int *x = <int*>malloc(sizeof(int) *n)
    #x[0] = n
    #cdef int i
    #for i in range(1, n):
    #    x[i] = 1
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
            t = m - h + 1 # 1 + "the number of 1s to the right of x[h] that belong to the partition"
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

