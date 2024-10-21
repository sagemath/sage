# sage_setup: distribution = sagemath-categories
"""
Permutations (Cython file)

This is a nearly-straightforward implementation of what
Knuth calls "Algorithm P" in TAOCP 7.2.1.2.  The intent is
to be able to enumerate permutation by "plain changes", or
multiplication by adjacent transpositions, as a generator.
This is useful when a class of objects is inherently
enumerated by permutations, but it is faster to swap items
in a permutation than construct the next object directly
from the next permutation in a list. The backtracking
algorithm in sage/graphs/genus.pyx is an example of this.

The lowest level is implemented as a struct with auxiliary
methods.  This is because Cython does not allow pointers to
class instances, so a list of these objects is inherently
slower than a list of structs.  The author prefers ugly code
to slow code.

For those willing to sacrifice a (very small) amount of
speed, we provide a class that wraps our struct.
"""
# ****************************************************************************
#       Copyright (C) 2010 Tom Boothby <tomas.boothby@gmail.com>
#       Copyright (C) 2017 Travis Scrimshaw <tscrim@ucdavis.edu>
#       Copyright (C) 2017 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cimport cython

from cysignals.memory cimport check_allocarray, sig_free


#########################################################
#
# The next two functions, reset_swap and next_swap do the
# real work.  They've been implemented separately because
# the target application in sage/graphs/genus.pyx is
# cleaner and faster if it can manage its own memory.
# Both take three arguments:
#
#    int  n: number of elements to permute
#    int *c: an int array of length at least n
#    int *o: an int array of length at least o
#
# The user is advised to call reset_swap before next_swap.
#
##########################################################

cdef void reset_swap(int n, int *c, int *o) noexcept:
    """
    Reset the plain_swapper to the initial state.
    """

    for i in range(n):
        c[i] = -1
        o[i] = 1

cdef int next_swap(int n, int *c, int *o) noexcept:
    """
    Here's the translation of Algorithm P.  We've modified
    it to

       a) work on zero-indexed lists
       b) operate as a generator
       c) yield the swap index rather than the current
          permutation.

    Note, Knuth's descriptions of algorithms tend to encourage
    one to think of finite state machines.  For convenience,
    we have added comments to show what state the machine is
    in at any given point in the algorithm. `plain_swap_reset`
    sets the state to 1, and this function begins and ends in
    state 2.

    Returns the index i such that the next permutation can be
    obtained by swapping P[i] <-> P[i+1]
    """

    cdef int j,s,q,offset

    #state 3
    j = n-1
    if j <= 0:
        return -1
    s = -1

    while True:
        #state 4
        q = c[j] + o[j]
        if q == j:
            #state 6
            if j == 1:
                return -1
            s = s+1
        elif q >= -1:
            break

        #state 7
        o[j] = -o[j]
        j = j-1

    #state 5
    offset = c[j]
    if q > offset:
        offset = q
    c[j] = q
    return j - offset + s


def permutation_iterator_transposition_list(int n):
    """
    Return a list of transposition indices to enumerate the
    permutations on `n` letters by adjacent transpositions.
    Assumes zero-based lists.  We artificially limit the
    argument to `n < 12` to avoid overflowing 32-bit pointers.
    While the algorithm works for larger `n`, the user is
    encouraged to avoid filling anything more than 4GB of
    memory with the output of this function.

    EXAMPLES::

        sage: import sage.combinat.permutation_cython
        sage: from sage.combinat.permutation_cython import permutation_iterator_transposition_list
        sage: permutation_iterator_transposition_list(4)
        [2, 1, 0, 2, 0, 1, 2, 0, 2, 1, 0, 2, 0, 1, 2, 0, 2, 1, 0, 2, 0, 1, 2]
        sage: permutation_iterator_transposition_list(200)
        Traceback (most recent call last):
        ...
        ValueError: Cowardly refusing to enumerate the permutations on more than 12 letters.
        sage: permutation_iterator_transposition_list(1)
        []

        sage: # Generate the permutations of [1,2,3,4] fixing 4.
        sage: Q = [1,2,3,4]
        sage: L = [copy(Q)]
        sage: for t in permutation_iterator_transposition_list(3):
        ....:     Q[t], Q[t+1] = Q[t+1], Q[t]
        ....:     L.append(copy(Q))
        sage: print(L)
        [[1, 2, 3, 4], [1, 3, 2, 4], [3, 1, 2, 4], [3, 2, 1, 4], [2, 3, 1, 4], [2, 1, 3, 4]]
    """
    if n <= 1:
        return []
    if n > 12:
        raise ValueError("Cowardly refusing to enumerate the permutations "
                         "on more than 12 letters.")
    # Compute N = n! - 1
    cdef Py_ssize_t N = n
    cdef Py_ssize_t i
    for i in range(2, n):
        N *= i
    N -= 1

    cdef int* c = <int *>check_allocarray(n, 2 * sizeof(int))
    cdef int* o = c + n
    try:
        reset_swap(n, c, o)
        return [next_swap(n, c, o) for i in range(N)]
    finally:
        sig_free(c)


#####################################################################
## iterator-type method for getting the next permutation

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef bint next_perm(array l) noexcept:
    """
    Obtain the next permutation under lex order of ``l``
    by mutating ``l``.

    Algorithm based on:
    http://marknelson.us/2002/03/01/next-permutation/

    INPUT:

    - ``l`` -- array of unsigned int (i.e., type ``'I'``)

    .. WARNING::

        This method mutates the array ``l``.

    OUTPUT: boolean; whether another permutation was obtained

    EXAMPLES::

        sage: from sage.combinat.permutation_cython import next_perm
        sage: from array import array
        sage: L = array('I', [1, 1, 2, 3])
        sage: while next_perm(L):
        ....:     print(L)
        array('I', [1, 1, 3, 2])
        array('I', [1, 2, 1, 3])
        array('I', [1, 2, 3, 1])
        array('I', [1, 3, 1, 2])
        array('I', [1, 3, 2, 1])
        array('I', [2, 1, 1, 3])
        array('I', [2, 1, 3, 1])
        array('I', [2, 3, 1, 1])
        array('I', [3, 1, 1, 2])
        array('I', [3, 1, 2, 1])
        array('I', [3, 2, 1, 1])
    """
    cdef Py_ssize_t n = len(l)

    if n <= 1:
        return False

    cdef Py_ssize_t one = n - 2
    cdef Py_ssize_t two = n - 1
    cdef Py_ssize_t j   = n - 1
    cdef unsigned int t

    # Starting from the end, find the first o such that
    #   l[o] < l[o+1]
    while two > 0 and l.data.as_uints[one] >= l.data.as_uints[two]:
        one -= 1
        two -= 1

    if two == 0:
        return False

    #starting from the end, find the first j such that
    #l[j] > l[one]
    while l.data.as_uints[j] <= l.data.as_uints[one]:
        j -= 1

    #Swap positions one and j
    t = l.data.as_uints[one]
    l.data.as_uints[one] = l.data.as_uints[j]
    l.data.as_uints[j] = t

    #Reverse the list between two and last
    #mset_list = mset_list[:two] + [x for x in reversed(mset_list[two:])]
    n -= 1 # In the loop, we only need n-1, so just do it once here
    cdef Py_ssize_t i
    for i in range((n + 1 - two) // 2 - 1, -1, -1):
        t = l.data.as_uints[i + two]
        l.data.as_uints[i + two] = l.data.as_uints[n - i]
        l.data.as_uints[n - i] = t

    return True


@cython.boundscheck(False)
cpdef map_to_list(array l, tuple values, int n):
    """
    Build a list by mapping the array ``l`` using ``values``.

    .. WARNING::

        There is no check of the input data at any point. Using wrong
        types or values with wrong length is likely to result in a Sage
        crash.

    INPUT:

    - ``l`` -- array of unsigned int (i.e., type ``'I'``)
    - ``values`` -- tuple; the values of the permutation
    - ``n`` -- integer; the length of the array ``l``

    OUTPUT: list representing the permutation

    EXAMPLES::

        sage: from array import array
        sage: from sage.combinat.permutation_cython import map_to_list
        sage: l = array('I', [0, 1, 0, 3, 3, 0, 1])
        sage: map_to_list(l, ('a', 'b', 'c', 'd'), 7)
        ['a', 'b', 'a', 'd', 'd', 'a', 'b']
    """
    cdef int i
    cdef unsigned int* ind = l.data.as_uints
    return [values[ind[i]] for i in range(n)]


#####################################################################
## Multiplication functions for permutations

cpdef list left_action_same_n(list S, list lp):
    r"""
    Return the permutation obtained by composing a permutation
    ``S`` with a permutation ``lp`` in such an order that ``lp``
    is applied first and ``S`` is applied afterwards and ``S``
    and ``lp`` are of the same length.

    .. SEEALSO::

        :meth:`sage.combinat.permutation.Permutation.left_action_product`

    EXAMPLES::

        sage: p = [2,1,3]
        sage: q = [3,1,2]
        sage: from sage.combinat.permutation_cython import left_action_same_n
        sage: left_action_same_n(p, q)
        [3, 2, 1]
        sage: left_action_same_n(q, p)
        [1, 3, 2]
    """
    cdef int i
    cdef list ret = []
    for i in lp:
        ret.append(S[i-1])
    return ret

cpdef list right_action_same_n(list S, list rp):
    """
    Return the permutation obtained by composing a permutation
    ``S`` with a permutation ``rp`` in such an order that ``S`` is
    applied first and ``rp`` is applied afterwards and ``S`` and
    ``rp`` are of the same length.

    .. SEEALSO::

        :meth:`sage.combinat.permutation.Permutation.right_action_product`

    EXAMPLES::

        sage: p = [2,1,3]
        sage: q = [3,1,2]
        sage: from sage.combinat.permutation_cython import right_action_same_n
        sage: right_action_same_n(p, q)
        [1, 3, 2]
        sage: right_action_same_n(q, p)
        [3, 2, 1]
    """
    cdef int i
    cdef list ret = []
    for i in S:
        ret.append(rp[i-1])
    return ret

cpdef list left_action_product(list S, list lp):
    r"""
    Return the permutation obtained by composing a permutation
    ``S`` with a permutation ``lp`` in such an order that ``lp`` is
    applied first and ``S`` is applied afterwards.

    .. SEEALSO::

        :meth:`sage.combinat.permutation.Permutation.left_action_product`

    EXAMPLES::

        sage: p = [2,1,3,4]
        sage: q = [3,1,2]
        sage: from sage.combinat.permutation_cython import left_action_product
        sage: left_action_product(p, q)
        [3, 2, 1, 4]
        sage: left_action_product(q, p)
        [1, 3, 2, 4]
        sage: q
        [3, 1, 2]
    """
    cdef int i

    # Pad the permutations if they are of
    # different sizes
    S = S[:]
    lp = lp[:]
    for i in range(len(S)+1, len(lp)+1):
        S.append(i)
    for i in range(len(lp)+1, len(S)+1):
        lp.append(i)
    return left_action_same_n(S, lp)

cpdef list right_action_product(list S, list rp):
    """
    Return the permutation obtained by composing a permutation
    ``S`` with a permutation ``rp`` in such an order that ``S`` is
    applied first and ``rp`` is applied afterwards.

    .. SEEALSO::

        :meth:`sage.combinat.permutation.Permutation.right_action_product`

    EXAMPLES::

        sage: p = [2,1,3,4]
        sage: q = [3,1,2]
        sage: from sage.combinat.permutation_cython import right_action_product
        sage: right_action_product(p, q)
        [1, 3, 2, 4]
        sage: right_action_product(q, p)
        [3, 2, 1, 4]
        sage: q
        [3, 1, 2]
    """
    cdef int i

    # Pad the permutations if they are of
    # different sizes
    S = S[:]
    rp = rp[:]
    for i in range(len(S)+1, len(rp)+1):
        S.append(i)
    for i in range(len(rp)+1, len(S)+1):
        rp.append(i)
    return right_action_same_n(S, rp)
