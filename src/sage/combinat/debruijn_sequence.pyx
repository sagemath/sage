r"""
De Bruijn sequences

A De Bruijn sequence is defined as the shortest cyclic sequence that
incorporates all substrings of a certain length of an alphabet.

For instance, the `2^3=8` binary strings of length 3 are all included in the
following string::

    sage: DeBruijnSequences(2,3).an_element()
    [0, 0, 0, 1, 0, 1, 1, 1]

They can be obtained as a subsequence of the *cyclic* De Bruijn sequence of
parameters `k=2` and `n=3`::

    sage: seq = DeBruijnSequences(2,3).an_element()
    sage: print(Word(seq).string_rep())
    00010111
    sage: shift = lambda i: [(i+j)%2**3 for j in range(3)]
    sage: for i in range(2**3):
    ....:    w = Word([b if j in shift(i) else '*' for j, b in enumerate(seq)])
    ....:    print(w.string_rep())
    000*****
    *001****
    **010***
    ***101**
    ****011*
    *****111
    0*****11
    00*****1

This sequence is of length `k^n`, which is best possible as it is the number of
`k`-ary strings of length `n`. One can equivalently define a De Bruijn sequence
of parameters `k` and `n` as a cyclic sequence of length `k^n` in which all
substring of length `n` are different.

See also :wikipedia:`De_Bruijn_sequence`.

TESTS:

Checking the sequences generated are indeed valid::

    sage: for n in range(1, 7):
    ....:    for k in range(1, 7):
    ....:       D = DeBruijnSequences(k, n)
    ....:       if not D.an_element() in D:
    ....:           print("Something's dead wrong (n=%s, k=%s)!" %(n,k))
    ....:           break

AUTHOR:

- Eviatar Bach (2011): initial version

- Nathann Cohen (2011): Some work on the documentation and defined the
  ``__contain__`` method
"""

# ******************************************************************************
#         Copyright (C) 2011 Eviatar Bach <eviatarbach@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         https://www.gnu.org/licenses/
# ******************************************************************************

from sage.data_structures.bitset_base cimport *


def debruijn_sequence(int k, int n):
    """
    The generating function for De Bruijn sequences. This avoids the object
    creation, so is significantly faster than accessing from DeBruijnSequence.
    For more information, see the documentation there. The algorithm used is
    from Frank Ruskey's "Combinatorial Generation".

    INPUT:

    - ``k`` -- arity; must be an integer

    - ``n`` -- substring length; must be an integer

    EXAMPLES::

        sage: from sage.combinat.debruijn_sequence import debruijn_sequence
        sage: debruijn_sequence(3, 1)
        [0, 1, 2]
    """
    global a, sequence
    if k == 1:
        return [0]
    a = [0] * k * n
    sequence = []
    gen(1, 1, k, n)
    return sequence


cdef gen(int t, int p, k, n):
    """
    The internal generation function. This should not be accessed by the
    user.
    """
    cdef int j
    if t > n:
        if n % p == 0:
            for j in range(1, p + 1):
                sequence.append(a[j])
    else:
        a[t] = a[t - p]
        gen(t + 1, p, k, n)
        for j in range((a[t - p] + 1), (k)):
            a[t] = j
            gen(t + 1, t, k, n)


def is_debruijn_sequence(seq, k, n):
    r"""
    Given a sequence of integer elements in `0, \ldots, k-1`, tests whether it
    corresponds to a De Bruijn sequence of parameters `k` and `n`.

    INPUT:

    - ``seq`` -- sequence of elements in `0, \ldots, k-1`

    - ``n``, ``k`` -- integers

    EXAMPLES::

        sage: from sage.combinat.debruijn_sequence import is_debruijn_sequence
        sage: s = DeBruijnSequences(2, 3).an_element()
        sage: is_debruijn_sequence(s, 2, 3)
        True
        sage: is_debruijn_sequence(s + [0], 2, 3)
        False
        sage: is_debruijn_sequence([1] + s[1:], 2, 3)
        False
    """

    if k == 1:
        return seq == [0]

    # The implementation is pretty straightforward.

    # The variable "current" is the integer representing the value of a
    # substring of length n. We iterate over all the possible substrings from
    # left to right, and keep track of the values met so far with the bitset
    # "seen".

    cdef int i
    cdef list s = seq
    cdef int nn = n
    cdef int kk = k

    cdef int k_p_n = kk ** nn

    cdef bitset_t seen

    # Checking if the length is correct
    if len(s) != kk ** nn:
        return False

    # Initializing the bitset
    bitset_init(seen, k_p_n)
    bitset_set_first_n(seen, 0)

    # We initialize "current" to correspond to the word formed by the (n-1) last elements
    cdef int current = 0

    for i in range(n - 1):
        current = kk * current + s[-n + i + 1]

    answer = True

    # Main loop, stopping if the same word has been met twice
    for i in s:
        current = (kk * current + i) % k_p_n

        if bitset_in(seen, current) or i < 0 or i >= k:
            answer = False
            break

        bitset_set(seen, current)

    bitset_free(seen)

    return answer


from sage.categories.finite_sets import FiniteSets
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ


class DeBruijnSequences(UniqueRepresentation, Parent):
    r"""
    Represent the De Bruijn sequences of given parameters `k` and `n`.

    A De Bruijn sequence of parameters `k` and `n` is defined as the shortest
    cyclic sequence that incorporates all substrings of length `n` a `k`-ary
    alphabet.

    This class can be used to generate the lexicographically smallest De Bruijn
    sequence, to count the number of existing De Bruijn sequences or to test
    whether a given sequence is De Bruijn.

    INPUT:

    - ``k`` -- a natural number to define arity; the letters used are the
      integers `0, \ldots, k-1`

    - ``n`` -- a natural number that defines the length of the substring

    EXAMPLES:

    Obtaining a De Bruijn sequence::

        sage: seq = DeBruijnSequences(2, 3).an_element()
        sage: seq
        [0, 0, 0, 1, 0, 1, 1, 1]

    Testing whether it is indeed one::

        sage: seq in DeBruijnSequences(2, 3)
        True

    The total number for these parameters::

        sage: DeBruijnSequences(2, 3).cardinality()
        2

    .. NOTE::

       This function only generates one De Bruijn sequence (the smallest
       lexicographically). Support for generating all possible ones may be
       added in the future.

    TESTS:

    Setting ``k`` to 1 will return 0:

    ::

        sage: DeBruijnSequences(1, 3).an_element()
        [0]

    Setting ``n`` to 1 will return the alphabet::

        sage: DeBruijnSequences(3, 1).an_element()
        [0, 1, 2]

    The test suite::

        sage: d = DeBruijnSequences(2, 3)
        sage: TestSuite(d).run()
    """
    def __init__(self, k, n):
        """
        Constructor.

        This checks the consistency of the given arguments.

        TESTS:

        Setting ``n`` or ``k`` to anything under 1 will return
        a :exc:`ValueError`::

            sage: DeBruijnSequences(3, 0).an_element()
            Traceback (most recent call last):
            ...
            ValueError: k and n cannot be under 1

        Setting ``n`` or ``k`` to any type except an integer will return a
        :exc:`TypeError`::

            sage: DeBruijnSequences(2.5, 3).an_element()
            Traceback (most recent call last):
            ...
            TypeError: k and n must be integers
        """
        Parent.__init__(self, category=FiniteSets())
        if n < 1 or k < 1:
            raise ValueError('k and n cannot be under 1')
        if (not isinstance(n, (Integer, int)) or
                not isinstance(k, (Integer, int))):
            raise TypeError('k and n must be integers')

        self.k = k
        self.n = n

    def _repr_(self):
        """
        Provides a string representation of the object's parameter.

        EXAMPLES::

            sage: repr(DeBruijnSequences(4, 50))
            'De Bruijn sequences with arity 4 and substring length 50'
        """
        return ("De Bruijn sequences with arity %s and substring length %s"
                % (self.k, self.n))

    def an_element(self):
        """
        Return the lexicographically smallest De Bruijn sequence with the given
        parameters.

        ALGORITHM:

        The algorithm is described in the book "Combinatorial Generation" by
        Frank Ruskey. This program is based on a Ruby implementation by Jonas
        Elfström, which is based on the C program by Joe Sadawa.

        EXAMPLES::

            sage: DeBruijnSequences(2, 3).an_element()
            [0, 0, 0, 1, 0, 1, 1, 1]
        """
        return debruijn_sequence(self.k, self.n)

    def __contains__(self, seq):
        r"""
        Test whether the given sequence is a De Bruijn sequence with
        the current object's parameters.

        INPUT:

        - ``seq`` -- a sequence of integers

        EXAMPLES::

           sage: Sequences =  DeBruijnSequences(2, 3)
           sage: Sequences.an_element() in Sequences
           True
        """
        return is_debruijn_sequence(seq, self.k, self.n)

    def cardinality(self):
        """
        Return the number of distinct De Bruijn sequences for the object's
        parameters.

        EXAMPLES::

            sage: DeBruijnSequences(2, 5).cardinality()
            2048

        ALGORITHM:

        The formula for cardinality is `k!^{k^{n-1}}/k^n` [Ros2002]_.
        """
        k = ZZ(self.k)
        n = ZZ(self.n)
        return (k.factorial() ** (k ** (n - 1))) // (k**n)
