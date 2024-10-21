# sage_setup: distribution = sagemath-categories
r"""
Bitsets

A Python interface to the fast bitsets in Sage.  Bitsets are fast
binary sets that store elements by toggling bits in an array of
numbers.  A bitset can store values between `0` and ``capacity - 1``,
inclusive (where ``capacity`` is finite, but arbitrary).  The storage cost is
linear in ``capacity``.

.. warning::

    This class is most likely to be useful as a way to store Cython
    bitsets in Python data structures, acting on them using the Cython
    inline functions.  If you want to use these classes for a Python
    set type, the Python ``set`` or ``frozenset`` data types may be
    faster.
"""

#*****************************************************************************
#       Copyright (C) 2009 Jason Grout <jason-sage@creativetrax.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.data_structures.bitset_base cimport *
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_NE, Py_GT, Py_GE


cdef class FrozenBitset:
    r"""
    A frozen bitset class which leverages inline Cython functions for
    creating and manipulating bitsets.

    A bitset can be thought of in two ways.  First, as a set of elements
    from the universe of the `n` natural numbers `0, 1, \dots, n-1` (where
    the capacity `n` can be specified), with typical set operations such as
    intersection, union, symmetric difference, etc.  Secondly, a bitset can
    be thought of as a binary vector with typical binary operations such as
    ``and``, ``or``, ``xor``, etc.  This class supports both interfaces.

    The interface in this class mirrors the interface in the ``frozenset``
    data type of Python. See the Python documentation on `set types
    <https://docs.python.org/library/stdtypes.html#set-types-set-frozenset>`_
    for more details on Python's ``set`` and ``frozenset`` classes.

    .. warning::

        This class is most likely to be useful as a way to store
        Cython bitsets in Python data structures, acting on them using
        the Cython inline functions.  If you want to use this class
        for a Python set type, the Python ``frozenset`` data type may
        be faster.

    INPUT:

    - ``iter`` -- initialization parameter (default: ``None``); valid inputs
      are:

      - :class:`Bitset` and :class:`FrozenBitset` -- If this is a
        :class:`Bitset` or :class:`FrozenBitset`, then it is copied

      - ``None`` -- if ``None``, then the bitset is set to the empty set

      - ``string`` -- if a nonempty string, then the bitset is initialized by
        including an element if the index of the string is ``1``. If the
        string is empty, then raise a :exc:`ValueError`.

      - ``iterable`` -- if an iterable, then it is assumed to contain a list of
        nonnegative integers and those integers are placed in the set

    - ``capacity`` -- (default: ``None``) the maximum capacity of the bitset.
      If this is not specified, then it is automatically calculated from the
      passed iterable.  It must be at least one.

    OUTPUT: none

    The string representation of a :class:`FrozenBitset` ``FB`` can be
    understood as follows. Let `B = b_0 b_1 b_2 \cdots b_k` be the string
    representation of the bitset ``FB``, where each `b_i \in \{0, 1\}`. We
    read the `b_i` from left to right. If `b_i = 1`, then the nonnegative
    integer `i` is in the bitset ``FB``. Similarly, if `b_i = 0`, then `i`
    is not in ``FB``. In other words, ``FB`` is a subset of
    `\{0, 1, 2, \dots, k\}` and the membership in ``FB`` of each `i` is
    determined by the binary value `b_i`.

    .. SEEALSO::

        - :class:`Bitset`

        - Python's `set types <https://docs.python.org/library/stdtypes.html#set-types-set-frozenset>`_

    EXAMPLES:

    The default bitset, which has capacity 1::

        sage: FrozenBitset()
        0
        sage: FrozenBitset(None)
        0

    Trying to create an empty bitset fails::

        sage: FrozenBitset([])
        Traceback (most recent call last):
        ...
        ValueError: Bitsets must not be empty
        sage: FrozenBitset(list())
        Traceback (most recent call last):
        ...
        ValueError: Bitsets must not be empty
        sage: FrozenBitset(())
        Traceback (most recent call last):
        ...
        ValueError: Bitsets must not be empty
        sage: FrozenBitset(tuple())
        Traceback (most recent call last):
        ...
        ValueError: Bitsets must not be empty
        sage: FrozenBitset("")
        Traceback (most recent call last):
        ...
        ValueError: Bitsets must not be empty

    We can create the all-zero bitset as follows::

        sage: FrozenBitset(capacity=10)
        0000000000
        sage: FrozenBitset([], capacity=10)
        0000000000

    We can initialize a :class:`FrozenBitset` with a :class:`Bitset` or
    another :class:`FrozenBitset`, and compare them for equality. As they
    are logically the same bitset, the equality test should return ``True``.
    Furthermore, each bitset is a subset of the other. ::

        sage: def bitcmp(a, b, c):  # custom function for comparing bitsets
        ....:     print(a == b == c)
        ....:     print((a <= b, b <= c, a <= c))
        ....:     print((a >= b, b >= c, a >= c))
        ....:     print((a != b, b != c, a != c))
        sage: a = Bitset("1010110"); b = FrozenBitset(a); c = FrozenBitset(b)
        sage: a; b; c
        1010110
        1010110
        1010110
        sage: a < b, b < c, a < c
        (False, False, False)
        sage: a > b, b > c, a > c
        (False, False, False)
        sage: bitcmp(a, b, c)
        True
        (True, True, True)
        (True, True, True)
        (False, False, False)

    Try a random bitset::

        sage: a = Bitset(randint(0, 1) for n in range(1, randint(1, 10^4)))
        sage: b = FrozenBitset(a); c = FrozenBitset(b)
        sage: bitcmp(a, b, c)
        True
        (True, True, True)
        (True, True, True)
        (False, False, False)

    A bitset with a hard-coded bitstring::

        sage: FrozenBitset('101')
        101

    For a string, only those positions with ``1`` would be initialized to ``1``
    in the corresponding position in the bitset. All other characters in the
    string, including ``0``, are set to ``0`` in the resulting bitset. ::

        sage: FrozenBitset('a')
        0
        sage: FrozenBitset('abc')
        000
        sage: FrozenBitset('abc1')
        0001
        sage: FrozenBitset('0abc1')
        00001
        sage: FrozenBitset('0abc10')
        000010
        sage: FrozenBitset('0a*c10')
        000010

    Represent the first 10 primes as a bitset. The primes are stored as a
    list and as a tuple. We then recover the primes from its bitset
    representation, and query the bitset for its length (how many elements
    it contains) and whether an element is in the bitset. Note that the
    length of a bitset is different from its capacity. The length counts
    the number of elements currently in the bitset, while the capacity
    is the number of elements that the bitset can hold. ::

        sage: p = primes_first_n(10); p                                                 # needs sage.libs.pari
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        sage: tuple(p)                                                                  # needs sage.libs.pari
        (2, 3, 5, 7, 11, 13, 17, 19, 23, 29)
        sage: F = FrozenBitset(p); F; FrozenBitset(tuple(p))                            # needs sage.libs.pari
        001101010001010001010001000001
        001101010001010001010001000001

    Recover the primes from the bitset::

        sage: for b in F:                                                               # needs sage.libs.pari
        ....:     print(b)
        2
        3
        ...
        29
        sage: list(F)                                                                   # needs sage.libs.pari
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

    Query the bitset::

        sage: # needs sage.libs.pari
        sage: len(F)
        10
        sage: len(list(F))
        10
        sage: F.capacity()
        30
        sage: s = str(F); len(s)
        30
        sage: 2 in F
        True
        sage: 1 in F
        False

    A random iterable, with all duplicate elements removed::

        sage: L = [randint(0, 100) for n in range(1, randint(1, 10^4))]
        sage: FrozenBitset(L) == FrozenBitset(list(set(L)))
        True
        sage: FrozenBitset(tuple(L)) == FrozenBitset(tuple(set(L)))
        True

    TESTS:

    Loading and dumping objects::

        sage: a = FrozenBitset('1101')
        sage: loads(dumps(a)) == a
        True
        sage: a = FrozenBitset('1101' * 64)
        sage: loads(dumps(a)) == a
        True

    If ``iter`` is a nonempty string and ``capacity`` is specified, then
    ``capacity`` must match the number of elements in ``iter``::

        sage: FrozenBitset("110110", capacity=3)
        Traceback (most recent call last):
        ...
        ValueError: bitset capacity does not match passed string
        sage: FrozenBitset("110110", capacity=100)
        Traceback (most recent call last):
        ...
        ValueError: bitset capacity does not match passed string

    The parameter ``capacity`` must be positive::

        sage: FrozenBitset("110110", capacity=0)
        Traceback (most recent call last):
        ...
        ValueError: bitset capacity must be greater than 0
        sage: FrozenBitset("110110", capacity=-2)
        Traceback (most recent call last):
        ...
        OverflowError: can...t convert negative value to mp_bitcnt_t
    """
    def __cinit__(self, iter=None, capacity=None):
        """
        Allocate the bitset, which is initially empty.

        See the class documentation of :class:`FrozenBitset` for details
        on the parameters.

        EXAMPLES::

            sage: FrozenBitset('1101')
            1101
            sage: FrozenBitset('1101' * 32)
            11011101110111011101110111011101110111011101110111011101110111011101110111011101110111011101110111011101110111011101110111011101
        """
        if capacity is None:
            bitset_init(self._bitset, 1)
        else:
            bitset_init(self._bitset, <mp_bitcnt_t> capacity)

    def __dealloc__(self):
        """
        Deallocate the C bitset data structure.

        EXAMPLES::

            sage: a = FrozenBitset('11010')
            sage: del a
            sage: b = FrozenBitset('11010' * 64)
            sage: del b
        """
        bitset_free(self._bitset)

    def __init__(self, iter=None, capacity=None):
        """
        Initialize a bitset.

        See the class documentation of ``FrozenBitset`` for details on the
        parameters.

        EXAMPLES::

            sage: FrozenBitset(capacity=3)
            000
            sage: FrozenBitset('11011')
            11011
            sage: FrozenBitset('110' * 32)
            110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110
            sage: FrozenBitset([0,3,2])
            1011
            sage: FrozenBitset(set([0,3,2]))
            1011
            sage: FrozenBitset(FrozenBitset('11011'))
            11011
            sage: FrozenBitset([0,3,2], capacity=10)
            1011000000
            sage: FrozenBitset([i for i in range(100) if i % 2 == 0])
            101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101

        TESTS:

        The capacity must be at least one::

            sage: FrozenBitset()
            0

        If capacity is specified, it must match up with the
        initialization items::

            sage: FrozenBitset('10', capacity=3)
            Traceback (most recent call last):
            ...
            ValueError: bitset capacity does not match passed string
            sage: FrozenBitset([0,3,2], capacity=2)
            Traceback (most recent call last):
            ...
            ValueError: bitset capacity does not allow storing the passed iterable
            sage: FrozenBitset(FrozenBitset('1101'), capacity=2)
            Traceback (most recent call last):
            ...
            ValueError: bitset capacity does not match passed bitset
            sage: FrozenBitset(FrozenBitset('1101'), capacity=5)
            Traceback (most recent call last):
            ...
            ValueError: bitset capacity does not match passed bitset
        """
        cdef unsigned long n
        cdef FrozenBitset b
        if iter is not None and not isinstance(iter, (str, tuple, list, FrozenBitset, Bitset)):
            iter = list(iter)

        if iter is None:
            pass  # Leave bitset empty
        elif isinstance(iter, (FrozenBitset, Bitset)):
            b = iter
            if capacity is None:
                bitset_realloc(self._bitset, b._bitset.size)
            elif self._bitset.size != b._bitset.size:
                raise ValueError("bitset capacity does not match passed bitset")
            bitset_copy(self._bitset, b._bitset)
        elif isinstance(iter, str):
            if len(iter) == 0:
                raise ValueError("Bitsets must not be empty")
            if capacity is None:
                bitset_realloc(self._bitset, len(iter))
            elif self._bitset.size != <mp_bitcnt_t>len(iter):
                raise ValueError("bitset capacity does not match passed string")
            bitset_from_str(self._bitset, iter)
        else:  # an iterable
            iter = list(iter)
            if iter:
                need_capacity = max(iter) + 1
            else:
                need_capacity = 0
            if capacity is None:
                if need_capacity == 0:
                    raise ValueError("Bitsets must not be empty")
                bitset_realloc(self._bitset, need_capacity)
            elif self._bitset.size < need_capacity:
                raise ValueError("bitset capacity does not allow storing the passed iterable")
            bitset_clear(self._bitset)
            for n in iter:
                bitset_add(self._bitset, n)

    cdef FrozenBitset _new(self, long int capacity):
        r"""
        Return an object of the same type as ``self``, initialized with a
        bitset of capacity ``capacity``.
        """
        cdef FrozenBitset b
        b = FrozenBitset.__new__(FrozenBitset, None, capacity)
        return b

    def __getstate__(self):
        """
        Return the current state of the object as a bytes string.

        EXAMPLES::

            sage: FrozenBitset('1101').__getstate__() == b'1101'
            True
            sage: FrozenBitset('110'*32).__getstate__() == (b'110' * 32)
            True
        """
        return bytes(self)

    def __setstate__(self, state):
        """
        Set the state of the object from the string state.

        EXAMPLES::

            sage: a = FrozenBitset()
            sage: a.__setstate__(b'1101')
            sage: a
            1101
            sage: a.__setstate__(b'110'*32)
            sage: a
            110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110
        """
        bitset_realloc(self._bitset, len(state))
        bitset_from_char(self._bitset, state)

    def __iter__(self):
        """
        Return an iterator over ``self``.

        EXAMPLES::

            sage: list(FrozenBitset('11011'))
            [0, 1, 3, 4]
            sage: list(FrozenBitset('00001' * 20))
            [4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 79, 84, 89, 94, 99]
            sage: set(FrozenBitset('11011'))
            {0, 1, 3, 4}
        """
        return iter(bitset_list(self._bitset))

    def __reversed__(self):
        """
        Return an iterator over ``self``, starting with the largest element.

        EXAMPLES::

            sage: list(reversed(FrozenBitset('11011')))
            [4, 3, 1, 0]
            sage: list(reversed(FrozenBitset('00001' * 20)))
            [99, 94, 89, 84, 79, 74, 69, 64, 59, 54, 49, 44, 39, 34, 29, 24, 19, 14, 9, 4]
        """
        return reversed(bitset_list(self._bitset))

    cpdef FrozenBitset _larger_capacity_(self, long capacity):
        """
        Return a copy of ``self`` where the bitset has the maximum of the
        current capacity and the capacity passed.  If no resizing is needed,
        return ``self``.

        This function is mainly used to satisfy the assumption of the
        underlying bitset functions that all bitsets are of the same
        capacity.

        INPUT:

        - ``capacity`` -- the underlying bitset of the returned bitset
          will have this capacity if it is bigger than the current
          capacity.

        EXAMPLES::

            sage: a = FrozenBitset('11010')
            sage: a.capacity()
            5
            sage: a._larger_capacity_(4) is a
            True
            sage: a._larger_capacity_(5) is a
            True
            sage: b = a._larger_capacity_(6)
            sage: b
            110100
            sage: b is a
            False
            sage: b.capacity()
            6
            sage: c = a._larger_capacity_(98)
            sage: c
            11010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: c.capacity()
            98
        """
        cdef FrozenBitset temp
        if self._bitset.size >= <mp_bitcnt_t>capacity:
            return self
        else:
            temp = self._new(self._bitset.size)
            bitset_copy(temp._bitset, self._bitset)
            bitset_realloc(temp._bitset, capacity)
            return temp

    cpdef long capacity(self) noexcept:
        """
        Return the size of the underlying bitset.

        The maximum value that can be stored in the current underlying
        bitset is ``self.capacity() - 1``.

        EXAMPLES::

            sage: FrozenBitset('11000').capacity()
            5
            sage: FrozenBitset('110' * 32).capacity()
            96
            sage: FrozenBitset(range(20), capacity=450).capacity()
            450
        """
        return self._bitset.size

    def __hash__(self):
        """
        Return a hash value for a bitset.

        EXAMPLES::

            sage: hash(FrozenBitset(capacity=5))
            0
            sage: hash(FrozenBitset('10110'))
            13
            sage: hash(FrozenBitset('10110' + '0' * 120, capacity=125))
            13
        """
        cdef long hash
        hash = bitset_hash(self._bitset)
        if hash == -1:
            return 0
        else:
            return hash

    cpdef bint isempty(self) noexcept:
        """
        Test if the bitset is empty.

        OUTPUT: boolean

        EXAMPLES::

            sage: FrozenBitset().isempty()
            True
            sage: FrozenBitset([1]).isempty()
            False
            sage: FrozenBitset([], capacity=110).isempty()
            True
            sage: FrozenBitset(range(99)).isempty()
            False
        """
        return bitset_isempty(self._bitset)

    def __richcmp__(FrozenBitset self, FrozenBitset other not None, int op):
        """
        Implement comparisons, using the Cython richcmp convention.
        Comparison is done by inclusion. That is, a set ``A`` is less than
        another set ``B``, written ``A < B``, if ``A`` is a subset of ``B``.

        EXAMPLES::

            sage: FrozenBitset('11') < FrozenBitset('101')
            False
            sage: FrozenBitset('11') <= FrozenBitset('110')
            True
            sage: FrozenBitset('11') != FrozenBitset('10')
            True
            sage: FrozenBitset('11') == FrozenBitset('10')
            False
            sage: FrozenBitset('11') == FrozenBitset('110')
            True
            sage: FrozenBitset('11') > FrozenBitset('10')
            True
            sage: FrozenBitset('11') >= FrozenBitset('10')
            True
            sage: FrozenBitset('11') < FrozenBitset('110' * 32)
            True

        TESTS:

        When performing comparison, ``other`` cannot be ``None``. ::

            sage: F = FrozenBitset('11')
            sage: F < None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: F <= None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: F > None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: F >= None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: F == None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: F != None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None < F
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None <= F
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None > F
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None >= F
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None == F
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None != F
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
        """
        cdef FrozenBitset left, right

        if self._bitset.size == other._bitset.size:
            left = self
            right = other
        elif self._bitset.size < other._bitset.size:
            left = self._larger_capacity_(other._bitset.size)
            right = other
        else:
            left = self
            right = other._larger_capacity_(self._bitset.size)

        if op == Py_EQ:
            return bitset_eq(left._bitset, right._bitset)
        elif op == Py_NE:
            return not bitset_eq(left._bitset, right._bitset)
        elif op == Py_LT:
            return bitset_issubset(left._bitset, right._bitset) and not bitset_eq(left._bitset, right._bitset)
        elif op == Py_LE:
            return bitset_issubset(left._bitset, right._bitset)
        elif op == Py_GT:
            return bitset_issuperset(left._bitset, right._bitset) and not bitset_eq(left._bitset, right._bitset)
        elif op == Py_GE:
            return bitset_issuperset(left._bitset, right._bitset)

    cpdef bint issubset(self, FrozenBitset other) except -1:
        """
        Test to see if ``self`` is a subset of ``other``.

        EXAMPLES::

            sage: FrozenBitset('11').issubset(FrozenBitset('01'))
            False
            sage: FrozenBitset('01').issubset(FrozenBitset('11'))
            True
            sage: FrozenBitset('01').issubset(FrozenBitset('01' * 45))
            True

        TESTS::

            sage: FrozenBitset('11').issubset(None)
            Traceback (most recent call last):
            ...
            ValueError: other cannot be None
        """
        if other is None:
            raise ValueError("other cannot be None")
        cdef FrozenBitset left, right
        if self._bitset.size <= other._bitset.size:
            left = self
            right = other
        else:
            left = self
            right = other._larger_capacity_(self._bitset.size)

        # Assumes ``left.size <= right.size``.
        return bitset_issubset(left._bitset, right._bitset)

    cpdef bint issuperset(self, FrozenBitset other) except -1:
        """
        Test to see if ``self`` is a superset of ``other``.

        EXAMPLES::

            sage: FrozenBitset('11').issuperset(FrozenBitset('01'))
            True
            sage: FrozenBitset('01').issuperset(FrozenBitset('11'))
            False
            sage: FrozenBitset('01').issuperset(FrozenBitset('10' * 45))
            False

        TESTS::

            sage: FrozenBitset('11').issuperset(None)
            Traceback (most recent call last):
            ...
            ValueError: other cannot be None
        """
        if other is None:
            raise ValueError("other cannot be None")
        cdef FrozenBitset left, right
        if self._bitset.size >= other._bitset.size:
            left = self
            right = other
        else:
            left = self._larger_capacity_(other._bitset.size)
            right = other

        # Assumes ``left.size >= right.size``.
        return bitset_issuperset(left._bitset, right._bitset)

    cpdef bint isdisjoint(self, FrozenBitset other) except -1:
        """
        Test to see if ``self`` is disjoint from ``other``.

        EXAMPLES::

            sage: FrozenBitset('11').isdisjoint(FrozenBitset('01'))
            False
            sage: FrozenBitset('01').isdisjoint(FrozenBitset('001'))
            True
            sage: FrozenBitset('00101').isdisjoint(FrozenBitset('110' * 35))
            False

        TESTS::

            sage: FrozenBitset('11').isdisjoint(None)
            Traceback (most recent call last):
            ...
            ValueError: other cannot be None
        """
        cdef bint retval
        if other is None:
            raise ValueError("other cannot be None")
        cdef FrozenBitset left, right

        if self._bitset.size <= other._bitset.size:
            left = self
            right = other
        else:
            left = other
            right = self

        # Assumes ``left.size <= right.size``.
        return bitset_are_disjoint(left._bitset, right._bitset)

    def __contains__(self, unsigned long n):
        """
        Test to see if ``n`` is in ``self``.

        EXAMPLES::

            sage: 0 in FrozenBitset([0,1])
            True
            sage: 0 in FrozenBitset([1,2])
            False
            sage: 10 in FrozenBitset([0,1])
            False
            sage: 121 in FrozenBitset('110' * 50)
            True
            sage: 122 in FrozenBitset('110' * 50)
            False

        TESTS::

            sage: None in FrozenBitset([0,1])
            Traceback (most recent call last):
            ...
            TypeError: an integer is required
        """
        if n < self._bitset.size:
            return bitset_in(self._bitset, n)
        else:
            return False

    def __len__(self):
        """
        Return the number of elements in the bitset (which may be
        different from the capacity of the bitset).

        EXAMPLES::

            sage: len(FrozenBitset([0,1], capacity=10))
            2
            sage: len(FrozenBitset(range(98)))
            98
        """
        return bitset_len(self._bitset)

    def __str__(self):
        """
        Return a string representing the bitset as a binary vector.

        EXAMPLES::

            sage: a = FrozenBitset('10110')
            sage: str(a)
            '10110'
            sage: str(FrozenBitset('110' * 32))
            '110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110'
        """
        return bitset_string(self._bitset)

    def __bytes__(self):
        """
        Return a bytes object representing the bitset as a binary vector.

        EXAMPLES::

            sage: a = FrozenBitset('10110')
            sage: bytes(a) == b'10110'
            True
            sage: bytes(FrozenBitset('110' * 32)) == b'110' * 32
            True
        """
        return bitset_bytes(self._bitset)

    def __repr__(self):
        """
        Return a string representing the bitset as a binary vector.

        EXAMPLES::

            sage: a = FrozenBitset('10110')
            sage: repr(a)
            '10110'
            sage: repr(FrozenBitset('110' * 32))
            '110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110'
        """
        return str(self)

    cpdef _union(self, FrozenBitset other):
        """
        Return the union of ``self`` and ``other``.

        In order to get a Cython "union" function, we have to use the
        underscore since "union" is a C keyword.

        EXAMPLES::

            sage: FrozenBitset('10101')._union(FrozenBitset('11100'))
            11101
            sage: FrozenBitset('10101' * 10)._union(FrozenBitset('01010' * 10))
            11111111111111111111111111111111111111111111111111

        TESTS::

            sage: set(FrozenBitset('10101' * 10)._union(FrozenBitset('01010' * 10))) == set(FrozenBitset('10101' * 10)).union(FrozenBitset('01010' * 10))
            True
            sage: set(FrozenBitset('10101')._union(FrozenBitset('01010' * 20))) == set(FrozenBitset('10101')).union(FrozenBitset('01010' * 20))
            True
            sage: set(FrozenBitset('10101' * 20)._union(FrozenBitset('01010'))) == set(FrozenBitset('10101' * 20)).union(FrozenBitset('01010'))
            True
            sage: FrozenBitset('10101' * 10)._union(None)
            Traceback (most recent call last):
            ...
            ValueError: other cannot be None
        """
        if other is None:
            raise ValueError("other cannot be None")
        cdef FrozenBitset temp, smaller, larger
        if self._bitset.size <= other._bitset.size:
            smaller = self
            larger = other
        else:
            smaller = other
            larger = self

        temp = self._new(smaller._bitset.size)
        bitset_copy(temp._bitset, smaller._bitset)
        bitset_realloc(temp._bitset, larger._bitset.size)
        bitset_union(temp._bitset, temp._bitset, larger._bitset)
        return temp

    def union(self, FrozenBitset other):
        """
        Return the union of ``self`` and ``other``.

        EXAMPLES::

            sage: FrozenBitset('10101').union(FrozenBitset('11100'))
            11101
            sage: FrozenBitset('10101' * 10).union(FrozenBitset('01010' * 10))
            11111111111111111111111111111111111111111111111111

        TESTS::

            sage: set(FrozenBitset('10101' * 10).union(FrozenBitset('01010' * 10))) == set(FrozenBitset('10101' * 10)).union(FrozenBitset('01010' * 10))
            True
            sage: set(FrozenBitset('10101').union(FrozenBitset('01010' * 20))) == set(FrozenBitset('10101')).union(FrozenBitset('01010' * 20))
            True
            sage: set(FrozenBitset('10101' * 20).union(FrozenBitset('01010'))) == set(FrozenBitset('10101' * 20)).union(FrozenBitset('01010'))
            True
            sage: FrozenBitset('10101' * 10).union(None)
            Traceback (most recent call last):
            ...
            ValueError: other cannot be None
        """
        return self._union(other)

    def __or__(self, FrozenBitset other not None):
        """
        Return the union of ``self`` and ``other``.

        EXAMPLES::

            sage: FrozenBitset('10101') | FrozenBitset('11100')
            11101
            sage: FrozenBitset('10101' * 10) | FrozenBitset('01010' * 10)
            11111111111111111111111111111111111111111111111111

       TESTS::

            sage: set(FrozenBitset('10101' * 10) | FrozenBitset('01010' * 10)) == set(FrozenBitset('10101' * 10)) | set(FrozenBitset('01010' * 10))
            True
            sage: set(FrozenBitset('10101') | FrozenBitset('01010' * 20)) == set(FrozenBitset('10101')) | set(FrozenBitset('01010' * 20))
            True
            sage: set(FrozenBitset('10101' * 20) | FrozenBitset('01010')) == set(FrozenBitset('10101' * 20)) | set(FrozenBitset('01010'))
            True
            sage: FrozenBitset('10101') | None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None | FrozenBitset('10101')
            Traceback (most recent call last):
            ...
            AttributeError: 'NoneType' object has no attribute '_union'...
        """
        return self._union(other)

    cpdef intersection(self, FrozenBitset other):
        """
        Return the intersection of ``self`` and ``other``.

        EXAMPLES::

            sage: FrozenBitset('10101').intersection(FrozenBitset('11100'))
            10100
            sage: FrozenBitset('11111' * 10).intersection(FrozenBitset('010101' * 10))
            010101010101010101010101010101010101010101010101010000000000

        TESTS::

            sage: set(FrozenBitset('11111' * 10).intersection(FrozenBitset('010101' * 10))) == set(FrozenBitset('11111' * 10)).intersection(FrozenBitset('010101' * 10))
            True
            sage: set(FrozenBitset('1' * 5).intersection(FrozenBitset('01010' * 20))) == set(FrozenBitset('1' * 5)).intersection(FrozenBitset('01010' * 20))
            True
            sage: set(FrozenBitset('10101' * 20).intersection(FrozenBitset('1' * 5))) == set(FrozenBitset('10101' * 20)).intersection(FrozenBitset('1' * 5))
            True
            sage: FrozenBitset("101011").intersection(None)
            Traceback (most recent call last):
            ...
            ValueError: other cannot be None
        """
        if other is None:
            raise ValueError("other cannot be None")
        cdef FrozenBitset temp, smaller, larger
        if self._bitset.size <= other._bitset.size:
            smaller = self
            larger = other
        else:
            smaller = other
            larger = self

        temp = self._new(smaller._bitset.size)
        bitset_copy(temp._bitset, smaller._bitset)
        bitset_realloc(temp._bitset, larger._bitset.size)
        bitset_intersection(temp._bitset, temp._bitset, larger._bitset)
        return temp

    def __and__(self, FrozenBitset other not None):
        """
        Return the intersection of ``self`` and ``other``.

        EXAMPLES::

            sage: FrozenBitset('10101') & FrozenBitset('11100')
            10100
            sage: FrozenBitset('11111' * 10) & FrozenBitset('010101' * 10)
            010101010101010101010101010101010101010101010101010000000000

        TESTS::

            sage: set(FrozenBitset('11111' * 10) & FrozenBitset('010101' * 10)) == set(FrozenBitset('11111' * 10)) & set(FrozenBitset('010101' * 10))
            True
            sage: set(FrozenBitset('1' * 5) & FrozenBitset('01010' * 20)) == set(FrozenBitset('1' * 5)) & set(FrozenBitset('01010' * 20))
            True
            sage: set(FrozenBitset('10101' * 20) & FrozenBitset('1' * 5)) == set(FrozenBitset('10101' * 20)) & set(FrozenBitset('1' * 5))
            True
            sage: FrozenBitset("101011") & None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None & FrozenBitset("101011")
            Traceback (most recent call last):
            ...
            AttributeError: 'NoneType' object has no attribute 'intersection'...
        """
        return self.intersection(other)

    cpdef difference(self, FrozenBitset other):
        """
        Return the difference of ``self`` and ``other``.

        EXAMPLES::

            sage: FrozenBitset('10101').difference(FrozenBitset('11100'))
            00001
            sage: FrozenBitset('11111' * 10).difference(FrozenBitset('010101' * 10))
            101010101010101010101010101010101010101010101010100000000000

        TESTS::

            sage: set(FrozenBitset('11111' * 10).difference(FrozenBitset('010101' * 10))) == set(FrozenBitset('11111' * 10)).difference(FrozenBitset('010101' * 10))
            True
            sage: set(FrozenBitset('1' * 5).difference(FrozenBitset('01010' * 20))) == set(FrozenBitset('1' * 5)).difference(FrozenBitset('01010' * 20))
            True
            sage: set(FrozenBitset('10101' * 20).difference(FrozenBitset('1' * 5))) == set(FrozenBitset('10101' * 20)).difference(FrozenBitset('1' * 5))
            True
            sage: FrozenBitset('10101').difference(None)
            Traceback (most recent call last):
            ...
            ValueError: other cannot be None
        """
        if other is None:
            raise ValueError("other cannot be None")
        cdef FrozenBitset temp = self._new(self._bitset.size)
        bitset_copy(temp._bitset, self._bitset)

        if temp._bitset.size == other._bitset.size:
            bitset_difference(temp._bitset, temp._bitset, other._bitset)
        elif temp._bitset.size < other._bitset.size:
            bitset_realloc(temp._bitset, other._bitset.size)
            bitset_difference(temp._bitset, temp._bitset, other._bitset)
        else:
            bitset_difference(temp._bitset, temp._bitset, other._larger_capacity_(temp._bitset.size)._bitset)

        return temp

    def __sub__(self, FrozenBitset other not None):
        """
        Return the difference of ``self`` and ``other``.

        EXAMPLES::

            sage: FrozenBitset('10101') - FrozenBitset('11100')
            00001
            sage: FrozenBitset('11111' * 10)-FrozenBitset('010101' * 10)
            101010101010101010101010101010101010101010101010100000000000

        TESTS::

            sage: set(FrozenBitset('11111' * 10)-FrozenBitset('010101' * 10)) == set(FrozenBitset('11111' * 10))-set(FrozenBitset('010101' * 10))
            True
            sage: set(FrozenBitset('1' * 5)-FrozenBitset('01010' * 20)) == set(FrozenBitset('1' * 5))-set(FrozenBitset('01010' * 20))
            True
            sage: set(FrozenBitset('10101' * 20)-FrozenBitset('1' * 5)) == set(FrozenBitset('10101' * 20))-set(FrozenBitset('1' * 5))
            True
            sage: FrozenBitset('10101') - None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None - FrozenBitset('10101')
            Traceback (most recent call last):
            ...
            AttributeError: 'NoneType' object has no attribute 'difference'...
        """
        return self.difference(other)

    cpdef symmetric_difference(self, FrozenBitset other):
        """
        Return the symmetric difference of ``self`` and ``other``.

        EXAMPLES::

            sage: FrozenBitset('10101').symmetric_difference(FrozenBitset('11100'))
            01001
            sage: FrozenBitset('11111' * 10).symmetric_difference(FrozenBitset('010101' * 10))
            101010101010101010101010101010101010101010101010100101010101

        TESTS::

            sage: set(FrozenBitset('11111' * 10).symmetric_difference(FrozenBitset('010101' * 10))) == set(FrozenBitset('11111' * 10)).symmetric_difference(FrozenBitset('010101' * 10))
            True
            sage: set(FrozenBitset('1' * 5).symmetric_difference(FrozenBitset('01010' * 20))) == set(FrozenBitset('1' * 5)).symmetric_difference(FrozenBitset('01010' * 20))
            True
            sage: set(FrozenBitset('10101' * 20).symmetric_difference(FrozenBitset('1' * 5))) == set(FrozenBitset('10101' * 20)).symmetric_difference(FrozenBitset('1' * 5))
            True
            sage: FrozenBitset('11111' * 10).symmetric_difference(None)
            Traceback (most recent call last):
            ...
            ValueError: other cannot be None
        """
        if other is None:
            raise ValueError("other cannot be None")
        cdef FrozenBitset temp, smaller, larger
        if self._bitset.size <= other._bitset.size:
            smaller = self
            larger = other
        else:
            smaller = other
            larger = self

        temp = self._new(smaller._bitset.size)
        bitset_copy(temp._bitset, smaller._bitset)
        bitset_realloc(temp._bitset, larger._bitset.size)
        bitset_symmetric_difference(temp._bitset, temp._bitset, larger._bitset)
        return temp

    def __xor__(self, FrozenBitset other not None):
        """
        Return the symmetric difference of ``self`` and ``other``.

        Note that because of the Sage preprocessor, in Sage, ``^^`` is the
        exclusive or, rather than ``^``.

        EXAMPLES::

            sage: FrozenBitset('10101') ^^ FrozenBitset('11100')
            01001
            sage: FrozenBitset('11111' * 10) ^^ FrozenBitset('010101' * 10)
            101010101010101010101010101010101010101010101010100101010101

        TESTS::

            sage: set(FrozenBitset('11111' * 10) ^^ FrozenBitset('010101' * 10)) == set(FrozenBitset('11111' * 10)) ^^ set(FrozenBitset('010101' * 10))
            True
            sage: set(FrozenBitset('1' * 5) ^^ FrozenBitset('01010' * 20)) == set(FrozenBitset('1' * 5)) ^^ set(FrozenBitset('01010' * 20))
            True
            sage: set(FrozenBitset('10101' * 20) ^^ FrozenBitset('1' * 5)) == set(FrozenBitset('10101' * 20)) ^^ set(FrozenBitset('1' * 5))
            True
            sage: FrozenBitset('11111' * 10) ^^ None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None ^^ FrozenBitset('11111' * 10)
            Traceback (most recent call last):
            ...
            AttributeError: 'NoneType' object has no attribute 'symmetric_difference'...
        """
        return self.symmetric_difference(other)

    cpdef complement(self):
        """
        Return the complement of ``self``.

        EXAMPLES::

            sage: ~FrozenBitset('10101')
            01010
            sage: ~FrozenBitset('11111'*10)
            00000000000000000000000000000000000000000000000000
            sage: x = FrozenBitset('10'*40)
            sage: x == ~x
            False
            sage: x == ~~x
            True
            sage: x|(~x) == FrozenBitset('11'*40)
            True
            sage: ~x == FrozenBitset('01'*40)
            True
        """
        cdef FrozenBitset temp = self._new(self._bitset.size)
        bitset_complement(temp._bitset, self._bitset)
        return temp

    def __invert__(self):
        """
        Return the complement of ``self``.

        EXAMPLES::

            sage: ~FrozenBitset('10101')
            01010
            sage: ~FrozenBitset('11111'*10)
            00000000000000000000000000000000000000000000000000
            sage: x = FrozenBitset('10'*40)
            sage: x == ~x
            False
            sage: x == ~~x
            True
            sage: x|(~x) == FrozenBitset('11'*40)
            True
            sage: ~x == FrozenBitset('01'*40)
            True
        """
        return self.complement()

    cpdef  __copy__(self):
        """
        Return ``self`` (since ``self`` is immutable).

        EXAMPLES::

            sage: a = FrozenBitset('10101')
            sage: from copy import copy
            sage: b = copy(a)
            sage: b is a
            True
            sage: c = FrozenBitset('1010' * 32)
            sage: d = copy(c)
            sage: d is c
            True
        """
        return self

cdef class Bitset(FrozenBitset):
    r"""
    A bitset class which leverages inline Cython functions for creating
    and manipulating bitsets. See the class documentation of
    :class:`FrozenBitset` for details on the parameters of the constructor
    and how to interpret the string representation of a :class:`Bitset`.

    A bitset can be thought of in two ways.  First, as a set of elements
    from the universe of the `n` natural numbers `0, 1, \dots, n-1` (where
    the capacity `n` can be specified), with typical set operations such as
    intersection, union, symmetric difference, etc.  Secondly, a bitset can
    be thought of as a binary vector with typical binary operations such as
    ``and``, ``or``, ``xor``, etc.  This class supports both interfaces.

    The interface in this class mirrors the interface in the ``set``
    data type of Python.

    .. warning::

        This class is most likely to be useful as a way to store
        Cython bitsets in Python data structures, acting on them using
        the Cython inline functions.  If you want to use this class
        for a Python set type, the Python ``set`` data type may be
        faster.

    .. SEEALSO::

        - :class:`FrozenBitset`
        - Python's `set types <https://docs.python.org/library/stdtypes.html#set-types-set-frozenset>`_

    EXAMPLES::

        sage: a = Bitset('1101')
        sage: loads(dumps(a)) == a
        True
        sage: a = Bitset('1101' * 32)
        sage: loads(dumps(a)) == a
        True
    """

    cpdef __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: a = Bitset('10101')
            sage: from copy import copy
            sage: b = copy(a)
            sage: b is a
            False
            sage: b == a
            True
            sage: c = Bitset('1010' * 32)
            sage: d = copy(c)
            sage: d is c
            False
            sage: d == c
            True
        """
        cdef FrozenBitset temp = self._new(self._bitset.size)
        bitset_copy(temp._bitset, self._bitset)
        return temp

    def __hash__(self):
        """
        Raise an error, since mutable ``Bitset``s are not hashable.

        EXAMPLES::

            sage: hash(Bitset('110'))
            Traceback (most recent call last):
            ...
            TypeError: Bitset objects are unhashable; use FrozenBitset
        """
        raise TypeError("Bitset objects are unhashable; use FrozenBitset")

    def __richcmp__(FrozenBitset self, FrozenBitset other not None, int op):
        """
        Implement comparisons, using the Cython richcmp convention.
        Comparison is done by inclusion. That is, a set ``A`` is less than
        another set ``B``, written ``A < B``, if ``A`` is a subset of ``B``.

        EXAMPLES::

            sage: Bitset('11') < Bitset('101')
            False
            sage: Bitset('11') <= Bitset('110')
            True
            sage: Bitset('11') != Bitset('10')
            True
            sage: Bitset('11') == Bitset('10')
            False
            sage: Bitset('11') == Bitset('110')
            True
            sage: Bitset('11') > Bitset('10')
            True
            sage: Bitset('11') >= Bitset('10')
            True
            sage: FrozenBitset('11') < FrozenBitset('110' * 32)
            True

        TESTS:

        During comparison, ``other`` cannot be ``None``. ::

            sage: Bitset('11') < None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: Bitset('11') <= None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: Bitset('11') > None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: Bitset('11') >= None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: Bitset('11') == None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: Bitset('11') != None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None < Bitset('11')
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None <= Bitset('11')
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None > Bitset('11')
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None >= Bitset('11')
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None == Bitset('11')
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
            sage: None != Bitset('11')
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
        """
        cdef FrozenBitset left, right

        if self._bitset.size == other._bitset.size:
            left = self
            right = other
        elif self._bitset.size < other._bitset.size:
            left = self._larger_capacity_(other._bitset.size)
            right = other
        else:
            left = self
            right = other._larger_capacity_(self._bitset.size)

        if op == Py_EQ:
            return bitset_eq(left._bitset, right._bitset)
        elif op == Py_NE:
            return not bitset_eq(left._bitset, right._bitset)
        elif op == Py_LT:
            return bitset_issubset(left._bitset, right._bitset) and not bitset_eq(left._bitset, right._bitset)
        elif op == Py_LE:
            return bitset_issubset(left._bitset, right._bitset)
        elif op == Py_GT:
            return bitset_issuperset(left._bitset, right._bitset) and not bitset_eq(left._bitset, right._bitset)
        elif op == Py_GE:
            return bitset_issuperset(left._bitset, right._bitset)

    cdef FrozenBitset _new(self, long int capacity):
        """
        Return an object of the same type as ``self``, initialized with a
        bitset of capacity ``capacity``.
        """
        cdef Bitset b
        b = Bitset.__new__(Bitset, None, capacity)
        return b

    cpdef update(self, FrozenBitset other):
        """
        Update the bitset to include items in ``other``.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a.update(Bitset('0101'))
            sage: a
            1101
            sage: a_set = set(a)
            sage: a.update(Bitset('00011' * 25))
            sage: a
            11011000110001100011000110001100011000110001100011000110001100011000110001100011000110001100011000110001100011000110001100011
            sage: a_set.update(Bitset('00011' * 25))
            sage: set(a) == a_set
            True

        TESTS:

        During update, ``other`` cannot be ``None``. ::

            sage: a = Bitset('1101')
            sage: a.update(None)
            Traceback (most recent call last):
            ...
            TypeError: other cannot be None
        """
        if other is None:
            raise TypeError("other cannot be None")
        cdef bitset_t temp
        if self._bitset.size == other._bitset.size:
            bitset_union(self._bitset, self._bitset, other._bitset)
        elif self._bitset.size < other._bitset.size:
            bitset_realloc(self._bitset, other._bitset.size)
            bitset_union(self._bitset, self._bitset, other._bitset)
        else:
            bitset_init(temp, other._bitset.size)
            bitset_copy(temp, other._bitset)
            bitset_realloc(temp, self._bitset.size)
            bitset_union(self._bitset, self._bitset, temp)
            bitset_free(temp)

    def __ior__(self, FrozenBitset other not None):
        """
        Update the bitset to include items in ``other``.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a |= Bitset('0101')
            sage: a
            1101
            sage: a_set = set(a)
            sage: a |= Bitset('00011' * 25)
            sage: a
            11011000110001100011000110001100011000110001100011000110001100011000110001100011000110001100011000110001100011000110001100011
            sage: a_set |= set(Bitset('00011' * 25))
            sage: set(a) == a_set
            True

        TESTS::

            sage: a = Bitset('110')
            sage: a |= None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
        """
        self.update(other)
        return self

    cpdef intersection_update(self, FrozenBitset other):
        """
        Update the bitset to the intersection of ``self`` and ``other``.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a.intersection_update(Bitset('0101'))
            sage: a
            0100
            sage: a_set = set(a)
            sage: a.intersection_update(Bitset('0110' * 25))
            sage: a
            0100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: a_set.intersection_update(Bitset('0110' * 25))
            sage: set(a) == a_set
            True

        TESTS::

            sage: Bitset('110').intersection_update(None)
            Traceback (most recent call last):
            ...
            TypeError: other cannot be None
        """
        if other is None:
            raise TypeError("other cannot be None")
        cdef bitset_t temp
        if self._bitset.size == other._bitset.size:
            bitset_intersection(self._bitset, self._bitset, other._bitset)
        elif self._bitset.size < other._bitset.size:
            bitset_realloc(self._bitset, other._bitset.size)
            bitset_intersection(self._bitset, self._bitset, other._bitset)
        else:
            bitset_init(temp, other._bitset.size)
            bitset_copy(temp, other._bitset)
            bitset_realloc(temp, self._bitset.size)
            bitset_intersection(self._bitset, self._bitset, temp)
            bitset_free(temp)

    def __iand__(self, FrozenBitset other not None):
        """
        Update the bitset to the intersection of ``self`` and ``other``.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a &= Bitset('0101')
            sage: a
            0100
            sage: a_set = set(a)
            sage: a &= Bitset('0110' * 25)
            sage: a
            0100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: a_set &= set(Bitset('0110' * 25))
            sage: a_set == set(a)
            True

        TESTS::

            sage: a = Bitset('110')
            sage: a &= None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
        """
        self.intersection_update(other)
        return self

    cpdef difference_update(self, FrozenBitset other):
        """
        Update the bitset to the difference of ``self`` and ``other``.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a.difference_update(Bitset('0101'))
            sage: a
            1000
            sage: a_set = set(a)
            sage: a.difference_update(FrozenBitset('010101' * 10)); a
            100000000000000000000000000000000000000000000000000000000000
            sage: a_set.difference_update(FrozenBitset('010101' * 10))
            sage: a_set == set(a)
            True
            sage: a.difference_update(FrozenBitset('110'))
            sage: a_set.difference_update(FrozenBitset('110'))
            sage: a_set == set(a)
            True
            sage: a.difference_update(FrozenBitset('01010' * 20)); a
            0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: a_set.difference_update(FrozenBitset('01010' * 20))
            sage: a_set == set(a)
            True
            sage: b = Bitset('10101' * 20)
            sage: b_set = set(b)
            sage: b.difference_update(FrozenBitset('1' * 5)); b
            0000010101101011010110101101011010110101101011010110101101011010110101101011010110101101011010110101
            sage: b_set.difference_update(FrozenBitset('1' * 5))
            sage: b_set == set(b)
            True

        TESTS::

            sage: Bitset('110').difference_update(None)
            Traceback (most recent call last):
            ...
            TypeError: other cannot be None
        """
        if other is None:
            raise TypeError("other cannot be None")
        cdef bitset_t temp
        if self._bitset.size == other._bitset.size:
            bitset_difference(self._bitset, self._bitset, other._bitset)
        elif self._bitset.size < other._bitset.size:
            bitset_realloc(self._bitset, other._bitset.size)
            bitset_difference(self._bitset, self._bitset, other._bitset)
        else:
            bitset_init(temp, other._bitset.size)
            bitset_copy(temp, other._bitset)
            bitset_realloc(temp, self._bitset.size)
            bitset_difference(self._bitset, self._bitset, temp)
            bitset_free(temp)

    def __isub__(self, FrozenBitset other not None):
        """
        Update the bitset to the difference of ``self`` and ``other``.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a -= Bitset('0101')
            sage: a
            1000
            sage: a_set = set(a)
            sage: a -= FrozenBitset('010101' * 10); a
            100000000000000000000000000000000000000000000000000000000000
            sage: a_set -= set(FrozenBitset('010101' * 10))
            sage: a_set == set(a)
            True
            sage: a = Bitset('110')
            sage: a_set = set(a)
            sage: a -= FrozenBitset('01010' * 20); a
            1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: a_set -= set(FrozenBitset('01010' * 20))
            sage: a_set == set(a)
            True
            sage: b = FrozenBitset('10101' * 20)
            sage: b_set = set(b)
            sage: b -= FrozenBitset('1' * 5); b
            0000010101101011010110101101011010110101101011010110101101011010110101101011010110101101011010110101
            sage: b_set -= FrozenBitset('1' * 5)
            sage: b_set == set(b)
            True

        TESTS::

            sage: a = Bitset('110')
            sage: a -= None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
        """
        self.difference_update(other)
        return self

    cpdef symmetric_difference_update(self, FrozenBitset other):
        """
        Update the bitset to the symmetric difference of ``self`` and
        ``other``.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a.symmetric_difference_update(Bitset('0101'))
            sage: a
            1001
            sage: a_set = set(a)
            sage: a.symmetric_difference_update(FrozenBitset('010101' * 10)); a
            110001010101010101010101010101010101010101010101010101010101
            sage: a_set.symmetric_difference_update(FrozenBitset('010101' * 10))
            sage: a_set == set(a)
            True
            sage: a.symmetric_difference_update(FrozenBitset('01010' * 20)); a
            1001011111000001111100000111110000011111000001111100000111110101001010010100101001010010100101001010
            sage: a_set.symmetric_difference_update(FrozenBitset('01010' * 20))
            sage: a_set == set(a)
            True
            sage: b = Bitset('10101' * 20)
            sage: b_set = set(b)
            sage: b.symmetric_difference_update( FrozenBitset('1' * 5)); b
            0101010101101011010110101101011010110101101011010110101101011010110101101011010110101101011010110101
            sage: b_set.symmetric_difference_update( FrozenBitset('1' * 5))
            sage: b_set == set(b)
            True

        TESTS::

            sage: Bitset('110').symmetric_difference_update(None)
            Traceback (most recent call last):
            ...
            TypeError: other cannot be None
        """
        if other is None:
            raise TypeError("other cannot be None")
        cdef bitset_t temp
        if self._bitset.size == other._bitset.size:
            bitset_symmetric_difference(self._bitset, self._bitset, other._bitset)
        elif self._bitset.size < other._bitset.size:
            bitset_realloc(self._bitset, other._bitset.size)
            bitset_symmetric_difference(self._bitset, self._bitset, other._bitset)
        else:
            bitset_init(temp, other._bitset.size)
            bitset_copy(temp, other._bitset)
            bitset_realloc(temp, self._bitset.size)
            bitset_symmetric_difference(self._bitset, self._bitset, temp)
            bitset_free(temp)

    def __ixor__(self, FrozenBitset other not None):
        """
        Update the bitset to the symmetric difference of ``self`` and
        ``other``.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a ^^=Bitset('0101')
            sage: a
            1001
            sage: a_set = set(a)
            sage: a ^^= FrozenBitset('010101' * 10); a
            110001010101010101010101010101010101010101010101010101010101
            sage: a_set ^^= set(FrozenBitset('010101' * 10))
            sage: a_set == set(a)
            True
            sage: a ^^= FrozenBitset('01010' * 20); a
            1001011111000001111100000111110000011111000001111100000111110101001010010100101001010010100101001010
            sage: a_set ^^= set(FrozenBitset('01010' * 20))
            sage: a_set == set(a)
            True
            sage: b = Bitset('10101' * 20)
            sage: b_set = set(b)
            sage: b ^^= FrozenBitset('1' * 5); b
            0101010101101011010110101101011010110101101011010110101101011010110101101011010110101101011010110101
            sage: b_set ^^= set(FrozenBitset('1' * 5))
            sage: b_set == set(b)
            True

        TESTS::

            sage: a = Bitset('110')
            sage: a ^^= None
            Traceback (most recent call last):
            ...
            TypeError: Argument 'other' has incorrect type (expected sage.data_structures.bitset.FrozenBitset, got NoneType)
        """
        self.symmetric_difference_update(other)
        return self

    cpdef add(self, unsigned long n):
        """
        Update the bitset by adding ``n``.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a.add(5)
            sage: a
            110001
            sage: a.add(100)
            sage: sorted(list(a))
            [0, 1, 5, 100]
            sage: a.capacity()
            101

        TESTS:

        The input ``n`` must be an integer. ::

            sage: Bitset('110').add(None)
            Traceback (most recent call last):
            ...
            TypeError: an integer is required
        """
        if n >= self._bitset.size:
            bitset_realloc(self._bitset, n + 1)
        bitset_add(self._bitset, n)

    cpdef remove(self, unsigned long n):
        """
        Update the bitset by removing ``n``.

        This raises a :exc:`KeyError` if ``n`` is not contained
        in the bitset.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a.remove(1)
            sage: a
            100
            sage: a.remove(2)
            Traceback (most recent call last):
            ...
            KeyError: 2
            sage: a.remove(4)
            Traceback (most recent call last):
            ...
            KeyError: 4
            sage: a
            100
            sage: a = Bitset('000001' * 15); sorted(list(a))
            [5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71, 77, 83, 89]
            sage: a.remove(83); sorted(list(a))
            [5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71, 77, 89]

        TESTS:

        The input ``n`` must be an integer. ::

            sage: Bitset('110').remove(None)
            Traceback (most recent call last):
            ...
            TypeError: an integer is required
        """
        if n >= self._bitset.size:
            raise KeyError(n)
        else:
            bitset_remove(self._bitset, n)

    cpdef discard(self, unsigned long n):
        """
        Update the bitset by removing ``n``.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a.discard(1)
            sage: a
            100
            sage: a.discard(2)
            sage: a.discard(4)
            sage: a
            100
            sage: a = Bitset('000001' * 15); sorted(list(a))
            [5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71, 77, 83, 89]
            sage: a.discard(83); sorted(list(a))
            [5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71, 77, 89]
            sage: a.discard(82); sorted(list(a))
            [5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71, 77, 89]

        TESTS:

        The input ``n`` must be an integer. ::

            sage: Bitset('110').discard(None)
            Traceback (most recent call last):
            ...
            TypeError: an integer is required
        """
        if n < self._bitset.size:
            bitset_discard(self._bitset, n)

    cpdef pop(self):
        """
        Remove and return an arbitrary element from the set.

        This raises a :exc:`KeyError` if the set is empty.

        EXAMPLES::

            sage: a = Bitset('011')
            sage: a.pop()
            1
            sage: a
            001
            sage: a.pop()
            2
            sage: a
            000
            sage: a.pop()
            Traceback (most recent call last):
            ...
            KeyError: 'pop from an empty set'
            sage: a = Bitset('0001'*32)
            sage: a.pop()
            3
            sage: [a.pop() for _ in range(20)]
            [7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83]
        """
        return bitset_pop(self._bitset)

    cpdef clear(self):
        """
        Remove all elements from the bitset.

        EXAMPLES::

            sage: a = Bitset('011')
            sage: a.clear()
            sage: a
            000
            sage: a = Bitset('011' * 32)
            sage: a.clear()
            sage: set(a)
            set()
        """
        bitset_clear(self._bitset)


#############################################################################
# Bitset Testing: test basic Cython bitsets
#############################################################################

def test_bitset(py_a, py_b, long n):
    """
    Test the Cython bitset functions so we can have some relevant doctests.

    TESTS::

        sage: from sage.data_structures.bitset import test_bitset
        sage: test_bitset('00101', '01110', 4)
        a 00101
        list a [2, 4]
        a.size 5
        len(a) 2
        a.limbs 1
        b 01110
        a.in(n)   True
        a.not_in(n)   False
        a.add(n)     00101
        a.discard(n)   00100
        a.set_to(n)  00101
        a.flip(n)    00100
        a.set_first_n(n)    11110
        a.first_in_complement()    4
        a.isempty()  False
        a.eq(b)      False
        a.cmp(b)     1
        a.lex_cmp(b) -1
        a.issubset(b) False
        a.issuperset(b) False
        a.copy()     00101
        r.clear()     00000
        complement a        11010
        a intersect b      00100
        a union b       01111
        a minus b      00001
        a symmetric_difference b      01011
        a.rshift(n)  10000
        a.lshift(n)  00000
        a.first()           2
        a.next(n)           4
        a.first_diff(b)     1
        a.next_diff(b, n)   4
        a.hamming_weight()  2
        a.map(m)  10100
        a == loads(dumps(a))  True
        reallocating a      00101
        to size 4          0010
        to size 8          00100000
        to original size    00100

    ::

        sage: test_bitset('11101', '11001', 2)
        a 11101
        list a [0, 1, 2, 4]
        a.size 5
        len(a) 4
        a.limbs 1
        b 11001
        a.in(n)   True
        a.not_in(n)   False
        a.add(n)     11101
        a.discard(n)   11001
        a.set_to(n)  11101
        a.flip(n)    11001
        a.set_first_n(n)    11000
        a.first_in_complement()    2
        a.isempty()  False
        a.eq(b)      False
        a.cmp(b)     1
        a.lex_cmp(b) 1
        a.issubset(b) False
        a.issuperset(b) True
        a.copy()     11101
        r.clear()     00000
        complement a        00010
        a intersect b      11001
        a union b       11101
        a minus b      00100
        a symmetric_difference b      00100
        a.rshift(n)  10100
        a.lshift(n)  00111
        a.first()           0
        a.next(n)           2
        a.first_diff(b)     2
        a.next_diff(b, n)   2
        a.hamming_weight()  4
        a.map(m)  10111
        a == loads(dumps(a))  True
        reallocating a      11101
        to size 2          11
        to size 4          1100
        to original size    11000

    Test a corner-case: a bitset that is a multiple of words::

        sage: test_bitset('00'*64, '01'*64, 127)
        a 00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        list a []
        a.size 128
        len(a) 0
        a.limbs ...
        b 01010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101
        a.in(n)   False
        a.not_in(n)   True
        a.add(n)     00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
        a.discard(n)   00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        a.set_to(n)  00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
        a.flip(n)    00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
        a.set_first_n(n)    11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111110
        a.first_in_complement()    127
        a.isempty()  True
        a.eq(b)      False
        a.cmp(b)     -1
        a.lex_cmp(b) -1
        a.issubset(b) True
        a.issuperset(b) False
        a.copy()     00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        r.clear()     00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        complement a        11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
        a intersect b      00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        a union b       01010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101
        a minus b      00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        a symmetric_difference b      01010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101
        a.rshift(n)  00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        a.lshift(n)  00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        a.first()           -1
        a.next(n)           -1
        a.first_diff(b)     1
        a.next_diff(b, n)   127
        a.hamming_weight()  0
        a.map(m)  00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        a == loads(dumps(a))  True
        rshifts add  True
        lshifts add  True
        intersection commutes True
        union commutes  True
        not not = id True
        flipped bit  127
        add bit      127
        discard bit    127
        lshift add unset ok True
        rshift set unset ok True
        reallocating a      00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        to size 127          0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        to size 254          00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        to original size    00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

    Large enough to span multiple limbs.  We don't explicitly check the number of limbs below because it will be different in the 32 bit versus 64 bit cases::

        sage: test_bitset('111001'*25, RealField(151)(pi).str(2)[2:], 69)               # needs sage.symbolic
        a 111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001
        list a [0, 1, 2, 5, 6, 7, 8, 11, 12, 13, 14, 17, 18, 19, 20, 23, 24, 25, 26, 29, 30, 31, 32, 35, 36, 37, 38, 41, 42, 43, 44, 47, 48, 49, 50, 53, 54, 55, 56, 59, 60, 61, 62, 65, 66, 67, 68, 71, 72, 73, 74, 77, 78, 79, 80, 83, 84, 85, 86, 89, 90, 91, 92, 95, 96, 97, 98, 101, 102, 103, 104, 107, 108, 109, 110, 113, 114, 115, 116, 119, 120, 121, 122, 125, 126, 127, 128, 131, 132, 133, 134, 137, 138, 139, 140, 143, 144, 145, 146, 149]
        a.size 150
        len(a) 100
        a.limbs ...
        b 000100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000000011011100000111001101000100101001000000100100111
        a.in(n)   False
        a.not_in(n)   True
        a.add(n)     111001111001111001111001111001111001111001111001111001111001111001111101111001111001111001111001111001111001111001111001111001111001111001111001111001
        a.discard(n)   111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001
        a.set_to(n)  111001111001111001111001111001111001111001111001111001111001111001111101111001111001111001111001111001111001111001111001111001111001111001111001111001
        a.flip(n)    111001111001111001111001111001111001111001111001111001111001111001111101111001111001111001111001111001111001111001111001111001111001111001111001111001
        a.set_first_n(n)    111111111111111111111111111111111111111111111111111111111111111111111000000000000000000000000000000000000000000000000000000000000000000000000000000000
        a.first_in_complement()    69
        a.isempty()  False
        a.eq(b)      False
        a.cmp(b)     -1
        a.lex_cmp(b) 1
        a.issubset(b) False
        a.issuperset(b) False
        a.copy()     111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001
        r.clear()     000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        complement a        000110000110000110000110000110000110000110000110000110000110000110000110000110000110000110000110000110000110000110000110000110000110000110000110000110
        a intersect b      000000100001111000110001010001000000001001010001100001000000100000001001100001001000010000010001000000011001100000111001101000100001001000000000100001
        a union b       111101111001111111111101111001111101111011111001111001111111111111111001111011111101111101111111111001111011111001111001111001111101111001111101111111
        a minus b      111001011000000001001000101000111001110000101000011000111001011001110000011000110001101001101000111001100000011001000000010001011000110001111001011000
        a symmetric_difference b      111101011000000111001100101000111101110010101000011000111111011111110000011010110101101101101110111001100010011001000000010001011100110001111101011110
        a.rshift(n)  001111001111001111001111001111001111001111001111001111001111001111001111001111001000000000000000000000000000000000000000000000000000000000000000000000
        a.lshift(n)  000000000000000000000000000000000000000000000000000000000000000000000111001111001111001111001111001111001111001111001111001111001111001111001111001111
        a.first()           0
        a.next(n)           71
        a.first_diff(b)     0
        a.next_diff(b, n)   73
        a.hamming_weight()  100
        a.map(m)  100111100111100111100111100111100111100111100111100111100111100111100111100111100111100111100111100111100111100111100111100111100111100111100111100111
        a == loads(dumps(a))  True
        rshifts add  True
        lshifts add  True
        intersection commutes True
        union commutes  True
        not not = id True
        flipped bit  69
        add bit      69
        discard bit    69
        lshift add unset ok True
        rshift set unset ok True
        reallocating a      111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001111001
        to size 69          111001111001111001111001111001111001111001111001111001111001111001111
        to size 138          111001111001111001111001111001111001111001111001111001111001111001111000000000000000000000000000000000000000000000000000000000000000000000
        to original size    111001111001111001111001111001111001111001111001111001111001111001111000000000000000000000000000000000000000000000000000000000000000000000000000000000
    """
    cdef bint bit = True
    cdef bitset_t a, b, r

    bitset_from_str(a, py_a)
    bitset_from_str(b, py_b)

    if a.size != b.size:
        raise ValueError("inputs must have same size")

    print("a", bitset_string(a))
    print("list a", bitset_list(a))
    print("a.size", a.size)
    print("len(a)", bitset_len(a))
    print("a.limbs", a.limbs)
    print("b", bitset_string(b))
    print("a.in(n)  ", bitset_in(a, n))
    print("a.not_in(n)  ", bitset_not_in(a, n))
    bitset_add(a, n)
    print("a.add(n)    ", bitset_string(a))
    bitset_from_str(a, py_a)
    bitset_discard(a, n)
    print("a.discard(n)  ", bitset_string(a))
    bitset_from_str(a, py_a)
    bitset_set_to(a, n, bit)
    print("a.set_to(n) ", bitset_string(a))
    bitset_from_str(a, py_a)
    bitset_flip(a, n)
    print("a.flip(n)   ", bitset_string(a))
    bitset_set_first_n(a, n)
    print("a.set_first_n(n)   ", bitset_string(a))
    print("a.first_in_complement()   ", bitset_first_in_complement(a))

    bitset_from_str(a, py_a)
    bitset_from_str(b, py_b)
    print("a.isempty() ", bitset_isempty(a))
    print("a.eq(b)     ", bitset_eq(a, b))
    print("a.cmp(b)    ", bitset_cmp(a, b))
    print("a.lex_cmp(b)", bitset_lex_cmp(a, b))
    print("a.issubset(b)", bitset_issubset(a, b))
    print("a.issuperset(b)", bitset_issuperset(a, b))

    bitset_from_str(a, py_a)
    bitset_from_str(b, py_b)

    bitset_init(r, a.size)
    bitset_copy(r, a)
    print("a.copy()    ", bitset_string(r))
    bitset_clear(r)
    print("r.clear()    ", bitset_string(r))
    bitset_complement(r, a)
    print("complement a       ", bitset_string(r))
    bitset_intersection(r, a, b)
    print("a intersect b     ", bitset_string(r))
    bitset_union(r, a, b)
    print("a union b      ", bitset_string(r))
    bitset_difference(r, a, b)
    print("a minus b     ", bitset_string(r))
    bitset_symmetric_difference(r, a, b)
    print("a symmetric_difference b     ", bitset_string(r))

    bitset_rshift(r, a, n)
    print("a.rshift(n) ", bitset_string(r))

    bitset_lshift(r, a, n)
    print("a.lshift(n) ", bitset_string(r))

    print("a.first()          ", bitset_first(a))
    print("a.next(n)          ", bitset_next(a, n))
    print("a.first_diff(b)    ", bitset_first_diff(a, b))
    print("a.next_diff(b, n)  ", bitset_next_diff(a, b, n))

    print("a.hamming_weight() ", bitset_hamming_weight(a))

    morphism = {i: a.size - i - 1 for i in range(a.size)}
    bitset_map(r, a, morphism)
    print("a.map(m) ", bitset_string(r))

    data = bitset_pickle(a)
    bitset_unpickle(r, data)
    print("a == loads(dumps(a)) ", bitset_eq(r, a))

    cdef bitset_t s
    bitset_init(s, a.size)

    if a.size > 100:
        bitset_rshift(r, b, 3)
        bitset_rshift(r, r, 77)
        bitset_rshift(s, b, 80)
        print("rshifts add ", bitset_eq(s, r))

        bitset_lshift(r, b, 69)
        bitset_lshift(r, r, 6)
        bitset_lshift(s, b, 75)
        print("lshifts add ", bitset_eq(s, r))

        bitset_intersection(r, a, b)
        bitset_intersection(s, b, a)
        print("intersection commutes", bitset_eq(s, r))

        bitset_union(r, a, b)
        bitset_union(s, b, a)
        print("union commutes ", bitset_eq(s, r))

        bitset_complement(r, b)
        bitset_complement(s, r)
        print("not not = id", bitset_eq(s, b))

        bitset_copy(r, b)
        bitset_flip(r, n)
        print("flipped bit ", bitset_first_diff(b, r))

        bitset_clear(r)
        bitset_add(r, n)
        print("add bit     ", bitset_first(r))

        bitset_clear(r)
        bitset_complement(r, r)
        bitset_discard(r, n)
        bitset_complement(r, r)
        print("discard bit   ", bitset_first(r))

        bitset_clear(r)
        bitset_add(r, 10)
        bitset_lshift(r, r, 68)
        bitset_flip(r, 78)
        print("lshift add unset ok", bitset_isempty(r))

        bitset_clear(r)
        bitset_add(r, 19)
        bitset_rshift(r, r, 8)
        bitset_discard(r, 11)
        print("rshift set unset ok", bitset_isempty(r))

    print("reallocating a     ", bitset_string(a))
    bitset_realloc(a, n)
    print("to size %d         " % n, bitset_string(a))
    bitset_realloc(a, 2 * n)
    print("to size %d         " % (2 * n), bitset_string(a))
    bitset_realloc(a, b.size)
    print("to original size   ", bitset_string(a))

    bitset_free(a)
    bitset_free(b)
    bitset_free(r)
    bitset_free(s)


def test_bitset_set_first_n(py_a, long n):
    """
    Test the bitset function set_first_n.

    TESTS::

        sage: from sage.data_structures.bitset import test_bitset_set_first_n
        sage: test_bitset_set_first_n('00'*64, 128)
        a.set_first_n(n)    11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
    """
    cdef bint bit = True
    cdef bitset_t a

    bitset_from_str(a, py_a)
    bitset_set_first_n(a, n)
    print("a.set_first_n(n)   ", bitset_string(a))
    bitset_free(a)


def test_bitset_remove(py_a, long n):
    """
    Test the bitset_remove function.

    TESTS::

        sage: from sage.data_structures.bitset import test_bitset_remove
        sage: test_bitset_remove('01', 0)
        Traceback (most recent call last):
        ...
        KeyError: 0
        sage: test_bitset_remove('01', 1)
        a 01
        a.size 2
        a.limbs 1
        n 1
        a.remove(n)   00
    """
    cdef bitset_t a
    bitset_from_str(a, py_a)

    print("a", bitset_string(a))
    print("a.size", a.size)
    print("a.limbs", a.limbs)
    print("n", n)

    bitset_remove(a, n)
    print("a.remove(n)  ", bitset_string(a))

    bitset_free(a)


def test_bitset_pop(py_a):
    """
    Test for the bitset_pop function.

    TESTS::

        sage: from sage.data_structures.bitset import test_bitset_pop
        sage: test_bitset_pop('0101')
        a.pop()   1
        new set:  0001
        sage: test_bitset_pop('0000')
        Traceback (most recent call last):
        ...
        KeyError: 'pop from an empty set'
    """
    cdef bitset_t a
    bitset_from_str(a, py_a)
    i = bitset_pop(a)
    print("a.pop()  ", i)
    print("new set: ", bitset_string(a))
    bitset_free(a)


def test_bitset_unpickle(data):
    """
    This (artificially) tests pickling of bitsets across systems.

    INPUT:

    - ``data`` -- tuple of data as would be produced by the internal,
      Cython-only, method ``bitset_pickle``

    OUTPUT: list form of the bitset corresponding to the pickled data

    EXAMPLES:

    We compare 64-bit and 32-bit encoding. Both should unpickle on any system::

        sage: from sage.data_structures.bitset import test_bitset_unpickle
        sage: test_bitset_unpickle((0, 100, 2, 8, (33, 6001)))
        [0, 5, 64, 68, 69, 70, 72, 73, 74, 76]
        sage: test_bitset_unpickle((0, 100, 4, 4, (33, 0, 6001, 0)))
        [0, 5, 64, 68, 69, 70, 72, 73, 74, 76]
    """
    cdef bitset_t bs
    bitset_init(bs, 1)
    bitset_unpickle(bs, data)
    L = bitset_list(bs)
    bitset_free(bs)
    return L


def test_bitset_copy_flex(py_a):
    """
    TESTS:

    Check that :issue:`33012` is fixed::

        sage: from sage.data_structures.bitset import test_bitset_copy_flex
        sage: test_bitset_copy_flex('0101'*100)
    """
    cdef bitset_t a, b

    bitset_from_str(a, py_a)
    bitset_init(b, a.size*2)

    bitset_copy_flex(b, a)
    if not bitset_list(b) == bitset_list(a):
        raise ValueError

    bitset_free(a)
    bitset_free(b)
