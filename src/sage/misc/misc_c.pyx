# sage_setup: distribution = sagemath-objects
"""
Miscellaneous functions (Cython)

This file contains support for products, running totals, balanced
sums, and also a function to flush output from external library calls.

AUTHORS:

- William Stein (2005)
- Joel B. Mohler (2007-10-03): Reimplemented in Cython and optimized
- Robert Bradshaw (2007-10-26): Balanced product tree, other optimizations, (lazy) generator support
- Robert Bradshaw (2008-03-26): Balanced product tree for generators and iterators
- Stefan van Zwam (2013-06-06): Added bitset tests, some docstring cleanup
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from copy import copy

from cpython.sequence cimport *
from cpython.list cimport *
from cpython.tuple cimport *
from cpython.number cimport *
from libc.stdio cimport fflush

cdef extern from *:
    bint PyGen_Check(x)


def running_total(L, start=None):
    """
    Return a list where the `i`-th entry is the sum of all entries up to (and
    including) `i`.

    INPUT:

    - ``L`` -- the list
    - ``start`` -- (optional) a default start value

    EXAMPLES::

        sage: running_total(range(5))
        [0, 1, 3, 6, 10]
        sage: running_total("abcdef")
        ['a', 'ab', 'abc', 'abcd', 'abcde', 'abcdef']
        sage: running_total("abcdef", start="%")
        ['%a', '%ab', '%abc', '%abcd', '%abcde', '%abcdef']
        sage: running_total([1..10], start=100)
        [101, 103, 106, 110, 115, 121, 128, 136, 145, 155]
        sage: running_total([])
        []
    """
    running = []
    total = start
    for x in L:
        if total is None:
            # This is the first entry
            total = x
        else:
            total += x
        PyList_Append(running, total)
    return running


def prod(x, z=None, Py_ssize_t recursion_cutoff=5):
    """
    Return the product of the elements in the list x.

    If optional argument z is not given, start the product with the first
    element of the list, otherwise use z.  The empty product is the int 1 if z
    is not specified, and is z if given.

    This assumes that your multiplication is associative; we do not promise
    which end of the list we start at.

    .. SEEALSO::

        For the symbolic product function, see :func:`sage.calculus.calculus.symbolic_product`.

    EXAMPLES::

        sage: prod([1,2,34])
        68
        sage: prod([2,3], 5)
        30
        sage: prod((1,2,3), 5)
        30
        sage: F = factor(-2006); F
        -1 * 2 * 17 * 59
        sage: prod(F)
        -2006

    AUTHORS:

    - Joel B. Mohler (2007-10-03): Reimplemented in Cython and optimized
    - Robert Bradshaw (2007-10-26): Balanced product tree, other optimizations, (lazy) generator support
    - Robert Bradshaw (2008-03-26): Balanced product tree for generators and iterators
    """
    cdef Py_ssize_t n

    if type(x) is not list and type(x) is not tuple:

        if not PyGen_Check(x):

            try:
                return x.prod()
            except AttributeError:
                pass

            try:
                return x.mul()
            except AttributeError:
                pass

            try:
                n = len(x)
                if n < 1000:   # arbitrary limit
                    x = list(x)
            except TypeError:
                pass

        if type(x) is not list:
            try:
                return iterator_prod(x, z)
            except StopIteration:
                x = []

    n = len(x)

    if n == 0:
        if z is None:
            import sage.rings.integer
            return sage.rings.integer.Integer(1)
        else:
            return z

    prod = balanced_list_prod(x, 0, n, recursion_cutoff)

    if z is not None:
        prod = z * prod

    return prod


cdef balanced_list_prod(L, Py_ssize_t offset, Py_ssize_t count, Py_ssize_t cutoff):
    """
    INPUT:

    - ``L`` -- the terms (MUST be a tuple or list)
    - ``off`` -- offset in the list from which to start
    - ``count`` -- how many terms in the product
    - ``cutoff`` -- the minimum count to recurse on

    OUTPUT:

    ``L[offset] * L[offset+1] * ... * L[offset+count-1]``

    .. NOTE::

        The parameter cutoff must be at least 1, and there is no reason to
        ever make it less than 3. However, there are at least two advantages
        to setting it higher (and consequently not recursing all the way
        down the tree). First, one avoids the overhead of the function
        calls at the base of the tree (which is the majority of them) and
        second, it allows one to save on object creation if inplace
        operations are used. The asymptotic gains should usually be at the
        top of the tree anyway.
    """
    cdef Py_ssize_t k
    if count <= cutoff:
        prod = <object>PySequence_Fast_GET_ITEM(L, offset)
        for k from offset < k < offset + count:
            prod *= <object>PySequence_Fast_GET_ITEM(L, k)
        return prod
    else:
        k = (1 + count) >> 1
        return balanced_list_prod(L, offset, k, cutoff) * balanced_list_prod(L, offset + k, count - k, cutoff)


cpdef iterator_prod(L, z=None, bint multiply=True):
    """
    Attempt to do a balanced product of an arbitrary and unknown length
    sequence (such as a generator). Intermediate multiplications are always
    done with subproducts of the same size (measured by the number of original
    factors) up until the iterator terminates. This is optimal when and only
    when there are exactly a power of two number of terms.

    A StopIteration is raised if the iterator is empty and z is not given.

    EXAMPLES::

        sage: from sage.misc.misc_c import iterator_prod
        sage: iterator_prod(1..5)
        120
        sage: iterator_prod([], z='anything')
        'anything'

        sage: from sage.misc.misc_c import NonAssociative
        sage: L = [NonAssociative(label) for label in 'abcdef']
        sage: iterator_prod(L)
        (((a*b)*(c*d))*(e*f))

    When ``multiply=False``, the items are added up instead (however this
    interface should not be used directly, use :func:`balanced_sum` instead)::

        sage: iterator_prod((1..5), multiply=False)
        15
    """
    cdef list sub_prods
    L = iter(L)
    if z is None:
        sub_prods = [next(L)] * 10  # only take one element from L, the rest are just placeholders
        # the list size can be dynamically increased later
    else:
        sub_prods = [z] * 10

    cdef Py_ssize_t j
    cdef Py_ssize_t i = 1
    cdef Py_ssize_t tip = 0

    for x in L:
        i += 1
        if i & 1:
            # for odd i we extend the stack
            tip += 1
            if len(sub_prods) == tip:
                sub_prods.append(x)
            else:
                sub_prods[tip] = x
            continue
        else:
            # for even i we multiply the stack down
            # by the number of factors of 2 in i
            if multiply:
                x = sub_prods[tip] * x
            else:
                x = sub_prods[tip] + x
            for j from 1 <= j < 64:
                if i & (1 << j):
                    break
                tip -= 1
                if multiply:
                    x = sub_prods[tip] * x
                else:
                    x = sub_prods[tip] + x
            sub_prods[tip] = x

    while tip > 0:
        tip -= 1
        if multiply:
            sub_prods[tip] *= sub_prods[tip + 1]
        else:
            sub_prods[tip] += sub_prods[tip + 1]

    return sub_prods[0]


class NonAssociative:
    """
    This class is to test the balance nature of prod.

    EXAMPLES::

        sage: from sage.misc.misc_c import NonAssociative
        sage: L = [NonAssociative(label) for label in 'abcdef']
        sage: prod(L)
        (((a*b)*c)*((d*e)*f))
        sage: L = [NonAssociative(label) for label in range(20)]
        sage: prod(L, recursion_cutoff=5)
        ((((((0*1)*2)*3)*4)*((((5*6)*7)*8)*9))*(((((10*11)*12)*13)*14)*((((15*16)*17)*18)*19)))
        sage: prod(L, recursion_cutoff=1)
        (((((0*1)*2)*(3*4))*(((5*6)*7)*(8*9)))*((((10*11)*12)*(13*14))*(((15*16)*17)*(18*19))))
        sage: L = [NonAssociative(label) for label in range(14)]
        sage: prod(L, recursion_cutoff=1)
        ((((0*1)*(2*3))*((4*5)*6))*(((7*8)*(9*10))*((11*12)*13)))
    """
    def __init__(self, left, right=None):
        """
        EXAMPLES::

            sage: from sage.misc.misc_c import NonAssociative
            sage: NonAssociative('a')
            a
            sage: NonAssociative('a','b')
            (a*b)
        """
        self.left = left
        self.right = right

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.misc.misc_c import NonAssociative
            sage: NonAssociative(1)
            1
            sage: NonAssociative(2,3)
            (2*3)
        """
        if self.right is None:
            return str(self.left)
        else:
            return "(%s*%s)" % (self.left, self.right)

    def __mul__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.misc_c import NonAssociative
            sage: a, b, c = [NonAssociative(label) for label in 'abc']
            sage: (a*b)*c
            ((a*b)*c)
            sage: a*(b*c)
            (a*(b*c))
        """
        return NonAssociative(self, other)


def balanced_sum(x, z=None, Py_ssize_t recursion_cutoff=5):
    """
    Return the sum of the elements in the list x.  If optional
    argument z is not given, start the sum with the first element of
    the list, otherwise use z.  The empty product is the int 0 if z is
    not specified, and is z if given.  The sum is computed
    recursively, where the sum is split up if the list is greater than
    recursion_cutoff.  recursion_cutoff must be at least 3.

    This assumes that your addition is associative; we do not promise
    which end of the list we start at.

    EXAMPLES::

        sage: balanced_sum([1,2,34])
        37
        sage: balanced_sum([2,3], 5)
        10
        sage: balanced_sum((1,2,3), 5)
        11

    Order should be preserved::

        sage: balanced_sum([[i] for i in range(10)], [], recursion_cutoff=3)
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    We make copies when appropriate so that we do not accidentally
    modify the arguments::

        sage: list(range(10^5))==balanced_sum([[i] for i in range(10^5)], [])
        True
        sage: list(range(10^5))==balanced_sum([[i] for i in range(10^5)], [])
        True

    TESTS::

        sage: balanced_sum((1..3)) # nonempty, z=None
        6
        sage: balanced_sum((1..-1)) # empty, z=None
        0
        sage: balanced_sum((1..3), 5) # nonempty, z is not None
        11
        sage: balanced_sum((1..-1), 5) # empty, z is not None
        5
        sage: balanced_sum([1])
        1

    AUTHORS:

    - Joel B. Mohler (2007-10-03): Reimplemented in Cython and optimized
    - Robert Bradshaw (2007-10-26): Balanced product tree, other optimizations, (lazy) generator support
    """
    if recursion_cutoff < 3:
        raise ValueError("recursion_cutoff must be at least 3")

    if type(x) is not list and type(x) is not tuple:

        if PyGen_Check(x):
            return iterator_prod(x, z, multiply=False)
        else:
            try:
                return x.sum()
            except AttributeError:
                pass

            x = list(x)

    cdef Py_ssize_t n = len(x)

    if n == 0:
        if z is None:
            import sage.rings.integer
            return sage.rings.integer.Integer(0)
        else:
            return z

    sum = balanced_list_sum(x, 0, n, recursion_cutoff)

    if z is not None:
        sum = z + sum

    return sum


cdef balanced_list_sum(L, Py_ssize_t offset, Py_ssize_t count, Py_ssize_t cutoff):
    """
    INPUT:

    - ``L`` -- the terms (MUST be a tuple or list)
    - ``off`` -- offset in the list from which to start
    - ``count`` -- how many terms in the sum; must be positive
    - ``cutoff`` -- the minimum count to recurse on; must be at least 2

    OUTPUT:

    ``L[offset] + L[offset+1] + ... + L[offset+count-1]``

    .. NOTE::

        The parameter cutoff must be at least 3. However, there are
        at least two advantages to setting it higher (and
        consequently not recursing all the way down the
        tree). First, one avoids the overhead of the function calls
        at the base of the tree (which is the majority of them) and
        second, it allows one to save on object creation if inplace
        operations are used. The asymptotic gains should usually be
        at the top of the tree anyway.
    """
    cdef Py_ssize_t k
    if count <= cutoff:
        sum = <object>PySequence_Fast_GET_ITEM(L, offset)
        for k in range(offset + 1, offset + count):
            sum += <object>PySequence_Fast_GET_ITEM(L, k)
        return sum
    else:
        k = (1 + count) >> 1
        return balanced_list_sum(L, offset, k, cutoff) + balanced_list_sum(L, offset + k, count - k, cutoff)


cpdef list normalize_index(object key, int size):
    """
    Normalize an index key and return a valid index or list of indices
    within the range(0, size).

    INPUT:

    - ``key`` -- the index key, which can be either an integer, a tuple/list of
      integers, or a slice
    - ``size`` -- the size of the collection

    OUTPUT:

    A tuple (SINGLE, VALUE), where SINGLE is True (i.e., 1) if VALUE
    is an integer and False (i.e., 0) if VALUE is a list.

    EXAMPLES::

        sage: from sage.misc.misc_c import normalize_index
        sage: normalize_index(-6,5)
        Traceback (most recent call last):
        ...
        IndexError: index out of range
        sage: normalize_index(-5,5)
        [0]
        sage: normalize_index(-4,5)
        [1]
        sage: normalize_index(-3,5)
        [2]
        sage: normalize_index(-2,5)
        [3]
        sage: normalize_index(-1,5)
        [4]
        sage: normalize_index(0,5)
        [0]
        sage: normalize_index(1,5)
        [1]
        sage: normalize_index(2,5)
        [2]
        sage: normalize_index(3,5)
        [3]
        sage: normalize_index(4,5)
        [4]
        sage: normalize_index(5,5)
        Traceback (most recent call last):
        ...
        IndexError: index out of range
        sage: normalize_index(6,5)
        Traceback (most recent call last):
        ...
        IndexError: index out of range
        sage: normalize_index((4,-6),5)
        Traceback (most recent call last):
        ...
        IndexError: index out of range
        sage: normalize_index((-2,3),5)
        [3, 3]
        sage: normalize_index((5,0),5)
        Traceback (most recent call last):
        ...
        IndexError: index out of range
        sage: normalize_index((-5,2),5)
        [0, 2]
        sage: normalize_index((0,-2),5)
        [0, 3]
        sage: normalize_index((2,-3),5)
        [2, 2]
        sage: normalize_index((3,3),5)
        [3, 3]
        sage: normalize_index((-2,-5),5)
        [3, 0]
        sage: normalize_index((-2,-4),5)
        [3, 1]
        sage: normalize_index([-2,-1,3],5)
        [3, 4, 3]
        sage: normalize_index([4,2,1],5)
        [4, 2, 1]
        sage: normalize_index([-2,-3,-4],5)
        [3, 2, 1]
        sage: normalize_index([3,-2,-3],5)
        [3, 3, 2]
        sage: normalize_index([-5,2,-3],5)
        [0, 2, 2]
        sage: normalize_index([4,4,-5],5)
        [4, 4, 0]
        sage: s=slice(None,None,None); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(None,None,-2); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(None,None,4); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(None,-2,None); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(None,-2,-2); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(None,-2,4); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(None,4,None); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(None,4,-2); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(None,4,4); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(-2,None,None); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(-2,None,-2); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(-2,None,4); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(-2,-2,None); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(-2,-2,-2); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(-2,-2,4); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(-2,4,None); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(-2,4,-2); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(-2,4,4); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(4,None,None); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(4,None,-2); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(4,None,4); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(4,-2,None); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(4,-2,-2); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(4,-2,4); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(4,4,None); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(4,4,-2); normalize_index(s,5)==list(range(5))[s]
        True
        sage: s=slice(4,4,4); normalize_index(s,5)==list(range(5))[s]
        True
    """
    cdef tuple index_tuple
    cdef list return_list = []
    cdef Py_ssize_t index
    cdef Py_ssize_t i
    cdef object index_obj

    if PyIndex_Check(key):
        index = key
        if index < 0:
            index += size
        if index < 0 or index >= size:
            raise IndexError("index out of range")
        return [index]
    elif isinstance(key, slice):
        return list(range(*key.indices(size)))
    elif type(key) is tuple:
        index_tuple = key
    elif type(key) is list:
        index_tuple = PyList_AsTuple(key)
    elif type(key) is range:
        index_tuple = tuple(key)
    else:
        raise TypeError("index must be an integer or slice or a tuple/list of integers and slices")

    # Cython does not automatically use PyTuple_GET_SIZE, even though
    # it knows that index_tuple is tuple
    for i in range(PyTuple_GET_SIZE(index_tuple)):
        index_obj = index_tuple[i]
        if PyIndex_Check(index_obj):
            index = index_obj
            if index < 0:
                index += size
            if index < 0 or index >= size:
                raise IndexError("index out of range")
            return_list.append(index)
        elif isinstance(index_obj, slice):
            return_list.extend(range(*index_obj.indices(size)))
        else:
            raise TypeError("index must be an integer or slice")
    return return_list


cdef class sized_iter:
    """
    Wrapper for an iterator to verify that it has a specified length.

    INPUT:

    - ``iterable`` -- object to be iterated over

    - ``length`` -- (optional) the required length; if this is not
      given, then ``len(iterable)`` will be used

    If the iterable does not have the given length, a :exc:`ValueError` is
    raised during iteration.

    EXAMPLES::

        sage: from sage.misc.misc_c import sized_iter
        sage: list(sized_iter(range(4)))
        [0, 1, 2, 3]
        sage: list(sized_iter(range(4), 4))
        [0, 1, 2, 3]
        sage: list(sized_iter(range(5), 4))
        Traceback (most recent call last):
        ...
        ValueError: sequence too long (expected length 4, got more)
        sage: list(sized_iter(range(3), 4))
        Traceback (most recent call last):
        ...
        ValueError: sequence too short (expected length 4, got 3)

    If the iterable is too long, we get the error on the last entry::

        sage: it = sized_iter(range(5), 2)
        sage: next(it)
        0
        sage: next(it)
        Traceback (most recent call last):
        ...
        ValueError: sequence too long (expected length 2, got more)

    When the expected length is zero, the iterator is checked on
    construction::

        sage: list(sized_iter([], 0))
        []
        sage: sized_iter([1], 0)
        Traceback (most recent call last):
        ...
        ValueError: sequence too long (expected length 0, got more)

    If no ``length`` is given, the iterable must implement ``__len__``::

        sage: sized_iter(x for x in range(4))
        Traceback (most recent call last):
        ...
        TypeError: object of type 'generator' has no len()
    """
    cdef iterator
    cdef Py_ssize_t index, size

    def __init__(self, iterable, length=None):
        self.iterator = iter(iterable)
        self.index = 0
        if length is None:
            self.size = len(iterable)
        else:
            self.size = length
        self.check()

    def __iter__(self):
        return self

    def __len__(self):
        """
        Number of entries remaining, assuming that the expected length
        is the actual length.

        EXAMPLES::

            sage: from sage.misc.misc_c import sized_iter
            sage: it = sized_iter(range(4), 4)
            sage: len(it)
            4
            sage: next(it)
            0
            sage: len(it)
            3
        """
        return self.size - self.index

    cdef inline int check(self) except -1:
        """
        If the iterator is supposed to be exhausted, check that it is.
        """
        if self.index < self.size:
            return 0
        try:
            next(self.iterator)
        except StopIteration:
            pass
        else:
            raise ValueError(f"sequence too long (expected length {self.size}, got more)")

    def __next__(self):
        if self.index >= self.size:
            raise StopIteration
        try:
            x = next(self.iterator)
        except StopIteration:
            raise ValueError(f"sequence too short (expected length {self.size}, got {self.index})")
        self.index += 1
        self.check()
        return x


def cyflush():
    """
    Flush any output left over from external library calls.

    Starting with Python 3, some output from external libraries (like
    FLINT) is not flushed, and so if a doctest produces such output,
    the output may not appear until a later doctest. See
    :issue:`28649`.

    Use this function after a doctest which produces potentially
    unflushed output to force it to be flushed.

    EXAMPLES::

        sage: R.<t> = QQ[]
        sage: t^(sys.maxsize//2)                                                        # needs sage.libs.flint
        Traceback (most recent call last):
        ...
        RuntimeError: FLINT exception
        sage: from sage.misc.misc_c import cyflush
        sage: cyflush()
        ...
    """
    fflush(NULL)
