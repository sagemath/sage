"""
Base class for matrices, part 0

.. NOTE::

   For design documentation see matrix/docs.py.

EXAMPLES::

    sage: matrix(2, [1,2,3,4])
    [1 2]
    [3 4]
"""
# ****************************************************************************
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython cimport *
from cysignals.signals cimport sig_check

import sage.modules.free_module
import sage.misc.latex
import sage.rings.integer

from sage.arith.power cimport generic_power
from sage.structure.sequence import Sequence
from sage.structure.parent cimport Parent

cimport sage.structure.element
from sage.structure.element cimport Element, Vector
from sage.misc.misc_c cimport normalize_index

from sage.categories.fields import Fields
from sage.categories.integral_domains import IntegralDomains
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.rings import Rings

import sage.rings.abc
from sage.rings.integer_ring import IntegerRing_class

import sage.modules.free_module

from sage.matrix.matrix_misc import row_iterator

_Fields = Fields()
_IntegralDomains = IntegralDomains()

cdef class Matrix(sage.structure.element.Matrix):
    r"""
    A generic matrix.

    The ``Matrix`` class is the base class for all matrix
    classes. To create a ``Matrix``, first create a
    ``MatrixSpace``, then coerce a list of elements into
    the ``MatrixSpace``. See the documentation of
    ``MatrixSpace`` for more details.

    EXAMPLES:

    We illustrate matrices and matrix spaces. Note that no actual
    matrix that you make should have class Matrix; the class should
    always be derived from Matrix.

    ::

        sage: M = MatrixSpace(CDF,2,3); M
        Full MatrixSpace of 2 by 3 dense matrices over Complex Double Field
        sage: a = M([1,2,3,  4,5,6]); a
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        sage: type(a)
        <class 'sage.matrix.matrix_complex_double_dense.Matrix_complex_double_dense'>
        sage: parent(a)
        Full MatrixSpace of 2 by 3 dense matrices over Complex Double Field

    ::

        sage: matrix(CDF, 2,3, [1,2,3, 4,5,6])
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        sage: Mat(CDF,2,3)(range(1,7))
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]

    ::

        sage: Q.<i,j,k> = QuaternionAlgebra(QQ, -1,-1)
        sage: matrix(Q,2,1,[1,2])
        [1]
        [2]
    """
    def __cinit__(self, parent, *args, **kwds):
        """
        The initialization routine of the ``Matrix`` base class ensures
        that it sets the attributes ``self._parent``, ``self._base_ring``,
        ``self._nrows``, ``self._ncols``.

        The private attributes ``self._is_immutable`` and ``self._cache``
        are implicitly initialized to valid values upon memory allocation.

        EXAMPLES::

            sage: import sage.matrix.matrix0
            sage: A = sage.matrix.matrix0.Matrix(MatrixSpace(QQ,2))
            sage: type(A)
            <class 'sage.matrix.matrix0.Matrix'>
        """
        P = <Parent?>parent
        self._parent = P
        self._base_ring = P._base
        self._nrows = P.nrows()
        self._ncols = P.ncols()
        self.hash = -1

    def list(self):
        """
        List of the elements of ``self`` ordered by elements in each
        row. It is safe to change the returned list.

        .. warning::

           This function returns a list of the entries in the matrix
           ``self``.  It does not return a list of the rows of ``self``,
           so it is different than the output of ``list(self)``, which
           returns ``[self[0],self[1],...]``.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: a = matrix(R,2,[x,y,x*y, y,x,2*x+y]); a
            [      x       y     x*y]
            [      y       x 2*x + y]
            sage: v = a.list(); v
            [x, y, x*y, y, x, 2*x + y]

        Note that list(a) is different than a.list()::

            sage: a.list()
            [x, y, x*y, y, x, 2*x + y]
            sage: list(a)
            [(x, y, x*y), (y, x, 2*x + y)]

        Notice that changing the returned list does not change a (the list
        is a copy)::

            sage: v[0] = 25
            sage: a
            [      x       y     x*y]
            [      y       x 2*x + y]
        """
        return list(self._list())

    def dense_coefficient_list(self, order=None) -> list:
        """
        Return a list of *all* coefficients of ``self``.

        By default, this is the same as :meth:`list`.

        INPUT:

        - ``order`` -- (optional) an ordering of the basis indexing set

        EXAMPLES::

            sage: A = matrix([[1,2,3], [4,5,6]])
            sage: A.dense_coefficient_list()
            [1, 2, 3, 4, 5, 6]
            sage: A.dense_coefficient_list([(1,2), (1,0), (0,1), (0,2), (0,0), (1,1)])
            [6, 4, 2, 3, 1, 5]
        """
        if order is None:
            return self.list()
        return [self[i] for i in order]

    def _list(self):
        """
        Unsafe version of the ``list`` method, mainly for internal use.
        This may return the list of elements, but as an *unsafe* reference
        to the underlying list of the object. It is dangerous to change
        entries of the returned list.

        EXAMPLES:

        Using ``_list`` is potentially fast and memory efficient,
        but very dangerous (at least for generic dense matrices).

        ::

            sage: a = matrix(QQ['x,y'],2,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: v = a._list(); v
            [0, 1, 2, 3, 4, 5]

        If you change an entry of the list, the corresponding entry of the
        matrix will be changed (but without clearing any caches of
        computing information about the matrix)::

            sage: v[0] = -2/3; v
            [-2/3, 1, 2, 3, 4, 5]
            sage: a._list()
            [-2/3, 1, 2, 3, 4, 5]

        Now the 0,0 entry of the matrix is `-2/3`, which is weird.

        ::

            sage: a[0,0]
            -2/3

        See::

            sage: a
            [-2/3    1    2]
            [   3    4    5]
        """
        cdef Py_ssize_t i, j

        x = self.fetch('list')
        if x is not None:
            return x
        x = []
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                x.append(self.get_unsafe(i, j))
        return x

    def dict(self, copy=True):
        r"""
        Dictionary of the elements of ``self`` with keys pairs ``(i,j)``
        and values the nonzero entries of ``self``.

        INPUT:

        - ``copy`` -- boolean (default: ``True``); make a copy of the ``dict``
          corresponding to ``self``

        If ``copy=True``, then is safe to change the returned dictionary.
        Otherwise, this can cause undesired behavior by mutating the ``dict``.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: a = matrix(R,2,[x,y,0, 0,0,2*x+y]); a
            [      x       y       0]
            [      0       0 2*x + y]
            sage: d = a.dict(); d
            {(0, 0): x, (0, 1): y, (1, 2): 2*x + y}

        Notice that changing the returned list does not change a (the list
        is a copy)::

            sage: d[0,0] = 25
            sage: a
            [      x       y       0]
            [      0       0 2*x + y]
        """
        if copy:
            return dict(self._dict())
        return self._dict()

    monomial_coefficients = dict

    def items(self):
        r"""
        Return an iterable of ``((i,j), value)`` elements.

        This may (but is not guaranteed to) suppress zero values.

        EXAMPLES::

            sage: a = matrix(QQ['x,y'], 2, range(6), sparse=True); a
            [0 1 2]
            [3 4 5]
            sage: list(a.items())
            [((0, 1), 1), ((0, 2), 2), ((1, 0), 3), ((1, 1), 4), ((1, 2), 5)]
        """
        return self._dict().items()

    def _dict(self):
        """
        Unsafe version of the dict method, mainly for internal use.
        This may return the dict of elements, but as an *unsafe*
        reference to the underlying dict of the object.  It might
        dangerous if you change entries of the returned dict.

        EXAMPLES: Using _dict is potentially fast and memory efficient,
        but very dangerous (at least for generic sparse matrices).

        ::

            sage: a = matrix(QQ['x,y'],2,range(6), sparse=True); a
            [0 1 2]
            [3 4 5]
            sage: v = a._dict(); v
            {(0, 1): 1, (0, 2): 2, (1, 0): 3, (1, 1): 4, (1, 2): 5}

        If you change a key of the dictionary, the corresponding entry of
        the matrix will be changed (but without clearing any caches of
        computing information about the matrix)::

            sage: v[0,1] = -2/3; v
            {(0, 1): -2/3, (0, 2): 2, (1, 0): 3, (1, 1): 4, (1, 2): 5}
            sage: a._dict()
            {(0, 1): -2/3, (0, 2): 2, (1, 0): 3, (1, 1): 4, (1, 2): 5}
            sage: a[0,1]
            -2/3

        But the matrix doesn't know the entry changed, so it returns the
        cached version of its print representation::

            sage: a
            [0 1 2]
            [3 4 5]

        If we change an entry, the cache is cleared, and the correct print
        representation appears::

            sage: a[1,2]=10
            sage: a
            [   0 -2/3    2]
            [   3    4   10]
        """
        d = self.fetch('dict')
        if d is not None:
            return d

        cdef Py_ssize_t i, j
        d = {}
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                x = self.get_unsafe(i, j)
                if x != 0:
                    d[(int(i),int(j))] = x
        self.cache('dict', d)
        return d

    ###########################################################
    # Cache
    ###########################################################
    def _clear_cache(self):
        """
        Clear anything cached about this matrix.

        EXAMPLES::

            sage: m = Matrix(QQ, 2, range(4))
            sage: m._clear_cache()
        """
        self.clear_cache()

    cdef void clear_cache(self) noexcept:
        """
        Clear the properties cache.
        """
        self._cache = None
        self.hash = -1

    cdef fetch(self, key):
        """
        Try to get an element from the cache; if there isn't anything
        there, return None.
        """
        if self._cache is None:
            return None
        try:
            return self._cache[key]
        except KeyError:
            return None

    cdef cache(self, key, x):
        """
        Record x in the cache with given key.
        """
        if self._cache is None:
            self._cache = {}
        self._cache[key] = x

    def _get_cache(self):
        """
        Return the cache.

        EXAMPLES::

            sage: m=Matrix(QQ,2,range(0,4))
            sage: m._get_cache()
            {}
        """
        if self._cache is None:
            self._cache = {}
        return self._cache

    ###########################################################
    # Mutability and bounds checking
    ###########################################################

    cdef check_bounds(self, Py_ssize_t i, Py_ssize_t j):
        """
        This function gets called when you're about to access the i,j entry
        of this matrix. If i, j are out of range, an :exc:`IndexError` is
        raised.
        """
        if i < 0 or i >= self._nrows or j < 0 or j >= self._ncols:
            raise IndexError("matrix index out of range")

    cdef check_mutability(self):
        """
        This function gets called when you're about to change this matrix.

        If ``self`` is immutable, a :exc:`ValueError` is raised, since you should
        never change a mutable matrix.

        If ``self`` is mutable, the cache of results about ``self`` is deleted.
        """
        if self._is_immutable:
            raise ValueError("matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).")
        else:
            self._cache = None

    cdef check_bounds_and_mutability(self, Py_ssize_t i, Py_ssize_t j):
        """
        This function gets called when you're about to set the i,j entry of
        this matrix. If i or j is out of range, an :exc:`IndexError`
        exception is raised.

        If ``self`` is immutable, a :exc:`ValueError` is raised, since you should
        never change a mutable matrix.

        If ``self`` is mutable, the cache of results about ``self`` is deleted.
        """
        if self._is_immutable:
            raise ValueError("matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).")
        else:
            self._cache = None

        if i < 0 or i >= self._nrows or j < 0 or j >= self._ncols:
            raise IndexError("matrix index out of range")

    def set_immutable(self):
        r"""
        Call this function to set the matrix as immutable.

        Matrices are always mutable by default, i.e., you can change their
        entries using ``A[i,j] = x``. However, mutable matrices
        aren't hashable, so can't be used as keys in dictionaries, etc.
        Also, often when implementing a class, you might compute a matrix
        associated to it, e.g., the matrix of a Hecke operator. If you
        return this matrix to the user you're really returning a reference
        and the user could then change an entry; this could be confusing.
        Thus you should set such a matrix immutable.

        EXAMPLES::

            sage: A = Matrix(QQ, 2, 2, range(4))
            sage: A.is_mutable()
            True
            sage: A[0,0] = 10
            sage: A
            [10   1]
            [ 2   3]

        Mutable matrices are not hashable, so can't be used as keys for
        dictionaries::

            sage: hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            sage: v = {A:1}
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable

        If we make A immutable it suddenly is hashable.

        ::

            sage: A.set_immutable()
            sage: A.is_mutable()
            False
            sage: A[0,0] = 10
            Traceback (most recent call last):
            ...
            ValueError: matrix is immutable; please change a copy instead
            (i.e., use copy(M) to change a copy of M).
            sage: hash(A) #random
            12
            sage: v = {A:1}; v
            {[10  1]
             [ 2  3]: 1}
        """
        self._is_immutable = True

    def is_immutable(self):
        """
        Return ``True`` if this matrix is immutable.

        See the documentation for self.set_immutable for more details
        about mutability.

        EXAMPLES::

            sage: A = Matrix(QQ['t','s'], 2, 2, range(4))
            sage: A.is_immutable()
            False
            sage: A.set_immutable()
            sage: A.is_immutable()
            True
        """
        return self._is_immutable

    def is_mutable(self):
        """
        Return ``True`` if this matrix is mutable.

        See the documentation for self.set_immutable for more details
        about mutability.

        EXAMPLES::

            sage: A = Matrix(QQ['t','s'], 2, 2, range(4))
            sage: A.is_mutable()
            True
            sage: A.set_immutable()
            sage: A.is_mutable()
            False
        """
        return not(self._is_immutable)

    ###########################################################
    # Entry access
    #    The first two must be overloaded in the derived class
    ###########################################################
    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x):
        """
        Set entry quickly without doing any bounds checking. Calling this
        with invalid arguments is allowed to produce a segmentation fault.

        This is fast since it is a cdef function and there is no bounds
        checking.
        """
        raise NotImplementedError("this must be defined in the derived class (type=%s)" % type(self))

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Entry access, but fast since it might be without bounds checking.

        This is fast since it is a cdef function and there is no bounds
        checking.
        """
        raise NotImplementedError("this must be defined in the derived type.")

    cdef bint get_is_zero_unsafe(self, Py_ssize_t i, Py_ssize_t j) except -1:
        """
        Return 1 if the entry ``(i, j)`` is zero, otherwise 0.

        Might/should be optimized for derived type.

        TESTS::

            sage: class MyAlgebraicNumber(sage.rings.qqbar.AlgebraicNumber):            # needs sage.rings.number_field
            ....:     def __bool__(self):
            ....:         raise ValueError
            sage: mat = matrix(1, 1, MyAlgebraicNumber(1))                              # needs sage.rings.number_field
            sage: bool(mat)                                                             # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            ValueError
        """
        if self.get_unsafe(i, j):
            return 0
        else:
            return 1

    def add_to_entry(self, Py_ssize_t i, Py_ssize_t j, elt):
        r"""
        Add ``elt`` to the entry at position ``(i, j)``.

        EXAMPLES::

            sage: m = matrix(QQ['x,y'], 2, 2)
            sage: m.add_to_entry(0, 1, 2)
            sage: m
            [0 2]
            [0 0]
        """
        elt = self.base_ring()(elt)
        if i < 0:
            i += self._nrows
        if i < 0 or i >= self._nrows:
            raise IndexError("row index out of range")
        if j < 0:
            j += self._ncols
        if j < 0 or j >= self._ncols:
            raise IndexError("column index out of range")

        self.set_unsafe(i, j, elt + self.get_unsafe(i, j))

    ##     def _get_very_unsafe(self, i, j):
    ##         r"""
    ##         Entry access, but potentially fast since it might be without
    ##         bounds checking.  (I know of no cases where this is actually
    ##         faster.)

    ##         This function it can very easily !! SEG FAULT !! if you call
    ##         it with invalid input.  Use with *extreme* caution.

    ##         EXAMPLES::
    ##
    ##             sage: a = matrix(ZZ,2,range(4))
    ##             sage: a._get_very_unsafe(0,1)
    ##             1

    ##         If you do \code{a.\_get\_very\_unsafe(0,10)} you'll very likely crash Sage
    ##         completely.
    ##         """
    ##         return self.get_unsafe(i, j)

    def __iter__(self):
        """
        Return an iterator for the rows of ``self``.

        EXAMPLES::

            sage: m = matrix(2,[1,2,3,4])
            sage: next(m.__iter__())
            (1, 2)
        """
        return row_iterator(self)

    def __getitem__(self, key):
        """
        Return element, row, or slice of ``self``.

        INPUT:

        - ``key``- tuple (i,j) where i, j can be integers, slices or lists

        USAGE:

        - ``A[i, j]`` -- the i,j element (or elements, if i or j are
          slices or lists) of A, or

        - ``A[i:j]`` -- rows of A, according to slice notation

        EXAMPLES::

            sage: A = Matrix(Integers(2006),2,2,[-1,2,3,4])
            sage: A[0,0]
            2005
            sage: A[0]
            (2005, 2)

        The returned row is immutable (mainly to avoid confusion)::

            sage: A[0][0] = 123
            Traceback (most recent call last):
            ...
            ValueError: vector is immutable; please change a copy instead (use copy())
            sage: A[0].is_immutable()
            True
            sage: a = matrix(ZZ,3,range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a[1,2]
            5
            sage: a[0]
            (0, 1, 2)
            sage: a[4,7]
            Traceback (most recent call last):
            ...
            IndexError: matrix index out of range
            sage: a[-1,0]
            6

        ::

            sage: a[2.7]
            Traceback (most recent call last):
            ...
            TypeError: index must be an integer
            sage: a[1, 2.7]
            Traceback (most recent call last):
            ...
            TypeError: index must be an integer
            sage: a[2.7, 1]
            Traceback (most recent call last):
            ...
            TypeError: index must be an integer

            sage: m = [(1, -2, -1, -1,9), (1, 8, 6, 2,2), (1, 1, -1, 1,4), (-1, 2, -2, -1,4)]; M = matrix(m)
            sage: M
            [ 1 -2 -1 -1  9]
            [ 1  8  6  2  2]
            [ 1  1 -1  1  4]
            [-1  2 -2 -1  4]

        Get the 2 x 2 submatrix of M, starting at row index and column
        index 1

        ::

            sage: M[1:3,1:3]
            [ 8  6]
            [ 1 -1]

        Get the 2 x 3 submatrix of M starting at row index and column index
        1::

            sage: M[1:3,[1..3]]
            [ 8  6  2]
            [ 1 -1  1]

        Get the second column of M::

            sage: M[:,1]
            [-2]
            [ 8]
            [ 1]
            [ 2]

        Get the first row of M::

            sage: M[0,:]
            [ 1 -2 -1 -1  9]

        More examples::

            sage: M[range(2),:]
            [ 1 -2 -1 -1  9]
            [ 1  8  6  2  2]
            sage: M[range(2),4]
            [9]
            [2]
            sage: M[range(3),range(5)]
            [ 1 -2 -1 -1  9]
            [ 1  8  6  2  2]
            [ 1  1 -1  1  4]

        ::

            sage: M[3,range(5)]
            [-1  2 -2 -1  4]
            sage: M[3,:]
            [-1  2 -2 -1  4]
            sage: M[3,4]
            4

            sage: M[-1,:]
            [-1  2 -2 -1  4]

            sage: A = matrix(ZZ,3,4, [3, 2, -5, 0, 1, -1, 1, -4, 1, 0, 1, -3]); A
            [ 3  2 -5  0]
            [ 1 -1  1 -4]
            [ 1  0  1 -3]

        ::

            sage: A[:,0:4:2]
            [ 3 -5]
            [ 1  1]
            [ 1  1]

        ::

            sage: A[1:,0:4:2]
            [1 1]
            [1 1]

            sage: A[2::-1,:]
            [ 1  0  1 -3]
            [ 1 -1  1 -4]
            [ 3  2 -5  0]

            sage: A[1:,3::-1]
            [-4  1 -1  1]
            [-3  1  0  1]

            sage: A[1:,3::-2]
            [-4 -1]
            [-3  0]

            sage: A[2::-1,3:1:-1]
            [-3  1]
            [-4  1]
            [ 0 -5]

        ::

            sage: A= matrix(3,4,[1, 0, -3, -1, 3, 0, -2, 1, -3, -5, -1, -5])
            sage: A[range(2,-1,-1),:]
            [-3 -5 -1 -5]
            [ 3  0 -2  1]
            [ 1  0 -3 -1]

        ::

            sage: A[range(2,-1,-1),range(3,-1,-1)]
            [-5 -1 -5 -3]
            [ 1 -2  0  3]
            [-1 -3  0  1]

        ::

            sage: A = matrix(2, [1, 2, 3, 4])
            sage: A[[0,0],[0,0]]
            [1 1]
            [1 1]

        ::

            sage: M = matrix(3, 4, range(12))
            sage: M[0:0, 0:0]
            []
            sage: M[0:0, 1:4]
            []
            sage: M[2:3, 3:3]
            []
            sage: M[range(2,2), :3]
            []
            sage: M[(1,2), 3]
            [ 7]
            [11]
            sage: M[(1,2),(0,1,1)]
            [4 5 5]
            [8 9 9]
            sage: m=[(1, -2, -1, -1), (1, 8, 6, 2), (1, 1, -1, 1), (-1, 2, -2, -1)]
            sage: M= matrix(m);M
            [ 1 -2 -1 -1]
            [ 1  8  6  2]
            [ 1  1 -1  1]
            [-1  2 -2 -1]

            sage: M[:2]
            [ 1 -2 -1 -1]
            [ 1  8  6  2]
            sage: M[:]
            [ 1 -2 -1 -1]
            [ 1  8  6  2]
            [ 1  1 -1  1]
            [-1  2 -2 -1]
            sage: M[1:3]
            [ 1  8  6  2]
            [ 1  1 -1  1]

            sage: A=matrix(QQ,10,range(100))
            sage: A[0:3]
            [ 0  1  2  3  4  5  6  7  8  9]
            [10 11 12 13 14 15 16 17 18 19]
            [20 21 22 23 24 25 26 27 28 29]
            sage: A[:2]
            [ 0  1  2  3  4  5  6  7  8  9]
            [10 11 12 13 14 15 16 17 18 19]
            sage: A[8:]
            [80 81 82 83 84 85 86 87 88 89]
            [90 91 92 93 94 95 96 97 98 99]
            sage: A[1:10:3]
            [10 11 12 13 14 15 16 17 18 19]
            [40 41 42 43 44 45 46 47 48 49]
            [70 71 72 73 74 75 76 77 78 79]
            sage: A[-1]
            (90, 91, 92, 93, 94, 95, 96, 97, 98, 99)
            sage: A[-1:-6:-2]
            [90 91 92 93 94 95 96 97 98 99]
            [70 71 72 73 74 75 76 77 78 79]
            [50 51 52 53 54 55 56 57 58 59]

            sage: A[3].is_immutable()
            True
            sage: A[1:3].is_immutable()
            True

        Slices that result in zero rows or zero columns are supported too::

            sage: m = identity_matrix(QQ, 4)[4:,:]
            sage: m.nrows(), m.ncols()
            (0, 4)
            sage: m * vector(QQ, 4)
            ()

        TESTS:

        If we're given lists as arguments, we should throw an
        appropriate error when those lists do not contain valid
        indices (:issue:`6569`)::

            sage: A = matrix(4, range(1,17))
            sage: A[[1.5], [1]]
            Traceback (most recent call last):
            ...
            IndexError: row indices must be integers
            sage: A[[1], [1.5]]
            Traceback (most recent call last):
            ...
            IndexError: column indices must be integers
            sage: A[[1.5]]
            Traceback (most recent call last):
            ...
            IndexError: row indices must be integers

        Before :issue:`6569` was fixed, sparse/dense matrices behaved
        differently due to implementation details. Given invalid
        indices, they should fail in the same manner. These tests
        just repeat the previous set with a sparse matrix::

            sage: A = matrix(4, range(1,17), sparse=True)
            sage: A[[1.5], [1]]
            Traceback (most recent call last):
            ...
            IndexError: row indices must be integers
            sage: A[[1], [1.5]]
            Traceback (most recent call last):
            ...
            IndexError: column indices must be integers
            sage: A[[1.5]]
            Traceback (most recent call last):
            ...
            IndexError: row indices must be integers

        Check that submatrices with a specified implementation have the
        same implementation::

            sage: # needs sage.libs.pari
            sage: M = MatrixSpace(GF(2), 3, 3, implementation='generic')
            sage: m = M(range(9))
            sage: type(m)
            <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
            sage: parent(m)
            Full MatrixSpace of 3 by 3 dense matrices
             over Finite Field of size 2 (using Matrix_generic_dense)
            sage: type(m[:2,:2])
            <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
            sage: parent(m[:2,:2])
            Full MatrixSpace of 2 by 2 dense matrices
             over Finite Field of size 2 (using Matrix_generic_dense)
        """
        cdef list row_list
        cdef list col_list
        cdef Py_ssize_t i
        cdef int row, col
        cdef int nrows = self._nrows
        cdef int ncols = self._ncols
        cdef tuple key_tuple
        cdef object row_index, col_index
        cdef int ind

        # used to keep track of when an index is a
        # single number
        cdef int single_row = 0, single_col = 0

        if type(key) is tuple:
            key_tuple = <tuple>key
            #if PyTuple_Size(key_tuple) != 2:
            if len(key_tuple) != 2:
                raise IndexError("index must be an integer or pair of integers")

            row_index = <object>PyTuple_GET_ITEM(key_tuple, 0)
            col_index = <object>PyTuple_GET_ITEM(key_tuple, 1)

            type_row = type(row_index)
            if type_row is list or type_row is tuple or type_row is range:
                if type_row is tuple or type_row is range:
                    row_list = list(row_index)
                else:
                    row_list = row_index

                for i from 0 <= i < len(row_list):
                    # The 'ind' variable is 'cdef int' and will
                    # truncate a float to a valid index. So, we have
                    # to test row_list[i] instead.
                    if not PyIndex_Check(row_list[i]):
                        raise IndexError('row indices must be integers')

                    ind = row_list[i]
                    if ind < 0:
                        ind += nrows
                        row_list[i] = ind

                    if ind < 0 or ind >= nrows:
                        raise IndexError("matrix index out of range")
            elif isinstance(row_index, slice):
                row_list = list(range(*row_index.indices(nrows)))
            else:
                if not PyIndex_Check(row_index):
                    raise TypeError("index must be an integer")
                row = row_index
                if row < 0:
                    row += nrows
                if row < 0 or row >= nrows:
                    raise IndexError("matrix index out of range")
                single_row = 1

            type_col = type(col_index)
            if type_col is list or type_col is tuple or type_col is range:
                if type_col is tuple or type_col is range:
                    col_list = list(col_index)
                else:
                    col_list = col_index

                for i from 0 <= i < len(col_list):
                    # The 'ind' variable is 'cdef int' and will
                    # truncate a float to a valid index. So, we have
                    # to test col_list[i] instead.
                    if not PyIndex_Check(col_list[i]):
                        raise IndexError('column indices must be integers')

                    ind = col_list[i]
                    if ind < 0:
                        ind += ncols
                        col_list[i] = ind

                    if ind < 0 or ind >= ncols:
                        raise IndexError("matrix index out of range")
            elif isinstance(col_index, slice):
                col_list =  list(range(*col_index.indices(ncols)))
            else:
                if not PyIndex_Check(col_index):
                    raise TypeError("index must be an integer")
                col = col_index
                if col < 0:
                    col += ncols
                if col < 0 or col >= ncols:
                    raise IndexError("matrix index out of range")
                single_col = 1

            # if we had a single row entry and a single column entry,
            # we want to just do a get_unsafe
            if single_row and single_col:
                return self.get_unsafe(row, col)

            # otherwise, prep these for the call to
            # matrix_from_rows_and_columns
            if single_row:
                row_list = [row]
            if single_col:
                col_list = [col]

            if not row_list or not col_list:
                return self.new_matrix(nrows=len(row_list), ncols=len(col_list))

            return self.matrix_from_rows_and_columns(row_list, col_list)

        row_index = key
        if type(row_index) is list or type(row_index) is tuple:
            if type(row_index) is tuple:
                row_list = list(row_index)
            else:
                row_list = row_index

            for i from 0 <= i < len(row_list):
                # The 'ind' variable is 'cdef int' and will
                # truncate a float to a valid index. So, we have
                # to test row_list[i] instead.
                if not PyIndex_Check(row_list[i]):
                    raise IndexError('row indices must be integers')

                ind = row_list[i]
                if ind < 0:
                    ind += nrows
                    row_list[i] = ind

                if ind < 0 or ind >= nrows:
                    raise IndexError("matrix index out of range")
            r = self.matrix_from_rows(row_list)
        elif isinstance(row_index, slice):
            row_list = list(range(*row_index.indices(nrows)))
            r = self.matrix_from_rows(row_list)
        else:
            if not PyIndex_Check(row_index):
                raise TypeError("index must be an integer")
            row = row_index
            if row < 0:
                row += nrows
            if row < 0 or row >= nrows:
                raise IndexError("matrix index out of range")
            r = self.row(row)

        r.set_immutable()
        return r

    def __setitem__(self, key, value):
        """
        Set elements of this matrix to values given in value.

        INPUT:

        - ``key`` -- any legal indexing (i.e., such that self[key] works)

        - ``value`` -- values that are used to set the elements indicated by key

        EXAMPLES::

            sage: A = Matrix(Integers(2006),2,2,[-1,2,3,4])
            sage: A[0,0]=43; A
            [43  2]
            [ 3  4]

            sage: A[0]=[10,20]; A
            [10 20]
            [ 3  4]

            sage: M=matrix([(1, -2, -1, -1,9), (1, 8, 6, 2,2), (1, 1, -1, 1,4), (-1, 2, -2, -1,4)]); M
            [ 1 -2 -1 -1  9]
            [ 1  8  6  2  2]
            [ 1  1 -1  1  4]
            [-1  2 -2 -1  4]

        Set the 2 x 2 submatrix of M, starting at row index and column
        index 1::

            sage: M[1:3,1:3] = [[1,0],[0,1]]; M
            [ 1 -2 -1 -1  9]
            [ 1  1  0  2  2]
            [ 1  0  1  1  4]
            [-1  2 -2 -1  4]

        Set the 2 x 3 submatrix of M starting at row index and column
        index 1::

            sage: M[1:3,[1..3]] = M[2:4,0:3]; M
            [ 1 -2 -1 -1  9]
            [ 1  1  0  1  2]
            [ 1 -1  2 -2  4]
            [-1  2 -2 -1  4]

        Set part of the first column of M::

            sage: M[1:,0]=[[2],[3],[4]]; M
            [ 1 -2 -1 -1  9]
            [ 2  1  0  1  2]
            [ 3 -1  2 -2  4]
            [ 4  2 -2 -1  4]

        Or do a similar thing with a vector::

            sage: M[1:,0]=vector([-2,-3,-4]); M
            [ 1 -2 -1 -1  9]
            [-2  1  0  1  2]
            [-3 -1  2 -2  4]
            [-4  2 -2 -1  4]

        Or a constant::

            sage: M[1:,0]=30; M
            [ 1 -2 -1 -1  9]
            [30  1  0  1  2]
            [30 -1  2 -2  4]
            [30  2 -2 -1  4]


        Set the first row of M::

            sage: M[0,:]=[[20,21,22,23,24]]; M
            [20 21 22 23 24]
            [30  1  0  1  2]
            [30 -1  2 -2  4]
            [30  2 -2 -1  4]
            sage: M[0,:]=vector([0,1,2,3,4]); M
            [ 0  1  2  3  4]
            [30  1  0  1  2]
            [30 -1  2 -2  4]
            [30  2 -2 -1  4]
            sage: M[0,:]=-3; M
            [-3 -3 -3 -3 -3]
            [30  1  0  1  2]
            [30 -1  2 -2  4]
            [30  2 -2 -1  4]


            sage: A = matrix(ZZ,3,4, [3, 2, -5, 0, 1, -1, 1, -4, 1, 0, 1, -3]); A
            [ 3  2 -5  0]
            [ 1 -1  1 -4]
            [ 1  0  1 -3]

        We can use the step feature of slices to set every other column::

            sage: A[:,0:3:2] = 5; A
            [ 5  2  5  0]
            [ 5 -1  5 -4]
            [ 5  0  5 -3]

            sage: A[1:,0:4:2] = [[100,200],[300,400]]; A
            [  5   2   5   0]
            [100  -1 200  -4]
            [300   0 400  -3]

        We can also count backwards to flip the matrix upside down.

        ::

            sage: A[::-1,:]=A; A
            [300   0 400  -3]
            [100  -1 200  -4]
            [  5   2   5   0]


            sage: A[1:,3::-1]=[[2,3,0,1],[9,8,7,6]]; A
            [300   0 400  -3]
            [  1   0   3   2]
            [  6   7   8   9]

            sage: A[1:,::-2] = A[1:,::2]; A
            [300   0 400  -3]
            [  1   3   3   1]
            [  6   8   8   6]

            sage: A[::-1,3:1:-1] = [[4,3],[1,2],[-1,-2]]; A
            [300   0  -2  -1]
            [  1   3   2   1]
            [  6   8   3   4]


        TESTS::

            sage: A = MatrixSpace(ZZ,3)(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A[1,2]=100; A
            [  0   1   2]
            [  3   4 100]
            [  6   7   8]
            sage: A[0]=(10,20,30); A
            [ 10  20  30]
            [  3   4 100]
            [  6   7   8]
            sage: A[4,7]=45
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: A[-1,0]=63; A[-1,0]
            63
            sage: A[2.7]=3
            Traceback (most recent call last):
            ...
            TypeError: index must be an integer or slice or a tuple/list of integers and slices
            sage: A[1, 2.7]=3
            Traceback (most recent call last):
            ...
            TypeError: index must be an integer or slice or a tuple/list of integers and slices
            sage: A[2.7, 1]=3
            Traceback (most recent call last):
            ...
            TypeError: index must be an integer or slice or a tuple/list of integers and slices
            sage: A.set_immutable()
            sage: A[0,0] = 7
            Traceback (most recent call last):
            ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).
            sage: A=matrix([[1,2],[3,4]]); B=matrix([[1,3],[5,7]])
            sage: A[1:2,1:2]=B[1:2,1:2]
            sage: A
            [1 2]
            [3 7]
            sage: A=matrix([[1,2],[3,4]]); B=matrix([[1,3],[5,7]])
            sage: A[1,0:1]=B[1,1:2]
            sage: A
            [1 2]
            [7 4]


        More examples::

            sage: M[range(2),:]=[[1..5], [6..10]]; M
            [ 1  2  3  4  5]
            [ 6  7  8  9 10]
            [30 -1  2 -2  4]
            [30  2 -2 -1  4]

            sage: M[range(2),4]=0; M
            [ 1  2  3  4  0]
            [ 6  7  8  9  0]
            [30 -1  2 -2  4]
            [30  2 -2 -1  4]

            sage: M[range(3),range(5)]=M[range(1,4), :]; M
            [ 6  7  8  9  0]
            [30 -1  2 -2  4]
            [30  2 -2 -1  4]
            [30  2 -2 -1  4]


            sage: M[3,range(5)]=vector([-2,3,4,-5,4]); M
            [ 6  7  8  9  0]
            [30 -1  2 -2  4]
            [30  2 -2 -1  4]
            [-2  3  4 -5  4]
            sage: M[3,:]=2*M[2,:]; M
            [ 6  7  8  9  0]
            [30 -1  2 -2  4]
            [30  2 -2 -1  4]
            [60  4 -4 -2  8]
            sage: M[3,4]=M[3,2]; M
            [ 6  7  8  9  0]
            [30 -1  2 -2  4]
            [30  2 -2 -1  4]
            [60  4 -4 -2 -4]

            sage: M[-1,:]=M[-3,:]; M
            [ 6  7  8  9  0]
            [30 -1  2 -2  4]
            [30  2 -2 -1  4]
            [30 -1  2 -2  4]


            sage: A= matrix(3,4,[1, 0, -3, -1, 3, 0, -2, 1, -3, -5, -1, -5]); A
            [ 1  0 -3 -1]
            [ 3  0 -2  1]
            [-3 -5 -1 -5]

            sage: A[range(2,-1,-1),:]=A; A
            [-3 -5 -1 -5]
            [ 3  0 -2  1]
            [ 1  0 -3 -1]

            sage: A[range(2,-1,-1),range(3,-1,-1)]=A; A
            [-1 -3  0  1]
            [ 1 -2  0  3]
            [-5 -1 -5 -3]

            sage: A = matrix(2, [1, 2, 3, 4])
            sage: A[[0,0],[0,0]]=10; A
            [10  2]
            [ 3  4]

            sage: M = matrix(3, 4, range(12))
            sage: M[0:0, 0:0]=20; M
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            sage: M[0:0, 1:4]=20; M
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            sage: M[2:3, 3:3]=20; M
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            sage: M[range(2,2), :3]=20; M
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            sage: M[(1,2), 3]=vector([-1,-2]); M
            [ 0  1  2  3]
            [ 4  5  6 -1]
            [ 8  9 10 -2]
            sage: M[(1,2),(0,1,1)]=[[-1,-2,-3],[-4,-5,-6]]; M
            [ 0  1  2  3]
            [-1 -3  6 -1]
            [-4 -6 10 -2]
            sage: M=matrix([(1, -2, -1, -1), (1, 8, 6, 2), (1, 1, -1, 1), (-1, 2, -2, -1)]); M
            [ 1 -2 -1 -1]
            [ 1  8  6  2]
            [ 1  1 -1  1]
            [-1  2 -2 -1]

            sage: M[:2]=M[2:]; M
            [ 1  1 -1  1]
            [-1  2 -2 -1]
            [ 1  1 -1  1]
            [-1  2 -2 -1]

            sage: M[:] = M.transpose(); M
            [ 1 -1  1 -1]
            [ 1  2  1  2]
            [-1 -2 -1 -2]
            [ 1 -1  1 -1]
            sage: M = matrix(ZZ,4,range(16)); M
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            [12 13 14 15]
            sage: M[::2]=M[::-2]; M
            [12 13 14 15]
            [ 4  5  6  7]
            [ 4  5  6  7]
            [12 13 14 15]
            sage: M[::2]=2; M
            [ 2  2  2  2]
            [ 4  5  6  7]
            [ 2  2  2  2]
            [12 13 14 15]

            sage: M[2:]=10; M
            [ 2  2  2  2]
            [ 4  5  6  7]
            [10 10 10 10]
            [10 10 10 10]

            sage: M=matrix(3,1,[1,2,3]); M
            [1]
            [2]
            [3]
            sage: M[1] = vector([20]); M
            [ 1]
            [20]
            [ 3]
            sage: M = matrix(3, 2, srange(6)); M[1] = 15; M
            [ 0  1]
            [15 15]
            [ 4  5]
            sage: M = matrix(3, 1, srange(3)); M[1] = 15; M
            [ 0]
            [15]
            [ 2]
            sage: M = matrix(3, 1, srange(3)); M[1] = [15]; M
            [ 0]
            [15]
            [ 2]
        """
        cdef list row_list
        cdef list col_list
        cdef Py_ssize_t row_list_len, col_list_len
        cdef list value_list
        cdef bint value_list_one_dimensional = 0
        cdef Py_ssize_t i
        cdef Py_ssize_t row, col
        cdef Py_ssize_t nrows = self._nrows
        cdef Py_ssize_t ncols = self._ncols
        cdef tuple key_tuple
        cdef object row_index, col_index
        cdef object value_row

        # used to keep track of when an index is a
        # single number
        cdef bint single_row = 0, single_col = 0
        cdef bint no_col_index = 0

        # If the matrix is immutable, check_mutability will raise an
        # exception.
        self.check_mutability()

        if type(key) is tuple:
            key_tuple = <tuple>key
            #if PyTuple_Size(key_tuple) != 2:
            if len(key_tuple) != 2:
                raise IndexError("index can't have more than two components")

            row_index = <object>PyTuple_GET_ITEM(key_tuple, 0)
            col_index = <object>PyTuple_GET_ITEM(key_tuple, 1)

            if PyIndex_Check(col_index):
                col = col_index
                if col < 0:
                    col += ncols
                if col < 0 or col >= ncols:
                    raise IndexError("index out of range")
                single_col = 1
                col_list_len = 1
            else:
                col_list = normalize_index(col_index, ncols)
                col_list_len = len(col_list)
                if col_list_len==0:
                    return

        else:
            no_col_index = 1
            row_index = key
            col_list_len = ncols
            if col_list_len==0:
                return

        # Special-case a single-row.
        if PyIndex_Check(row_index):
            row = row_index
            if row < 0:
                row += nrows
            if row < 0 or row >= nrows:
                raise IndexError("index out of range")
            single_row = 1
            row_list_len = 1
        else:
            row_list = normalize_index(row_index, nrows)
            row_list_len = len(row_list)
            if row_list_len==0:
                return

        if single_row and single_col and not no_col_index:
            self.set_unsafe(row, col, self._coerce_element(value))
            return

        if type(value) is list:
            if single_row and no_col_index:
                # A convenience addition, so we can set a row by
                # M[1] = [1,2,3] or M[1,:]=[1,2,3]
                value_list_one_dimensional = 1
            value_list = value
        elif type(value) is tuple:
            if single_row and no_col_index:
                # A convenience addition, so we can set a row by
                # M[1] = [1,2,3] or M[1,:]=[1,2,3]
                value_list_one_dimensional = 1
            value_list = list(value)
        elif isinstance(value, Matrix):
            value_list = list(value)
        elif isinstance(value, Vector):
            if single_row or single_col:
                value_list_one_dimensional = 1
                value_list = list(value)
            else:
                raise IndexError("value does not have the right dimensions")
        else:
            # If value is not a list, tuple, matrix, or vector, try
            # broadcasting the element to all positions.
            value_element = self._coerce_element(value)
            if single_row:
                if no_col_index:
                    for col in range(col_list_len):
                        self.set_unsafe(row, col, value_element)
                else:
                    for col in col_list:
                        self.set_unsafe(row, col, value_element)
            elif single_col:
                for row in row_list:
                    self.set_unsafe(row, col, value_element)
            else:
                if no_col_index:
                    for row in row_list:
                        for col in range(col_list_len):
                            self.set_unsafe(row, col, value_element)
                else:
                    for row in row_list:
                        for col in col_list:
                            self.set_unsafe(row, col, value_element)
            return

        if value_list_one_dimensional:
            # This will break when assigning a vector to a column
            if single_row and col_list_len != len(value_list):
                raise IndexError("value does not have the right number of columns")
            elif single_col and row_list_len != len(value_list):
                raise IndexError("value does not have the right number of rows")
        else:
            if row_list_len != len(value_list):
                raise IndexError("value does not have the right number of rows")
            for value_row in value_list:
                if col_list_len != len(value_row):
                    raise IndexError("value does not have the right number of columns")

        if single_row:
            if value_list_one_dimensional:
                value_row = value_list
            else:
                value_row = value_list[0]

            if no_col_index:
                for col in range(col_list_len):
                    self.set_unsafe(row, col, self._coerce_element(value_row[col]))
            else:
                for col in range(col_list_len):
                    self.set_unsafe(row, col_list[col], self._coerce_element(value_row[col]))
        elif single_col:
            if value_list_one_dimensional:
                for row in range(row_list_len):
                    self.set_unsafe(row_list[row], col, self._coerce_element(value_list[row]))
            else:
                for row in range(row_list_len):
                    self.set_unsafe(row_list[row], col, self._coerce_element(value_list[row][0]))
        else:
            if no_col_index:
                for i in range(row_list_len):
                    row = row_list[i]
                    value_row = value_list[i]
                    for col in range(col_list_len):
                        self.set_unsafe(row, col, self._coerce_element(value_row[col]))
            else:
                for i in range(row_list_len):
                    row = row_list[i]
                    value_row = value_list[i]
                    for col in range(col_list_len):
                        self.set_unsafe(row, col_list[col], self._coerce_element(value_row[col]))
        return

    cdef _coerce_element(self, x):
        """
        Return coercion of x into the base ring of ``self``.
        """
        if isinstance(x, Element) and (<Element> x)._parent is self._base_ring:
            return x
        return self._base_ring(x)

    ###########################################################
    # Pickling
    ###########################################################

    def __reduce__(self):
        """
        EXAMPLES::

            sage: a = matrix(Integers(8),3,range(9))
            sage: a == loads(dumps(a))
            True
        """
        data, version = self._pickle()
        return unpickle, (self.__class__, self._parent, self._is_immutable,
                                          self._cache, data, version)

    def _pickle(self):
        """
        Not yet implemented!

        EXAMPLES::

            sage: m=matrix(QQ,2,range(0,4))
            sage: m._pickle() # todo: not implemented
        """
        raise NotImplementedError

    def _test_reduce(self, **options):
        """
        Check that the pickling function works.

        EXAMPLES::

            sage: a=matrix([[1,2],[3,4]])
            sage: a._test_reduce()
        """
        tester = self._tester(**options)
        a, b = self.__reduce__()
        tester.assertEqual(a(*b),self)

    ###########################################################
    # Base Change
    ###########################################################
    def base_ring(self):
        """
        Return the base ring of the matrix.

        EXAMPLES::

            sage: m = matrix(QQ, 2, [1,2,3,4])
            sage: m.base_ring()
            Rational Field
        """
        return self._base_ring

    def change_ring(self, ring):
        """
        Return the matrix obtained by coercing the entries of this matrix
        into the given ring.

        Always returns a copy (unless ``self`` is immutable, in which case
        returns ``self``).

        EXAMPLES::

            sage: A = Matrix(QQ, 2, 2, [1/2, 1/3, 1/3, 1/4])
            sage: A.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: A.change_ring(GF(25,'a'))                                             # needs sage.rings.finite_rings
            [3 2]
            [2 4]
            sage: A.change_ring(GF(25,'a')).parent()                                    # needs sage.rings.finite_rings
            Full MatrixSpace of 2 by 2 dense matrices
             over Finite Field in a of size 5^2
            sage: A.change_ring(ZZ)                                                     # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            TypeError: matrix has denominators so can...t change to ZZ

        Changing rings preserves subdivisions::

            sage: A.subdivide([1], []); A
            [1/2 1/3]
            [-------]
            [1/3 1/4]
            sage: A.change_ring(GF(25,'a'))                                             # needs sage.rings.finite_rings
            [3 2]
            [---]
            [2 4]
        """
        if ring not in Rings():
            raise TypeError("ring must be a ring")

        if ring is self._base_ring:
            if self._is_immutable:
                return self
            return self.__copy__()

        try:
            return self._change_ring(ring)
        except (AttributeError, NotImplementedError):
            M = sage.matrix.matrix_space.MatrixSpace(ring, self._nrows, self._ncols, sparse=self.is_sparse())
            mat = M(self.list(), coerce=True, copy=False)
            if self._subdivisions is not None:
                mat.subdivide(self.subdivisions())
            return mat

    def _test_change_ring(self, **options):
        """
        Check that :meth:`change_ring` works.

        EXAMPLES::

            sage: a = matrix([[1,2],[3,4]])
            sage: a._test_change_ring()
        """
        tester = self._tester(**options)
        # Test to make sure the returned matrix is a copy
        tester.assertIsNot(self.change_ring(self.base_ring()), self)

    def _matrix_(self, R=None):
        """
        Return ``self`` as a matrix over the ring ``R``. If ``R`` is ``None``,
        then return ``self``.

        EXAMPLES::

            sage: A = Matrix(ZZ[['t']], 2, 2, range(4))
            sage: A.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Power Series Ring in t over Integer Ring
            sage: A._matrix_(QQ[['t']])
            [0 1]
            [2 3]
            sage: A._matrix_(QQ[['t']]).parent()
            Full MatrixSpace of 2 by 2 dense matrices over Power Series Ring in t over Rational Field

        Check that :issue:`14314` is fixed::

            sage: m = Matrix({(1,2):2})
            sage: matrix(m) == m
            True
        """
        if R is None:
            return self
        return self.change_ring(R)

    ###########################################################
    # Representation -- string, latex, etc.
    ###########################################################
    def __repr__(self):
        r"""
        EXAMPLES::

            sage: A = matrix([[1,2], [3,4], [5,6]])
            sage: A.__repr__()
            '[1 2]\n[3 4]\n[5 6]'
            sage: print(A)
            [1 2]
            [3 4]
            [5 6]

        If the matrix is too big, don't print all of the elements::

            sage: A = random_matrix(ZZ, 100)
            sage: A.__repr__()
            '100 x 100 dense matrix over Integer Ring'

        When a big matrix returned, include a hint on how to get the entries.
        This is a feature of the sage command-line::

            sage: A
            100 x 100 dense matrix over Integer Ring (use the '.str()' method to see the entries)

        But don't do that when the matrix is part of a larger data structure::

            sage: [A]*2
            [100 x 100 dense matrix over Integer Ring,
             100 x 100 dense matrix over Integer Ring]
        """
        from sage.matrix.constructor import options
        if self._nrows <= options.max_rows() and self._ncols <= options.max_cols():
            return self.str()
        if self.is_sparse():
            s = 'sparse'
        else:
            s = 'dense'
        return "{} x {} {} matrix over {}".format(self._nrows, self._ncols, s, self.base_ring())

    def __str__(self):
        r"""
        Return a string representation of this matrix. Unlike
        ``__repr__`` (used by interactive sessions), this always prints
        the matrix entries.

        EXAMPLES::

            sage: A = zero_matrix(ZZ, 20)
            sage: A
            20 x 20 dense matrix over Integer Ring (use the '.str()' method to see the entries)
            sage: print(A)
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        """
        return self.str()

    def str(self, rep_mapping=None, zero=None, plus_one=None, minus_one=None,
            *, unicode=False, shape=None, character_art=False,
            left_border=None, right_border=None,
            top_border=None, bottom_border=None):
        r"""
        Return a nice string representation of the matrix.

        INPUT:

        - ``rep_mapping`` -- dictionary or callable used to override
          the usual representation of elements

          If ``rep_mapping`` is a dictionary then keys should be
          elements of the base ring and values the desired string
          representation.  Values sent in via the other keyword
          arguments will override values in the dictionary.
          Use of a dictionary can potentially take a very long time
          due to the need to hash entries of the matrix.  Matrices
          with entries from ``QQbar`` are one example.

          If ``rep_mapping`` is callable then it will be called with
          elements of the matrix and must return a string.  Simply
          call :func:`repr` on elements which should have the default
          representation.

        - ``zero`` -- string (default: ``None``); if not ``None`` use
          the value of ``zero`` as the representation of the zero
          element.

        - ``plus_one`` -- string (default: ``None``); if not ``None``
          use the value of ``plus_one`` as the representation of the
          one element.

        - ``minus_one`` -- string (default: ``None``); if not ``None``
          use the value of ``minus_one`` as the representation of the
          negative of the one element.

        - ``unicode`` -- boolean (default: ``False``);
          whether to use Unicode symbols instead of ASCII symbols
          for brackets and subdivision lines

        - ``shape`` -- one of ``'square'`` or ``'round'`` (default: ``None``).
          Switches between round and square brackets.
          The default depends on the setting of the ``unicode`` keyword
          argument. For Unicode symbols, the default is round brackets
          in accordance with the TeX rendering,
          while the ASCII rendering defaults to square brackets.

        - ``character_art`` -- boolean (default: ``False``); if ``True``, the
          result will be of type :class:`~sage.typeset.ascii_art.AsciiArt` or
          :class:`~sage.typeset.unicode_art.UnicodeArt` which support line
          breaking of wide matrices that exceed the window width

        - ``left_border``, ``right_border`` -- sequence (default: ``None``);
          if not ``None``, call :func:`str` on the elements and use the
          results as labels for the rows of the matrix. The labels appear
          outside of the parentheses.

        - ``top_border``, ``bottom_border`` -- sequence (default: ``None``);
          if not ``None``, call :func:`str` on the elements and use the
          results as labels for the columns of the matrix. The labels appear
          outside of the parentheses.

        EXAMPLES::

            sage: R = PolynomialRing(QQ,6,'z')
            sage: a = matrix(2,3, R.gens())
            sage: a.__repr__()
            '[z0 z1 z2]\n[z3 z4 z5]'

            sage: M = matrix([[1,0],[2,-1]])
            sage: M.str()
            '[ 1  0]\n[ 2 -1]'
            sage: M.str(plus_one='+',minus_one='-',zero='.')
            '[+ .]\n[2 -]'
            sage: M.str({1:"not this one",2:"II"},minus_one='*',plus_one='I')
            '[ I  0]\n[II  *]'

            sage: def print_entry(x):
            ....:   if x>0:
            ....:       return '+'
            ....:   elif x<0:
            ....:       return '-'
            ....:   else: return '.'
            ...
            sage: M.str(print_entry)
            '[+ .]\n[+ -]'
            sage: M.str(repr)
            '[ 1  0]\n[ 2 -1]'

            sage: M = matrix([[1,2,3],[4,5,6],[7,8,9]])
            sage: M.subdivide(None, 2)
            sage: print(M.str(unicode=True))
            ⎛1 2│3⎞
            ⎜4 5│6⎟
            ⎝7 8│9⎠
            sage: M.subdivide([0,1,1,3], [0,2,3,3])
            sage: print(M.str(unicode=True, shape='square'))
            ⎡┼───┼─┼┼⎤
            ⎢│1 2│3││⎥
            ⎢┼───┼─┼┼⎥
            ⎢┼───┼─┼┼⎥
            ⎢│4 5│6││⎥
            ⎢│7 8│9││⎥
            ⎣┼───┼─┼┼⎦

        If ``character_art`` is set, the lines of large matrices are wrapped in
        a readable way::

            sage: set_random_seed(0)
            sage: matrix.random(RDF, 3, 5).str(unicode=True, character_art=True)
            ⎛ -0.27440062056807446    0.5031965950979831 -0.001975438590219314
            ⎜ -0.05461130074681608 -0.033673314214051286   -0.9401270875197381
            ⎝  0.19906256610645512    0.3242250183948632    0.6026443545751128
            <BLANKLINE>
               -0.9467802263760512    0.5056889961514748⎞
              -0.35104242112828943    0.5084492941557279⎟
               -0.9541798283979341   -0.8948790563276592⎠

        The number of floating point digits to display is controlled by
        :obj:`matrix.options.precision <.constructor.options>` and can also be
        set by the `IPython magic
        <https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-precision>`_
        ``%precision``. This does not affect the internal precision of the
        represented data, but only the textual display of matrices::

            sage: matrix.options.precision = 4
            sage: A = matrix(RR, [[1/3, 200/3], [-3, 1e6]]); A
            [  0.3333    66.67]
            [  -3.000 1.000E+6]
            sage: unicode_art(A)
            ⎛  0.3333    66.67⎞
            ⎝  -3.000 1.000E+6⎠
            sage: matrix.options.precision = None
            sage: A
            [ 0.333333333333333   66.6666666666667]
            [ -3.00000000000000 1.00000000000000e6]

        Matrices with borders::

            sage: M = matrix([[1,2,3], [4,5,6], [7,8,9]])
            sage: M.subdivide(None, 2)
            sage: print(M.str(unicode=True,
            ....:             top_border=['ab', 'cde', 'f'],
            ....:             bottom_border=['*', '', ''],
            ....:             left_border=[1, 10, 100],
            ....:             right_border=['', ' <', '']))
                 ab cde   f
              1⎛  1   2│  3⎞
             10⎜  4   5│  6⎟ <
            100⎝  7   8│  9⎠
                  *

        TESTS:

        Prior to :issue:`11544` this could take a full minute to run (2011). ::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQ, 4, 4, [1, 2, -2, 2, 1, 0, -1, -1, 0, -1, 1, 1, -1, 2, 1/2, 0])
            sage: e = A.eigenvalues()[3]
            sage: K = (A - e).kernel()
            sage: P = K.basis_matrix()
            sage: P.str()
            '[              1.000000000000000? + 0.?e-17*I  -2.116651487479748? + 0.0255565807096352?*I -0.2585224251020429? + 0.2886023409047535?*I  -0.4847545623533090? - 1.871890760086142?*I]'

        Use single-row delimiters where appropriate::

            sage: print(matrix([[1]]).str(unicode=True))
            (1)
            sage: print(matrix([[],[]]).str(unicode=True))
            ()
            sage: M = matrix([[1]])
            sage: M.subdivide([0,1], [])
            sage: print(M.str(unicode=True))
            ⎛─⎞
            ⎜1⎟
            ⎝─⎠

        Check that exact number types are not affected by the precision
        option::

            sage: matrix.options.precision = 4
            sage: matrix(ZZ, [[10^10]])
            [10000000000]
            sage: matrix(QQ, [[2/3, 10^6]])
            [    2/3 1000000]
            sage: R.<x,y> = QQ[[]]
            sage: matrix(R, [[2/3 - 10^6 * x^3 + 3 * y + O(x, y)^4]])
            [2/3 + 3*y - 1000000*x^3 + O(x, y)^4]
            sage: matrix.options._reset()

        Edge cases of matrices with borders::

            sage: print(matrix(ZZ, 0, 0).str(
            ....:     top_border=[], bottom_border=[], left_border=[], right_border=[]))
            []
            sage: print(matrix(ZZ, 0, 4).str(
            ....:     unicode=True,
            ....:     top_border='abcd', bottom_border=range(4)))
            ()
            sage: print(matrix(ZZ, 1, 4).str(
            ....:     unicode=True,
            ....:     top_border='abcd', bottom_border=range(4)))
             a b c d
            (0 0 0 0)
             0 1 2 3
            sage: print(matrix(ZZ, 2, 4).str(
            ....:     unicode=True,
            ....:     top_border='abcd', bottom_border=range(4), left_border='uv'))
              a b c d
            u⎛0 0 0 0⎞
            v⎝0 0 0 0⎠
              0 1 2 3
            sage: print(matrix(ZZ, 2, 0).str(
            ....:     top_border='', left_border='uv', right_border=['*', '']))
              []
        """
        cdef Py_ssize_t nr, nc, r, c
        nr = self._nrows
        nc = self._ncols

        # symbols is a string with 11 elements:
        # - top left bracket         (tlb)
        # - middle left bracket      (mlb)
        # - bottom left bracket      (blb)
        # - single-row left bracket  (slb)
        # - top right bracket        (trb)
        # - middle right bracket     (mrb)
        # - bottom right bracket     (brb)
        # - single-row right bracket (srb)
        # - vertical line            (vl)
        # - horizontal line          (hl)
        # - crossing lines           (cl)
        if shape is None:
            shape = "round" if unicode else "square"
        if unicode:
            import unicodedata
            hl = unicodedata.lookup('BOX DRAWINGS LIGHT HORIZONTAL')
            vl = unicodedata.lookup('BOX DRAWINGS LIGHT VERTICAL')
            cl = unicodedata.lookup('BOX DRAWINGS LIGHT VERTICAL AND HORIZONTAL')
        else:
            hl = '-'        # - horizontal line
            vl = '|'        # - vertical line
            cl = '+'        # - crossing lines
        if shape == "square":
            if unicode:
                from sage.typeset.symbols import (
                    unicode_left_square_bracket as left,
                    unicode_right_square_bracket as right
                )
            else:
                from sage.typeset.symbols import (
                    ascii_left_square_bracket as left,
                    ascii_right_square_bracket as right
                )
        elif shape == "round":
            if unicode:
                from sage.typeset.symbols import (
                    unicode_left_parenthesis as left,
                    unicode_right_parenthesis as right
                )
            else:
                from sage.typeset.symbols import (
                    ascii_left_parenthesis as left,
                    ascii_right_parenthesis as right
                )
        else:
            raise ValueError("No such shape")
        tlb = left.top              # - top left bracket
        mlb = left.extension        # - extension piece left bracket
        blb = left.bottom           # - bottom left bracket
        slb = left.character        # - single-row left bracket
        trb = right.top             # - top right bracket
        mrb = right.extension       # - extension piece right bracket
        brb = right.bottom          # - bottom right bracket
        srb = right.character       # - single-row right bracket

        if character_art:
            if unicode:
                from sage.typeset.unicode_art import UnicodeArt as CharacterArt
            else:
                from sage.typeset.ascii_art import AsciiArt as CharacterArt

        if nr == 0 or nc == 0:
            result = slb + srb
            return CharacterArt([result]) if character_art else result

        row_divs, col_divs = self.subdivisions()
        row_div_counts = [0] * (nr + 1)
        for r in row_divs:
            row_div_counts[r] += 1
        col_div_counts = [0] * (nc + 1)
        for c in col_divs:
            col_div_counts[c] += 1

        # Set the mapping based on keyword arguments
        if rep_mapping is None:
            rep_mapping = {}
        if isinstance(rep_mapping, dict):
            if zero is not None:
                rep_mapping[self.base_ring().zero()] = zero
            if plus_one is not None:
                rep_mapping[self.base_ring().one()] = plus_one
            if minus_one is not None:
                rep_mapping[-self.base_ring().one()] = minus_one

        entries = self.list()

        # only use floating point formatting for inexact types that have
        # custom implementation of __format__
        from sage.matrix.constructor import options
        prec = options.precision()
        if prec is None or callable(rep_mapping) or not entries \
                or type(entries[0]).__format__ is Element.__format__ \
                or self._base_ring.is_exact():
            fmt_numeric = None
        else:
            fmt_numeric = options.format_numeric()

        # compute column widths
        S = []
        if top_border is not None:
            for x in top_border:
                S.append(str(x))
            top_count = 1
        else:
            top_count = 0
        for x in entries:
            # Override the usual representations with those specified
            if callable(rep_mapping):
                rep = rep_mapping(x)
            # avoid hashing entries, especially algebraic numbers
            elif rep_mapping and x in rep_mapping:
                rep = rep_mapping.get(x)
            elif fmt_numeric is not None:
                rep = fmt_numeric.format(x, prec=prec)
            else:
                rep = repr(x)
            S.append(rep)
        if bottom_border is not None:
            for x in bottom_border:
                S.append(str(x))
            bottom_count = 1
        else:
            bottom_count = 0

        width = max(map(len, S))
        left = []
        rows = []
        right = []

        hline = cl.join(hl * ((width + 1)*(b - a) - 1)
                        for a,b in zip([0] + col_divs, col_divs + [nc]))

        # compute rows
        for r in range(-top_count, nr + bottom_count):
            if 0 <= r < nr:
                n = row_divs.count(r)
                if n:
                    left.extend([""] * n)
                    rows.extend([hline] * n)
                    right.extend([""] * n)
            if left_border is not None and 0 <= r < nr:
                left.append(str(left_border[r]))
            else:
                left.append("")
            s = ""
            for c from 0 <= c < nc:
                if col_div_counts[c]:
                    if 0 <= r < nr:
                        sep = vl * col_div_counts[c]
                    else:
                        sep = " " * col_div_counts[c]
                elif c == 0:
                    sep = ""
                else:
                    sep = " "
                entry = S[(r + top_count) * nc + c]
                entry = " " * (width - len(entry)) + entry
                s = s + sep + entry
            else:
                if 0 <= r < nr:
                    s = s + vl * col_div_counts[nc]
                else:
                    s = s + " " * col_div_counts[nc]
            rows.append(s)
            if right_border is not None and 0 <= r < nr:
                right.append(str(right_border[r]))
            else:
                right.append("")
        else:
            if nr == nr + bottom_count:
                n = row_divs.count(nr)
                if n:
                    left.extend([""] * n)
                    rows.extend([hline] * n)
                    right.extend([""] * n)

        # left and right brackets
        for i in range(top_count):
            rows[i] = " "*len(slb) + rows[i] + " "*len(srb)
        if len(rows) == top_count + 1 + bottom_count:
            rows[top_count] = slb + rows[top_count] + srb
        else:
            rows[top_count] = tlb + rows[top_count] + trb
            for i in range(top_count + 1, len(rows) - bottom_count - 1):
                rows[i] = mlb + rows[i] + mrb
            rows[-1 - bottom_count] = blb + rows[-1 - bottom_count] + brb
        for i in range(bottom_count):
            rows[-1 - i] = " "*len(slb) + rows[-1 - i] + " "*len(srb)

        # left and right border
        left_width = max(len(s) for s in left)
        right_width = max(len(s) for s in right)
        for i in range(len(rows)):
            rows[i] = left[i].rjust(left_width) + rows[i] + right[i].rjust(right_width)

        if character_art:
            breakpoints = []
            idx = len(tlb) + (col_div_counts[0] if nc > 0 else 0) + width
            for c from 1 <= c < nc:
                breakpoints.append(idx)
                len_sep = max(col_div_counts[c], 1)
                idx += len_sep + width
            return CharacterArt(rows, breakpoints=breakpoints)
        else:
            return "\n".join(rows)

    def _ascii_art_(self):
        """
        Return an ASCII art representation of this matrix.

        EXAMPLES::

            sage: set_random_seed(0)
            sage: ascii_art(matrix.random(RDF, 3, 5))  # indirect doctest
            [ -0.27440062056807446    0.5031965950979831 -0.001975438590219314
            [ -0.05461130074681608 -0.033673314214051286   -0.9401270875197381
            [  0.19906256610645512    0.3242250183948632    0.6026443545751128
            <BLANKLINE>
               -0.9467802263760512    0.5056889961514748]
              -0.35104242112828943    0.5084492941557279]
               -0.9541798283979341   -0.8948790563276592]
        """
        from sage.matrix.constructor import options
        if self._nrows <= options.max_rows() and self._ncols <= options.max_cols():
            return self.str(character_art=True)
        else:
            from sage.typeset.ascii_art import AsciiArt
            return AsciiArt(repr(self).splitlines())

    def _unicode_art_(self):
        """
        Return a unicode art representation of this matrix.

        EXAMPLES::

            sage: A = matrix([[1,2], [3,4], [5,6]])
            sage: A._unicode_art_()
            ⎛1 2⎞
            ⎜3 4⎟
            ⎝5 6⎠
            sage: unicode_art(A)    # indirect doctest
            ⎛1 2⎞
            ⎜3 4⎟
            ⎝5 6⎠

        If the matrix is too big, don't print all of the elements::

            sage: A = random_matrix(ZZ, 100)
            sage: unicode_art(A)
            100 x 100 dense matrix over Integer Ring
        """
        from sage.matrix.constructor import options
        if self._nrows <= options.max_rows() and self._ncols <= options.max_cols():
            return self.str(unicode=True, character_art=True)
        else:
            from sage.typeset.unicode_art import UnicodeArt
            return UnicodeArt(repr(self).splitlines())

    def _latex_(self):
        r"""
        Return latex representation of this matrix.  The matrix is
        enclosed in parentheses by default, but the delimiters can be
        changed using the command
        ``latex.matrix_delimiters(...)``.

        EXAMPLES::

            sage: R = PolynomialRing(QQ,4,'z')
            sage: a = matrix(2,2, R.gens())
            sage: b = a*a
            sage: latex(b) # indirect doctest
            \left(\begin{array}{rr}
            z_{0}^{2} + z_{1} z_{2} & z_{0} z_{1} + z_{1} z_{3} \\
            z_{0} z_{2} + z_{2} z_{3} & z_{1} z_{2} + z_{3}^{2}
            \end{array}\right)

        Latex representation for block matrices::

            sage: B = matrix(3,4)
            sage: B.subdivide([2,2], [3])
            sage: latex(B)
            \left(\begin{array}{rrr|r}
            0 & 0 & 0 & 0 \\
            0 & 0 & 0 & 0 \\
            \hline\hline
            0 & 0 & 0 & 0
            \end{array}\right)
        """
        latex = sage.misc.latex.latex
        matrix_delimiters = latex.matrix_delimiters()
        align = latex.matrix_column_alignment()
        cdef Py_ssize_t nr, nc, r, c
        nr = self._nrows
        nc = self._ncols
        if nr == 0 or nc == 0:
            return matrix_delimiters[0] + matrix_delimiters[1]

        S = self.list()
        rows = []

        row_divs, col_divs = self.subdivisions()

        # construct one large array, using \hline and vertical
        # bars | in the array descriptor to indicate subdivisions.
        for r from 0 <= r < nr:
            if r in row_divs:
                s = "\\hline"*row_divs.count(r) + "\n"
            else:
                s = ""
            for c from 0 <= c < nc:
                if c == nc-1:
                    sep=""
                else:
                    sep=" & "
                entry = latex(S[r*nc+c])
                s = s + entry + sep
            rows.append(s)

        # Put brackets around in a single string
        tmp = []
        for row in rows:
            tmp.append(str(row))
        s = " \\\\\n".join(tmp)

        tmp = [align*(b-a) for a,b in zip([0] + col_divs, col_divs + [nc])]
        format = '|'.join(tmp)

        return "\\left" + matrix_delimiters[0] + "\\begin{array}{%s}\n" % format + s + "\n\\end{array}\\right" + matrix_delimiters[1]

    ###################################################
    ## Basic Properties
    ###################################################

    def ncols(self):
        """
        Return the number of columns of this matrix.

        EXAMPLES::

            sage: M = MatrixSpace(QQ, 2, 3)
            sage: A = M([1,2,3, 4,5,6])
            sage: A
            [1 2 3]
            [4 5 6]
            sage: A.ncols()
            3
            sage: A.nrows()
            2

        AUTHORS:

        - Naqi Jaffery (2006-01-24): examples
        """
        return self._ncols

    def nrows(self):
        r"""
        Return the number of rows of this matrix.

        EXAMPLES::

            sage: M = MatrixSpace(QQ,6,7)
            sage: A = M([1,2,3,4,5,6,7, 22,3/4,34,11,7,5,3, 99,65,1/2,2/3,3/5,4/5,5/6, 9,8/9, 9/8,7/6,6/7,76,4, 0,9,8,7,6,5,4, 123,99,91,28,6,1024,1])
            sage: A
            [   1    2    3    4    5    6    7]
            [  22  3/4   34   11    7    5    3]
            [  99   65  1/2  2/3  3/5  4/5  5/6]
            [   9  8/9  9/8  7/6  6/7   76    4]
            [   0    9    8    7    6    5    4]
            [ 123   99   91   28    6 1024    1]
            sage: A.ncols()
            7
            sage: A.nrows()
            6

        AUTHORS:

        - Naqi Jaffery (2006-01-24): examples
        """
        return self._nrows

    def dimensions(self):
        r"""
        Return the dimensions of this matrix as the tuple (nrows, ncols).

        EXAMPLES::

            sage: M = matrix([[1,2,3],[4,5,6]])
            sage: N = M.transpose()
            sage: M.dimensions()
            (2, 3)
            sage: N.dimensions()
            (3, 2)

        AUTHORS:

        - Benjamin Lundell (2012-02-09): examples
        """
        return (self._nrows,self._ncols)

    ###################################################
    # Functions
    ###################################################

    def act_on_polynomial(self, f):
        r"""
        Return the polynomial ``f(self*x)``.

        INPUT:

        - ``self`` -- an nxn matrix

        - ``f`` -- a polynomial in n variables x=(x1,...,xn)

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: x, y = R.gens()
            sage: f = x**2 - y**2
            sage: M = MatrixSpace(QQ, 2)
            sage: A = M([1,2,3,4])
            sage: A.act_on_polynomial(f)
            -8*x^2 - 20*x*y - 12*y^2
        """
        cdef Py_ssize_t i, j, n

        if self._nrows != self._ncols:
            raise ArithmeticError("self must be a square matrix")

        vars = f.parent().gens()
        n = len(self.rows())
        ans = []
        for i from 0 <= i < n:
            tmp = []
            for j from 0 <= j < n:
                tmp.append(self.get_unsafe(i, j)*vars[j])
            ans.append( sum(tmp) )
        return f(tuple(ans))

    def __call__(self, *args, **kwargs):
        """
        Calling a matrix returns the result of calling each component.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: f(x,y) = x^2 + y
            sage: m = matrix([[f, f*f], [f^3, f^4]]); m
            [    (x, y) |--> x^2 + y (x, y) |--> (x^2 + y)^2]
            [(x, y) |--> (x^2 + y)^3 (x, y) |--> (x^2 + y)^4]
            sage: m(1, 2)
            [ 3  9]
            [27 81]
            sage: m(y=2, x=1)
            [ 3  9]
            [27 81]
            sage: m(2, 1)
            [  5  25]
            [125 625]
        """
        from sage.matrix.constructor import matrix
        return matrix(self.nrows(), self.ncols(), [e(*args, **kwargs) for e in self.list()])

    ###################################################
    # Arithmetic
    ###################################################
    def commutator(self, other):
        r"""
        Return the commutator self\*other - other\*self.

        EXAMPLES::

            sage: A = Matrix(ZZ, 2, 2, range(4))
            sage: B = Matrix(ZZ, 2, 2, [0, 1, 0, 0])
            sage: A.commutator(B)
            [-2 -3]
            [ 0  2]
            sage: A.commutator(B) == -B.commutator(A)
            True
        """
        return self*other - other*self

    def anticommutator(self, other):
        r"""
        Return the anticommutator ``self`` and ``other``.

        The *anticommutator* of two `n \times n` matrices `A` and `B`
        is defined as `\{A, B\} := AB + BA` (sometimes this is written as
        `[A, B]_+`).

        EXAMPLES::

            sage: A = Matrix(ZZ, 2, 2, range(4))
            sage: B = Matrix(ZZ, 2, 2, [0, 1, 0, 0])
            sage: A.anticommutator(B)
            [2 3]
            [0 2]
            sage: A.anticommutator(B) == B.anticommutator(A)
            True
            sage: A.commutator(B) + B.anticommutator(A) == 2*A*B
            True
        """
        return self*other + other*self

    ###################################################
    # Row and column operations
    # The _c versions do no bounds checking.
    # The with_ versions do not change the input matrix.
    # Some of the functions assume that input values
    # have parent that is self._base_ring.
    # AUTHORS:
    #     -- Karl-Dieter Crisman (June 2008):
    # Improved examples and error messages for methods which could
    # involve multiplication outside base ring, including
    # with_ versions of these methods for this situation
    ###################################################
    cdef check_row_bounds(self, Py_ssize_t r1, Py_ssize_t r2):
        if r1 < 0 or r1 >= self._nrows or r2 < 0 or r2 >= self._nrows:
            raise IndexError("matrix row index out of range")

    cdef check_row_bounds_and_mutability(self, Py_ssize_t r1, Py_ssize_t r2):
        if self._is_immutable:
            raise ValueError("Matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).")
        else:
            self._cache = None
        if r1 < 0 or r1 >= self._nrows or r2 < 0 or r2 >= self._nrows:
            raise IndexError("matrix row index out of range")

    cdef check_column_bounds(self, Py_ssize_t c1, Py_ssize_t c2):
        if c1 < 0 or c1 >= self._ncols or c2 < 0 or c2 >= self._ncols:
            raise IndexError("matrix column index out of range")

    cdef check_column_bounds_and_mutability(self, Py_ssize_t c1, Py_ssize_t c2):
        if self._is_immutable:
            raise ValueError("Matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).")
        else:
            self._cache = None
        if c1 < 0 or c1 >= self._ncols or c2 < 0 or c2 >= self._ncols:
            raise IndexError("matrix column index out of range")

    def swap_columns(self, Py_ssize_t c1, Py_ssize_t c2):
        """
        Swap columns c1 and c2 of ``self``.

        EXAMPLES: We create a rational matrix::

            sage: M = MatrixSpace(QQ,3,3)
            sage: A = M([1,9,-7,4/5,4,3,6,4,3])
            sage: A
            [  1   9  -7]
            [4/5   4   3]
            [  6   4   3]

        Since the first column is numbered zero, this swaps the second and
        third columns::

            sage: A.swap_columns(1,2); A
            [  1  -7   9]
            [4/5   3   4]
            [  6   3   4]
        """
        self.check_column_bounds_and_mutability(c1, c2)
        if c1 != c2:
            self.swap_columns_c(c1, c2)

    def with_swapped_columns(self, c1, c2):
        r"""
        Swap columns ``c1`` and ``c2`` of ``self`` and return a new matrix.

        INPUT:

        - ``c1``, ``c2`` -- integers specifying columns of ``self`` to interchange

        OUTPUT:

        A new matrix, identical to ``self`` except that columns ``c1`` and ``c2``
        are swapped.

        EXAMPLES:

        Remember that columns are numbered starting from zero. ::

            sage: A = matrix(QQ, 4, range(20))
            sage: A.with_swapped_columns(1, 2)
            [ 0  2  1  3  4]
            [ 5  7  6  8  9]
            [10 12 11 13 14]
            [15 17 16 18 19]

        Trying to swap a column with itself will succeed, but still return
        a new matrix. ::

            sage: A = matrix(QQ, 4, range(20))
            sage: B = A.with_swapped_columns(2, 2)
            sage: A == B
            True
            sage: A is B
            False

        The column specifications are checked. ::

            sage: A = matrix(4, range(20))
            sage: A.with_swapped_columns(-1, 2)
            Traceback (most recent call last):
            ...
            IndexError: matrix column index out of range

            sage: A.with_swapped_columns(2, 5)
            Traceback (most recent call last):
            ...
            IndexError: matrix column index out of range
        """
        cdef Matrix temp
        self.check_column_bounds_and_mutability(c1,c2)
        temp = self.__copy__()
        if c1 != c2:
            temp.swap_columns_c(c1,c2)
        return temp

    def permute_columns(self, permutation):
        r"""
        Permute the columns of ``self`` by applying the permutation
        group element ``permutation``.

        As permutation group elements act on integers `\{1,\dots,n\}`,
        columns are considered numbered from 1 for this operation.

        INPUT:

        - ``permutation`` -- a ``PermutationGroupElement``

        EXAMPLES: We create a matrix::

            sage: M = matrix(ZZ, [[1,0,0,0,0],[0,2,0,0,0],[0,0,3,0,0],[0,0,0,4,0],[0,0,0,0,5]])
            sage: M
            [1 0 0 0 0]
            [0 2 0 0 0]
            [0 0 3 0 0]
            [0 0 0 4 0]
            [0 0 0 0 5]

        Next of all, create a permutation group element and act
        on ``M`` with it::

            sage: # needs sage.groups
            sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
            sage: sigma, tau = G.gens()
            sage: sigma
            (1,2,3)(4,5)
            sage: M.permute_columns(sigma)
            sage: M
            [0 0 1 0 0]
            [2 0 0 0 0]
            [0 3 0 0 0]
            [0 0 0 0 4]
            [0 0 0 5 0]
        """
        self.check_mutability()
        for cycle in permutation.cycle_tuples():
            cycle = [elt-1 for elt in reversed(cycle)]
            for elt in cycle:
                self.check_column_bounds(cycle[0], elt)
                if cycle[0] != elt:
                    self.swap_columns_c(cycle[0], elt)

    def with_permuted_columns(self, permutation):
        r"""
        Return the matrix obtained from permuting the columns
        of ``self`` by applying the permutation group element
        ``permutation``.

        As permutation group elements act on integers `\{1,\dots,n\}`,
        columns are considered numbered from 1 for this operation.

        INPUT:

        - ``permutation`` -- a ``PermutationGroupElement``

        OUTPUT: a matrix

        EXAMPLES: We create some matrix::

            sage: M = matrix(ZZ, [[1,0,0,0,0],[0,2,0,0,0],[0,0,3,0,0],[0,0,0,4,0],[0,0,0,0,5]])
            sage: M
            [1 0 0 0 0]
            [0 2 0 0 0]
            [0 0 3 0 0]
            [0 0 0 4 0]
            [0 0 0 0 5]

        Next of all, create a permutation group element and
        act on ``M``::

            sage: # needs sage.groups
            sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
            sage: sigma, tau = G.gens()
            sage: sigma
            (1,2,3)(4,5)
            sage: M.with_permuted_columns(sigma)
            [0 0 1 0 0]
            [2 0 0 0 0]
            [0 3 0 0 0]
            [0 0 0 0 4]
            [0 0 0 5 0]
        """
        cdef Matrix temp
        temp = self.__copy__()
        for cycle in permutation.cycle_tuples():
            cycle = [(elt - 1) for elt in reversed(cycle)]
            for elt in cycle:
                self.check_column_bounds(cycle[0], elt)
                if cycle[0] != elt:
                    temp.swap_columns_c(cycle[0], elt)
        return temp

    cdef swap_columns_c(self, Py_ssize_t c1, Py_ssize_t c2):
        cdef Py_ssize_t r
        for r from 0 <= r < self._nrows:
            a = self.get_unsafe(r, c2)
            self.set_unsafe(r, c2, self.get_unsafe(r,c1))
            self.set_unsafe(r, c1, a)

    def swap_rows(self, r1, r2):
        """
        Swap rows r1 and r2 of ``self``.

        EXAMPLES: We create a rational matrix::

            sage: M = MatrixSpace(QQ, 3, 3)
            sage: A = M([1,9,-7, 4/5,4,3, 6,4,3])
            sage: A
            [  1   9  -7]
            [4/5   4   3]
            [  6   4   3]

        Since the first row is numbered zero, this swaps the first and
        third rows::

            sage: A.swap_rows(0, 2); A
            [  6   4   3]
            [4/5   4   3]
            [  1   9  -7]
        """
        self.check_row_bounds_and_mutability(r1, r2)
        if r1 != r2:
            self.swap_rows_c(r1, r2)

    def with_swapped_rows(self, r1, r2):
        r"""
        Swap rows ``r1`` and ``r2`` of ``self`` and return a new matrix.

        INPUT:

        - ``r1``, ``r2`` -- integers specifying rows of ``self`` to interchange

        OUTPUT:

        A new matrix, identical to ``self`` except that rows ``r1`` and ``r2``
        are swapped.

        EXAMPLES:

        Remember that rows are numbered starting from zero. ::

            sage: A = matrix(QQ, 4, range(20))
            sage: A.with_swapped_rows(1, 2)
            [ 0  1  2  3  4]
            [10 11 12 13 14]
            [ 5  6  7  8  9]
            [15 16 17 18 19]

        Trying to swap a row with itself will succeed, but still return
        a new matrix. ::

            sage: A = matrix(QQ, 4, range(20))
            sage: B = A.with_swapped_rows(2, 2)
            sage: A == B
            True
            sage: A is B
            False

        The row specifications are checked. ::

            sage: A = matrix(4, range(20))
            sage: A.with_swapped_rows(-1, 2)
            Traceback (most recent call last):
            ...
            IndexError: matrix row index out of range

            sage: A.with_swapped_rows(2, 5)
            Traceback (most recent call last):
            ...
            IndexError: matrix row index out of range
        """
        cdef Matrix temp
        self.check_row_bounds_and_mutability(r1,r2)
        temp = self.__copy__()
        if r1 != r2:
            temp.swap_rows_c(r1,r2)
        return temp

    def permute_rows(self, permutation):
        r"""
        Permute the rows of ``self`` by applying the permutation
        group element ``permutation``.

        As permutation group elements act on integers `\{1,\dots,n\}`,
        rows are considered numbered from 1 for this operation.

        INPUT:

        - ``permutation`` -- a ``PermutationGroupElement``

        EXAMPLES: We create a matrix::

            sage: M = matrix(ZZ, [[1,0,0,0,0],[0,2,0,0,0],[0,0,3,0,0],[0,0,0,4,0],[0,0,0,0,5]])
            sage: M
            [1 0 0 0 0]
            [0 2 0 0 0]
            [0 0 3 0 0]
            [0 0 0 4 0]
            [0 0 0 0 5]

        Next of all, create a permutation group element and act on ``M``::

            sage: # needs sage.groups
            sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
            sage: sigma, tau = G.gens()
            sage: sigma
            (1,2,3)(4,5)
            sage: M.permute_rows(sigma)
            sage: M
            [0 2 0 0 0]
            [0 0 3 0 0]
            [1 0 0 0 0]
            [0 0 0 0 5]
            [0 0 0 4 0]
        """
        self.check_mutability()
        for cycle in permutation.cycle_tuples():
            cycle = [elt - 1 for elt in reversed(cycle)]
            for elt in cycle:
                self.check_row_bounds(cycle[0], elt)
                if cycle[0] != elt:
                    self.swap_rows_c(cycle[0], elt)

    def with_permuted_rows(self, permutation):
        r"""
        Return the matrix obtained from permuting the rows
        of ``self`` by applying the permutation group element
        ``permutation``.

        As permutation group elements act on integers `\{1,\dots,n\}`,
        rows are considered numbered from 1 for this operation.

        INPUT:

        - ``permutation`` -- a ``PermutationGroupElement``

        OUTPUT: a matrix

        EXAMPLES: We create a matrix::

            sage: M = matrix(ZZ, [[1,0,0,0,0],[0,2,0,0,0],[0,0,3,0,0],[0,0,0,4,0],[0,0,0,0,5]])
            sage: M
            [1 0 0 0 0]
            [0 2 0 0 0]
            [0 0 3 0 0]
            [0 0 0 4 0]
            [0 0 0 0 5]

        Next of all, create a permutation group element and act on ``M``::

            sage: # needs sage.groups
            sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
            sage: sigma, tau = G.gens()
            sage: sigma
            (1,2,3)(4,5)
            sage: M.with_permuted_rows(sigma)
            [0 2 0 0 0]
            [0 0 3 0 0]
            [1 0 0 0 0]
            [0 0 0 0 5]
            [0 0 0 4 0]
        """
        cdef Matrix temp
        temp = self.__copy__()
        for cycle in permutation.cycle_tuples():
            cycle = [elt - 1 for elt in reversed(cycle)]
            for elt in cycle:
                self.check_row_bounds(cycle[0], elt)
                if cycle[0] != elt:
                    temp.swap_rows_c(cycle[0], elt)
        return temp

    cdef swap_rows_c(self, Py_ssize_t r1, Py_ssize_t r2):
        cdef Py_ssize_t c
        for c from 0 <= c < self._ncols:
            a = self.get_unsafe(r2, c)
            self.set_unsafe(r2, c, self.get_unsafe(r1, c))
            self.set_unsafe(r1, c, a)

    def permute_rows_and_columns(self, row_permutation, column_permutation):
        r"""
        Permute the rows and columns of ``self`` by applying the permutation
        group elements ``row_permutation`` and ``column_permutation``
        respectively.

        As permutation group elements act on integers `\{1,\dots,n\}`,
        rows and columns are considered numbered from 1 for this operation.

        INPUT:

        - ``row_permutation`` -- a ``PermutationGroupElement``
        - ``column_permutation`` -- a ``PermutationGroupElement``

        OUTPUT: a matrix

        EXAMPLES: We create a matrix::

            sage: M = matrix(ZZ, [[1,0,0,0,0],[0,2,0,0,0],[0,0,3,0,0],[0,0,0,4,0],[0,0,0,0,5]])
            sage: M
            [1 0 0 0 0]
            [0 2 0 0 0]
            [0 0 3 0 0]
            [0 0 0 4 0]
            [0 0 0 0 5]

        Next of all, create a permutation group element and act on ``M``::

            sage: # needs sage.groups
            sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
            sage: sigma, tau = G.gens()
            sage: sigma
            (1,2,3)(4,5)
            sage: M.permute_rows_and_columns(sigma,tau)
            sage: M
            [2 0 0 0 0]
            [0 3 0 0 0]
            [0 0 0 0 1]
            [0 0 0 5 0]
            [0 0 4 0 0]
        """
        self.permute_rows(row_permutation)
        self.permute_columns(column_permutation)

    def with_permuted_rows_and_columns(self,row_permutation,column_permutation):
        r"""
        Return the matrix obtained from permuting the rows and
        columns of ``self`` by applying the permutation group
        elements ``row_permutation`` and ``column_permutation``.

        As permutation group elements act on integers `\{1,\dots,n\}`,
        rows and columns are considered numbered from 1 for this operation.

        INPUT:

        - ``row_permutation`` -- a ``PermutationGroupElement``
        - ``column_permutation`` -- a ``PermutationGroupElement``

        OUTPUT: a matrix

        EXAMPLES: We create a matrix::

            sage: M = matrix(ZZ, [[1,0,0,0,0],[0,2,0,0,0],[0,0,3,0,0],[0,0,0,4,0],[0,0,0,0,5]])
            sage: M
            [1 0 0 0 0]
            [0 2 0 0 0]
            [0 0 3 0 0]
            [0 0 0 4 0]
            [0 0 0 0 5]

        Next of all, create a permutation group element and act on ``M``::

            sage: # needs sage.groups
            sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
            sage: sigma, tau = G.gens()
            sage: sigma
            (1,2,3)(4,5)
            sage: M.with_permuted_rows_and_columns(sigma,tau)
            [2 0 0 0 0]
            [0 3 0 0 0]
            [0 0 0 0 1]
            [0 0 0 5 0]
            [0 0 4 0 0]
        """
        return self.with_permuted_rows(row_permutation).with_permuted_columns(column_permutation)

    def add_multiple_of_row(self, Py_ssize_t i, Py_ssize_t j, s, Py_ssize_t start_col=0):
        """
        Add s times row j to row i.

        EXAMPLES: We add -3 times the first row to the second row of an
        integer matrix, remembering to start numbering rows at zero::

            sage: a = matrix(ZZ,2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.add_multiple_of_row(1,0,-3)
            sage: a
            [ 0  1  2]
            [ 3  1 -1]

        To add a rational multiple, we first need to change the base ring::

            sage: a = a.change_ring(QQ)
            sage: a.add_multiple_of_row(1,0,1/3)
            sage: a
            [   0    1    2]
            [   3  4/3 -1/3]

        If not, we get an error message::

            sage: a.add_multiple_of_row(1, 0, SR.I())                                   # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: Multiplying row by Symbolic Ring element cannot be done over
            Rational Field, use change_ring or with_added_multiple_of_row instead.
        """
        self.check_row_bounds_and_mutability(i, j)
        try:
            s = self._coerce_element(s)
            self.add_multiple_of_row_c(i, j, s, start_col)
        except TypeError:
            raise TypeError('Multiplying row by %s element cannot be done over %s, use change_ring or with_added_multiple_of_row instead.' % (s.parent(), self.base_ring()))

    cdef add_multiple_of_row_c(self, Py_ssize_t i, Py_ssize_t j, s, Py_ssize_t start_col):
        cdef Py_ssize_t c
        for c from start_col <= c < self._ncols:
            self.set_unsafe(i, c, self.get_unsafe(i, c) + s*self.get_unsafe(j, c))

    def with_added_multiple_of_row(self, Py_ssize_t i, Py_ssize_t j, s, Py_ssize_t start_col=0):
        """
        Add s times row j to row i, returning new matrix.

        EXAMPLES: We add -3 times the first row to the second row of an
        integer matrix, remembering to start numbering rows at zero::

            sage: a = matrix(ZZ,2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: b = a.with_added_multiple_of_row(1,0,-3); b
            [ 0  1  2]
            [ 3  1 -1]

        The original matrix is unchanged::

            sage: a
            [0 1 2]
            [3 4 5]

        Adding a rational multiple is okay, and reassigning a variable is
        okay::

            sage: a = a.with_added_multiple_of_row(0,1,1/3); a
            [   1  7/3 11/3]
            [   3    4    5]
        """
        cdef Matrix temp
        self.check_row_bounds_and_mutability(i, j)
        try:
            s = self._coerce_element(s)
            temp = self.__copy__()
            temp.add_multiple_of_row_c(i, j, s, start_col)
            return temp
        # If scaling factor cannot be coerced, change the base ring to
        # one acceptable to both the original base ring and the scaling factor.
        except TypeError:
            temp = self.change_ring(Sequence([s,self.base_ring()(0)]).universe())
            s = temp._coerce_element(s)
            temp.add_multiple_of_row_c(i, j, s, start_col)
            return temp

    def add_multiple_of_column(self, Py_ssize_t i, Py_ssize_t j, s, Py_ssize_t start_row=0):
        """
        Add s times column j to column i.

        EXAMPLES: We add -1 times the third column to the second column of
        an integer matrix, remembering to start numbering cols at zero::

            sage: a = matrix(ZZ,2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.add_multiple_of_column(1,2,-1)
            sage: a
            [ 0 -1  2]
            [ 3 -1  5]

        To add a rational multiple, we first need to change the base ring::

            sage: a = a.change_ring(QQ)
            sage: a.add_multiple_of_column(1,0,1/3)
            sage: a
            [ 0 -1  2]
            [ 3  0  5]

        If not, we get an error message::

            sage: a.add_multiple_of_column(1, 0, SR.I())                                # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: Multiplying column by Symbolic Ring element cannot be done over
            Rational Field, use change_ring or with_added_multiple_of_column instead.
        """
        self.check_column_bounds_and_mutability(i, j)
        try:
            s = self._coerce_element(s)
            self.add_multiple_of_column_c(i, j, s, start_row)
        except TypeError:
            raise TypeError('Multiplying column by %s element cannot be done over %s, use change_ring or with_added_multiple_of_column instead.' % (s.parent(), self.base_ring()))

    cdef add_multiple_of_column_c(self, Py_ssize_t i, Py_ssize_t j, s, Py_ssize_t start_row):
        cdef Py_ssize_t r
        for r from start_row <= r < self._nrows:
            self.set_unsafe(r, i, self.get_unsafe(r, i) + s*self.get_unsafe(r, j))

    def with_added_multiple_of_column(self, Py_ssize_t i, Py_ssize_t j, s, Py_ssize_t start_row=0):
        """
        Add s times column j to column i, returning new matrix.

        EXAMPLES: We add -1 times the third column to the second column of
        an integer matrix, remembering to start numbering cols at zero::

            sage: a = matrix(ZZ, 2, 3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: b = a.with_added_multiple_of_column(1, 2, -1); b
            [ 0 -1  2]
            [ 3 -1  5]

        The original matrix is unchanged::

            sage: a
            [0 1 2]
            [3 4 5]

        Adding a rational multiple is okay, and reassigning a variable is
        okay::

            sage: a = a.with_added_multiple_of_column(0, 1, 1/3); a
            [ 1/3    1    2]
            [13/3    4    5]
        """
        cdef Matrix temp
        self.check_column_bounds_and_mutability(i, j)
        try:
            s = self._coerce_element(s)
            temp = self.__copy__()
            temp.add_multiple_of_column_c(i, j, s, start_row)
            return temp
        # If scaling factor cannot be coerced, change the base ring to
        # one acceptable to both the original base ring and the scaling factor.
        except TypeError:
            temp = self.change_ring(Sequence([s,self.base_ring()(0)]).universe())
            s = temp._coerce_element(s)
            temp.add_multiple_of_column_c(i, j, s, start_row)
            return temp

    def rescale_row(self, Py_ssize_t i, s, Py_ssize_t start_col=0):
        """
        Replace `i`-th row of ``self`` by `s` times `i`-th row of ``self``.

        INPUT:

        - ``i`` -- `i`-th row

        - ``s`` -- scalar

        - ``start_col`` -- only rescale entries at this column
          and to the right

        EXAMPLES: We rescale the second row of a matrix over the rational
        numbers::

            sage: a = matrix(QQ, 3, range(6)); a
            [0 1]
            [2 3]
            [4 5]
            sage: a.rescale_row(1, 1/2); a
            [ 0   1]
            [ 1 3/2]
            [ 4   5]

        We rescale the second row of a matrix over a polynomial ring::

            sage: R.<x> = QQ[]
            sage: a = matrix(R, 3, [1,x,x^2,x^3,x^4,x^5]); a
            [  1   x]
            [x^2 x^3]
            [x^4 x^5]
            sage: a.rescale_row(1, 1/2); a
            [      1       x]
            [1/2*x^2 1/2*x^3]
            [    x^4     x^5]

        We try and fail to rescale a matrix over the integers by a
        non-integer::

            sage: a = matrix(ZZ, 2, 3, [0,1,2, 3,4,4]); a
            [0 1 2]
            [3 4 4]
            sage: a.rescale_row(1, 1/2)
            Traceback (most recent call last):
            ...
            TypeError: Rescaling row by Rational Field element cannot be done
            over Integer Ring, use change_ring or with_rescaled_row instead.

        To rescale the matrix by 1/2, you must change the base ring to the
        rationals::

            sage: a = a.change_ring(QQ); a
            [0 1 2]
            [3 4 4]
            sage: a.rescale_col(1, 1/2); a
            [  0 1/2   2]
            [  3   2   4]
        """
        self.check_row_bounds_and_mutability(i, i)
        try:
            s = self._coerce_element(s)
            self.rescale_row_c(i, s, start_col)
        except TypeError:
            raise TypeError('Rescaling row by %s element cannot be done over %s, use change_ring or with_rescaled_row instead.' % (s.parent(), self.base_ring()))

    cdef rescale_row_c(self, Py_ssize_t i, s, Py_ssize_t start_col):
        cdef Py_ssize_t j
        for j from start_col <= j < self._ncols:
            self.set_unsafe(i, j, self.get_unsafe(i, j)*s)

    def with_rescaled_row(self, Py_ssize_t i, s, Py_ssize_t start_col=0):
        """
        Replace `i`-th row of ``self`` by s times `i`-th row of self, returning
        new matrix.

        EXAMPLES: We rescale the second row of a matrix over the integers::

            sage: a = matrix(ZZ, 3, 2, range(6)); a
            [0 1]
            [2 3]
            [4 5]
            sage: b = a.with_rescaled_row(1, -2); b
            [ 0  1]
            [-4 -6]
            [ 4  5]

        The original matrix is unchanged::

            sage: a
            [0 1]
            [2 3]
            [4 5]

        Adding a rational multiple is okay, and reassigning a variable is
        okay::

            sage: a = a.with_rescaled_row(2, 1/3); a
            [  0   1]
            [  2   3]
            [4/3 5/3]
        """
        cdef Matrix temp
        self.check_row_bounds_and_mutability(i,i)
        try:
            s = self._coerce_element(s)
            temp = self.__copy__()
            temp.rescale_row_c(i, s, start_col)
            return temp
        # If scaling factor cannot be coerced, change the base ring to
        # one acceptable to both the original base ring and the scaling factor.
        except TypeError:
            temp = self.change_ring(Sequence([s,self.base_ring()(0)]).universe())
            s = temp._coerce_element(s)
            temp.rescale_row_c(i, s, start_col)
            return temp

    def rescale_col(self, Py_ssize_t i, s, Py_ssize_t start_row=0):
        """
        Replace `i`-th col of ``self`` by `s` times `i`-th col of ``self``.

        INPUT:

        - ``i`` -- `i`-th column

        - ``s`` -- scalar

        - ``start_row`` -- only rescale entries at this row
          and lower

        EXAMPLES: We rescale the last column of a matrix over the rational
        numbers::

            sage: a = matrix(QQ, 2, 3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.rescale_col(2, 1/2); a
            [  0   1   1]
            [  3   4 5/2]
            sage: R.<x> = QQ[]

        We rescale the last column of a matrix over a polynomial ring::

            sage: a = matrix(R, 2, 3, [1,x,x^2,x^3,x^4,x^5]); a
            [  1   x x^2]
            [x^3 x^4 x^5]
            sage: a.rescale_col(2, 1/2); a
            [      1       x 1/2*x^2]
            [    x^3     x^4 1/2*x^5]

        We try and fail to rescale a matrix over the integers by a
        non-integer::

            sage: a = matrix(ZZ, 2, 3, [0,1,2, 3,4,4]); a
            [0 1 2]
            [3 4 4]
            sage: a.rescale_col(2, 1/2)
            Traceback (most recent call last):
            ...
            TypeError: Rescaling column by Rational Field element cannot be done
            over Integer Ring, use change_ring or with_rescaled_col instead.

        To rescale the matrix by 1/2, you must change the base ring to the
        rationals::

            sage: a = a.change_ring(QQ); a
            [0 1 2]
            [3 4 4]
            sage: a.rescale_col(2,1/2); a
            [0 1 1]
            [3 4 2]
        """
        self.check_column_bounds_and_mutability(i, i)
        try:
            s = self._coerce_element(s)
            self.rescale_col_c(i, s, start_row)
        except TypeError:
            raise TypeError('Rescaling column by %s element cannot be done over %s, use change_ring or with_rescaled_col instead.' % (s.parent(), self.base_ring()))

    cdef rescale_col_c(self, Py_ssize_t i, s, Py_ssize_t start_row):
        cdef Py_ssize_t j
        for j from start_row <= j < self._nrows:
            self.set_unsafe(j, i, self.get_unsafe(j, i)*s)

    def with_rescaled_col(self, Py_ssize_t i, s, Py_ssize_t start_row=0):
        """
        Replaces `i`-th col of ``self`` by `s` times `i`-th col of self, returning
        new matrix.

        EXAMPLES: We rescale the last column of a matrix over the
        integers::

            sage: a = matrix(ZZ, 2, 3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: b = a.with_rescaled_col(2, -2); b
            [  0   1  -4]
            [  3   4 -10]

        The original matrix is unchanged::

            sage: a
            [0 1 2]
            [3 4 5]

        Adding a rational multiple is okay, and reassigning a variable is
        okay::

            sage: a = a.with_rescaled_col(1, 1/3); a
            [  0 1/3   2]
            [  3 4/3   5]
        """
        cdef Matrix temp
        self.check_column_bounds_and_mutability(i,i)
        try:
            s = self._coerce_element(s)
            temp = self.__copy__()
            temp.rescale_col_c(i, s, start_row)
            return temp
        # If scaling factor cannot be coerced, change the base ring to
        # one acceptable to both the original base ring and the scaling factor.
        except TypeError:
            temp = self.change_ring(Sequence([s,self.base_ring()(0)]).universe())
            s = temp._coerce_element(s)
            temp.rescale_col_c(i, s, start_row)
            return temp

    def set_row_to_multiple_of_row(self, Py_ssize_t i, Py_ssize_t j, s):
        """
        Set row i equal to s times row j.

        EXAMPLES: We change the second row to -3 times the first row::

            sage: a = matrix(ZZ, 2, 3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.set_row_to_multiple_of_row(1, 0, -3)
            sage: a
            [ 0  1  2]
            [ 0 -3 -6]

        If we try to multiply a row by a rational number, we get an error
        message::

            sage: a.set_row_to_multiple_of_row(1, 0, 1/2)
            Traceback (most recent call last):
            ...
            TypeError: Multiplying row by Rational Field element cannot be done over
            Integer Ring, use change_ring or with_row_set_to_multiple_of_row instead.
        """
        self.check_row_bounds_and_mutability(i, j)
        cdef Py_ssize_t n
        try:
            s = self._coerce_element(s)
            for n from 0 <= n < self._ncols:
                self.set_unsafe(i, n, s * self.get_unsafe(j, n))  # self[i] = s*self[j]
        except TypeError:
            raise TypeError('Multiplying row by %s element cannot be done over %s, use change_ring or with_row_set_to_multiple_of_row instead.' % (s.parent(), self.base_ring()))

    def with_row_set_to_multiple_of_row(self, Py_ssize_t i, Py_ssize_t j, s):
        """
        Set row i equal to s times row j, returning a new matrix.

        EXAMPLES: We change the second row to -3 times the first row::

            sage: a = matrix(ZZ, 2, 3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: b = a.with_row_set_to_multiple_of_row(1, 0, -3); b
            [ 0  1  2]
            [ 0 -3 -6]

        Note that the original matrix is unchanged::

            sage: a
            [0 1 2]
            [3 4 5]

        Adding a rational multiple is okay, and reassigning a variable is
        okay::

            sage: a = a.with_row_set_to_multiple_of_row(1, 0, 1/2); a
            [  0   1   2]
            [  0 1/2   1]
        """
        self.check_row_bounds_and_mutability(i, j)
        cdef Matrix temp
        cdef Py_ssize_t n
        try:
            s = self._coerce_element(s)
            temp = self.__copy__()
            for n from 0 <= n < temp._ncols:
                temp.set_unsafe(i, n, s * temp.get_unsafe(j, n))  # temp[i] = s*temp[j]
            return temp
        # If scaling factor cannot be coerced, change the base ring to
        # one acceptable to both the original base ring and the scaling factor.
        except TypeError:
            temp = self.change_ring(Sequence([s,self.base_ring()(0)]).universe())
            s = temp._coerce_element(s)
            for n from 0 <= n < temp._ncols:
                temp.set_unsafe(i, n, s * temp.get_unsafe(j, n))  # temp[i] = s*temp[j]
            return temp

    def set_col_to_multiple_of_col(self, Py_ssize_t i, Py_ssize_t j, s):
        """
        Set column i equal to s times column j.

        EXAMPLES: We change the second column to -3 times the first
        column.

        ::

            sage: a = matrix(ZZ, 2, 3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.set_col_to_multiple_of_col(1, 0, -3)
            sage: a
            [ 0  0  2]
            [ 3 -9  5]

        If we try to multiply a column by a rational number, we get an
        error message::

            sage: a.set_col_to_multiple_of_col(1, 0, 1/2)
            Traceback (most recent call last):
            ...
            TypeError: Multiplying column by Rational Field element cannot be done over Integer Ring, use change_ring or with_col_set_to_multiple_of_col instead.
        """
        self.check_column_bounds_and_mutability(i, j)
        cdef Py_ssize_t n
        try:
            s = self._coerce_element(s)
            for n from 0 <= n < self._nrows:
                self.set_unsafe(n, i, s * self.get_unsafe(n, j))
        # If scaling factor cannot be coerced, change the base ring to
        # one acceptable to both the original base ring and the scaling factor.
        except TypeError:
            raise TypeError('Multiplying column by %s element cannot be done over %s, use change_ring or with_col_set_to_multiple_of_col instead.' % (s.parent(), self.base_ring()))

    def with_col_set_to_multiple_of_col(self, Py_ssize_t i, Py_ssize_t j, s):
        """
        Set column i equal to s times column j, returning a new matrix.

        EXAMPLES: We change the second column to -3 times the first
        column.

        ::

            sage: a = matrix(ZZ, 2, 3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: b = a.with_col_set_to_multiple_of_col(1, 0, -3); b
            [ 0  0  2]
            [ 3 -9  5]

        Note that the original matrix is unchanged::

            sage: a
            [0 1 2]
            [3 4 5]

        Adding a rational multiple is okay, and reassigning a variable is
        okay::

            sage: a = a.with_col_set_to_multiple_of_col(1, 0, 1/2); a
            [  0   0   2]
            [  3 3/2   5]
        """
        self.check_column_bounds_and_mutability(i, j)
        cdef Py_ssize_t n
        cdef Matrix temp
        try:
            s = self._coerce_element(s)
            temp = self.__copy__()
            for n from 0 <= n < temp._nrows:
                temp.set_unsafe(n, i, s * temp.get_unsafe(n, j))
            return temp
        # If scaling factor cannot be coerced, change the base ring to
        # one acceptable to both the original base ring and the scaling factor.
        except TypeError:
            temp = self.change_ring(Sequence([s,self.base_ring()(0)]).universe())
            s = temp._coerce_element(s)
            for n from 0 <= n < temp._nrows:
                temp.set_unsafe(n, i, s * temp.get_unsafe(n, j))
            return temp

    def _set_row_to_negative_of_row_of_A_using_subset_of_columns(self, Py_ssize_t i, Matrix A,
                                                                 Py_ssize_t r, cols,
                                                                 cols_index=None):
        """
        Set row i of ``self`` to -(row r of A), but where we only take the
        given column positions in that row of A. We do not zero out the
        other entries of ``self``'s row i either.

        INPUT:

        - ``i`` -- integer, index into the rows of self

        - ``A`` -- a matrix

        - ``r`` -- integer, index into rows of A

        - ``cols`` -- a *sorted* list of integers

        - ``(cols_index`` -- ignored)

        EXAMPLES::

            sage: a = matrix(QQ,2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a._set_row_to_negative_of_row_of_A_using_subset_of_columns(0,a,1,[1,2])
            sage: a
            [-4 -5  2]
            [ 3  4  5]
        """
        self.check_row_bounds_and_mutability(i, i)
        if r < 0 or r >= A.nrows():
            raise IndexError("invalid row")
        # this function exists just because it is useful for modular symbols presentations.
        cdef Py_ssize_t l
        l = 0
        for k in cols:
            self.set_unsafe(i, l, -A.get_unsafe(r, k))  # self[i,l] = -A[r,k]
            l += 1

    def reverse_rows_and_columns(self):
        r"""
        Reverse the row order and column order of this matrix.

        This method transforms a matrix `m_{i,j}` with `0 \leq i < nrows` and
        `0 \leq j < ncols` into `m_{nrows - i - 1, ncols - j - 1}`.

        EXAMPLES::

            sage: m = matrix(ZZ, 2, 2, range(4))
            sage: m.reverse_rows_and_columns()
            sage: m
            [3 2]
            [1 0]

            sage: m = matrix(ZZ, 2, 3, range(6), sparse=True)
            sage: m.reverse_rows_and_columns()
            sage: m
            [5 4 3]
            [2 1 0]
            sage: m = matrix(ZZ, 3, 2, range(6), sparse=True)
            sage: m.reverse_rows_and_columns()
            sage: m
            [5 4]
            [3 2]
            [1 0]
            sage: m.reverse_rows_and_columns()
            sage: m
            [0 1]
            [2 3]
            [4 5]

            sage: m = matrix(QQ, 3, 2, [1/i for i in range(1,7)])
            sage: m.reverse_rows_and_columns()
            sage: m
            [1/6 1/5]
            [1/4 1/3]
            [1/2   1]

            sage: R.<x,y> = ZZ['x,y']
            sage: m = matrix(R, 3, 3, lambda i,j: x**i*y**j, sparse=True)
            sage: m.reverse_rows_and_columns()
            sage: m
            [x^2*y^2   x^2*y     x^2]
            [  x*y^2     x*y       x]
            [    y^2       y       1]

        If the matrix is immutable, the method raises an error::

            sage: m = matrix(ZZ, 2, [1, 3, -2, 4])
            sage: m.set_immutable()
            sage: m.reverse_rows_and_columns()
            Traceback (most recent call last):
            ...
            ValueError: matrix is immutable; please change a copy
            instead (i.e., use copy(M) to change a copy of M).
        """
        self.check_mutability()
        self.clear_cache()
        self._reverse_unsafe()

    ###################################################
    # Methods needed for quiver and cluster mutations
    # - mutate
    # - _travel_column
    # - is_symmetrizable
    # - is_skew_symmetrizable
    # - _check_symmetrizability
    #
    # AUTHORS:
    #     -- Christian Stump (Jun 2011)
    ###################################################

    def mutate(self, Py_ssize_t k ):
        """
        Mutates ``self`` at row and column index ``k``.

        .. warning:: Only makes sense if ``self`` is skew-symmetrizable.

        INPUT:

        - ``k`` -- integer at which row/column ``self`` is mutated

        EXAMPLES:

        Mutation of the B-matrix of the quiver of type `A_3`::

            sage: M = matrix(ZZ, 3, [0,1,0,-1,0,-1,0,1,0]); M
            [ 0  1  0]
            [-1  0 -1]
            [ 0  1  0]

            sage: M.mutate(0); M
            [ 0 -1  0]
            [ 1  0 -1]
            [ 0  1  0]

            sage: M.mutate(1); M
            [ 0  1 -1]
            [-1  0  1]
            [ 1 -1  0]

            sage: M = matrix(ZZ, 6, [0,1,0,-1,0,-1,0,1,0,1,0,0,0,1,0,0,0,1]); M
            [ 0  1  0]
            [-1  0 -1]
            [ 0  1  0]
            [ 1  0  0]
            [ 0  1  0]
            [ 0  0  1]

            sage: M.mutate(0); M
            [ 0 -1  0]
            [ 1  0 -1]
            [ 0  1  0]
            [-1  1  0]
            [ 0  1  0]
            [ 0  0  1]

        REFERENCES:

        - [FZ2001]_
        """
        cdef Py_ssize_t i, j, _
        cdef list pairs, k0_pairs, k1_pairs
        cdef bint ineg, jneg

        if k < 0 or k >= self._nrows or k >= self._ncols:
            raise IndexError("The mutation index is invalid")

        pairs = self.nonzero_positions()
        k0_pairs = [pair for pair in pairs if pair[0] == k]
        k1_pairs = [pair for pair in pairs if pair[1] == k]
        for _, j in k0_pairs:
            self[k, j] = -self.get_unsafe(k, j)
        for i,_ in k1_pairs:
            self[i, k] = -self.get_unsafe(i, k)

        for i,_ in k1_pairs:
            ik = self.get_unsafe(i, k)
            ineg = (ik < 0)
            for _, j in k0_pairs:
                kj = self.get_unsafe(k, j)
                jneg = (kj < 0)
                if ineg and jneg:
                    self[i, j] = self.get_unsafe(i, j) + self.get_unsafe(i, k)*self.get_unsafe(k, j)
                elif not ineg and not jneg:
                    self[i, j] = self.get_unsafe(i, j) - self.get_unsafe(i, k)*self.get_unsafe(k, j)

    def _travel_column( self, dict d, int k, int sign, positive ):
        r"""
        Helper function for testing symmetrizability. Tests dependencies within entries in ``self`` and entries in the dictionary ``d``.

        .. warning:: the dictionary ``d`` gets new values for keys in L.

        INPUT:

        - ``d`` -- dictionary modelling partial entries of a diagonal matrix

        - ``k`` -- integer for which row and column of ``self`` should be tested with the dictionary d

        - ``sign`` -- `\pm 1`, depending on symmetric or skew-symmetric is tested

        - ``positive`` -- if ``True``, only positive entries for the values of the dictionary are allowed

        OUTPUT: ``L`` -- list of new keys in d

        EXAMPLES::

            sage: M = matrix(ZZ, 3, [0,1,0,-1,0,-1,0,1,0]); M
            [ 0  1  0]
            [-1  0 -1]
            [ 0  1  0]

            sage: M._travel_column({0: 1}, 0, -1, True)
            [1]
        """
        cdef list L = []
        cdef int i

        for i from 0 <= i < self._ncols:
            if i not in d:
                self_ik = self.get_unsafe(i,k)
                self_ki = self.get_unsafe(k,i)
                if bool(self_ik) != bool(self_ki):
                    return False
                if self_ik != 0:
                    L.append(i)
                    d[i] = sign * d[k] * self_ki / self_ik
                    if positive and not d[i] > 0:
                        return False
                    for j in d:
                        if d[i] * self.get_unsafe(i, j) != sign * d[j] * self.get_unsafe(j, i):
                            return False
        return L

    def _check_symmetrizability(self, return_diag=False, skew=False, positive=True):
        r"""
        This function takes a square matrix over an *ordered integral domain*
        and checks if it is (skew-)symmetrizable.

        A matrix `B` is (skew-)symmetrizable iff there exists an invertible
        diagonal matrix `D` such that `DB` is (skew-)symmetric.

        INPUT:

        - ``return_diag`` -- boolean (default: ``False``); if ``True`` and
          ``self`` is (skew)-symmetrizable the diagonal entries of the matrix
          `D` are returned
        - ``skew`` -- boolean (default: ``False``); if ``True``,
          (skew-)symmetrizability is checked
        - ``positive`` -- boolean (default: ``True``); if ``True``, the
          condition that `D` has positive entries is added

        OUTPUT:

        - ``True`` -- if ``self`` is (skew-)symmetrizable and ``return_diag``
          is ``False``
        - the diagonal entries of the matrix `D` such that `DB` is
          (skew-)symmetric -- iff ``self`` is (skew-)symmetrizable and
          ``return_diag`` is ``True``
        - ``False`` -- iff ``self`` is not (skew-)symmetrizable

        EXAMPLES::

            sage: matrix([[0,6],[3,0]])._check_symmetrizability(positive=False)
            True
            sage: matrix([[0,6],[3,0]])._check_symmetrizability(positive=True)
            True
            sage: matrix([[0,6],[3,0]])._check_symmetrizability(skew=True, positive=False)
            True
            sage: matrix([[0,6],[3,0]])._check_symmetrizability(skew=True, positive=True)
            False

        REFERENCES:

        - [FZ2001]_
        """
        cdef dict d = {}
        cdef list queue = list(range(self._ncols))
        cdef int l, sign, i

        if skew:
            # testing the diagonal entries to be zero
            for i from 0 <= i < self._nrows:
                if not self.get_is_zero_unsafe(i,i):
                    return False
            sign = -1
        else:
            sign = 1

        while queue:
            i = queue.pop(0)
            d[i] = 1
            L = self._travel_column( d, i, sign, positive )
            if L is False:
                return False
            while L:
                l = L.pop(0)
                queue.remove( l )
                L_prime = self._travel_column( d, l, sign, positive )
                if L_prime is False:
                    return False
                else:
                    L.extend( L_prime )
        if return_diag:
            return [d[i] for i in range(self._nrows)]
        else:
            return True

    ###################################################
    # Matrix-vector multiply
    ###################################################
    def linear_combination_of_rows(self, v):
        r"""
        Return the linear combination of the rows of ``self`` given by the
        coefficients in the list ``v``.

        INPUT:

        - ``v`` -- a list of scalars.  The length can be less than
          the number of rows of ``self`` but not greater.

        OUTPUT:

        The vector (or free module element) that is a linear
        combination of the rows of ``self``. If the list of
        scalars has fewer entries than the number of rows,
        additional zeros are appended to the list until it
        has as many entries as the number of rows.

        EXAMPLES::

            sage: a = matrix(ZZ, 2, 3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.linear_combination_of_rows([1,2])
            (6, 9, 12)

            sage: a.linear_combination_of_rows([0,0])
            (0, 0, 0)

            sage: a.linear_combination_of_rows([1/2,2/3])
            (2, 19/6, 13/3)

        The list ``v`` can be anything that is iterable.  Perhaps most
        naturally, a vector may be used. ::

            sage: v = vector(ZZ, [1,2])
            sage: a.linear_combination_of_rows(v)
            (6, 9, 12)

        We check that a matrix with no rows behaves properly. ::

            sage: matrix(QQ, 0, 2).linear_combination_of_rows([])
            (0, 0)

        The object returned is a vector, or a free module element. ::

            sage: B = matrix(ZZ, 4, 3, range(12))
            sage: w = B.linear_combination_of_rows([-1,2,-3,4])
            sage: w
            (24, 26, 28)
            sage: w.parent()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: x = B.linear_combination_of_rows([1/2,1/3,1/4,1/5])
            sage: x
            (43/10, 67/12, 103/15)
            sage: x.parent()
            Vector space of dimension 3 over Rational Field

        The length of v can be less than the number of rows, but not
        greater. ::

            sage: A = matrix(QQ, 3, 4, range(12))
            sage: A.linear_combination_of_rows([2,3])
            (12, 17, 22, 27)
            sage: A.linear_combination_of_rows([1,2,3,4])
            Traceback (most recent call last):
            ...
            ValueError: length of v must be at most the number of rows of self
        """
        if len(v) > self._nrows:
            raise ValueError("length of v must be at most the number of rows of self")
        if not self._nrows:
            return self.parent().row_space().zero_vector()
        from sage.matrix.constructor import matrix
        v = matrix(list(v)+[0]*(self._nrows-len(v)))
        return (v * self)[0]

    def linear_combination_of_columns(self, v):
        r"""
        Return the linear combination of the columns of ``self`` given by the
        coefficients in the list ``v``.

        INPUT:

        - ``v`` -- a list of scalars.  The length can be less than
          the number of columns of ``self`` but not greater.

        OUTPUT:

        The vector (or free module element) that is a linear
        combination of the columns of ``self``. If the list of
        scalars has fewer entries than the number of columns,
        additional zeros are appended to the list until it
        has as many entries as the number of columns.

        EXAMPLES::

            sage: a = matrix(ZZ, 2, 3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.linear_combination_of_columns([1,1,1])
            (3, 12)

            sage: a.linear_combination_of_columns([0,0,0])
            (0, 0)

            sage: a.linear_combination_of_columns([1/2,2/3,3/4])
            (13/6, 95/12)

        The list ``v`` can be anything that is iterable.  Perhaps most
        naturally, a vector may be used. ::

            sage: v = vector(ZZ, [1,2,3])
            sage: a.linear_combination_of_columns(v)
            (8, 26)

        We check that a matrix with no columns behaves properly. ::

            sage: matrix(QQ, 2, 0).linear_combination_of_columns([])
            (0, 0)

        The object returned is a vector, or a free module element. ::

            sage: B = matrix(ZZ, 4, 3, range(12))
            sage: w = B.linear_combination_of_columns([-1,2,-3])
            sage: w
            (-4, -10, -16, -22)
            sage: w.parent()
            Ambient free module of rank 4 over the principal ideal domain Integer Ring
            sage: x = B.linear_combination_of_columns([1/2,1/3,1/4])
            sage: x
            (5/6, 49/12, 22/3, 127/12)
            sage: x.parent()
            Vector space of dimension 4 over Rational Field

        The length of v can be less than the number of columns, but not
        greater. ::

            sage: A = matrix(QQ, 3, 5, range(15))
            sage: A.linear_combination_of_columns([1,-2,3,-4])
            (-8, -18, -28)
            sage: A.linear_combination_of_columns([1,2,3,4,5,6])
            Traceback (most recent call last):
            ...
            ValueError: length of v must be at most the number of columns of self
        """
        if len(v) > self._ncols:
            raise ValueError("length of v must be at most the number of columns of self")
        if not self._ncols:
            return self.parent().column_space().zero_vector()
        from sage.matrix.constructor import matrix
        v = matrix(self._ncols, 1, list(v)+[0]*(self._ncols-len(v)))
        return (self * v).column(0)

    ###################################################
    # Predicates
    ###################################################

    def is_symmetric(self):
        """
        Return ``True`` if this is a symmetric matrix.

        A symmetric matrix is necessarily square.

        EXAMPLES::

            sage: m = Matrix(QQ, 2, range(0,4))
            sage: m.is_symmetric()
            False

            sage: m = Matrix(QQ, 2, (1,1,1,1,1,1))
            sage: m.is_symmetric()
            False

            sage: m = Matrix(QQ, 1, (2,))
            sage: m.is_symmetric()
            True
        """
        if self._ncols != self._nrows: return False
        # could be bigger than an int on a 64-bit platform, this
        #  is the type used for indexing.
        cdef Py_ssize_t i, j

        for i from 0 <= i < self._nrows:
            for j from 0 <= j < i:
                if self.get_unsafe(i, j) != self.get_unsafe(j, i):
                    return False
        return True

    def _is_hermitian(self, skew, tolerance):
        r"""
        Return ``True`` if the matrix is (skew-)Hermitian up to the
        entry-wise ``tolerance``.

        For internal purposes. This function is used to implement both
        the :meth:`is_hermitian` and :meth:`is_skew_hermitian` methods.

        INPUT:

        - ``skew`` -- boolean (default: ``False``); set to ``True`` to
          check if the matrix is skew-Hermitian instead of Hermitian

        - ``tolerance`` -- a real number; the maximum difference we'll
          tolerate between entries of the given matrix and its conjugate-
          transpose.

        OUTPUT:

        ``True`` if the matrix is square and (skew-)Hermitian, and
        ``False`` otherwise.

        Note that if conjugation has no effect on elements of the base
        ring (such as for integers), then the :meth:`is_(skew_)symmetric`
        method is equivalent and faster.

        The result is cached.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [[ 1 + I,  1 - 6*I, -1 - I],
            ....:                    [-3 - I,     -4*I,     -2],
            ....:                    [-1 + I, -2 - 8*I,  2 + I]])
            sage: A._is_hermitian(skew=False, tolerance=0)
            False
            sage: B = A*A.conjugate_transpose()
            sage: B._is_hermitian(skew=False, tolerance=0)
            True

        Sage has several fields besides the entire complex numbers
        where conjugation is non-trivial::

            sage: # needs sage.rings.number_field
            sage: F.<b> = QuadraticField(-7)
            sage: C = matrix(F, [[-2*b - 3,  7*b - 6, -b + 3],
            ....:                [-2*b - 3, -3*b + 2,   -2*b],
            ....:                [   b + 1,        0,     -2]])
            sage: C._is_hermitian(skew=False, tolerance=0)
            False
            sage: C = C*C.conjugate_transpose()
            sage: C._is_hermitian(skew=False, tolerance=0)
            True

        A matrix that is nearly Hermitian, but for a non-real
        diagonal entry::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [[    2,   2-I, 1+4*I],
            ....:                    [  2+I,   3+I, 2-6*I],
            ....:                    [1-4*I, 2+6*I,     5]])
            sage: A._is_hermitian(skew=False, tolerance=0)
            False
            sage: A[1, 1] = 132
            sage: A._is_hermitian(skew=False, tolerance=0)
            True

        Rectangular matrices are never Hermitian::

            sage: A = matrix(QQbar, 3, 4)                                               # needs sage.rings.number_field
            sage: A._is_hermitian(skew=False, tolerance=0)                              # needs sage.rings.number_field
            False

        A square, empty matrix is trivially Hermitian::

            sage: A = matrix(QQ, 0, 0)
            sage: A._is_hermitian(skew=False, tolerance=0)
            True

        A matrix that is skew-Hermitian::

            sage: A = matrix(QQbar, [[-I, 2+I], [-2+I, 0]])                             # needs sage.rings.number_field
            sage: A._is_hermitian(skew=False, tolerance=0)                              # needs sage.rings.number_field
            False
            sage: A._is_hermitian(skew=True, tolerance=0)
            True
        """
        key = ("_is_hermitian", skew, tolerance)

        cached = self.fetch(key)
        if cached is not None:
            return cached
        if not self.is_square():
            self.cache(key, False)
            return False
        if self._nrows == 0:
            self.cache(key, True)
            return True

        s = 1
        if skew:
            s = -1

        tolerance = self.base_ring()(tolerance)
        cdef bint tolerance_is_zero = tolerance.is_zero()
        cdef Py_ssize_t i, j

        if self.is_sparse_c():
            # The dense algorithm checks all of the on-or-below-diagonal
            # entries, of which there are (n^2 + n)/2. If the matrix
            # is sparse, however, we can get away with checking only
            # the nonzero positions. This will be faster if the matrix
            # is truly sparse (if there are not so many of those positions)
            # even after taking numerical issues into account.
            #
            # We access this list of entries directly, without making a
            # copy, so it's important that we don't modify it.
            entries = self._nonzero_positions_by_row(copy=False)
        else:
            entries = ((i, j) for i in range(self._nrows)
                       for j in range(i + 1))

        for (i, j) in entries:
            entry_a = self.get_unsafe(i, j)
            entry_b = s*self.get_unsafe(j, i).conjugate()

            if tolerance_is_zero:
                # When the tolerance is exactly zero, as will
                # usually be the case for exact rings, testing for
                # literal equality provides a simple answer to the
                # question of how we should test against the
                # tolerance in rings such as finite fields and
                # polynomials where abs/norm support is spotty and
                # an ordering may not be intelligently defined.
                if entry_a != entry_b:
                    self.cache(key, False)
                    return False
            else:
                d = entry_a - entry_b
                # sqrt() can have a different parent, and doesn't
                # preserve order in e.g. finite fields, so we
                # square both sides of the usual test here.
                if (d*d.conjugate()) > tolerance**2:
                    self.cache(key, False)
                    return False

        self.cache(key, True)
        return True

    def is_hermitian(self):
        r"""
        Return ``True`` if the matrix is equal to its conjugate-transpose.

        OUTPUT:

        ``True`` if the matrix is square and equal to the transpose with
        every entry conjugated, and ``False`` otherwise.

        Note that if conjugation has no effect on elements of the base
        ring (such as for integers), then the :meth:`is_symmetric`
        method is equivalent and faster.

        This routine is for matrices over exact rings and so may not
        work properly for matrices over ``RR`` or ``CC``.  For matrices with
        approximate entries, the rings of double-precision floating-point
        numbers, ``RDF`` and ``CDF``, are a better choice since the
        :meth:`sage.matrix.matrix_double_dense.Matrix_double_dense.is_hermitian`
        method has a tolerance parameter.  This provides control over
        allowing for minor discrepancies between entries when checking
        equality.

        The result is cached.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [[ 1 + I,  1 - 6*I, -1 - I],
            ....:                    [-3 - I,     -4*I,     -2],
            ....:                    [-1 + I, -2 - 8*I,  2 + I]])
            sage: A.is_hermitian()
            False
            sage: B = A * A.conjugate_transpose()
            sage: B.is_hermitian()
            True

        Sage has several fields besides the entire complex numbers
        where conjugation is non-trivial. ::

            sage: # needs sage.rings.number_field
            sage: F.<b> = QuadraticField(-7)
            sage: C = matrix(F, [[-2*b - 3,  7*b - 6, -b + 3],
            ....:                [-2*b - 3, -3*b + 2,   -2*b],
            ....:                [   b + 1,        0,     -2]])
            sage: C.is_hermitian()
            False
            sage: C = C*C.conjugate_transpose()
            sage: C.is_hermitian()
            True

        A matrix that is nearly Hermitian, but for a non-real
        diagonal entry. ::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [[    2,   2-I, 1+4*I],
            ....:                    [  2+I,   3+I, 2-6*I],
            ....:                    [1-4*I, 2+6*I,     5]])
            sage: A.is_hermitian()
            False
            sage: A[1, 1] = 132
            sage: A.is_hermitian()
            True

        Rectangular matrices are never Hermitian.  ::

            sage: A = matrix(QQbar, 3, 4)                                               # needs sage.rings.number_field
            sage: A.is_hermitian()                                                      # needs sage.rings.number_field
            False

        A square, empty matrix is trivially Hermitian.  ::

            sage: A = matrix(QQ, 0, 0)
            sage: A.is_hermitian()
            True
        """
        return self._is_hermitian(skew=False, tolerance=0)

    def is_skew_hermitian(self):
        r"""
        Return ``True`` if the matrix is equal to the negative of its
        conjugate transpose.

        OUTPUT:

        ``True`` if the matrix is square and equal to the negative of
        its conjugate transpose, and ``False`` otherwise.

        Note that if conjugation has no effect on elements of the base
        ring (such as for integers), then the :meth:`is_skew_symmetric`
        method is equivalent and faster.

        This routine is for matrices over exact rings and so may not
        work properly for matrices over ``RR`` or ``CC``.  For matrices with
        approximate entries, the rings of double-precision floating-point
        numbers, ``RDF`` and ``CDF``, are a better choice since the
        :meth:`sage.matrix.matrix_double_dense.Matrix_double_dense.is_skew_hermitian`
        method has a tolerance parameter.  This provides control over
        allowing for minor discrepancies between entries when checking
        equality.

        The result is cached.

        EXAMPLES::

            sage: A = matrix(QQbar, [[0, -1],                                           # needs sage.rings.number_field
            ....:                    [1,  0]])
            sage: A.is_skew_hermitian()                                                 # needs sage.rings.number_field
            True

        A matrix that is nearly skew-Hermitian, but for a non-real
        diagonal entry. ::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [[  -I, -1, 1-I],
            ....:                    [   1,  1,  -1],
            ....:                    [-1-I,  1,  -I]])
            sage: A.is_skew_hermitian()
            False
            sage: A[1, 1] = -I
            sage: A.is_skew_hermitian()
            True

        Rectangular matrices are never skew-Hermitian. ::

            sage: A = matrix(QQbar, 3, 4)                                               # needs sage.rings.number_field
            sage: A.is_skew_hermitian()                                                 # needs sage.rings.number_field
            False

        A square, empty matrix is trivially Hermitian. ::

            sage: A = matrix(QQ, 0, 0)
            sage: A.is_skew_hermitian()
            True
        """
        return self._is_hermitian(skew=True, tolerance=0)

    def is_skew_symmetric(self):
        """
        Return ``True`` if ``self`` is a skew-symmetric matrix.

        Here, "skew-symmetric matrix" means a square matrix `A`
        satisfying `A^T = -A`. It does not require that the
        diagonal entries of `A` are `0` (although this
        automatically follows from `A^T = -A` when `2` is
        invertible in the ground ring over which the matrix is
        considered). Skew-symmetric matrices `A` whose diagonal
        entries are `0` are said to be "alternating", and this
        property is checked by the :meth:`is_alternating`
        method.

        EXAMPLES::

            sage: m = matrix(QQ, [[0,2], [-2,0]])
            sage: m.is_skew_symmetric()
            True
            sage: m = matrix(QQ, [[1,2], [2,1]])
            sage: m.is_skew_symmetric()
            False

        Skew-symmetric is not the same as alternating when
        `2` is a zero-divisor in the ground ring::

            sage: n = matrix(Zmod(4), [[0, 1], [-1, 2]])
            sage: n.is_skew_symmetric()
            True

        but yet the diagonal cannot be completely
        arbitrary in this case::

            sage: n = matrix(Zmod(4), [[0, 1], [-1, 3]])
            sage: n.is_skew_symmetric()
            False
        """
        if self._ncols != self._nrows: return False
        # could be bigger than an int on a 64-bit platform, this
        #  is the type used for indexing.
        cdef Py_ssize_t i, j

        for i from 0 <= i < self._nrows:
            for j from 0 <= j <= i:
                if self.get_unsafe(i, j) != -self.get_unsafe(j, i):
                    return False
        return True

    def is_alternating(self):
        """
        Return ``True`` if ``self`` is an alternating matrix.

        Here, "alternating matrix" means a square matrix `A`
        satisfying `A^T = -A` and such that the diagonal entries
        of `A` are `0`. Notice that the condition that the
        diagonal entries be `0` is not redundant for matrices over
        arbitrary ground rings (but it is redundant when `2` is
        invertible in the ground ring). A square matrix `A` only
        required to satisfy `A^T = -A` is said to be
        "skew-symmetric", and this property is checked by the
        :meth:`is_skew_symmetric` method.

        EXAMPLES::

            sage: m = matrix(QQ, [[0,2], [-2,0]])
            sage: m.is_alternating()
            True
            sage: m = matrix(QQ, [[1,2], [2,1]])
            sage: m.is_alternating()
            False

        In contrast to the property of being skew-symmetric, the
        property of being alternating does not tolerate nonzero
        entries on the diagonal even if they are their own
        negatives::

            sage: n = matrix(Zmod(4), [[0, 1], [-1, 2]])
            sage: n.is_alternating()
            False
        """
        if self._ncols != self._nrows: return False
        # could be bigger than an int on a 64-bit platform, this
        #  is the type used for indexing.
        cdef Py_ssize_t i, j

        zero = self._base_ring.zero()
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < i:
                if self.get_unsafe(i, j) != -self.get_unsafe(j, i):
                    return False
            if not self.get_unsafe(i, i) == zero:
                return False
        return True

    def is_symmetrizable(self, return_diag=False, positive=True):
        r"""
        This function takes a square matrix over an *ordered integral domain*
        and checks if it is symmetrizable.

        A matrix `B` is symmetrizable iff there exists an invertible diagonal
        matrix `D` such that `DB` is symmetric.

        .. warning:: Expects ``self`` to be a matrix over an *ordered integral domain*.

        INPUT:

        - ``return_diag`` -- boolean (default: ``False``); if ``True`` and
          ``self`` is symmetrizable the diagonal entries of the matrix `D` are
          returned
        - ``positive`` -- boolean (default: ``True``); if ``True``, the
          condition that `D` has positive entries is added

        OUTPUT:

        - ``True`` -- if ``self`` is symmetrizable and ``return_diag`` is
          ``False``
        - the diagonal entries of a matrix `D` such that `DB` is symmetric --
          iff ``self`` is symmetrizable and ``return_diag`` is ``True``
        - ``False`` -- iff ``self`` is not symmetrizable

        EXAMPLES::

            sage: matrix([[0,6],[3,0]]).is_symmetrizable(positive=False)
            True

            sage: matrix([[0,6],[3,0]]).is_symmetrizable(positive=True)
            True

            sage: matrix([[0,6],[0,0]]).is_symmetrizable(return_diag=True)
            False

            sage: matrix([2]).is_symmetrizable(positive=True)
            True

            sage: matrix([[1,2],[3,4]]).is_symmetrizable(return_diag=true)
            [1, 2/3]

        REFERENCES:

        - [FZ2001]_
        """
        if self._ncols != self._nrows:
            raise ValueError("The matrix is not a square matrix")
        return self._check_symmetrizability(return_diag=return_diag, skew=False, positive=positive)

    def is_skew_symmetrizable(self, return_diag=False, positive=True):
        r"""
        This function takes a square matrix over an *ordered integral domain*
        and checks if it is skew-symmetrizable.
        A matrix `B` is skew-symmetrizable iff there exists an invertible
        diagonal matrix `D` such that `DB` is skew-symmetric.

        .. warning:: Expects ``self`` to be a matrix over an *ordered integral domain*.

        INPUT:

        - ``return_diag`` -- boolean (default: ``False``); if ``True`` and
          ``self`` is skew-symmetrizable the diagonal entries of the matrix `D`
          are returned
        - ``positive`` -- boolean (default: ``True``); if ``True``, the
          condition that `D` has positive entries is added

        OUTPUT:

        - ``True`` -- if ``self`` is skew-symmetrizable and ``return_diag`` is
          ``False``
        - the diagonal entries of a matrix `D` such that `DB` is
          skew-symmetric -- iff ``self`` is skew-symmetrizable and ``return_diag``
          is ``True``
        - ``False`` -- iff ``self`` is not skew-symmetrizable

        EXAMPLES::

            sage: matrix([[0,6],[3,0]]).is_skew_symmetrizable(positive=False)
            True
            sage: matrix([[0,6],[3,0]]).is_skew_symmetrizable(positive=True)
            False

            sage: M = matrix(4, [0,1,0,0, -1,0,-1,0, 0,2,0,1, 0,0,-1,0]); M
            [ 0  1  0  0]
            [-1  0 -1  0]
            [ 0  2  0  1]
            [ 0  0 -1  0]

            sage: M.is_skew_symmetrizable(return_diag=True)
            [1, 1, 1/2, 1/2]

            sage: M2 = diagonal_matrix([1,1,1/2,1/2]) * M; M2
            [   0    1    0    0]
            [  -1    0   -1    0]
            [   0    1    0  1/2]
            [   0    0 -1/2    0]

            sage: M2.is_skew_symmetric()
            True

        REFERENCES:

        - [FZ2001]_
        """
        if self._ncols != self._nrows:
            raise ValueError("The matrix is not a square matrix")
        return self._check_symmetrizability(return_diag=return_diag, skew=True, positive=positive)

    def is_dense(self):
        """
        Return ``True`` if this is a dense matrix.

        In Sage, being dense is a property of the underlying
        representation, not the number of nonzero entries.

        EXAMPLES::

            sage: matrix(QQ, 2, 2, range(4)).is_dense()
            True
            sage: matrix(QQ, 2, 2, range(4), sparse=True).is_dense()
            False
        """
        return self.is_dense_c()

    def is_sparse(self):
        """
        Return ``True`` if this is a sparse matrix.

        In Sage, being sparse is a property of the underlying
        representation, not the number of nonzero entries.

        EXAMPLES::

            sage: matrix(QQ, 2, 2, range(4)).is_sparse()
            False
            sage: matrix(QQ, 2, 2, range(4), sparse=True).is_sparse()
            True
        """
        return self.is_sparse_c()

    def is_square(self):
        """
        Return ``True`` precisely if this matrix is square, i.e., has the same
        number of rows and columns.

        EXAMPLES::

            sage: matrix(QQ, 2, 2, range(4)).is_square()
            True
            sage: matrix(QQ, 2, 3, range(6)).is_square()
            False
        """
        return self._nrows == self._ncols

    def is_invertible(self):
        r"""
        Return ``True`` if this matrix is invertible.

        EXAMPLES: The following matrix is invertible over
        `\QQ` but not over `\ZZ`.

        ::

            sage: A = MatrixSpace(ZZ, 2)(range(4))
            sage: A.is_invertible()
            False
            sage: A.matrix_over_field().is_invertible()
            True

        The inverse function is a constructor for matrices over the
        fraction field, so it can work even if A is not invertible.

        ::

            sage: ~A   # inverse of A
            [-3/2  1/2]
            [   1    0]

        The next matrix is invertible over `\ZZ`.

        ::

            sage: A = MatrixSpace(IntegerRing(), 2)([1,10,0,-1])
            sage: A.is_invertible()
            True
            sage: ~A                # compute the inverse
            [ 1 10]
            [ 0 -1]

        The following nontrivial matrix is invertible over
        `\ZZ[x]`.

        ::

            sage: R.<x> = PolynomialRing(IntegerRing())
            sage: A = MatrixSpace(R, 2)([1,x,0,-1])
            sage: A.is_invertible()
            True
            sage: ~A
            [ 1  x]
            [ 0 -1]
        """
        return self.is_square() and self.determinant().is_unit()

    is_unit = is_invertible

    def is_singular(self):
        r"""
        Return ``True`` if ``self`` is singular.

        OUTPUT:

        A square matrix is singular if it has a zero
        determinant and this method will return ``True``
        in exactly this case. When the entries of the
        matrix come from a field, this is equivalent
        to having a nontrivial kernel, or lacking an
        inverse, or having linearly dependent rows,
        or having linearly dependent columns.

        For square matrices over a field the methods
        :meth:`is_invertible` and :meth:`is_singular`
        are logical opposites.  However, it is an error
        to apply :meth:`is_singular` to a matrix that
        is not square, while :meth:`is_invertible` will
        always return ``False`` for a matrix that is not
        square.

        EXAMPLES:

        A singular matrix over the field ``QQ``. ::

            sage: A = matrix(QQ, 4, [-1,2,-3,6, 0,-1,-1,0, -1,1,-5,7, -1,6,5,2])
            sage: A.is_singular()
            True
            sage: A.right_kernel().dimension()
            1

        A matrix that is not singular, i.e. nonsingular, over a field. ::

            sage: B = matrix(QQ, 4, [1,-3,-1,-5, 2,-5,-2,-7, -2,5,3,4, -1,4,2,6])
            sage: B.is_singular()
            False
            sage: B.left_kernel().dimension()
            0

        For *rectangular* matrices, invertibility is always
        ``False``, but asking about singularity will give an error. ::

            sage: C = matrix(QQ, 5, range(30))
            sage: C.is_invertible()
            False
            sage: C.is_singular()
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix

        When the base ring is not a field, then a matrix
        may be both not invertible and not singular. ::

            sage: D = matrix(ZZ, 4, [2,0,-4,8, 2,1,-2,7, 2,5,7,0, 0,1,4,-6])
            sage: D.is_invertible()
            False
            sage: D.is_singular()
            False
            sage: d = D.determinant(); d
            2
            sage: d.is_unit()
            False
        """
        if self._ncols == self._nrows:
            return self.rank() != self._nrows
        else:
            raise ValueError("self must be a square matrix")

    ###################################################
    # Invariants of a matrix
    ###################################################

    def pivots(self):
        """
        Return the pivot column positions of this matrix.

        OUTPUT: a tuple of Python integers: the position of the
        first nonzero entry in each row of the echelon form.

        This returns a tuple so it is immutable; see :issue:`10752`.

        EXAMPLES::

            sage: A = matrix(QQ, 2, 2, range(4))
            sage: A.pivots()
            (0, 1)
        """
        x = self.fetch('pivots')
        if x is not None:
            return tuple(x)
        self.echelon_form()
        x = self.fetch('pivots')
        if x is None:
            print(self)
            print(self.nrows())
            print(self.dict())
            raise RuntimeError("BUG: matrix pivots should have been set but weren't, matrix parent = '%s'" % self.parent())
        return tuple(x)

    def rank(self):
        """
        Return the rank of this matrix.

        EXAMPLES::

            sage: m = matrix(GF(7), 5, range(25))
            sage: m.rank()
            2

        Rank is not implemented over the integers modulo a composite yet.::

            sage: m = matrix(Integers(4), 2, [2,2,2,2])
            sage: m.rank()
            Traceback (most recent call last):
            ...
            NotImplementedError: Echelon form not implemented over 'Ring of integers modulo 4'.

        TESTS:

        We should be able to compute the rank of a matrix whose
        entries are polynomials over a finite field (:issue:`5014`)::

            sage: P.<x> = PolynomialRing(GF(17))
            sage: m = matrix(P, [[ 6*x^2 + 8*x + 12, 10*x^2 + 4*x + 11],
            ....:                [8*x^2 + 12*x + 15,  8*x^2 + 9*x + 16]])
            sage: m.rank()
            2
        """
        x = self.fetch('rank')
        if x is not None:
            return x
        if self._nrows == 0 or self._ncols == 0:
            return 0
        r = len(self.pivots())
        self.cache('rank', r)
        return r

    def nonpivots(self):
        """
        Return the list of `i` such that the `i`-th column of ``self`` is NOT a
        pivot column of the reduced row echelon form of ``self``.

        OUTPUT: sorted tuple of (Python) integers

        EXAMPLES::

            sage: a = matrix(QQ, 3, 3, range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a.echelon_form()
            [ 1  0 -1]
            [ 0  1  2]
            [ 0  0  0]
            sage: a.nonpivots()
            (2,)
        """
        x = self.fetch('nonpivots')
        if x is not None:
            return tuple(x)

        X = set(self.pivots())
        np = []
        for j in range(self.ncols()):
            if j not in X:
                np.append(j)
        np = tuple(np)
        self.cache('nonpivots',np)
        return np

    def nonzero_positions(self, copy=True, column_order=False):
        r"""
        Return the sorted list of pairs ``(i,j)`` such that ``self[i,j] != 0``.

        INPUT:

        - ``copy`` -- boolean (default: ``True``); it is safe to change the
          resulting list (unless you give the option ``copy=False``)

        - ``column_order`` -- boolean (default: ``False``); if ``True``,
          returns the list of pairs ``(i,j)`` such that ``self[i,j] != 0``, but
          sorted by columns, i.e., column ``j=0`` entries occur first, then
          column ``j=1`` entries, etc.

        EXAMPLES::

            sage: a = matrix(QQ, 2,3, [1,2,0,2,0,0]); a
            [1 2 0]
            [2 0 0]
            sage: a.nonzero_positions()
            [(0, 0), (0, 1), (1, 0)]
            sage: a.nonzero_positions(copy=False)
            [(0, 0), (0, 1), (1, 0)]
            sage: a.nonzero_positions(column_order=True)
            [(0, 0), (1, 0), (0, 1)]
            sage: a = matrix(QQ, 2,3, [1,2,0,2,0,0], sparse=True); a
            [1 2 0]
            [2 0 0]
            sage: a.nonzero_positions()
            [(0, 0), (0, 1), (1, 0)]
            sage: a.nonzero_positions(copy=False)
            [(0, 0), (0, 1), (1, 0)]
            sage: a.nonzero_positions(column_order=True)
            [(0, 0), (1, 0), (0, 1)]
        """
        if column_order:
            return self._nonzero_positions_by_column(copy)
        else:
            return self._nonzero_positions_by_row(copy)

    def _nonzero_positions_by_row(self, copy=True):
        """
        Return the list of pairs ``(i,j)`` such that ``self[i,j] != 0``.

        It is safe to change the resulting list (unless you give the
        option ``copy=False``).

        EXAMPLES::

            sage: M = Matrix(CC, [[1,0],[0,1]], sparse=True)
            sage: M._nonzero_positions_by_row()
            [(0, 0), (1, 1)]
        """
        x = self.fetch('nonzero_positions')
        if x is not None:
            if copy:
                return list(x)
            return x
        cdef Py_ssize_t i, j
        nzp = []
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                if not self.get_is_zero_unsafe(i, j):
                    nzp.append((i, j))
        self.cache('nonzero_positions', nzp)
        if copy:
            return list(nzp)
        return nzp

    def _nonzero_positions_by_column(self, copy=True):
        """
        Return the list of pairs ``(i,j)`` such that ``self[i,j] != 0``, but
        sorted by columns, i.e., column ``j=0`` entries occur first, then
        column ``j=1`` entries, etc.

        It is safe to change the resulting list (unless you give the option
        ``copy=False``).

        EXAMPLES::

            sage: m=matrix(QQ,2,[1,0,1,1,1,0])
            sage: m._nonzero_positions_by_column()
            [(0, 0), (1, 0), (1, 1), (0, 2)]
        """
        x = self.fetch('nonzero_positions_by_column')
        if x is not None:
            if copy:
                return list(x)
            return x
        cdef Py_ssize_t i, j
        nzp = []
        for j from 0 <= j < self._ncols:
            for i from 0 <= i < self._nrows:
                if not self.get_is_zero_unsafe(i, j):
                    nzp.append((i, j))
        self.cache('nonzero_positions_by_column', nzp)
        if copy:
            return list(nzp)
        return nzp

    def nonzero_positions_in_column(self, Py_ssize_t i):
        """
        Return a sorted list of the integers ``j`` such that ``self[j,i]`` is
        nonzero, i.e., such that the ``j``-th position of the ``i``-th column
        is nonzero.

        INPUT:

        - ``i`` -- integer

        OUTPUT: list

        EXAMPLES::

            sage: a = matrix(QQ, 3,2, [1,2,0,2,0,0]); a
            [1 2]
            [0 2]
            [0 0]
            sage: a.nonzero_positions_in_column(0)
            [0]
            sage: a.nonzero_positions_in_column(1)
            [0, 1]

        You will get an :exc:`IndexError` if you select an invalid column::

            sage: a.nonzero_positions_in_column(2)
            Traceback (most recent call last):
            ...
            IndexError: matrix column index out of range
        """
        cdef Py_ssize_t j
        tmp = []

        if i < 0 or i >= self._ncols:
            raise IndexError("matrix column index out of range")
        for j from 0 <= j < self._nrows:
            if not self.get_is_zero_unsafe(j, i):
                tmp.append(j)
        return tmp

    def nonzero_positions_in_row(self, Py_ssize_t i):
        """
        Return the integers ``j`` such that ``self[i,j]`` is nonzero, i.e.,
        such that the ``j``-th position of the ``i``-th row is nonzero.

        INPUT:

        - ``i`` -- integer

        OUTPUT: list

        EXAMPLES::

            sage: a = matrix(QQ, 3,2, [1,2,0,2,0,0]); a
            [1 2]
            [0 2]
            [0 0]
            sage: a.nonzero_positions_in_row(0)
            [0, 1]
            sage: a.nonzero_positions_in_row(1)
            [1]
            sage: a.nonzero_positions_in_row(2)
            []
        """
        cdef Py_ssize_t j

        if i < 0 or i >= self._nrows:
            raise IndexError("matrix row index out of range")

        tmp = []

        for j from 0 <= j < self._ncols:
            if not self.get_is_zero_unsafe(i, j):
                tmp.append(j)
        return tmp

    def multiplicative_order(self):
        r"""
        Return the multiplicative order of this matrix, which must
        therefore be invertible.

        Only implemented over finite fields and over `\ZZ`.

        EXAMPLES:

        Over finite fields::

            sage: A = matrix(GF(59), 3, [10,56,39,53,56,33,58,24,55])
            sage: A.multiplicative_order()                                              # needs sage.libs.pari
            580
            sage: (A^580).is_one()
            True

            sage: B = matrix(GF(10007^3, 'b'), 0)                                       # needs sage.rings.finite_rings
            sage: B.multiplicative_order()                                              # needs sage.rings.finite_rings
            1

            sage: # needs sage.rings.finite_rings
            sage: M = MatrixSpace(GF(11^2, 'e'), 5)
            sage: E = M.random_element()
            sage: while E.det() == 0:
            ....:     E = M.random_element()
            sage: (E^E.multiplicative_order()).is_one()
            True

        Over `\ZZ`::

            sage: m = matrix(ZZ, 2, 2, [-1,1,-1,0])
            sage: m.multiplicative_order()                                              # needs sage.libs.pari
            3

            sage: m = posets.ChainPoset(6).coxeter_transformation()                     # needs sage.combinat sage.graphs
            sage: m.multiplicative_order()                                              # needs sage.combinat sage.graphs sage.groups
            7

            sage: P = posets.TamariLattice(4).coxeter_transformation()                  # needs sage.combinat sage.graphs
            sage: P.multiplicative_order()                                              # needs sage.combinat sage.graphs sage.groups
            10

            sage: M = matrix(ZZ, 2, 2, [1, 1, 0, 1])
            sage: M.multiplicative_order()                                              # needs sage.libs.pari
            +Infinity

            sage: for k in range(600):                                                  # needs sage.groups sage.modular
            ....:     m = SL2Z.random_element()
            ....:     o = m.multiplicative_order()
            ....:     if o != Infinity and m**o != SL2Z.one():
            ....:         raise RuntimeError

            sage: m24 = matrix.companion(cyclotomic_polynomial(24))
            sage: def val(i, j):
            ....:     if i < j:
            ....:         return 0
            ....:     elif i == j:
            ....:         return 1
            ....:     else:
            ....:         return ZZ.random_element(-100,100)
            sage: rnd = matrix(ZZ, 8, 8, val)
            sage: (rnd * m24 * rnd.inverse_of_unit()).multiplicative_order()            # needs sage.libs.pari
            24

        TESTS::

            sage: C = matrix(GF(2^10, 'c'), 2, 3, [1]*6)                                # needs sage.rings.finite_rings
            sage: C.multiplicative_order()                                              # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ArithmeticError: self must be invertible ...

            sage: D = matrix(IntegerModRing(6), 3, [5,5,3,0,2,5,5,4,0])
            sage: D.multiplicative_order()
            Traceback (most recent call last):
            ...
            NotImplementedError: ... only ... over finite fields or ZZ

        REFERENCES:

        - [CLG1997]_

        - [KP2002b]_
        """
        from sage.rings.integer import Integer
        from sage.rings.integer_ring import ZZ
        from sage.categories.fields import Fields

        n = self.ncols()
        if not n:
            return Integer(1)

        if not self.is_invertible():
            raise ArithmeticError("self must be invertible to have a multiplicative order")

        K = self.base_ring()

        if K in Fields().Finite():
            from sage.groups.generic import order_from_multiple
            P = self.minimal_polynomial()
            R = P.parent()
            P = P.factor()
            q = K.cardinality()
            p = K.characteristic()
            a = 0
            res = Integer(1)
            for f, m in P:
                a = max(a, m)
                S = R.quotient(f, 'y')
                res = res._lcm(order_from_multiple(S.gen(),
                                                   q**f.degree() - 1,
                                                   operation='*'))
            ppart = p**Integer(a).exact_log(p)
            if ppart < a:
                ppart *= p
            return res * ppart
        elif K is ZZ:
            from sage.rings.infinity import Infinity

            # two small odd prime numbers
            p1 = Integer(3)
            p2 = Integer(5)
            o1 = self.mod(p1).multiplicative_order()

            # Test if o1 cannot be the order of a matrix in GL_n(QQ)
            # Uses Thm 2.7 [KuPa2002]
            fac = o1.factor()
            S = sum((pi - 1) * pi**(ei - 1) for pi, ei in fac)
            if fac[0] == (2, 1):
                impossible_order = not(S <= n + 1)
            else:
                impossible_order = not(S <= n)
            if impossible_order:
                return Infinity

            o2 = self.mod(p2).multiplicative_order()
            if o1 != o2:
                return Infinity
            P = self.minimal_polynomial()
            x = P.parent().gen()
            if x**o1 % P == 1:  # or (x % P)**o1 == 1 ? maybe faster
                return o1
            else:
                return Infinity
        else:
            raise NotImplementedError("multiplicative order is only implemented"
                                      " for matrices over finite fields or ZZ")

    ###################################################
    # Arithmetic
    ###################################################
    cdef _vector_times_matrix_(self, Vector v):
        r"""
        Return the vector times matrix product.

        INPUT:

        - ``v`` -- a free module element

        OUTPUT: the vector times matrix product ``v*A``

        EXAMPLES::

            sage: B = matrix(QQ, 2, [1,2,3,4])
            sage: V = VectorSpace(QQ, 2)
            sage: v = V([-1,5])
            sage: v*B
            (14, 18)
            sage: -1*B.row(0) + 5*B.row(1)
            (14, 18)
            sage: B*v    # computes B*v, where v is interpreted as a column vector.
            (9, 17)
            sage: -1*B.column(0) + 5*B.column(1)
            (9, 17)

        We mix dense and sparse over different rings::

            sage: v = FreeModule(ZZ, 3, sparse=True)([1, 2, 3])
            sage: m = matrix(QQ, 3, 4, range(12))
            sage: v * m
            (32, 38, 44, 50)
            sage: v = FreeModule(ZZ, 3, sparse=False)([1, 2, 3])
            sage: m = matrix(QQ, 3, 4, range(12), sparse=True)
            sage: v * m
            (32, 38, 44, 50)
            sage: (v * m).parent() is m.row(0).parent()
            True

        TESTS:

        Check that :issue:`8198` is fixed::

            sage: # needs sage.rings.padics
            sage: R = Qp(5, 5)
            sage: x = R(5).add_bigoh(1)
            sage: I = matrix(R, [[1, 0], [0, 1]])
            sage: v = vector(R, [1, x])
            sage: v*I
            (1 + O(5^5), O(5))
        """
        M = self.row_ambient_module()
        if self._nrows != v._degree:
            raise ArithmeticError("number of rows of matrix must equal degree of vector")
        cdef Py_ssize_t i
        return sum([v[i] * self.row(i, from_list=True)
                    for i in range(self._nrows)], M(0))

    cdef _matrix_times_vector_(self, Vector v):
        """
        EXAMPLES::

            sage: v = FreeModule(ZZ, 3, sparse=True)([1, 2, 3])
            sage: m = matrix(QQ, 4, 3, range(12))
            sage: m * v
            (8, 26, 44, 62)
            sage: v = FreeModule(ZZ, 3, sparse=False)([1, 2, 3])
            sage: m = matrix(QQ, 4, 3, range(12), sparse=True)
            sage: m * v
            (8, 26, 44, 62)
            sage: (m * v).parent() is m.column(0).parent()
            True

        TESTS:

        Check that :issue:`8198` is fixed::

            sage: # needs sage.rings.padics
            sage: R = Qp(5, 5)
            sage: x = R(5).add_bigoh(1)
            sage: I = matrix(R, [[1, 0], [0, 1]])
            sage: v = vector(R, [1, x])
            sage: I*v
            (1 + O(5^5), O(5))
        """
        M = self.column_ambient_module()
        if self._ncols != v._degree:
            raise ArithmeticError("number of columns of matrix must equal degree of vector")
        cdef Py_ssize_t i
        return sum([self.column(i, from_list=True) * v[i]
                    for i in range(self._ncols)], M(0))

    def iterates(self, v, n, rows=True):
        r"""
        Let `A` be this matrix and `v` be a free module
        element. If rows is True, return a matrix whose rows are the
        entries of the following vectors:

        .. MATH::

                       v, v A, v A^2, \dots, v A^{n-1}.

        If rows is False, return a matrix whose columns are the entries of
        the following vectors:

        .. MATH::

                       v, Av, A^2 v, \dots, A^{n-1} v.

        INPUT:

        - ``v`` -- free module element

        - ``n`` -- nonnegative integer

        EXAMPLES::

            sage: A = matrix(ZZ, 2, [1,1,3,5]); A
            [1 1]
            [3 5]
            sage: v = vector([1,0])
            sage: A.iterates(v, 0)
            []
            sage: A.iterates(v, 5)
            [  1   0]
            [  1   1]
            [  4   6]
            [ 22  34]
            [124 192]

        Another example::

            sage: a = matrix(ZZ, 3, range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: v = vector([1,0,0])
            sage: a.iterates(v, 4)
            [  1   0   0]
            [  0   1   2]
            [ 15  18  21]
            [180 234 288]
            sage: a.iterates(v, 4, rows=False)
            [  1   0  15 180]
            [  0   3  42 558]
            [  0   6  69 936]
        """
        n = int(n)
        if n >= 2 and self.nrows() != self.ncols():
            raise ArithmeticError("matrix must be square if n >= 2.")
        if n == 0:
            return self.matrix_space(n, self.ncols())(0)
        m = self.nrows()
        M = sage.modules.free_module.FreeModule(self._base_ring, m, sparse=self.is_sparse())
        v = M(v)
        X = [v]

        if rows:
            for _ in range(n-1):
                X.append(X[len(X)-1]*self)
            MS = self.matrix_space(n, m)
            return MS(X)
        else:
            for _ in range(n-1):
                X.append(self*X[len(X)-1])
            MS = self.matrix_space(n, m)
            return MS(X).transpose()

    cpdef _add_(self, _right):
        """
        Add two matrices with the same parent.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: R.<x,y> = FreeAlgebra(QQ, 2)
            sage: a = matrix(2, 2, [1,2,x*y,y*x])
            sage: b = matrix(2, 2, [1,2,y*x,y*x])
            sage: a + b  # indirect doctest
            [        2         4]
            [x*y + y*x     2*y*x]
        """
        cdef Py_ssize_t i, j
        cdef Matrix A
        cdef Matrix right = _right
        A = self.new_matrix()
        for i in range(self._nrows):
            for j in range(self._ncols):
                A.set_unsafe(i, j, self.get_unsafe(i,j)._add_(right.get_unsafe(i,j)))
        return A

    cpdef _sub_(self, _right):
        """
        Subtract two matrices with the same parent.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: R.<x,y> = FreeAlgebra(QQ,2)
            sage: a = matrix(2, 2, [1,2,x*y,y*x])
            sage: b = matrix(2, 2, [1,2,y*x,y*x])
            sage: a - b  # indirect doctest
            [        0         0]
            [x*y - y*x         0]
        """
        cdef Py_ssize_t i, j
        cdef Matrix A
        cdef Matrix right = _right
        A = self.new_matrix()
        for i in range(self._nrows):
            for j in range(self._ncols):
                A.set_unsafe(i,j,self.get_unsafe(i,j)._sub_(right.get_unsafe(i,j)))
        return A

    def __mod__(self, p):
        r"""
        Return matrix mod `p`, returning again a matrix over the
        same base ring.

        .. NOTE::

           Use :meth:`mod` to obtain a matrix over the residue class ring
           modulo `p`.

        EXAMPLES::

            sage: M = Matrix(ZZ, 2, 2, [5, 9, 13, 15])
            sage: M % 7
            [5 2]
            [6 1]
            sage: parent(M % 7)
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        """
        cdef Py_ssize_t i, j
        cdef Matrix s = self
        cdef Matrix A = s.new_matrix()
        for i in range(A._nrows):
            for j in range(A._ncols):
                A[i,j] = s.get_unsafe(i,j) % p
        return A

    def mod(self, p):
        """
        Return matrix mod `p`, over the reduced ring.

        EXAMPLES::

            sage: M = matrix(ZZ, 2, 2, [5, 9, 13, 15])
            sage: M.mod(7)
            [5 2]
            [6 1]
            sage: parent(M.mod(7))
            Full MatrixSpace of 2 by 2 dense matrices over Ring of integers modulo 7
        """
        return self.change_ring(self._base_ring.quotient_ring(p))

    cpdef _rmul_(self, Element left):
        """
        EXAMPLES::

            sage: a = matrix(QQ['x'], 2, range(6))
            sage: (3/4) * a
            [   0  3/4  3/2]
            [ 9/4    3 15/4]

            sage: R.<x,y> = QQ[]
            sage: a = matrix(R, 2, 3, [1,x,y, -x*y,x+y,x-y]); a
            [    1     x     y]
            [ -x*y x + y x - y]
            sage: (x*y) * a
            [          x*y         x^2*y         x*y^2]
            [     -x^2*y^2 x^2*y + x*y^2 x^2*y - x*y^2]

            sage: # needs sage.combinat
            sage: R.<x,y> = FreeAlgebra(ZZ,2)
            sage: a = matrix(R, 2, 3, [1,x,y, -x*y,x+y,x-y]); a
            [    1     x     y]
            [ -x*y x + y x - y]
            sage: (x*y) * a  # indirect doctest
            [          x*y         x*y*x         x*y^2]
            [     -x*y*x*y x*y*x + x*y^2 x*y*x - x*y^2]
        """
        # derived classes over a commutative base *just* overload _lmul_ (!!)
        if self._base_ring in CommutativeRings():
            return self._lmul_(left)
        cdef Py_ssize_t r,c
        x = self._base_ring(left)
        cdef Matrix ans
        ans = self._parent.zero_matrix().__copy__()
        for r from 0 <= r < self._nrows:
            for c from 0 <= c < self._ncols:
                ans.set_unsafe(r, c, x * self.get_unsafe(r, c))
        return ans

    cpdef _lmul_(self, Element right):
        """
        EXAMPLES:

        A simple example in which the base ring is commutative::

            sage: a = matrix(QQ['x'],2,range(6))
            sage: a*(3/4)
            [   0  3/4  3/2]
            [ 9/4    3 15/4]

        An example in which the base ring is not commutative::

            sage: # needs sage.combinat
            sage: F.<x,y> = FreeAlgebra(QQ,2)
            sage: a = matrix(2, [x,y, x^2,y^2]); a
            [  x   y]
            [x^2 y^2]
            sage: x * a # indirect doctest
            [  x^2   x*y]
            [  x^3 x*y^2]
            sage: a * y
            [  x*y   y^2]
            [x^2*y   y^3]

            sage: # needs sage.combinat
            sage: R.<x,y> = FreeAlgebra(ZZ,2)
            sage: a = matrix(R, 2, 3, [1,x,y, -x*y,x+y,x-y]); a
            [    1     x     y]
            [ -x*y x + y x - y]
            sage: a * (x*y)
            [          x*y         x^2*y         y*x*y]
            [     -x*y*x*y x^2*y + y*x*y x^2*y - y*x*y]
        """
        # derived classes over a commutative base *just* overload this and not _rmul_
        cdef Py_ssize_t r,c
        x = self._base_ring(right)
        cdef Matrix ans
        ans = self._parent.zero_matrix().__copy__()
        for r from 0 <= r < self._nrows:
            for c from 0 <= c < self._ncols:
                ans.set_unsafe(r, c, self.get_unsafe(r, c) * x)
        return ans

    cdef sage.structure.element.Matrix _matrix_times_matrix_(self, sage.structure.element.Matrix right):
        r"""
        Return the product of two matrices.

        EXAMPLE of matrix times matrix over same base ring: We multiply
        matrices over `\QQ[x,y]`.

        ::

            sage: R.<x,y> = QQ[]
            sage: a = matrix(R, 2, 3, [1,x,y, -x*y,x+y,x-y]); a
            [    1     x     y]
            [ -x*y x + y x - y]
            sage: b = a.transpose(); b
            [    1  -x*y]
            [    x x + y]
            [    y x - y]
            sage: a*b  # indirect doctest
            [          x^2 + y^2 + 1         x^2 + x*y - y^2]
            [        x^2 + x*y - y^2 x^2*y^2 + 2*x^2 + 2*y^2]
            sage: b*a  # indirect doctest
            [        x^2*y^2 + 1  -x^2*y - x*y^2 + x  -x^2*y + x*y^2 + y]
            [ -x^2*y - x*y^2 + x 2*x^2 + 2*x*y + y^2     x^2 + x*y - y^2]
            [ -x^2*y + x*y^2 + y     x^2 + x*y - y^2 x^2 - 2*x*y + 2*y^2]

        We verify that the matrix multiplies are correct by comparing them
        with what PARI gets::

            sage: gp(a)*gp(b) - gp(a*b)                                                 # needs sage.libs.pari
            [0, 0; 0, 0]
            sage: gp(b)*gp(a) - gp(b*a)                                                 # needs sage.libs.pari
            [0, 0, 0; 0, 0, 0; 0, 0, 0]

        EXAMPLE of matrix times matrix over different base rings::

            sage: a = matrix(ZZ, 2, 2, range(4))
            sage: b = matrix(GF(7), 2, 2, range(4))
            sage: c = matrix(QQ, 2, 2, range(4))
            sage: d = a * b; d
            [2 3]
            [6 4]
            sage: parent(d)
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7
            sage: parent(b * a)
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7
            sage: d = a * c; d
            [ 2  3]
            [ 6 11]
            sage: parent(d)
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: d = b + c
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for +:
             'Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7' and
             'Full MatrixSpace of 2 by 2 dense matrices over Rational Field'
            sage: d = b + c.change_ring(GF(7)); d
            [0 2]
            [4 6]

        EXAMPLE of matrix times matrix where one matrix is sparse and the
        other is dense (in such mixed cases, the result is always dense)::

            sage: a = matrix(ZZ, 2, 2, range(4), sparse=True)
            sage: b = matrix(GF(7), 2, 2, range(4), sparse=False)
            sage: c = a * b; c
            [2 3]
            [6 4]
            sage: parent(c)
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7
            sage: c = b * a; c
            [2 3]
            [6 4]
            sage: parent(c)
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7

        EXAMPLE of matrix multiplication over a noncommutative base ring::

            sage: # needs sage.combinat
            sage: R.<x,y> = FreeAlgebra(QQ, 2)
            sage: x*y - y*x
            x*y - y*x
            sage: a = matrix(2, 2, [1,2, x,y])
            sage: b = matrix(2, 2, [x,y, x^2,y^2])
            sage: a*b
            [  x + 2*x^2   y + 2*y^2]
            [x^2 + y*x^2   x*y + y^3]
            sage: b*a
            [    x + y*x   2*x + y^2]
            [x^2 + y^2*x 2*x^2 + y^3]

        EXAMPLE of row vector times matrix (vectors are row vectors, so
        matrices act from the right)::

            sage: a = matrix(2,3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: V = ZZ^2
            sage: v = V([-2,3]); v
            (-2, 3)
            sage: v*a
            (9, 10, 11)

        This is not allowed, since v is a *row* vector::

            sage: a*v
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *:
            'Full MatrixSpace of 2 by 3 dense matrices over Integer Ring' and
            'Ambient free module of rank 2 over the principal ideal domain Integer Ring'

        This illustrates how coercion works::

            sage: V = QQ^2
            sage: v = V([-2,3]); v
            (-2, 3)
            sage: parent(v*a)
            Vector space of dimension 3 over Rational Field

        EXAMPLE of matrix times column vector: (column vectors are not
        implemented yet) TODO TODO

        EXAMPLE of scalar times matrix::

            sage: a = matrix(2,3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: b = 3*a; b
            [ 0  3  6]
            [ 9 12 15]
            sage: parent(b)
            Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
            sage: b = (2/3)*a; b
            [   0  2/3  4/3]
            [   2  8/3 10/3]
            sage: parent(b)
            Full MatrixSpace of 2 by 3 dense matrices over Rational Field

        EXAMPLE of matrix times scalar::

            sage: a = matrix(2,3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: b = a*3; b
            [ 0  3  6]
            [ 9 12 15]
            sage: parent(b)
            Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
            sage: b = a*(2/3); b
            [   0  2/3  4/3]
            [   2  8/3 10/3]
            sage: parent(b)
            Full MatrixSpace of 2 by 3 dense matrices over Rational Field

        EXAMPLE of scalar multiplication in the noncommutative case::

            sage: # needs sage.combinat
            sage: R.<x,y> = FreeAlgebra(ZZ, 2)
            sage: a = matrix(2, [x,y, x^2,y^2])
            sage: a * x
            [  x^2   y*x]
            [  x^3 y^2*x]
            sage: x * a
            [  x^2   x*y]
            [  x^3 x*y^2]
            sage: a*x - x*a
            [             0     -x*y + y*x]
            [             0 -x*y^2 + y^2*x]
        """
        # Both self and right are matrices with compatible dimensions and base ring.
        if self._will_use_strassen(right):
            return self._multiply_strassen(right)
        else:
            return self._multiply_classical(right)

    cdef bint _will_use_strassen(self, Matrix right) except -2:
        """
        Whether or not matrix multiplication of ``self`` by ``right`` should be
        done using Strassen.

        Overload _strassen_default_cutoff to return -1 to not use
        Strassen by default.
        """
        cdef int n
        n = self._strassen_default_cutoff(right)
        if n == -1:
            return 0  # do not use Strassen
        if self._nrows > n and self._ncols > n and \
               right._nrows > n and right._ncols > n:
            return 1
        return 0

    cdef bint _will_use_strassen_echelon(self) except -2:
        """
        Whether or not matrix multiplication of ``self`` by ``right`` should be
        done using Strassen.

        Overload this in derived classes to not use Strassen by default.
        """
        cdef int n
        n = self._strassen_default_echelon_cutoff()
        if n == -1:
            return 0  # do not use Strassen
        if self._nrows > n and self._ncols > n:
            return 1
        return 0

    def __neg__(self):
        """
        Return the negative of ``self``.

        EXAMPLES::

            sage: a = matrix(ZZ,2,range(4))
            sage: a.__neg__()
            [ 0 -1]
            [-2 -3]
            sage: -a
            [ 0 -1]
            [-2 -3]
        """
        return self._lmul_(self._base_ring(-1))

    def __invert__(self):
        r"""
        Return the inverse of this matrix, as a matrix over the fraction
        field.

        Raises a :exc:`ZeroDivisionError` if the matrix has zero
        determinant, and raises an :exc:`ArithmeticError`, if the
        inverse doesn't exist because the matrix is nonsquare. Also, note,
        e.g., that the inverse of a matrix over `\ZZ` is
        always a matrix defined over `\QQ` (even if the
        entries are integers).

        EXAMPLES::

            sage: A = MatrixSpace(ZZ, 2)([1,1,3,5])
            sage: ~A
            [ 5/2 -1/2]
            [-3/2  1/2]
            sage: A.__invert__()
            [ 5/2 -1/2]
            [-3/2  1/2]

        Even if the inverse lies in the base field, the result is still a
        matrix over the fraction field.

        ::

            sage: I = MatrixSpace(ZZ, 2)(1)  # identity matrix
            sage: ~I
            [1 0]
            [0 1]
            sage: (~I).parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

        This is analogous to the situation for ring elements, e.g., for
        `\ZZ` we have::

            sage: parent(~1)
            Rational Field

        A matrix with 0 rows and 0 columns is invertible (see :issue:`3734`)::

            sage: M = MatrixSpace(RR, 0, 0)(0); M
            []
            sage: M.determinant()
            1.00000000000000
            sage: M.is_invertible()
            True
            sage: M.inverse() == M
            True

        Matrices over the integers modulo a composite modulus::

            sage: m = matrix(Zmod(49), 2, [2,1,3,3])
            sage: type(m)
            <class 'sage.matrix.matrix_modn_dense_float.Matrix_modn_dense_float'>
            sage: ~m
            [ 1 16]
            [48 17]
            sage: m = matrix(Zmod(2^100), 2, [2,1,3,3])
            sage: type(m)
            <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
            sage: (~m)*m                                                                # needs sage.libs.pari
            [1 0]
            [0 1]
            sage: ~m                                                                    # needs sage.libs.pari
            [                              1  422550200076076467165567735125]
            [1267650600228229401496703205375  422550200076076467165567735126]

        Matrices over `p`-adics. See :issue:`17272` ::

            sage: # needs sage.rings.padics
            sage: R = ZpCA(5, 5, print_mode='val-unit')
            sage: A = matrix(R, 3, 3, [250,2369,1147,106,927,362,90,398,2483])
            sage: A
            [5^3 * 2 + O(5^5)    2369 + O(5^5)    1147 + O(5^5)]
            [    106 + O(5^5)     927 + O(5^5)     362 + O(5^5)]
            [ 5 * 18 + O(5^5)     398 + O(5^5)    2483 + O(5^5)]
            sage: ~A
            [5 * 212 + O(5^5)    3031 + O(5^5)    2201 + O(5^5)]
            [   1348 + O(5^5) 5 * 306 + O(5^5)    2648 + O(5^5)]
            [   1987 + O(5^5) 5 * 263 + O(5^5)     154 + O(5^5)]

        This matrix is not invertible::

            sage: m = matrix(Zmod(9), 2, [2,1,3,3])
            sage: ~m
            Traceback (most recent call last):
            ...
            ZeroDivisionError: input matrix must be nonsingular

        Check to make sure that :issue:`2256` is still fixed::

            sage: M = MatrixSpace(CC, 2)(-1.10220440881763)
            sage: N = ~M
            sage: (N*M).norm()
            0.9999999999999999

        Check that :issue:`28402` is fixed::

            sage: B = matrix(RR, [[1/6, -1/24, -1/30, 1/120,1/12, 0, 0, 0, 0],
            ....:                 [-1/24,1/60,1/60, 1/420, -1/24, 0, 0, 0, 0],
            ....:                 [-1/30,1/60, 2/105, 1/140, -1/20, 0, 0, 0, 0],
            ....:                 [1/120, 1/420, 1/140, 13/1260, -1/40, 0, 0, 0, 0],
            ....:                 [1/12, -1/24, -1/20, -1/40, 1/3, -1/24, -1/30, 1/120,1/12],
            ....:                 [0, 0, 0, 0, -1/24,1/60,1/60, 1/420, -1/24],
            ....:                 [0, 0, 0, 0, -1/30,1/60, 2/105, 1/140, -1/20],
            ....:                 [0, 0, 0, 0, 1/120, 1/420, 1/140, 13/1260, -1/40],
            ....:                 [0, 0, 0, 0,1/12, -1/24, -1/20, -1/40, 1/6]],
            ....:           sparse=True)
            sage: (B.inverse()*B).norm(1)  # rel tol 2e-12
            1.0
            sage: B = matrix(QQ, [[1/6, -1/24, -1/30, 1/120,1/12, 0, 0, 0, 0],
            ....:                 [-1/24,1/60,1/60, 1/420, -1/24, 0, 0, 0, 0],
            ....:                 [-1/30,1/60, 2/105, 1/140, -1/20, 0, 0, 0, 0],
            ....:                 [1/120, 1/420, 1/140, 13/1260, -1/40, 0, 0, 0, 0],
            ....:                 [1/12, -1/24, -1/20, -1/40, 1/3, -1/24, -1/30, 1/120,1/12],
            ....:                 [0, 0, 0, 0, -1/24,1/60,1/60, 1/420, -1/24],
            ....:                 [0, 0, 0, 0, -1/30,1/60, 2/105, 1/140, -1/20],
            ....:                 [0, 0, 0, 0, 1/120, 1/420, 1/140, 13/1260, -1/40],
            ....:                 [0, 0, 0, 0,1/12, -1/24, -1/20, -1/40, 1/6]],
            ....:           sparse=True)
            sage: (B.inverse()*B).norm(1)
            1.0
        """
        if not self.is_square():
            raise ArithmeticError("self must be a square matrix")
        if not self.nrows():
            return self

        R = self.base_ring()
        if R not in _Fields:
            if R in _IntegralDomains:
                return ~self.matrix_over_field()
            else:
                return self.inverse_of_unit()
        else:
            A = self.augment(self.parent().identity_matrix())
            A.echelonize()

            # Now we want to make sure that B is of the form [I|X], in
            # which case X is the inverse of self. We can simply look at
            # the lower right entry of the left half of B, and make sure
            # that it's 1.
            #
            # However, doing this naively causes trouble over inexact
            # fields -- see Issue #2256. The *right* thing to do would
            # probably be to make sure that self.det() is nonzero. That
            # doesn't work here, because our det over an arbitrary field
            # just does expansion by minors and is unusable for even 10x10
            # matrices over CC. Instead, we choose a different band-aid:
            # we check to make sure that the lower right entry isn't
            # 0. Since we're over a field, we know that it *should* be
            # either 1 or 0. This can still cause trouble, but it's
            # significantly better than it was before.
            #
            # Over exact rings, of course, we still want the old
            # behavior.

            if R.is_exact():
                if not A[self._nrows-1, self._ncols-1].is_one():
                    raise ZeroDivisionError("input matrix must be nonsingular")
                if self.is_sparse():
                    return self.build_inverse_from_augmented_sparse(A)
            else:
                if not A[self._nrows-1, self._ncols-1]:
                    raise ZeroDivisionError("input matrix must be nonsingular")
            return A.matrix_from_columns(list(range(self._ncols, 2 * self._ncols)))

    cdef build_inverse_from_augmented_sparse(self, A):
        # We can directly use the dict entries of A
        cdef Py_ssize_t i, nrows
        cdef dict data = <dict> A._dict()
        nrows = self._nrows
        # We can modify data because A is local to this function
        for i in range(nrows):
            del data[i,i]
        data = {(r,c-nrows): data[r,c] for (r,c) in data}
        return self._parent(data)

    def inverse_of_unit(self, algorithm=None):
        r"""
        Return the inverse of this matrix in the same matrix space.

        The matrix must be invertible on the base ring. Otherwise, an
        :exc:`ArithmeticError` is raised.

        The computation goes through the matrix of cofactors and avoids
        division. In particular the base ring does not need to have a
        fraction field.

        INPUT:

        - ``algorithm`` -- (default: ``None``) either ``None`` or ``'df'`` (for
          division free)

        EXAMPLES::

            sage: R.<a,b,c,d> = ZZ[]
            sage: RR = R.quotient(a*d - b*c - 1)
            sage: a,b,c,d = RR.gens()                                                   # needs sage.libs.singular
            sage: m = matrix(2, [a,b, c,d])
            sage: n = m.inverse_of_unit()                                               # needs sage.libs.singular
            sage: m * n                                                                 # needs sage.libs.singular
            [1 0]
            [0 1]

            sage: matrix(RR, 2, 1, [a,b]).inverse_of_unit()                             # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            ArithmeticError: self must be a square matrix

            sage: matrix(RR, 1, 1, [2]).inverse_of_unit()                               # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            ArithmeticError: non-invertible matrix

            sage: R = ZZ.cartesian_product(ZZ)
            sage: m = matrix(R, 2, [R((2,1)), R((1,1)), R((1,1)), R((1,2))])
            sage: m * m.inverse_of_unit()
            [(1, 1) (0, 0)]
            [(0, 0) (1, 1)]

        Tests for :issue:`28570`::

            sage: P = posets.TamariLattice(7)                                           # needs sage.graphs
            sage: M = P._hasse_diagram._leq_matrix                                      # needs sage.graphs
            sage: M.inverse_of_unit()   # this was very slow, now 1s                    # needs sage.graphs
            429 x 429 sparse matrix over Integer Ring...

            sage: m = matrix(Zmod(2**2), 1, 1, [1], sparse=True)
            sage: mi = ~m; mi
            [1]
            sage: mi.parent()
            Full MatrixSpace of 1 by 1 sparse matrices over Ring of integers modulo 4
        """
        n = self.nrows()
        if n != self.ncols():
            raise ArithmeticError("self must be a square matrix")

        R = self.base_ring()
        if algorithm is None and R in _Fields:
            return ~self
        elif algorithm is None and isinstance(R, sage.rings.abc.IntegerModRing):
            # Finite fields are handled above.
            # This is "easy" in that we either get an error or
            # the right answer. Note that of course there
            # could be a much faster algorithm, e.g., using
            # CRT or p-adic lifting.
            try:
                return (~self.lift_centered()).change_ring(R)
            except (TypeError, ZeroDivisionError):
                raise ZeroDivisionError("input matrix must be nonsingular")
        elif algorithm is None and isinstance(R, IntegerRing_class):
            try:
                return (~self).change_ring(R)
            except (TypeError, ZeroDivisionError):
                raise ZeroDivisionError("input matrix must be nonsingular")
        elif algorithm is None or algorithm == "df":
            d = self.det()
            if d.is_one():
                dinv = d
            elif not d.is_unit():
                raise ArithmeticError("non-invertible matrix")
            else:
                dinv = d.inverse_of_unit()

            return dinv * self.adjugate()
        else:
            raise ValueError('algorithm can only be "df"')

    def __pos__(self):
        """
        Return +self, which is just self, of course.

        EXAMPLES::

            sage: a = matrix(ZZ,2,range(4))
            sage: +a
            [0 1]
            [2 3]
            sage: a.__pos__()
            [0 1]
            [2 3]
        """
        return self

    def __pow__(self, n, ignored):
        """
        EXAMPLES::

            sage: MS = MatrixSpace(QQ, 3, 3)
            sage: A = MS([0, 0, 1, 1, 0, '-2/11', 0, 1, '-3/11'])
            sage: A * A^(-1) == 1
            True
            sage: A^4
            [      -3/11     -13/121   1436/1331]
            [    127/121   -337/1331 -4445/14641]
            [    -13/121   1436/1331 -8015/14641]
            sage: A.__pow__(4)
            [      -3/11     -13/121   1436/1331]
            [    127/121   -337/1331 -4445/14641]
            [    -13/121   1436/1331 -8015/14641]

        Sage follows Python's convention 0^0 = 1, as each of the following
        examples show::

            sage: a = Matrix([[1,0],[0,0]]); a
            [1 0]
            [0 0]
            sage: a^0 # lower right entry is 0^0
            [1 0]
            [0 1]
            sage: Matrix([[0]])^0
            [1]
            sage: 0^0
            1

        Non-integer (symbolic) exponents are also supported::

            sage: k = var('k')                                                          # needs sage.symbolic
            sage: A = matrix([[2, -1], [1,  0]])
            sage: A^(2*k+1)                                                             # needs sage.symbolic
            [ 2*k + 2 -2*k - 1]
            [ 2*k + 1     -2*k]
        """
        from sage.structure.element import Expression

        if not self.is_square():
            raise ArithmeticError("self must be a square matrix")
        if ignored is not None:
            raise RuntimeError("__pow__ third argument not used")
        if isinstance(n, Expression):
            from sage.matrix.matrix2 import _matrix_power_symbolic
            return _matrix_power_symbolic(self, n)
        return generic_power(self, n)

    ###################################################
    # Comparison
    ###################################################
    def __hash__(self):
        """
        Return the hash of this matrix.

        Equal matrices should have equal hashes, even if one is sparse
        and the other is dense. We also ensure that zero matrices hash
        to zero and that scalar matrices have the same hash as the
        scalar.

        EXAMPLES::

            sage: m = matrix(2, range(24), sparse=True)
            sage: m.set_immutable()
            sage: hash(m)
            3327233128576517516  # 64-bit
            -373881460           # 32-bit

        ::

            sage: d = m.dense_matrix()
            sage: d.set_immutable()
            sage: hash(m) == hash(d)
            True

        ::

            sage: R.<x> = ZZ[]
            sage: M = matrix(R, 10, 20); M.set_immutable()
            sage: hash(M)
            0
            sage: M = matrix(R, 10, 10, x); M.set_immutable()
            sage: hash(M) == hash(x)
            True
        """
        if not self._is_immutable:
            raise TypeError("mutable matrices are unhashable")
        if self.hash != -1:
            return self.hash
        cdef long h = self._hash_()
        self.hash = h
        return h

    cdef long _hash_(self) except -1:
        """
        Implementation of hash function.

        AUTHOR: Jeroen Demeyer
        """
        cdef long C[5]
        self.get_hash_constants(C)

        # The hash is of the form
        #
        #     sum_{i,j} F(i,j) * hash(M[i,j])
        #
        # The fact that it is a sum means that it can be computed in
        # any order, which is useful for sparse matrices or other
        # matrix implementations.
        #
        # Entries which have zero hash do not contribute to the matrix
        # hash. This is again useful for sparse matrices, where we can
        # safely skip the zero entries (assuming that the hash of the
        # zero element is zero, which should be the case)
        #
        # To get a predictable hash for scalar matrices, some tricks
        # are needed. First of all, we compute F(i,j) as k xor l, where
        # l is zero for diagonal entries and where k (which depends only
        # on the row) is called the row multiplier.
        #
        # So the hash of a scalar matrix is the sum of all row
        # multipliers times the hash of the scalar. Therefore, the hash
        # of a scalar matrix equals the hash of the scalar if the sum of
        # all row multipliers is 1. We choose the constants (see
        # get_hash_constants()) such that this in indeed the case.
        # Actually, this is not the complete story: we additionally
        # multiply all row constants with some random value C[2] and
        # then at the end we multiply with C[2]^-1 = C[4].
        # This gives better mixing.
        #
        # The value for l in the loop below is not so important: it
        # must be zero if i == j and sufficiently complicated to avoid
        # hash collisions.
        cdef long h = 0, k, l
        cdef Py_ssize_t i, j
        for i in range(self._nrows):
            k = C[0] if i == 0 else C[1] + C[2] * i
            for j in range(self._ncols):
                sig_check()
                l = C[3] * (i - j) * (i ^ j)
                h += (k ^ l) * hash(self.get_unsafe(i, j))
        h *= C[4]

        if h == -1:
            return -2
        return h

    cdef void get_hash_constants(self, long C[5]) noexcept:
        """
        Get constants for the hash algorithm.
        """
        cdef long m = self._nrows
        cdef long n = self._ncols

        # XKCD-221 compliant random numbers
        C[1] = 0x6951766c055d2c0a
        C[2] = 0x1155b61baeb88b61  # must be odd
        C[3] = 0x0d58d3c0539376c1  # should be odd

        # Multiplicative inverse of C[2] mod 2^64
        C[4] = 0x7c7067f7da6758a1

        # Change some of these constants such that matrices with the
        # same entries but a different size have different hashes.
        # C[2] must be the same for all square matrices though.
        C[1] *= n + C[1] * m
        C[2] += (n - m) * ((C[3] * m) ^ n)

        # The k in the hashing loop is called the row multiplier. For
        # the i-th row with i > 0, this is (C[1] + C[2]*i). We choose
        # C[0] (the row multiplier for the 0-th row) such that the sum
        # of all row multipliers is C[2].
        #
        # This way, the row multiplier is never a small number in
        # absolute value and the row multipliers depend on the size of
        # the matrix.

        # mm = m * (m - 1)/2 computed correctly mod 2^wordsize.
        cdef long mm = (m // 2) * ((m - 1) | 1)

        # C[0] = (1 - m * (m - 1)/2) * C[2] - (m - 1) * C[1]
        C[0] = (1 - mm) * C[2] - (m - 1) * C[1]

    cpdef _richcmp_(left, right, int op):
        """
        Compare two matrices.

        Matrices are compared in lexicographic order on the underlying list
        of coefficients. A dense matrix and a sparse matrix are equal if
        their coefficients are the same.

        EXAMPLES: EXAMPLE comparing sparse and dense matrices::

            sage: matrix(QQ,2,range(4)) == matrix(QQ,2,range(4),sparse=True)
            True
            sage: matrix(QQ,2,range(4)) == matrix(QQ,2,range(4),sparse=True)
            True

        Dictionary order::

            sage: matrix(ZZ,2,[1,2,3,4]) < matrix(ZZ,2,[3,2,3,4])
            True
            sage: matrix(ZZ,2,[1,2,3,4]) > matrix(ZZ,2,[3,2,3,4])
            False
            sage: matrix(ZZ,2,[0,2,3,4]) < matrix(ZZ,2,[0,3,3,4], sparse=True)
            True
        """
        raise NotImplementedError  # this is defined in the derived classes

    def __bool__(self):
        """
        EXAMPLES::

            sage: M = Matrix(ZZ, 2, 2, [0, 0, 0, 0])
            sage: bool(M)
            False
            sage: M = Matrix(ZZ, 2, 2, [1, 2, 3, 5])
            sage: bool(M)
            True
        """
        cdef Py_ssize_t i, j
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                if not self.get_is_zero_unsafe(i,j):
                    return True
        return False

    cdef int _strassen_default_cutoff(self, Matrix right) except -2:
        return -1

    cdef int _strassen_default_echelon_cutoff(self) except -2:
        return -1


#######################
# Unpickling
#######################

def unpickle(cls, parent, immutability, cache, data, version):
    r"""
    Unpickle a matrix. This is only used internally by Sage. Users
    should never call this function directly.

    EXAMPLES: We illustrating saving and loading several different
    types of matrices.

    OVER `\ZZ`::

        sage: A = matrix(ZZ, 2, range(4))
        sage: loads(dumps(A))  # indirect doctest
        [0 1]
        [2 3]

    Sparse OVER `\QQ`:

    Dense over `\QQ[x,y]`:

    Dense over finite field.
    """
    cdef Matrix A
    A = cls.__new__(cls, parent, 0,0,0)
    A._parent = parent  # make sure -- __new__ doesn't have to set it, but unpickle may need to know.
    A._nrows = parent.nrows()
    A._ncols = parent.ncols()
    A._is_immutable = immutability
    A._base_ring = parent.base_ring()
    A._cache = cache
    if version >= 0:
        A._unpickle(data, version)
    else:
        A._unpickle_generic(data, version)
    return A


def set_max_rows(n):
    """
    Set the global variable ``max_rows`` (which is used in deciding how to
    output a matrix).

    EXAMPLES::

        sage: from sage.matrix.matrix0 import set_max_rows
        sage: set_max_rows(20)
        doctest:...: DeprecationWarning: 'set_max_rows' is replaced by 'matrix.options.max_rows'
        See https://github.com/sagemath/sage/issues/30552 for details.
    """
    from sage.misc.superseded import deprecation
    deprecation(30552, "'set_max_rows' is replaced by 'matrix.options.max_rows'")
    from sage.matrix.constructor import options
    options.max_rows = n-1


def set_max_cols(n):
    """
    Set the global variable ``max_cols`` (which is used in deciding how to
    output a matrix).

    EXAMPLES::

        sage: from sage.matrix.matrix0 import set_max_cols
        sage: set_max_cols(50)
        doctest:...: DeprecationWarning: 'set_max_cols' is replaced by 'matrix.options.max_cols'
        See https://github.com/sagemath/sage/issues/30552 for details.
    """
    from sage.misc.superseded import deprecation
    deprecation(30552, "'set_max_cols' is replaced by 'matrix.options.max_cols'")
    from sage.matrix.constructor import options
    options.max_cols = n-1
