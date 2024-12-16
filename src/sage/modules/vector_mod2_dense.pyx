# distutils: libraries = M4RI_LIBRARIES GDLIB_LIBRARIES LIBPNG_LIBRARIES
# distutils: library_dirs = M4RI_LIBDIR GDLIB_LIBDIR LIBPNG_LIBDIR
# distutils: include_dirs = M4RI_INCDIR GDLIB_INCDIR LIBPNG_INCDIR
# distutils: extra_compile_args = M4RI_CFLAGS

r"""
Vectors with elements in `\GF{2}`

AUTHOR:

- Martin Albrecht (2009-12): initial implementation
- Thomas Feulner (2012-11): added :meth:`Vector_mod2_dense.hamming_weight`

EXAMPLES::

    sage: VS = GF(2)^3
    sage: e = VS.random_element()
    sage: e.parent() is VS
    True
    sage: S = set(vector(v, immutable=True) for v in VS)
    sage: S1 = set()
    sage: while S != S1:
    ....:     S1.add(vector(VS.random_element(), immutable=True))

TESTS::

    sage: w = vector(GF(2), [-1,0,0,0])
    sage: w.set_immutable()
    sage: isinstance(hash(w), int)
    True
"""

# ****************************************************************************
#       Copyright (C) 2009 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.finite_rings.integer_mod cimport IntegerMod_int, IntegerMod_abstract
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.structure.element cimport Element, Vector
from sage.structure.richcmp cimport rich_to_bool
cimport sage.modules.free_module_element as free_module_element
from libc.stdint cimport uintptr_t

from sage.libs.m4ri cimport mzd_add, mzd_copy, mzd_cmp, mzd_free, mzd_init, mzd_set_ui, mzd_read_bit, mzd_row, mzd_write_bit, m4ri_word

cdef class Vector_mod2_dense(free_module_element.FreeModuleElement):
    cdef _new_c(self):
        """
        EXAMPLES::

            sage: VS = VectorSpace(GF(2),3)
            sage: VS([0,0,1])
            (0, 0, 1)
            sage: type(_)
            <class 'sage.modules.vector_mod2_dense.Vector_mod2_dense'>
        """
        cdef Vector_mod2_dense y
        y = Vector_mod2_dense.__new__(Vector_mod2_dense)
        y._init(self._degree, self._parent)
        return y

    cdef bint is_dense_c(self) noexcept:
        """
        EXAMPLES::

            sage: VS = VectorSpace(GF(2),3)
            sage: VS([0,0,1]).is_dense()
            True
        """
        return 1

    cdef bint is_sparse_c(self) noexcept:
        """
        EXAMPLES::

            sage: VS = VectorSpace(GF(2),3)
            sage: VS([0,0,1]).is_sparse()
            False
        """
        return 0

    def __copy__(self):
        """
        EXAMPLES::

            sage: VS = VectorSpace(GF(2),10^4)
            sage: v = VS.random_element()
            sage: w = copy(v)
            sage: w == v
            True
            sage: v[:10] == w[:10]
            True
            sage: v[5] += 1
            sage: v == w
            False
        """
        cdef Vector_mod2_dense y = self._new_c()
        if self._degree:
            mzd_copy(y._entries, self._entries)
        return y

    cdef _init(self, Py_ssize_t degree, parent):
        """
        EXAMPLES::

            sage: VS = VectorSpace(GF(2),3)
            sage: VS([0,0,1])
            (0, 0, 1)
            sage: type(_)
            <class 'sage.modules.vector_mod2_dense.Vector_mod2_dense'>
        """
        self._degree = degree
        self._parent = parent
        self._base_ring = parent.base_ring()
        self._entries = mzd_init(1, degree)
        if self._entries == NULL:
            raise MemoryError("Allocation of Vector_mod2_dense failed.")

    def __cinit__(self, parent=None, x=None, coerce=True, copy=True):
        """
        EXAMPLES::

            sage: VS = VectorSpace(GF(2),3)
            sage: VS((0,0,1/3))
            (0, 0, 1)
            sage: type(_)
            <class 'sage.modules.vector_mod2_dense.Vector_mod2_dense'>
        """
        self._entries = NULL
        self._is_immutable = 0
        if parent is not None:
            self._init(parent.degree(), parent)

    def __init__(self, parent, x, coerce=True, copy=True):
        """
        EXAMPLES::

            sage: VS = VectorSpace(GF(2),3)
            sage: VS((0,0,1/3))
            (0, 0, 1)
            sage: type(_)
            <class 'sage.modules.vector_mod2_dense.Vector_mod2_dense'>
            sage: VS((0,0,int(3)))
            (0, 0, 1)
            sage: VS((0,0,3))
            (0, 0, 1)
            sage: VS((0,0,GF(2)(1)))
            (0, 0, 1)

        TESTS:

        Check that issue :issue:`8601` is fixed::

            sage: VS = VectorSpace(GF(2), 3)
            sage: VS((-1,-2,-3))
            (1, 0, 1)
            sage: V = VectorSpace(GF(2), 2)
            sage: V([1,3])
            (1, 1)
            sage: V([1,-3])
            (1, 1)

        Check integer overflow prior to :issue:`21746`::

            sage: VS = VectorSpace(GF(2),1)
            sage: VS([2**64])
            (0)
            sage: VS([3**100/5**100])
            (1)

        Check division error over rationals::

            sage: V = VectorSpace(GF(2), 2)
            sage: V([1/3, 3/4])
            Traceback (most recent call last):
            ...
            ZeroDivisionError: inverse does not exist

        Check zero initialization::

            sage: for _ in range(1,100):
            ....:     assert VectorSpace(GF(2), randint(1,5000))(0).is_zero()
            sage: (GF(2)**5)(1)
            Traceback (most recent call last):
            ...
            TypeError: can...t initialize vector from nonzero non-list
            sage: (GF(2)**0).zero_vector()
            ()

        Check construction from numpy arrays::

            sage: # needs numpy
            sage: import numpy
            sage: VS = VectorSpace(GF(2),3)
            sage: VS(numpy.array([0,-3,7], dtype=numpy.int8))
            (0, 1, 1)
            sage: VS(numpy.array([0,-3,7], dtype=numpy.int32))
            (0, 1, 1)
            sage: VS(numpy.array([0,-3,7], dtype=numpy.int64))
            (0, 1, 1)
            sage: VS(numpy.array([False,True,False], dtype=bool))
            (0, 1, 0)
            sage: VS(numpy.array([[1]]))
            Traceback (most recent call last):
            ...
            ValueError: numpy array must have dimension 1
            sage: VS(numpy.array([1,2,3,4]))
            Traceback (most recent call last):
            ...
            ValueError: numpy array must have the right length

        Make sure it's reasonably fast::

            sage: # needs numpy
            sage: import numpy
            sage: VS = VectorSpace(GF(2),2*10^7)
            sage: v = VS(numpy.random.randint(0, 1, size=VS.dimension()))  # around 300ms
        """
        try:
            import numpy
        except ImportError:
            pass
        else:
            from .numpy_util import set_mzd_from_numpy
            if set_mzd_from_numpy(<uintptr_t>self._entries, self._degree, x):
                return
        if isinstance(x, (list, tuple)):
            if len(x) != self._degree:
                raise TypeError("x must be a list of the right length")
            for i in range(len(x)):
                xi = x[i]
                if isinstance(xi, (IntegerMod_int, int, Integer)):
                    # the if/else statement is because in some compilers, (-1)%2 is -1
                    mzd_write_bit(self._entries, 0, i, 1 if xi % 2 else 0)
                elif isinstance(xi, Rational):
                    if not (xi.denominator() % 2):
                        raise ZeroDivisionError("inverse does not exist")
                    mzd_write_bit(self._entries, 0, i, 1 if (xi.numerator() % 2) else 0)
                else:
                    mzd_write_bit(self._entries, 0, i, xi % 2)
        elif x != 0:
            raise TypeError("can't initialize vector from nonzero non-list")
        elif self._degree:
            mzd_set_ui(self._entries, 0)

    def __dealloc__(self):
        """
        EXAMPLES::

            sage: VS = VectorSpace(GF(2),10^3)
            sage: import gc
            sage: for i in range(10):
            ....:     v = VS.random_element()
            ....:     del v
            ....:     _ = gc.collect()
        """
        if self._entries:
            mzd_free(self._entries)

    cpdef _richcmp_(left, right, int op):
        """
        EXAMPLES::

            sage: v = vector(GF(2), [0,0,0,0])
            sage: v == 0
            True
            sage: v == 1
            False
            sage: v == v
            True
            sage: w = vector(GF(2), [1,0,0,0])
            sage: w < v
            False
            sage: w > v
            True
            sage: w = vector(GF(2), [-1,0,0,0])
            sage: w == w
            True
        """
        cdef int c
        if left._degree == 0:
            return rich_to_bool(op, 0)
        c = mzd_cmp(left._entries, (<Vector_mod2_dense>right)._entries)
        return rich_to_bool(op, c)

    cdef get_unsafe(self, Py_ssize_t i):
        """
        EXAMPLES::

            sage: v = vector(GF(2), [1,2,3]); v
            (1, 0, 1)
            sage: v[0]
            1
            sage: v[2]
            1
            sage: v[-2]
            0
            sage: v[0:2]
            (1, 0)
        """
        return self._base_ring(mzd_read_bit(self._entries, 0, i))

    cdef int set_unsafe(self, Py_ssize_t i, value) except -1:
        """
        EXAMPLES::

            sage: VS = VectorSpace(GF(2),4)
            sage: v = VS.random_element()
            sage: v[0] = 0; v[0]
            0
            sage: v[1:3] = [1, 1]; v[1:3]
            (1, 1)
            sage: v[3] = 0; v
            (0, 1, 1, 0)
            sage: v[4] = 0
            Traceback (most recent call last):
            ...
            IndexError: vector index out of range
        """
        mzd_write_bit(self._entries, 0, i, value)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: VS = VectorSpace(GF(2),10^4)
            sage: e = VS.random_element()
            sage: loads(dumps(e)) == e
            True
        """
        return unpickle_v0, (self._parent, self.list(), self._degree,
                             self._is_immutable)

    cpdef _add_(self, right):
        """
        EXAMPLES::

            sage: VS = VectorSpace(GF(2),10)
            sage: e = VS([0,0,1,1,0,0,1,1,0,0])
            sage: f = VS([0,1,0,1,0,1,0,1,0,1])
            sage: e + f #indirect doctest
            (0, 1, 1, 0, 0, 1, 1, 0, 0, 1)
        """
        cdef Vector_mod2_dense z = self._new_c()
        if self._degree:
            mzd_add(z._entries, self._entries, (<Vector_mod2_dense>right)._entries)
        return z

    cpdef _sub_(self, right):
        """
        EXAMPLES::

            sage: VS = VectorSpace(GF(2),10)
            sage: e = VS([0,0,1,1,0,0,1,1,0,0])
            sage: f = VS([0,1,0,1,0,1,0,1,0,1])
            sage: e - f #indirect doctest
            (0, 1, 1, 0, 0, 1, 1, 0, 0, 1)
        """
        cdef Vector_mod2_dense z = self._new_c()
        if self._degree:
            mzd_add(z._entries, self._entries, (<Vector_mod2_dense>right)._entries)
        return z

    cpdef int hamming_weight(self) noexcept:
        """
        Return the number of positions ``i`` such that ``self[i] != 0``.

        EXAMPLES::

            sage: vector(GF(2), [1,1,0]).hamming_weight()
            2
        """
        cdef int i
        cdef int res = 0
        cdef m4ri_word *row = mzd_row(self._entries, 0)
        for i from 0 <= i < self._entries.width:
            res += Integer(row[i]).popcount()
        return res

    cpdef _dot_product_(self, Vector right):
        """
        EXAMPLES::

           sage: VS = VectorSpace(GF(2),3)
           sage: v = VS([1,1,1]); w = VS([0,0,0])
           sage: v * w, w * v #indirect doctest
           (0, 0)
           sage: v = VS([1,1,1]); w = VS([0,1,0])
           sage: v * w, w * v
           (1, 1)
           sage: v = VS([1,1,1]); w = VS([0,1,1])
           sage: v * w, w * v
           (0, 0)
           sage: v = VS([1,1,1]); w = VS([1,1,1])
           sage: v * w, w * v
           (1, 1)

           sage: VS = VectorSpace(GF(2),10^4)
           sage: v = VS(0); w = VS(0)
           sage: v[1337] = 1; w[1337] = 1
           sage: v * w, w * v
           (1, 1)
           sage: v[9881] = 1; w[9881] = 1
           sage: v * w, w * v
           (0, 0)
           sage: v[5172] = 1; w[6178] = 1
           sage: v * w, w * v
           (0, 0)
        """
        cdef Py_ssize_t i
        cdef IntegerMod_int n
        cdef Vector_mod2_dense r = right
        cdef m4ri_word tmp = 0
        n =  IntegerMod_int.__new__(IntegerMod_int)
        IntegerMod_abstract.__init__(n, self.base_ring())
        n.ivalue = 0
        cdef m4ri_word *lrow = mzd_row(self._entries, 0)
        cdef m4ri_word *rrow = mzd_row(r._entries, 0)
        for i from 0 <= i < self._entries.width:
            tmp ^= lrow[i] & rrow[i]

        for i in range(64):
            n.ivalue ^= <int>(tmp & 1)
            tmp = tmp >> 1

        return n

    cpdef _pairwise_product_(self, Vector right):
        """
        EXAMPLES::

            sage: VS = VectorSpace(GF(2),10)
            sage: e = VS.random_element()
            sage: f = VS.random_element()
            sage: g = e.pairwise_product(f) #indirect doctest
            sage: all(g[i] == e[i]*f[i] for i in range(10))
            True
        """
        cdef Vector_mod2_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        cdef m4ri_word *lrow = mzd_row(self._entries, 0)
        cdef m4ri_word *rrow = mzd_row(r._entries, 0)
        cdef m4ri_word *zrow = mzd_row(z._entries, 0)
        for i from 0 <= i < self._entries.width:
            zrow[i] = (lrow[i] & rrow[i])
        return z

    cpdef _lmul_(self, Element left):
        """
        EXAMPLES::

            sage: VS = VectorSpace(GF(2),10)
            sage: e = VS.random_element()
            sage: 0 * e
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            sage: 1 * e == e
            True
            sage: 2 * e  # indirect doctest
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

        ::

            sage: VS = VectorSpace(GF(2), 100)
            sage: e = VS.random_element()
            sage: e * 0 == 0  # indirect doctest
            True
            sage: e * 1 == e
            True
            sage: e * 2 == 0
            True
        """
        if left:
            return self.__copy__()
        return self._new_c()

    cpdef _neg_(self):
        """
        EXAMPLES::

            sage: VS = VectorSpace(GF(2),10)
            sage: e = VS.random_element()
            sage: -e == e
            True
        """
        return self.__copy__()

    def list(self, copy=True):
        """
        Return a list of entries in ``self``.

        INPUT:

        - ``copy`` -- always ``True``

        EXAMPLES::

            sage: VS = VectorSpace(GF(2), 10)
            sage: entries = [GF(2).random_element() for _ in range(10)]
            sage: e = VS(entries)
            sage: e.list() == entries
            True
        """
        cdef Py_ssize_t d = self._degree
        cdef Py_ssize_t i
        cdef list v = [0]*d
        K = self.base_ring()
        z = K.zero()
        o = K.one()
        cdef list switch = [z,o]
        for i in range(d):
            v[i] = switch[mzd_read_bit(self._entries, 0, i)]
        return v


def unpickle_v0(parent, entries, degree, is_immutable):
    """
    EXAMPLES::

        sage: from sage.modules.vector_mod2_dense import unpickle_v0
        sage: VS = VectorSpace(GF(2),10)
        sage: unpickle_v0(VS, [0,1,2,3,4,5,6,7,8,9], 10, 0)
        (0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
    """
    # If you think you want to change this function, don't.
    cdef Vector_mod2_dense v
    v = Vector_mod2_dense.__new__(Vector_mod2_dense)
    v._init(degree, parent)
    cdef int xi

    for i in range(degree):
        if isinstance(entries[i], (IntegerMod_int, int, Integer)):
            xi = entries[i]
            mzd_write_bit(v._entries, 0, i, xi % 2)
        else:
            mzd_write_bit(v._entries, 0, i, entries[i] % 2)
    v._is_immutable = int(is_immutable)
    return v
