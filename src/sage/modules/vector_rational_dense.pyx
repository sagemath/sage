# sage_setup: distribution = sagemath-modules
"""
Vectors with rational entries

AUTHOR:

- William Stein (2007)
- Soroosh Yazdani (2007)

EXAMPLES::

    sage: v = vector(QQ,[1,2,3,4,5])
    sage: v
    (1, 2, 3, 4, 5)
    sage: 3*v
    (3, 6, 9, 12, 15)
    sage: v/2
    (1/2, 1, 3/2, 2, 5/2)
    sage: -v
    (-1, -2, -3, -4, -5)
    sage: v - v
    (0, 0, 0, 0, 0)
    sage: v + v
    (2, 4, 6, 8, 10)
    sage: v * v
    55

We make a large zero vector::

    sage: k = QQ^100000; k
    Vector space of dimension 100000 over Rational Field
    sage: v = k(0)
    sage: v[:10]
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

TESTS::

    sage: v = vector(QQ, [1,2/5,-3/8,4])
    sage: loads(dumps(v)) == v
    True

    sage: w = vector(QQ, [-1,0,0,0])
    sage: w.set_immutable()
    sage: isinstance(hash(w), int)
    True
"""

# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport check_allocarray, sig_free
from cysignals.signals cimport sig_on, sig_off

from sage.structure.element cimport Element, Vector
from sage.structure.richcmp cimport rich_to_bool
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational

cimport sage.modules.free_module_element as free_module_element

from sage.libs.gmp.mpq cimport *


cdef inline _Rational_from_mpq(mpq_t e):
    cdef Rational z = Rational.__new__(Rational)
    mpq_set(z.value, e)
    return z

cdef class Vector_rational_dense(free_module_element.FreeModuleElement):
    cdef bint is_dense_c(self) noexcept:
        return 1
    cdef bint is_sparse_c(self) noexcept:
        return 0

    def __copy__(self):
        cdef Vector_rational_dense y
        y = self._new_c()
        cdef Py_ssize_t i
        for i in range(self._degree):
            mpq_set(y._entries[i], self._entries[i])
        return y

    cdef int _init(self, Py_ssize_t degree, Parent parent) except -1:
        """
        Initialize the C data structures in this vector. After calling
        this, ``self`` is a zero vector of degree ``degree`` with
        parent ``parent``.

        Only if you call ``__new__`` without a ``parent`` argument, it
        is needed to call this function manually. The only reason to do
        that is for efficiency: calling ``__new__`` without any
        additional arguments besides the type and then calling ``_init``
        is (slightly) more efficient than calling ``__new__`` with a
        ``parent`` argument.

        TESTS:

        Check implicitly that :issue:`10257` works::

            sage: from sage.modules.vector_rational_dense import Vector_rational_dense
            sage: Vector_rational_dense(QQ^(sys.maxsize))
            Traceback (most recent call last):
            ...
            MemoryError: failed to allocate ... bytes
            sage: try:
            ....:     # Note: some malloc() implementations (on OS X
            ....:     # for example) print stuff when an allocation
            ....:     # fails.  # We catch this with the ... in the
            ....:     # doctest result. The * is needed because a
            ....:     # result cannot start with ...
            ....:     print("*")
            ....:     Vector_rational_dense(QQ^(2^56))
            ....: except (MemoryError, OverflowError):
            ....:     print("allocation failed")
            *...
            allocation failed
        """
        # Assign variables only when the array is fully initialized
        cdef mpq_t* entries = <mpq_t*>check_allocarray(degree, sizeof(mpq_t))
        cdef Py_ssize_t i
        sig_on()
        for i in range(degree):
            mpq_init(entries[i])
        sig_off()
        self._entries = entries
        self._degree = degree
        self._parent = parent

    def __cinit__(self, parent=None, x=None, coerce=True, copy=True):
        self._entries = NULL
        self._is_immutable = 0
        if parent is None:
            self._degree = 0
            return
        self._init(parent.degree(), <Parent?>parent)

    def __init__(self, parent, x, coerce=True, copy=True):
        cdef Py_ssize_t i
        cdef Rational z
        if isinstance(x, (list, tuple)):
            if len(x) != self._degree:
                raise TypeError("entries must be a list of length %s" % self._degree)
            for i from 0 <= i < self._degree:
                z = Rational(x[i])
                mpq_set(self._entries[i], z.value)
            return
        if x != 0:
            raise TypeError("can't initialize vector from nonzero non-list")

    def __dealloc__(self):
        cdef Py_ssize_t i
        if self._entries:
            # Do *not* use sig_on() here, since __dealloc__
            # cannot raise exceptions!
            for i in range(self._degree):
                mpq_clear(self._entries[i])
            sig_free(self._entries)

    cpdef _richcmp_(left, right, int op):
        """
        EXAMPLES::

            sage: v = vector(QQ, [0,0,0,0])
            sage: v == 0
            True
            sage: v == 1
            False
            sage: v == v
            True
            sage: w = vector(QQ, [-1,3/2,0,0])
            sage: w < v
            True
            sage: w > v
            False
            sage: w = vector(QQ, [-1,0,0,0])
            sage: w == w
            True
        """
        cdef Py_ssize_t i
        cdef int c
        for i in range(left.degree()):
            c = mpq_cmp(left._entries[i], (<Vector_rational_dense>right)._entries[i])
            if c < 0:
                return rich_to_bool(op, -1)
            elif c > 0:
                return rich_to_bool(op, 1)
        return rich_to_bool(op, 0)

    cdef get_unsafe(self, Py_ssize_t i):
        """
        EXAMPLES::

            sage: v = vector([1/2,2/3,3/4]); v
            (1/2, 2/3, 3/4)
            sage: v[0]
            1/2
            sage: v[2]
            3/4
            sage: v[-2]
            2/3
            sage: v[5]
            Traceback (most recent call last):
            ...
            IndexError: vector index out of range
            sage: v[-5]
            Traceback (most recent call last):
            ...
            IndexError: vector index out of range
        """
        cdef Rational z = Rational.__new__(Rational)
        mpq_set(z.value, self._entries[i])
        return z

    cdef int set_unsafe(self, Py_ssize_t i, value) except -1:
        """
        EXAMPLES::

            sage: v = vector(QQ, [1/2,2/5,0]); v
            (1/2, 2/5, 0)
            sage: v.set(2, -15/17); v
            (1/2, 2/5, -15/17)
        """
        mpq_set(self._entries[i], (<Rational>value).value)

    def list(self, copy=True):
        """
        The list of entries of the vector.

        INPUT:

        - ``copy`` -- ignored optional argument

        EXAMPLES::

            sage: v = vector(QQ,[1,2,3,4])
            sage: a = v.list(copy=False); a
            [1, 2, 3, 4]
            sage: a[0] = 0
            sage: v
            (1, 2, 3, 4)
        """
        cdef int i
        return [_Rational_from_mpq(self._entries[i])
                for i in range(self._degree)]

    def __reduce__(self):
        return (unpickle_v1, (self._parent, self.list(), self._degree,
                              not self._is_immutable))

    cpdef _add_(self, right):
        cdef Vector_rational_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i in range(self._degree):
            mpq_add(z._entries[i], self._entries[i], r._entries[i])
        return z

    cpdef _sub_(self, right):
        cdef Vector_rational_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i in range(self._degree):
            mpq_sub(z._entries[i], self._entries[i], r._entries[i])
        return z

    cpdef _dot_product_(self, Vector right):
        """
        Dot product of dense vectors over the rationals.

        EXAMPLES::

            sage: v = vector(QQ, [1,2,-3]); w = vector(QQ,[4,3,2])
            sage: v*w
            4
            sage: w*v
            4
        """
        cdef Vector_rational_dense r = right
        cdef Rational z
        z = Rational.__new__(Rational)
        cdef mpq_t t
        mpq_init(t)
        mpq_set_si(z.value, 0, 1)
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            mpq_mul(t, self._entries[i], r._entries[i])
            mpq_add(z.value, z.value, t)
        mpq_clear(t)
        return z

    cpdef _pairwise_product_(self, Vector right):
        """
        EXAMPLES::

            sage: v = vector(QQ, [1,2,-3]); w = vector(QQ,[4,3,2])
            sage: v.pairwise_product(w)
            (4, 6, -6)
        """
        cdef Vector_rational_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i  in range(self._degree):
            mpq_mul(z._entries[i], self._entries[i], r._entries[i])
        return z

    cpdef _rmul_(self, Element left):
        cdef Vector_rational_dense z
        cdef Rational a
        if isinstance(left, Rational):
            a = <Rational>left
        elif isinstance(left, Integer):
            a = <Rational>Rational.__new__(Rational)
            mpq_set_z(a.value, (<Integer>left).value)
        else:
            # should not happen
            raise TypeError("cannot convert %s to %s" % (type(left).__name__,
                                                         Rational.__name__))
        z = self._new_c()
        cdef Py_ssize_t i
        for i in range(self._degree):
            mpq_mul(z._entries[i], self._entries[i], a.value)
        return z

    cpdef _lmul_(self, Element right):
        cdef Vector_rational_dense z
        cdef Rational a
        if isinstance(right, Rational):
            a = <Rational>right
        elif isinstance(right, Integer):
            a = <Rational>Rational.__new__(Rational)
            mpq_set_z(a.value, (<Integer>right).value)
        else:
            # should not happen
            raise TypeError("cannot convert %s to %s" % (type(right).__name__,
                                                         Rational.__name__))
        z = self._new_c()
        cdef Py_ssize_t i
        for i in range(self._degree):
            mpq_mul(z._entries[i], self._entries[i], a.value)
        return z

    cpdef _neg_(self):
        cdef Vector_rational_dense z
        z = self._new_c()
        cdef Py_ssize_t i
        for i in range(self._degree):
            mpq_neg(z._entries[i], self._entries[i])
        return z


def unpickle_v0(parent, entries, degree):
    # If you think you want to change this function, don't.
    # Instead make a new version with a name like
    #    make_FreeModuleElement_generic_dense_v1
    # and changed the reduce method below.
    cdef Vector_rational_dense v
    v = Vector_rational_dense.__new__(Vector_rational_dense)
    v._init(degree, parent)
    cdef Rational z
    cdef Py_ssize_t i
    for i in range(degree):
        z = Rational(entries[i])
        mpq_set(v._entries[i], z.value)
    return v


def unpickle_v1(parent, entries, degree, is_mutable):
    cdef Vector_rational_dense v
    v = Vector_rational_dense.__new__(Vector_rational_dense)
    v._init(degree, parent)
    cdef Rational z
    cdef Py_ssize_t i
    for i in range(degree):
        z = Rational(entries[i])
        mpq_set(v._entries[i], z.value)
    v._is_immutable = not is_mutable
    return v
