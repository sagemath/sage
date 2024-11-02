"""
Vectors with integer entries

AUTHOR:

- William Stein (2007)

EXAMPLES::

    sage: v = vector(ZZ, [1,2,3,4,5])
    sage: v
    (1, 2, 3, 4, 5)
    sage: 3*v
    (3, 6, 9, 12, 15)
    sage: v*7
    (7, 14, 21, 28, 35)
    sage: -v
    (-1, -2, -3, -4, -5)
    sage: v - v
    (0, 0, 0, 0, 0)
    sage: v + v
    (2, 4, 6, 8, 10)
    sage: v * v   # dot product.
    55

We make a large zero vector::

    sage: k = ZZ^100000; k
    Ambient free module of rank 100000 over the principal ideal domain Integer Ring
    sage: v = k(0)
    sage: v[:10]
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

TESTS::

    sage: v = vector(ZZ, [1,2,3,4])
    sage: loads(dumps(v)) == v
    True

    sage: w = vector(ZZ, [-1,0,0,0])
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
from sage.rings.integer cimport Integer, _Integer_from_mpz

cimport sage.modules.free_module_element as free_module_element

from sage.libs.gmp.mpz cimport *

cdef class Vector_integer_dense(free_module_element.FreeModuleElement):
    cdef bint is_dense_c(self) noexcept:
        return 1
    cdef bint is_sparse_c(self) noexcept:
        return 0

    def __copy__(self):
        cdef Vector_integer_dense y
        y = self._new_c()
        cdef Py_ssize_t i
        for i in range(self._degree):
            mpz_set(y._entries[i], self._entries[i])
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
        """
        # Assign variables only when the array is fully initialized
        cdef mpz_t* entries = <mpz_t*>check_allocarray(degree, sizeof(mpz_t))
        cdef Py_ssize_t i
        sig_on()
        for i in range(degree):
            mpz_init(entries[i])
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
        cdef Integer z
        if isinstance(x, (list, tuple)):
            if len(x) != self._degree:
                raise TypeError("x must be a list of the right length")
            for i from 0 <= i < self._degree:
                z = Integer(x[i])
                mpz_set(self._entries[i], z.value)
            return
        if x != 0:
            raise TypeError("can't initialize vector from nonzero non-list")

    def __dealloc__(self):
        cdef Py_ssize_t i
        if self._entries:
            for i in range(self._degree):
                mpz_clear(self._entries[i])
            sig_free(self._entries)

    cpdef _richcmp_(left, right, int op):
        """
        EXAMPLES::

            sage: v = vector(ZZ, [0,0,0,0])
            sage: v == 0
            True
            sage: v == 1
            False
            sage: v == v
            True
            sage: w = vector(ZZ, [-1,0,0,0])
            sage: w == w
            True
            sage: w < v
            True
            sage: w > v
            False
        """
        cdef Py_ssize_t i
        cdef int c
        for i in range(left.degree()):
            c = mpz_cmp(left._entries[i], (<Vector_integer_dense>right)._entries[i])
            if c < 0:
                return rich_to_bool(op, -1)
            elif c > 0:
                return rich_to_bool(op, 1)
        return rich_to_bool(op, 0)

    cdef get_unsafe(self, Py_ssize_t i):
        """
        EXAMPLES::

            sage: v = vector([1,2,3]); v
            (1, 2, 3)
            sage: v[0]
            1
            sage: v[-2]
            2
            sage: v[0:2]
            (1, 2)
            sage: v[::-1]
            (3, 2, 1)
        """
        cdef Integer z = Integer.__new__(Integer)
        mpz_set(z.value, self._entries[i])
        return z

    cdef int set_unsafe(self, Py_ssize_t i, value) except -1:
        """
        EXAMPLES::

            sage: v = vector([1,2,3]); v
            (1, 2, 3)
            sage: v[0] = 2
            sage: v[1:3] = [1, 4]; v
            (2, 1, 4)
        """
        mpz_set(self._entries[i], (<Integer>value).value)

    def list(self, copy=True):
        """
        The list of entries of the vector.

        INPUT:

        - ``copy`` -- ignored optional argument

        EXAMPLES::

            sage: v = vector([1,2,3,4])
            sage: a = v.list(copy=False); a
            [1, 2, 3, 4]
            sage: a[0] = 0
            sage: v
            (1, 2, 3, 4)
        """
        cdef int i
        return [_Integer_from_mpz(self._entries[i])
                for i in range(self._degree)]

    def __reduce__(self):
        return (unpickle_v1, (self._parent, self.list(), self._degree,
                              not self._is_immutable))

    cpdef _add_(self, right):
        cdef Vector_integer_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i in range(self._degree):
            mpz_add(z._entries[i], self._entries[i], r._entries[i])
        return z

    cpdef _sub_(self, right):
        cdef Vector_integer_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i in range(self._degree):
            mpz_sub(z._entries[i], self._entries[i], r._entries[i])
        return z

    cpdef _dot_product_(self, Vector right):
        """
        Dot product of dense vectors over the integers.

        EXAMPLES::

            sage: v = vector(ZZ, [1,2,-3]); w = vector(ZZ,[4,3,2])
            sage: v*w
            4
            sage: w*v
            4
        """
        cdef Vector_integer_dense r = right
        cdef Integer z = Integer.__new__(Integer)
        cdef mpz_t t
        mpz_init(t)
        mpz_set_si(z.value, 0)
        cdef Py_ssize_t i
        for i in range(self._degree):
            mpz_mul(t, self._entries[i], r._entries[i])
            mpz_add(z.value, z.value, t)
        mpz_clear(t)
        return z

    cpdef _pairwise_product_(self, Vector right):
        """
        EXAMPLES::

            sage: v = vector(ZZ, [1,2,-3]); w = vector(ZZ,[4,3,2])
            sage: v.pairwise_product(w)
            (4, 6, -6)
        """
        cdef Vector_integer_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i in range(self._degree):
            mpz_mul(z._entries[i], self._entries[i], r._entries[i])
        return z

    cpdef _rmul_(self, Element left):
        cdef Vector_integer_dense z
        cdef Integer a
        a = left
        z = self._new_c()
        cdef Py_ssize_t i
        for i in range(self._degree):
            mpz_mul(z._entries[i], self._entries[i], a.value)
        return z

    cpdef _lmul_(self, Element right):
        cdef Vector_integer_dense z
        cdef Integer a
        a = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i in range(self._degree):
            mpz_mul(z._entries[i], self._entries[i], a.value)
        return z

    cpdef _neg_(self):
        cdef Vector_integer_dense z
        z = self._new_c()
        cdef Py_ssize_t i
        for i in range(self._degree):
            mpz_neg(z._entries[i], self._entries[i])
        return z

    def _singular_(self, singular=None):
        r"""
        Return \Singular representation of this integer vector.

        INPUT:

        - ``singular`` -- \Singular interface instance (default: ``None``)

        EXAMPLES::

            sage: A = random_matrix(ZZ,1,3)
            sage: v = A.row(0)
            sage: vs = singular(v)                                                      # needs sage.libs.singular
            sage: vs._repr_() == '{},\n{},\n{}'.format(*v)                              # needs sage.libs.singular
            True
            sage: vs.type()                                                             # needs sage.libs.singular
            'intvec'
        """
        if singular is None:
            from sage.interfaces.singular import singular as singular_default
            singular = singular_default

        name = singular._next_var_name()
        values = str(self.list())[1:-1]
        singular.eval("intvec %s = %s"%(name, values))

        from sage.interfaces.singular import SingularElement
        return SingularElement(singular, 'foobar', name, True)


def unpickle_v0(parent, entries, degree):
    # If you think you want to change this function, don't.
    # Instead make a new version with a name like
    #    make_FreeModuleElement_generic_dense_v1
    # and changed the reduce method below.
    cdef Vector_integer_dense v
    v = Vector_integer_dense.__new__(Vector_integer_dense)
    v._init(degree, parent)
    cdef Integer z
    cdef Py_ssize_t i
    for i in range(degree):
        z = Integer(entries[i])
        mpz_set(v._entries[i], z.value)
    return v


def unpickle_v1(parent, entries, degree, is_mutable):
    cdef Vector_integer_dense v
    v = Vector_integer_dense.__new__(Vector_integer_dense)
    v._init(degree, parent)
    cdef Integer z
    cdef Py_ssize_t i
    for i in range(degree):
        z = Integer(entries[i])
        mpz_set(v._entries[i], z.value)
    v._is_immutable = not is_mutable
    return v
