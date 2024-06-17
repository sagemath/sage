# sage_setup: distribution = sagemath-objects
cimport cython
from cpython.object cimport PyObject


cdef struct mono_cell:
    PyObject* key_id
    PyObject* key_weakref
    PyObject* value


@cython.final
cdef class MonoDict:
    cdef __weakref__
    cdef size_t mask  # size - 1 with 'size' the length of the table array
    cdef size_t used  # number of valid entries
    cdef size_t fill  # number of non-NULL entries (including deleted entries)
    cdef mono_cell* table
    cdef bint weak_values
    cdef eraser

    cdef mono_cell* lookup(self, PyObject* key) noexcept
    cdef get(self, k)
    cdef int set(self, k, value) except -1
    cdef int resize(self) except -1


cdef struct triple_cell:
    PyObject* key_id1
    PyObject* key_id2
    PyObject* key_id3
    PyObject* key_weakref1
    PyObject* key_weakref2
    PyObject* key_weakref3
    PyObject* value


@cython.final
cdef class TripleDict:
    cdef __weakref__
    cdef size_t mask  # size - 1 with 'size' the length of the table array
    cdef size_t used  # number of valid entries
    cdef size_t fill  # number of non-NULL entries (including deleted entries)
    cdef triple_cell* table
    cdef bint weak_values
    cdef eraser

    cdef triple_cell* lookup(self, PyObject* key1, PyObject* key2, PyObject* key3) noexcept
    cdef get(self, k1, k2, k3)
    cdef int set(self, k1, k2, k3, value) except -1
    cdef int resize(self) except -1
