cimport cython
from cpython cimport PyObject

cdef struct mono_cell

@cython.final
cdef class MonoDict:
    cdef __weakref__
    cdef size_t mask
    cdef size_t used
    cdef size_t fill
    cdef mono_cell* table
    cdef bint weak_values
    cdef eraser
    cdef mono_cell* lookup(self,PyObject* key)
    cdef get(self, object k)
    cdef set(self, object k, value)
    cdef int resize(self) except -1

cdef struct triple_cell

@cython.final
cdef class TripleDict:
    cdef __weakref__
    cdef size_t mask
    cdef size_t used
    cdef size_t fill
    cdef triple_cell* table
    cdef bint weak_values
    cdef eraser
    cdef triple_cell* lookup(self, PyObject* key1, PyObject* key2, PyObject* key3)
    cdef get(self, object k1, object k2, object k3)
    cdef set(self, object k1, object k2, object k3, value)
    cdef int resize(self) except -1
