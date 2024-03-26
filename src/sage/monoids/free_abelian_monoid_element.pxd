from sage.structure.element cimport MonoidElement
from sage.libs.gmp.types cimport mpz_t
from sage.structure.parent cimport Parent

cdef class FreeAbelianMonoidElement(MonoidElement):
    cdef mpz_t *_element_vector
    cdef Py_ssize_t _n

    cdef int _init(self, Py_ssize_t n, Parent parent) except -1

    cdef inline FreeAbelianMonoidElement _new_c(self):
        cdef type t = type(self)
        cdef FreeAbelianMonoidElement x = <FreeAbelianMonoidElement>(t.__new__(t))
        x._init(self._n, self._parent)
        return x
