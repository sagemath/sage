# sage_setup: distribution = sagemath-modules

from sage.libs.gmp.types cimport mpz_t
from sage.modules.free_module_element cimport FreeModuleElement
from sage.structure.parent cimport Parent


cdef class Vector_integer_dense(FreeModuleElement):
    cdef mpz_t* _entries
    cdef int _init(self, Py_ssize_t degree, Parent parent) except -1

    cdef inline Vector_integer_dense _new_c(self):
        cdef type t = type(self)
        cdef Vector_integer_dense x = <Vector_integer_dense>(t.__new__(t))
        x._init(self._degree, self._parent)
        return x
