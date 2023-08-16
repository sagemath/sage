"""
Octonions
"""

from sage.structure.element cimport AlgebraElement
from sage.modules.free_module_element cimport FreeModuleElement

cdef class Octonion_generic(AlgebraElement):
    cdef FreeModuleElement vec

    cpdef Octonion_generic conjugate(self)
    cpdef quadratic_form(self)
    cpdef norm(self)
    cpdef abs(self)
    cpdef real_part(self)
    cpdef Octonion_generic imag_part(self)

cdef class Octonion(Octonion_generic):
    pass
