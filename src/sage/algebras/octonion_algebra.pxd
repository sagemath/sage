"""
Octonions
"""

from sage.structure.element cimport AlgebraElement
from sage.modules.free_module_element cimport FreeModuleElement

cdef class Octonion_generic(AlgebraElement):
    cdef FreeModuleElement vec

    cpdef Octonion_generic conjugate(self) noexcept
    cpdef quadratic_form(self) noexcept
    cpdef norm(self) noexcept
    cpdef abs(self) noexcept
    cpdef real_part(self) noexcept
    cpdef Octonion_generic imag_part(self) noexcept

cdef class Octonion(Octonion_generic):
    pass
