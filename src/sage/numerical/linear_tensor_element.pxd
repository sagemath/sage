# sage_setup: distribution = sagemath-polyhedra
from sage.structure.element cimport Element, ModuleElement

cdef class LinearTensor(ModuleElement):
    cdef dict _f
    cpdef _add_(self, other)
