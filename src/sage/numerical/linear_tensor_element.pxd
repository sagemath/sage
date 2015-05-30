from sage.structure.element cimport Element, ModuleElement

cdef class LinearTensor(ModuleElement):
    cpdef dict _f
