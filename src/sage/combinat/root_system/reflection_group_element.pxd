from sage.groups.perm_gps.permgroup_element cimport PermutationGroupElement

cdef class ComplexReflectionGroupElement(PermutationGroupElement):
    cpdef action(self, vec, on_space=*)
    cpdef action_on_root_indices(self, i)

cdef class RealReflectionGroupElement(ComplexReflectionGroupElement):
    cpdef bint has_left_descent(self, i) noexcept
    cpdef bint has_descent(self, i, side=*, positive=*) noexcept
    cpdef action(self, vec, side=*, on_space=*)
    cpdef action_on_root_indices(self, i, side=*)
