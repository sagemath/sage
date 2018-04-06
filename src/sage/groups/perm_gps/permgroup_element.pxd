from sage.structure.element cimport MultiplicativeGroupElement, MonoidElement, Element
from sage.structure.list_clone cimport ClonableIntArray

cdef class PermutationGroupElement(MultiplicativeGroupElement):
    cdef int* perm
    cdef int n
    cdef int perm_buf[15] # to avoid malloc for small elements
    cdef __gap
    cdef Element _gap_element
    cdef __tuple
    cdef PermutationGroupElement _new_c(self)
    cpdef _mul_(self, other)
    cpdef _gap_list(self)
    cpdef domain(self)
    cdef public __custom_name
    cpdef list _act_on_list_on_position(self, list x)
    cpdef ClonableIntArray _act_on_array_on_position(self, ClonableIntArray x)
