# sage_setup: distribution = sagemath-gap
from sage.structure.element cimport MultiplicativeGroupElement, MonoidElement, Element
from sage.structure.list_clone cimport ClonableIntArray
from sage.rings.polynomial.polydict cimport ETuple
from sage.libs.gap.element cimport GapElement

cdef class PermutationGroupElement(MultiplicativeGroupElement):
    cdef int* perm
    cdef int n
    cdef int perm_buf[15] # to avoid malloc for small elements
    cdef GapElement _libgap
    cdef PermutationGroupElement _new_c(self)
    cdef _alloc(self, int)
    cpdef _set_identity(self)
    cpdef _set_list_images(self, v, bint convert)
    cpdef _set_libgap(self, GapElement p)
    cpdef _set_list_cycles(self, c, bint convert)
    cpdef _set_string(self, str s)
    cpdef _set_permutation_group_element(self, PermutationGroupElement p, bint convert)

    cpdef _mul_(self, other)
    cpdef PermutationGroupElement _transpose_left(self, j, k)
    cpdef PermutationGroupElement _generate_new(self, list new_list)
    cpdef PermutationGroupElement _generate_new_GAP(self, old)
    cpdef _gap_list(self)
    cpdef domain(self)
    cdef public _SageObject__custom_name
    cpdef list _act_on_list_on_position(self, list x)
    cpdef ClonableIntArray _act_on_array_on_position(self, ClonableIntArray x)
    cpdef ETuple _act_on_etuple_on_position(self, ETuple x)
