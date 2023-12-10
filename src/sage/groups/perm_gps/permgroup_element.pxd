from sage.structure.element cimport MultiplicativeGroupElement, MonoidElement, Element
from sage.structure.list_clone cimport ClonableIntArray
from sage.rings.polynomial.polydict cimport ETuple
from sage.libs.gap.element cimport GapElement

cdef class PermutationGroupElement(MultiplicativeGroupElement):
    cdef int* perm
    cdef int n
    cdef int perm_buf[15] # to avoid malloc for small elements
    cdef GapElement _libgap
    cdef PermutationGroupElement _new_c(self) noexcept
    cdef _alloc(self, int) noexcept
    cpdef _set_identity(self) noexcept
    cpdef _set_list_images(self, v, bint convert) noexcept
    cpdef _set_libgap(self, GapElement p) noexcept
    cpdef _set_list_cycles(self, c, bint convert) noexcept
    cpdef _set_string(self, str s) noexcept
    cpdef _set_permutation_group_element(self, PermutationGroupElement p, bint convert) noexcept

    cpdef _mul_(self, other) noexcept
    cpdef PermutationGroupElement _transpose_left(self, j, k) noexcept
    cpdef PermutationGroupElement _generate_new(self, list new_list) noexcept
    cpdef PermutationGroupElement _generate_new_GAP(self, old) noexcept
    cpdef _gap_list(self) noexcept
    cpdef domain(self) noexcept
    cdef public _SageObject__custom_name
    cpdef list _act_on_list_on_position(self, list x) noexcept
    cpdef ClonableIntArray _act_on_array_on_position(self, ClonableIntArray x) noexcept
    cpdef ETuple _act_on_etuple_on_position(self, ETuple x) noexcept
