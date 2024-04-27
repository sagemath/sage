from sage.groups.libgap_wrapper cimport ElementLibGAP

cdef class MatrixGroupElement_gap(ElementLibGAP):
    cpdef _act_on_(self, x, bint self_on_left)
    cpdef list list(self)
