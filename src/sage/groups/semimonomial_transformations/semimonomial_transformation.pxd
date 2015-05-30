from sage.structure.element cimport Element, MonoidElement, MultiplicativeGroupElement
from sage.rings.finite_rings.integer_mod cimport IntegerMod_int

cdef class SemimonomialTransformation(MultiplicativeGroupElement):
    cdef tuple v
    cdef object perm, alpha

    cdef _new_c(self)
