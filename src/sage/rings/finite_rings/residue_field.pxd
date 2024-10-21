# sage_setup: distribution = sagemath-categories
from sage.rings.morphism cimport RingHomomorphism
from sage.structure.element cimport Element
from sage.rings.finite_rings.integer_mod cimport NativeIntStruct

from sage.categories.map cimport Map, Section

cdef class ReductionMap(Map):
    cdef public object _K
    cdef public object _F
    cdef public object _to_vs
    cdef public object _PBinv
    cdef object _to_order
    cdef object _PB
    cdef object _section

cdef class ResidueFieldHomomorphism_global(RingHomomorphism):
    cdef public object _K
    cdef public object _F
    cdef public object _to_vs
    cdef public object _PBinv
    cdef object _to_order
    cdef object _PB
    cdef object _section

cdef class LiftingMap(Section):
    cdef public object _K
    cdef public object _F
    cdef public object _to_order
    cdef public object _PB
