# sage_setup: distribution = sagemath-modules
from sage.structure.element cimport ModuleElement

cdef class FiniteZZsubmodule_iterator:
    #### Global Data
    cdef FiniteZZsubmodule_iterator _other_ZZ
    cdef ModuleElement _basis
    cdef ModuleElement _cw
    cdef ModuleElement _coset_rep
    cdef ModuleElement _other
    cdef list _basis_all
    cdef list _plus
    cdef int _basis_length
    cdef int _count
    cdef int _order
    cdef bint _immutable
    cdef ModuleElement _iteration(FiniteZZsubmodule_iterator self)

cdef class FiniteFieldsubspace_iterator(FiniteZZsubmodule_iterator):
    pass

cdef class FiniteFieldsubspace_projPoint_iterator:
    cdef int _basis_length, _normalized_pos
    cdef int _one_dimensional_case  # takes the values 0, 1, 2
    cdef list _basis
    cdef bint _immutable
    cdef FiniteFieldsubspace_iterator _it
