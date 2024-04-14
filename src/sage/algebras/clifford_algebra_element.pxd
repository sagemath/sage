# sage_setup: distribution = sagemath-modules
"""
Clifford algebra elements
"""
from sage.modules.with_basis.indexed_element cimport IndexedFreeModuleElement
from sage.data_structures.bitset cimport FrozenBitset

cdef class CliffordAlgebraElement(IndexedFreeModuleElement):
    cdef CliffordAlgebraElement _mul_self_term(self, FrozenBitset supp, coeff)
    cdef CliffordAlgebraElement _mul_term_self(self, FrozenBitset supp, coeff)

cdef class ExteriorAlgebraElement(CliffordAlgebraElement):
    pass

cdef class CohomologyRAAGElement(CliffordAlgebraElement):
    pass
