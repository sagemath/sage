from sage.matroids.matroid cimport Matroid
from sage.matroids.basis_exchange_matroid cimport BasisExchangeMatroid

cdef class TransversalMatroid(BasisExchangeMatroid):
    cdef dict _matching
    cdef object _sets
    cdef object _D
    cdef list _set_labels, _sets_input, _set_labels_input

    cpdef list sets(self)
    cdef dict _translate_matching(self)
    cpdef reduce_presentation(self)
    cpdef transversal_extension(self, element=*, newset=*, sets=*)
