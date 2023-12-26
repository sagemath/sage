from sage.combinat.integer_lists.base cimport IntegerListsBackend
cdef class IntegerListsBackend_invlex(IntegerListsBackend):
    cdef public bint check
