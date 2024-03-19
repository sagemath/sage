# sage_setup: distribution = sagemath-cmr
from sage.libs.cmr.cmr cimport CMR_MATROID_DEC, CMR_ELEMENT
from sage.structure.sage_object cimport SageObject


cdef class DecompositionNode(SageObject):
    cdef object _matrix
    cdef CMR_MATROID_DEC *_dec
    cdef object _row_keys
    cdef object _column_keys
    cdef public object _parent_indices
    cdef object _child_nodes
    cdef object _child_row_keys
    cdef object _child_column_keys

    cdef _set_dec(self, CMR_MATROID_DEC *dec)
    cdef _set_row_keys(self, row_keys)
    cdef _set_column_keys(self, column_keys)

    cdef _CMRelement_to_key(self, CMR_ELEMENT element)


cdef class SymbolicNode(DecompositionNode):

    cdef object _symbol


cdef create_DecompositionNode(CMR_MATROID_DEC *dec, row_keys=?, column_keys=?)
