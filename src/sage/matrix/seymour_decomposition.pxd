# sage_setup: distribution = sagemath-cmr
from sage.libs.cmr.cmr cimport CMR_MATROID_DEC
from sage.structure.sage_object cimport SageObject


cdef class DecompositionNode(SageObject):

    cdef CMR_MATROID_DEC *_dec
    cdef object _row_keys
    cdef object _column_keys
    cdef DecompositionNode _root   # my CMR_MATROID_DEC is owned by this

    cdef _set_dec(self, CMR_MATROID_DEC *dec, root)


cdef create_DecompositionNode(CMR_MATROID_DEC *dec, root=?)
