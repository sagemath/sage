# sage_setup: distribution = sagemath-cmr
from sage.libs.cmr.cmr cimport CMR_DEC
from sage.structure.sage_object cimport SageObject


cdef class DecompositionNode(SageObject):

    cdef CMR_DEC *_dec
    cdef DecompositionNode _root   # my CMR_DEC is owned by this

    cdef _set_dec(self, CMR_DEC *dec, root)


cdef create_DecompositionNode(CMR_DEC *dec, root=?)
