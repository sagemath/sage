# sage_setup: distribution = sagemath-plot
from sage.structure.sage_object cimport SageObject


cdef class Graphics3d(SageObject):
    cdef public object texture
    cdef object _aspect_ratio
    cdef object _frame_aspect_ratio
    cdef public dict _extra_kwds

cdef class PrimitiveObject(Graphics3d):
    pass
