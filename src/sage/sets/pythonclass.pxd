# sage_setup: distribution = sagemath-objects
from sage.structure.parent cimport Set_generic


cdef class Set_PythonType_class(Set_generic):
    cdef type _type


cpdef Set_PythonType(typ)
