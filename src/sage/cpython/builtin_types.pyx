# sage_setup: distribution = sagemath-objects
from cpython.object cimport PyTypeObject

cdef extern from *:
    PyTypeObject PyWrapperDescr_Type

wrapper_descriptor = <type>(&PyWrapperDescr_Type)
