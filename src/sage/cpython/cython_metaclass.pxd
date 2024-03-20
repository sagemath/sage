# sage_setup: distribution = sagemath-objects
cdef extern from "cython_metaclass.h":
    PyMethodDescr_CallSelf(desc, self)
