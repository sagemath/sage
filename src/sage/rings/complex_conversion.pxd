# sage_setup: distribution = sagemath-modules
from sage.structure.element cimport Element
from sage.categories.map cimport Map


cdef class CCtoCDF(Map):

    cpdef Element _call_(self, x)
