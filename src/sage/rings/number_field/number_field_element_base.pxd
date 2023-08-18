# sage_setup: distribution = sagemath-flint
from sage.structure.element cimport FieldElement


cdef class NumberFieldElement_base(FieldElement):
    pass
