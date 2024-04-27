from sage.libs.ntl.GF2X cimport GF2X_c

ctypedef GF2X_c celement
ctypedef long cparent

include "polynomial_template_header.pxi"

cdef class Polynomial_GF2X(Polynomial_template):
    pass

