"""
Polynomial Template for C/C++ Library Interfaces
"""

from sage.rings.polynomial.polynomial_element cimport Polynomial

cdef class Polynomial_template(Polynomial):
    cdef celement x
    cdef cparent _cparent
    cpdef _mod_(self, right)
