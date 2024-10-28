"""
Direct Access to Singular's Functions via libSingular

AUTHORS:

- Michael Brickenstein (2009-07): initial implementation, overall design
- Martin Albrecht (2009-07): clean up, enhancements, etc.
"""
#*****************************************************************************
#       Copyright (C) 2009 Michael Brickenstein <brickenstein@mfo.de>
#       Copyright (C) 2009 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object cimport SageObject
from sage.libs.singular.decl cimport leftv, idhdl, syStrategy, matrix, poly, ideal, intvec
from sage.libs.singular.decl cimport ring as singular_ring
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular, MPolynomial_libsingular

cdef new_sage_polynomial(ring,  poly *p)
cdef poly* access_singular_poly(p) except <poly*> -1
cdef singular_ring* access_singular_ring(r) except <singular_ring*> -1

cdef class RingWrap:
    cdef singular_ring *_ring

cdef class Resolution:
    cdef syStrategy *_resolution
    cdef object base_ring

cdef class Converter(SageObject):
    cdef leftv *args
    cdef object _sage_ring
    cdef singular_ring* _singular_ring
    cdef leftv* pop_front(self) except NULL
    cdef leftv * _append_leftv(self, leftv *v) noexcept
    cdef leftv * _append(self, void* data, int res_type) noexcept
    cdef leftv * append_polynomial(self, p) except NULL
    cdef leftv * append_ideal(self,  i) except NULL
    cdef leftv * append_number(self, n) except NULL
    cdef leftv * append_int(self, n) except NULL
    cdef leftv * append_str(self, n) except NULL
    cdef leftv * append_intmat(self, a) except NULL
    cdef leftv * append_resolution(self, Resolution resolution) except NULL
    cdef leftv * append_vector(self, v) except NULL
    cdef leftv * append_intvec(self, v) except NULL
    cdef leftv * append_list(self, l) except NULL
    cdef leftv * append_matrix(self, a) except NULL
    cdef leftv * append_ring(self, r) except NULL
    cdef leftv * append_module(self, m) except NULL
    cdef to_sage_integer_matrix(self, intvec *mat)
    cdef object to_sage_module_element_sequence_destructive(self, ideal *i)
    cdef to_sage_vector_destructive(self, poly *p, free_module = ?)
    cdef to_sage_matrix(self, matrix* mat)
    cdef to_python(self, leftv* to_convert)

cdef class BaseCallHandler:
    cdef leftv* handle_call(self, Converter argument_list, singular_ring *_ring=?) noexcept
    cdef bint free_res(self) noexcept

cdef class LibraryCallHandler(BaseCallHandler):
    cdef idhdl * proc_idhdl

cdef class KernelCallHandler(BaseCallHandler):
    cdef long arity
    cdef long cmd_n

cdef class SingularFunction(SageObject):
    cdef object _name
    cdef MPolynomialRing_libsingular _ring
    cdef BaseCallHandler call_handler

    cdef BaseCallHandler get_call_handler(self)
    cdef bint function_exists(self) noexcept
    cdef common_ring(self, tuple args, ring=?)

cdef class SingularLibraryFunction(SingularFunction):
    pass

cdef class SingularKernelFunction(SingularFunction):
    pass

# the most direct function call interface
cdef call_function(SingularFunction self, tuple args, object R, bint signal_handler=?, object attributes=?)
