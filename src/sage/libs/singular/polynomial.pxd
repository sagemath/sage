"""
Wrapper for Singular's Polynomial Arithmetic

AUTHOR:

- Martin Albrecht (2009-07): initial implementation
"""
#*****************************************************************************
#       Copyright (C) 2009 Martin Albrecht <malb@informatik.uni-bremen.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element cimport RingElement

from sage.libs.singular.decl cimport poly, ring

cdef int singular_polynomial_check(poly *p, ring *r) except -1
cdef int singular_polynomial_add (poly **ret, poly *p, poly *q, ring *r) noexcept
cdef int singular_polynomial_call (poly **ret, poly *p, ring *r, list args,
                                   poly *(*get_element)(object) noexcept) noexcept
cdef int singular_polynomial_cmp (poly *p, poly *q, ring *r) noexcept
cdef int singular_polynomial_rmul (poly **ret, poly *p, RingElement q, ring *r) noexcept
cdef int singular_polynomial_mul (poly **ret, poly *p, poly *q, ring *r) except -1
cdef int singular_polynomial_sub (poly **ret, poly *p, poly *q, ring *r) noexcept
cdef int singular_polynomial_div_coeff (poly **ret, poly *p, poly *q, ring *r) except -1
cdef int singular_polynomial_pow (poly **ret, poly *p, unsigned long exp, ring *r) except -1
cdef int singular_polynomial_neg(poly **ret, poly *p, ring *r) noexcept

cdef object singular_polynomial_latex(poly *p, ring *r, object base, object latex_gens)
cdef object singular_polynomial_str(poly *p, ring *r)
cdef object singular_polynomial_str_with_changed_varnames(poly *p, ring *r, object varnames)
cdef long singular_polynomial_deg(poly *p, poly *x, ring *r) noexcept

cdef int singular_polynomial_length_bounded(poly *p, int bound) noexcept
cdef int singular_vector_maximal_component(poly *v, ring *r) except -1
cdef int singular_polynomial_subst(poly **p, int var_index, poly *value, ring *r) except -1
