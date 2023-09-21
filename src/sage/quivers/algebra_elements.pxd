# sage_setup: distribution = sagemath-flint
"""
Cython types for elements of path algebras
"""
# ****************************************************************************
#     Copyright (C) 2015 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# This file declares the types. The implementation of the basic on these types
# is in algebra_elements.pxi, the implementation of the Python class is in
# algebra_elements.pyx. The latter file also contains all doctests.

from cpython.object cimport PyObject
from sage.data_structures.bounded_integer_sequences cimport *
from sage.structure.element cimport RingElement, ModuleElement, Element
from sage.quivers.paths cimport QuiverPath

# Type definitions

cdef struct path_mon_s:
    # The monomial is of the form "a*I_i*b" with "a" a path
    # of length "l_len", and I_i the generator of the i-th component of a free module.
    mp_size_t l_len
    # In a sub-module of a direct sum, "i=pos" denotes the direct summand that
    # this monomial belongs to. If pos==-1 then the monomial is supposed to be
    # element of a path algebra, and should be formed by a single path. In particular,
    # l_len has to be zero.
    long pos
    # In the Schreyer order, monomials of the form a*I_i*b are not
    # compared directly, but to each position is associated a monomial s_i,
    # and a*I_i*b is compared with c*I_j*d by first comparing a*s_i*b with
    # c*s_j*d. s_length is the length of s_i.
    mp_size_t s_len
    # paths are encoded as lists of integers. We store a*s_i*b if the monomial is
    # a*I_i*b.
    biseq_t path

ctypedef path_mon_s path_mon_t[1]

cdef struct path_term_t:
    # A term is given by a monomial "mon" and a coefficient "coef".
    path_mon_t mon
    # We need to manually take care of the reference count for the
    # coefficient!
    PyObject *coef
    # In a polynomial, the terms are arranged in a pointered list
    path_term_t *nxt

# Type of monomial ordering functions.
# Returns -1, 0 or 1, depending on whether the first argument is
# smaller (wrt. the chosen ordering function), equal to, or greater
# than the second argument.
ctypedef int (*path_order_t)(path_mon_t, path_mon_t) except -2

# Polynomials are decreasingly sorted lists of terms, wrt. some fixed
# monomial ordering. The list starts with pointer .lead. For convenience,
# the number of terms is directly available by .nterms.
cdef struct path_poly_t:
    path_term_t *lead
    size_t nterms

# In path_poly_t, the terms need not to have all the same start and end
# points. path_homog_poly_t points to a path_poly_t whose terms are all
# guaranteed to start and end at the given vertex labels (which are integers).
# We will work with lists of path_homog_poly_t, and thus have a pointer
# to the next start and end point homogeneous polynomial.
cdef struct path_homog_poly_t:
    path_poly_t *poly
    int start, end
    path_homog_poly_t *nxt

cdef class PathAlgebraElement(RingElement):
    # The terms of path algebra element are stored as a list of start and
    # end point homogeneous polynomials. These are sorted increasingly
    # according to the (start, end) point pairs.
    cdef path_homog_poly_t *data
    # A fixed term ordering is stored along with the path algebra element.
    # This ordering has to be passed as an argument to all boilerplate
    # functions.
    cdef path_order_t cmp_terms
    cdef long _hash
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef ssize_t degree(self) except -2
    cpdef dict monomial_coefficients(self)
    cpdef list coefficients(self)
    cpdef list monomials(self)
    cpdef list support(self)
    cpdef list terms(self)
    cpdef object coefficient(self, QuiverPath P)
    cdef list _sorted_items_for_printing(self)
    cdef inline PathAlgebraElement _new_(self, path_homog_poly_t *h)
