###############################################################################
#
#       Copyright (C) 2011 Simon King <simon.king@uni-jena.de>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#
###############################################################################

cdef class FreeAlgebra_letterplace

from sage.rings.ring cimport Algebra
from sage.structure.element cimport AlgebraElement, ModuleElement, RingElement, Element
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular, MPolynomial_libsingular
from sage.algebras.letterplace.free_algebra_element_letterplace cimport FreeAlgebraElement_letterplace


cdef class FreeAlgebra_letterplace(Algebra):
    cdef MPolynomialRing_libsingular _commutative_ring
    cdef MPolynomialRing_libsingular _current_ring
    cdef int _degbound
    cdef int __ngens
    cdef int _nb_slackvars
    cdef object __monoid
    cdef public object __custom_name
    cdef str exponents_to_string(self, E)
    cdef str exponents_to_latex(self, E)
    cdef tuple _degrees
