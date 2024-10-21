# sage_setup: distribution = sagemath-singular
###############################################################################
#
#       Copyright (C) 2011 Simon King <simon.king@uni-jena.de>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
#
###############################################################################

from sage.structure.parent cimport Parent
from sage.structure.element cimport AlgebraElement, ModuleElement, RingElement, Element
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular, MPolynomial_libsingular
from sage.algebras.letterplace.free_algebra_element_letterplace cimport FreeAlgebraElement_letterplace
from sage.libs.singular.decl cimport ring


cdef class FreeAlgebra_letterplace_libsingular():
    cdef ring* _lp_ring
    cdef MPolynomialRing_libsingular _commutative_ring
    cdef MPolynomialRing_libsingular _lp_ring_internal
    cdef object _ngens


cdef class FreeAlgebra_letterplace(Parent):
    cdef MPolynomialRing_libsingular _commutative_ring
    cdef MPolynomialRing_libsingular _current_ring
    cdef int _degbound
    cdef int _ngens
    cdef int _nb_slackvars
    cdef object __monoid
    cdef str exponents_to_string(self, E)
    cdef str exponents_to_latex(self, E)
    cdef tuple _degrees
    cdef public object _latex_names
