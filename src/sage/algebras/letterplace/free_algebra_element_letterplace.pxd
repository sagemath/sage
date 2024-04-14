# sage_setup: distribution = sagemath-singular
###############################################################################
#
#       Copyright (C) 2011 Simon King <simon.king@uni-jena.de>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#
###############################################################################

from sage.structure.element cimport AlgebraElement, ModuleElement, Element
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular, MPolynomial_libsingular
from sage.algebras.letterplace.free_algebra_letterplace cimport FreeAlgebra_letterplace

cdef class FreeAlgebraElement_letterplace(AlgebraElement):
    cdef MPolynomial_libsingular _poly
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
