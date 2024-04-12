# sage_setup: distribution = sagemath-flint
###############################################################################
#          Copyright (C) 2010 Sebastian Pancratz <sfp@pancratz.org>           #
#                                                                             #
#     Distributed under the terms of the GNU General Public License (GPL)     #
#                                                                             #
#                        http://www.gnu.org/licenses/                         #
###############################################################################

from sage.libs.flint.types cimport fmpq_poly_t

from sage.rings.polynomial.polynomial_element cimport Polynomial

cdef class Polynomial_rational_flint(Polynomial):
    cdef fmpq_poly_t _poly

    cdef Polynomial_rational_flint _new(self)
    cpdef _mod_(self, right)
    cpdef _unsafe_mutate(self, unsigned long n, value)
    cpdef Polynomial truncate(self, long n)

