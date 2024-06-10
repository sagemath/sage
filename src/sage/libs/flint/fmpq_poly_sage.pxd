# distutils: libraries = flint
# distutils: depends = flint/fmpq_poly.h

#*****************************************************************************
#          Copyright (C) 2010 Sebastian Pancratz <sfp@pancratz.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.gmp.types cimport mpz_t, mpq_t
from sage.libs.flint.types cimport *
from sage.libs.flint.fmpz_vec cimport _fmpz_vec_max_limbs
from sage.libs.flint.fmpq_poly cimport fmpq_poly_numref, fmpq_poly_length


# since the fmpq_poly header seems to be lacking this inline function
cdef inline sage_fmpq_poly_max_limbs(const fmpq_poly_t poly):
    return _fmpz_vec_max_limbs(fmpq_poly_numref(poly), fmpq_poly_length(poly))

# functions removed from flint but still needed in sage
cdef void fmpq_poly_scalar_mul_mpz(fmpq_poly_t, const fmpq_poly_t, const mpz_t) noexcept
cdef void fmpq_poly_scalar_mul_mpq(fmpq_poly_t, const fmpq_poly_t, const mpq_t) noexcept
cdef void fmpq_poly_set_coeff_mpq(fmpq_poly_t, slong, const mpq_t) noexcept
cdef void fmpq_poly_get_coeff_mpq(mpq_t, const fmpq_poly_t, slong) noexcept
cdef void fmpq_poly_set_mpz(fmpq_poly_t, const mpz_t) noexcept
cdef void fmpq_poly_set_mpq(fmpq_poly_t, const mpq_t) noexcept
