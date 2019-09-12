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

# flint/fmpq_poly.h
cdef extern from "flint_wrap.h":
    # Memory management
    void fmpq_poly_init(fmpq_poly_t)

    void fmpq_poly_init2(fmpq_poly_t, slong)
    void fmpq_poly_realloc(fmpq_poly_t, slong)

    void fmpq_poly_fit_length(fmpq_poly_t, slong)

    void fmpq_poly_clear(fmpq_poly_t)

    void fmpq_poly_canonicalise(fmpq_poly_t)
    int fmpq_poly_is_canonical(const fmpq_poly_t)

    # Polynomial parameters
    slong fmpq_poly_degree(const fmpq_poly_t)
    ulong fmpq_poly_length(const fmpq_poly_t)

    # Accessing the numerator and denominator
    fmpz *fmpq_poly_numref(fmpq_poly_t)
    fmpz *fmpq_poly_denref(fmpq_poly_t)

    void fmpq_poly_get_numerator(fmpz_poly_t, const fmpq_poly_t)

    # Assignment, swap, negation
    void fmpq_poly_set(fmpq_poly_t, const fmpq_poly_t)
    void fmpq_poly_set_si(fmpq_poly_t, slong)
    void fmpq_poly_set_ui(fmpq_poly_t, ulong)
    void fmpq_poly_set_fmpz(fmpq_poly_t, const fmpz_t)
    void fmpq_poly_set_fmpq(fmpq_poly_t, const fmpq_t)
    void fmpq_poly_set_mpz(fmpq_poly_t, const mpz_t)
    void fmpq_poly_set_mpq(fmpq_poly_t, const mpq_t)
    void fmpq_poly_set_fmpz_poly(fmpq_poly_t, const fmpz_poly_t)
    void fmpq_poly_set_array_mpq(fmpq_poly_t, const mpq_t *, slong)

    void fmpq_poly_set_str(fmpq_poly_t, const char *)
    char *fmpq_poly_get_str(const fmpq_poly_t)
    char *fmpq_poly_get_str_pretty(const fmpq_poly_t, const char *)

    void fmpq_poly_zero(fmpq_poly_t)
    void fmpq_poly_one(fmpq_poly_t)

    void fmpq_poly_neg(fmpq_poly_t, const fmpq_poly_t)
    void fmpq_poly_inv(fmpq_poly_t, const fmpq_poly_t)

    void fmpq_poly_swap(fmpq_poly_t, fmpq_poly_t)
    void fmpq_poly_truncate(fmpq_poly_t, slong)
    void fmpq_poly_get_slice(fmpq_poly_t, const fmpq_poly_t, slong, slong)
    void fmpq_poly_reverse(fmpq_poly_t, const fmpq_poly_t, slong)

    void fmpq_poly_get_coeff_fmpq(fmpq_t, const fmpq_poly_t, slong)
    void fmpq_poly_get_coeff_mpq(mpq_t, const fmpq_poly_t, slong)
    void fmpq_poly_get_coeff_si(slong, const fmpq_poly_t, slong)
    void fmpq_poly_get_coeff_ui(ulong, const fmpq_poly_t, slong)

    void fmpq_poly_set_coeff_si(fmpq_poly_t, slong, slong)
    void fmpq_poly_set_coeff_ui(fmpq_poly_t, slong, ulong)
    void fmpq_poly_set_coeff_fmpz(fmpq_poly_t, slong, const fmpz_t)
    void fmpq_poly_set_coeff_fmpq(fmpq_poly_t, slong, const fmpq_t)
    void fmpq_poly_set_coeff_mpz(fmpq_poly_t, slong, const mpz_t)
    void fmpq_poly_set_coeff_mpq(fmpq_poly_t, slong, const mpq_t)

    # Comparison
    int fmpq_poly_equal(const fmpq_poly_t, const fmpq_poly_t)
    int fmpq_poly_cmp(const fmpq_poly_t, const fmpq_poly_t)
    int fmpq_poly_is_one(const fmpq_poly_t)
    int fmpq_poly_is_zero(const fmpq_poly_t)

    # Addition and subtraction
    void fmpq_poly_add(fmpq_poly_t, const fmpq_poly_t, const fmpq_poly_t)
    void fmpq_poly_sub(fmpq_poly_t, const fmpq_poly_t, const fmpq_poly_t)

    void fmpq_poly_add_can(
            fmpq_poly_t, const fmpq_poly_t, const fmpq_poly_t, int)
    void fmpq_poly_sub_can(
            fmpq_poly_t, const fmpq_poly_t, const fmpq_poly_t, int)

    # Scalar multiplication and division
    void fmpq_poly_scalar_mul_si(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_scalar_mul_ui(fmpq_poly_t, const fmpq_poly_t, ulong)
    void fmpq_poly_scalar_mul_fmpz(
            fmpq_poly_t, const fmpq_poly_t, const fmpz_t)
    void fmpq_poly_scalar_mul_fmpq(
            fmpq_poly_t, const fmpq_poly_t, const fmpq_t)
    void fmpq_poly_scalar_mul_mpz(fmpq_poly_t, const fmpq_poly_t, const mpz_t)
    void fmpq_poly_scalar_mul_mpq(fmpq_poly_t, const fmpq_poly_t, const mpq_t)

    void fmpq_poly_scalar_div_si(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_scalar_div_ui(fmpq_poly_t, const fmpq_poly_t, ulong)
    void fmpq_poly_scalar_div_fmpz(
            fmpq_poly_t, const fmpq_poly_t, const fmpz_t)
    void fmpq_poly_scalar_div_fmpq(
            fmpq_poly_t, const fmpq_poly_t, const fmpq_t)
    void fmpq_poly_scalar_div_mpz(fmpq_poly_t, const fmpq_poly_t, const mpz_t)
    void fmpq_poly_scalar_div_mpq(fmpq_poly_t, const fmpq_poly_t, const mpq_t)

    # Multiplication
    void fmpq_poly_mul(fmpq_poly_t, const fmpq_poly_t, const fmpq_poly_t)
    void fmpq_poly_mullow(fmpq_poly_t, const fmpq_poly_t, const fmpq_poly_t, slong)

    void fmpq_poly_addmul(fmpq_poly_t, const fmpq_poly_t, const fmpq_poly_t)
    void fmpq_poly_submul(fmpq_poly_t, const fmpq_poly_t, const fmpq_poly_t)

    # Powering
    void fmpq_poly_pow(fmpq_poly_t, const fmpq_poly_t, ulong)

    # Shifting
    void fmpq_poly_shift_left(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_shift_right(fmpq_poly_t, const fmpq_poly_t, slong)

    # Euclidean division
    void fmpq_poly_divrem(
            fmpq_poly_t, fmpq_poly_t, const fmpq_poly_t, const fmpq_poly_t)
    void fmpq_poly_div(fmpq_poly_t, const fmpq_poly_t, const fmpq_poly_t)
    void fmpq_poly_rem(fmpq_poly_t, const fmpq_poly_t, const fmpq_poly_t)

    # Greatest common divisor
    void fmpq_poly_gcd(fmpq_poly_t, const fmpq_poly_t, const fmpq_poly_t)
    void fmpq_poly_xgcd(
            fmpq_poly_t, fmpq_poly_t, fmpq_poly_t,
                    const fmpq_poly_t, const fmpq_poly_t)

    void fmpq_poly_lcm(fmpq_poly_t, const fmpq_poly_t, const fmpq_poly_t)

    void fmpq_poly_resultant(fmpq_t, const fmpq_poly_t, const fmpq_poly_t)

    # Power series division
    void fmpq_poly_inv_series_newton(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_inv_series(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_div_series(
            fmpq_poly_t, const fmpq_poly_t, const fmpq_poly_t, slong)

    # Derivative and integral
    void fmpq_poly_derivative(fmpq_poly_t, const fmpq_poly_t)
    void fmpq_poly_integral(fmpq_poly_t, const fmpq_poly_t)

    # Evaluation
    void fmpq_poly_evaluate_fmpz(fmpq_t, const fmpq_poly_t, const fmpz_t)
    void fmpq_poly_evaluate_fmpq(fmpq_t, const fmpq_poly_t, const fmpq_t)
    void fmpq_poly_evaluate_mpz(mpq_t, const fmpq_poly_t, const mpz_t)
    void fmpq_poly_evaluate_mpq(mpq_t, const fmpq_poly_t, const mpq_t)

    # Composition
    void fmpq_poly_compose(fmpq_poly_t, const fmpq_poly_t, const fmpq_poly_t)
    void fmpq_poly_rescale(fmpq_poly_t, const fmpq_poly_t, const fmpq_t)

    # Revert
    void fmpq_poly_revert_series(fmpq_poly_t, fmpq_poly_t, unsigned long)

    # Gaussian content
    void fmpq_poly_content(fmpq_t, const fmpq_poly_t)
    void fmpq_poly_primitive_part(fmpq_poly_t, const fmpq_poly_t)

    int fmpq_poly_is_monic(const fmpq_poly_t)
    void fmpq_poly_make_monic(fmpq_poly_t, const fmpq_poly_t)

    # Transcendental functions
    void fmpq_poly_log_series(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_exp_series(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_atan_series(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_atanh_series(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_asin_series(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_asinh_series(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_tan_series(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_sin_series(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_cos_series(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_sinh_series(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_cosh_series(fmpq_poly_t, const fmpq_poly_t, slong)
    void fmpq_poly_tanh_series(fmpq_poly_t, const fmpq_poly_t, slong)

# since the fmpq_poly header seems to be lacking this inline function
cdef inline sage_fmpq_poly_max_limbs(const fmpq_poly_t poly):
    return _fmpz_vec_max_limbs(fmpq_poly_numref(poly), fmpq_poly_length(poly))
