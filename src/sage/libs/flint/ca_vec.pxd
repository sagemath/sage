# distutils: libraries = flint
# distutils: depends = flint/ca_vec.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    ca_ptr _ca_vec_init(slong len, ca_ctx_t ctx) noexcept
    # Returns a pointer to an array of *len* coefficients
    # initialized to zero.

    void ca_vec_init(ca_vec_t vec, slong len, ca_ctx_t ctx) noexcept
    # Initializes *vec* to a length *len* vector. All entries
    # are set to zero.

    void _ca_vec_clear(ca_ptr vec, slong len, ca_ctx_t ctx) noexcept
    # Clears all *len* entries in *vec* and frees the pointer
    # *vec* itself.

    void ca_vec_clear(ca_vec_t vec, ca_ctx_t ctx) noexcept
    # Clears the vector *vec*.

    void _ca_vec_swap(ca_ptr vec1, ca_ptr vec2, slong len, ca_ctx_t ctx) noexcept
    # Swaps the entries in *vec1* and *vec2* efficiently.

    void ca_vec_swap(ca_vec_t vec1, ca_vec_t vec2, ca_ctx_t ctx) noexcept
    # Swaps the vectors *vec1* and *vec2* efficiently.

    slong ca_vec_length(const ca_vec_t vec, ca_ctx_t ctx) noexcept
    # Returns the length of *vec*.

    void _ca_vec_fit_length(ca_vec_t vec, slong len, ca_ctx_t ctx) noexcept
    # Allocates space in *vec* for *len* elements.

    void ca_vec_set_length(ca_vec_t vec, slong len, ca_ctx_t ctx) noexcept
    # Sets the length of *vec* to *len*.
    # If *vec* is shorter on input, it will be zero-extended.
    # If *vec* is longer on input, it will be truncated.

    void _ca_vec_set(ca_ptr res, ca_srcptr src, slong len, ca_ctx_t ctx) noexcept
    # Sets *res* to a copy of *src* of length *len*.

    void ca_vec_set(ca_vec_t res, const ca_vec_t src, ca_ctx_t ctx) noexcept
    # Sets *res* to a copy of *src*.

    void _ca_vec_zero(ca_ptr res, slong len, ca_ctx_t ctx) noexcept
    # Sets the *len* entries in *res* to zeros.

    void ca_vec_zero(ca_vec_t res, slong len, ca_ctx_t ctx) noexcept
    # Sets *res* to the length *len* zero vector.

    void ca_vec_print(const ca_vec_t vec, ca_ctx_t ctx) noexcept
    # Prints *vec* to standard output. The coefficients are printed on separate lines.

    void ca_vec_printn(const ca_vec_t poly, slong digits, ca_ctx_t ctx) noexcept
    # Prints a decimal representation of *vec* with precision specified by *digits*.
    # The coefficients are comma-separated and the whole list is enclosed in square brackets.

    void ca_vec_append(ca_vec_t vec, const ca_t f, ca_ctx_t ctx) noexcept
    # Appends *f* to the end of *vec*.

    void _ca_vec_neg(ca_ptr res, ca_srcptr src, slong len, ca_ctx_t ctx) noexcept

    void ca_vec_neg(ca_vec_t res, const ca_vec_t src, ca_ctx_t ctx) noexcept
    # Sets *res* to the negation of *src*.

    void _ca_vec_add(ca_ptr res, ca_srcptr vec1, ca_srcptr vec2, slong len, ca_ctx_t ctx) noexcept

    void _ca_vec_sub(ca_ptr res, ca_srcptr vec1, ca_srcptr vec2, slong len, ca_ctx_t ctx) noexcept
    # Sets *res* to the sum or difference of *vec1* and *vec2*,
    # all vectors having length *len*.

    void _ca_vec_scalar_mul_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx) noexcept
    # Sets *res* to *src* multiplied by *c*, all vectors having
    # length *len*.

    void _ca_vec_scalar_div_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx) noexcept
    # Sets *res* to *src* divided by *c*, all vectors having
    # length *len*.

    void _ca_vec_scalar_addmul_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx) noexcept
    # Adds *src* multiplied by *c* to the vector *res*, all vectors having
    # length *len*.

    void _ca_vec_scalar_submul_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx) noexcept
    # Subtracts *src* multiplied by *c* from the vector *res*, all vectors having
    # length *len*.

    truth_t _ca_vec_check_is_zero(ca_srcptr vec, slong len, ca_ctx_t ctx) noexcept
    # Returns whether *vec* is the zero vector.

    bint _ca_vec_is_fmpq_vec(ca_srcptr vec, slong len, ca_ctx_t ctx) noexcept
    # Checks if all elements of *vec* are structurally rational numbers.

    bint _ca_vec_fmpq_vec_is_fmpz_vec(ca_srcptr vec, slong len, ca_ctx_t ctx) noexcept
    # Assuming that all elements of *vec* are structurally rational numbers,
    # checks if all elements are integers.

    void _ca_vec_fmpq_vec_get_fmpz_vec_den(fmpz * c, fmpz_t den, ca_srcptr vec, slong len, ca_ctx_t ctx) noexcept
    # Assuming that all elements of *vec* are structurally rational numbers,
    # converts them to a vector of integers *c* on a common denominator
    # *den*.

    void _ca_vec_set_fmpz_vec_div_fmpz(ca_ptr res, const fmpz * v, const fmpz_t den, slong len, ca_ctx_t ctx) noexcept
    # Sets *res* to the rational vector given by numerators *v*
    # and the common denominator *den*.
