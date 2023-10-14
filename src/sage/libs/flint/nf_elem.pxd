# distutils: libraries = flint
# distutils: depends = flint/nf_elem.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void nf_elem_init(nf_elem_t a, const nf_t nf)
    # Initialise a number field element to belong to the given number field
    # ``nf``. The element is set to zero.

    void nf_elem_clear(nf_elem_t a, const nf_t nf)
    # Clear resources allocated by the given number field element in the given
    # number field.

    void nf_elem_randtest(nf_elem_t a, flint_rand_t state, mp_bitcnt_t bits, const nf_t nf)
    # Generate a random number field element `a` in the number field ``nf``
    # whose coefficients have up to the given number of bits.

    void nf_elem_canonicalise(nf_elem_t a, const nf_t nf)
    # Canonicalise a number field element, i.e. reduce numerator and denominator
    # to lowest terms. If the numerator is `0`, set the denominator to `1`.

    void _nf_elem_reduce(nf_elem_t a, const nf_t nf)
    # Reduce a number field element modulo the defining polynomial. This is used
    # with functions such as ``nf_elem_mul_red`` which allow reduction to be
    # delayed. Does not canonicalise.

    void nf_elem_reduce(nf_elem_t a, const nf_t nf)
    # Reduce a number field element modulo the defining polynomial. This is used
    # with functions such as ``nf_elem_mul_red`` which allow reduction to be
    # delayed.

    int _nf_elem_invertible_check(nf_elem_t a, const nf_t nf)
    # Whilst the defining polynomial for a number field should by definition be
    # irreducible, it is not enforced. Thus in test code, it is convenient to be
    # able to check that a given number field element is invertible modulo the
    # defining polynomial of the number field. This function does precisely this.
    # If `a` is invertible modulo the defining polynomial of ``nf`` the value
    # `1` is returned, otherwise `0` is returned.
    # The function is only intended to be used in test code.

    void nf_elem_set_fmpz_mat_row(nf_elem_t b, const fmpz_mat_t M, const slong i, fmpz_t den, const nf_t nf)
    # Set `b` to the element specified by row `i` of the matrix `M` and with the
    # given denominator `d`. Column `0` of the matrix corresponds to the constant
    # coefficient of the number field element.

    void nf_elem_get_fmpz_mat_row(fmpz_mat_t M, const slong i, fmpz_t den, const nf_elem_t b, const nf_t nf)
    # Set the row `i` of the matrix `M` to the coefficients of the numerator of
    # the element `b` and `d` to the denominator of `b`. Column `0` of the matrix
    # corresponds to the constant coefficient of the number field element.

    void nf_elem_set_fmpq_poly(nf_elem_t a, const fmpq_poly_t pol, const nf_t nf)
    # Set `a` to the element corresponding to the polynomial ``pol``.

    void nf_elem_get_fmpq_poly(fmpq_poly_t pol, const nf_elem_t a, const nf_t nf)
    # Set ``pol`` to a polynomial corresponding to `a`, reduced modulo the
    # defining polynomial of ``nf``.

    void nf_elem_get_nmod_poly_den(nmod_poly_t pol, const nf_elem_t a, const nf_t nf, int den)
    # Set ``pol`` to the reduction of the polynomial corresponding to the
    # numerator of `a`. If ``den == 1``, the result is multiplied by the
    # inverse of the denominator of `a`. In this case it is assumed that the
    # reduction of the denominator of `a` is invertible.

    void nf_elem_get_nmod_poly(nmod_poly_t pol, const nf_elem_t a, const nf_t nf)
    # Set ``pol`` to the reduction of the polynomial corresponding to the
    # numerator of `a`. The result is multiplied by the inverse of the
    # denominator of `a`. It is assumed that the reduction of the denominator of
    # `a` is invertible.

    void nf_elem_get_fmpz_mod_poly_den(fmpz_mod_poly_t pol, const nf_elem_t a, const nf_t nf, int den, const fmpz_mod_ctx_t ctx)
    # Set ``pol`` to the reduction of the polynomial corresponding to the
    # numerator of `a`. If ``den == 1``, the result is multiplied by the
    # inverse of the denominator of `a`. In this case it is assumed that the
    # reduction of the denominator of `a` is invertible.

    void nf_elem_get_fmpz_mod_poly(fmpz_mod_poly_t pol, const nf_elem_t a, const nf_t nf, const fmpz_mod_ctx_t ctx)
    # Set ``pol`` to the reduction of the polynomial corresponding to the
    # numerator of `a`. The result is multiplied by the inverse of the
    # denominator of `a`. It is assumed that the reduction of the denominator of
    # `a` is invertible.

    void nf_elem_set_den(nf_elem_t b, fmpz_t d, const nf_t nf)
    # Set the denominator of the ``nf_elem_t b`` to the given integer `d`.
    # Assumes `d > 0`.

    void nf_elem_get_den(fmpz_t d, const nf_elem_t b, const nf_t nf)
    # Set `d` to the denominator of the ``nf_elem_t b``.

    void _nf_elem_set_coeff_num_fmpz(nf_elem_t a, slong i, const fmpz_t d, const nf_t nf)
    # Set the `i`-th coefficient of the denominator of `a` to the given integer
    # `d`.

    bint _nf_elem_equal(const nf_elem_t a, const nf_elem_t b, const nf_t nf)
    # Return `1` if the given number field elements are equal in the given
    # number field ``nf``. This function does \emph{not} assume `a` and `b`
    # are canonicalised.

    bint nf_elem_equal(const nf_elem_t a, const nf_elem_t b, const nf_t nf)
    # Return `1` if the given number field elements are equal in the given
    # number field ``nf``. This function assumes `a` and `b` \emph{are}
    # canonicalised.

    bint nf_elem_is_zero(const nf_elem_t a, const nf_t nf)
    # Return `1` if the given number field element is equal to zero,
    # otherwise return `0`.

    bint nf_elem_is_one(const nf_elem_t a, const nf_t nf)
    # Return `1` if the given number field element is equal to one,
    # otherwise return `0`.

    void nf_elem_print_pretty(const nf_elem_t a, const nf_t nf, const char * var)
    # Print the given number field element to ``stdout`` using the
    # null-terminated string ``var`` not equal to ``"\0"`` as the
    # name of the primitive element.

    void nf_elem_zero(nf_elem_t a, const nf_t nf)

    void nf_elem_one(nf_elem_t a, const nf_t nf)

    void nf_elem_set(nf_elem_t a, const nf_elem_t b, const nf_t nf)
    # Set the number field element `a` to equal the number field element `b`,
    # i.e. set `a = b`.

    void nf_elem_neg(nf_elem_t a, const nf_elem_t b, const nf_t nf)
    # Set the number field element `a` to minus the number field element `b`,
    # i.e. set `a = -b`.

    void nf_elem_swap(nf_elem_t a, nf_elem_t b, const nf_t nf)
    # Efficiently swap the two number field elements `a` and `b`.

    void nf_elem_mul_gen(nf_elem_t a, const nf_elem_t b, const nf_t nf)
    # Multiply the element `b` with the generator of the number field.

    void _nf_elem_add(nf_elem_t r, const nf_elem_t a, const nf_elem_t b, const nf_t nf)
    # Add two elements of a number field ``nf``, i.e. set `r = a + b`.
    # Canonicalisation is not performed.

    void nf_elem_add(nf_elem_t r, const nf_elem_t a, const nf_elem_t b, const nf_t nf)
    # Add two elements of a number field ``nf``, i.e. set `r = a + b`.

    void _nf_elem_sub(nf_elem_t r, const nf_elem_t a, const nf_elem_t b, const nf_t nf)
    # Subtract two elements of a number field ``nf``, i.e. set `r = a - b`.
    # Canonicalisation is not performed.

    void nf_elem_sub(nf_elem_t r, const nf_elem_t a, const nf_elem_t b, const nf_t nf)
    # Subtract two elements of a number field ``nf``, i.e. set `r = a - b`.

    void _nf_elem_mul(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)
    # Multiply two elements of a number field ``nf``, i.e. set `r = a * b`.
    # Does not canonicalise. Aliasing of inputs with output is not supported.

    void _nf_elem_mul_red(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf, int red)
    # As per ``_nf_elem_mul``, but reduction modulo the defining polynomial
    # of the number field is only carried out if ``red == 1``. Assumes both
    # inputs are reduced.

    void nf_elem_mul(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)
    # Multiply two elements of a number field ``nf``, i.e. set `r = a * b`.

    void nf_elem_mul_red(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf, int red)
    # As per ``nf_elem_mul``, but reduction modulo the defining polynomial
    # of the number field is only carried out if ``red == 1``. Assumes both
    # inputs are reduced.

    void _nf_elem_inv(nf_elem_t r, const nf_elem_t a, const nf_t nf)
    # Invert an element of a number field ``nf``, i.e. set `r = a^{-1}`.
    # Aliasing of the input with the output is not supported.

    void nf_elem_inv(nf_elem_t r, const nf_elem_t a, const nf_t nf)
    # Invert an element of a number field ``nf``, i.e. set `r = a^{-1}`.

    void _nf_elem_div(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)
    # Set `a` to `b/c` in the given number field. Aliasing of `a` and `b` is not
    # permitted.

    void nf_elem_div(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)
    # Set `a` to `b/c` in the given number field.

    void _nf_elem_pow(nf_elem_t res, const nf_elem_t a, ulong e, const nf_t nf)
    # Set ``res`` to `a^e` using left-to-right binary exponentiation as
    # described on p. 461 of [Knu1997]_.
    # Assumes that `a \neq 0` and `e > 1`. Does not support aliasing.

    void nf_elem_pow(nf_elem_t res, const nf_elem_t a, ulong e, const nf_t nf)
    # Set ``res`` = ``a^e`` using the binary exponentiation algorithm.
    # If `e` is zero, returns one, so that in particular ``0^0 = 1``.

    void _nf_elem_norm(fmpz_t rnum, fmpz_t rden, const nf_elem_t a, const nf_t nf)
    # Set ``rnum, rden`` to the absolute norm of the given number field
    # element `a`.

    void nf_elem_norm(fmpq_t res, const nf_elem_t a, const nf_t nf)
    # Set ``res`` to the absolute norm of the given number field
    # element `a`.

    void nf_elem_norm_div(fmpq_t res, const nf_elem_t a, const nf_t nf, const fmpz_t div, slong nbits)
    # Set ``res`` to the absolute norm of the given number field element `a`,
    # divided by ``div`` . Assumes the result to be an integer and having
    # at most ``nbits`` bits.

    void _nf_elem_norm_div(fmpz_t rnum, fmpz_t rden, const nf_elem_t a, const nf_t nf, const fmpz_t divisor, slong nbits)
    # Set ``rnum, rden`` to the absolute norm of the given number field element `a`,
    # divided by ``div`` . Assumes the result to be an integer and having
    # at most ``nbits`` bits.

    void _nf_elem_trace(fmpz_t rnum, fmpz_t rden, const nf_elem_t a, const nf_t nf)
    # Set ``rnum, rden`` to the absolute trace of the given number field
    # element `a`.

    void nf_elem_trace(fmpq_t res, const nf_elem_t a, const nf_t nf)
    # Set ``res`` to the absolute trace of the given number field
    # element `a`.

    void nf_elem_rep_mat(fmpq_mat_t res, const nf_elem_t a, const nf_t nf)
    # Set ``res`` to the matrix representing the multiplication with `a` with
    # respect to the basis `1, a, \dotsc, a^{d - 1}`, where `a` is the generator
    # of the number field of `d` is its degree.

    void nf_elem_rep_mat_fmpz_mat_den(fmpz_mat_t res, fmpz_t den, const nf_elem_t a, const nf_t nf)
    # Return a tuple `M, d` such that `M/d` is the matrix representing the
    # multiplication with `a` with respect to the basis `1, a, \dotsc, a^{d - 1}`,
    # where `a` is the generator of the number field of `d` is its degree.
    # The integral matrix `M` is primitive.

    void nf_elem_mod_fmpz_den(nf_elem_t z, const nf_elem_t a, const fmpz_t mod, const nf_t nf, int den)
    # If ``den == 0``, return an element `z` with denominator `1`, such that
    # the coefficients of `z - da` are divisble by ``mod``, where `d` is the
    # denominator of `a`. The coefficients of `z` are reduced modulo ``mod``.
    # If ``den == 1``, return an element `z`, such that `z - a` has
    # denominator `1` and the coefficients of `z - a` are divisible by ``mod``.
    # The coefficients of `z` are reduced modulo `\mathtt{mod} \cdot d`, where `d` is the
    # denominator of `a`.
    # Reduction takes place with respect to the positive residue system.

    void nf_elem_smod_fmpz_den(nf_elem_t z, const nf_elem_t a, const fmpz_t mod, const nf_t nf, int den)
    # If ``den == 0``, return an element `z` with denominator `1`, such that
    # the coefficients of `z - da` are divisble by ``mod``, where `d` is the
    # denominator of `a`. The coefficients of `z` are reduced modulo ``mod``.
    # If ``den == 1``, return an element `z`, such that `z - a` has
    # denominator `1` and the coefficients of `z - a` are divisible by ``mod``.
    # The coefficients of `z` are reduced modulo `\mathtt{mod} \cdot d`, where `d` is the
    # denominator of `a`.
    # Reduction takes place with respect to the symmetric residue system.

    void nf_elem_mod_fmpz(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf)
    # Return an element `z` such that `z - a` has denominator `1` and the
    # coefficients of `z - a` are divisible by ``mod``. The coefficients of
    # `z` are reduced modulo `\mathtt{mod} \cdot d`, where `d` is the denominator of `b`.
    # Reduction takes place with respect to the positive residue system.

    void nf_elem_smod_fmpz(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf)
    # Return an element `z` such that `z - a` has denominator `1` and the
    # coefficients of `z - a` are divisible by ``mod``. The coefficients of
    # `z` are reduced modulo `\mathtt{mod} \cdot d`, where `d` is the denominator of `b`.
    # Reduction takes place with respect to the symmetric residue system.

    void nf_elem_coprime_den(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf)
    # Return an element `z` such that the denominator of `z - a` is coprime to
    # ``mod``.
    # Reduction takes place with respect to the positive residue system.

    void nf_elem_coprime_den_signed(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf)
    # Return an element `z` such that the denominator of `z - a` is coprime to
    # ``mod``.
    # Reduction takes place with respect to the symmetric residue system.
