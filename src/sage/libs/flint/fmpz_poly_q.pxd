# distutils: libraries = flint
# distutils: depends = flint/fmpz_poly_q.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fmpz_poly_q_init(fmpz_poly_q_t rop)
    # Initialises ``rop``.

    void fmpz_poly_q_clear(fmpz_poly_q_t rop)
    # Clears the object ``rop``.

    fmpz_poly_struct * fmpz_poly_q_numref(const fmpz_poly_q_t op)
    # Returns a reference to the numerator of ``op``.

    fmpz_poly_struct * fmpz_poly_q_denref(const fmpz_poly_q_t op)
    # Returns a reference to the denominator of ``op``.

    void fmpz_poly_q_canonicalise(fmpz_poly_q_t rop)
    # Brings ``rop`` into canonical form, only assuming that
    # the denominator is non-zero.

    int fmpz_poly_q_is_canonical(const fmpz_poly_q_t op)
    # Checks whether the rational function ``op`` is in
    # canonical form.

    void fmpz_poly_q_randtest(fmpz_poly_q_t poly, flint_rand_t state, slong len1, flint_bitcnt_t bits1, slong len2, flint_bitcnt_t bits2)
    # Sets ``poly`` to a random rational function.

    void fmpz_poly_q_randtest_not_zero(fmpz_poly_q_t poly, flint_rand_t state, slong len1, flint_bitcnt_t bits1, slong len2, flint_bitcnt_t bits2)
    # Sets ``poly`` to a random non-zero rational function.

    void fmpz_poly_q_set(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
    # Sets the element ``rop`` to the same value as the element ``op``.

    void fmpz_poly_q_set_si(fmpz_poly_q_t rop, slong op)
    # Sets the element ``rop`` to the value given by the ``slong``
    # ``op``.

    void fmpz_poly_q_swap(fmpz_poly_q_t op1, fmpz_poly_q_t op2)
    # Swaps the elements ``op1`` and ``op2``.
    # This is done efficiently by swapping pointers.

    void fmpz_poly_q_zero(fmpz_poly_q_t rop)
    # Sets ``rop`` to zero.

    void fmpz_poly_q_one(fmpz_poly_q_t rop)
    # Sets ``rop`` to one.

    void fmpz_poly_q_neg(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
    # Sets the element ``rop`` to the additive inverse of ``op``.

    void fmpz_poly_q_inv(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
    # Sets the element ``rop`` to the multiplicative inverse of ``op``.
    # Assumes that the element ``op`` is non-zero.

    int fmpz_poly_q_is_zero(const fmpz_poly_q_t op)
    # Returns whether the element ``op`` is zero.

    int fmpz_poly_q_is_one(const fmpz_poly_q_t op)
    # Returns whether the element ``rop`` is equal to the constant
    # polynomial `1`.

    int fmpz_poly_q_equal(const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
    # Returns whether the two elements ``op1`` and ``op2`` are equal.

    void fmpz_poly_q_add(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
    # Sets ``rop`` to the sum of ``op1`` and ``op2``.

    void fmpz_poly_q_sub(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
    # Sets ``rop`` to the difference of ``op1`` and ``op2``.

    void fmpz_poly_q_addmul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
    # Adds the product of ``op1`` and ``op2`` to ``rop``.

    void fmpz_poly_q_submul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
    # Subtracts the product of ``op1`` and ``op2`` from ``rop``.

    void fmpz_poly_q_scalar_mul_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, slong x)
    # Sets ``rop`` to the product of the rational function ``op``
    # and the ``slong`` integer `x`.

    void fmpz_poly_q_scalar_mul_fmpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpz_t x)
    # Sets ``rop`` to the product of the rational function ``op``
    # and the ``fmpz_t`` integer `x`.

    void fmpz_poly_q_scalar_mul_fmpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpq_t x)
    # Sets ``rop`` to the product of the rational function ``op``
    # and the ``fmpq_t`` rational `x`.

    void fmpz_poly_q_scalar_div_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, slong x)
    # Sets ``rop`` to the quotient of the rational function ``op``
    # and the ``slong`` integer `x`.

    void fmpz_poly_q_scalar_div_fmpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpz_t x)
    # Sets ``rop`` to the quotient of the rational function ``op``
    # and the ``fmpz_t`` integer `x`.

    void fmpz_poly_q_scalar_div_fmpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpq_t x)
    # Sets ``rop`` to the quotient of the rational function ``op``
    # and the ``fmpq_t`` rational `x`.

    void fmpz_poly_q_mul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
    # Sets ``rop`` to the product of ``op1`` and ``op2``.

    void fmpz_poly_q_div(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
    # Sets ``rop`` to the quotient of ``op1`` and ``op2``.

    void fmpz_poly_q_pow(fmpz_poly_q_t rop, const fmpz_poly_q_t op, ulong exp)
    # Sets ``rop`` to the ``exp``-th power of ``op``.
    # The corner case of ``exp == 0`` is handled by setting ``rop`` to
    # the constant function `1`.  Note that this includes the case `0^0 = 1`.

    void fmpz_poly_q_derivative(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
    # Sets ``rop`` to the derivative of ``op``.

    int fmpz_poly_q_evaluate_fmpq(fmpq_t rop, const fmpz_poly_q_t f, const fmpq_t a)
    # Sets ``rop`` to `f` evaluated at the rational `a`.
    # If the denominator evaluates to zero at `a`, returns non-zero and
    # does not modify any of the variables.  Otherwise, returns `0` and
    # sets ``rop`` to the rational `f(a)`.

    int fmpz_poly_q_set_str(fmpz_poly_q_t rop, const char *s)
    # Sets ``rop`` to the rational function given
    # by the string ``s``.

    char * fmpz_poly_q_get_str(const fmpz_poly_q_t op)
    # Returns the string representation of
    # the rational function ``op``.

    char * fmpz_poly_q_get_str_pretty(const fmpz_poly_q_t op, const char *x)
    # Returns the pretty string representation of
    # the rational function ``op``.

    int fmpz_poly_q_print(const fmpz_poly_q_t op)
    # Prints the representation of the rational
    # function ``op`` to ``stdout``.

    int fmpz_poly_q_print_pretty(const fmpz_poly_q_t op, const char *x)
    # Prints the pretty representation of the rational
    # function ``op`` to ``stdout``.
