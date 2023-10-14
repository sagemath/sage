# distutils: libraries = flint
# distutils: depends = flint/fq_nmod_poly.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fq_nmod_poly_init(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    # Initialises ``poly`` for use, with context ctx, and setting its
    # length to zero. A corresponding call to :func:`fq_nmod_poly_clear`
    # must be made after finishing with the ``fq_nmod_poly_t`` to free the
    # memory used by the polynomial.

    void fq_nmod_poly_init2(fq_nmod_poly_t poly, slong alloc, const fq_nmod_ctx_t ctx)
    # Initialises ``poly`` with space for at least ``alloc``
    # coefficients and sets the length to zero.  The allocated
    # coefficients are all set to zero.  A corresponding call to
    # :func:`fq_nmod_poly_clear` must be made after finishing with the
    # ``fq_nmod_poly_t`` to free the memory used by the polynomial.

    void fq_nmod_poly_realloc(fq_nmod_poly_t poly, slong alloc, const fq_nmod_ctx_t ctx)
    # Reallocates the given polynomial to have space for ``alloc``
    # coefficients.  If ``alloc`` is zero the polynomial is cleared
    # and then reinitialised.  If the current length is greater than
    # ``alloc`` the polynomial is first truncated to length
    # ``alloc``.

    void fq_nmod_poly_fit_length(fq_nmod_poly_t poly, slong len, const fq_nmod_ctx_t ctx)
    # If ``len`` is greater than the number of coefficients currently
    # allocated, then the polynomial is reallocated to have space for at
    # least ``len`` coefficients.  No data is lost when calling this
    # function.
    # The function efficiently deals with the case where
    # ``fit_length`` is called many times in small increments by at
    # least doubling the number of allocated coefficients when length is
    # larger than the number of coefficients currently allocated.

    void _fq_nmod_poly_set_length(fq_nmod_poly_t poly, slong newlen, const fq_nmod_ctx_t ctx)
    # Sets the coefficients of ``poly`` beyond ``len`` to zero and
    # sets the length of ``poly`` to ``len``.

    void fq_nmod_poly_clear(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    # Clears the given polynomial, releasing any memory used.  It must
    # be reinitialised in order to be used again.

    void _fq_nmod_poly_normalise(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    # Sets the length of ``poly`` so that the top coefficient is
    # non-zero.  If all coefficients are zero, the length is set to
    # zero.  This function is mainly used internally, as all functions
    # guarantee normalisation.

    void _fq_nmod_poly_normalise2(const fq_nmod_struct *poly, slong *length, const fq_nmod_ctx_t ctx)
    # Sets the length ``length`` of ``(poly,length)`` so that the
    # top coefficient is non-zero. If all coefficients are zero, the
    # length is set to zero. This function is mainly used internally, as
    # all functions guarantee normalisation.

    void fq_nmod_poly_truncate(fq_nmod_poly_t poly, slong newlen, const fq_nmod_ctx_t ctx)
    # Truncates the polynomial to length at most ``n``.

    void fq_nmod_poly_set_trunc(fq_nmod_poly_t poly1, fq_nmod_poly_t poly2, slong newlen, const fq_nmod_ctx_t ctx)
    # Sets ``poly1`` to ``poly2`` truncated to length `n`.

    void _fq_nmod_poly_reverse(fq_nmod_struct* output, const fq_nmod_struct* input, slong len, slong m, const fq_nmod_ctx_t ctx)
    # Sets ``output`` to the reverse of ``input``, which is of
    # length ``len``, but thinking of it as a polynomial of
    # length ``m``, notionally zero-padded if necessary. The
    # length ``m`` must be non-negative, but there are no other
    # restrictions. The polynomial ``output`` must have space for
    # ``m`` coefficients.

    void fq_nmod_poly_reverse(fq_nmod_poly_t output, const fq_nmod_poly_t input, slong m, const fq_nmod_ctx_t ctx)
    # Sets ``output`` to the reverse of ``input``, thinking of it
    # as a polynomial of length ``m``, notionally zero-padded if
    # necessary).  The length ``m`` must be non-negative, but there
    # are no other restrictions. The output polynomial will be set to
    # length ``m`` and then normalised.

    slong fq_nmod_poly_degree(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    # Returns the degree of the polynomial ``poly``.

    slong fq_nmod_poly_length(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    # Returns the length of the polynomial ``poly``.

    fq_nmod_struct * fq_nmod_poly_lead(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    # Returns a pointer to the leading coefficient of ``poly``, or
    # ``NULL`` if ``poly`` is the zero polynomial.

    void fq_nmod_poly_randtest(fq_nmod_poly_t f, flint_rand_t state, slong len, const fq_nmod_ctx_t ctx)
    # Sets `f` to a random polynomial of length at most ``len``
    # with entries in the field described by ``ctx``.

    void fq_nmod_poly_randtest_not_zero(fq_nmod_poly_t f, flint_rand_t state, slong len, const fq_nmod_ctx_t ctx)
    # Same as ``fq_nmod_poly_randtest`` but guarantees that the polynomial
    # is not zero.

    void fq_nmod_poly_randtest_monic(fq_nmod_poly_t f, flint_rand_t state, slong len, const fq_nmod_ctx_t ctx)
    # Sets `f` to a random monic polynomial of length ``len`` with
    # entries in the field described by ``ctx``.

    void fq_nmod_poly_randtest_irreducible(fq_nmod_poly_t f, flint_rand_t state, slong len, const fq_nmod_ctx_t ctx)
    # Sets `f` to a random monic, irreducible polynomial of length
    # ``len`` with entries in the field described by ``ctx``.

    void _fq_nmod_poly_set(fq_nmod_struct *rop, const fq_nmod_struct *op, slong len, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, len``) to ``(op, len)``.

    void fq_nmod_poly_set(fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_ctx_t ctx)
    # Sets the polynomial ``poly1`` to the polynomial ``poly2``.

    void fq_nmod_poly_set_fq_nmod(fq_nmod_poly_t poly, const fq_nmod_t c, const fq_nmod_ctx_t ctx)
    # Sets the polynomial ``poly`` to ``c``.

    void fq_nmod_poly_set_fmpz_mod_poly(fq_nmod_poly_t rop, const fmpz_mod_poly_t op, const fq_nmod_ctx_t ctx)
    # Sets the polynomial ``rop`` to the polynomial ``op``

    void fq_nmod_poly_set_nmod_poly(fq_nmod_poly_t rop, const nmod_poly_t op, const fq_nmod_ctx_t ctx)
    # Sets the polynomial ``rop`` to the polynomial ``op``

    void fq_nmod_poly_swap(fq_nmod_poly_t op1, fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    # Swaps the two polynomials ``op1`` and ``op2``.

    void _fq_nmod_poly_zero(fq_nmod_struct *rop, slong len, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, len)`` to the zero polynomial.

    void fq_nmod_poly_zero(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    # Sets ``poly`` to the zero polynomial.

    void fq_nmod_poly_one(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    # Sets ``poly`` to the constant polynomial `1`.

    void fq_nmod_poly_gen(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    # Sets ``poly`` to the polynomial `x`.

    void fq_nmod_poly_make_monic(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to ``op``, normed to have leading coefficient 1.

    void _fq_nmod_poly_make_monic(fq_nmod_struct *rop, const fq_nmod_struct *op, slong length, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to ``(op,length)``, normed to have leading coefficient 1.
    # Assumes that ``rop`` has enough space for the polynomial, assumes that
    # ``op`` is not zero (and thus has an invertible leading coefficient).

    void fq_nmod_poly_get_coeff(fq_nmod_t x, const fq_nmod_poly_t poly, slong n, const fq_nmod_ctx_t ctx)
    # Sets `x` to the coefficient of `X^n` in ``poly``.

    void fq_nmod_poly_set_coeff(fq_nmod_poly_t poly, slong n, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    # Sets the coefficient of `X^n` in ``poly`` to `x`.

    void fq_nmod_poly_set_coeff_fmpz(fq_nmod_poly_t poly, slong n, const fmpz_t x, const fq_nmod_ctx_t ctx)
    # Sets the coefficient of `X^n` in the polynomial to `x`,
    # assuming `n \geq 0`.

    int fq_nmod_poly_equal(const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_ctx_t ctx)
    # Returns nonzero if the two polynomials ``poly1`` and ``poly2``
    # are equal, otherwise return zero.

    int fq_nmod_poly_equal_trunc(const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, slong n, const fq_nmod_ctx_t ctx)
    # Notionally truncate ``poly1`` and ``poly2`` to length `n` and
    # return nonzero if they are equal, otherwise return zero.

    int fq_nmod_poly_is_zero(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    # Returns whether the polynomial ``poly`` is the zero polynomial.

    int fq_nmod_poly_is_one(const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    # Returns whether the polynomial ``poly`` is equal
    # to the constant polynomial `1`.

    int fq_nmod_poly_is_gen(const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    # Returns whether the polynomial ``poly`` is equal
    # to the polynomial `x`.

    int fq_nmod_poly_is_unit(const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    # Returns whether the polynomial ``poly`` is a unit in the polynomial
    # ring `\mathbf{F}_q[X]`, i.e. if it has degree `0` and is non-zero.

    int fq_nmod_poly_equal_fq_nmod(const fq_nmod_poly_t poly, const fq_nmod_t c, const fq_nmod_ctx_t ctx)
    # Returns whether the polynomial ``poly`` is equal the (constant)
    # `\mathbf{F}_q` element ``c``

    void _fq_nmod_poly_add(fq_nmod_struct *res, const fq_nmod_struct *poly1, slong len1, const fq_nmod_struct *poly2, slong len2, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the sum of ``(poly1,len1)`` and ``(poly2,len2)``.

    void fq_nmod_poly_add(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the sum of ``poly1`` and ``poly2``.

    void fq_nmod_poly_add_si(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, slong c, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the sum of ``poly1`` and ``c``.

    void fq_nmod_poly_add_series(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, slong n, const fq_nmod_ctx_t ctx)
    # Notionally truncate ``poly1`` and ``poly2`` to length ``n`` and set
    # ``res`` to the sum.

    void _fq_nmod_poly_sub(fq_nmod_struct *res, const fq_nmod_struct *poly1, slong len1, const fq_nmod_struct *poly2, slong len2, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the difference of ``(poly1,len1)`` and
    # ``(poly2,len2)``.

    void fq_nmod_poly_sub(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the difference of ``poly1`` and ``poly2``.

    void fq_nmod_poly_sub_series(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, slong n, const fq_nmod_ctx_t ctx)
    # Notionally truncate ``poly1`` and ``poly2`` to length ``n`` and set
    # ``res`` to the difference.

    void _fq_nmod_poly_neg(fq_nmod_struct *rop, const fq_nmod_struct *op, slong len, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the additive inverse of ``(poly,len)``.

    void fq_nmod_poly_neg(fq_nmod_poly_t res, const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the additive inverse of ``poly``.

    void _fq_nmod_poly_scalar_mul_fq_nmod(fq_nmod_struct *rop, const fq_nmod_struct *op, slong len, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    # Sets ``(rop,len)`` to the product of ``(op,len)`` by the
    # scalar ``x``, in the context defined by ``ctx``.

    void fq_nmod_poly_scalar_mul_fq_nmod(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the product of ``op`` by the scalar ``x``, in the context
    # defined by ``ctx``.

    void _fq_nmod_poly_scalar_addmul_fq_nmod(fq_nmod_struct *rop, const fq_nmod_struct *op, slong len, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    # Adds to ``(rop,len)`` the product of ``(op,len)`` by the
    # scalar ``x``, in the context defined by ``ctx``.
    # In particular, assumes the same length for ``op`` and
    # ``rop``.

    void fq_nmod_poly_scalar_addmul_fq_nmod(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    # Adds to ``rop`` the product of ``op`` by the
    # scalar ``x``, in the context defined by ``ctx``.

    void _fq_nmod_poly_scalar_submul_fq_nmod(fq_nmod_struct *rop, const fq_nmod_struct *op, slong len, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    # Subtracts from ``(rop,len)`` the product of ``(op,len)`` by the
    # scalar ``x``, in the context defined by ``ctx``.
    # In particular, assumes the same length for ``op`` and
    # ``rop``.

    void fq_nmod_poly_scalar_submul_fq_nmod(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    # Subtracts from ``rop`` the product of ``op`` by the
    # scalar ``x``, in the context defined by ``ctx``.

    void _fq_nmod_poly_scalar_div_fq(fq_nmod_struct *rop, const fq_nmod_struct *op, slong len, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    # Sets ``(rop,len)`` to the quotient of ``(op,len)`` by the
    # scalar ``x``, in the context defined by ``ctx``. An exception is raised
    # if ``x`` is zero.

    void fq_nmod_poly_scalar_div_fq(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the quotient of ``op`` by the scalar ``x``, in the context
    # defined by ``ctx``. An exception is raised if ``x`` is zero.

    void _fq_nmod_poly_mul_classical(fq_nmod_struct *rop, const fq_nmod_struct *op1, slong len1, const fq_nmod_struct *op2, slong len2, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, len1 + len2 - 1)`` to the product of ``(op1, len1)``
    # and ``(op2, len2)``, assuming that ``len1`` is at least ``len2``
    # and neither is zero.
    # Permits zero padding.  Does not support aliasing of ``rop``
    # with either ``op1`` or ``op2``.

    void fq_nmod_poly_mul_classical(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the product of ``op1`` and ``op2``
    # using classical polynomial multiplication.

    void _fq_nmod_poly_mul_reorder(fq_nmod_struct *rop, const fq_nmod_struct *op1, slong len1, const fq_nmod_struct *op2, slong len2, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, len1 + len2 - 1)`` to the product of ``(op1, len1)``
    # and ``(op2, len2)``, assuming that ``len1`` and ``len2`` are
    # non-zero.
    # Permits zero padding.  Supports aliasing.

    void fq_nmod_poly_mul_reorder(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the product of ``op1`` and ``op2``,
    # reordering the two indeterminates `X` and `Y` when viewing
    # the polynomials as elements of `\mathbf{F}_p[X,Y]`.
    # Suppose `\mathbf{F}_q = \mathbf{F}_p[X]/ (f(X))` and recall
    # that elements of `\mathbf{F}_q` are internally represented
    # by elements of type ``fmpz_poly``.  For small degree extensions
    # but polynomials in `\mathbf{F}_q[Y]` of large degree `n`, we
    # change the representation to
    # .. math ::
    # \begin{split}
    # g(Y) & = \sum_{i=0}^{n} a_i(X) Y^i \\
    # & = \sum_{j=0}^{d} \sum_{i=0}^{n} \text{Coeff}(a_i(X), j) Y^i.
    # \end{split}
    # This allows us to use a poor algorithm (such as classical multiplication)
    # in the `X`-direction and leverage the existing fast integer
    # multiplication routines in the `Y`-direction where the polynomial
    # degree `n` is large.

    void _fq_nmod_poly_mul_univariate(fq_nmod_struct *rop, const fq_nmod_struct *op1, slong len1, const fq_nmod_struct *op2, slong len2, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, len1 + len2 - 1)`` to the product of ``(op1, len1)``
    # and ``(op2, len2)``.
    # Permits zero padding and makes no assumptions on ``len1`` and ``len2``.
    # Supports aliasing.

    void fq_nmod_poly_mul_univariate(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the product of ``op1`` and ``op2``
    # using a bivariate to univariate transformation and reducing
    # this problem to multiplying two univariate polynomials.

    void _fq_nmod_poly_mul_KS(fq_nmod_struct *rop, const fq_nmod_struct *op1, slong len1, const fq_nmod_struct *op2, slong len2, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, len1 + len2 - 1)`` to the product of ``(op1, len1)``
    # and ``(op2, len2)``.
    # Permits zero padding and places no assumptions on the
    # lengths ``len1`` and ``len2``.  Supports aliasing.

    void fq_nmod_poly_mul_KS(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the product of ``op1`` and ``op2``
    # using Kronecker substitution, that is, by encoding each
    # coefficient in `\mathbf{F}_{q}` as an integer and reducing
    # this problem to multiplying two polynomials over the integers.

    void _fq_nmod_poly_mul(fq_nmod_struct *rop, const fq_nmod_struct *op1, slong len1, const fq_nmod_struct *op2, slong len2, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, len1 + len2 - 1)`` to the product of ``(op1, len1)``
    # and ``(op2, len2)``, choosing an appropriate algorithm.
    # Permits zero padding.  Does not support aliasing.

    void fq_nmod_poly_mul(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the product of ``op1`` and ``op2``,
    # choosing an appropriate algorithm.

    void _fq_nmod_poly_mullow_classical(fq_nmod_struct *rop, const fq_nmod_struct *op1, slong len1, const fq_nmod_struct *op2, slong len2, slong n, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, n)`` to the first `n` coefficients of
    # ``(op1, len1)`` multiplied by ``(op2, len2)``.
    # Assumes ``0 < n <= len1 + len2 - 1``.  Assumes neither
    # ``len1`` nor ``len2`` is zero.

    void fq_nmod_poly_mullow_classical(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong n, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the product of ``op1`` and ``op2``,
    # computed using the classical or schoolbook method.

    void _fq_nmod_poly_mullow_univariate(fq_nmod_struct *rop, const fq_nmod_struct *op1, slong len1, const fq_nmod_struct *op2, slong len2, slong n, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, n)`` to the lowest `n` coefficients of the product of
    # ``(op1, len1)`` and ``(op2, len2)``, computed using a
    # bivariate to univariate transformation.
    # Assumes that ``len1`` and ``len2`` are positive, but does allow
    # for the polynomials to be zero-padded.  The polynomials may be zero,
    # too.  Assumes `n` is positive.  Supports aliasing between ``rop``,
    # ``op1`` and ``op2``.

    void fq_nmod_poly_mullow_univariate(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong n, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the lowest `n` coefficients of the product of
    # ``poly1`` and ``poly2``, computed using a bivariate to
    # univariate transformation.

    void _fq_nmod_poly_mullow_KS(fq_nmod_struct *rop, const fq_nmod_struct *op1, slong len1, const fq_nmod_struct *op2, slong len2, slong n, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, n)`` to the lowest `n` coefficients of the product of
    # ``(op1, len1)`` and ``(op2, len2)``.
    # Assumes that ``len1`` and ``len2`` are positive, but does allow
    # for the polynomials to be zero-padded.  The polynomials may be zero,
    # too.  Assumes `n` is positive.  Supports aliasing between ``rop``,
    # ``op1`` and ``op2``.

    void fq_nmod_poly_mullow_KS(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong n, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the product of ``op1`` and ``op2``.

    void _fq_nmod_poly_mullow(fq_nmod_struct *rop, const fq_nmod_struct *op1, slong len1, const fq_nmod_struct *op2, slong len2, slong n, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, n)`` to the lowest `n` coefficients of the product of
    # ``(op1, len1)`` and ``(op2, len2)``.
    # Assumes ``0 < n <= len1 + len2 - 1``.  Allows for zero-padding in
    # the inputs.  Does not support aliasing between the inputs and the output.

    void fq_nmod_poly_mullow(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong n, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the lowest `n` coefficients of the product of
    # ``op1`` and ``op2``.

    void _fq_nmod_poly_mulhigh_classical(fq_nmod_struct *res, const fq_nmod_struct *poly1, slong len1, const fq_nmod_struct *poly2, slong len2, slong start, const fq_nmod_ctx_t ctx)
    # Computes the product of ``(poly1, len1)`` and ``(poly2, len2)``
    # and writes the coefficients from ``start`` onwards into the high
    # coefficients of ``res``, the remaining coefficients being arbitrary
    # but reduced.  Assumes that ``len1 >= len2 > 0``. Aliasing of inputs
    # and output is not permitted.  Algorithm is classical multiplication.

    void fq_nmod_poly_mulhigh_classical(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, slong start, const fq_nmod_ctx_t ctx)
    # Computes the product of ``poly1`` and ``poly2`` and writes the
    # coefficients from ``start`` onwards into the high coefficients of
    # ``res``, the remaining coefficients being arbitrary but reduced.
    # Algorithm is classical multiplication.

    void _fq_nmod_poly_mulhigh(fq_nmod_struct *res, const fq_nmod_struct *poly1, slong len1, const fq_nmod_struct *poly2, slong len2, slong start, fq_nmod_ctx_t ctx)
    # Computes the product of ``(poly1, len1)`` and ``(poly2, len2)``
    # and writes the coefficients from ``start`` onwards into the high
    # coefficients of ``res``, the remaining coefficients being arbitrary
    # but reduced.  Assumes that ``len1 >= len2 > 0``. Aliasing of inputs
    # and output is not permitted.

    void fq_nmod_poly_mulhigh(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, slong start, const fq_nmod_ctx_t ctx)
    # Computes the product of ``poly1`` and ``poly2`` and writes the
    # coefficients from ``start`` onwards into the high coefficients of
    # ``res``, the remaining coefficients being arbitrary but reduced.

    void _fq_nmod_poly_mulmod(fq_nmod_struct* res, const fq_nmod_struct* poly1, slong len1, const fq_nmod_struct* poly2, slong len2, const fq_nmod_struct* f, slong lenf, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the remainder of the product of ``poly1``
    # and ``poly2`` upon polynomial division by ``f``.
    # It is required that ``len1 + len2 - lenf > 0``, which is
    # equivalent to requiring that the result will actually be
    # reduced. Otherwise, simply use ``_fq_nmod_poly_mul`` instead.
    # Aliasing of ``f`` and ``res`` is not permitted.

    void fq_nmod_poly_mulmod(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the remainder of the product of ``poly1``
    # and ``poly2`` upon polynomial division by ``f``.

    void _fq_nmod_poly_mulmod_preinv(fq_nmod_struct* res, const fq_nmod_struct* poly1, slong len1, const fq_nmod_struct* poly2, slong len2, const fq_nmod_struct* f, slong lenf, const fq_nmod_struct* finv, slong lenfinv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the remainder of the product of ``poly1``
    # and ``poly2`` upon polynomial division by ``f``.
    # It is required that ``finv`` is the inverse of the reverse of
    # ``f`` mod ``x^lenf``.
    # Aliasing of ``res`` with any of the inputs is not permitted.

    void fq_nmod_poly_mulmod_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_poly_t f, const fq_nmod_poly_t finv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the remainder of the product of ``poly1``
    # and ``poly2`` upon polynomial division by ``f``. ``finv``
    # is the inverse of the reverse of ``f``.

    void _fq_nmod_poly_sqr_classical(fq_nmod_struct *rop, const fq_nmod_struct *op, slong len, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, 2*len - 1)`` to the square of ``(op, len)``,
    # assuming that ``(op,len)`` is not zero and using classical
    # polynomial multiplication.
    # Permits zero padding.  Does not support aliasing of ``rop``
    # with either ``op1`` or ``op2``.

    void fq_nmod_poly_sqr_classical(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the square of ``op`` using classical
    # polynomial multiplication.

    void _fq_nmod_poly_sqr_KS(fq_nmod_struct *rop, const fq_nmod_struct *op, slong len, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, 2*len - 1)`` to the square of ``(op, len)``.
    # Permits zero padding and places no assumptions on the
    # lengths ``len1`` and ``len2``.  Supports aliasing.

    void fq_nmod_poly_sqr_KS(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the square ``op`` using Kronecker substitution,
    # that is, by encoding each coefficient in `\mathbf{F}_{q}` as an integer
    # and reducing this problem to multiplying two polynomials over the integers.

    void _fq_nmod_poly_sqr(fq_nmod_struct *rop, const fq_nmod_struct *op, slong len, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, 2* len - 1)`` to the square of ``(op, len)``,
    # choosing an appropriate algorithm.
    # Permits zero padding.  Does not support aliasing.

    void fq_nmod_poly_sqr(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the square of ``op``,
    # choosing an appropriate algorithm.

    void _fq_nmod_poly_pow(fq_nmod_struct *rop, const fq_nmod_struct *op, slong len, ulong e, const fq_nmod_ctx_t ctx)
    # Sets ``rop = op^e``, assuming that ``e, len > 0`` and that
    # ``rop`` has space for ``e*(len - 1) + 1`` coefficients.  Does
    # not support aliasing.

    void fq_nmod_poly_pow(fq_nmod_poly_t rop, const fq_nmod_poly_t op, ulong e, const fq_nmod_ctx_t ctx)
    # Computes ``rop = op^e``.  If `e` is zero, returns one,
    # so that in particular ``0^0 = 1``.

    void _fq_nmod_poly_powmod_ui_binexp(fq_nmod_struct* res, const fq_nmod_struct* poly, ulong e, const fq_nmod_struct* f, slong lenf, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    # ``f``, using binary exponentiation. We require ``e > 0``.
    # We require ``lenf > 1``. It is assumed that ``poly`` is
    # already reduced modulo ``f`` and zero-padded as necessary to
    # have length exactly ``lenf - 1``. The output ``res`` must
    # have room for ``lenf - 1`` coefficients.

    void fq_nmod_poly_powmod_ui_binexp(fq_nmod_poly_t res, const fq_nmod_poly_t poly, ulong e, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    # ``f``, using binary exponentiation. We require ``e >= 0``.

    void _fq_nmod_poly_powmod_ui_binexp_preinv(fq_nmod_struct* res, const fq_nmod_struct* poly, ulong e, const fq_nmod_struct* f, slong lenf, const fq_nmod_struct* finv, slong lenfinv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    # ``f``, using binary exponentiation. We require ``e > 0``.
    # We require ``finv`` to be the inverse of the reverse of
    # ``f``.
    # We require ``lenf > 1``. It is assumed that ``poly`` is
    # already reduced modulo ``f`` and zero-padded as necessary to
    # have length exactly ``lenf - 1``. The output ``res`` must
    # have room for ``lenf - 1`` coefficients.

    void fq_nmod_poly_powmod_ui_binexp_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly, ulong e, const fq_nmod_poly_t f, const fq_nmod_poly_t finv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    # ``f``, using binary exponentiation. We require ``e >= 0``.
    # We require ``finv`` to be the inverse of the reverse of
    # ``f``.

    void _fq_nmod_poly_powmod_fmpz_binexp(fq_nmod_struct* res, const fq_nmod_struct* poly, const fmpz_t e, const fq_nmod_struct* f, slong lenf, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    # ``f``, using binary exponentiation. We require ``e > 0``.
    # We require ``lenf > 1``. It is assumed that ``poly`` is
    # already reduced modulo ``f`` and zero-padded as necessary to
    # have length exactly ``lenf - 1``. The output ``res`` must
    # have room for ``lenf - 1`` coefficients.

    void fq_nmod_poly_powmod_fmpz_binexp(fq_nmod_poly_t res, const fq_nmod_poly_t poly, const fmpz_t e, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    # ``f``, using binary exponentiation. We require ``e >= 0``.

    void _fq_nmod_poly_powmod_fmpz_binexp_preinv(fq_nmod_struct* res, const fq_nmod_struct* poly, const fmpz_t e, const fq_nmod_struct* f, slong lenf, const fq_nmod_struct* finv, slong lenfinv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    # ``f``, using binary exponentiation. We require ``e > 0``.
    # We require ``finv`` to be the inverse of the reverse of
    # ``f``.
    # We require ``lenf > 1``. It is assumed that ``poly`` is
    # already reduced modulo ``f`` and zero-padded as necessary to
    # have length exactly ``lenf - 1``. The output ``res`` must
    # have room for ``lenf - 1`` coefficients.

    void fq_nmod_poly_powmod_fmpz_binexp_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly, const fmpz_t e, const fq_nmod_poly_t f, const fq_nmod_poly_t finv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    # ``f``, using binary exponentiation. We require ``e >= 0``.
    # We require ``finv`` to be the inverse of the reverse of
    # ``f``.

    void _fq_nmod_poly_powmod_fmpz_sliding_preinv(fq_nmod_struct* res, const fq_nmod_struct* poly, const fmpz_t e, ulong k, const fq_nmod_struct* f, slong lenf, const fq_nmod_struct* finv, slong lenfinv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    # ``f``, using sliding-window exponentiation with window size
    # ``k``. We require ``e > 0``.  We require ``finv`` to be
    # the inverse of the reverse of ``f``. If ``k`` is set to
    # zero, then an "optimum" size will be selected automatically base
    # on ``e``.
    # We require ``lenf > 1``. It is assumed that ``poly`` is
    # already reduced modulo ``f`` and zero-padded as necessary to
    # have length exactly ``lenf - 1``. The output ``res`` must
    # have room for ``lenf - 1`` coefficients.

    void fq_nmod_poly_powmod_fmpz_sliding_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly, const fmpz_t e, ulong k, const fq_nmod_poly_t f, const fq_nmod_poly_t finv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    # ``f``, using sliding-window exponentiation with window size
    # ``k``. We require ``e >= 0``.  We require ``finv`` to be
    # the inverse of the reverse of ``f``.  If ``k`` is set to
    # zero, then an "optimum" size will be selected automatically base
    # on ``e``.

    void _fq_nmod_poly_powmod_x_fmpz_preinv(fq_nmod_struct * res, const fmpz_t e, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * finv, slong lenfinv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to ``x`` raised to the power ``e`` modulo ``f``,
    # using sliding window exponentiation. We require ``e > 0``.
    # We require ``finv`` to be the inverse of the reverse of ``f``.
    # We require ``lenf > 2``. The output ``res`` must have room for
    # ``lenf - 1`` coefficients.

    void fq_nmod_poly_powmod_x_fmpz_preinv(fq_nmod_poly_t res, const fmpz_t e, const fq_nmod_poly_t f, const fq_nmod_poly_t finv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to ``x`` raised to the power ``e``
    # modulo ``f``, using sliding window exponentiation. We require
    # ``e >= 0``. We require ``finv`` to be the inverse of the reverse of
    # ``f``.

    void _fq_nmod_poly_pow_trunc_binexp(fq_nmod_struct * res, const fq_nmod_struct * poly, ulong e, slong trunc, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the low ``trunc`` coefficients of ``poly``
    # (assumed to be zero padded if necessary to length ``trunc``) to
    # the power ``e``. This is equivalent to doing a powering followed
    # by a truncation. We require that ``res`` has enough space for
    # ``trunc`` coefficients, that ``trunc > 0`` and that
    # ``e > 1``. Aliasing is not permitted. Uses the binary
    # exponentiation method.

    void fq_nmod_poly_pow_trunc_binexp(fq_nmod_poly_t res, const fq_nmod_poly_t poly, ulong e, slong trunc, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the low ``trunc`` coefficients of ``poly``
    # to the power ``e``. This is equivalent to doing a powering
    # followed by a truncation. Uses the binary exponentiation method.

    void _fq_nmod_poly_pow_trunc(fq_nmod_struct * res, const fq_nmod_struct * poly, ulong e, slong trunc, const fq_nmod_ctx_t mod)
    # Sets ``res`` to the low ``trunc`` coefficients of ``poly``
    # (assumed to be zero padded if necessary to length ``trunc``) to
    # the power ``e``. This is equivalent to doing a powering followed
    # by a truncation. We require that ``res`` has enough space for
    # ``trunc`` coefficients, that ``trunc > 0`` and that
    # ``e > 1``. Aliasing is not permitted.

    void fq_nmod_poly_pow_trunc(fq_nmod_poly_t res, const fq_nmod_poly_t poly, ulong e, slong trunc, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the low ``trunc`` coefficients of ``poly``
    # to the power ``e``. This is equivalent to doing a powering
    # followed by a truncation.

    void _fq_nmod_poly_shift_left(fq_nmod_struct *rop, const fq_nmod_struct *op, slong len, slong n, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, len + n)`` to ``(op, len)`` shifted left by
    # `n` coefficients.
    # Inserts zero coefficients at the lower end.  Assumes that
    # ``len`` and `n` are positive, and that ``rop`` fits
    # ``len + n`` elements.  Supports aliasing between ``rop`` and
    # ``op``.

    void fq_nmod_poly_shift_left(fq_nmod_poly_t rop, const fq_nmod_poly_t op, slong n, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to ``op`` shifted left by `n` coeffs.  Zero
    # coefficients are inserted.

    void _fq_nmod_poly_shift_right(fq_nmod_struct *rop, const fq_nmod_struct *op, slong len, slong n, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, len - n)`` to ``(op, len)`` shifted right by
    # `n` coefficients.
    # Assumes that ``len`` and `n` are positive, that ``len > n``,
    # and that ``rop`` fits ``len - n`` elements.  Supports
    # aliasing between ``rop`` and ``op``, although in this case
    # the top coefficients of ``op`` are not set to zero.

    void fq_nmod_poly_shift_right(fq_nmod_poly_t rop, const fq_nmod_poly_t op, slong n, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to ``op`` shifted right by `n` coefficients.
    # If `n` is equal to or greater than the current length of
    # ``op``, ``rop`` is set to the zero polynomial.

    slong _fq_nmod_poly_hamming_weight(const fq_nmod_struct *op, slong len, const fq_nmod_ctx_t ctx)
    # Returns the number of non-zero entries in ``(op, len)``.

    slong fq_nmod_poly_hamming_weight(const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    # Returns the number of non-zero entries in the polynomial ``op``.

    void _fq_nmod_poly_divrem(fq_nmod_struct *Q, fq_nmod_struct *R, const fq_nmod_struct *A, slong lenA, const fq_nmod_struct *B, slong lenB, const fq_nmod_t invB, const fq_nmod_ctx_t ctx)
    # Computes ``(Q, lenA - lenB + 1)``, ``(R, lenA)`` such that
    # `A = B Q + R` with `0 \leq \operatorname{len}(R) < \operatorname{len}(B)`.
    # Assumes that the leading coefficient of `B` is invertible
    # and that ``invB`` is its inverse.
    # Assumes that `\operatorname{len}(A), \operatorname{len}(B) > 0`.  Allows zero-padding in
    # ``(A, lenA)``.  `R` and `A` may be aliased, but apart from
    # this no aliasing of input and output operands is allowed.

    void fq_nmod_poly_divrem(fq_nmod_poly_t Q, fq_nmod_poly_t R, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    # Computes `Q`, `R` such that `A = B Q + R` with
    # `0 \leq \operatorname{len}(R) < \operatorname{len}(B)`.
    # Assumes that the leading coefficient of `B` is invertible.  This can
    # be taken for granted the context is for a finite field, that is, when
    # `p` is prime and `f(X)` is irreducible.

    void fq_nmod_poly_divrem_f(fq_nmod_t f, fq_nmod_poly_t Q, fq_nmod_poly_t R, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    # Either finds a non-trivial factor `f` of the modulus of
    # ``ctx``, or computes `Q`, `R` such that `A = B Q + R` and
    # `0 \leq \operatorname{len}(R) < \operatorname{len}(B)`.
    # If the leading coefficient of `B` is invertible, the division with
    # remainder operation is carried out, `Q` and `R` are computed
    # correctly, and `f` is set to `1`.  Otherwise, `f` is set to a
    # non-trivial factor of the modulus and `Q` and `R` are not touched.
    # Assumes that `B` is non-zero.

    void _fq_nmod_poly_rem(fq_nmod_struct *R, const fq_nmod_struct *A, slong lenA, const fq_nmod_struct *B, slong lenB, const fq_nmod_t invB, const fq_nmod_ctx_t ctx)
    # Sets ``R`` to the remainder of the division of ``(A,lenA)`` by
    # ``(B,lenB)``. Assumes that the leading coefficient of ``(B,lenB)``
    # is invertible and that ``invB`` is its inverse.

    void fq_nmod_poly_rem(fq_nmod_poly_t R, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    # Sets ``R`` to the remainder of the division of ``A`` by
    # ``B`` in the context described by ``ctx``.

    void _fq_nmod_poly_div(fq_nmod_struct *Q, const fq_nmod_struct *A, slong lenA, const fq_nmod_struct *B, slong lenB, const fq_nmod_t invB, const fq_nmod_ctx_t ctx)
    # Notationally, computes `Q`, `R` such that `A = B Q + R` with `0
    # \leq \operatorname{len}(R) < \operatorname{len}(B)` but only sets ``(Q, lenA - lenB + 1)``.
    # Allows zero-padding in `A` but not in `B`.  Assumes that the leading coefficient of `B` is a
    # unit.

    void fq_nmod_poly_div(fq_nmod_poly_t Q, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    # Notionally finds polynomials `Q` and `R` such that `A = B Q + R` with
    # `\operatorname{len}(R) < \operatorname{len}(B)`, but returns only ``Q``. If `\operatorname{len}(B) = 0` an
    # exception is raised.

    void _fq_nmod_poly_div_newton_n_preinv(fq_nmod_struct* Q, const fq_nmod_struct* A, slong lenA, const fq_nmod_struct* B, slong lenB, const fq_nmod_struct* Binv, slong lenBinv, const fq_nmod_ctx_t ctx)
    # Notionally computes polynomials `Q` and `R` such that `A = BQ + R` with
    # `\operatorname{len}(R)` less than ``lenB``, where ``A`` is of length ``lenA``
    # and ``B`` is of length ``lenB``, but return only `Q`.
    # We require that `Q` have space for ``lenA - lenB + 1`` coefficients
    # and assume that the leading coefficient of `B` is a unit. Furthermore, we
    # assume that `Binv` is the inverse of the reverse of `B` mod `x^{\operatorname{len}(B)}`.
    # The algorithm used is to reverse the polynomials and divide the
    # resulting power series, then reverse the result.

    void fq_nmod_poly_div_newton_n_preinv(fq_nmod_poly_t Q, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_poly_t Binv, const fq_nmod_ctx_t ctx)
    # Notionally computes `Q` and `R` such that `A = BQ + R` with
    # `\operatorname{len}(R) < \operatorname{len}(B)`, but returns only `Q`.
    # We assume that the leading coefficient of `B` is a unit and that `Binv` is
    # the inverse of the reverse of `B` mod `x^{\operatorname{len}(B)}`.
    # It is required that the length of `A` is less than or equal to
    # 2*the length of `B` - 2.
    # The algorithm used is to reverse the polynomials and divide the
    # resulting power series, then reverse the result.

    void _fq_nmod_poly_divrem_newton_n_preinv(fq_nmod_struct* Q, fq_nmod_struct* R, const fq_nmod_struct* A, slong lenA, const fq_nmod_struct* B, slong lenB, const fq_nmod_struct* Binv, slong lenBinv, const fq_nmod_ctx_t ctx)
    # Computes `Q` and `R` such that `A = BQ + R` with `\operatorname{len}(R)` less
    # than ``lenB``, where `A` is of length ``lenA`` and `B` is of
    # length ``lenB``. We require that `Q` have space for
    # ``lenA - lenB + 1`` coefficients. Furthermore, we assume that `Binv` is
    # the inverse of the reverse of `B` mod `x^{\operatorname{len}(B)}`. The algorithm
    # used is to call :func:`div_newton_preinv` and then multiply out
    # and compute the remainder.

    void fq_nmod_poly_divrem_newton_n_preinv(fq_nmod_poly_t Q, fq_nmod_poly_t R, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_poly_t Binv, const fq_nmod_ctx_t ctx)
    # Computes `Q` and `R` such that `A = BQ + R` with `\operatorname{len}(R) <
    # \operatorname{len}(B)`.  We assume `Binv` is the inverse of the reverse of `B`
    # mod `x^{\operatorname{len}(B)}`.
    # It is required that the length of `A` is less than or equal to
    # 2*the length of `B` - 2.
    # The algorithm used is to call :func:`div_newton` and then
    # multiply out and compute the remainder.

    void _fq_nmod_poly_inv_series_newton(fq_nmod_struct* Qinv, const fq_nmod_struct* Q, slong n, const fq_nmod_t cinv, const fq_nmod_ctx_t ctx)
    # Given ``Q`` of length ``n`` whose constant coefficient is
    # invertible modulo the given modulus, find a polynomial ``Qinv``
    # of length ``n`` such that ``Q * Qinv`` is ``1`` modulo
    # `x^n`. Requires ``n > 0``.  This function can be viewed as
    # inverting a power series via Newton iteration.

    void fq_nmod_poly_inv_series_newton(fq_nmod_poly_t Qinv, const fq_nmod_poly_t Q, slong n, const fq_nmod_ctx_t ctx)
    # Given ``Q`` find ``Qinv`` such that ``Q * Qinv`` is
    # ``1`` modulo `x^n`. The constant coefficient of ``Q`` must
    # be invertible modulo the modulus of ``Q``. An exception is
    # raised if this is not the case or if ``n = 0``. This function
    # can be viewed as inverting a power series via Newton iteration.

    void _fq_nmod_poly_inv_series(fq_nmod_struct* Qinv, const fq_nmod_struct* Q, slong n, const fq_nmod_t cinv, const fq_nmod_ctx_t ctx)
    # Given ``Q`` of length ``n`` whose constant coefficient is
    # invertible modulo the given modulus, find a polynomial ``Qinv``
    # of length ``n`` such that ``Q * Qinv`` is ``1`` modulo
    # `x^n`. Requires ``n > 0``.

    void fq_nmod_poly_inv_series(fq_nmod_poly_t Qinv, const fq_nmod_poly_t Q, slong n, const fq_nmod_ctx_t ctx)
    # Given ``Q`` find ``Qinv`` such that ``Q * Qinv`` is
    # ``1`` modulo `x^n`. The constant coefficient of ``Q`` must
    # be invertible modulo the modulus of ``Q``. An exception is
    # raised if this is not the case or if ``n = 0``.

    void _fq_nmod_poly_div_series(fq_nmod_struct *Q, const fq_nmod_struct *A, mp_limb_signed_t Alen, const fq_nmod_struct *B, mp_limb_signed_t Blen, mp_limb_signed_t n, const fq_nmod_ctx_t ctx)
    # Set ``(Q, n)`` to the quotient of the series ``(A, Alen``) and
    # ``(B, Blen)`` assuming ``Alen, Blen <= n``. We assume the bottom
    # coefficient of ``B`` is invertible.

    void fq_nmod_poly_div_series(fq_nmod_poly_t Q, const fq_nmod_poly_t A, const fq_nmod_poly_t B, slong n, fq_nmod_ctx_t ctx)
    # Set `Q` to the quotient of the series `A` by `B`, thinking of the series as
    # though they were of length `n`. We assume that the bottom coefficient of
    # `B` is invertible.

    void fq_nmod_poly_gcd(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the greatest common divisor of ``op1`` and
    # ``op2``, using the either the Euclidean or HGCD algorithm. The
    # GCD of zero polynomials is defined to be zero, whereas the GCD of
    # the zero polynomial and some other polynomial `P` is defined to be
    # `P`. Except in the case where the GCD is zero, the GCD `G` is made
    # monic.

    slong _fq_nmod_poly_gcd(fq_nmod_struct* G, const fq_nmod_struct* A, slong lenA, const fq_nmod_struct* B, slong lenB, const fq_nmod_ctx_t ctx)
    # Computes the GCD of `A` of length ``lenA`` and `B` of length
    # ``lenB``, where ``lenA >= lenB > 0`` and sets `G` to it. The
    # length of the GCD `G` is returned by the function. No attempt is
    # made to make the GCD monic. It is required that `G` have space for
    # ``lenB`` coefficients.

    slong _fq_nmod_poly_gcd_euclidean_f(fq_nmod_t f, fq_nmod_struct *G, const fq_nmod_struct *A, slong lenA, const fq_nmod_struct *B, slong lenB, const fq_nmod_ctx_t ctx)
    # Either sets `f = 1` and `G` to the greatest common divisor of
    # `(A,\operatorname{len}(A))` and `(B, \operatorname{len}(B))` and returns its length, or sets
    # `f` to a non-trivial factor of the modulus of ``ctx`` and leaves
    # the contents of the vector `(G, lenB)` undefined.
    # Assumes that `\operatorname{len}(A) \geq \operatorname{len}(B) > 0` and that the vector `G`
    # has space for sufficiently many coefficients.

    void fq_nmod_poly_gcd_euclidean_f(fq_nmod_t f, fq_nmod_poly_t G, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    # Either sets `f = 1` and `G` to the greatest common divisor of `A`
    # and `B` or sets `f` to a factor of the modulus of ``ctx``.

    slong _fq_nmod_poly_xgcd(fq_nmod_struct *G, fq_nmod_struct *S, fq_nmod_struct *T, const fq_nmod_struct *A, slong lenA, const fq_nmod_struct *B, slong lenB, const fq_nmod_ctx_t ctx)
    # Computes the GCD of `A` and `B` together with cofactors `S` and `T`
    # such that `S A + T B = G`.  Returns the length of `G`.
    # Assumes that `\operatorname{len}(A) \geq \operatorname{len}(B) \geq 1` and
    # `(\operatorname{len}(A),\operatorname{len}(B)) \neq (1,1)`.
    # No attempt is made to make the GCD monic.
    # Requires that `G` have space for `\operatorname{len}(B)` coefficients.  Writes
    # `\operatorname{len}(B)-1` and `\operatorname{len}(A)-1` coefficients to `S` and `T`, respectively.
    # Note that, in fact, `\operatorname{len}(S) \leq \max(\operatorname{len}(B) - \operatorname{len}(G), 1)` and
    # `\operatorname{len}(T) \leq \max(\operatorname{len}(A) - \operatorname{len}(G), 1)`.
    # No aliasing of input and output operands is permitted.

    void fq_nmod_poly_xgcd(fq_nmod_poly_t G, fq_nmod_poly_t S, fq_nmod_poly_t T, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    # Computes the GCD of `A` and `B`. The GCD of zero polynomials is
    # defined to be zero, whereas the GCD of the zero polynomial and some other
    # polynomial `P` is defined to be `P`. Except in the case where
    # the GCD is zero, the GCD `G` is made monic.
    # Polynomials ``S`` and ``T`` are computed such that
    # ``S*A + T*B = G``. The length of ``S`` will be at most
    # ``lenB`` and the length of ``T`` will be at most ``lenA``.

    slong _fq_nmod_poly_xgcd_euclidean_f(fq_nmod_t f, fq_nmod_struct *G, fq_nmod_struct *S, fq_nmod_struct *T, const fq_nmod_struct *A, slong lenA, const fq_nmod_struct *B, slong lenB, const fq_nmod_ctx_t ctx)
    # Either sets `f = 1` and computes the GCD of `A` and `B` together
    # with cofactors `S` and `T` such that `S A + T B = G`; otherwise,
    # sets `f` to a non-trivial factor of the modulus of ``ctx`` and
    # leaves `G`, `S`, and `T` undefined.  Returns the length of `G`.
    # Assumes that `\operatorname{len}(A) \geq \operatorname{len}(B) \geq 1` and
    # `(\operatorname{len}(A),\operatorname{len}(B)) \neq (1,1)`.
    # No attempt is made to make the GCD monic.
    # Requires that `G` have space for `\operatorname{len}(B)` coefficients.  Writes
    # `\operatorname{len}(B)-1` and `\operatorname{len}(A)-1` coefficients to `S` and `T`, respectively.
    # Note that, in fact, `\operatorname{len}(S) \leq \max(\operatorname{len}(B) - \operatorname{len}(G), 1)` and
    # `\operatorname{len}(T) \leq \max(\operatorname{len}(A) - \operatorname{len}(G), 1)`.
    # No aliasing of input and output operands is permitted.

    void fq_nmod_poly_xgcd_euclidean_f(fq_nmod_t f, fq_nmod_poly_t G, fq_nmod_poly_t S, fq_nmod_poly_t T, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    # Either sets `f = 1` and computes the GCD of `A` and `B` or sets
    # `f` to a non-trivial factor of the modulus of ``ctx``.
    # If the GCD is computed, polynomials ``S`` and ``T`` are
    # computed such that ``S*A + T*B = G``; otherwise, they are
    # undefined.  The length of ``S`` will be at most ``lenB`` and
    # the length of ``T`` will be at most ``lenA``.
    # The GCD of zero polynomials is defined to be zero, whereas the GCD
    # of the zero polynomial and some other polynomial `P` is defined to
    # be `P`. Except in the case where the GCD is zero, the GCD `G` is
    # made monic.

    int _fq_nmod_poly_divides(fq_nmod_struct *Q, const fq_nmod_struct *A, slong lenA, const fq_nmod_struct *B, slong lenB, const fq_nmod_t invB, const fq_nmod_ctx_t ctx)
    # Returns `1` if ``(B, lenB)`` divides ``(A, lenA)`` exactly and
    # sets `Q` to the quotient, otherwise returns `0`.
    # It is assumed that `\operatorname{len}(A) \geq \operatorname{len}(B) > 0` and that `Q` has space
    # for `\operatorname{len}(A) - \operatorname{len}(B) + 1` coefficients.
    # Aliasing of `Q` with either of the inputs is not permitted.
    # This function is currently unoptimised and provided for convenience
    # only.

    int fq_nmod_poly_divides(fq_nmod_poly_t Q, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    # Returns `1` if `B` divides `A` exactly and sets `Q` to the quotient,
    # otherwise returns `0`.
    # This function is currently unoptimised and provided for convenience
    # only.

    void _fq_nmod_poly_derivative(fq_nmod_struct *rop, const fq_nmod_struct *op, slong len, const fq_nmod_ctx_t ctx)
    # Sets ``(rop, len - 1)`` to the derivative of ``(op, len)``.
    # Also handles the cases where ``len`` is `0` or `1` correctly.
    # Supports aliasing of ``rop`` and ``op``.

    void fq_nmod_poly_derivative(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the derivative of ``op``.

    void _fq_nmod_poly_invsqrt_series(fq_nmod_struct * g, const fq_nmod_struct * h, slong n, fq_nmod_ctx_t mod)
    # Set the first `n` terms of `g` to the series expansion of `1/\sqrt{h}`.
    # It is assumed that `n > 0`, that `h` has constant term 1 and that `h`
    # is zero-padded as necessary to length `n`. Aliasing is not permitted.

    void fq_nmod_poly_invsqrt_series(fq_nmod_poly_t g, const fq_nmod_poly_t h, slong n, fq_nmod_ctx_t ctx)
    # Set `g` to the series expansion of `1/\sqrt{h}` to order `O(x^n)`.
    # It is assumed that `h` has constant term 1.

    void _fq_nmod_poly_sqrt_series(fq_nmod_struct * g, const fq_nmod_struct * h, slong n, fq_nmod_ctx_t ctx)
    # Set the first `n` terms of `g` to the series expansion of `\sqrt{h}`.
    # It is assumed that `n > 0`, that `h` has constant term 1 and that `h`
    # is zero-padded as necessary to length `n`. Aliasing is not permitted.

    void fq_nmod_poly_sqrt_series(fq_nmod_poly_t g, const fq_nmod_poly_t h, slong n, fq_nmod_ctx_t ctx)
    # Set `g` to the series expansion of `\sqrt{h}` to order `O(x^n)`.
    # It is assumed that `h` has constant term 1.

    int _fq_nmod_poly_sqrt(fq_nmod_struct * s, const fq_nmod_struct * p, slong n, fq_nmod_ctx_t mod)
    # If ``(p, n)`` is a perfect square, sets ``(s, n / 2 + 1)``
    # to a square root of `p` and returns 1. Otherwise returns 0.

    int fq_nmod_poly_sqrt(fq_nmod_poly_t s, const fq_nmod_poly_t p, fq_nmod_ctx_t mod)
    # If `p` is a perfect square, sets `s` to a square root of `p`
    # and returns 1. Otherwise returns 0.

    void _fq_nmod_poly_evaluate_fq_nmod(fq_nmod_t rop, const fq_nmod_struct *op, slong len, const fq_nmod_t a, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to ``(op, len)`` evaluated at `a`.
    # Supports zero padding.  There are no restrictions on ``len``, that
    # is, ``len`` is allowed to be zero, too.

    void fq_nmod_poly_evaluate_fq_nmod(fq_nmod_t rop, const fq_nmod_poly_t f, const fq_nmod_t a, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the value of `f(a)`.
    # As the coefficient ring `\mathbf{F}_q` is finite, Horner's method
    # is sufficient.

    void _fq_nmod_poly_compose(fq_nmod_struct *rop, const fq_nmod_struct *op1, slong len1, const fq_nmod_struct *op2, slong len2, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the composition of ``(op1, len1)`` and
    # ``(op2, len2)``.
    # Assumes that ``rop`` has space for ``(len1-1)*(len2-1) + 1``
    # coefficients.  Assumes that ``op1`` and ``op2`` are non-zero
    # polynomials.  Does not support aliasing between any of the inputs and
    # the output.

    void fq_nmod_poly_compose(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    # Sets ``rop`` to the composition of ``op1`` and ``op2``.
    # To be precise about the order of composition, denoting ``rop``,
    # ``op1``, and ``op2`` by `f`, `g`, and `h`, respectively,
    # sets `f(t) = g(h(t))`.

    void _fq_nmod_poly_compose_mod_horner(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the composition `f(g)` modulo `h`. We require that
    # `h` is nonzero and that the length of `g` is one less than the
    # length of `h` (possibly with zero padding). The output is not allowed
    # to be aliased with any of the inputs.
    # The algorithm used is Horner's rule.

    void fq_nmod_poly_compose_mod_horner(fq_nmod_poly_t res, const fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_poly_t h, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the composition `f(g)` modulo `h`. We require that
    # `h` is nonzero. The algorithm used is Horner's rule.

    void _fq_nmod_poly_compose_mod_horner_preinv(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_struct * hinv, slong lenhiv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the composition `f(g)` modulo `h`. We require
    # that `h` is nonzero and that the length of `g` is one less than
    # the length of `h` (possibly with zero padding). We also require
    # that the length of `f` is less than the length of
    # `h`. Furthermore, we require ``hinv`` to be the inverse of the
    # reverse of ``h``.  The output is not allowed to be aliased with
    # any of the inputs.
    # The algorithm used is Horner's rule.

    void fq_nmod_poly_compose_mod_horner_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_poly_t h, const fq_nmod_poly_t hinv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the composition `f(g)` modulo `h`. We require
    # that `h` is nonzero and that `f` has smaller degree than
    # `h`. Furthermore, we require ``hinv`` to be the inverse of the
    # reverse of ``h``.  The algorithm used is Horner's rule.

    void _fq_nmod_poly_compose_mod_brent_kung(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the composition `f(g)` modulo `h`. We require
    # that `h` is nonzero and that the length of `g` is one less than
    # the length of `h` (possibly with zero padding). We also require
    # that the length of `f` is less than the length of `h`. The output
    # is not allowed to be aliased with any of the inputs.
    # The algorithm used is the Brent-Kung matrix algorithm.

    void fq_nmod_poly_compose_mod_brent_kung(fq_nmod_poly_t res, const fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_poly_t h, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the composition `f(g)` modulo `h`. We require
    # that `h` is nonzero and that `f` has smaller degree than `h`.  The
    # algorithm used is the Brent-Kung matrix algorithm.

    void _fq_nmod_poly_compose_mod_brent_kung_preinv(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_struct * hinv, slong lenhiv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the composition `f(g)` modulo `h`. We require
    # that `h` is nonzero and that the length of `g` is one less than
    # the length of `h` (possibly with zero padding). We also require
    # that the length of `f` is less than the length of
    # `h`. Furthermore, we require ``hinv`` to be the inverse of the
    # reverse of ``h``.  The output is not allowed to be aliased with
    # any of the inputs.
    # The algorithm used is the Brent-Kung matrix algorithm.

    void fq_nmod_poly_compose_mod_brent_kung_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_poly_t h, const fq_nmod_poly_t hinv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the composition `f(g)` modulo `h`. We require
    # that `h` is nonzero and that `f` has smaller degree than
    # `h`. Furthermore, we require ``hinv`` to be the inverse of the
    # reverse of ``h``.  The algorithm used is the Brent-Kung matrix
    # algorithm.

    void _fq_nmod_poly_compose_mod(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the composition `f(g)` modulo `h`. We require
    # that `h` is nonzero and that the length of `g` is one less than
    # the length of `h` (possibly with zero padding). The output is not
    # allowed to be aliased with any of the inputs.

    void fq_nmod_poly_compose_mod(fq_nmod_poly_t res, const fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_poly_t h, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the composition `f(g)` modulo `h`. We require
    # that `h` is nonzero.

    void _fq_nmod_poly_compose_mod_preinv(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_struct * hinv, slong lenhiv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the composition `f(g)` modulo `h`. We require
    # that `h` is nonzero and that the length of `g` is one less than
    # the length of `h` (possibly with zero padding). We also require
    # that the length of `f` is less than the length of
    # `h`. Furthermore, we require ``hinv`` to be the inverse of the
    # reverse of ``h``.  The output is not allowed to be aliased with
    # any of the inputs.

    void fq_nmod_poly_compose_mod_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_poly_t h, const fq_nmod_poly_t hinv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the composition `f(g)` modulo `h`. We require
    # that `h` is nonzero and that `f` has smaller degree than
    # `h`. Furthermore, we require ``hinv`` to be the inverse of the
    # reverse of ``h``.

    void _fq_nmod_poly_reduce_matrix_mod_poly (fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)
    # Sets the ith row of ``A`` to the reduction of the ith row of `B` modulo
    # `f` for `i=1,\ldots,\sqrt{\deg(f)}`. We require `B` to be at least
    # a `\sqrt{\deg(f)}\times \deg(f)` matrix and `f` to be nonzero.

    void _fq_nmod_poly_precompute_matrix (fq_nmod_mat_t A, const fq_nmod_struct* f, const fq_nmod_struct* g, slong leng, const fq_nmod_struct* ginv, slong lenginv, const fq_nmod_ctx_t ctx)
    # Sets the ith row of ``A`` to `f^i` modulo `g` for
    # `i=1,\ldots,\sqrt{\deg(g)}`. We require `A` to be a
    # `\sqrt{\deg(g)}\times \deg(g)` matrix. We require ``ginv`` to
    # be the inverse of the reverse of ``g`` and `g` to be nonzero.

    void fq_nmod_poly_precompute_matrix (fq_nmod_mat_t A, const fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_poly_t ginv, const fq_nmod_ctx_t ctx)
    # Sets the ith row of ``A`` to `f^i` modulo `g` for
    # `i=1,\ldots,\sqrt{\deg(g)}`. We require `A` to be a
    # `\sqrt{\deg(g)}\times \deg(g)` matrix. We require ``ginv`` to
    # be the inverse of the reverse of ``g``.

    void _fq_nmod_poly_compose_mod_brent_kung_precomp_preinv(fq_nmod_struct* res, const fq_nmod_struct* f, slong lenf, const fq_nmod_mat_t A, const fq_nmod_struct* h, slong lenh, const fq_nmod_struct* hinv, slong lenhinv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the composition `f(g)` modulo `h`. We require
    # that `h` is nonzero. We require that the ith row of `A` contains
    # `g^i` for `i=1,\ldots,\sqrt{\deg(h)}`, i.e. `A` is a
    # `\sqrt{\deg(h)}\times \deg(h)` matrix. We also require that the
    # length of `f` is less than the length of `h`. Furthermore, we
    # require ``hinv`` to be the inverse of the reverse of ``h``.
    # The output is not allowed to be aliased with any of the inputs.
    # The algorithm used is the Brent-Kung matrix algorithm.

    void fq_nmod_poly_compose_mod_brent_kung_precomp_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t f, const fq_nmod_mat_t A, const fq_nmod_poly_t h, const fq_nmod_poly_t hinv, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to the composition `f(g)` modulo `h`. We require
    # that the ith row of `A` contains `g^i` for
    # `i=1,\ldots,\sqrt{\deg(h)}`, i.e. `A` is a `\sqrt{\deg(h)}\times
    # \deg(h)` matrix. We require that `h` is nonzero and that `f` has
    # smaller degree than `h`. Furthermore, we require ``hinv`` to be
    # the inverse of the reverse of ``h``. This version of Brent-Kung
    # modular composition is particularly useful if one has to perform
    # several modular composition of the form `f(g)` modulo `h` for
    # fixed `g` and `h`.

    int _fq_nmod_poly_fprint_pretty(FILE *file, const fq_nmod_struct *poly, slong len, const char *x, const fq_nmod_ctx_t ctx)
    # Prints the pretty representation of ``(poly, len)`` to the stream
    # ``file``, using the string ``x`` to represent the indeterminate.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int fq_nmod_poly_fprint_pretty(FILE * file, const fq_nmod_poly_t poly, const char *x, const fq_nmod_ctx_t ctx)
    # Prints the pretty representation of ``poly`` to the stream
    # ``file``, using the string ``x`` to represent the indeterminate.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int _fq_nmod_poly_print_pretty(const fq_nmod_struct *poly, slong len, const char *x, const fq_nmod_ctx_t ctx)
    # Prints the pretty representation of ``(poly, len)`` to ``stdout``,
    # using the string ``x`` to represent the indeterminate.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int fq_nmod_poly_print_pretty(const fq_nmod_poly_t poly, const char *x, const fq_nmod_ctx_t ctx)
    # Prints the pretty representation of ``poly`` to ``stdout``,
    # using the string ``x`` to represent the indeterminate.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int _fq_nmod_poly_fprint(FILE *file, const fq_nmod_struct *poly, slong len, const fq_nmod_ctx_t ctx)
    # Prints the pretty representation of ``(poly, len)`` to the stream
    # ``file``.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int fq_nmod_poly_fprint(FILE * file, const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    # Prints the pretty representation of ``poly`` to the stream
    # ``file``.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int _fq_nmod_poly_print(const fq_nmod_struct *poly, slong len, const fq_nmod_ctx_t ctx)
    # Prints the pretty representation of ``(poly, len)`` to ``stdout``.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int fq_nmod_poly_print(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    # Prints the representation of ``poly`` to ``stdout``.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    char * _fq_nmod_poly_get_str(const fq_nmod_struct * poly, slong len, const fq_nmod_ctx_t ctx)
    # Returns the plain FLINT string representation of the polynomial
    # ``(poly, len)``.

    char * fq_nmod_poly_get_str(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    # Returns the plain FLINT string representation of the polynomial
    # ``poly``.

    char * _fq_nmod_poly_get_str_pretty(const fq_nmod_struct * poly, slong len, const char * x, const fq_nmod_ctx_t ctx)
    # Returns a pretty representation of the polynomial
    # ``(poly, len)`` using the null-terminated string ``x`` as the
    # variable name.

    char * fq_nmod_poly_get_str_pretty(const fq_nmod_poly_t poly, const char * x, const fq_nmod_ctx_t ctx)
    # Returns a pretty representation of the polynomial ``poly`` using the
    # null-terminated string ``x`` as the variable name

    void fq_nmod_poly_inflate(fq_nmod_poly_t result, const fq_nmod_poly_t input, ulong inflation, const fq_nmod_ctx_t ctx)
    # Sets ``result`` to the inflated polynomial `p(x^n)` where
    # `p` is given by ``input`` and `n` is given by ``inflation``.

    void fq_nmod_poly_deflate(fq_nmod_poly_t result, const fq_nmod_poly_t input, ulong deflation, const fq_nmod_ctx_t ctx)
    # Sets ``result`` to the deflated polynomial `p(x^{1/n})` where
    # `p` is given by ``input`` and `n` is given by ``deflation``.
    # Requires `n > 0`.

    ulong fq_nmod_poly_deflation(const fq_nmod_poly_t input, const fq_nmod_ctx_t ctx)
    # Returns the largest integer by which ``input`` can be deflated.
    # As special cases, returns 0 if ``input`` is the zero polynomial
    # and 1 of ``input`` is a constant polynomial.
