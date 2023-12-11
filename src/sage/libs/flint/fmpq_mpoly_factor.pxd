# distutils: libraries = flint
# distutils: depends = flint/fmpq_mpoly_factor.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fmpq_mpoly_factor_init(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx) noexcept
    # Initialise *f*.

    void fmpq_mpoly_factor_clear(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx) noexcept
    # Clear *f*.

    void fmpq_mpoly_factor_swap(fmpq_mpoly_factor_t f, fmpq_mpoly_factor_t g, const fmpq_mpoly_ctx_t ctx) noexcept
    # Efficiently swap *f* and *g*.

    slong fmpq_mpoly_factor_length(const fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx) noexcept
    # Return the length of the product in *f*.

    void fmpq_mpoly_factor_get_constant_fmpq(fmpq_t c, const fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx) noexcept
    # Set *c* to the constant of *f*.

    void fmpq_mpoly_factor_get_base(fmpq_mpoly_t B, const fmpq_mpoly_factor_t f, slong i, const fmpq_mpoly_ctx_t ctx) noexcept
    void fmpq_mpoly_factor_swap_base(fmpq_mpoly_t B, fmpq_mpoly_factor_t f, slong i, const fmpq_mpoly_ctx_t ctx) noexcept
    # Set (resp. swap) *B* to (resp. with) the base of the term of index *i* in  *A*.

    slong fmpq_mpoly_factor_get_exp_si(fmpq_mpoly_factor_t f, slong i, const fmpq_mpoly_ctx_t ctx) noexcept
    # Return the exponent of the term of index *i* in *A*. It is assumed to fit an ``slong``.

    void fmpq_mpoly_factor_sort(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx) noexcept
    # Sort the product of *f* first by exponent and then by base.

    int fmpq_mpoly_factor_make_monic(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx) noexcept
    int fmpq_mpoly_factor_make_integral(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx) noexcept
    # Make the bases in *f* monic (resp. integral and primitive with positive leading coefficient).
    # Return `1` for success, `0` for failure.

    int fmpq_mpoly_factor_squarefree(fmpq_mpoly_factor_t f, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx) noexcept
    # Set *f* to a factorization of *A* where the bases are primitive and
    # pairwise relatively prime. If the product of all irreducible factors with
    # a given exponent is desired, it is recommended to call :func:`fmpq_mpoly_factor_sort`
    # and then multiply the bases with the desired exponent.

    int fmpq_mpoly_factor(fmpq_mpoly_factor_t f, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx) noexcept
    # Set *f* to a factorization of *A* where the bases are irreducible.
