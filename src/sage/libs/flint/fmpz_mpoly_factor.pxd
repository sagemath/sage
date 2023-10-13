# distutils: libraries = flint
# distutils: depends = flint/fmpz_mpoly_factor.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fmpz_mpoly_factor_init(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
    # Initialise *f*.

    void fmpz_mpoly_factor_clear(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
    # Clear *f*.

    void fmpz_mpoly_factor_swap(fmpz_mpoly_factor_t f, fmpz_mpoly_factor_t g, const fmpz_mpoly_ctx_t ctx)
    # Efficiently swap *f* and *g*.

    long fmpz_mpoly_factor_length(const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
    # Return the length of the product in *f*.

    void fmpz_mpoly_factor_get_constant_fmpz(fmpz_t c, const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_factor_get_constant_fmpq(fmpq_t c, const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
    # Set `c` to the constant of *f*.

    void fmpz_mpoly_factor_get_base(fmpz_mpoly_t B, const fmpz_mpoly_factor_t f, long i, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_factor_swap_base(fmpz_mpoly_t B, fmpz_mpoly_factor_t f, long i, const fmpz_mpoly_ctx_t ctx)
    # Set (resp. swap) *B* to (resp. with) the base of the term of index `i` in  *A*.

    long fmpz_mpoly_factor_get_exp_si(fmpz_mpoly_factor_t f, long i, const fmpz_mpoly_ctx_t ctx)
    # Return the exponent of the term of index `i` in *A*. It is assumed to fit an ``slong``.

    void fmpz_mpoly_factor_sort(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
    # Sort the product of *f* first by exponent and then by base.

    int fmpz_mpoly_factor_squarefree(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    # Set *f* to a factorization of *A* where the bases are primitive and
    # pairwise relatively prime. If the product of all irreducible factors with
    # a given exponent is desired, it is recommended to call :func:`fmpz_mpoly_factor_sort`
    # and then multiply the bases with the desired exponent.

    int fmpz_mpoly_factor(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    # Set *f* to a factorization of *A* where the bases are irreducible.
