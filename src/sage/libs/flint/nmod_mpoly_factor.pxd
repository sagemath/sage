# distutils: libraries = flint
# distutils: depends = flint/nmod_mpoly_factor.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void nmod_mpoly_factor_init(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
    # Initialise *f*.

    void nmod_mpoly_factor_clear(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
    # Clear *f*.

    void nmod_mpoly_factor_swap(nmod_mpoly_factor_t f, nmod_mpoly_factor_t g, const nmod_mpoly_ctx_t ctx)
    # Efficiently swap *f* and *g*.

    long nmod_mpoly_factor_length(const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
    # Return the length of the product in *f*.

    unsigned long nmod_mpoly_factor_get_constant_ui(const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
    # Return the constant of *f*.

    void nmod_mpoly_factor_get_base(nmod_mpoly_t p, const nmod_mpoly_factor_t f, long i, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_factor_swap_base(nmod_mpoly_t p, nmod_mpoly_factor_t f, long i, const nmod_mpoly_ctx_t ctx)
    # Set (resp. swap) *B* to (resp. with) the base of the term of index `i` in  *A*.

    long nmod_mpoly_factor_get_exp_si(nmod_mpoly_factor_t f, long i, const nmod_mpoly_ctx_t ctx)
    # Return the exponent of the term of index `i` in *A*. It is assumed to fit an ``slong``.

    void nmod_mpoly_factor_sort(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
    # Sort the product of *f* first by exponent and then by base.

    int nmod_mpoly_factor_squarefree(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Set *f* to a factorization of *A* where the bases are primitive and
    # pairwise relatively prime. If the product of all irreducible factors with
    # a given exponent is desired, it is recommended to call :func:`nmod_mpoly_factor_sort`
    # and then multiply the bases with the desired exponent.

    int nmod_mpoly_factor(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Set *f* to a factorization of *A* where the bases are irreducible.
