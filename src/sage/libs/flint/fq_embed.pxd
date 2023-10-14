# distutils: libraries = flint
# distutils: depends = flint/fq_embed.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fq_embed_gens(fq_t gen_sub, fq_t gen_sup, fmpz_mod_poly_t minpoly, const fq_ctx_t sub_ctx, const fq_ctx_t sup_ctx)
    # Given two contexts ``sub_ctx`` and ``sup_ctx``, such that
    # ``degree(sub_ctx)`` divides ``degree(sup_ctx)``, compute:
    # * an element ``gen_sub`` in ``sub_ctx`` such that
    # ``gen_sub`` generates the finite field defined by
    # ``sub_ctx``,
    # * its minimal polynomial ``minpoly``,
    # * a root ``gen_sup`` of ``minpoly`` inside the field
    # defined by ``sup_ctx``.
    # These data uniquely define an embedding of ``sub_ctx`` into
    # ``sup_ctx``.

    void _fq_embed_gens_naive(fq_t gen_sub, fq_t gen_sup, fmpz_mod_poly_t minpoly, const fq_ctx_t sub_ctx, const fq_ctx_t sup_ctx)
    # Given two contexts ``sub_ctx`` and ``sup_ctx``, such that
    # ``degree(sub_ctx)`` divides ``degree(sup_ctx)``, compute an
    # embedding of ``sub_ctx`` into ``sup_ctx`` defined as follows:
    # * ``gen_sub`` is the canonical generator of ``sup_ctx``
    # (i.e., the class of `X`),
    # * ``minpoly`` is the defining polynomial of ``sub_ctx``,
    # * ``gen_sup`` is a root of ``minpoly`` inside the field
    # defined by ``sup_ctx``.

    void fq_embed_matrices(fmpz_mod_mat_t embed, fmpz_mod_mat_t project, const fq_t gen_sub, const fq_ctx_t sub_ctx, const fq_t gen_sup, const fq_ctx_t sup_ctx, const fmpz_mod_poly_t gen_minpoly)
    # Given:
    # * two contexts ``sub_ctx`` and ``sup_ctx``, of
    # respective degrees `m` and `n`, such that `m` divides `n`;
    # * a generator ``gen_sub`` of ``sub_ctx``, its minimal
    # polynomial ``gen_minpoly``, and a root ``gen_sup`` of
    # ``gen_minpoly`` in ``sup_ctx``, as returned by
    # :func:`fq_embed_gens`;
    # Compute:
    # * the `n\times m` matrix ``embed`` mapping ``gen_sub``
    # to ``gen_sup``, and all their powers accordingly;
    # * an `m\times n` matrix ``project`` such that
    # ``project`` `\times` ``embed`` is the `m\times m` identity
    # matrix.

    void fq_embed_trace_matrix(fmpz_mod_mat_t res, const fmpz_mod_mat_t basis, const fq_ctx_t sub_ctx, const fq_ctx_t sup_ctx)
    # Given:
    # * two contexts ``sub_ctx`` and ``sup_ctx``, of degrees
    # `m` and `n`, such that `m` divides `n`;
    # * an `n\times m` matrix ``basis`` that maps ``sub_ctx``
    # to an isomorphic subfield in ``sup_ctx``;
    # Compute the `m\times n` matrix of the trace from ``sup_ctx`` to
    # ``sub_ctx``.
    # This matrix is computed as
    # ``embed_dual_to_mono_matrix(_, sub_ctx)``
    # `\times` ``basis``:sup:`t` `\times`
    # ``embed_mono_to_dual_matrix(_, sup_ctx)``.
    # **Note:** if
    # `m=n`, ``basis`` represents a Frobenius, and the result is its
    # inverse matrix.

    void fq_embed_composition_matrix(fmpz_mod_mat_t matrix, const fq_t gen, const fq_ctx_t ctx)
    # Compute the *composition matrix* of ``gen``.
    # For an element `a\in\mathbf{F}_{p^n}`, its composition matrix is the
    # matrix whose columns are `a^0, a^1, \ldots, a^{n-1}`.

    void fq_embed_composition_matrix_sub(fmpz_mod_mat_t matrix, const fq_t gen, const fq_ctx_t ctx, slong trunc)
    # Compute the *composition matrix* of ``gen``, truncated to
    # ``trunc`` columns.

    void fq_embed_mul_matrix(fmpz_mod_mat_t matrix, const fq_t gen, const fq_ctx_t ctx)
    # Compute the *multiplication matrix* of ``gen``.
    # For an element `a` in `\mathbf{F}_{p^n}=\mathbf{F}_p[x]`, its
    # multiplication matrix is the matrix whose columns are `a, ax,
    # \dots, ax^{n-1}`.

    void fq_embed_mono_to_dual_matrix(fmpz_mod_mat_t res, const fq_ctx_t ctx)
    # Compute the change of basis matrix from the monomial basis of
    # ``ctx`` to its dual basis.

    void fq_embed_dual_to_mono_matrix(fmpz_mod_mat_t res, const fq_ctx_t ctx)
    # Compute the change of basis matrix from the dual basis of
    # ``ctx`` to its monomial basis.

    void fq_modulus_pow_series_inv(fmpz_mod_poly_t res, const fq_ctx_t ctx, slong trunc)
    # Compute the power series inverse of the reverse of the modulus of
    # ``ctx`` up to `O(x^\texttt{trunc})`.

    void fq_modulus_derivative_inv(fq_t m_prime, fq_t m_prime_inv, const fq_ctx_t ctx)
    # Compute the derivative ``m_prime`` of the modulus of ``ctx``
    # as an element of ``ctx``, and its inverse ``m_prime_inv``.
