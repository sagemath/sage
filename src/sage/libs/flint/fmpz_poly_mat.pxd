# distutils: libraries = flint
# distutils: depends = flint/fmpz_poly_mat.h

from sage.libs.flint.types cimport fmpz_poly_mat_t, fmpz_poly_t, fmpz_t, slong, fmpz_mat_t

# flint/fmpz_poly_mat.h
cdef extern from "flint_wrap.h":

    void fmpz_poly_mat_init(fmpz_poly_mat_t mat, slong rows, slong cols)
    void fmpz_poly_mat_init_set(fmpz_poly_mat_t mat, const fmpz_poly_mat_t src)
    void fmpz_poly_mat_clear(fmpz_poly_mat_t mat)
    fmpz_poly_t fmpz_poly_mat_entry(fmpz_poly_mat_t mat, long i, long j)
    slong fmpz_poly_mat_nrows(const fmpz_poly_mat_t mat)
    slong fmpz_poly_mat_ncols(const fmpz_poly_mat_t mat)

    void fmpz_poly_mat_set(fmpz_poly_mat_t mat1, const fmpz_poly_mat_t mat2)

    void fmpz_poly_mat_swap(fmpz_poly_mat_t mat1, fmpz_poly_mat_t mat2)

    void fmpz_poly_mat_transpose(fmpz_poly_mat_t B, const fmpz_poly_mat_t A)

    void fmpz_poly_mat_evaluate_fmpz(fmpz_mat_t B, const fmpz_poly_mat_t A,
                                                               const fmpz_t x)

    void fmpz_poly_mat_trace(fmpz_poly_t trace, const fmpz_poly_mat_t mat)

    void fmpz_poly_mat_det(fmpz_poly_t det, const fmpz_poly_mat_t A)
