##############################################################################
#       Copyright (C) 2010 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.libs.m4ri cimport mzd_t, m4ri_word


cdef extern from "m4rie/m4rie.h":
    ctypedef struct gf2e:
        int degree
        m4ri_word minpoly

        m4ri_word *pow_gen
        m4ri_word *red

        m4ri_word (*inv)(gf2e *ff, m4ri_word a)
        m4ri_word (*mul)(gf2e *ff, m4ri_word a, m4ri_word b)

    gf2e *gf2e_init(m4ri_word minpoly)
    void gf2e_free(gf2e *ff)

#cdef extern from "m4rie/mzed.h":
    ctypedef struct mzed_t:
        mzd_t *x
        gf2e *finite_field
        int nrows
        int ncols
        int w

    ctypedef int const_int "const int"
    ctypedef size_t const_size_t "const size_t"
    ctypedef mzed_t const_mzed_t "const mzed_t"

    mzed_t *mzed_init(gf2e *, size_t m, size_t n)

    void mzed_free(mzed_t *)

    int mzed_read_elem(const_mzed_t *M, const_size_t row, const_size_t col)

    void mzed_write_elem(mzed_t *, const_size_t row, const_size_t col, const_int elem)

    void mzed_set_ui(mzed_t *A, m4ri_word value)

    mzed_t *mzed_copy(mzed_t *o, const_mzed_t *i)

    int mzed_cmp(mzed_t *l, mzed_t *r)

    mzed_t *mzed_randomize(mzed_t *)

    mzed_t *mzed_add(mzed_t *, mzed_t *, mzed_t *)

    size_t mzed_echelonize_naive(mzed_t *, size_t)

    void mzed_add_elem(mzed_t *a, const_size_t row, const_size_t col, const_int elem)

    void mzed_add_multiple_of_row(mzed_t *A, size_t ar, mzed_t *B, size_t br, m4ri_word x, size_t start_col)

    void mzed_rescale_row(mzed_t *A, size_t r, size_t c, m4ri_word x)

    void mzed_row_swap(mzed_t *M, const_size_t rowa, const_size_t rowb)

    void mzed_copy_row(mzed_t* B, size_t i, const_mzed_t* A, size_t j)

    void mzed_col_swap(mzed_t *M, const_size_t cola, const_size_t colb)

    void mzed_row_add(mzed_t *M, const_size_t sourcerow, const_size_t destrow)

    size_t mzed_first_zero_row(mzed_t *A)

    int mzed_is_zero(mzed_t *A)

    void mzed_row_clear_offset(mzed_t *M, const_size_t row, const_size_t coloffset)

    mzed_t *mzed_concat(mzed_t *C, const_mzed_t *A, const_mzed_t *B)

    mzed_t *mzed_stack(mzed_t *C, const_mzed_t *A, const_mzed_t *B)

    mzed_t *mzed_submatrix(mzed_t *S, const_mzed_t *M, size_t lowr, size_t lowc, size_t highr, size_t highc)

    mzed_t *mzed_mul(mzed_t *C, const_mzed_t *A, const_mzed_t *B)

    mzed_t *mzed_mul_naive(mzed_t *C, const_mzed_t *A, const_mzed_t *B)

    mzed_t *mzed_mul_scalar(mzed_t *C, m4ri_word a, const_mzed_t *B)

    # TODO: not implemented yet in m4rie
    mzed_t *mzed_transpose(mzed_t *DST, const_mzed_t *A )

    void mzed_print(const_mzed_t *M)

    mzed_t *mzed_invert_newton_john(mzed_t *A, mzed_t *B)

    # TODO: not implemented yet in m4rie
    double mzed_density(mzed_t *A, int res)

    # TODO: not implemented yet in m4rie
    double _mzed_density(mzed_t *A, int res, size_t r, size_t c)

#cdef extern from "m4rie/newton_john.h":
    size_t mzed_echelonize_newton_john(mzed_t *, size_t)

    mzed_t *mzed_mul_newton_john(mzed_t *, mzed_t *, mzed_t *)

#cdef extern from "m4rie/echelonform.h":
    size_t mzed_echelonize(mzed_t *, size_t)

    size_t mzed_echelonize_ple(mzed_t *, size_t)

#cdef extern from "m4rie/strassen.h":
    mzed_t *mzed_mul_strassen(mzed_t *, mzed_t *, mzed_t *, size_t cutoff)

    size_t _mzed_strassen_cutoff(mzed_t *C, mzed_t *A, mzed_t *B)

#cdef extern from "m4rie/mzd_slice.h":

    int __M4RIE_MAX_KARATSUBA_DEGREE

    ctypedef struct mzd_slice_t:
        mzd_t *x[16]
        gf2e *finite_field
        int nrows
        int ncols
        int depth

    mzd_slice_t *mzd_slice_init(gf2e *ff, size_t m, size_t n)

    void mzd_slice_free(mzd_slice_t *A)

    mzed_t *mzed_cling(mzed_t *A, mzd_slice_t *Z)

    mzd_slice_t *mzed_slice(mzd_slice_t *A, mzed_t *Z)

    mzd_slice_t *mzd_slice_concat(mzd_slice_t *C, mzd_slice_t *A, mzd_slice_t *B)

    mzd_slice_t *mzd_slice_stack(mzd_slice_t *C, mzd_slice_t *A, mzd_slice_t *B)

    mzd_slice_t *mzd_slice_submatrix(mzd_slice_t *S, mzd_slice_t *A, size_t lowr, size_t lowc, size_t highr, size_t highc)

    mzd_slice_t *mzd_slice_init_window(mzd_slice_t *A, size_t lowr, size_t lowc, size_t highr, size_t highc)

    void mzd_slice_free_window(mzd_slice_t *A)

    mzd_slice_t *_mzd_slice_add(mzd_slice_t *C, mzd_slice_t *A, mzd_slice_t *B)

    mzd_slice_t *mzd_slice_add(mzd_slice_t *C, mzd_slice_t *A, mzd_slice_t *B)

    void mzd_slice_randomize(mzd_slice_t *A)

    mzd_slice_t *mzd_slice_copy(mzd_slice_t *B, mzd_slice_t *A)

    void mzd_slice_set_ui(mzd_slice_t *A, m4ri_word value)

    m4ri_word mzd_slice_read_elem(mzd_slice_t *A, size_t row, size_t col)

    void mzd_slice_add_elem(mzd_slice_t *A, size_t row, size_t col, m4ri_word elem)

    void mzd_slice_write_elem(mzd_slice_t *A, size_t row, size_t col, m4ri_word elem)

    int mzd_slice_cmp(mzd_slice_t *A, mzd_slice_t *B)

    int mzd_slice_is_zero(mzd_slice_t *A)

    void mzd_slice_rescale_row(mzd_slice_t *A, size_t r, size_t c, m4ri_word x)

    void mzd_slice_row_swap(mzd_slice_t *A, size_t rowa, size_t rowb)

    void mzd_slice_copy_row(mzd_slice_t* B, size_t i, mzd_slice_t* A, size_t j)

    void mzd_slice_col_swap(mzd_slice_t *A, size_t cola, size_t colb)

    void mzd_slice_row_add(mzd_slice_t *A, size_t sourcerow, size_t destrow)

    void mzd_slice_row_clear_offset(mzd_slice_t *A, size_t row, size_t coloffset)

    void mzd_slice_print(mzd_slice_t *A)

    mzd_slice_t *_mzed_slice2(mzd_slice_t *A, mzed_t *Z)

    mzd_slice_t *mzed_slice2(mzd_slice_t *A, mzed_t *Z)

    mzd_slice_t *_mzed_slice4(mzd_slice_t *A, mzed_t *Z)

    mzed_t *_mzed_cling2(mzed_t *A, mzd_slice_t *Z)

    mzed_t* mzed_cling2(mzed_t *A, mzd_slice_t *Z)

    mzed_t *_mzed_cling4(mzed_t *A, mzd_slice_t *Z)

    mzed_t* mzed_cling4(mzed_t *A, mzd_slice_t *Z)

    mzed_t *mzed_mul_karatsuba(mzed_t *C, mzed_t *A, mzed_t *B)

    mzed_t *mzed_addmul_karatsuba(mzed_t *C, mzed_t *A, mzed_t *B)

    mzed_t *_mzed_mul_karatsuba(mzed_t *C, mzed_t *A, mzed_t *B)

    mzd_slice_t *_mzd_slice_mul_karatsuba2(mzd_slice_t *C, mzd_slice_t *A, mzd_slice_t *B)

    mzd_slice_t *_mzd_slice_mul_karatsuba3(mzd_slice_t *C, mzd_slice_t *A, mzd_slice_t *B)
