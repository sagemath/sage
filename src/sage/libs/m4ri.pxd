# sage_setup: distribution = sagemath-linbox
# distutils: extra_compile_args = -std=c++11
# distutils: language = c++

cdef extern from "m4ri/m4ri.h":
    ctypedef int rci_t
    ctypedef int wi_t
    ctypedef unsigned long long m4ri_word "word"
    ctypedef int BIT

    ctypedef struct mzd_t:
        rci_t nrows
        rci_t ncols
        wi_t width

    ctypedef struct mzp_t:
        rci_t *values
        rci_t size

    cdef int m4ri_radix

    ##############
    # Maintenance
    ##############

    # builds all gray codes up to a certain size
    cdef void m4ri_build_all_codes()
    cdef void m4ri_destroy_all_codes()

    ##############
    # Constructors
    ##############

    # create empty matrix
    cdef mzd_t *mzd_init(rci_t , rci_t)

    # create the identity permutation
    cdef mzp_t *mzp_init(rci_t)

    # free memory for the matrix
    cdef void mzd_free(mzd_t *)

    # free memory for the permutation
    cdef void mzp_free(mzp_t *)

    # filled uniformly random
    cdef void mzd_randomize(mzd_t *)

    # identity matrix if i%2
    cdef void mzd_set_ui(mzd_t *, unsigned int )

    # [A],[B] -> [AB]
    cdef mzd_t *mzd_concat(mzd_t *, mzd_t *, mzd_t *)

    # [A],[B] -> | A |
    #            | B |
    cdef mzd_t *mzd_stack(mzd_t *, mzd_t *, mzd_t *)

    # returns a submatrix from a
    cdef mzd_t *mzd_submatrix(mzd_t *, mzd_t *, rci_t lowr, rci_t lowc, rci_t highr, rci_t highc)

    # return a matrix window to A
    cdef mzd_t *mzd_init_window(mzd_t *, rci_t lowr, rci_t lowc, rci_t highr, rci_t highc)

    cdef void mzd_free_window(mzd_t *)

    # deep copy
    cdef mzd_t *mzd_copy(mzd_t *, mzd_t *)

    # printing
    cdef void mzd_print(mzd_t *)

    ##############
    # Bit Level IO
    ##############

    # set BIT
    cdef void mzd_write_bit( mzd_t *m, rci_t row, rci_t col, BIT value)

    # get BIT
    cdef BIT mzd_read_bit( mzd_t *m, rci_t row, rci_t col )

    # get BITs (n<=64)
    cdef m4ri_word mzd_read_bits( mzd_t *m, rci_t row, rci_t col, int n)

    #####################
    # Row/Column Based IO
    #####################

    cdef m4ri_word *mzd_row(mzd_t *, rci_t)

    cdef void mzd_row_swap(mzd_t *, rci_t, rci_t)

    cdef void mzd_col_swap(mzd_t *, rci_t, rci_t)

    cdef void mzd_row_clear_offset(mzd_t *m, rci_t, rci_t)

    cdef void mzd_row_add_offset(mzd_t *m, rci_t, rci_t, rci_t)

    ############
    # Arithmetic
    ############

    # matrix addition
    cdef mzd_t *mzd_add(mzd_t *, mzd_t *, mzd_t *)

    # naive cubic matrix multiply
    cdef mzd_t *mzd_mul_naive(mzd_t *, mzd_t *, mzd_t *)

    # naive cubic matrix multiply (b is pre-transposed)
    cdef mzd_t *_mzd_mul_naive(mzd_t *, mzd_t *, mzd_t *, int)

    # matrix multiply using Gray codes
    cdef mzd_t *mzd_mul_m4rm(mzd_t *, mzd_t *, mzd_t *, int k)

    # matrix multiply using Gray codes (transposed)
    cdef mzd_t *mzd_mul_m4rm_t(mzd_t *, mzd_t *, mzd_t *, int k)

    # matrix multiply and addition using Gray codes: C = C + AB
    cdef mzd_t *mzd_addmul_m4rm(mzd_t *, mzd_t *, mzd_t *, int k)

    # matrix multiplication via Strassen's formula
    cdef mzd_t *mzd_mul(mzd_t *, mzd_t *, mzd_t *, int cutoff)

    # C = C + AB via Strassen's formula
    cdef mzd_t *mzd_addmul(mzd_t *, mzd_t *, mzd_t *, int cutoff)

    # equality testing
    cdef int mzd_equal(mzd_t *, mzd_t *)

    # returns -1,0,1
    cdef int mzd_cmp(mzd_t *, mzd_t *)

    # transpose
    cdef mzd_t *mzd_transpose(mzd_t *, mzd_t *)

    # density with given resolution
    cdef double mzd_density(mzd_t *, int resolution)

    ########################
    # LAPACK Level Functions
    ########################

    # cubic Gaussian elimination
    cdef int mzd_echelonize_naive(mzd_t *, int full)

    # row echelon form using Gray codes
    cdef int mzd_echelonize_m4ri(mzd_t *m, int full, int k)

    # reduced row echelon form from upper triangular form
    cdef void mzd_top_echelonize_m4ri(mzd_t *m, int k)

    # heuristic choice of algorithms
    cdef int mzd_echelonize(mzd_t *m, int full)

    # matrix inversion using Gray codes
    cdef mzd_t *mzd_inv_m4ri(mzd_t *dst, mzd_t *src, int k)

    # asymptotically fast PLUQ factorization
    cdef long mzd_pluq(mzd_t *A, mzp_t *P, mzp_t *Q, int cutoff)

    # PLUQ factorization using Gray codes
    cdef long _mzd_pluq_russian(mzd_t *A, mzp_t *P, mzp_t *Q, int k)

    # cubic PLUQ factorization
    cdef long _mzd_pluq_naive(mzd_t *A, mzp_t *P, mzp_t *Q)

    # asymptotically fast PLS factorization
    cdef long mzd_ple(mzd_t *A, mzp_t *P, mzp_t *Q, int cutoff)

    # PLS factorization using Gray codes
    cdef long _mzd_ple_russian(mzd_t *A, mzp_t *P, mzp_t *Q, int k)

    # cubic PLS factorization
    cdef long _mzd_ple_naive(mzd_t *A, mzp_t *P, mzp_t *Q)

    # reduced row echelon form using PLUQ factorization
    cdef long mzd_echelonize_pluq(mzd_t *A, int full)

    # reduced row echelon form using PLUQ factorization
    cdef mzd_t *mzd_kernel_left_pluq(mzd_t *A, int cutoff)

    ########################
    # Bit operations
    ########################

    cdef m4ri_word __M4RI_LEFT_BITMASK(int)

    cdef m4ri_word m4ri_swap_bits(m4ri_word)

    ##################################
    # Internal Functions for Debugging
    ##################################

    cdef void mzd_clear_bits(mzd_t *m, int x, int y, int n)

