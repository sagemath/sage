# Issue #33153: fflas-ffpack-2.4.3 is missing a return value in one of
# its functions and runs afoul of -Werror=return-type. Compounding the
# problem on openSUSE tumbleweed, the CFLAGS in python's sysconfig
# contain -Werror=return-type and wind up being used to compile this
# extension. To avoid a compilation failure on that platform, we
# temporarily append "-Wno-error=return-type" to those flags.
#
# distutils: extra_compile_args = FFLASFFPACK_CFLAGS -Wno-error=return-type
#
# distutils: include_dirs = FFLASFFPACK_INCDIR
# distutils: libraries = FFLASFFPACK_LIBRARIES
# distutils: library_dirs = FFLASFFPACK_LIBDIR
# distutils: extra_link_args = FFLASFFPACK_LIBEXTRA
# distutils: language = c++ 

from sage.libs.linbox.givaro cimport Modular_double, Modular_float, Dense, Sparse
from sage.libs.linbox.givaro cimport givvector, Poly1Dom
from libcpp.vector cimport vector
from libcpp cimport bool
ctypedef Poly1Dom[Modular_double, Dense] PolynomialRing_Modular_double
ctypedef Poly1Dom[Modular_float, Dense] PolynomialRing_Modular_float

ctypedef givvector[Modular_double.Element] ModDoubleDensePolynomial
ctypedef givvector[Modular_float.Element] ModFloatDensePolynomial

cdef extern from "fflas-ffpack/fflas-ffpack.h" namespace "FFLAS":
    ctypedef enum FFLAS_TRANSPOSE:
        FflasNoTrans
        FflasTrans

    ctypedef enum FFLAS_SIDE:
        FflasRight

    # double
    Modular_double.Element* fgemv (Modular_double F, FFLAS_TRANSPOSE transA,
             size_t nrows, size_t ncols,
             Modular_double.Element alpha, Modular_double.Element* A,
             size_t lda, Modular_double.Element* X, size_t incX,
             Modular_double.Element beta, Modular_double.Element* Y,
             size_t incY)

    Modular_double.Element* fgemm (Modular_double F,
             FFLAS_TRANSPOSE transA, FFLAS_TRANSPOSE transB,
             size_t nrowsA, size_t ncolsB, size_t ncolsA,
             Modular_double.Element alpha, Modular_double.Element* A,
             size_t A_stride, Modular_double.Element* B, int B_stride,
             Modular_double.Element beta, Modular_double.Element* C,
             size_t C_stride)

    Modular_double.Element* pfgemm (Modular_double F,
             FFLAS_TRANSPOSE transA, FFLAS_TRANSPOSE transB,
             size_t nrowsA, size_t ncolsB, size_t ncolsA,
             Modular_double.Element alpha, Modular_double.Element* A,
             size_t A_stride, Modular_double.Element* B, int B_stride,
             Modular_double.Element beta, Modular_double.Element* C,
             size_t C_stride, size_t numthreads)

    # float
    Modular_float.Element* fgemv (Modular_float F, FFLAS_TRANSPOSE transA,
             size_t nrows, size_t ncols,
             Modular_float.Element alpha, Modular_float.Element* A,
             size_t lda, Modular_float.Element* X, size_t incX,
             Modular_float.Element beta, Modular_float.Element* Y,
             size_t incY)

    Modular_float.Element* fgemm (Modular_float F,
             FFLAS_TRANSPOSE transA, FFLAS_TRANSPOSE transB,
             size_t nrowsA, size_t ncolsB, size_t ncolsA,
             Modular_float.Element alpha, Modular_float.Element* A,
             size_t A_stride, Modular_float.Element* B, int B_stride,
             Modular_float.Element beta, Modular_float.Element* C,
             size_t C_stride)

    Modular_float.Element* pfgemm (Modular_float F,
             FFLAS_TRANSPOSE transA, FFLAS_TRANSPOSE transB,
             size_t nrowsA, size_t ncolsB, size_t ncolsA,
             Modular_float.Element alpha, Modular_float.Element* A,
             size_t A_stride, Modular_float.Element* B, int B_stride,
             Modular_float.Element beta, Modular_float.Element* C,
             size_t C_stride, size_t numthreads)

cdef extern from "fflas-ffpack/fflas-ffpack.h" namespace "FFPACK":
    # double
    bint IsSingular (Modular_double F,
                     size_t nrows, size_t ncols, Modular_double.Element* A,
                     size_t A_stride)

    Modular_double.Element* Invert (Modular_double F, size_t order,
                                    Modular_double.Element* A, size_t A_stride, int nullity)

    Modular_double.Element Det (Modular_double F, Modular_double.Element& d, size_t n,
                                Modular_double.Element* A, size_t A_stride)

    Modular_double.Element pDet (Modular_double F, Modular_double.Element& d, size_t n,
                                Modular_double.Element* A, size_t A_stride, size_t numthreads)

    int Rank (Modular_double,
              size_t nrows, size_t ncols,
              Modular_double.Element *A, size_t lda)

    int pRank (Modular_double,
              size_t nrows, size_t ncols,
              Modular_double.Element *A, size_t lda, size_t numthreads)

    size_t ReducedRowEchelonForm (Modular_double F, size_t a, size_t b,
                                  Modular_double.Element* matrix,
                                  size_t s, size_t* P, size_t* Q)

    size_t pReducedRowEchelonForm (Modular_double F, size_t a, size_t b,
                                  Modular_double.Element* matrix,
                                  size_t s, size_t* P, size_t* Q, bool transform, size_t numthreads)

    Modular_double.Element* Solve (Modular_double F, size_t M,
           Modular_double.Element* A, size_t lda,
           Modular_double.Element* x, int incx,
           Modular_double.Element* b, int incb)

    Modular_double.Element* pSolve (Modular_double F, size_t M,
           Modular_double.Element* A, size_t lda,
           Modular_double.Element* x, int incx,
           Modular_double.Element* b, int incb, size_t numthreads)

    void applyP (Modular_double F,
                 FFLAS_SIDE s, FFLAS_TRANSPOSE tr,
                 size_t nr, size_t foo, size_t r,
                 Modular_double.Element* matrix, size_t nc, size_t* Q)

    void MinPoly ( Modular_double& F,
                   vector[Modular_double.Element] minP, size_t N,
                   Modular_double.Element*A, size_t lda)

    void CharPoly ( PolynomialRing_Modular_double& R,
                    ModDoubleDensePolynomial& charp, size_t N,
                    Modular_double.Element* A, size_t lda)

    # float

    bint IsSingular (Modular_float F,
                     size_t nrows, size_t ncols, Modular_float.Element* A,
                     size_t A_stride)

    Modular_float.Element* Invert (Modular_float F, size_t order,
                                   Modular_float.Element* A, size_t A_stride, int nullity)

    Modular_float.Element Det (Modular_float F, Modular_float.Element& d, size_t n,
                               Modular_float.Element* A, size_t A_stride)

    Modular_float.Element pDet (Modular_float F, Modular_float.Element& d, size_t n,
                               Modular_float.Element* A, size_t A_stride, size_t numthreads)

    int Rank (Modular_float,
              size_t nrows, size_t ncols,
              Modular_float.Element *A, size_t lda)

    int pRank (Modular_float,
              size_t nrows, size_t ncols,
              Modular_float.Element *A, size_t lda, size_t numthreads)

    size_t ReducedRowEchelonForm (Modular_float F, size_t a, size_t b,
                                  Modular_float.Element* matrix,
                                  size_t s, size_t* P, size_t* Q)

    size_t pReducedRowEchelonForm (Modular_float F, size_t a, size_t b,
                                  Modular_float.Element* matrix,
                                  size_t s, size_t* P, size_t* Q, bool transform, size_t numthreads)

    Modular_float.Element* Solve (Modular_float F, size_t M,
           Modular_float.Element* A, size_t lda,
           Modular_float.Element* x, int incx,
           Modular_float.Element* b, int incb)

    Modular_float.Element* pSolve (Modular_float F, size_t M,
           Modular_float.Element* A, size_t lda,
           Modular_float.Element* x, int incx,
           Modular_float.Element* b, int incb, size_t numthreads)

    void applyP (Modular_float F,
                 FFLAS_SIDE s, FFLAS_TRANSPOSE tr,
                 size_t nr, size_t foo, size_t r,
                 Modular_float.Element* matrix, size_t nc, size_t* Q)

    void MinPoly ( Modular_float F,
                   vector[Modular_float.Element] minP, size_t N,
                   Modular_float.Element* A, size_t lda)

    void CharPoly ( PolynomialRing_Modular_float& F,
                    ModFloatDensePolynomial& charp, size_t N,
                    Modular_float.Element* A, size_t lda )
