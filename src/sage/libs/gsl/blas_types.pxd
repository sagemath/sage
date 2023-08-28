# sage_setup: distribution = sagemath-modules
# distutils: include_dirs = GSL_INCDIR
cdef extern from "gsl/gsl_cblas.h":
  cdef enum CBLAS_ORDER:
    CblasRowMajor=101
    CblasColMajor=102

  cdef enum CBLAS_TRANSPOSE:
    CblasNoTrans=111
    CblasTrans=112
    CblasConjTrans=113

  cdef enum CBLAS_UPLO:
    CblasUpper=121
    CblasLower=122

  cdef enum CBLAS_DIAG:
    CblasNonUnit=131
    CblasUnit=132

  cdef enum CBLAS_SIDE:
    CblasLeft=141
    CblasRight=142

ctypedef size_t CBLAS_INDEX

cdef extern from "gsl/gsl_blas_types.h":
  ctypedef  CBLAS_INDEX  CBLAS_INDEX_t
  ctypedef  CBLAS_ORDER       CBLAS_ORDER_t
  ctypedef  CBLAS_TRANSPOSE   CBLAS_TRANSPOSE_t
  ctypedef  CBLAS_UPLO        CBLAS_UPLO_t
  ctypedef  CBLAS_DIAG        CBLAS_DIAG_t
  ctypedef  CBLAS_SIDE        CBLAS_SIDE_t

