# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_blas.h":

  # Level 1
  int  gsl_blas_ddot(gsl_vector * x, gsl_vector * y, double * result)

  int  gsl_blas_zdotu(gsl_vector_complex * x, gsl_vector_complex * y, gsl_complex * dotu)

  int  gsl_blas_zdotc(gsl_vector_complex * x, gsl_vector_complex * y, gsl_complex * dotc)

  double  gsl_blas_dnrm2(gsl_vector * x)

  double  gsl_blas_dasum(gsl_vector * x)

  double  gsl_blas_dznrm2(gsl_vector_complex * x)

  double  gsl_blas_dzasum(gsl_vector_complex * x)

  CBLAS_INDEX_t  gsl_blas_idamax(gsl_vector * x)

  CBLAS_INDEX_t  gsl_blas_izamax(gsl_vector_complex * x)

  int  gsl_blas_dswap(gsl_vector * x, gsl_vector * y)

  int  gsl_blas_dcopy(gsl_vector * x, gsl_vector * y)

  int  gsl_blas_daxpy(double alpha, gsl_vector * x, gsl_vector * y)

  int  gsl_blas_zswap(gsl_vector_complex * x, gsl_vector_complex * y)

  int  gsl_blas_zcopy(gsl_vector_complex * x, gsl_vector_complex * y)

  int  gsl_blas_zaxpy(gsl_complex alpha, gsl_vector_complex * x, gsl_vector_complex * y)

  int  gsl_blas_drotg(double a[], double b[], double c[], double s[])

  int  gsl_blas_drotmg(double d1[], double d2[], double b1[], double b2, double P[])

  int  gsl_blas_drot(gsl_vector * x, gsl_vector * y, double c, double s)

  int  gsl_blas_drotm(gsl_vector * x, gsl_vector * y, double P[])

  void  gsl_blas_dscal(double alpha, gsl_vector * x)

  void  gsl_blas_zscal(gsl_complex alpha, gsl_vector_complex * x)

  void  gsl_blas_zdscal(double alpha, gsl_vector_complex * x)


  # Level 2
  int  gsl_blas_dgemv(CBLAS_TRANSPOSE_t TransA, double alpha, gsl_matrix * A, gsl_vector * x, double beta, gsl_vector * y)

  int  gsl_blas_dtrmv(CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, gsl_matrix * A, gsl_vector * x)

  int  gsl_blas_dtrsv(CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, gsl_matrix * A, gsl_vector * x)

  int  gsl_blas_zgemv (CBLAS_TRANSPOSE_t TransA,  gsl_complex alpha,  gsl_matrix_complex * A,  gsl_vector_complex * X,  gsl_complex beta, gsl_vector_complex * Y)

  int  gsl_blas_ztrmv (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,  gsl_matrix_complex * A, gsl_vector_complex * X)

  int  gsl_blas_ztrsv (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,  gsl_matrix_complex * A, gsl_vector_complex *X)

  int  gsl_blas_dsymv(CBLAS_UPLO_t Uplo, double alpha, gsl_matrix * A, gsl_vector * x, double beta, gsl_vector * y)

  int  gsl_blas_dger(double alpha, gsl_vector * x, gsl_vector * y, gsl_matrix * A)

  int  gsl_blas_dsyr(CBLAS_UPLO_t Uplo, double alpha, gsl_vector * x, gsl_matrix * A)

  int  gsl_blas_dsyr2(CBLAS_UPLO_t Uplo, double alpha, gsl_vector * x, gsl_vector * y, gsl_matrix * A)

  int  gsl_blas_zhemv (CBLAS_UPLO_t Uplo,  gsl_complex alpha,  gsl_matrix_complex * A,  gsl_vector_complex * X,  gsl_complex beta, gsl_vector_complex * Y)

  int  gsl_blas_zgeru ( gsl_complex alpha,  gsl_vector_complex * X,  gsl_vector_complex * Y, gsl_matrix_complex * A)

  int  gsl_blas_zgerc ( gsl_complex alpha,  gsl_vector_complex * X,  gsl_vector_complex * Y, gsl_matrix_complex * A)

  int  gsl_blas_zher (CBLAS_UPLO_t Uplo, double alpha,  gsl_vector_complex * X, gsl_matrix_complex * A)

  int  gsl_blas_zher2 (CBLAS_UPLO_t Uplo,  gsl_complex alpha,  gsl_vector_complex * X,  gsl_vector_complex * Y, gsl_matrix_complex * A)


  # Level 3
  int  gsl_blas_dgemm(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, gsl_matrix * A, gsl_matrix * B, double beta, gsl_matrix * C)

  int  gsl_blas_dsymm(CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, double alpha, gsl_matrix * A, gsl_matrix * B, double beta, gsl_matrix * C)

  int  gsl_blas_dsyrk(CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, double alpha, gsl_matrix * A, double beta, gsl_matrix * C)

  int  gsl_blas_dsyr2k(CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, double alpha, gsl_matrix * A, gsl_matrix * B, double beta, gsl_matrix * C)

  int  gsl_blas_dtrmm(CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, double alpha, gsl_matrix * A, gsl_matrix * B)

  int  gsl_blas_dtrsm(CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, double alpha, gsl_matrix * A, gsl_matrix * B)

  int  gsl_blas_zgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB,  gsl_complex alpha,  gsl_matrix_complex * A,  gsl_matrix_complex * B,  gsl_complex beta, gsl_matrix_complex * C)

  int  gsl_blas_zsymm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo,  gsl_complex alpha,  gsl_matrix_complex * A,  gsl_matrix_complex * B,  gsl_complex beta, gsl_matrix_complex * C)

  int  gsl_blas_zsyrk (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans,  gsl_complex alpha,  gsl_matrix_complex * A,  gsl_complex beta, gsl_matrix_complex * C)

  int  gsl_blas_zsyr2k (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans,  gsl_complex alpha,  gsl_matrix_complex * A,  gsl_matrix_complex * B,  gsl_complex beta, gsl_matrix_complex *C)

  int  gsl_blas_ztrmm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,  gsl_complex alpha,  gsl_matrix_complex * A, gsl_matrix_complex * B)

  int  gsl_blas_ztrsm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,  gsl_complex alpha,  gsl_matrix_complex * A, gsl_matrix_complex * B)

  int  gsl_blas_zhemm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo,  gsl_complex alpha,  gsl_matrix_complex * A,  gsl_matrix_complex * B,  gsl_complex beta, gsl_matrix_complex * C)

  int  gsl_blas_zherk (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, double alpha,  gsl_matrix_complex * A, double beta, gsl_matrix_complex * C)

  int  gsl_blas_zher2k (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans,  gsl_complex alpha,  gsl_matrix_complex * A,  gsl_matrix_complex * B, double beta, gsl_matrix_complex * C)
