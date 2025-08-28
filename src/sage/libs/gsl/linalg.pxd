# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_linalg.h":

  cdef enum gsl_linalg_matrix_mod_t:
    GSL_LINALG_MOD_NONE = 0
    GSL_LINALG_MOD_TRANSPOSE = 1
    GSL_LINALG_MOD_CONJUGATE = 2

  int gsl_linalg_matmult ( gsl_matrix * A,  gsl_matrix * B, gsl_matrix * C)

  int gsl_linalg_matmult_mod ( gsl_matrix * A, gsl_linalg_matrix_mod_t modA,  gsl_matrix * B, gsl_linalg_matrix_mod_t modB, gsl_matrix * C)

  int gsl_linalg_exponential_ss(gsl_matrix * A, gsl_matrix * eA, gsl_mode_t mode)

  # Householder Transformations
  double gsl_linalg_householder_transform (gsl_vector * v)

  gsl_complex gsl_linalg_complex_householder_transform (gsl_vector_complex * v)

  int gsl_linalg_householder_hm (double tau,  gsl_vector * v, gsl_matrix * A)

  int gsl_linalg_householder_mh (double tau,  gsl_vector * v, gsl_matrix * A)

  int gsl_linalg_householder_hv (double tau,  gsl_vector * v, gsl_vector * w)

  int gsl_linalg_householder_hm1 (double tau, gsl_matrix * A)

  int gsl_linalg_complex_householder_hm (gsl_complex tau,  gsl_vector_complex * v, gsl_matrix_complex * A)

  int gsl_linalg_complex_householder_hv (gsl_complex tau,  gsl_vector_complex * v, gsl_vector_complex * w)

  # Singular Value Decomposition
  int  gsl_linalg_SV_decomp(gsl_matrix * A, gsl_matrix * V, gsl_vector * S, gsl_vector * work)

  int  gsl_linalg_SV_decomp_mod(gsl_matrix * A, gsl_matrix * X, gsl_matrix * V, gsl_vector * S, gsl_vector * work)

  int  gsl_linalg_SV_decomp_jacobi(gsl_matrix * A, gsl_matrix * V, gsl_vector * S)

  int  gsl_linalg_SV_solve(gsl_matrix * U, gsl_matrix * V, gsl_vector * S, gsl_vector * b, gsl_vector * x)


  # LU Decomposition, Gaussian elimination with partial pivoting
  int  gsl_linalg_LU_decomp(gsl_matrix * A, gsl_permutation * p, int *signum)

  int  gsl_linalg_LU_solve(gsl_matrix * LU, gsl_permutation * p, gsl_vector * b, gsl_vector * x)

  int  gsl_linalg_LU_svx(gsl_matrix * LU, gsl_permutation * p, gsl_vector * x)

  int  gsl_linalg_LU_refine(gsl_matrix * A, gsl_matrix * LU, gsl_permutation * p, gsl_vector * b, gsl_vector * x, gsl_vector * residual)

  int  gsl_linalg_LU_invert(gsl_matrix * LU, gsl_permutation * p, gsl_matrix * inverse)

  double  gsl_linalg_LU_det(gsl_matrix * LU, int signum)

  double  gsl_linalg_LU_lndet(gsl_matrix * LU)

  int  gsl_linalg_LU_sgndet(gsl_matrix * LU, int signum)


  # Complex LU Decomposition

  int  gsl_linalg_complex_LU_decomp(gsl_matrix_complex * A, gsl_permutation * p, int *signum)

  int  gsl_linalg_complex_LU_solve(gsl_matrix_complex * LU, gsl_permutation * p, gsl_vector_complex * b, gsl_vector_complex * x)

  int  gsl_linalg_complex_LU_svx(gsl_matrix_complex * LU, gsl_permutation * p, gsl_vector_complex * x)

  int  gsl_linalg_complex_LU_refine(gsl_matrix_complex * A, gsl_matrix_complex * LU, gsl_permutation * p, gsl_vector_complex * b, gsl_vector_complex * x, gsl_vector_complex * residual)

#  int  gsl_complex_linalg_LU_invert(gsl_matrix_complex * LU, gsl_permutation * p, gsl_matrix_complex * inverse)
  int  gsl_linalg_complex_LU_invert(gsl_matrix_complex * LU, gsl_permutation * p, gsl_matrix_complex * inverse)

  gsl_complex  gsl_linalg_complex_LU_det(gsl_matrix_complex * LU, int signum)

  double  gsl_linalg_complex_LU_lndet(gsl_matrix_complex * LU)

  gsl_complex  gsl_linalg_complex_LU_sgndet(gsl_matrix_complex * LU, int signum)


  # QR decomposition

  int  gsl_linalg_QR_decomp(gsl_matrix * A, gsl_vector * tau)

  int  gsl_linalg_QR_solve(gsl_matrix * QR, gsl_vector * tau, gsl_vector * b, gsl_vector * x)

  int  gsl_linalg_QR_svx(gsl_matrix * QR, gsl_vector * tau, gsl_vector * x)

  int  gsl_linalg_QR_lssolve(gsl_matrix * QR, gsl_vector * tau, gsl_vector * b, gsl_vector * x, gsl_vector * residual)

  int  gsl_linalg_QR_QRsolve(gsl_matrix * Q, gsl_matrix * R, gsl_vector * b, gsl_vector * x)

  int  gsl_linalg_QR_Rsolve(gsl_matrix * QR, gsl_vector * b, gsl_vector * x)

  int  gsl_linalg_QR_Rsvx(gsl_matrix * QR, gsl_vector * x)

  int  gsl_linalg_QR_update(gsl_matrix * Q, gsl_matrix * R, gsl_vector * w, gsl_vector * v)

  int  gsl_linalg_QR_QTvec(gsl_matrix * QR, gsl_vector * tau, gsl_vector * v)

  int  gsl_linalg_QR_Qvec(gsl_matrix * QR, gsl_vector * tau, gsl_vector * v)

  int  gsl_linalg_QR_unpack(gsl_matrix * QR, gsl_vector * tau, gsl_matrix * Q, gsl_matrix * R)

  int  gsl_linalg_R_solve(gsl_matrix * R, gsl_vector * b, gsl_vector * x)

  int  gsl_linalg_R_svx(gsl_matrix * R, gsl_vector * x)


  # Q R P^T decomposition

  int  gsl_linalg_QRPT_decomp(gsl_matrix * A, gsl_vector * tau, gsl_permutation * p, int *signum, gsl_vector * norm)

  int  gsl_linalg_QRPT_decomp2(gsl_matrix * A, gsl_matrix * q, gsl_matrix * r, gsl_vector * tau, gsl_permutation * p, int *signum, gsl_vector * norm)

  int  gsl_linalg_QRPT_solve(gsl_matrix * QR, gsl_vector * tau, gsl_permutation * p, gsl_vector * b, gsl_vector * x)

  int  gsl_linalg_QRPT_svx(gsl_matrix * QR, gsl_vector * tau, gsl_permutation * p, gsl_vector * x)

  int  gsl_linalg_QRPT_QRsolve(gsl_matrix * Q, gsl_matrix * R, gsl_permutation * p, gsl_vector * b, gsl_vector * x)

  int  gsl_linalg_QRPT_Rsolve(gsl_matrix * QR, gsl_permutation * p, gsl_vector * b, gsl_vector * x)

  int  gsl_linalg_QRPT_Rsvx(gsl_matrix * QR, gsl_permutation * p, gsl_vector * x)

  int  gsl_linalg_QRPT_update(gsl_matrix * Q, gsl_matrix * R, gsl_permutation * p, gsl_vector * u, gsl_vector * v)


  # Cholesky Decomposition
  int  gsl_linalg_cholesky_decomp(gsl_matrix * A)

  int  gsl_linalg_cholesky_solve(gsl_matrix * cholesky, gsl_vector * b, gsl_vector * x)

  int  gsl_linalg_cholesky_svx(gsl_matrix * cholesky, gsl_vector * x)

  # Symmetric to symmetric tridiagonal decomposition
  int  gsl_linalg_symmtd_decomp(gsl_matrix * A, gsl_vector * tau)

  int  gsl_linalg_symmtd_unpack(gsl_matrix * A, gsl_vector * tau, gsl_matrix * Q, gsl_vector * diag, gsl_vector * subdiag)

  int  gsl_linalg_symmtd_unpack_T(gsl_matrix * A, gsl_vector * diag, gsl_vector * subdiag)


  # Hermitian to symmetric tridiagonal decomposition
  int  gsl_linalg_hermtd_decomp(gsl_matrix_complex * A, gsl_vector_complex * tau)

  int  gsl_linalg_hermtd_unpack(gsl_matrix_complex * A, gsl_vector_complex * tau, gsl_matrix_complex * Q, gsl_vector * diag, gsl_vector * subdiag)

  int  gsl_linalg_hermtd_unpack_T(gsl_matrix_complex * A, gsl_vector * diag, gsl_vector * subdiag)

  # Linear Solve Using Householder Transformations
  int  gsl_linalg_HH_solve(gsl_matrix * A, gsl_vector * b, gsl_vector * x)

  int  gsl_linalg_HH_svx(gsl_matrix * A, gsl_vector * x)

  # Linear solve for a symmetric tridiagonal system
  int  gsl_linalg_solve_symm_tridiag(gsl_vector * diag, gsl_vector * e, gsl_vector * b, gsl_vector * x)

  # Linear solve for a nonsymmetric tridiagonal system
  int gsl_linalg_solve_tridiag ( gsl_vector * diag,  gsl_vector * abovediag,  gsl_vector * belowdiag,  gsl_vector * b, gsl_vector * x)

  # Linear solve for a symmetric cyclic tridiagonal system
  int  gsl_linalg_solve_symm_cyc_tridiag(gsl_vector * diag, gsl_vector * e, gsl_vector * b, gsl_vector * x)


  # Linear solve for a nonsymmetric cyclic tridiagonal system
  int gsl_linalg_solve_cyc_tridiag ( gsl_vector * diag,  gsl_vector * abovediag,  gsl_vector * belowdiag,  gsl_vector * b, gsl_vector * x)

  # Bidiagonal decomposition
  int  gsl_linalg_bidiag_decomp(gsl_matrix * A, gsl_vector * tau_U, gsl_vector * tau_V)

  int  gsl_linalg_bidiag_unpack(gsl_matrix * A, gsl_vector * tau_U, gsl_matrix * U, gsl_vector * tau_V, gsl_matrix * V, gsl_vector * diag, gsl_vector * superdiag)

  int  gsl_linalg_bidiag_unpack2(gsl_matrix * A, gsl_vector * tau_U, gsl_vector * tau_V, gsl_matrix * V)

  int  gsl_linalg_bidiag_unpack_B(gsl_matrix * A, gsl_vector * diag, gsl_vector * superdiag)


  # Balancing
  int gsl_linalg_balance_columns (gsl_matrix * A, gsl_vector * D)

