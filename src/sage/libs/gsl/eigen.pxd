# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from .types cimport *

cdef extern from "gsl/gsl_eigen.h":

  ctypedef enum gsl_eigen_sort_t:
    GSL_EIGEN_SORT_VAL_ASC,
    GSL_EIGEN_SORT_VAL_DESC,
    GSL_EIGEN_SORT_ABS_ASC,
    GSL_EIGEN_SORT_ABS_DESC
  ctypedef struct gsl_eigen_symm_workspace
  ctypedef struct gsl_eigen_symmv_workspace
  ctypedef struct gsl_eigen_herm_workspace
  ctypedef struct gsl_eigen_hermv_workspace

  gsl_eigen_symm_workspace *  gsl_eigen_symm_alloc(size_t n)

  void  gsl_eigen_symm_free(gsl_eigen_symm_workspace * w)

  int  gsl_eigen_symm(gsl_matrix * A, gsl_vector * eval, gsl_eigen_symm_workspace * w)

  gsl_eigen_symmv_workspace *  gsl_eigen_symmv_alloc(size_t n)

  void  gsl_eigen_symmv_free(gsl_eigen_symmv_workspace * w)

  int  gsl_eigen_symmv(gsl_matrix * A, gsl_vector * eval, gsl_matrix * evec, gsl_eigen_symmv_workspace * w)

  gsl_eigen_herm_workspace *  gsl_eigen_herm_alloc(size_t n)

  void  gsl_eigen_herm_free(gsl_eigen_herm_workspace * w)

  int  gsl_eigen_herm(gsl_matrix_complex * A, gsl_vector * eval, gsl_eigen_herm_workspace * w)

  gsl_eigen_hermv_workspace *  gsl_eigen_hermv_alloc(size_t n)

  void  gsl_eigen_hermv_free(gsl_eigen_hermv_workspace * w)

  int  gsl_eigen_hermv(gsl_matrix_complex * A, gsl_vector * eval, gsl_matrix_complex * evec, gsl_eigen_hermv_workspace * w)

  int  gsl_eigen_symmv_sort(gsl_vector * eval, gsl_matrix * evec, gsl_eigen_sort_t sort_type)

  int  gsl_eigen_hermv_sort(gsl_vector * eval, gsl_matrix_complex * evec, gsl_eigen_sort_t sort_type)

