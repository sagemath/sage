# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_sf_dilog.h":

  double  gsl_sf_dilog(double x)

  int  gsl_sf_dilog_e(double x, gsl_sf_result * result)

  int  gsl_sf_complex_dilog_e(double r, double theta, gsl_sf_result * result_re, gsl_sf_result * result_im)

