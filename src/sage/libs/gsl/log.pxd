# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_sf_log.h":

  double  gsl_sf_log(double x)

  int  gsl_sf_log_e(double x, gsl_sf_result * result)

  double  gsl_sf_log_abs(double x)

  int  gsl_sf_log_abs_e(double x, gsl_sf_result * result)

  int  gsl_sf_complex_log_e(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * theta)

  double  gsl_sf_log_1plusx(double x)

  int  gsl_sf_log_1plusx_e(double x, gsl_sf_result * result)

  double  gsl_sf_log_1plusx_mx(double x)

  int  gsl_sf_log_1plusx_mx_e(double x, gsl_sf_result * result)

