# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_sf_laguerre.h":

  double  gsl_sf_laguerre_1(double a, double x)

  double  gsl_sf_laguerre_2(double a, double x)

  double  gsl_sf_laguerre_3(double a, double x)

  int  gsl_sf_laguerre_1_e(double a, double x, gsl_sf_result * result)

  int  gsl_sf_laguerre_2_e(double a, double x, gsl_sf_result * result)

  int  gsl_sf_laguerre_3_e(double a, double x, gsl_sf_result * result)

  double  gsl_sf_laguerre_n(int n, double a, double x)

  int  gsl_sf_laguerre_n_e(int n, double a, double x, gsl_sf_result * result)

