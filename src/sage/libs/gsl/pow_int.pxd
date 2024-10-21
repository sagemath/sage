# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_sf_pow_int.h":

  double  gsl_sf_pow_int(double x, int n)

  int  gsl_sf_pow_int_e(double x, int n, gsl_sf_result * result)

