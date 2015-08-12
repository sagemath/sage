# distutils: libraries = GSL_LIBRARIES
from .types cimport *

cdef extern from "gsl/gsl_sf_pow_int.h":

  double  gsl_sf_pow_int(double x, int n)

  int  gsl_sf_pow_int_e(double x, int n, gsl_sf_result * result)

