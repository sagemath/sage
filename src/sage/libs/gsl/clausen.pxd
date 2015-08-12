# distutils: libraries = GSL_LIBRARIES
from .types cimport *

cdef extern from "gsl/gsl_sf_clausen.h":

  double  gsl_sf_clausen(double x)

  int  gsl_sf_clausen_e(double x, gsl_sf_result * result)

