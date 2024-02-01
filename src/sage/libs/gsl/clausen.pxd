# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_sf_clausen.h":

  double  gsl_sf_clausen(double x)

  int  gsl_sf_clausen_e(double x, gsl_sf_result * result)

