# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from .types cimport *

cdef extern from "gsl/gsl_sf_lambert.h":

  double  gsl_sf_lambert_W0(double x)

  int  gsl_sf_lambert_W0_e(double x, gsl_sf_result * result)

  double  gsl_sf_lambert_Wm1(double x)

  int  gsl_sf_lambert_Wm1_e(double x, gsl_sf_result * result)

