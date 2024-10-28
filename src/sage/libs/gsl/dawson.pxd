# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_sf_dawson.h":

  double  gsl_sf_dawson(double x)

  int  gsl_sf_dawson_e(double x, gsl_sf_result * result)

