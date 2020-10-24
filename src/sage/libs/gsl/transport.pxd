# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from .types cimport *

cdef extern from "gsl/gsl_sf_transport.h":

  double  gsl_sf_transport_2(double x)

  int  gsl_sf_transport_2_e(double x, gsl_sf_result * result)

  double  gsl_sf_transport_3(double x)

  int  gsl_sf_transport_3_e(double x, gsl_sf_result * result)

  double  gsl_sf_transport_4(double x)

  int  gsl_sf_transport_4_e(double x, gsl_sf_result * result)

  double  gsl_sf_transport_5(double x)

  int  gsl_sf_transport_5_e(double x, gsl_sf_result * result)

