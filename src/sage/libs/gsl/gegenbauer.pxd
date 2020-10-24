# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from .types cimport *

cdef extern from "gsl/gsl_sf_gegenbauer.h":

  double  gsl_sf_gegenpoly_1(double lambd, double x)

  double  gsl_sf_gegenpoly_2(double lambd, double x)

  double  gsl_sf_gegenpoly_3(double lambd, double x)

  int  gsl_sf_gegenpoly_1_e(double lambd, double x, gsl_sf_result * result)

  int  gsl_sf_gegenpoly_2_e(double lambd, double x, gsl_sf_result * result)

  int  gsl_sf_gegenpoly_3_e(double lambd, double x, gsl_sf_result * result)

  double  gsl_sf_gegenpoly_n(int n, double lambd, double x)

  int  gsl_sf_gegenpoly_n_e(int n, double lambd, double x, gsl_sf_result * result)

  int  gsl_sf_gegenpoly_array(int nmax, double lambd, double x, double result_array[])

