# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_sf_debye.h":

  double  gsl_sf_debye_1(double x)

  int  gsl_sf_debye_1_e(double x, gsl_sf_result * result)

  double  gsl_sf_debye_2(double x)

  int  gsl_sf_debye_2_e(double x, gsl_sf_result * result)

  double  gsl_sf_debye_3(double x)

  int  gsl_sf_debye_3_e(double x, gsl_sf_result * result)

  double  gsl_sf_debye_4(double x)

  int  gsl_sf_debye_4_e(double x, gsl_sf_result * result)

