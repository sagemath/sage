# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_sf_expint.h":

  double  gsl_sf_expint_E1(double x)

  int  gsl_sf_expint_E1_e(double x, gsl_sf_result * result)

  double  gsl_sf_expint_E2(double x)

  int  gsl_sf_expint_E2_e(double x, gsl_sf_result * result)

  double  gsl_sf_expint_Ei(double x)

  int  gsl_sf_expint_Ei_e(double x, gsl_sf_result * result)

  double  gsl_sf_Shi(double x)

  int  gsl_sf_Shi_e(double x, gsl_sf_result * result)

  double  gsl_sf_Chi(double x)

  int  gsl_sf_Chi_e(double x, gsl_sf_result * result)

  double  gsl_sf_expint_3(double x)

  int  gsl_sf_expint_3_e(double x, gsl_sf_result * result)

  double  gsl_sf_Si(double x)

  int  gsl_sf_Si_e(double x, gsl_sf_result * result)

  double  gsl_sf_Ci(double x)

  int  gsl_sf_Ci_e(double x, gsl_sf_result * result)

  double  gsl_sf_atanint(double x)

  int  gsl_sf_atanint_e(double x, gsl_sf_result * result)

