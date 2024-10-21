# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_sf_exp.h":

  double  gsl_sf_exp(double x)

  int  gsl_sf_exp_e(double x, gsl_sf_result * result)

  int  gsl_sf_exp_e10_e(double x, gsl_sf_result_e10 * result)

  double  gsl_sf_exp_mult(double x, double y)

  int  gsl_sf_exp_mult_e(double x, double y, gsl_sf_result * result)

  int  gsl_sf_exp_mult_e10_e(double x, double y, gsl_sf_result_e10 * result)

  double  gsl_sf_expm1(double x)

  int  gsl_sf_expm1_e(double x, gsl_sf_result * result)

  double  gsl_sf_exprel(double x)

  int  gsl_sf_exprel_e(double x, gsl_sf_result * result)

  double  gsl_sf_exprel_2(double x)

  int  gsl_sf_exprel_2_e(double x, gsl_sf_result * result)

  double  gsl_sf_exprel_n(int n, double x)

  int  gsl_sf_exprel_n_e(int n, double x, gsl_sf_result * result)

  int  gsl_sf_exp_err_e(double x, double dx, gsl_sf_result * result)

  int  gsl_sf_exp_err_e10_e(double x, double dx, gsl_sf_result_e10 * result)

  int  gsl_sf_exp_mult_err_e(double x, double dx, double y, double dy, gsl_sf_result * result)

  int  gsl_sf_exp_mult_err_e10_e(double x, double dx, double y, double dy, gsl_sf_result_e10 * result)

