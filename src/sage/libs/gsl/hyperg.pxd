# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from .types cimport *

cdef extern from "gsl/gsl_sf_hyperg.h":

  double  gsl_sf_hyperg_0F1(double c, double x)

  int  gsl_sf_hyperg_0F1_e(double c, double x, gsl_sf_result * result)

  double  gsl_sf_hyperg_1F1_int(int m, int n, double x)

  int  gsl_sf_hyperg_1F1_int_e(int m, int n, double x, gsl_sf_result * result)

  double  gsl_sf_hyperg_1F1(double a, double b, double x)

  int  gsl_sf_hyperg_1F1_e(double a, double b, double x, gsl_sf_result * result)

  double  gsl_sf_hyperg_U_int(int m, int n, double x)

  int  gsl_sf_hyperg_U_int_e(int m, int n, double x, gsl_sf_result * result)

  int  gsl_sf_hyperg_U_int_e10_e(int m, int n, double x, gsl_sf_result_e10 * result)

  double  gsl_sf_hyperg_U(double a, double b, double x)

  int  gsl_sf_hyperg_U_e(double a, double b, double x)

  int  gsl_sf_hyperg_U_e10_e(double a, double b, double x, gsl_sf_result_e10 * result)

  double  gsl_sf_hyperg_2F1(double a, double b, double c, double x)

  int  gsl_sf_hyperg_2F1_e(double a, double b, double c, double x, gsl_sf_result * result)

  double  gsl_sf_hyperg_2F1_conj(double aR, double aI, double c, double x)

  int  gsl_sf_hyperg_2F1_conj_e(double aR, double aI, double c, double x, gsl_sf_result * result)

  double  gsl_sf_hyperg_2F1_renorm(double a, double b, double c, double x)

  int  gsl_sf_hyperg_2F1_renorm_e(double a, double b, double c, double x, gsl_sf_result * result)

  double  gsl_sf_hyperg_2F1_conj_renorm(double aR, double aI, double c, double x)

  int  gsl_sf_hyperg_2F1_conj_renorm_e(double aR, double aI, double c, double x, gsl_sf_result * result)

  double  gsl_sf_hyperg_2F0(double a, double b, double x)

  int  gsl_sf_hyperg_2F0_e(double a, double b, double x, gsl_sf_result * result)

