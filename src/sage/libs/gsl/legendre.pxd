# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_sf_legendre.h":

  double  gsl_sf_legendre_P1(double x)

  double  gsl_sf_legendre_P2(double x)

  double  gsl_sf_legendre_P3(double x)

  int  gsl_sf_legendre_P1_e(double x, gsl_sf_result * result)

  int  gsl_sf_legendre_P2_e(double x, gsl_sf_result * result)

  int  gsl_sf_legendre_P3_e(double x, gsl_sf_result * result)

  double  gsl_sf_legendre_Pl(int l, double x)

  int  gsl_sf_legendre_Pl_e(int l, double x, gsl_sf_result * result)

  int  gsl_sf_legendre_Pl_array(int lmax, double x, double result_array[])

  double  gsl_sf_legendre_Q0(double x)

  int  gsl_sf_legendre_Q0_e(double x, gsl_sf_result * result)

  double  gsl_sf_legendre_Q1(double x)

  int  gsl_sf_legendre_Q1_e(double x, gsl_sf_result * result)

  double  gsl_sf_legendre_Ql(int l, double x)

  int  gsl_sf_legendre_Ql_e(int l, double x, gsl_sf_result * result)

  double  gsl_sf_legendre_Plm(int l, int m, double x)

  int  gsl_sf_legendre_Plm_e(int l, int m, double x, gsl_sf_result * result)

  int  gsl_sf_legendre_Plm_array(int lmax, int m, double x, double result_array[])

  double  gsl_sf_legendre_sphPlm(int l, int m, double x)

  int  gsl_sf_legendre_sphPlm_e(int l, int m, double x, gsl_sf_result * result)

  int  gsl_sf_legendre_sphPlm_array(int lmax, int m, double x, double result_array[])

  int  gsl_sf_legendre_array_size(int lmax, int m)

  double  gsl_sf_conicalP_half(double lambd, double x)

  int  gsl_sf_conicalP_half_e(double lambd, double x, gsl_sf_result * result)

  double  gsl_sf_conicalP_mhalf(double lambd, double x)

  int  gsl_sf_conicalP_mhalf_e(double lambd, double x, gsl_sf_result * result)

  double  gsl_sf_conicalP_0(double lambd, double x)

  int  gsl_sf_conicalP_0_e(double lambd, double x, gsl_sf_result * result)

  double  gsl_sf_conicalP_1(double lambd, double x)

  int  gsl_sf_conicalP_1_e(double lambd, double x, gsl_sf_result * result)

  double  gsl_sf_conicalP_sph_reg(int l, double lambd, double x)

  int  gsl_sf_conicalP_sph_reg_e(int l, double lambd, double x, gsl_sf_result * result)

  double  gsl_sf_conicalP_cyl_reg(int m, double lambd, double x)

  int  gsl_sf_conicalP_cyl_reg_e(int m, double lambd, double x, gsl_sf_result * result)

  double  gsl_sf_legendre_H3d_0(double lambd, double eta)

  int  gsl_sf_legendre_H3d_0_e(double lambd, double eta, gsl_sf_result * result)

  double  gsl_sf_legendre_H3d_1(double lambd, double eta)

  int  gsl_sf_legendre_H3d_1_e(double lambd, double eta, gsl_sf_result * result)

  double  gsl_sf_legendre_H3d(int l, double lambd, double eta)

  int  gsl_sf_legendre_H3d_e(int l, double lambd, double eta, gsl_sf_result * result)

  int  gsl_sf_legendre_H3d_array(int lmax, double lambd, double eta, double result_array[])

