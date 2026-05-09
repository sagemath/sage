# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_sf_bessel.h":

  double  gsl_sf_bessel_J0(double x)

  int  gsl_sf_bessel_J0_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_J1(double x)

  int  gsl_sf_bessel_J1_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_Jn(int n, double x)

  int  gsl_sf_bessel_Jn_e(int n, double x, gsl_sf_result * result)

  int  gsl_sf_bessel_Jn_array(int nmin, int nmax, double x, double result_array[])

  double  gsl_sf_bessel_Y0(double x)

  int  gsl_sf_bessel_Y0_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_Y1(double x)

  int  gsl_sf_bessel_Y1_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_Yn(int n, double x)

  int  gsl_sf_bessel_Yn_e(int n, double x, gsl_sf_result * result)

  int  gsl_sf_bessel_Yn_array(int nmin, int nmax, double x, double result_array[])

  double  gsl_sf_bessel_I0(double x)

  int  gsl_sf_bessel_I0_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_I1(double x)

  int  gsl_sf_bessel_I1_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_In(int n, double x)

  int  gsl_sf_bessel_In_e(int n, double x, gsl_sf_result * result)

  int  gsl_sf_bessel_In_array(int nmin, int nmax, double x, double result_array[])

  double  gsl_sf_bessel_I0_scaled(double x)

  int  gsl_sf_bessel_I0_scaled_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_I1_scaled(double x)

  int  gsl_sf_bessel_I1_scaled_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_In_scaled(int n, double x)

  int  gsl_sf_bessel_In_scaled_e(int n, double x, gsl_sf_result * result)

  int  gsl_sf_bessel_In_scaled_array(int nmin, int nmax, double x, double result_array[])

  double  gsl_sf_bessel_K0(double x)

  int  gsl_sf_bessel_K0_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_K1(double x)

  int  gsl_sf_bessel_K1_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_Kn(int n, double x)

  int  gsl_sf_bessel_Kn_e(int n, double x, gsl_sf_result * result)

  int  gsl_sf_bessel_Kn_array(int nmin, int nmax, double x, double result_array[])

  double  gsl_sf_bessel_K0_scaled(double x)

  int  gsl_sf_bessel_K0_scaled_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_K1_scaled(double x)

  int  gsl_sf_bessel_K1_scaled_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_Kn_scaled(int n, double x)

  int  gsl_sf_bessel_Kn_scaled_e(int n, double x, gsl_sf_result * result)

  int  gsl_sf_bessel_Kn_scaled_array(int nmin, int nmax, double x, double result_array[])

  double  gsl_sf_bessel_j0(double x)

  int  gsl_sf_bessel_j0_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_j1(double x)

  int  gsl_sf_bessel_j1_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_j2(double x)

  int  gsl_sf_bessel_j2_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_jl(int l, double x)

  int  gsl_sf_bessel_jl_e(int l, double x, gsl_sf_result * result)

  int  gsl_sf_bessel_jl_array(int lmax, double x, double result_array[])

  int  gsl_sf_bessel_jl_steed_array(int lmax, double x, double * jl_x_array)

  double  gsl_sf_bessel_y0(double x)

  int  gsl_sf_bessel_y0_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_y1(double x)

  int  gsl_sf_bessel_y1_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_y2(double x)

  int  gsl_sf_bessel_y2_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_yl(int l, double x)

  int  gsl_sf_bessel_yl_e(int l, double x, gsl_sf_result * result)

  int  gsl_sf_bessel_yl_array(int lmax, double x, double result_array[])

  double  gsl_sf_bessel_i0_scaled(double x)

  int  gsl_sf_bessel_i0_scaled_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_i1_scaled(double x)

  int  gsl_sf_bessel_i1_scaled_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_i2_scaled(double x)

  int  gsl_sf_bessel_i2_scaled_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_il_scaled(int l, double x)

  int  gsl_sf_bessel_il_scaled_e(int l, double x, gsl_sf_result * result)

  int  gsl_sf_bessel_il_scaled_array(int lmax, double x, double result_array[])

  double  gsl_sf_bessel_k0_scaled(double x)

  int  gsl_sf_bessel_k0_scaled_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_k1_scaled(double x)

  int  gsl_sf_bessel_k1_scaled_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_k2_scaled(double x)

  int  gsl_sf_bessel_k2_scaled_e(double x, gsl_sf_result * result)

  double  gsl_sf_bessel_kl_scaled(int l, double x)

  int  gsl_sf_bessel_kl_scaled_e(int l, double x, gsl_sf_result * result)

  int  gsl_sf_bessel_kl_scaled_array(int lmax, double x, double result_array[])

  double  gsl_sf_bessel_Jnu(double nu, double x)

  int  gsl_sf_bessel_Jnu_e(double nu, double x, gsl_sf_result * result)

  int  gsl_sf_bessel_sequence_Jnu_e(double nu, gsl_mode_t mode, size_t size, double v[])

  double  gsl_sf_bessel_Ynu(double nu, double x)

  int  gsl_sf_bessel_Ynu_e(double nu, double x, gsl_sf_result * result)

  double  gsl_sf_bessel_Inu(double nu, double x)

  int  gsl_sf_bessel_Inu_e(double nu, double x, gsl_sf_result * result)

  double  gsl_sf_bessel_Inu_scaled(double nu, double x)

  int  gsl_sf_bessel_Inu_scaled_e(double nu, double x, gsl_sf_result * result)

  double  gsl_sf_bessel_Knu(double nu, double x)

  int  gsl_sf_bessel_Knu_e(double nu, double x, gsl_sf_result * result)

  double  gsl_sf_bessel_lnKnu(double nu, double x)

  int  gsl_sf_bessel_lnKnu_e(double nu, double x, gsl_sf_result * result)

  double  gsl_sf_bessel_Knu_scaled(double nu, double x)

  int  gsl_sf_bessel_Knu_scaled_e(double nu, double x, gsl_sf_result * result)

  double  gsl_sf_bessel_zero_J0(unsigned int s)

  int  gsl_sf_bessel_zero_J0_e(unsigned int s, gsl_sf_result * result)

  double  gsl_sf_bessel_zero_J1(unsigned int s)

  int  gsl_sf_bessel_zero_J1_e(unsigned int s, gsl_sf_result * result)

  double  gsl_sf_bessel_zero_Jnu(double nu, unsigned int s)

  int  gsl_sf_bessel_zero_Jnu_e(double nu, unsigned int s, gsl_sf_result * result)

