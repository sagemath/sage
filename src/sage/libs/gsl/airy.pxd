# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from .types cimport *

cdef extern from "gsl/gsl_sf_airy.h":

  double  gsl_sf_airy_Ai(double x, gsl_mode_t mode)

  int  gsl_sf_airy_Ai_e(double x, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_airy_Bi(double x, gsl_mode_t mode)

  int  gsl_sf_airy_Bi_e(double x, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_airy_Ai_scaled(double x, gsl_mode_t mode)

  int  gsl_sf_airy_Ai_scaled_e(double x, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_airy_Bi_scaled(double x, gsl_mode_t mode)

  int  gsl_sf_airy_Bi_scaled_e(double x, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_airy_Ai_deriv(double x, gsl_mode_t mode)

  int  gsl_sf_airy_Ai_deriv_e(double x, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_airy_Bi_deriv(double x, gsl_mode_t mode)

  int  gsl_sf_airy_Bi_deriv_e(double x, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_airy_Ai_deriv_scaled(double x, gsl_mode_t mode)

  int  gsl_sf_airy_Ai_deriv_scaled_e(double x, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_airy_Bi_deriv_scaled(double x, gsl_mode_t mode)

  int  gsl_sf_airy_Bi_deriv_scaled_e(double x, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_airy_zero_Ai(unsigned int s)

  int  gsl_sf_airy_zero_Ai_e(unsigned int s, gsl_sf_result * result)

  double  gsl_sf_airy_zero_Bi(unsigned int s)

  int  gsl_sf_airy_zero_Bi_e(unsigned int s, gsl_sf_result * result)

  double  gsl_sf_airy_zero_Ai_deriv(unsigned int s)

  int  gsl_sf_airy_zero_Ai_deriv_e(unsigned int s, gsl_sf_result * result)

  double  gsl_sf_airy_zero_Bi_deriv(unsigned int s)

  int  gsl_sf_airy_zero_Bi_deriv_e(unsigned int s, gsl_sf_result * result)

