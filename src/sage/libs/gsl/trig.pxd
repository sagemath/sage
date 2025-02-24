# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_sf_trig.h":

  double  gsl_sf_sin(double x)

  int  gsl_sf_sin_e(double x, gsl_sf_result * result)

  double  gsl_sf_cos(double x)

  int  gsl_sf_cos_e(double x, gsl_sf_result * result)

  double  gsl_sf_hypot(double x, double y)

  int  gsl_sf_hypot_e(double x, double y, gsl_sf_result * result)

  double  gsl_sf_sinc(double x)

  int  gsl_sf_sinc_e(double x, gsl_sf_result * result)

  int  gsl_sf_complex_sin_e(double zr, double zi, gsl_sf_result * szr, gsl_sf_result * szi)

  int  gsl_sf_complex_cos_e(double zr, double zi, gsl_sf_result * czr, gsl_sf_result * czi)

  int  gsl_sf_complex_logsin_e(double zr, double zi, gsl_sf_result * lszr, gsl_sf_result * lszi)

  double  gsl_sf_lnsinh(double x)

  int  gsl_sf_lnsinh_e(double x, gsl_sf_result * result)

  double  gsl_sf_lncosh(double x)

  int  gsl_sf_lncosh_e(double x, gsl_sf_result * result)

  int  gsl_sf_polar_to_rect(double r, double theta, gsl_sf_result * x, gsl_sf_result * y)

  int  gsl_sf_rect_to_polar(double x, double y, gsl_sf_result * r, gsl_sf_result * theta)

  double  gsl_sf_angle_restrict_symm(double theta)

  int  gsl_sf_angle_restrict_symm_e(double * theta)

  double  gsl_sf_angle_restrict_pos(double theta)

  int  gsl_sf_angle_restrict_pos_e(double * theta)

  double  gsl_sf_sin_err(double x, double dx)

  int  gsl_sf_sin_err_e(double x, double dx, gsl_sf_result * result)

  double  gsl_sf_cos_err(double x, double dx)

  int  gsl_sf_cos_err_e(double x, double dx, gsl_sf_result * result)

