# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_sf_ellint.h":

  double  gsl_sf_ellint_Kcomp(double k, gsl_mode_t mode)

  int  gsl_sf_ellint_Kcomp_e(double k, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_ellint_Ecomp(double k, gsl_mode_t mode)

  int  gsl_sf_ellint_Ecomp_e(double k, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_ellint_F(double phi, double k, gsl_mode_t mode)

  int  gsl_sf_ellint_F_e(double phi, double k, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_ellint_E(double phi, double k, gsl_mode_t mode)

  int  gsl_sf_ellint_E_e(double phi, double k, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_ellint_P(double phi, double k, double n, gsl_mode_t mode)

  int  gsl_sf_ellint_P_e(double phi, double k, double n, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_ellint_D(double phi, double k, double n, gsl_mode_t mode)

  int  gsl_sf_ellint_D_e(double phi, double k, double n, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_ellint_RC(double x, double y, gsl_mode_t mode)

  int  gsl_sf_ellint_RC_e(double x, double y, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_ellint_RD(double x, double y, double z, gsl_mode_t mode)

  int  gsl_sf_ellint_RD_e(double x, double y, double z, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_ellint_RF(double x, double y, double z, gsl_mode_t mode)

  int  gsl_sf_ellint_RF_e(double x, double y, double z, gsl_mode_t mode, gsl_sf_result * result)

  double  gsl_sf_ellint_RJ(double x, double y, double z, double p, gsl_mode_t mode)

  int  gsl_sf_ellint_RJ_e(double x, double y, double z, double p, gsl_mode_t mode, gsl_sf_result * result)

