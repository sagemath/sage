# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_sf_coulomb.h":

  double  gsl_sf_hydrogenicR_1(double Z, double r)

  int  gsl_sf_hydrogenicR_1_e(double Z, double r, gsl_sf_result * result)

  double  gsl_sf_hydrogenicR(int n, int l, double Z, double r)

  int  gsl_sf_hydrogenicR_e(int n, int l, double Z, double r, gsl_sf_result * result)

  int  gsl_sf_coulomb_wave_FG_e(double eta, double x, double L_F, int k, gsl_sf_result * F, gsl_sf_result * Fp, gsl_sf_result * G, gsl_sf_result * Gp, double * exp_F, double * exp_G)

  int  gsl_sf_coulomb_wave_F_array(double L_min, int kmax, double eta, double x, double fc_array[], double * F_exponent)

  int  gsl_sf_coulomb_wave_FG_array(double L_min, int kmax, double eta, double x, double fc_array[], double gc_array[], double * F_exponent, double * G_exponent)

  int  gsl_sf_coulomb_wave_FGp_array(double L_min, int kmax, double eta, double x, double fc_array[], double fcp_array[], double gc_array[], double gcp_array[], double * F_exponent, double * G_exponent)

  int  gsl_sf_coulomb_wave_sphF_array(double L_min, int kmax, double eta, double x, double fc_array[], double F_exponent[])

  int  gsl_sf_coulomb_CL_e(double L, double eta, gsl_sf_result * result)

  int  gsl_sf_coulomb_CL_array(double Lmin, int kmax, double eta, double cl[])

