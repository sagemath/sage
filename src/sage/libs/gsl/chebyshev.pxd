# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_chebyshev.h":

  ctypedef struct gsl_cheb_series


  gsl_cheb_series * gsl_cheb_alloc( size_t order)

  void gsl_cheb_free(gsl_cheb_series * cs)

  int gsl_cheb_init(gsl_cheb_series * cs,  gsl_function * func,
                     double a,  double b)


  double gsl_cheb_eval( gsl_cheb_series * cs,  double x)
  int gsl_cheb_eval_err( gsl_cheb_series * cs,  double x,
                        double * result, double * abserr)


  double gsl_cheb_eval_n( gsl_cheb_series * cs,  size_t order, double x)
  int gsl_cheb_eval_n_err( gsl_cheb_series * cs,  size_t order,
                           double x, double * result, double * abserr)


  double gsl_cheb_eval_mode( gsl_cheb_series * cs, double x, gsl_mode_t mode)
  int gsl_cheb_eval_mode_e( gsl_cheb_series * cs,  double x, gsl_mode_t mode, double * result, double * abserr)


  int gsl_cheb_calc_deriv(gsl_cheb_series * deriv,  gsl_cheb_series * cs)

  int gsl_cheb_calc_integ(gsl_cheb_series * integ,  gsl_cheb_series * cs)

