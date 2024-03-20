# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_min.h":

  ctypedef struct gsl_min_fminimizer_type:
      char *name
      size_t size
      int (*set) (void *state, gsl_function * f, double x_minimum, double f_minimum, double x_lower, double f_lower, double x_upper, double f_upper)
      int (*iterate) (void *state, gsl_function * f, double * x_minimum, double * f_minimum, double * x_lower, double * f_lower, double * x_upper, double * f_upper)

  ctypedef struct gsl_min_fminimizer:
      gsl_min_fminimizer_type * type
      gsl_function * function
      double x_minimum
      double x_lower
      double x_upper
      double f_minimum, f_lower, f_upper
      void *state

  gsl_min_fminimizer * gsl_min_fminimizer_alloc ( gsl_min_fminimizer_type * T)

  void gsl_min_fminimizer_free (gsl_min_fminimizer * s)

  int gsl_min_fminimizer_set (gsl_min_fminimizer * s,
                              gsl_function * f, double x_minimum,
                              double x_lower, double x_upper)

  int gsl_min_fminimizer_set_with_values (gsl_min_fminimizer * s,
                                          gsl_function * f,
                                          double x_minimum, double f_minimum,
                                          double x_lower, double f_lower,
                                          double x_upper, double f_upper)

  int gsl_min_fminimizer_iterate (gsl_min_fminimizer * s)

  char * gsl_min_fminimizer_name ( gsl_min_fminimizer * s)

  double gsl_min_fminimizer_x_minimum ( gsl_min_fminimizer * s)
  double gsl_min_fminimizer_x_lower ( gsl_min_fminimizer * s)
  double gsl_min_fminimizer_x_upper ( gsl_min_fminimizer * s)
  double gsl_min_fminimizer_f_minimum ( gsl_min_fminimizer * s)
  double gsl_min_fminimizer_f_lower ( gsl_min_fminimizer * s)
  double gsl_min_fminimizer_f_upper ( gsl_min_fminimizer * s)

  #Deprecated, use x_minimum instead
  #double gsl_min_fminimizer_minimum ( gsl_min_fminimizer * s)

  int gsl_min_test_interval (double x_lower, double x_upper, double epsabs, double epsrel)

  gsl_min_fminimizer_type  * gsl_min_fminimizer_goldensection
  gsl_min_fminimizer_type  * gsl_min_fminimizer_brent

  ctypedef int (*gsl_min_bracketing_function)(gsl_function *f,
                                     double *x_minimum,double * f_minimum,
                                     double *x_lower, double * f_lower,
                                     double *x_upper, double * f_upper,
                                     size_t eval_max)

  int gsl_min_find_bracket(gsl_function *f,double *x_minimum,double * f_minimum,
                       double *x_lower, double * f_lower,
                       double *x_upper, double * f_upper,
                       size_t eval_max)

