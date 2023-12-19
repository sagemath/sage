# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_roots.h":

  ctypedef struct gsl_root_fsolver_type:
      char *name
      size_t size
      int (*set) (void *state, gsl_function * f, double * root, double x_lower, double x_upper)
      int (*iterate) (void *state, gsl_function * f, double * root, double * x_lower, double * x_upper)

  ctypedef struct gsl_root_fsolver:
      gsl_root_fsolver_type * type
      gsl_function * function
      double root
      double x_lower
      double x_upper
      void *state

  ctypedef struct gsl_root_fdfsolver_type:
      char *name
      size_t size
      int (*set) (void *state, gsl_function_fdf * f, double * root)
      int (*iterate) (void *state, gsl_function_fdf * f, double * root)

  ctypedef struct gsl_root_fdfsolver:
      gsl_root_fdfsolver_type * type
      gsl_function_fdf * fdf
      double root
      void *state

  gsl_root_fsolver * gsl_root_fsolver_alloc ( gsl_root_fsolver_type * T)
  void gsl_root_fsolver_free (gsl_root_fsolver * s)

  int gsl_root_fsolver_set (gsl_root_fsolver * s, gsl_function * f,
                            double x_lower, double x_upper)

  int gsl_root_fsolver_iterate (gsl_root_fsolver * s)

  char * gsl_root_fsolver_name ( gsl_root_fsolver * s)
  double gsl_root_fsolver_root ( gsl_root_fsolver * s)
  double gsl_root_fsolver_x_lower ( gsl_root_fsolver * s)
  double gsl_root_fsolver_x_upper ( gsl_root_fsolver * s)


  gsl_root_fdfsolver * gsl_root_fdfsolver_alloc ( gsl_root_fdfsolver_type * T)

  int gsl_root_fdfsolver_set (gsl_root_fdfsolver * s,
                           gsl_function_fdf * fdf, double root)

  int gsl_root_fdfsolver_iterate (gsl_root_fdfsolver * s)

  void gsl_root_fdfsolver_free (gsl_root_fdfsolver * s)

  char * gsl_root_fdfsolver_name ( gsl_root_fdfsolver * s)
  double gsl_root_fdfsolver_root ( gsl_root_fdfsolver * s)

  int gsl_root_test_interval (double x_lower, double x_upper, double epsabs, double epsrel)

  int gsl_root_test_residual (double f, double epsabs)

  int gsl_root_test_delta (double x1, double x0, double epsabs, double epsrel)

  gsl_root_fsolver_type  * gsl_root_fsolver_bisection
  gsl_root_fsolver_type  * gsl_root_fsolver_brent
  gsl_root_fsolver_type  * gsl_root_fsolver_falsepos
  gsl_root_fdfsolver_type  * gsl_root_fdfsolver_newton
  gsl_root_fdfsolver_type  * gsl_root_fdfsolver_secant
  gsl_root_fdfsolver_type  * gsl_root_fdfsolver_steffenson
