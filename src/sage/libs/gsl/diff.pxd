# distutils: libraries = GSL_LIBRARIES
from .types cimport *

cdef extern from "gsl/gsl_diff.h":
  int gsl_diff_central ( gsl_function *f, double x,
                        double *result, double *abserr)

  int gsl_diff_backward ( gsl_function *f, double x,
                         double *result, double *abserr)

  int gsl_diff_forward ( gsl_function *f, double x,
                        double *result, double *abserr)

