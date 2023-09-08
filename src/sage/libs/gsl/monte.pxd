# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_monte.h":
  ctypedef struct gsl_monte_function:
    double (*f)(double * x_array, size_t dim, void * params)
    size_t dim
    void * params

cdef extern from "gsl/gsl_monte_plain.h":
  ctypedef struct gsl_monte_plain_state

  int gsl_monte_plain_integrate ( gsl_monte_function * f, \
                              double xl[],  double xu[], \
                              size_t dim, size_t calls, gsl_rng * r,  \
                             gsl_monte_plain_state * state, \
                             double *result, double *abserr)

  gsl_monte_plain_state* gsl_monte_plain_alloc(size_t dim)

  int gsl_monte_plain_init(gsl_monte_plain_state* state)

  void gsl_monte_plain_free (gsl_monte_plain_state* state)

cdef extern from "gsl/gsl_monte_miser.h":

  ctypedef struct gsl_monte_miser_state

  int gsl_monte_miser_integrate(gsl_monte_function * f, double xl[],  double xh[],  size_t dim, size_t calls, gsl_rng *r,  gsl_monte_miser_state* state, double *result, double *abserr)

  gsl_monte_miser_state* gsl_monte_miser_alloc(size_t dim)

  int gsl_monte_miser_init(gsl_monte_miser_state* state)

  void gsl_monte_miser_free(gsl_monte_miser_state* state)

  # vegas
cdef extern from "gsl/gsl_monte_vegas.h":
  cdef enum:
    GSL_VEGAS_MODE_IMPORTANCE = 1
    GSL_VEGAS_MODE_IMPORTANCE_ONLY = 0
    GSL_VEGAS_MODE_STRATIFIED = -1

  ctypedef struct gsl_monte_vegas_state:
    double chisq

  int gsl_monte_vegas_integrate(gsl_monte_function * f, double xl[], double xu[], size_t dim, size_t calls, gsl_rng * r, gsl_monte_vegas_state *state, double* result, double* abserr)

  gsl_monte_vegas_state* gsl_monte_vegas_alloc(size_t dim)

  int gsl_monte_vegas_init(gsl_monte_vegas_state* state)

  void gsl_monte_vegas_free (gsl_monte_vegas_state* state)
