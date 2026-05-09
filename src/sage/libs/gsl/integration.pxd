# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_integration.h":

  ctypedef struct gsl_integration_workspace
  ctypedef struct gsl_integration_qaws_table
  ctypedef struct  gsl_integration_qawo_table
  cdef enum:
    GSL_INTEG_GAUSS15 = 1
    GSL_INTEG_GAUSS21 = 2
    GSL_INTEG_GAUSS31 = 3
    GSL_INTEG_GAUSS41 = 4
    GSL_INTEG_GAUSS51 = 5
    GSL_INTEG_GAUSS61 = 6
  cdef enum gsl_integration_qawo_enum:
    GSL_INTEG_COSINE, GSL_INTEG_SINE

  int  gsl_integration_qng(gsl_function *f, double a, double b, double epsabs, double epsrel, double * result, double * abserr, size_t * neval)

  gsl_integration_workspace *  gsl_integration_workspace_alloc(size_t n)

  void  gsl_integration_workspace_free(gsl_integration_workspace * w)

  int  gsl_integration_qag(gsl_function *f, double a, double b, double epsabs, double epsrel, size_t limit, int key, gsl_integration_workspace * workspace, double * result, double * abserr)

  int  gsl_integration_qags(gsl_function * f, double a, double b, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr)

  int  gsl_integration_qagp(gsl_function * f, double *pts, size_t npts, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr)

  int  gsl_integration_qagi(gsl_function * f, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr)

  int  gsl_integration_qagiu(gsl_function * f, double a, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr)

  int  gsl_integration_qagil(gsl_function * f, double b, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr)

  int  gsl_integration_qawc(gsl_function *f, double a, double b, double c, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr)

  gsl_integration_qaws_table *  gsl_integration_qaws_table_alloc(double alpha, double beta, int mu, int nu)

  int  gsl_integration_qaws_table_set(gsl_integration_qaws_table * t, double alpha, double beta, int mu, int nu)

  void  gsl_integration_qaws_table_free(gsl_integration_qaws_table * t)

  int  gsl_integration_qaws(gsl_function * f, double a, double b, gsl_integration_qaws_table * t, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr)

  gsl_integration_qawo_table *  gsl_integration_qawo_table_alloc(double omega, double L,  gsl_integration_qawo_enum sine, size_t n)

  int  gsl_integration_qawo_table_set(gsl_integration_qawo_table * t, double omega, double L,  gsl_integration_qawo_enum sine)

  int  gsl_integration_qawo_table_set_length(gsl_integration_qawo_table * t, double L)

  void  gsl_integration_qawo_table_free(gsl_integration_qawo_table * t)

  int  gsl_integration_qawo(gsl_function * f, double a, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, gsl_integration_qawo_table * wf, double *result, double *abserr)

  int  gsl_integration_qawf(gsl_function * f, double a, double epsabs, size_t limit, gsl_integration_workspace * workspace, gsl_integration_workspace * cycle_workspace, gsl_integration_qawo_table * wf, double *result, double *abserr)

  double GSL_EMAXITER

  double GSL_EROUND

  double GSL_ESING

  double GSL_EDIVERGE

