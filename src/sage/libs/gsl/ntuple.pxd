# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from .types cimport *

cdef extern from "gsl/gsl_ntuple.h":
  ctypedef struct gsl_ntuple

  ctypedef struct gsl_ntuple_select_fn:
    int (* function) (void * ntuple_data, void * params)
    void * params

  ctypedef struct gsl_ntuple_value_fn:
    double (* function) (void * ntuple_data, void * params)
    void * params

  gsl_ntuple * gsl_ntuple_open (char * filename, void * ntuple_data, size_t size)

  gsl_ntuple * gsl_ntuple_create (char * filename, void * ntuple_data, size_t size)

  int gsl_ntuple_write (gsl_ntuple * ntuple)
  int gsl_ntuple_read (gsl_ntuple * ntuple)

  int gsl_ntuple_bookdata (gsl_ntuple * ntuple)

  int gsl_ntuple_project (gsl_histogram * h, gsl_ntuple * ntuple,
                          gsl_ntuple_value_fn *value_func,
                          gsl_ntuple_select_fn *select_func)

  int gsl_ntuple_close (gsl_ntuple * ntuple)

