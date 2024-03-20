# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_block_double.h":

  ctypedef struct gsl_block:
    size_t size
    double * data

  gsl_block *  gsl_block_alloc(size_t n)

  gsl_block *  gsl_block_calloc(size_t n)

  void  gsl_block_free(gsl_block * b)

  int  gsl_block_fread(FILE * stream, gsl_block * b)

  int  gsl_block_fwrite(FILE * stream, gsl_block * b)

  int  gsl_block_fscanf(FILE * stream, gsl_block * b)

  int  gsl_block_fprintf(FILE * stream, gsl_block * b, char * format)

  size_t gsl_block_size (gsl_block * b)
  double * gsl_block_data (gsl_block * b)


cdef extern from "gsl/gsl_block_complex_double.h":

  ctypedef struct gsl_block_complex:
    size_t size
    double * data

  gsl_block_complex *  gsl_block_complex_alloc(size_t n)

  gsl_block_complex *  gsl_block_complex_calloc(size_t n)

  void  gsl_block_complex_free(gsl_block_complex * b)

  int  gsl_block_complex_fread(FILE * stream, gsl_block_complex * b)

  int  gsl_block_complex_fwrite(FILE * stream, gsl_block_complex * b)

  int  gsl_block_complex_fscanf(FILE * stream, gsl_block_complex * b)

  int  gsl_block_complex_fprintf(FILE * stream, gsl_block_complex * b, char * format)

  size_t gsl_block_complex_size (gsl_block_complex * b)
  double * gsl_block_complex_data (gsl_block_complex * b)
