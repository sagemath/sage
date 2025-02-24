# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_combination.h":

  ctypedef struct gsl_combination:
    size_t n
    size_t k
    size_t *data

  gsl_combination *  gsl_combination_alloc(size_t n, size_t k)

  gsl_combination *  gsl_combination_calloc(size_t n, size_t k)

  void  gsl_combination_init_first(gsl_combination * c)

  void  gsl_combination_init_last(gsl_combination * c)

  void  gsl_combination_free(gsl_combination * c)

  int  gsl_combination_memcpy(gsl_combination * dest, gsl_combination * src)

  int  gsl_combination_fread(FILE * stream, gsl_combination * c)

  int  gsl_combination_fwrite(FILE * stream, gsl_combination * c)

  int  gsl_combination_fscanf(FILE * stream, gsl_combination * c)

  int  gsl_combination_fprintf(FILE * stream, gsl_combination * c, char *format)

  size_t  gsl_combination_n(gsl_combination * c)

  size_t  gsl_combination_k(gsl_combination * c)

  size_t *  gsl_combination_data(gsl_combination * c)

  size_t  gsl_combination_get(gsl_combination * c, size_t i)

  int  gsl_combination_valid(gsl_combination * c)

  int  gsl_combination_next(gsl_combination * c)

  int  gsl_combination_prev(gsl_combination * c)

