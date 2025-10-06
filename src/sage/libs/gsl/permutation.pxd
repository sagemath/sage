# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_permutation.h":
  # Allocation
  gsl_permutation *  gsl_permutation_alloc(size_t n)

  gsl_permutation *  gsl_permutation_calloc(size_t n)

  void  gsl_permutation_init(gsl_permutation * p)

  void  gsl_permutation_free(gsl_permutation * p)

  # Reading and writing permutations
  int  gsl_permutation_fread(FILE * stream, gsl_permutation * p)

  int  gsl_permutation_fwrite(FILE * stream, gsl_permutation * p)

  int  gsl_permutation_fscanf(FILE * stream, gsl_permutation * p)

  int  gsl_permutation_fprintf(FILE * stream, gsl_permutation * p, char *format)

  # Permutation properties
  size_t  gsl_permutation_size(gsl_permutation * p)

  size_t *  gsl_permutation_data(gsl_permutation * p)

  int  gsl_permutation_memcpy(gsl_permutation * dest, gsl_permutation * src)

  size_t  gsl_permutation_get(gsl_permutation * p, size_t i)

  int  gsl_permutation_swap(gsl_permutation * p, size_t i, size_t j)

  int  gsl_permutation_valid(gsl_permutation * p)

  # Permutation functions
  void  gsl_permutation_reverse(gsl_permutation * p)

  int  gsl_permutation_inverse(gsl_permutation * inv, gsl_permutation * p)

  int  gsl_permutation_next(gsl_permutation * p)

  int  gsl_permutation_prev(gsl_permutation * p)

  int  gsl_permutation_mul(gsl_permutation * p, gsl_permutation * pa, gsl_permutation * pb)

  # Permutations in Cyclic Form
  int  gsl_permutation_linear_to_canonical(gsl_permutation * q, gsl_permutation * p)

  int  gsl_permutation_canonical_to_linear(gsl_permutation * p, gsl_permutation * q)

  size_t  gsl_permutation_inversions(gsl_permutation * p)

  size_t  gsl_permutation_linear_cycles(gsl_permutation * p)

  size_t  gsl_permutation_canonical_cycles(gsl_permutation * q)


# Applying Permutations
cdef extern from "gsl/gsl_permute_double.h":
  int  gsl_permute(size_t * p, double * data, size_t stride, size_t n)

  int  gsl_permute_inverse(size_t * p, double * data, size_t stride, size_t n)


cdef extern from "gsl/gsl_permute_vector_double.h":
  int  gsl_permute_vector(gsl_permutation * p, gsl_vector * v)

  int  gsl_permute_vector_inverse(gsl_permutation * p, gsl_vector * v)

cdef extern from "gsl/gsl_permute_vector_complex_double.h":
  int gsl_permute_vector_complex ( gsl_permutation * p, gsl_vector_complex * v)

  int gsl_permute_vector_complex_inverse ( gsl_permutation * p, gsl_vector_complex * v)

