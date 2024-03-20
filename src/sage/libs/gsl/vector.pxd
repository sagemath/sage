# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_vector.h":
  # Allocation
  gsl_vector *  gsl_vector_alloc(size_t n)

  gsl_vector *  gsl_vector_calloc(size_t n)

  gsl_vector_alloc_from_block(gsl_block * b, size_t offset,
                              size_t n, size_t stride)

  gsl_vector *gsl_vector_alloc_from_vector(gsl_vector * v,
                         size_t offset, size_t n, size_t stride)

  void  gsl_vector_free(gsl_vector * v)

  # Views
  gsl_vector_view  gsl_vector_view_array(double *base, size_t n)

  gsl_vector_view  gsl_vector_subvector(gsl_vector *v, size_t offset, size_t n)

  gsl_vector_view  gsl_vector_view_array_with_stride(double * base, size_t stride, size_t n)

  gsl_vector_const_view  gsl_vector_const_view_array(double *base, size_t n)

  gsl_vector_const_view  gsl_vector_const_view_array_with_stride(double * base, size_t stride, size_t n)

  gsl_vector_const_view  gsl_vector_const_subvector(gsl_vector * v, size_t offset, size_t n)

  gsl_vector_view  gsl_vector_subvector_with_stride(gsl_vector *v, size_t offset, size_t stride, size_t n)

  gsl_vector_const_view  gsl_vector_const_subvector_with_stride(gsl_vector * v, size_t offset, size_t stride, size_t n)


  # Operations
  double  gsl_vector_get(gsl_vector * v, size_t i)

  void  gsl_vector_set(gsl_vector * v, size_t i, double x)

  double *  gsl_vector_ptr(gsl_vector * v, size_t i)

  double *  gsl_vector_const_ptr(gsl_vector * v, size_t i)

  void  gsl_vector_set_zero(gsl_vector * v)

  void  gsl_vector_set_all(gsl_vector * v, double x)

  int  gsl_vector_set_basis(gsl_vector * v, size_t i)

  # Reading and writing vectors
  int  gsl_vector_fread(FILE * stream, gsl_vector * v)

  int  gsl_vector_fwrite(FILE * stream, gsl_vector * v)

  int  gsl_vector_fscanf(FILE * stream, gsl_vector * v)

  int  gsl_vector_fprintf(FILE * stream, gsl_vector * v, char * format)

  # Copying or exchanging elements
  int  gsl_vector_memcpy(gsl_vector * dest, gsl_vector * src)

  int  gsl_vector_reverse(gsl_vector * v)

  int  gsl_vector_swap(gsl_vector * v, gsl_vector * w)

  int  gsl_vector_swap_elements(gsl_vector * v, size_t i, size_t j)

  # Finding maximum and minimum elements of vectors

  double  gsl_vector_max(gsl_vector * v)

  double  gsl_vector_min(gsl_vector * v)

  void  gsl_vector_minmax(gsl_vector * v, double * min_out, double * max_out)

  size_t  gsl_vector_max_index(gsl_vector * v)

  size_t  gsl_vector_min_index(gsl_vector * v)

  void  gsl_vector_minmax_index(gsl_vector * v, size_t * imin, size_t * imax)

  # Vector operations
  int  gsl_vector_add(gsl_vector * a, gsl_vector * b)

  int  gsl_vector_sub(gsl_vector * a, gsl_vector * b)

  int  gsl_vector_mul(gsl_vector * a, gsl_vector * b)

  int  gsl_vector_div(gsl_vector * a, gsl_vector * b)

  int  gsl_vector_scale(gsl_vector * a, double x)

  int  gsl_vector_add_constant(gsl_vector * a, double x)

  int  gsl_vector_isnull(gsl_vector * v)


