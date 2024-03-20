# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_vector_complex_double.h":
  # Allocation
  gsl_vector_complex *  gsl_vector_complex_alloc(size_t n)

  gsl_vector_complex *  gsl_vector_complex_calloc(size_t n)

  gsl_vector_complex_alloc_from_block(gsl_block_complex * b, size_t offset,
                              size_t n, size_t stride)

  gsl_vector_complex *gsl_vector_complex_alloc_from_vector(gsl_vector_complex * v,
                         size_t offset, size_t n, size_t stride)

  void  gsl_vector_complex_free(gsl_vector_complex * v)

  # Views
  gsl_vector_complex_view  gsl_vector_complex_view_array(double *base, size_t n)

  gsl_vector_complex_view  gsl_vector_complex_view_array_with_stride(double * base, size_t stride, size_t n)

  gsl_vector_complex_const_view  gsl_vector_complex_const_view_array(double *base, size_t n)

  gsl_vector_complex_const_view  gsl_vector_complex_const_view_array_with_stride(double * base, size_t stride, size_t n)

  gsl_vector_complex_view  gsl_vector_complex_subvector(gsl_vector_complex *v, size_t offset, size_t n)

  gsl_vector_complex_const_view  gsl_vector_complex_const_subvector(gsl_vector_complex * v, size_t offset, size_t n)

  gsl_vector_complex_view  gsl_vector_complex_subvector_with_stride(gsl_vector_complex *v, size_t offset, size_t stride, size_t n)

  gsl_vector_complex_const_view  gsl_vector_complex_const_subvector_with_stride(gsl_vector_complex * v, size_t offset, size_t stride, size_t n)


  gsl_vector_view  gsl_vector_complex_real(gsl_vector_complex *v)

  gsl_vector_view  gsl_vector_complex_imag(gsl_vector_complex *v)

  gsl_vector_const_view  gsl_vector_complex_const_real(gsl_vector_complex *v)

  gsl_vector_const_view  gsl_vector_complex_const_imag(gsl_vector_complex *v)


  # Operations
  gsl_complex  gsl_vector_complex_get(gsl_vector_complex * v, size_t i)

  void  gsl_vector_complex_set(gsl_vector_complex * v, size_t i, gsl_complex x)

  gsl_complex *  gsl_vector_complex_ptr(gsl_vector_complex * v, size_t i)

  gsl_complex *  gsl_vector_complex_const_ptr(gsl_vector_complex * v, size_t i)

  void  gsl_vector_complex_set_zero(gsl_vector_complex * v)

  void  gsl_vector_complex_set_all(gsl_vector_complex * v, gsl_complex x)

  int  gsl_vector_complex_set_basis(gsl_vector_complex * v, size_t i)

  # Reading and writing vector_complexs
  int  gsl_vector_complex_fread(FILE * stream, gsl_vector_complex * v)

  int  gsl_vector_complex_fwrite(FILE * stream, gsl_vector_complex * v)

  int  gsl_vector_complex_fscanf(FILE * stream, gsl_vector_complex * v)

  int  gsl_vector_complex_fprintf(FILE * stream, gsl_vector_complex * v, char * format)

  # Copying or exchanging elements
  int  gsl_vector_complex_memcpy(gsl_vector_complex * dest, gsl_vector_complex * src)

  int  gsl_vector_complex_reverse(gsl_vector_complex * v)

  int  gsl_vector_complex_swap(gsl_vector_complex * v, gsl_vector_complex * w)

  int  gsl_vector_complex_swap_elements(gsl_vector_complex * v, size_t i, size_t j)

  int  gsl_vector_complex_isnull(gsl_vector_complex * v)


