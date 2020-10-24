# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from .types cimport *

cdef extern from "gsl/gsl_matrix_complex_double.h":
  # Allocation
  gsl_matrix_complex *  gsl_matrix_complex_alloc(size_t n1, size_t n2)

  gsl_matrix_complex *  gsl_matrix_complex_calloc(size_t n1, size_t n2)

  gsl_matrix_complex * gsl_matrix_complex_alloc_from_block (gsl_block_complex * b,  size_t offset,  size_t n1,  size_t n2,  size_t d2)

  gsl_matrix_complex * gsl_matrix_complex_alloc_from_matrix (gsl_matrix_complex * b,  size_t k1,  size_t k2,  size_t n1,  size_t n2)

  gsl_vector_complex * gsl_vector_complex_alloc_row_from_matrix (gsl_matrix_complex * m,  size_t i)

  gsl_vector_complex * gsl_vector_complex_alloc_col_from_matrix (gsl_matrix_complex * m,  size_t j)

  void  gsl_matrix_complex_free(gsl_matrix_complex * m)

  # Views
  gsl_matrix_complex_view  gsl_matrix_complex_submatrix(gsl_matrix_complex * m, size_t k1, size_t k2, size_t n1, size_t n2)

  gsl_vector_complex_view  gsl_matrix_complex_row(gsl_matrix_complex * m, size_t i)

  gsl_vector_complex_view  gsl_matrix_complex_column(gsl_matrix_complex * m, size_t j)

  gsl_vector_complex_view  gsl_matrix_complex_diagonal(gsl_matrix_complex * m)

  gsl_vector_complex_view  gsl_matrix_complex_subdiagonal(gsl_matrix_complex * m, size_t k)

  gsl_vector_complex_view  gsl_matrix_complex_superdiagonal(gsl_matrix_complex * m, size_t k)

  gsl_matrix_complex_view gsl_matrix_complex_view_array (double * base,  size_t n1,  size_t n2)

  gsl_matrix_complex_view gsl_matrix_complex_view_array_with_tda (double * base,  size_t n1,  size_t n2,  size_t tda)

  gsl_matrix_complex_view  gsl_matrix_complex_view_vector(gsl_vector_complex * v, size_t n1, size_t n2)

  gsl_matrix_complex_view  gsl_matrix_complex_view_vector_with_tda(gsl_vector_complex* v, size_t n1, size_t n2, size_t tda)

  gsl_matrix_complex_const_view  gsl_matrix_complex_const_submatrix(gsl_matrix_complex * m, size_t k1, size_t k2, size_t n1, size_t n2)

  gsl_vector_complex_const_view  gsl_matrix_complex_const_row(gsl_matrix_complex * m, size_t i)

  gsl_vector_complex_const_view  gsl_matrix_complex_const_column(gsl_matrix_complex * m, size_t j)

  gsl_vector_complex_const_view  gsl_matrix_complex_const_diagonal(gsl_matrix_complex * m)

  gsl_vector_complex_const_view  gsl_matrix_complex_const_subdiagonal(gsl_matrix_complex * m, size_t k)

  gsl_vector_complex_const_view  gsl_matrix_complex_const_superdiagonal(gsl_matrix_complex * m, size_t k)

  gsl_matrix_complex_const_view gsl_matrix_complex_const_view_array ( double * base,  size_t n1,  size_t n2)

  gsl_matrix_complex_const_view gsl_matrix_complex_const_view_array_with_tda ( double * base,  size_t n1,  size_t n2,  size_t tda)

  gsl_matrix_complex_const_view  gsl_matrix_complex_const_view_vector(gsl_vector_complex * v, size_t n1, size_t n2)

  gsl_matrix_complex_const_view  gsl_matrix_complex_const_view_vector_with_tda(gsl_vector_complex * v, size_t n1, size_t n2, size_t tda)

  # Operations
  gsl_complex  gsl_matrix_complex_get(gsl_matrix_complex * m, size_t i, size_t j)

  void  gsl_matrix_complex_set(gsl_matrix_complex * m, size_t i, size_t j, gsl_complex x)

  gsl_complex *  gsl_matrix_complex_ptr(gsl_matrix_complex * m, size_t i, size_t j)

  gsl_complex *  gsl_matrix_complex_const_ptr(gsl_matrix_complex * m, size_t i, size_t j)

  void  gsl_matrix_complex_set_zero(gsl_matrix_complex * m)

  void  gsl_matrix_complex_set_identity(gsl_matrix_complex * m)

  void  gsl_matrix_complex_set_all(gsl_matrix_complex * m, gsl_complex x)

  # Reading and writing matrices
  int  gsl_matrix_complex_fread(FILE * stream, gsl_matrix_complex * m)

  int  gsl_matrix_complex_fwrite(FILE * stream, gsl_matrix_complex * m)

  int  gsl_matrix_complex_fscanf(FILE * stream, gsl_matrix_complex * m)

  int  gsl_matrix_complex_fprintf(FILE * stream, gsl_matrix_complex * m, char * format)



  # Copying or exchanging elements
  int  gsl_matrix_complex_memcpy(gsl_matrix_complex * dest, gsl_matrix_complex * src)

  int  gsl_matrix_complex_swap(gsl_matrix_complex * m1, gsl_matrix_complex * m2)

  int  gsl_matrix_complex_swap_rows(gsl_matrix_complex * m, size_t i, size_t j)

  int  gsl_matrix_complex_swap_columns(gsl_matrix_complex * m, size_t i, size_t j)

  int  gsl_matrix_complex_swap_rowcol(gsl_matrix_complex * m, size_t i, size_t j)

  int  gsl_matrix_complex_transpose(gsl_matrix_complex * m)

  int  gsl_matrix_complex_transpose_memcpy(gsl_matrix_complex * dest, gsl_matrix_complex * src)

  int  gsl_matrix_complex_isnull(gsl_matrix_complex * m)

  # Matrix operations
  int  gsl_matrix_complex_add(gsl_matrix_complex * a, gsl_matrix_complex * b)

  int  gsl_matrix_complex_sub(gsl_matrix_complex * a, gsl_matrix_complex * b)

  int  gsl_matrix_complex_mul_elements(gsl_matrix_complex * a, gsl_matrix_complex * b)

  int  gsl_matrix_complex_div_elements(gsl_matrix_complex * a, gsl_matrix_complex * b)

  int  gsl_matrix_complex_scale(gsl_matrix_complex * a, gsl_complex x)

  int  gsl_matrix_complex_add_constant(gsl_matrix_complex * a, gsl_complex x)

  int gsl_matrix_complex_add_diagonal (gsl_matrix_complex * a,  gsl_complex x)

  # The functions below are obsolete
  int  gsl_matrix_complex_get_row(gsl_vector * v, gsl_matrix_complex * m, size_t i)

  int  gsl_matrix_complex_get_col(gsl_vector * v, gsl_matrix_complex * m, size_t j)

  int  gsl_matrix_complex_set_row(gsl_matrix_complex * m, size_t i, gsl_vector * v)

  int  gsl_matrix_complex_set_col(gsl_matrix_complex * m, size_t j, gsl_vector * v)


