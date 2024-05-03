# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_wavelet.h":

  cdef enum gsl_wavelet_direction:
    gsl_wavelet_forward  =  1
    gsl_wavelet_backward = -1

  ctypedef struct gsl_wavelet_type:
    char *name
    int (*init)(double **h1, double **g1, double **h2, double **g2, size_t * nc, size_t * offset, size_t member)


  ctypedef struct gsl_wavelet:
    gsl_wavelet_type *type
    double *h1
    double *g1
    double *h2
    double *g2
    size_t nc
    size_t offset

  ctypedef struct gsl_wavelet_workspace:
    double *scratch
    size_t n

  gsl_wavelet_type *gsl_wavelet_daubechies
  gsl_wavelet_type *gsl_wavelet_daubechies_centered
  gsl_wavelet_type *gsl_wavelet_haar
  gsl_wavelet_type *gsl_wavelet_haar_centered
  gsl_wavelet_type *gsl_wavelet_bspline
  gsl_wavelet_type *gsl_wavelet_bspline_centered

  gsl_wavelet *gsl_wavelet_alloc (gsl_wavelet_type * T, size_t k)
  void gsl_wavelet_free (gsl_wavelet * w)
  char *gsl_wavelet_name (gsl_wavelet * w)

  gsl_wavelet_workspace *gsl_wavelet_workspace_alloc (size_t n)
  void gsl_wavelet_workspace_free (gsl_wavelet_workspace * work)

  int gsl_wavelet_transform (gsl_wavelet * w, double *data, size_t stride, size_t n, gsl_wavelet_direction dir, gsl_wavelet_workspace * work)

  int gsl_wavelet_transform_forward (gsl_wavelet * w, double *data, size_t stride, size_t n, gsl_wavelet_workspace * work)

  int gsl_wavelet_transform_inverse (gsl_wavelet * w, double *data, size_t stride, size_t n, gsl_wavelet_workspace * work)
