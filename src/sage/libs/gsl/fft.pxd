# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from .types cimport *

cdef extern from "gsl/gsl_fft.h":
  cdef enum gsl_fft_direction:
    forward = -1
    backward = 1

cdef extern from "gsl/gsl_fft_complex.h":

  # Power of 2 routines
  int  gsl_fft_complex_radix2_forward(gsl_complex_packed_array data, size_t stride, size_t n)

  int  gsl_fft_complex_radix2_backward(gsl_complex_packed_array data, size_t stride, size_t n)

  int  gsl_fft_complex_radix2_inverse(gsl_complex_packed_array data, size_t stride, size_t n)

  int  gsl_fft_complex_radix2_transform(gsl_complex_packed_array data, size_t stride, size_t n, gsl_fft_direction sign)

  int  gsl_fft_complex_radix2_dif_forward(gsl_complex_packed_array data, size_t stride, size_t n)

  int  gsl_fft_complex_radix2_dif_backward(gsl_complex_packed_array data, size_t stride, size_t n)

  int  gsl_fft_complex_radix2_dif_inverse(gsl_complex_packed_array data, size_t stride, size_t n)

  int  gsl_fft_complex_radix2_dif_transform(gsl_complex_packed_array data, size_t stride, size_t n, gsl_fft_direction sign)

  int gsl_fft_complex_bitreverse_order (gsl_complex_packed_array data, size_t stride, size_t n, size_t n_bits)


  # Mixed Radix general-N routines
  # or ctypedef struct gsl_fft_complex_wavetable
  ctypedef struct gsl_fft_complex_wavetable:
    size_t n
    size_t nf
    size_t factor[64]
    gsl_complex * trig
    gsl_complex * twiddle[64]

  ctypedef struct gsl_fft_complex_workspace

  gsl_fft_complex_wavetable *  gsl_fft_complex_wavetable_alloc(size_t n)

  void  gsl_fft_complex_wavetable_free(gsl_fft_complex_wavetable * wavetable)

  gsl_fft_complex_workspace *  gsl_fft_complex_workspace_alloc(size_t n)

  void  gsl_fft_complex_workspace_free(gsl_fft_complex_workspace * workspace)

  int gsl_fft_complex_memcpy (gsl_fft_complex_wavetable * dest, gsl_fft_complex_wavetable * src)

  int  gsl_fft_complex_forward(gsl_complex_packed_array data, size_t stride, size_t n, gsl_fft_complex_wavetable * wavetable, gsl_fft_complex_workspace * work)

  int  gsl_fft_complex_backward(gsl_complex_packed_array data, size_t stride, size_t n, gsl_fft_complex_wavetable * wavetable, gsl_fft_complex_workspace * work)

  int  gsl_fft_complex_inverse(gsl_complex_packed_array data, size_t stride, size_t n, gsl_fft_complex_wavetable * wavetable, gsl_fft_complex_workspace * work)

  int  gsl_fft_complex_transform(gsl_complex_packed_array data, size_t stride, size_t n, gsl_fft_complex_wavetable * wavetable, gsl_fft_complex_workspace * work, gsl_fft_direction sign)

  # end of gsl_fft_complex.h

cdef extern from "gsl/gsl_fft_real.h":

  int  gsl_fft_real_radix2_transform(double data[], size_t stride, size_t n)

  # or ctypedef struct gsl_fft_real_wavetable
  cdef struct gsl_fft_real_wavetable:
    size_t n
    size_t nf
    size_t factor[64]
    gsl_complex *twiddle[64]
    gsl_complex *trig

  cdef struct gsl_fft_real_workspace

  gsl_fft_real_wavetable *  gsl_fft_real_wavetable_alloc(size_t n)

  void  gsl_fft_real_wavetable_free(gsl_fft_real_wavetable * wavetable)

  gsl_fft_real_workspace *  gsl_fft_real_workspace_alloc(size_t n)

  void  gsl_fft_real_workspace_free(gsl_fft_real_workspace * workspace)

  int  gsl_fft_real_transform(double data[], size_t stride, size_t n, gsl_fft_real_wavetable * wavetable, gsl_fft_real_workspace * work)

  int  gsl_fft_real_unpack(double real_coefficient[], gsl_complex_packed_array complex_coefficient[], size_t stride, size_t n)

  # end of gsl_fft_real.h

cdef extern from "gsl/gsl_fft_halfcomplex.h":

  int  gsl_fft_halfcomplex_radix2_backward(double data[], size_t stride, size_t n)

  int  gsl_fft_halfcomplex_radix2_inverse(double data[], size_t stride, size_t n)

  int gsl_fft_halfcomplex_radix2_transform (double data[],  size_t stride,  size_t n)

  # or ctypedef struct gsl_fft_halfcomplex_wavetable
  cdef struct gsl_fft_halfcomplex_wavetable:
    size_t n
    size_t nf
    size_t factor[64]
    gsl_complex *twiddle[64]
    gsl_complex *trig

  gsl_fft_halfcomplex_wavetable *  gsl_fft_halfcomplex_wavetable_alloc(size_t n)

  void  gsl_fft_halfcomplex_wavetable_free(gsl_fft_halfcomplex_wavetable * wavetable)

  int gsl_fft_halfcomplex_backward (double data[],  size_t stride,  size_t n,  gsl_fft_halfcomplex_wavetable * wavetable, gsl_fft_real_workspace * work)

  int gsl_fft_halfcomplex_inverse (double data[],  size_t stride,  size_t n,  gsl_fft_halfcomplex_wavetable * wavetable, gsl_fft_real_workspace * work)

  int gsl_fft_halfcomplex_transform (double data[],  size_t stride,  size_t n,  gsl_fft_halfcomplex_wavetable * wavetable, gsl_fft_real_workspace * work)

  int  gsl_fft_halfcomplex_unpack(double halfcomplex_coefficient[], gsl_complex_packed_array complex_coefficient, size_t stride, size_t n)

  int gsl_fft_halfcomplex_radix2_unpack ( double halfcomplex_coefficient[], double complex_coefficient[],  size_t stride,  size_t n)
