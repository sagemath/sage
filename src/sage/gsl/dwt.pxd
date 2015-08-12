from sage.libs.gsl.wavelet cimport *
cimport gsl_array


cdef class DiscreteWaveletTransform(gsl_array.GSLDoubleArray):
  cdef gsl_wavelet* wavelet
  cdef gsl_wavelet_workspace* workspace
