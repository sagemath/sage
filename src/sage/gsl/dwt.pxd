from sage.libs.gsl.wavelet cimport *
from .gsl_array cimport GSLDoubleArray


cdef class DiscreteWaveletTransform(GSLDoubleArray):
    cdef gsl_wavelet* wavelet
    cdef gsl_wavelet_workspace* workspace
