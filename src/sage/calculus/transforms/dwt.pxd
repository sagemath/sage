# sage_setup: distribution = sagemath-modules
from sage.libs.gsl.wavelet cimport *
from sage.libs.gsl.array cimport GSLDoubleArray

cdef class DiscreteWaveletTransform(GSLDoubleArray):
    cdef gsl_wavelet* wavelet
    cdef gsl_wavelet_workspace* workspace
