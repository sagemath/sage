# sage_setup: distribution = sagemath-modules
from sage.libs.gsl.fft cimport *

cdef class FastFourierTransform_base:
    pass

cdef class FastFourierTransform_complex(FastFourierTransform_base):
    cdef double * data
    cdef size_t n
    cdef size_t stride
    cdef gsl_fft_complex_wavetable * wt
    cdef gsl_fft_complex_workspace * mem
