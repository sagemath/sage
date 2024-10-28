"""
Discrete Wavelet Transform

Wraps ``GSL's gsl_wavelet_transform_forward()``,
and ``gsl_wavelet_transform_inverse()`` and creates plot methods.

AUTHOR:

- Josh Kantor (2006-10-07)  - initial version
- David Joyner (2006-10-09) - minor changes to docstrings and examples.
"""

# ****************************************************************************
#       Copyright (C) 2006 Joshua Kantor <jkantor@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


def WaveletTransform(n, wavelet_type, wavelet_k):
    r"""
    This function initializes an GSLDoubleArray of length n which
    can perform a discrete wavelet transform.

    INPUT:

    - ``n`` -- a power of 2
    - ``T`` -- the data in the GSLDoubleArray must be real
    - ``wavelet_type`` -- the name of the type of wavelet, valid choices are:

      * ``'daubechies'``
      * ``'daubechies_centered'``
      * ``'haar'``
      * ``'haar_centered'``
      * ``'bspline'``
      * ``'bspline_centered'``

    For daubechies wavelets, ``wavelet_k`` specifies a daubechie wavelet
    with `k/2` vanishing moments. `k = 4,6,...,20` for `k` even are the
    only ones implemented.

    For Haar wavelets, ``wavelet_k`` must be 2.

    For bspline wavelets, ``wavelet_k`` of `103,105,202,204,206,208,301,
    305,307,309` will give biorthogonal B-spline wavelets of order `(i,j)`
    where ``wavelet_k`` is `100*i+j`.
    The wavelet transform uses `J = \log_2(n)` levels.

    OUTPUT:

    An array of the form
    `(s_{-1,0}, d_{0,0}, d_{1,0}, d_{1,1}, d_{2,0}, \ldots, d_{J-1,2^{J-1}-1})`
    for `d_{j,k}` the detail coefficients of level `j`.
    The centered forms align the coefficients of the sub-bands on edges.

    EXAMPLES::

        sage: a = WaveletTransform(128,'daubechies',4)
        sage: for i in range(1, 11):
        ....:     a[i] = 1
        ....:     a[128-i] = 1
        sage: a.plot().show(ymin=0)                                                     # needs sage.plot
        sage: a.forward_transform()
        sage: a.plot().show()                                                           # needs sage.plot
        sage: a = WaveletTransform(128,'haar',2)
        sage: for i in range(1, 11): a[i] = 1; a[128-i] = 1
        sage: a.forward_transform()
        sage: a.plot().show(ymin=0)                                                     # needs sage.plot
        sage: a = WaveletTransform(128,'bspline_centered',103)
        sage: for i in range(1, 11): a[i] = 1; a[100+i] = 1
        sage: a.forward_transform()
        sage: a.plot().show(ymin=0)                                                     # needs sage.plot

    This example gives a simple example of wavelet compression::

        sage: # needs sage.symbolic
        sage: a = DWT(2048,'daubechies',6)
        sage: for i in range(2048): a[i]=float(sin((i*5/2048)**2))
        sage: a.plot().show()                   # long time (7s on sage.math, 2011), needs sage.plot
        sage: a.forward_transform()
        sage: for i in range(1800): a[2048-i-1] = 0
        sage: a.backward_transform()
        sage: a.plot().show()                   # long time (7s on sage.math, 2011), needs sage.plot
    """
    cdef size_t _n, _k
    _n = int(n)
    if _n < 0:
        raise ValueError("n must be nonnegative.")
    _k = int(wavelet_k)
    if not is2pow(_n):
        raise NotImplementedError("discrete wavelet transform only implemented when n is a 2-power")
    return DiscreteWaveletTransform(_n,1,wavelet_type,_k)


DWT = WaveletTransform


cdef class DiscreteWaveletTransform(GSLDoubleArray):
    """
    Discrete wavelet transform class.
    """
    def __cinit__(self,size_t n,size_t stride, wavelet_type, size_t wavelet_k):
        self.wavelet = NULL
        self.workspace = NULL

    def __init__(self,size_t n,size_t stride, wavelet_type, size_t wavelet_k):
        if not is2pow(n):
            raise NotImplementedError("discrete wavelet transform only implemented when n is a 2-power")
        GSLDoubleArray.__init__(self,n,stride)
        if wavelet_type=="daubechies":
            self.wavelet = <gsl_wavelet*> gsl_wavelet_alloc(gsl_wavelet_daubechies, wavelet_k)
        elif wavelet_type == "daubechies_centered":
            self.wavelet = <gsl_wavelet*> gsl_wavelet_alloc(gsl_wavelet_daubechies_centered,wavelet_k)
        elif wavelet_type == "haar":
            self.wavelet = <gsl_wavelet *> gsl_wavelet_alloc(gsl_wavelet_haar,wavelet_k)
        elif wavelet_type == "haar_centered":
            self.wavelet = <gsl_wavelet*> gsl_wavelet_alloc(gsl_wavelet_haar_centered,wavelet_k)
        elif wavelet_type == "bspline":
            self.wavelet = <gsl_wavelet*> gsl_wavelet_alloc(gsl_wavelet_bspline,wavelet_k)
        elif wavelet_type == "bspline_centered":
            self.wavelet = <gsl_wavelet*> gsl_wavelet_alloc(gsl_wavelet_bspline_centered,wavelet_k)
        self.workspace = <gsl_wavelet_workspace*> gsl_wavelet_workspace_alloc(n)

    def __dealloc__(self):
        if self.wavelet != NULL:
            gsl_wavelet_free(self.wavelet)
            gsl_wavelet_workspace_free(self.workspace)

    def forward_transform(self):
        gsl_wavelet_transform_forward(self.wavelet,self.data,self.stride,self.n,self.workspace)

    def backward_transform(self):
        gsl_wavelet_transform_inverse(self.wavelet,self.data,self.stride,self.n,self.workspace)

    def plot(self, xmin=None, xmax=None, **args):
        from sage.plot.point import point

        cdef int i
        cdef double x
        v = []
        if xmin is None:
            x_min = 0
        if xmax is None:
            x_max = self.n
        for i in range(x_min, x_max):
            x = self.data[i]
            if i > 0:
                v.append(point([(i, x)], hue=(1, 1, 1), **args))
        return sum(v)


def is2pow(unsigned int n):
    while n and not n % 2:
        n = n >> 1
    return n == 1
