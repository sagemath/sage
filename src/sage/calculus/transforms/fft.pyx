"""
Fast Fourier Transforms Using GSL

AUTHORS:

- William Stein (2006-9): initial file (radix2)
- D. Joyner (2006-10): Minor modifications (from radix2 to general case\
                        and some documentation).
- M. Hansen (2013-3): Fix radix2 backwards transformation
- L.F. Tabera Alonso (2013-3): Documentation
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport sig_malloc, sig_free

from sage.rings.integer import Integer
from sage.rings.complex_mpfr import ComplexNumber


def FastFourierTransform(size, base_ring=None):
    """
    Create an array for fast Fourier transform conversion using gsl.

    INPUT:

    - ``size`` -- the size of the array
    - ``base_ring`` -- unused (2013-03)

    EXAMPLES:

    We create an array of the desired size::

        sage: a = FastFourierTransform(8)
        sage: a
        [(0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)]

    Now, set the values of the array::

        sage: for i in range(8): a[i] = i + 1
        sage: a
        [(1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0), (5.0, 0.0), (6.0, 0.0), (7.0, 0.0), (8.0, 0.0)]

    We can perform the forward Fourier transform on the array::

        sage: a.forward_transform()
        sage: a                       #abs tol 1e-2
        [(36.0, 0.0), (-4.00, 9.65), (-4.0, 4.0), (-4.0, 1.65), (-4.0, 0.0), (-4.0, -1.65), (-4.0, -4.0), (-4.0, -9.65)]

    And backwards::

        sage: a.backward_transform()
        sage: a                       #abs tol 1e-2
        [(8.0, 0.0), (16.0, 0.0), (24.0, 0.0), (32.0, 0.0), (40.0, 0.0), (48.0, 0.0), (56.0, 0.0), (64.0, 0.0)]

    Other example::

        sage: a = FastFourierTransform(128)
        sage: for i in range(1, 11):
        ....:    a[i] = 1
        ....:    a[128-i] = 1
        sage: a[:6:2]
        [(0.0, 0.0), (1.0, 0.0), (1.0, 0.0)]
        sage: a.plot().show(ymin=0)                                                     # needs sage.plot
        sage: a.forward_transform()
        sage: a.plot().show()                                                           # needs sage.plot
    """
    return FastFourierTransform_complex(int(size))


FFT = FastFourierTransform


cdef class FastFourierTransform_base:
    pass

cdef class FastFourierTransform_complex(FastFourierTransform_base):
    """
    Wrapper class for GSL's fast Fourier transform.
    """

    def __init__(self, size_t n, size_t stride=1):
        """
        Create an array-like object of fixed size that will contain the vector to
        apply the Fast Fourier Transform.

        INPUT:

        - ``n`` -- integer, the size of the array
        - ``stride`` -- the stride to be applied when manipulating the array

        EXAMPLES::

            sage: a = FastFourierTransform(1) # indirect doctest
            sage: a
            [(0.0, 0.0)]
        """
        self.n = n
        self.stride = stride
        self.data = <double*>sig_malloc(sizeof(double)*(2*n))
        cdef int i
        for i in range(2 * n):
            self.data[i] = 0

    def __dealloc__(self):
        """
        Frees allocated memory.

        EXAMPLES::

            sage: a = FastFourierTransform(128)
            sage: del a
        """
        sig_free(self.data)

    def __len__(self):
        """
        Return the size of the array.

        OUTPUT: the size of the array

        EXAMPLES::

            sage: a = FastFourierTransform(48)
            sage: len(a)
            48
        """
        return self.n

    def __setitem__(self, size_t i, xy):
        """
        Assign a value to an index of the array.

        Currently the input has to be en element that can be coerced
        to ``float`` or a ``ComplexNumber`` element.

        INPUT:

        - ``i`` -- integer; the index
        - ``xy`` -- an object to store as `i`-th element of the array ``self[i]``

        EXAMPLES::

            sage: # needs sage.rings.mpfr sage.symbolic
            sage: I = CC(I)
            sage: a = FastFourierTransform(4)
            sage: a[0] = 1
            sage: a[1] = I
            sage: a[2] = 1+I
            sage: a[3] = (2,2)
            sage: a
            [(1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (2.0, 2.0)]
            sage: I = CDF(I)
            sage: a[1] = I
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 1.0*I to float; use abs() or real_part() as desired

        TESTS:

        Verify that :issue:`10758` is fixed. ::

            sage: F = FFT(1)
            sage: F[0] = (1,1)
            sage: F[0] = 1
            sage: F[0]
            (1.0, 0.0)
        """
        # just set real for now
        if i < 0 or i >= self.n:
            raise IndexError
        if isinstance(xy, (tuple, ComplexNumber)):
            self.data[2*i] = xy[0]
            self.data[2*i+1] = xy[1]
        else:
            self.data[2*i] = xy
            self.data[2*i+1] = 0

    def __getitem__(self, i):
        """
        Get the `i`-th element of the array.

        INPUT:

        - ``i`` -- integer

        OUTPUT: the `i`-th element of the array ``self[i]``

        EXAMPLES::

            sage: a = FastFourierTransform(4)
            sage: a[0]
            (0.0, 0.0)
            sage: a[0] = 1
            sage: a[0] == (1,0)
            True
        """
        if isinstance(i, slice):
            start, stop, step = i.indices(self.n)
            return list(self)[start:stop:step]
        else:
            if i < 0 or i >= self.n:
                raise IndexError
            return self.data[2*i], self.data[2*i+1]

    def __repr__(self):
        """
        String representation of the array.

        OUTPUT:

        - A string representing this array. The complex numbers are
            presented as a tuple of two float elements.

        EXAMPLES::

            sage: a = FastFourierTransform(4)
            sage: for i in range(4): a[i] = i
            sage: a
            [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0)]
        """
        return str(list(self))

    def _plot_polar(self, xmin, xmax, **args):
        """
        Plot a slice of the array using polar coordinates.

        INPUT:

        - ``xmin`` -- the lower bound of the slice to plot
        - ``xmax`` -- the upper bound of the slice to plot
        - ``**args`` -- passed on to the line plotting function

        OUTPUT: a plot of the array interpreting each element as polar coordinates

        This method should not be called directly. See :meth:`plot` for the details.

        EXAMPLES::

            sage: a = FastFourierTransform(4)
            sage: a._plot_polar(0,2)                                                    # needs sage.plot sage.symbolic
            Graphics object consisting of 2 graphics primitives
        """
        from sage.plot.point import point
        from sage.symbolic.constants import pi, I

        cdef int i
        v = []

        pi = pi.n()
        I = I.n()
        s = 1/(3*pi)   # so arg gets scaled between -1/3 and 1/3.

        for i in range(xmin, xmax):
            z = self.data[2*i] + I*self.data[2*i+1]
            mag = z.abs()
            arg = z.arg()*s
            v.append(point((i, mag), hue=arg, **args))
        return sum(v)

    def _plot_rect(self, xmin, xmax, **args):
        """
        Plot a slice of the array.

        INPUT:

        - ``xmin`` -- the lower bound of the slice to plot
        - ``xmax`` -- the upper bound of the slice to plot
        - ``**args`` -- passed on to the line plotting function

        OUTPUT: a plot of the array

        This method should not be called directly. See :meth:`plot` for the details.

        EXAMPLES::

            sage: a = FastFourierTransform(4)
            sage: a._plot_rect(0,3)                                                     # needs sage.plot
            Graphics object consisting of 3 graphics primitives
        """
        from sage.plot.point import point

        cdef int i
        cdef double x, h
        v = []

        for i in range(xmin, xmax):
            x = self.data[2*i]
            h = self.data[2*i+1]
            v.append(point((i, x), hue=h, **args))
        return sum(v)

    def plot(self, style='rect', xmin=None, xmax=None, **args):
        """
        Plot a slice of the array.

        - ``style`` -- style of the plot, options are ``'rect'`` or ``'polar'``

          - ``'rect'`` -- height represents real part, color represents
            imaginary part
          - ``'polar'`` -- height represents absolute value, color
            represents argument

        - ``xmin`` -- the lower bound of the slice to plot; 0 by default
        - ``xmax`` -- the upper bound of the slice to plot; ``len(self)`` by default
        - ``**args`` -- passed on to the line plotting function

        OUTPUT: a plot of the array

        EXAMPLES::

            sage: # needs sage.plot
            sage: a = FastFourierTransform(16)
            sage: for i in range(16): a[i] = (random(),random())
            sage: A = plot(a)
            sage: B = plot(a, style='polar')                                            # needs sage.symbolic
            sage: type(A)
            <class 'sage.plot.graphics.Graphics'>
            sage: type(B)                                                               # needs sage.symbolic
            <class 'sage.plot.graphics.Graphics'>
            sage: a = FastFourierTransform(125)
            sage: b = FastFourierTransform(125)
            sage: for i in range(1, 60): a[i]=1
            sage: for i in range(1, 60): b[i]=1
            sage: a.forward_transform()
            sage: a.inverse_transform()
            sage: a.plot() + b.plot()
            Graphics object consisting of 250 graphics primitives
        """
        if xmin is None:
            xmin = 0
        else:
            xmin = int(xmin)
        if xmax is None:
            xmax = self.n
        else:
            xmax = int(xmax)
        if style == 'rect':
            return self._plot_rect(xmin, xmax, **args)
        elif style == 'polar':
            return self._plot_polar(xmin, xmax, **args)
        else:
            raise ValueError("unknown style '%s'" % style)

    def forward_transform(self):
        """
        Compute the in-place forward Fourier transform of this data
        using the Cooley-Tukey algorithm.

        OUTPUT: none, the transformation is done in-place

        If the number of sample points in the input is a power of 2 then the
        gsl function ``gsl_fft_complex_radix2_forward`` is automatically called.
        Otherwise, ``gsl_fft_complex_forward`` is called.

        EXAMPLES::

            sage: a = FastFourierTransform(4)
            sage: for i in range(4): a[i] = i
            sage: a.forward_transform()
            sage: a #abs tol 1e-2
            [(6.0, 0.0), (-2.0, 2.0), (-2.0, 0.0), (-2.0, -2.0)]
        """
        cdef gsl_fft_complex_wavetable * wt
        cdef gsl_fft_complex_workspace * mem
        N = Integer(self.n)
        e = N.exact_log(2)
        if N==2**e:
            gsl_fft_complex_radix2_forward(self.data, self.stride, self.n)
        else:
            mem = gsl_fft_complex_workspace_alloc(self.n)
            wt = gsl_fft_complex_wavetable_alloc(self.n)
            gsl_fft_complex_forward(self.data, self.stride, self.n, wt, mem)
            gsl_fft_complex_workspace_free(mem)
            gsl_fft_complex_wavetable_free(wt)

    def inverse_transform(self):
        """
        Compute the in-place inverse Fourier transform of this data
        using the Cooley-Tukey algorithm.

        OUTPUT: none, the transformation is done in-place

        If the number of sample points in the input is a power of 2 then the
        function ``gsl_fft_complex_radix2_inverse`` is automatically called.
        Otherwise, ``gsl_fft_complex_inverse`` is called.

        This transform is normalized so ``f.forward_transform().inverse_transform() == f``
        modulo round-off errors. See also :meth:`backward_transform`.

        EXAMPLES::

            sage: a = FastFourierTransform(125)
            sage: b = FastFourierTransform(125)
            sage: for i in range(1, 60): a[i]=1
            sage: for i in range(1, 60): b[i]=1
            sage: a.forward_transform()
            sage: a.inverse_transform()
            sage: a.plot() + b.plot()                                                   # needs sage.plot
            Graphics object consisting of 250 graphics primitives
            sage: abs(sum([CDF(a[i])-CDF(b[i]) for i in range(125)])) < 2**-16
            True

        Here we check it with a power of two::

            sage: a = FastFourierTransform(128)
            sage: b = FastFourierTransform(128)
            sage: for i in range(1, 60): a[i]=1
            sage: for i in range(1, 60): b[i]=1
            sage: a.forward_transform()
            sage: a.inverse_transform()
            sage: a.plot() + b.plot()                                                   # needs sage.plot
            Graphics object consisting of 256 graphics primitives
        """
        cdef gsl_fft_complex_wavetable * wt
        cdef gsl_fft_complex_workspace * mem
        N = Integer(self.n)
        e = N.exact_log(2)
        if N==2**e:
            gsl_fft_complex_radix2_inverse(self.data, self.stride, self.n)
        else:
            mem = gsl_fft_complex_workspace_alloc(self.n)
            wt = gsl_fft_complex_wavetable_alloc(self.n)
            gsl_fft_complex_inverse(self.data, self.stride, self.n, wt, mem)
            gsl_fft_complex_workspace_free(mem)
            gsl_fft_complex_wavetable_free(wt)

    def backward_transform(self):
        """
        Compute the in-place backwards Fourier transform of this data
        using the Cooley-Tukey algorithm.

        OUTPUT: none, the transformation is done in-place

        This is the same as :meth:`inverse_transform` but lacks normalization
        so that ``f.forward_transform().backward_transform() == n*f``. Where
        ``n`` is the size of the array.

        EXAMPLES::

            sage: a = FastFourierTransform(125)
            sage: b = FastFourierTransform(125)
            sage: for i in range(1, 60): a[i]=1
            sage: for i in range(1, 60): b[i]=1
            sage: a.forward_transform()
            sage: a.backward_transform()
            sage: (a.plot() + b.plot()).show(ymin=0)    # long time (2s on sage.math, 2011), needs sage.plot
            sage: abs(sum([CDF(a[i])/125-CDF(b[i]) for i in range(125)])) < 2**-16
            True

        Here we check it with a power of two::

            sage: a = FastFourierTransform(128)
            sage: b = FastFourierTransform(128)
            sage: for i in range(1, 60): a[i]=1
            sage: for i in range(1, 60): b[i]=1
            sage: a.forward_transform()
            sage: a.backward_transform()
            sage: (a.plot() + b.plot()).show(ymin=0)                                    # needs sage.plot
        """
        cdef gsl_fft_complex_wavetable * wt
        cdef gsl_fft_complex_workspace * mem
        N = Integer(self.n)
        e = N.exact_log(2)
        if N==2**e:
            gsl_fft_complex_radix2_backward(self.data, self.stride, self.n)
        else:
            mem = gsl_fft_complex_workspace_alloc(self.n)
            wt = gsl_fft_complex_wavetable_alloc(self.n)
            gsl_fft_complex_backward(self.data, self.stride, self.n, wt, mem)
            gsl_fft_complex_workspace_free(mem)
            gsl_fft_complex_wavetable_free(wt)

cdef class FourierTransform_complex:
    pass

cdef class FourierTransform_real:
    pass
