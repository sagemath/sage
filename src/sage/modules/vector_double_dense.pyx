# sage_setup: distribution = sagemath-modules
# sage.doctest: optional - numpy
r"""
Dense vectors using a NumPy backend

This serves as a base class for dense vectors over Real Double Field and
Complex Double Field

EXAMPLES::

    sage: # needs sage.symbolic
    sage: v = vector(CDF,[(1,-1), (2,pi), (3,5)]); v
    (1.0 - 1.0*I, 2.0 + 3.141592653589793*I, 3.0 + 5.0*I)
    sage: type(v)
    <class 'sage.modules.vector_complex_double_dense.Vector_complex_double_dense'>
    sage: parent(v)
    Vector space of dimension 3 over Complex Double Field
    sage: v[0] = 5
    sage: v
    (5.0, 2.0 + 3.141592653589793*I, 3.0 + 5.0*I)
    sage: loads(dumps(v)) == v
    True

    sage: v = vector(RDF, [1,2,3,4]); v
    (1.0, 2.0, 3.0, 4.0)
    sage: loads(dumps(v)) == v
    True

AUTHORS:

- Jason Grout, Oct 2008: switch to numpy backend, factored out
  ``Vector_double_dense`` class
- Josh Kantor
- William Stein
"""

#*****************************************************************************
#       Copyright (C) 2006-2010 William Stein <wstein@gmail.com>
#       Copyright (C) 2009      Alexandru Ghitza
#       Copyright (C) 2020      Antonio Rojas
#       Copyright (C) 2017      Frédéric Chapoton
#       Copyright (C) 2008-2009 Jason Grout
#       Copyright (C) 2014-2016 Jeroen Demeyer
#       Copyright (C) 2011      Mike Hansen
#       Copyright (C) 2011      Rob Beezer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

cimport numpy
import numpy

from sage.structure.element cimport Element, Vector

from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF

# This is for the NumPy C API (the PyArray... functions) to work
numpy.import_array()


cdef class Vector_double_dense(Vector_numpy_dense):
    """
    Base class for vectors over the Real Double Field and the Complex
    Double Field.  These are supposed to be fast vector operations
    using C doubles. Most operations are implemented using numpy which
    will call the underlying BLAS, if needed, on the system.

    This class cannot be instantiated on its own.  The numpy vector
    creation depends on several variables that are set in the
    subclasses.

    EXAMPLES::

        sage: v = vector(RDF, [1,2,3,4]); v
        (1.0, 2.0, 3.0, 4.0)
        sage: v*v
        30.0
    """

    cpdef _add_(self, right):
        """
        Add two vectors together.

        EXAMPLES::

            sage: A = vector(RDF, range(3))
            sage: A+A
            (0.0, 2.0, 4.0)
        """
        if self._degree == 0:
            from copy import copy
            return copy(self)

        cdef Vector_double_dense _right, _left
        _right = right
        _left = self

        return self._new(_left._vector_numpy + _right._vector_numpy)

    cpdef _sub_(self, right):
        """
        Return ``self - right``.

        EXAMPLES::

            sage: A = vector(RDF, range(3))
            sage: (A-A).is_zero()
            True
        """
        if self._degree == 0:
            from copy import copy
            return copy(self)

        cdef Vector_double_dense _right, _left
        _right = right
        _left = self

        return self._new(_left._vector_numpy - _right._vector_numpy)

    cpdef _dot_product_(self, Vector right):
        """
        Dot product of ``self`` and ``right``.

        EXAMPLES::

            sage: v = vector(RDF, [1,2,3]); w = vector(RDF, [2, 4, -3])
            sage: v*w
            1.0
            sage: w*v
            1.0

        This works correctly for zero-dimensional vectors::

            sage: v = vector(RDF, [])
            sage: v._dot_product_(v)
            0.0
        """
        cdef Vector_double_dense _right, _left
        _right = right
        _left = self

        return self._sage_dtype(numpy.dot(_left._vector_numpy, _right._vector_numpy))

    cpdef _pairwise_product_(self, Vector right):
        """
        Return the component-wise product of ``self`` and ``right``.

        EXAMPLES::

            sage: v = vector(CDF, [1,2,3]); w = vector(CDF, [2, 4, -3])
            sage: v.pairwise_product(w)
            (2.0, 8.0, -9.0)
        """
        if not right.parent() == self.parent():
            right = self.parent().ambient_module()(right)

        if self._degree == 0:
            from copy import copy
            return copy(self)

        cdef Vector_double_dense _right, _left
        _right = right
        _left = self

        return self._new(_left._vector_numpy * _right._vector_numpy)

    cpdef _rmul_(self, Element left):
        """
        Multiply a scalar and vector.

        EXAMPLES::

            sage: v = vector(CDF, range(3))
            sage: 3*v
            (0.0, 3.0, 6.0)
        """
        if self._degree == 0:
            from copy import copy
            return copy(self)

        return self._new(self._python_dtype(left)*self._vector_numpy)

    cpdef _lmul_(self, Element right):
        """
        Multiply a scalar and vector.

        EXAMPLES::

            sage: v = vector(CDF, range(3))
            sage: v*3
            (0.0, 3.0, 6.0)
        """
        if self._degree == 0:
            from copy import copy
            return copy(self)

        return self._new(self._vector_numpy*self._python_dtype(right))

    def inv_fft(self, algorithm='radix2', inplace=False):
        """
        This performs the inverse fast Fourier transform on the vector.

        The Fourier transform can be done in place using the keyword
        inplace=True

        This will be fastest if the vector's length is a power of 2.

        EXAMPLES::

            sage: # needs scipy
            sage: v = vector(CDF,[1,2,3,4])
            sage: w = v.fft()
            sage: max(v - w.inv_fft()) < 1e-12
            True
        """
        return self.fft(direction='backward', algorithm=algorithm, inplace=inplace)

    def fft(self, direction='forward', algorithm='radix2', inplace=False):
        """
        This performs a fast Fourier transform on the vector.

        INPUT:

        - ``direction`` -- string; ``'forward'`` (default) or ``'backward'``

        The algorithm and inplace arguments are ignored.

        This function is fastest if the vector's length is a power of 2.

        EXAMPLES::

            sage: # needs scipy
            sage: v = vector(CDF,[1+2*I,2,3*I,4])
            sage: v.fft()
            (7.0 + 5.0*I, 1.0 + 1.0*I, -5.0 + 5.0*I, 1.0 - 3.0*I)
            sage: v.fft(direction='backward')
            (1.75 + 1.25*I, 0.25 - 0.75*I, -1.25 + 1.25*I, 0.25 + 0.25*I)
            sage: v.fft().fft(direction='backward')
            (1.0 + 2.0*I, 2.0, 3.0*I, 4.0)
            sage: v.fft().parent()
            Vector space of dimension 4 over Complex Double Field
            sage: v.fft(inplace=True)
            sage: v
            (7.0 + 5.0*I, 1.0 + 1.0*I, -5.0 + 5.0*I, 1.0 - 3.0*I)

            sage: # needs scipy
            sage: v = vector(RDF,4,range(4)); v
            (0.0, 1.0, 2.0, 3.0)
            sage: v.fft()
            (6.0, -2.0 + 2.0*I, -2.0, -2.0 - 2.0*I)
            sage: v.fft(direction='backward')
            (1.5, -0.5 - 0.5*I, -0.5, -0.5 + 0.5*I)
            sage: v.fft().fft(direction='backward')
            (0.0, 1.0, 2.0, 3.0)
            sage: v.fft().parent()
            Vector space of dimension 4 over Complex Double Field
            sage: v.fft(inplace=True)
            Traceback (most recent call last):
            ...
            ValueError: inplace can only be True for CDF vectors
        """
        if direction not in ('forward', 'backward'):
            raise ValueError("direction must be either 'forward' or 'backward'")

        if self._degree == 0:
            return self

        import scipy.fftpack

        if inplace:
            if self._sage_dtype is not CDF:
                raise ValueError("inplace can only be True for CDF vectors")
            if direction == 'forward':
                self._vector_numpy = scipy.fftpack.fft(self._vector_numpy, overwrite_x = True)
            else:
                self._vector_numpy = scipy.fftpack.ifft(self._vector_numpy, overwrite_x = True)
        else:
            try:
                fft = scipy.fft.fft
                ifft = scipy.fft.ifft
            except AttributeError:
                fft = scipy.fft
                ifft = scipy.ifft
            V = CDF ** self._degree
            from sage.modules.vector_complex_double_dense import Vector_complex_double_dense
            if direction == 'forward':
                return Vector_complex_double_dense(V, fft(self._vector_numpy))
            else:
                return Vector_complex_double_dense(V, ifft(self._vector_numpy))

    def complex_vector(self):
        """
        Return the associated complex vector, i.e., this vector but with
        coefficients viewed as complex numbers.

        EXAMPLES::

            sage: v = vector(RDF,4,range(4)); v
            (0.0, 1.0, 2.0, 3.0)
            sage: v.complex_vector()
            (0.0, 1.0, 2.0, 3.0)
            sage: v = vector(RDF,0)
            sage: v.complex_vector()
            ()
        """
        return self.change_ring(CDF)

    def zero_at(self, eps):
        r"""
        Return a copy with small entries replaced by zeros.

        This is useful for modifying output from algorithms
        which have large relative errors when producing zero
        elements, e.g. to create reliable doctests.

        INPUT:

        - ``eps`` -- cutoff value

        OUTPUT:

        A modified copy of the vector.  Elements smaller than
        or equal to ``eps`` are replaced with zeroes.  For
        complex vectors, the real and imaginary parts are
        considered individually.

        EXAMPLES::

            sage: v = vector(RDF, [1.0, 2.0, 10^-10, 3.0])
            sage: v.zero_at(1e-8)
            (1.0, 2.0, 0.0, 3.0)
            sage: v.zero_at(1e-12)
            (1.0, 2.0, 1e-10, 3.0)

        For complex numbers the real and imaginary parts are considered
        separately.  ::

            sage: w = vector(CDF, [10^-6 + 5*I, 5 + 10^-6*I, 5 + 5*I, 10^-6 + 10^-6*I])
            sage: w.zero_at(1.0e-4)
            (5.0*I, 5.0, 5.0 + 5.0*I, 0.0)
            sage: w.zero_at(1.0e-8)
            (1e-06 + 5.0*I, 5.0 + 1e-06*I, 5.0 + 5.0*I, 1e-06 + 1e-06*I)
        """
        import sage.rings.complex_double
        cdef Vector_double_dense v
        eps = float(eps)
        out = self._vector_numpy.copy()
        if self._sage_dtype is sage.rings.complex_double.CDF:
            out.real[numpy.abs(out.real) <= eps] = 0
            out.imag[numpy.abs(out.imag) <= eps] = 0
        else:
            out[numpy.abs(out) <= eps] = 0
        v = self._new(out)
        return v

    def norm(self, p=2):
        r"""
        Return the norm (or related computations) of the vector.

        INPUT:

        - ``p`` -- (default: 2) controls which norm is computed,
          allowable values are any real number and positive and
          negative infinity.  See output discussion for specifics.

        OUTPUT:

        Returned value is a double precision floating point value
        in ``RDF`` (or an integer when ``p=0``).  The default value
        of ``p = 2`` is the "usual" Euclidean norm.  For other values:

        - ``p = Infinity`` or ``p = oo``: the maximum of the
          absolute values of the entries, where the absolute value
          of the complex number `a+bi` is `\sqrt{a^2+b^2}`
        - ``p = -Infinity`` or ``p = -oo``: the minimum of the
          absolute values of the entries
        - ``p = 0``: the number of nonzero entries in the vector
        - ``p`` any other real number: for a vector `\vec{x}` this method computes

          .. MATH::

                \left(\sum_i x_i^p\right)^{1/p}

          For ``p < 0`` this function is not a norm, but the above
          computation may be useful for other purposes.

        ALGORITHM:

        Computation is performed by the :func:`~scipy:scipy.linalg.norm`
        function of the SciPy/NumPy library.

        EXAMPLES:

        First over the reals.  ::

            sage: v = vector(RDF, range(9))
            sage: v.norm()
            14.28285685...
            sage: v.norm(p=2)
            14.28285685...
            sage: v.norm(p=6)
            8.744039097...
            sage: v.norm(p=Infinity)
            8.0
            sage: v.norm(p=-oo)
            0.0
            sage: v.norm(p=0)
            8.0
            sage: v.norm(p=0.3)
            4099.153615...

        And over the complex numbers.  ::

            sage: w = vector(CDF, [3-4*I, 0, 5+12*I])
            sage: w.norm()
            13.9283882...
            sage: w.norm(p=2)
            13.9283882...
            sage: w.norm(p=0)
            2.0
            sage: w.norm(p=4.2)
            13.0555695...
            sage: w.norm(p=oo)
            13.0

        Negative values of ``p`` are allowed and will
        provide the same computation as for positive values.
        A zero entry in the vector will raise a warning and return
        zero. ::

            sage: v = vector(CDF, range(1,10))
            sage: v.norm(p=-3.2)
            0.953760808...
            sage: w = vector(CDF, [-1,0,1])
            sage: w.norm(p=-1.6)
            doctest:...: RuntimeWarning: divide by zero encountered in power
            0.0

        Return values are in ``RDF``, or an integer when ``p = 0``.  ::

            sage: v = vector(RDF, [1,2,4,8])
            sage: v.norm() in RDF
            True
            sage: v.norm(p=0) in ZZ
            True

        Improper values of ``p`` are caught.  ::

            sage: w = vector(CDF, [-1,0,1])
            sage: w.norm(p='junk')
            Traceback (most recent call last):
            ...
            ValueError: vector norm 'p' must be +/- infinity or a real number, not junk
        """
        import sage.rings.infinity
        import sage.rings.integer
        if p == sage.rings.infinity.Infinity:
            p = numpy.inf
        elif p == -sage.rings.infinity.Infinity:
            p = -numpy.inf
        else:
            try:
                p = RDF(p)
            except Exception:
                raise ValueError("vector norm 'p' must be +/- infinity or a real number, not %s" % p)
        n = numpy.linalg.norm(self._vector_numpy, ord=p)
        # p = 0 returns integer *count* of nonzero entries
        return RDF(n)

    #############################
    # statistics
    #############################

    def mean(self):
        """
        Calculate the arithmetic mean of the vector.

        EXAMPLES::

            sage: v = vector(RDF, range(9))
            sage: w = vector(CDF, [k+(9-k)*I for k in range(9)])
            sage: v.mean()
            4.0
            sage: w.mean()
            4.0 + 5.0*I
        """
        return self._sage_dtype(numpy.mean(self._vector_numpy))

    def variance(self, population=True):
        """
        Calculate the variance of entries of the vector.

        INPUT:

        - ``population`` -- if ``False``, calculate the sample variance

        EXAMPLES::

            sage: v = vector(RDF, range(9))
            sage: w = vector(CDF, [k+(9-k)*I for k in range(9)])
            sage: v.variance()
            7.5
            sage: v.variance(population=False)
            6.666666666666667
            sage: w.variance()
            15.0
            sage: w.variance(population=False)
            13.333333333333334
        """
        if population is True:
            return self._sage_dtype(numpy.var(self._vector_numpy, ddof=1))
        else:
            return self._sage_dtype(numpy.var(self._vector_numpy, ddof=0))

    def standard_deviation(self, population=True):
        """
        Calculate the standard deviation of entries of the vector.

        INPUT:

        - ``population`` -- If ``False``, calculate the sample standard
          deviation

        EXAMPLES::

            sage: v = vector(RDF, range(9))
            sage: w = vector(CDF, [k+(9-k)*I for k in range(9)])
            sage: v.standard_deviation()
            2.7386127875258306
            sage: v.standard_deviation(population=False)
            2.581988897471611
            sage: w.standard_deviation()
            3.872983346207417
            sage: w.standard_deviation(population=False)
            3.6514837167011076
        """
        if population is True:
            return self._sage_dtype(numpy.std(self._vector_numpy, ddof=1))
        else:
            return self._sage_dtype(numpy.std(self._vector_numpy, ddof=0))

    def stats_kurtosis(self):
        """
        Compute the kurtosis of a dataset.

        Kurtosis is the fourth central moment divided by the square of
        the variance. Since we use Fisher's definition, 3.0 is
        subtracted from the result to give 0.0 for a normal
        distribution. (Paragraph from the scipy.stats docstring.)

        EXAMPLES::

            sage: # needs scipy
            sage: v = vector(RDF, range(9))
            sage: w = vector(CDF, [k+(9-k)*I for k in range(9)])
            sage: v.stats_kurtosis()  # rel tol 5e-15
            -1.2300000000000000
            sage: w.stats_kurtosis()  # rel tol 5e-15
            -1.2300000000000000
        """
        import scipy.stats
        return self._sage_dtype(scipy.stats.kurtosis(self._vector_numpy))

    def prod(self):
        """
        Return the product of the entries of ``self``.

        EXAMPLES::

            sage: v = vector(RDF, range(9))
            sage: w = vector(CDF, [k+(9-k)*I for k in range(9)])
            sage: v.prod()
            0.0
            sage: w.prod()
            57204225.0*I
        """
        return self._sage_dtype(self._vector_numpy.prod())

    def sum(self):
        """
        Return the sum of the entries of ``self``.

        EXAMPLES::

            sage: v = vector(RDF, range(9))
            sage: w = vector(CDF, [k+(9-k)*I for k in range(9)])
            sage: v.sum()
            36.0
            sage: w.sum()
            36.0 + 45.0*I
        """
        return self._sage_dtype(self._vector_numpy.sum())
