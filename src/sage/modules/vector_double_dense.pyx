r"""
Dense vectors using a NumPy backend.

This serves as a base class for dense vectors over Real Double Field and
Complex Double Field

EXAMPLES::

    sage: v = vector(CDF,[(1,-1), (2,pi), (3,5)])
    sage: v
    (1.0 - 1.0*I, 2.0 + 3.141592653589793*I, 3.0 + 5.0*I)
    sage: type(v)
    <type 'sage.modules.vector_complex_double_dense.Vector_complex_double_dense'>
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
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cimport numpy
import numpy
from .free_module_element import FreeModuleElement

from sage.structure.element cimport Element, ModuleElement, RingElement, Vector

from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF

# This is for the NumPy C API (the PyArray... functions) to work
numpy.import_array()


cdef class Vector_double_dense(FreeModuleElement):
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
    def __cinit__(self, parent, entries, coerce=True, copy=True):
        """
        Set up a new vector

        EXAMPLES::

            sage: v = vector(RDF, range(3))
            sage: v.__new__(v.__class__, v.parent(), [1,1,1])  # random output
            (3.00713073107e-261, 3.06320700422e-322, 2.86558074588e-322)
            sage: v = vector(CDF, range(3))
            sage: v.__new__(v.__class__, v.parent(), [1,1,1])  # random output
            (-2.26770549592e-39, 5.1698223615e-252*I, -5.9147262602e-62 + 4.63145528786e-258*I)
        """
        self._is_immutable = 0
        self._degree = parent.degree()
        self._parent = parent

    cdef Vector_double_dense _new(self, numpy.ndarray vector_numpy):
        """
        Return a new vector with same parent as self.
        """
        cdef Vector_double_dense v
        v = self.__class__.__new__(self.__class__,self._parent,None,None,None)
        v._is_immutable = 0
        v._parent = self._parent
        v._degree = self._parent.degree()

        v._vector_numpy = vector_numpy
        return v

    def __create_vector__(self):
        """
        Create a new uninitialized numpy array to hold the data for the class.

        This function assumes that self._numpy_dtypeint and
        self._nrows and self._ncols have already been initialized.

        EXAMPLES:

        In this example, we throw away the current array and make a
        new uninitialized array representing the data for the class. ::

            sage: a=vector(RDF, range(9))
            sage: a.__create_vector__()
        """
        cdef numpy.npy_intp dims[1]
        dims[0] = self._degree
        self._vector_numpy = numpy.PyArray_SimpleNew(1, dims, self._numpy_dtypeint)
        return

    cdef bint is_dense_c(self):
        """
        Return True (i.e., 1) if self is dense.
        """
        return 1

    cdef bint is_sparse_c(self):
        """
        Return True (i.e., 1) if self is sparse.
        """
        return 0


    def __copy__(self, copy=True):
        """
        Return a copy of the vector

        EXAMPLES::

            sage: a = vector(RDF, range(9))
            sage: a == copy(a)
            True
        """
        if self._degree == 0:
            return self
        from copy import copy
        return self._new(copy(self._vector_numpy))

    def __init__(self, parent, entries, coerce = True, copy = True):
        """
        Fill the vector with entries.

        The numpy array must have already been allocated.

        EXAMPLES::

            sage: vector(RDF, range(9))
            (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)
            sage: vector(CDF, 5)
            (0.0, 0.0, 0.0, 0.0, 0.0)

        TESTS::

            sage: vector(CDF, 0)
            ()
            sage: vector(RDF, 0)
            ()
            sage: vector(CDF, 4)
            (0.0, 0.0, 0.0, 0.0)
            sage: vector(RDF, 4)
            (0.0, 0.0, 0.0, 0.0)
            sage: vector(CDF, [CDF(1+I)*j for j in range(4)])
            (0.0, 1.0 + 1.0*I, 2.0 + 2.0*I, 3.0 + 3.0*I)
            sage: vector(RDF, 4, range(4))
            (0.0, 1.0, 2.0, 3.0)

            sage: V = RDF^2
            sage: V.element_class(V, 5)
            Traceback (most recent call last):
            ...
            TypeError: entries must be a list or 0
            sage: V.element_class(V, 0)
            (0.0, 0.0)
        """
        cdef Py_ssize_t i,j
        if isinstance(entries,(tuple, list)):
            if len(entries)!=self._degree:
                    raise TypeError("entries has wrong length")

            if coerce:
                for i from 0<=i<self._degree:
                    self.set_unsafe(i,self._python_dtype(entries[i]))
            else:
                for i from 0<=i<self._degree:
                    self.set_unsafe(i,entries[i])

        elif isinstance(entries, numpy.ndarray):
            self._replace_self_with_numpy(entries)
        else:
            numpy.PyArray_FILLWBYTE(self._vector_numpy, 0)
            if entries is None:
                return
            else:
                try:
                    z = self._python_dtype(entries)
                except TypeError:
                    raise TypeError("unable to convert {!r} to {}".format(entries, self._python_dtype))
                if z != 0:
                    raise TypeError("entries must be a list or 0")
                else:
                    # Set all entries to z=0.
                    for i from 0<=i<self._degree:
                        self.set_unsafe(i,z)

    def __len__(self):
        """
        Return the length of the vector.

        EXAMPLES::

            sage: v = vector(RDF, 5); v
            (0.0, 0.0, 0.0, 0.0, 0.0)
            sage: len(v)
            5
        """
        return self._degree

    cdef int set_unsafe(self, Py_ssize_t i, value) except -1:
        """
        EXAMPLES::

            sage: v = vector(CDF, [1,CDF(3,2), -1]); v
            (1.0, 3.0 + 2.0*I, -1.0)
            sage: v[1] = 2
            sage: v[-1] = I
            sage: v
            (1.0, 2.0, 1.0*I)
            sage: v[1:3] = [1, 1]; v
            (1.0, 1.0, 1.0)
        """
        # We assume that Py_ssize_t is the same as npy_intp

        # We call the self._python_dtype function on the value since
        # numpy does not know how to deal with complex numbers other
        # than the built-in complex number type.
        cdef int status
        status = numpy.PyArray_SETITEM(self._vector_numpy,
                        numpy.PyArray_GETPTR1(self._vector_numpy, i),
                        self._python_dtype(value))
        #TODO: Throw an error if status == -1

    cdef get_unsafe(self, Py_ssize_t i):
        """
        EXAMPLES::

            sage: v = vector(CDF, [1,CDF(3,2), -1]); v
            (1.0, 3.0 + 2.0*I, -1.0)
            sage: v[1]
            3.0 + 2.0*I
            sage: v[-1]
            -1.0
            sage: v[1:3]
            (3.0 + 2.0*I, -1.0)
        """
        # We assume that Py_ssize_t is the same as npy_intp
        return self._sage_dtype(numpy.PyArray_GETITEM(self._vector_numpy,
                                                numpy.PyArray_GETPTR1(self._vector_numpy, i)))


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
        Return self - right

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
        Dot product of self and right.

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
        Return the component-wise product of self and right.

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
        Multiply a scalar and vector

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
        Multiply a scalar and vector

        EXAMPLES::

            sage: v = vector(CDF, range(3))
            sage: v*3
            (0.0, 3.0, 6.0)
        """
        if self._degree == 0:
            from copy import copy
            return copy(self)

        return self._new(self._vector_numpy*self._python_dtype(right))


    def inv_fft(self,algorithm="radix2", inplace=False):
        """
        This performs the inverse fast Fourier transform on the vector.

        The Fourier transform can be done in place using the keyword
        inplace=True

        This will be fastest if the vector's length is a power of 2.

        EXAMPLES::

            sage: v = vector(CDF,[1,2,3,4])
            sage: w = v.fft()
            sage: max(v - w.inv_fft()) < 1e-12
            True
        """
        return self.fft(direction="backward",algorithm=algorithm,inplace=inplace)

    def fft(self, direction = "forward", algorithm = "radix2", inplace=False):
        """
        This performs a fast Fourier transform on the vector.

        INPUT:

        - direction -- 'forward' (default) or 'backward'

        The algorithm and inplace arguments are ignored.

        This function is fastest if the vector's length is a power of 2.

        EXAMPLES::

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
            from .vector_complex_double_dense import Vector_complex_double_dense
            if direction == 'forward':
                return Vector_complex_double_dense(V, fft(self._vector_numpy))
            else:
                return Vector_complex_double_dense(V, ifft(self._vector_numpy))


    cdef _replace_self_with_numpy(self, numpy.ndarray numpy_array):
        """
        Replace the underlying numpy array with numpy_array.
        """
        if self._degree == 0:
            return
        if numpy_array.ndim != 1 or len(self._vector_numpy) != numpy_array.shape[0]:
            raise ValueError("vector lengths are not the same")

        self._vector_numpy = numpy_array.astype(self._numpy_dtype)


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
        Returns a copy with small entries replaced by zeros.

        This is useful for modifying output from algorithms
        which have large relative errors when producing zero
        elements, e.g. to create reliable doctests.

        INPUT:

        - ``eps`` - cutoff value

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
        Returns the norm (or related computations) of the vector.

        INPUT:

        - ``p`` - default: 2 - controls which norm is computed,
          allowable values are any real number and positive and
          negative infinity.  See output discussion for specifics.

        OUTPUT:

        Returned value is a double precision floating point value
        in ``RDF`` (or an integer when ``p=0``).  The default value
        of ``p = 2`` is the "usual" Euclidean norm.  For other values:

        - ``p = Infinity`` or ``p = oo``: the maximum of the
          absolute values of the entries, where the absolute value
          of the complex number `a+bi` is `\sqrt{a^2+b^2}`.
        - ``p = -Infinity`` or ``p = -oo``: the minimum of the
          absolute values of the entries.
        - ``p = 0`` : the number of nonzero entries in the vector.
        - ``p`` is any other real number: for a vector `\vec{x}`
          this method computes

          .. MATH::

                \left(\sum_i x_i^p\right)^{1/p}

          For ``p < 0`` this function is not a norm, but the above
          computation may be useful for other purposes.

        ALGORITHM:

        Computation is performed by the ``norm()`` function of
        the SciPy/NumPy library.

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
        # p = 0 returns integer *count* of non-zero entries
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

        - ``population`` -- If False, calculate the sample variance.

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
            population -- If False, calculate the sample standard deviation.

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
        Return the product of the entries of self.

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
        Return the sum of the entries of self.

        EXAMPLES::

            sage: v = vector(RDF, range(9))
            sage: w = vector(CDF, [k+(9-k)*I for k in range(9)])
            sage: v.sum()
            36.0
            sage: w.sum()
            36.0 + 45.0*I
        """
        return self._sage_dtype(self._vector_numpy.sum())

    # Put this method last, otherwise it overrides the "numpy" cimport
    def numpy(self, dtype=None):
        """
        Return numpy array corresponding to this vector.

        INPUT:

        - ``dtype`` -- if specified, the `numpy dtype <http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html>`_
                       of the returned array.

        EXAMPLES::

            sage: v = vector(CDF,4,range(4))
            sage: v.numpy()
            array([0.+0.j, 1.+0.j, 2.+0.j, 3.+0.j])
            sage: v = vector(CDF,0)
            sage: v.numpy()
            array([], dtype=complex128)
            sage: v = vector(RDF,4,range(4))
            sage: v.numpy()
            array([0., 1., 2., 3.])
            sage: v = vector(RDF,0)
            sage: v.numpy()
            array([], dtype=float64)

        A numpy dtype may be requested manually::

            sage: import numpy
            sage: v = vector(CDF, 3, range(3))
            sage: v.numpy()
            array([0.+0.j, 1.+0.j, 2.+0.j])
            sage: v.numpy(dtype=numpy.float64)
            array([0., 1., 2.])
            sage: v.numpy(dtype=numpy.float32)
            array([0., 1., 2.], dtype=float32)
        """
        if dtype is None or dtype is self._vector_numpy.dtype:
            from copy import copy
            return copy(self._vector_numpy)
        else:
            return super(Vector_double_dense, self).numpy(dtype)
