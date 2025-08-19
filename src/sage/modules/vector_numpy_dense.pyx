# sage.doctest: optional - numpy
r"""
Dense vectors using a NumPy backend.

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
from sage.modules.free_module_element import FreeModuleElement

# This is for the NumPy C API (the PyArray... functions) to work
numpy.import_array()


cdef class Vector_numpy_dense(FreeModuleElement):
    """
    Base class for vectors implemented using numpy arrays.

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
        Set up a new vector.

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

    cdef Vector_numpy_dense _new(self, numpy.ndarray vector_numpy):
        """
        Return a new vector with same parent as ``self``.
        """
        cdef Vector_numpy_dense v
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

    cdef bint is_dense_c(self) noexcept:
        """
        Return ``True`` (i.e., 1) if ``self`` is dense.
        """
        return 1

    cdef bint is_sparse_c(self) noexcept:
        """
        Return ``True`` (i.e., 1) if ``self`` is sparse.
        """
        return 0

    def __copy__(self, copy=True):
        """
        Return a copy of the vector.

        EXAMPLES::

            sage: a = vector(RDF, range(9))
            sage: a == copy(a)
            True
        """
        if self._degree == 0:
            return self
        from copy import copy
        return self._new(copy(self._vector_numpy))

    def __init__(self, parent, entries, coerce=True, copy=True):
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
            sage: vector(CDF, [CDF(1+I)*j for j in range(4)])                           # needs sage.symbolic
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
        cdef Py_ssize_t i
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

            sage: v = vector(CDF, [1, CDF(3,2), -1]); v
            (1.0, 3.0 + 2.0*I, -1.0)
            sage: v[1] = 2
            sage: v[-1] = I                                                             # needs sage.symbolic
            sage: v                                                                     # needs sage.symbolic
            (1.0, 2.0, 1.0*I)
            sage: v[1:3] = [1, 1]; v                                                    # needs sage.symbolic
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
        # TODO: Throw an error if status == -1

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

    cdef _replace_self_with_numpy(self, numpy.ndarray numpy_array):
        """
        Replace the underlying numpy array with numpy_array.
        """
        if self._degree == 0:
            return
        if numpy_array.ndim != 1 or len(self._vector_numpy) != numpy_array.shape[0]:
            raise ValueError("vector lengths are not the same")

        self._vector_numpy = numpy_array.astype(self._numpy_dtype)

    # Put this method last, otherwise it overrides the "numpy" cimport
    def numpy(self, dtype=None):
        """
        Return numpy array corresponding to this vector.

        INPUT:

        - ``dtype`` -- if specified, the `numpy dtype
          <http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html>`_ of
          the returned array

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
            return super().numpy(dtype)
