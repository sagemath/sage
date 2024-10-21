# sage_setup: distribution = sagemath-modules
# sage.doctest: optional - numpy
"""
Dense matrices using a NumPy backend

AUTHORS:

- Jason Grout, Sep 2008: switch to NumPy backend, factored out the Matrix_double_dense class

- Josh Kantor

- William Stein: many bug fixes and touch ups.
"""

# ****************************************************************************
#       Copyright (C) 2004-2006 Joshua Kantor <kantor.jm@gmail.com>
#       Copyright (C) 2008      Georg S. Weber
#       Copyright (C) 2008-2011 Mike Hansen
#       Copyright (C) 2008-2012 Jason Grout
#       Copyright (C) 2009      Dag Sverre Seljebotn
#       Copyright (C) 2009      Yann Laigle-Chapuy
#       Copyright (C) 2009-2010 Florent Hivert
#       Copyright (C) 2010-2012 Rob Beezer
#       Copyright (C) 2011      Martin Raum
#       Copyright (C) 2011-2012 J. H. Palmieri
#       Copyright (C) 2011-2014 André Apitzsch
#       Copyright (C) 2011-2018 Jeroen Demeyer
#       Copyright (C) 2012      Kenneth Smith
#       Copyright (C) 2016-2019 Frédéric Chapoton
#       Copyright (C) 2017      Kiran Kedlaya
#       Copyright (C) 2019      Chaman Agrawal
#       Copyright (C) 2019-2021 Markus Wageringel
#       Copyright (C) 2020      Michael Orlitzky
#       Copyright (C) 2020      Victor Santos
#       Copyright (C) 2021      Jonathan Kliem
#       Copyright (C) 2021      Travis Scrimshaw
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.args cimport MatrixArgs_init
cimport sage.structure.element

cimport numpy as cnumpy

numpy = None
scipy = None

# This is for the Numpy C API to work
cnumpy.import_array()


cdef class Matrix_numpy_dense(Matrix_dense):

    def __create_matrix__(self):
        """
        Create a new uninitialized numpy matrix to hold the data for the class.

        This function assumes that self._numpy_dtypeint and
        self._nrows and self._ncols have already been initialized.

        EXAMPLES:

        In this example, we throw away the current matrix and make a
        new uninitialized matrix representing the data for the class::

            sage: a = matrix(RDF, 3, range(9))
            sage: a.__create_matrix__()
        """
        cdef cnumpy.npy_intp dims[2]
        dims[0] = self._nrows
        dims[1] = self._ncols
        self._matrix_numpy = cnumpy.PyArray_SimpleNew(2, dims, self._numpy_dtypeint)
        return

    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        r"""
        Fill the matrix with entries.

        The numpy matrix must have already been allocated.

        INPUT:

        - ``parent`` -- a matrix space over ``RDF``

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- if ``True`` (the default), convert elements to the
          base ring before passing them to NumPy. If ``False``, pass the
          elements to NumPy as given.

        EXAMPLES::

            sage: matrix(RDF,3,range(9))
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]
            [6.0 7.0 8.0]
            sage: matrix(CDF,3,3,2)
            [2.0 0.0 0.0]
            [0.0 2.0 0.0]
            [0.0 0.0 2.0]

        TESTS::

            sage: matrix(RDF,3,0)
            []
            sage: matrix(RDF,3,3,0)
            [0.0 0.0 0.0]
            [0.0 0.0 0.0]
            [0.0 0.0 0.0]
            sage: matrix(RDF,3,3,1)
            [1.0 0.0 0.0]
            [0.0 1.0 0.0]
            [0.0 0.0 1.0]
            sage: matrix(RDF,3,3,2)
            [2.0 0.0 0.0]
            [0.0 2.0 0.0]
            [0.0 0.0 2.0]
            sage: matrix(CDF,3,0)
            []
            sage: matrix(CDF,3,3,0)
            [0.0 0.0 0.0]
            [0.0 0.0 0.0]
            [0.0 0.0 0.0]
            sage: matrix(CDF,3,3,1)
            [1.0 0.0 0.0]
            [0.0 1.0 0.0]
            [0.0 0.0 1.0]
            sage: matrix(CDF,3,3,range(9))
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]
            [6.0 7.0 8.0]
            sage: matrix(CDF,2,2,[CDF(1+I)*j for j in range(4)])                        # needs sage.symbolic
            [        0.0 1.0 + 1.0*I]
            [2.0 + 2.0*I 3.0 + 3.0*I]
        """
        ma = MatrixArgs_init(parent, entries)
        cdef long i, j
        it = ma.iter(coerce)
        for i in range(ma.nrows):
            for j in range(ma.ncols):
                self.set_unsafe(i, j, next(it))

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object value):
        """
        Set the (i,j) entry to value without any bounds checking,
        mutability checking, etc.
        """
        # We assume that Py_ssize_t is the same as cnumpy.npy_intp

        # We must patch the ndarrayobject.h file so that the SETITEM
        # macro does not have a semicolon at the end for this to work.
        # Cython wraps the macro in a function that converts the
        # returned int to a python object, which leads to compilation
        # errors because after preprocessing you get something that
        # looks like "););".  This is bug
        # http://scipy.org/scipy/numpy/ticket/918

        # We call the self._python_dtype function on the value since
        # numpy does not know how to deal with complex numbers other
        # than the built-in complex number type.
        cdef int status
        status = cnumpy.PyArray_SETITEM(self._matrix_numpy,
                        cnumpy.PyArray_GETPTR2(self._matrix_numpy, i, j),
                        self._python_dtype(value))
        #TODO: Throw an error if status == -1

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Get the (i,j) entry without any bounds checking, etc.
        """
        # We assume that Py_ssize_t is the same as cnumpy.npy_intp
        return self._sage_dtype(cnumpy.PyArray_GETITEM(self._matrix_numpy,
                                                cnumpy.PyArray_GETPTR2(self._matrix_numpy, i, j)))

    cdef Matrix_numpy_dense _new(self, int nrows=-1, int ncols=-1):
        """
        Return a new uninitialized matrix with same parent as ``self``.

        INPUT:

        - ``nrows`` -- (default: ``self._nrows``) number of rows in returned matrix
        - ``ncols`` -- (default: ``self._ncols``) number of columns in returned matrix
        """
        cdef Matrix_numpy_dense m
        if nrows == -1 and ncols == -1:
            nrows = self._nrows
            ncols = self._ncols
            parent = self._parent
        else:
            if nrows == -1:
                nrows = self._nrows
            if ncols == -1:
                ncols = self._ncols
            parent = self.matrix_space(nrows, ncols)
        m = self.__class__.__new__(self.__class__, parent, None, None, None)
        return m

    def __copy__(self):
        r"""
        Return a new copy of this matrix.

        EXAMPLES::

            sage: a = matrix(RDF,1,3, [1,2,-3])
            sage: a
            [ 1.0  2.0 -3.0]
            sage: b = a.__copy__()
            sage: b
            [ 1.0  2.0 -3.0]
            sage: b is a
            False
            sage: b == a
            True
            sage: b[0,0] = 3
            sage: a[0,0] # note that a hasn't changed
            1.0

        ::

            sage: copy(MatrixSpace(RDF,0,0,sparse=False).zero_matrix())
            []
        """
        if self._nrows == 0 or self._ncols == 0:
            # Create a brand new empty matrix. This is needed to prevent a
            # recursive loop: a copy of zero_matrix is asked otherwise.
            return self.__class__(self.parent(), [], self._nrows, self._ncols)

        cdef Matrix_numpy_dense A
        A = self._new(self._nrows, self._ncols)
        A._matrix_numpy = self._matrix_numpy.copy()
        if self._subdivisions is not None:
            A.subdivide(*self.subdivisions())
        return A

    def transpose(self):
        """
        Return the transpose of this matrix, without changing ``self``.

        EXAMPLES::

            sage: m = matrix(RDF,2,3,range(6)); m
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]
            sage: m2 = m.transpose()
            sage: m[0,0] = 2
            sage: m2           #note that m2 hasn't changed
            [0.0 3.0]
            [1.0 4.0]
            [2.0 5.0]

        ``.T`` is a convenient shortcut for the transpose::

            sage: m.T
            [2.0 3.0]
            [1.0 4.0]
            [2.0 5.0]

            sage: m = matrix(RDF,0,3); m
            []
            sage: m.transpose()
            []
            sage: m.transpose().parent()
            Full MatrixSpace of 3 by 0 dense matrices over Real Double Field
        """
        if self._nrows == 0 or self._ncols == 0:
            return self.new_matrix(self._ncols, self._nrows)

        cdef Matrix_numpy_dense trans
        trans = self._new(self._ncols, self._nrows)
        trans._matrix_numpy = self._matrix_numpy.transpose().copy()
        if self._subdivisions is not None:
            row_divs, col_divs = self.subdivisions()
            trans.subdivide(col_divs, row_divs)
        return trans

    def is_symmetric(self, tol=1e-12):
        """
        Return whether this matrix is symmetric, to the given tolerance.

        EXAMPLES::

            sage: m = matrix(RDF,2,2,range(4)); m
            [0.0 1.0]
            [2.0 3.0]
            sage: m.is_symmetric()
            False
            sage: m[1,0] = 1.1; m
            [0.0 1.0]
            [1.1 3.0]
            sage: m.is_symmetric()
            False

        The tolerance inequality is strict::

            sage: m.is_symmetric(tol=0.1)
            False
            sage: m.is_symmetric(tol=0.11)
            True

        TESTS:

        Complex entries are supported (:issue:`27831`).  ::

            sage: a = matrix(CDF, [(21, 0.6 + 18.5*i), (0.6 - 18.5*i, 21)])             # needs sage.symbolic
            sage: a.is_symmetric()                                                      # needs sage.symbolic
            False
        """
        cdef Py_ssize_t i, j
        tol = float(tol)
        key = 'symmetric_%s' % tol
        b = self.fetch(key)
        if b is not None:
            return b
        if self._nrows != self._ncols:
            self.cache(key, False)
            return False
        b = True
        for i from 0 < i < self._nrows:
            for j from 0 <= j < i:
                if abs(self.get_unsafe(i,j) - self.get_unsafe(j,i)) > tol:
                    b = False
                    break
        self.cache(key, b)
        return b

    def _is_lower_triangular(self, tol):
        r"""
        Return ``True`` if the entries above the diagonal are all zero.

        INPUT:

        - ``tol`` -- the largest value of the absolute value of the
          difference between two matrix entries for which they will
          still be considered equal.

        OUTPUT:

        Return ``True`` if each entry above the diagonal (entries
        with a row index less than the column index) is zero.

        EXAMPLES::

            sage: A = matrix(RDF, [[ 2.0, 0.0,  0.0],
            ....:                  [ 1.0, 3.0,  0.0],
            ....:                  [-4.0, 2.0, -1.0]])
            sage: A._is_lower_triangular(1.0e-17)
            True
            sage: A[1,2] = 10^-13
            sage: A._is_lower_triangular(1.0e-14)
            False
        """
        global numpy
        if numpy is None:
            import numpy
        cdef Py_ssize_t i, j
        for i in range(self._nrows):
            for j in range(i+1, self._ncols):
                if abs(self.get_unsafe(i,j)) > tol:
                    return False
        return True

    def numpy(self, dtype=None):
        """
        Return a copy of the matrix as a numpy array.

        It uses the numpy C/api so is very fast.

        INPUT:

        - ``dtype`` -- the desired data-type for the array. If not given,
          then the type will be determined as the minimum type required
          to hold the objects in the sequence.

        EXAMPLES::

            sage: m = matrix(RDF,[[1,2],[3,4]])
            sage: n = m.numpy()
            sage: import numpy
            sage: tuple(numpy.linalg.eig(n))
            (array([-0.37228132,  5.37228132]),
             array([[-0.82456484, -0.41597356],
                   [ 0.56576746, -0.90937671]]))
            sage: m = matrix(RDF, 2, range(6)); m
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]
            sage: m.numpy()
            array([[0., 1., 2.],
                   [3., 4., 5.]])

        Alternatively, numpy automatically calls this function (via
        the magic :meth:`__array__` method) to convert Sage matrices
        to numpy arrays::

            sage: import numpy
            sage: m = matrix(RDF, 2, range(6)); m
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]
            sage: numpy.array(m)
            array([[0., 1., 2.],
                   [3., 4., 5.]])
            sage: numpy.array(m).dtype
            dtype('float64')
            sage: m = matrix(CDF, 2, range(6)); m
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]
            sage: numpy.array(m)
            array([[0.+0.j, 1.+0.j, 2.+0.j],
                   [3.+0.j, 4.+0.j, 5.+0.j]])
            sage: numpy.array(m).dtype
            dtype('complex128')

        TESTS::

            sage: m = matrix(RDF,0,5,[]); m
            []
            sage: m.numpy()
            array([], shape=(0, 5), dtype=float64)
            sage: m = matrix(RDF,5,0,[]); m
            []
            sage: m.numpy()
            array([], shape=(5, 0), dtype=float64)
        """
        import numpy as np
        if dtype is None or self._numpy_dtype == np.dtype(dtype):
            return self._matrix_numpy.copy()
        else:
            return Matrix_dense.numpy(self, dtype=dtype)

    def _replace_self_with_numpy(self, numpy_matrix):
        """

        EXAMPLES::

            sage: import numpy
            sage: a = numpy.array([[1,2],[3,4]], 'float64')
            sage: m = matrix(RDF,2,2,0)
            sage: m._replace_self_with_numpy(a)
            sage: m
            [1.0 2.0]
            [3.0 4.0]
        """
        if (<object>self._matrix_numpy).shape != (<object>numpy_matrix).shape:
            raise ValueError("matrix shapes are not the same")
        self._matrix_numpy = numpy_matrix.astype(self._numpy_dtype)
