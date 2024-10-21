# sage_setup: distribution = sagemath-modules
# sage.doctest: optional - numpy
"""
Dense matrices over the Complex Double Field using NumPy

EXAMPLES::

    sage: b = Mat(CDF,2,3).basis()
    sage: b[0,0]
    [1.0 0.0 0.0]
    [0.0 0.0 0.0]

We deal with the case of zero rows or zero columns::

    sage: m = MatrixSpace(CDF,0,3)
    sage: m.zero_matrix()
    []

TESTS::

    sage: # needs sage.symbolic
    sage: a = matrix(CDF,2,[i+(4-i)*I for i in range(4)], sparse=False)
    sage: TestSuite(a).run()
    sage: Mat(CDF,0,0).zero_matrix().inverse()
    []

AUTHORS:

- Jason Grout (2008-09): switch to NumPy backend

- Josh Kantor

- William Stein: many bug fixes and touch ups.
"""

##############################################################################
#       Copyright (C) 2004,2005,2006 Joshua Kantor <kantor.jm@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.rings.complex_double import CDF

cimport numpy as cnumpy

numpy=None

cdef class Matrix_complex_double_dense(Matrix_double_dense):
    """
    Class that implements matrices over the real double field. These
    are supposed to be fast matrix operations using C doubles. Most
    operations are implemented using numpy which will call the
    underlying BLAS on the system.

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: m = Matrix(CDF, [[1,2*I],[3+I,4]])
        sage: m**2
        [-1.0 + 6.0*I       10.0*I]
        [15.0 + 5.0*I 14.0 + 6.0*I]
        sage: n= m^(-1); n  # abs tol 1e-15
        [  0.3333333333333333 + 0.3333333333333333*I 0.16666666666666669 - 0.16666666666666666*I]
        [-0.16666666666666666 - 0.3333333333333333*I 0.08333333333333331 + 0.08333333333333333*I]

    To compute eigenvalues, use the methods
    :meth:`~.Matrix_double_dense.left_eigenvectors` or
    :meth:`~.Matrix_double_dense.right_eigenvectors`::

        sage: p,e = m.right_eigenvectors()                                              # needs sage.symbolic

    The result is a pair ``(p,e)``, where ``p`` is a diagonal matrix of
    eigenvalues and ``e`` is a matrix whose columns are the
    eigenvectors.

    To solve a linear system `Ax = b` where ``A = [[1,2*I],[3+I,4]]`` and
    ``b = [5,6]``::

        sage: b = vector(CDF,[5,6])                                                     # needs sage.symbolic
        sage: m.solve_right(b)  # abs tol 1e-14                                         # needs sage.symbolic
        (2.6666666666666665 + 0.6666666666666669*I, -0.3333333333333333 - 1.1666666666666667*I)

    See the methods :meth:`~.Matrix_double_dense.QR`,
    :meth:`~.Matrix_double_dense.LU`, and :meth:`.SVD` for QR, LU, and singular
    value decomposition.
    """
    def __cinit__(self):
        global numpy
        if numpy is None:
            import numpy
        self._numpy_dtype = numpy.dtype('complex128')
        self._python_dtype = complex
        self._numpy_dtypeint = cnumpy.NPY_CDOUBLE
        # TODO: Make ComplexDoubleElement instead of CDF for speed
        self._sage_dtype = CDF
        self.__create_matrix__()
        return
