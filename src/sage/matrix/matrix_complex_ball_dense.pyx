r"""
Arbitrary precision complex ball matrices using Arb

AUTHORS:

- Clemens Heuberger (2014-10-25): Initial version.

This is a rudimentary binding to the `Arb library
<http://fredrikj.net/arb/>`_; it may be useful to refer to its
documentation for more details.

TESTS::

    sage: mat = matrix(CBF, 2, 2, range(4))
    sage: x = polygen(QQ)
    sage: pol = x^3 + 2
    sage: pol(mat)
    [8.000000000000000 11.00000000000000]
    [22.00000000000000 41.00000000000000]
"""
#*****************************************************************************
# Copyright (C) 2014 Clemens Heuberger <clemens.heuberger@aau.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                http://www.gnu.org/licenses/
#*****************************************************************************

include 'cysignals/signals.pxi'


from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_mat cimport *
from sage.matrix.matrix cimport Matrix
from sage.matrix.constructor import matrix
from sage.matrix.matrix_generic_sparse cimport Matrix_generic_sparse
from sage.rings.complex_interval_field import ComplexIntervalField_class, ComplexIntervalField
from sage.rings.complex_interval cimport ComplexIntervalFieldElement
from sage.rings.complex_arb cimport (
    ComplexBall,
    ComplexIntervalFieldElement_to_acb,
    acb_to_ComplexIntervalFieldElement)
from sage.structure.element cimport Element


cdef void matrix_to_acb_mat(acb_mat_t target, source):
    """
    Convert a matrix containing :class:`ComplexIntervalFieldElement` to an ``acb_mat_t``.

    INPUT:

    - ``target`` -- an ``acb_mat_t``

    - ``source`` -- a matrix consisting of :class:`ComplexIntervalFieldElement`

    OUTPUT:

    None.
    """
    cdef unsigned long nrows, ncols, r, c, precision

    nrows = acb_mat_nrows(target)
    ncols = acb_mat_ncols(target)

    for r in range(nrows):
        for c in range(ncols):
            ComplexIntervalFieldElement_to_acb(acb_mat_entry(target, r, c),
                                               source[r][c])

cdef ComplexIntervalFieldElement _to_CIF(acb_t source, ComplexIntervalFieldElement template):
    cdef ComplexIntervalFieldElement result
    result = template._new()
    acb_to_ComplexIntervalFieldElement(
        result, source)
    return result

cdef Matrix_generic_dense acb_mat_to_matrix(
    acb_mat_t source, Parent CIF):
    """
    Convert an ``acb_mat_t`` to a matrix containing :class:`ComplexIntervalFieldElement`.

    INPUT:

    - ``source`` -- an ``acb_mat_t``

    - ``precision`` -- a positive integer.

    OUTPUT:

    A :class:`~sage.matrix.matrix_generic_dense.Matrix_generic_dense`
    containing :class:`ComplexIntervalFieldElement`.
    """
    cdef unsigned long nrows, ncols, r, c
    cdef ComplexIntervalFieldElement template

    nrows = acb_mat_nrows(source)
    ncols = acb_mat_ncols(source)
    template = CIF(0)

    return matrix(
                  [[_to_CIF(acb_mat_entry(source, r, c), template)
                    for c in range(ncols)]
                   for r in range(nrows)])


cdef class Matrix_complex_ball_dense(matrix_dense.Matrix_dense):
    """
    Matrix over a complex ball field. Implemented using the
    ``acb_mat`` type of the Arb library.

    EXAMPLES::

        sage: MatrixSpace(CBF, 3)(2)
        [2.000000000000000                 0                 0]
        [                0 2.000000000000000                 0]
        [                0                 0 2.000000000000000]
        sage: matrix(CBF, 1, 3, [1, 2, -3])
        [ 1.000000000000000  2.000000000000000 -3.000000000000000]
    """
    #################################################################
    # LEVEL 1 functionality
    # * __cinit__
    # * __init__
    # * __dealloc__
    # * set_unsafe(self, size_t i, size_t j, x)
    # * get_unsafe(self, size_t i, size_t j)
    ################################################################

    def __cinit__(self,
                  parent,
                  entries,
                  coerce,
                  copy):
        """
        Create and allocate memory for the matrix.

        INPUT:

        -  ``parent, entries, coerce, copy`` - as for
           ``__init__``.

        EXAMPLES::

            sage: from sage.matrix.matrix_complex_ball_dense import Matrix_complex_ball_dense
            sage: a = Matrix_complex_ball_dense.__new__( # indirect doctest
            ....:     Matrix_complex_ball_dense, Mat(CBF, 2), 0, 0, 0)
            sage: type(a)
            <type 'sage.matrix.matrix_complex_ball_dense.Matrix_complex_ball_dense'>
        """
        self._parent = parent
        self._base_ring = parent.base_ring()
        self._nrows = parent.nrows()
        self._ncols = parent.ncols()
        sig_str("Arb exception")
        acb_mat_init(self.value, self._nrows, self._ncols)
        sig_off()

    def __dealloc__(self):
        """
        Free all the memory allocated for this matrix.

        EXAMPLES::

            sage: a = Matrix(CBF, 2, [1, 2, 3, 4]) # indirect doctest
            sage: del a
        """
        acb_mat_clear(self.value)


    def __init__(self,
                 parent,
                 entries,
                 copy,
                 coerce):
        r"""
        Initialize a dense matrix over the complex ball field.

        INPUT:

        -  ``parent`` - a matrix space

        -  ``entries`` - list - create the matrix with those
           entries along the rows.

        -  ``other`` - a scalar; entries is coerced to a complex ball
           and the diagonal entries of this matrix are set to that
           complex ball.

        -  ``coerce`` - whether need to coerce entries to the
           complex ball field (program may crash if you get this wrong)

        -  ``copy`` - ignored (since complex balls are immutable)

        EXAMPLES:

        The ``__init__`` function is called implicitly in each of the
        examples below to actually fill in the values of the matrix.

        We create a `2 \times 2` and a `1\times 4` matrix::

            sage: matrix(CBF, 2, 2, range(4))
            [                0 1.000000000000000]
            [2.000000000000000 3.000000000000000]
            sage: Matrix(CBF, 1, 4, range(4))
            [                0 1.000000000000000 2.000000000000000 3.000000000000000]

        If the number of columns isn't given, it is determined from the
        number of elements in the list. ::

            sage: matrix(CBF, 2, range(4))
            [                0 1.000000000000000]
            [2.000000000000000 3.000000000000000]
            sage: matrix(CBF, 2, range(6))
            [                0 1.000000000000000 2.000000000000000]
            [3.000000000000000 4.000000000000000 5.000000000000000]

        Another way to make a matrix is to create the space of matrices and
        convert lists into it. ::

            sage: A = Mat(CBF, 2); A
            Full MatrixSpace of 2 by 2 dense matrices over
            Complex ball field with 53 bits precision
            sage: A(range(4))
            [                0 1.000000000000000]
            [2.000000000000000 3.000000000000000]

        Actually it is only necessary that the input can be converted to a
        list, so the following also works::

            sage: v = reversed(range(4)); type(v)
            <... 'listreverseiterator'>
            sage: A(v)
            [3.000000000000000 2.000000000000000]
            [1.000000000000000                 0]

        Matrices can have many rows or columns (in fact, on a 64-bit
        machine they could have up to `2^{63}-1` rows or columns)::

            sage: v = matrix(CBF, 1, 10^5, range(10^5))
            sage: v.parent()
            Full MatrixSpace of 1 by 100000 dense matrices over
            Complex ball field with 53 bits precision

        TESTS::

            sage: MatrixSpace(CBF, 0, 0).one()
            []
            sage: Matrix(CBF, 0, 100)
            0 x 100 dense matrix over Complex ball field with 53 bits precision
            (use the '.str()' method to see the entries)
            sage: Matrix(CBF, 100, 0)
            100 x 0 dense matrix over Complex ball field with 53 bits precision
            (use the '.str()' method to see the entries)
        """
        cdef Py_ssize_t i, j, k
        cdef bint is_list
        cdef ComplexBall x

        if entries is None:
            x = self._base_ring.zero()
            is_list = False
        elif isinstance(entries, (int, long, Element)):
            try:
                x = self._base_ring(entries)
            except TypeError:
                raise TypeError("unable to convert entry to a complex ball")
            is_list = False
        else:
            entries = list(entries)
            is_list = True

        if is_list:
            # Create the matrix whose entries are in the given entry list.
            if len(entries) != self._nrows * self._ncols:
                raise TypeError("entries has the wrong length")
            if coerce:
                k = 0
                for i in range(self._nrows):
                    for j in range(self._ncols):
                        x = self._base_ring(entries[k])
                        acb_set(acb_mat_entry(self.value, i, j),
                                x.value)
                        k += 1
            else:
                k = 0
                for i in range(self._nrows):
                    for j in range(self._ncols):
                        acb_set(acb_mat_entry(self.value, i, j),
                                (<ComplexBall> entries[k]).value)
                        k += 1
        else:
            # If x is zero, make the zero matrix and be done.
            if acb_is_zero(x.value):
                acb_mat_zero(self.value)
                return

            # the matrix must be square:
            if self._nrows != self._ncols:
                raise TypeError("nonzero scalar matrix must be square")

            # Now we set all the diagonal entries to x and all other entries to 0.
            acb_mat_zero(self.value)
            for i in range(self._nrows):
                acb_set(acb_mat_entry(self.value, i, i), x.value)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x):
        """
        Set position ``i``, ``j`` of this matrix to ``x``.

        The object ``x`` must be of type ``ComplexBall``.

        INPUT:

        - ``i`` -- row

        - ``j`` -- column

        - ``x`` -- must be ComplexBall! The value to set self[i,j] to.

        EXAMPLES::

            sage: a = matrix(CBF, 2, 3, range(6)); a
            [                0 1.000000000000000 2.000000000000000]
            [3.000000000000000 4.000000000000000 5.000000000000000]
            sage: a[0, 0] = 10
            sage: a
            [10.00000000000000  1.000000000000000  2.000000000000000]
            [3.000000000000000  4.000000000000000  5.000000000000000]
        """
        acb_set(acb_mat_entry(self.value, i, j), (<ComplexBall> x).value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Return ``(i, j)`` entry of this matrix as a new ComplexBall.

        .. warning::

           This is very unsafe; it assumes ``i`` and ``j`` are in the right
           range.

        EXAMPLES::

            sage: a = MatrixSpace(CBF, 3)(range(9)); a
            [                0 1.000000000000000 2.000000000000000]
            [3.000000000000000 4.000000000000000 5.000000000000000]
            [6.000000000000000 7.000000000000000 8.000000000000000]
            sage: a[1, 2]
            5.000000000000000
            sage: a[4, 7]
            Traceback (most recent call last):
            ...
            IndexError: matrix index out of range
            sage: a[-1, 0]
            6.000000000000000
        """
        cdef ComplexBall z = ComplexBall.__new__(ComplexBall)
        z._parent = self._base_ring
        acb_set(z.value, acb_mat_entry(self.value, i, j))
        return z
