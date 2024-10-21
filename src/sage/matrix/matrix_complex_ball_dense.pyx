# sage_setup: distribution = sagemath-flint
# distutils: libraries = flint
r"""
Arbitrary precision complex ball matrices

AUTHORS:

- Clemens Heuberger (2014-10-25): Initial version.

This is an incomplete interface to the `acb_mat module
<https://flintlib.org/doc/acb_mat.html>`_ of FLINT; it may be useful to refer
to its documentation for more details.

TESTS::

    sage: mat = matrix(CBF, 2, 2, range(4))
    sage: x = polygen(QQ)
    sage: pol = x^3 + 2
    sage: pol(mat)
    [8.000000000000000 11.00000000000000]
    [22.00000000000000 41.00000000000000]

    sage: mat = matrix(ComplexBallField(20), 2, 2, list(range(4)))*i/3                  # needs sage.symbolic
    sage: loads(dumps(mat)).identical(mat)                                              # needs sage.symbolic
    True
"""
# ****************************************************************************
#       Copyright (C) 2014 Clemens Heuberger <clemens.heuberger@aau.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************
from cpython.object cimport Py_EQ, Py_NE
from cysignals.signals cimport sig_on, sig_str, sig_off

from sage.arith.power cimport generic_power_pos
from sage.libs.flint.acb cimport *
from sage.libs.flint.acb_mat cimport *
from sage.libs.gmp.mpz cimport mpz_fits_ulong_p, mpz_get_ui
from sage.matrix.constructor import matrix
from sage.matrix.args cimport SparseEntry, MatrixArgs_init
from sage.rings.complex_interval cimport ComplexIntervalFieldElement
from sage.rings.complex_arb cimport (
    ComplexBall,
    ComplexIntervalFieldElement_to_acb,
    acb_to_ComplexIntervalFieldElement)
from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_complex_arb cimport Polynomial_complex_arb
from sage.structure.element cimport Element, Matrix
from sage.structure.parent cimport Parent
from sage.structure.sequence import Sequence

from sage.misc.superseded import experimental
from sage.rings.polynomial import polynomial_ring_constructor


cdef void matrix_to_acb_mat(acb_mat_t target, source) noexcept:
    """
    Convert a matrix containing :class:`ComplexIntervalFieldElement` to an ``acb_mat_t``.

    INPUT:

    - ``target`` -- an ``acb_mat_t``

    - ``source`` -- a matrix consisting of :class:`ComplexIntervalFieldElement`

    OUTPUT:

    None.
    """
    cdef unsigned long nrows, ncols, r, c

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

cdef Matrix_generic_dense acb_mat_to_matrix(acb_mat_t source, Parent CIF):
    """
    Convert an ``acb_mat_t`` to a matrix containing :class:`ComplexIntervalFieldElement`.

    INPUT:

    - ``source`` -- an ``acb_mat_t``

    - ``precision`` -- positive integer

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

cdef inline long prec(Matrix_complex_ball_dense mat) noexcept:
    return mat._base_ring._prec

cdef class Matrix_complex_ball_dense(Matrix_dense):
    """
    Matrix over a complex ball field. Implemented using the
    ``acb_mat`` type of the FLINT library.

    EXAMPLES::

        sage: MatrixSpace(CBF, 3)(2)
        [2.000000000000000                 0                 0]
        [                0 2.000000000000000                 0]
        [                0                 0 2.000000000000000]
        sage: matrix(CBF, 1, 3, [1, 2, -3])
        [ 1.000000000000000  2.000000000000000 -3.000000000000000]
    """
    def __cinit__(self):
        """
        Create and allocate memory for the matrix.

        EXAMPLES::

            sage: from sage.matrix.matrix_complex_ball_dense import Matrix_complex_ball_dense
            sage: a = Matrix_complex_ball_dense.__new__( # indirect doctest
            ....:     Matrix_complex_ball_dense, Mat(CBF, 2), 0, 0, 0)
            sage: type(a)
            <class 'sage.matrix.matrix_complex_ball_dense.Matrix_complex_ball_dense'>
        """
        sig_str("FLINT exception")
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

    cdef Matrix_complex_ball_dense _new(self, Py_ssize_t nrows, Py_ssize_t ncols):
        r"""
        Return a new matrix over the same base ring.
        """
        cdef Parent P
        if nrows == self._nrows and ncols == self._ncols:
            P = self._parent
        else:
            P = self.matrix_space(nrows, ncols)
        return Matrix_complex_ball_dense.__new__(Matrix_complex_ball_dense, P, None, None, None)

    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        r"""
        Initialize a dense matrix over the complex ball field.

        INPUT:

        - ``parent`` -- a matrix space over a complex ball field

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- if ``False``, assume without checking that the
          entries lie in the base ring

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
            Complex ball field with 53 bits of precision
            sage: A(range(4))
            [                0 1.000000000000000]
            [2.000000000000000 3.000000000000000]

        Actually it is only necessary that the input can be converted to a
        list, so the following also works::

            sage: v = reversed(range(4)); type(v)
            <...iterator'>
            sage: A(v)
            [3.000000000000000 2.000000000000000]
            [1.000000000000000                 0]

        Matrices can have many rows or columns (in fact, on a 64-bit
        machine they could have up to `2^{63}-1` rows or columns)::

            sage: v = matrix(CBF, 1, 10^5, range(10^5))
            sage: v.parent()
            Full MatrixSpace of 1 by 100000 dense matrices over
            Complex ball field with 53 bits of precision

        TESTS::

            sage: MatrixSpace(CBF, 0, 0).one()
            []
            sage: Matrix(CBF, 0, 100)
            0 x 100 dense matrix over Complex ball field with 53 bits
            of precision (use the '.str()' method to see the entries)
            sage: Matrix(CBF, 100, 0)
            100 x 0 dense matrix over Complex ball field with 53 bits
            of precision (use the '.str()' method to see the entries)
        """
        ma = MatrixArgs_init(parent, entries)
        cdef ComplexBall z
        for t in ma.iter(coerce, True):
            se = <SparseEntry>t
            z = <ComplexBall>se.entry
            acb_set(acb_mat_entry(self.value, se.i, se.j), z.value)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x):
        """
        Set position ``i``, ``j`` of this matrix to ``x``.

        The object ``x`` must be of type ``ComplexBall``.

        INPUT:

        - ``i`` -- row

        - ``j`` -- column

        - ``x`` -- must be ComplexBall! The value to set ``self[i,j]`` to.

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

    cpdef _richcmp_(left, right, int op):
        r"""
        EXAMPLES::

            sage: a = matrix(CBF, [[1,2],[3,4]])
            sage: b = matrix(CBF, [[1,2],[3,4]])
            sage: a == b
            True
            sage: a + 1/3 == b + 1/3
            False
            sage: a < b
            Traceback (most recent call last):
            ...
            TypeError: no order is defined on complex ball matrices

        TESTS::

            sage: a = matrix(CBF, [1/3])
            sage: b = matrix(CBF, [1/3])
            sage: a == a or b == b or a[0,0] == a[0,0] or a[0,0] == b[0,0]
            False
        """
        cdef Matrix_complex_ball_dense lt = <Matrix_complex_ball_dense> left
        cdef Matrix_complex_ball_dense rt = <Matrix_complex_ball_dense> right
        if op == Py_EQ:
            return acb_mat_eq(lt.value, rt.value)
        elif op == Py_NE:
            return acb_mat_ne(lt.value, rt.value)
        else:
            raise TypeError("no order is defined on complex ball matrices")

    def identical(self, Matrix_complex_ball_dense other):
        r"""
        Test if the corresponding entries of two complex ball matrices
        represent the same balls.

        EXAMPLES::

            sage: a = matrix(CBF, [[1/3,2],[3,4]])
            sage: b = matrix(CBF, [[1/3,2],[3,4]])
            sage: a == b
            False
            sage: a.identical(b)
            True
        """
        return acb_mat_equal(self.value, other.value)

    def overlaps(self, Matrix_complex_ball_dense other):
        r"""
        Test if two matrices with complex ball entries represent overlapping
        sets of complex matrices.

        EXAMPLES::

            sage: b = CBF(0, RBF(0, rad=0.1r)); b
            [+/- 0.101]*I
            sage: matrix(CBF, [0, b]).overlaps(matrix(CBF, [b, 0]))
            True
            sage: matrix(CBF, [1, 0]).overlaps(matrix(CBF, [b, 0]))
            False
        """
        return acb_mat_overlaps(self.value, other.value)

    def contains(self, Matrix_complex_ball_dense other):
        r"""
        Test if the set of complex matrices represented by ``self`` is
        contained in that represented by ``other``.

        EXAMPLES::

            sage: b = CBF(0, RBF(0, rad=.1r)); b
            [+/- 0.101]*I
            sage: matrix(CBF, [0, b]).contains(matrix(CBF, [0, 0]))
            True
            sage: matrix(CBF, [0, b]).contains(matrix(CBF, [b, 0]))
            False
            sage: matrix(CBF, [b, b]).contains(matrix(CBF, [b, 0]))
            True
        """
        return acb_mat_contains(self.value, other.value)

    def __neg__(self):
        r"""
        TESTS::

            sage: -matrix(CBF, [[1,2]])
            [-1.000000000000000 -2.000000000000000]
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, self._ncols)
        sig_on()
        acb_mat_neg(res.value, self.value)
        sig_off()
        return res

    cpdef _add_(self, other):
        r"""
        TESTS::

            sage: matrix(CBF, [[1,2]])._add_(matrix(CBF, [3,4]))
            [4.000000000000000 6.000000000000000]
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, self._ncols)
        sig_on()
        acb_mat_add(res.value, self.value, (<Matrix_complex_ball_dense> other).value, prec(self))
        sig_off()
        return res

    cpdef _sub_(self, other):
        r"""
        TESTS::

            sage: matrix(CBF, [[1,2]])._sub_(matrix(CBF, [3,4]))
            [-2.000000000000000 -2.000000000000000]
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, self._ncols)
        sig_on()
        acb_mat_sub(res.value, self.value, (<Matrix_complex_ball_dense> other).value, prec(self))
        sig_off()
        return res

    cpdef _lmul_(self, Element a):
        r"""
        TESTS::

            sage: matrix(CBF, [[1,2]])._lmul_(CBF(I))
            [1.000000000000000*I 2.000000000000000*I]
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, self._ncols)
        sig_on()
        acb_mat_scalar_mul_acb(res.value, self.value, (<ComplexBall> a).value, prec(self))
        sig_off()
        return res

    cpdef _rmul_(self, Element a):
        r"""
        TESTS::

            sage: matrix(CBF, [[1,2]])._rmul_(CBF(I))
            [1.000000000000000*I 2.000000000000000*I]
        """
        return self._lmul_(a)

    cdef _matrix_times_matrix_(self, Matrix other):
        r"""
        TESTS::

            sage: matrix(CBF, [[1,2]])*matrix([[3], [4]]) # indirect doctest
            [11.00000000000000]
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, other._ncols)
        sig_on()
        acb_mat_mul(res.value, self.value, (<Matrix_complex_ball_dense> other).value, prec(self))
        sig_off()
        return res

    cpdef _pow_int(self, n):
        r"""
        Return the ``n``-th power of this matrix.

        EXAMPLES::

            sage: mat = matrix(CBF, [[1/2, 1/3], [1, 1]])
            sage: mat**2
            [[0.5833333333333...] [0.500000000000000 +/- ...e-16]]
            [               1.500000000000000 [1.333333333333333 +/- ...e-16]]
            sage: mat**(-2)
            [ [48.00000000000...] [-18.00000000000...]]
            [[-54.0000000000...]  [21.000000000000...]]

        TESTS::

            sage: mat**(0r)
            [1.000000000000000                 0]
            [                0 1.000000000000000]

            sage: mat**(1/2)
            Traceback (most recent call last):
                ...
            NotImplementedError: non-integral exponents not supported

            sage: (-(matrix(CBF, [2])**(-2**100))[0,0].log(2)).log(2)
            [100.000000000000 +/- ...e-14]
            sage: (-(matrix(CBF, [2])**(-2**64+1))[0,0].log(2)).log(2)
            [64.0000000000000 +/- ...e-14]
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, self._ncols)
        cdef Matrix_complex_ball_dense tmp
        cdef unsigned long expo
        n = Integer(n)
        if self._nrows != self._ncols:
            raise ArithmeticError("self must be a square matrix")

        neg = (n < 0)
        if neg:
            n = -n
        if mpz_fits_ulong_p((<Integer>n).value):
            expo = mpz_get_ui((<Integer>n).value)
            sig_on()
            acb_mat_pow_ui(res.value, self.value, expo, prec(self))
            sig_off()
        else:
            tmp = generic_power_pos(self, n)
            acb_mat_set(res.value, tmp.value)
        if neg:
            sig_on()
            acb_mat_inv(res.value, res.value, prec(self))
            sig_off()

        return res

    def __invert__(self):
        r"""
        TESTS::

            sage: ~matrix(CBF, [[1/2, 1/3], [1, 1]])
            [ [6.00000000000000 +/- ...e-15] [-2.00000000000000 +/- ...e-15]]
            [[-6.00000000000000 +/- ...e-15]  [3.00000000000000 +/- ...e-15]]
            sage: ~matrix(CBF, [[1/2, 1/3]])
            Traceback (most recent call last):
            ...
            ArithmeticError: self must be a square matrix
            sage: mat = matrix(CBF, [[1/3, 1/2], [0, 1]]) - 1/3
            sage: ~mat
            Traceback (most recent call last):
            ...
            ZeroDivisionError: unable to compute the inverse, is the matrix singular?
        """
        if not self.is_square():
            raise ArithmeticError("self must be a square matrix")
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, self._ncols)
        sig_on()
        cdef bint success = acb_mat_inv(res.value, self.value, prec(self))
        sig_off()
        if success:
            return res
        else:
            raise ZeroDivisionError("unable to compute the inverse, is the matrix singular?")

    def transpose(self):
        r"""
        Return the transpose of ``self``.

        EXAMPLES::

            sage: m = matrix(CBF, 2, 3, [1, 2, 3, 4, 5, 6])
            sage: m.transpose()
            [1.000000000000000 4.000000000000000]
            [2.000000000000000 5.000000000000000]
            [3.000000000000000 6.000000000000000]
            sage: m.transpose().parent()
            Full MatrixSpace of 3 by 2 dense matrices over Complex ball field with 53 bits of precision
        """
        cdef Py_ssize_t nc = self._ncols
        cdef Py_ssize_t nr = self._nrows
        cdef Matrix_complex_ball_dense trans = self._new(nc, nr)
        acb_mat_transpose(trans.value, self.value)
        return trans

    def _solve_right_nonsingular_square(self, Matrix_complex_ball_dense rhs, check_rank=None):
        r"""
        TESTS::

            sage: matrix(CBF, [[1/2, 1/3], [1, 1]]).solve_right(vector([-1, 1]))
            ([-8.00000000000000 +/- ...], [9.00000000000000 +/- ...])
            sage: matrix(CBF, 2, 2, 0).solve_right(vector([-1, 1]))
            Traceback (most recent call last):
            ...
            ValueError: unable to invert this matrix
            sage: b = CBF(0, RBF(0, rad=.1r))
            sage: matrix(CBF, [[1, 1], [0, b]]).solve_right(vector([-1, 1]))
            Traceback (most recent call last):
            ...
            ValueError: unable to invert this matrix
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, rhs._ncols)
        sig_on()
        success = acb_mat_solve(res.value, self.value, rhs.value, min(prec(self), prec(rhs)))
        sig_off()
        if success:
            return res
        else:
            raise ValueError("unable to invert this matrix")

    def determinant(self):
        r"""
        Compute the determinant of this matrix.

        EXAMPLES::

            sage: matrix(CBF, [[1/2, 1/3], [1, 1]]).determinant()
            [0.1666666666666667 +/- ...e-17]
            sage: matrix(CBF, [[1/2, 1/3], [1, 1]]).det()
            [0.1666666666666667 +/- ...e-17]
            sage: matrix(CBF, [[1/2, 1/3]]).determinant()
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix
        """
        cdef ComplexBall res = ComplexBall.__new__(ComplexBall)
        res._parent = self._base_ring
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        sig_on()
        acb_mat_det(res.value, self.value, prec(self))
        sig_off()
        return res

    def trace(self):
        r"""
        Compute the trace of this matrix.

        EXAMPLES::

            sage: matrix(CBF, [[1/3, 1/3], [1, 1]]).trace()
            [1.333333333333333 +/- ...e-16]
            sage: matrix(CBF, [[1/2, 1/3]]).trace()
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix
        """
        cdef ComplexBall res = ComplexBall.__new__(ComplexBall)
        res._parent = self._base_ring
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        sig_on()
        acb_mat_trace(res.value, self.value, prec(self))
        sig_off()
        return res

    def charpoly(self, var='x', algorithm=None):
        r"""
        Compute the characteristic polynomial of this matrix.

        EXAMPLES::

            sage: from sage.matrix.benchmark import hilbert_matrix
            sage: mat = hilbert_matrix(5).change_ring(ComplexBallField(10))
            sage: mat.charpoly()
            x^5 + ([-1.8 +/- 0.0258])*x^4 + ([0.3 +/- 0.05...)*x^3 +
            ([+/- 0.0...])*x^2 + ([+/- 0.0...])*x + [+/- 0.0...]

        TESTS::

            sage: mat.charpoly(algorithm='hessenberg')
            x^5 + ([-1.8 +/- 0.04...])*x^4 + ([0.3 +/- 0.08...])*x^3
            + ([+/- 0.0...])*x^2 + ([+/- ...e-4])*x + [+/- ...e-6]
            sage: mat.charpoly('y')
            y^5 + ([-1.8 +/- 0.02...])*y^4 + ([0.3 +/- 0.05...])*y^3 +
            ([+/- 0.0...])*y^2 + ([+/- 0.0...])*y + [+/- 0.0...]
        """
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        if algorithm is not None:
            return super(Matrix_dense, self).charpoly(var=var, algorithm=algorithm)
        Pol = polynomial_ring_constructor._single_variate(self.base_ring(), var)
        cdef Polynomial_complex_arb res = Polynomial_complex_arb(Pol)
        sig_on()
        acb_mat_charpoly(res._poly, self.value, prec(self))
        sig_off()
        return res

    @experimental(issue_number=30393)
    def eigenvalues(self, other=None, *, extend=None):
        r"""
        (Experimental.) Compute rigorous enclosures of the eigenvalues of this matrix.

        INPUT:

        - ``self`` -- an `n \times n` matrix
        - ``other`` -- unsupported (generalized eigenvalue problem), should be ``None``
        - ``extend`` -- ignored

        OUTPUT:

        A :class:`~sage.structure.sequence.Sequence` of complex balls of
        length equal to the size of the matrix.

        Each element represents one eigenvalue with the correct multiplicities
        in case of overlap. The output intervals are either disjoint or
        identical, and identical intervals are guaranteed to be grouped
        consecutively. Each complete run of `k` identical balls thus represents
        a cluster of exactly `k` eigenvalues which could not be separated from
        each other at the current precision, but which could be isolated from
        the other eigenvalues.

        There is currently no guarantee that the algorithm converges as the
        working precision is increased.

        See the `FLINT documentation <https://flintlib.org/doc/acb_mat.html#c.acb_mat_eig_multiple>`__
        for more information.

        EXAMPLES::

            sage: from sage.matrix.benchmark import hilbert_matrix
            sage: mat = hilbert_matrix(5).change_ring(CBF)
            sage: mat.eigenvalues()
            doctest:...: FutureWarning: This class/method/function is marked as experimental.
            ...
            [[1.567050691098...] + [+/- ...]*I, [0.208534218611...] + [+/- ...]*I,
            [3.287928...e-6...] + [+/- ...]*I, [0.000305898040...] + [+/- ...]*I,
            [0.011407491623...] + [+/- ...]*I]

            sage: mat = Permutation([2, 1, 4, 5, 3]).to_matrix().dense_matrix().change_ring(CBF)
            sage: mat.eigenvalues()
            Traceback (most recent call last):
            ...
            ValueError: unable to certify the eigenvalues
            sage: precond = matrix(ZZ, [[-1, -2, 2, 2, -2], [2, -2, -2, -2, 2],
            ....:     [-2, 2, -1, 2, 1], [2, 1, -1, 0, 2], [-2, 0, 1, -1, 1]])
            sage: (~precond*mat*precond).eigenvalues()
            [[-0.5000000000000...] + [-0.8660254037844...]*I, [-1.000000000000...] + [+/- ...]*I,
             [-0.5000000000000...] + [0.8660254037844...]*I,
             [1.000000000000...] + [+/- ...]*I, [1.000000000000...] + [+/- ...]*I]

        .. SEEALSO:: :meth:`eigenvectors_right`
        """
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        cdef long n = self._ncols
        cdef acb_ptr eigval_approx, eigval
        cdef acb_mat_t eigvec_approx
        if other is not None:
            raise NotImplementedError
        try:
            eigval_approx = _acb_vec_init(n)
            acb_mat_init(eigvec_approx, n, n)
            acb_mat_approx_eig_qr(eigval_approx, NULL, eigvec_approx, self.value, NULL, 0, prec(self))
            eigval = _acb_vec_init(n)
            if not acb_mat_eig_multiple(eigval, self.value, eigval_approx, eigvec_approx, prec(self)):
                raise ValueError("unable to certify the eigenvalues")
            res = _acb_vec_to_list(eigval, n, self._parent._base)
        finally:
            acb_mat_clear(eigvec_approx)
            _acb_vec_clear(eigval, n)
            _acb_vec_clear(eigval_approx, n)
        return Sequence(res)

    @experimental(issue_number=30393)
    def eigenvectors_right_approx(self, other=None, *, extend=None):
        r"""
        (Experimental.) Compute *non-rigorous* approximations of the
        eigenvalues and eigenvectors of this matrix.

        INPUT:

        - ``self`` -- an `n \times n` matrix
        - ``other`` -- unsupported (generalized eigenvalue problem), should be ``None``
        - ``extend`` -- ignored

        OUTPUT:

        A list of triples of the form ``(eigenvalue, [eigenvector], 1)``. The
        eigenvalue and the entries of the eigenvector are complex balls with
        zero radius.

        No guarantees are made about the accuracy of the output.

        See the `FLINT documentation <https://flintlib.org/doc/acb_mat.html#c.acb_mat_approx_eig_qr>`__
        for more information.

        EXAMPLES::

            sage: from sage.matrix.benchmark import hilbert_matrix
            sage: mat = hilbert_matrix(3).change_ring(CBF)
            sage: eigval, eigvec, _ = mat.eigenvectors_right_approx()[0]
            doctest:...: FutureWarning: This class/method/function is marked as experimental.
            ...
            sage: eigval
            [1.40831892712...]
            sage: eigval.rad()
            0.00000000
            sage: eigvec
            [([0.8270449269720...], [0.4598639043655...], [0.3232984352444...])]
            sage: (mat - eigval)*eigvec[0]
            ([1e-15 +/- ...], [2e-15 +/- ...], [+/- ...])

        .. SEEALSO:: :meth:`eigenvectors_right`
        """
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        cdef long n = self._ncols
        cdef Matrix_complex_ball_dense eigvec = self._new(n, n)
        cdef acb_ptr _eigval
        if other is not None:
            raise NotImplementedError
        try:
            _eigval = _acb_vec_init(n)
            acb_mat_approx_eig_qr(_eigval, NULL, eigvec.value, self.value, NULL, 0, prec(self))
            eigval = _acb_vec_to_list(_eigval, n, self._parent._base)
        finally:
            _acb_vec_clear(_eigval, n)
        return [(val, [vec], 1) for val, vec in zip(eigval, eigvec.columns())]

    @experimental(issue_number=30393)
    def eigenvectors_right(self, other=None, *, extend=None):
        r"""
        (Experimental.) Compute rigorous enclosures of the eigenvalues and
        eigenvectors of this matrix.

        INPUT:

        - ``self`` -- an `n \times n` matrix
        - ``other`` -- unsupported (generalized eigenvalue problem), should be ``None``
        - ``extend`` -- ignored

        OUTPUT:

        A list of triples of the form ``(eigenvalue, [eigenvector], 1)``.

        Unlike :meth:`eigenvalues` and :meth:`eigenvectors_right_approx`, this
        method currently fails in the presence of multiple eigenvalues.

        Additionally, there is currently no guarantee that the algorithm
        converges as the working precision is increased.

        See the `FLINT documentation <https://flintlib.org/doc/acb_mat.html#c.acb_mat_eig_simple>`__
        for more information.

        EXAMPLES::

            sage: from sage.matrix.benchmark import hilbert_matrix
            sage: mat = hilbert_matrix(3).change_ring(CBF)
            sage: eigval, eigvec, _ = mat.eigenvectors_right()[0]
            doctest:...: FutureWarning: This class/method/function is marked as experimental.
            ...
            sage: eigval
            [1.40831892712...] + [+/- ...]*I
            sage: eigvec
            [([0.82704492697...] + [+/- ...]*I, [0.45986390436...] + [+/- ...]*I, [0.32329843524...] + [+/- ...]*I)]
            sage: (mat - eigval)*eigvec[0]
            ([+/- ...] + [+/- ...]*I, [+/- ...] + [+/- ...]*I, [+/- ...] + [+/- ...]*I)

        .. SEEALSO:: :meth:`eigenvectors_right_approx`, :meth:`eigenvalues`
        """
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        cdef long n = self._ncols
        cdef acb_ptr eigval_approx, _eigval
        cdef acb_mat_t eigvec_approx
        cdef Matrix_complex_ball_dense eigvec = self._new(n, n)
        if other is not None:
            raise NotImplementedError
        try:
            _eigval = _acb_vec_init(n)
            eigval_approx = _acb_vec_init(n)
            acb_mat_init(eigvec_approx, n, n)
            acb_mat_approx_eig_qr(eigval_approx, NULL, eigvec_approx, self.value, NULL, 0, prec(self))
            if not acb_mat_eig_simple(_eigval, NULL, eigvec.value, self.value, eigval_approx, eigvec_approx, prec(self)):
                raise ValueError("unable to isolate the eigenvalues (multiple eigenvalues?)")
            eigval = _acb_vec_to_list(_eigval, n, self._parent._base)
        finally:
            acb_mat_clear(eigvec_approx)
            _acb_vec_clear(_eigval, n)
            _acb_vec_clear(eigval_approx, n)
        return [(val, [vec], 1) for val, vec in zip(eigval, eigvec.columns())]

    def eigenvectors_left_approx(self, other=None, *, extend=None):
        r"""
        (Experimental.) Compute *non-rigorous* approximations of the
        left eigenvalues and eigenvectors of this matrix.

        INPUT:

        - ``self`` -- an `n \times n` matrix
        - ``other`` -- unsupported (generalized eigenvalue problem), should be ``None``
        - ``extend`` -- ignored

        OUTPUT:

        A list of triples of the form ``(eigenvalue, [eigenvector], 1)``. The
        eigenvalue and the entries of the eigenvector are complex balls with
        zero radius.

        No guarantees are made about the accuracy of the output.

        See the `FLINT documentation <https://flintlib.org/doc/acb_mat.html#c.acb_mat_approx_eig_qr>`__
        for more information.

        EXAMPLES::

            sage: mat = matrix(CBF, 3, [2, 3, 5, 7, 11, 13, 17, 19, 23])
            sage: eigval, eigvec, _ = mat.eigenvectors_left_approx()[0]
            sage: eigval
            [1.1052996349... +/- ...]
            sage: eigvec[0]
            ([0.69817246751...], [-0.67419514369...], [0.240865343781...])
            sage: eigvec[0] * (mat - eigval)
            ([+/- ...], [+/- ...], [+/- ...])

        .. SEEALSO:: :meth:`eigenvectors_left`
        """
        return self.transpose().eigenvectors_right_approx(other=None, extend=extend)

    def eigenvectors_left(self, other=None, *, extend=True):
        r"""
        (Experimental.) Compute rigorous enclosures of the eigenvalues and
        left eigenvectors of this matrix.

        INPUT:

        - ``self`` -- an `n \times n` matrix
        - ``other`` -- unsupported (generalized eigenvalue problem), should be ``None``
        - ``extend`` -- ignored

        OUTPUT:

        A list of triples of the form ``(eigenvalue, [eigenvector], 1)``.

        Unlike :meth:`eigenvalues` and :meth:`eigenvectors_left_approx`, this
        method currently fails in the presence of multiple eigenvalues.

        Additionally, there is currently no guarantee that the algorithm
        converges as the working precision is increased.

        See the `FLINT documentation <https://flintlib.org/doc/acb_mat.html#c.acb_mat_eig_simple>`__
        for more information.

        EXAMPLES::

            sage: mat = matrix(CBF, 3, [2, 3, 5, 7, 11, 13, 17, 19, 23])
            sage: eigval, eigvec, _ = mat.eigenvectors_left()[0]
            sage: eigval
            [1.1052996349...] + [+/- ...]*I
            sage: eigvec[0]
            ([0.69817246751...] + [+/- ...]*I, [-0.67419514369...] + [+/- ...]*I, [0.240865343781...] + [+/- ...]*I)
            sage: eigvec[0] * (mat - eigval)
            ([+/- ...] + [+/- ...]*I, [+/- ...] + [+/- ...]*I, [+/- ...] + [+/- ...]*I)

        .. SEEALSO:: :meth:`eigenvectors_right`, :meth:`eigenvalues`, :meth:`eigenvectors_left_approx`
        """
        return self.transpose().eigenvectors_right(other=other, extend=extend)

    def exp(self):
        r"""
        Compute the exponential of this matrix.

        EXAMPLES::

            sage: matrix(CBF, [[i*pi, 1], [0, i*pi]]).exp()                             # needs sage.symbolic
            [[-1.00000000000000 +/- ...e-16] + [+/- ...e-16]*I [-1.00000000000000 +/- ...e-16] + [+/- ...e-16]*I]
            [                                                0 [-1.00000000000000 +/- ...e-16] + [+/- ...e-16]*I]
            sage: matrix(CBF, [[1/2, 1/3]]).exp()
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, self._ncols)
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        sig_on()
        acb_mat_exp(res.value, self.value, prec(self))
        sig_off()
        return res

cdef _acb_vec_to_list(acb_ptr vec, long n, Parent parent):
    cdef ComplexBall b
    res = []
    for i in range(n):
        b = ComplexBall.__new__(ComplexBall)
        b._parent = parent
        acb_set(b.value, &vec[i])
        res.append(b)
    return res
