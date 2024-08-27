"""
Dense matrices over univariate polynomials over fields

The implementation inherits from Matrix_generic_dense but some algorithms
are optimized for polynomial matrices.

AUTHORS:

- Kwankyu Lee (2016-12-15): initial version with code moved from other files.

- Johan Rosenkilde (2017-02-07): added weak_popov_form()

- Vincent Neiger (2018-06-13): added basic functions (row/column degrees,
  leading positions, leading matrix, testing reduced and canonical forms)

- Vincent Neiger (2018-09-29): added functions for computing and for verifying
  minimal approximant bases

- Vincent Neiger (2020-04-01): added functions for computing and for verifying
  minimal kernel bases

- Vincent Neiger (2021-03-11): added matrix-wise basic functions for univariate
  polynomials (shifts, reverse, truncate, get coefficient of specified degree)

- Vincent Neiger (2021-07-29): added popov_form(). Added more options to
  weak_popov_form() (column-wise, ordered, zero rows).

- Vincent Neiger (2021-08-07): added inverse_series_trunc(),
  solve_{left/right}_series_trunc(), {left/right}_quo_rem(), reduce().

- Vincent Neiger (2024-02-13): added basis_completion(), _is_basis_completion(),
  _basis_completion_via_reversed_approx().
"""
# ****************************************************************************
#       Copyright (C) 2016 Kwankyu Lee <ekwankyu@gmail.com>
#       Copyright (C) 2017 Johan Rosenkilde
#       Copyright (C) 2018,2020,2021,2024 Vincent Neiger
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.matrix.matrix2 cimport Matrix
from sage.rings.integer_ring import ZZ


cdef class Matrix_polynomial_dense(Matrix_generic_dense):
    r"""
    Dense matrix over a univariate polynomial ring over a field.

    For a field `\Bold{K}`, we consider matrices over the univariate
    polynomial ring `\Bold{K}[x]`.

    They are often used to represent bases of some `\Bold{K}[x]`-modules. In
    this context, there are two possible representations which are both
    commonly used in the literature.

    - Working column-wise: each column of the matrix is a vector in the basis;
      then, a `\Bold{K}[x]`-submodule of `\Bold{K}[x]^{m}` of rank `n` is
      represented by an `m \times n` matrix, whose columns span the module
      (via `\Bold{K}[x]`-linear combinations). This matrix has full rank,
      and `n \leq m`.

    - Working row-wise: each row of the matrix is a vector in the basis; then,
      a `\Bold{K}[x]`-submodule of `\Bold{K}[x]^{n}` of rank `m` is
      represented by an `m \times n` matrix, whose rows span the module (via
      `\Bold{K}[x]`-linear combinations). This matrix has full rank, and `m
      \leq n`.

    For the rest of this class description, we assume that one is working
    row-wise. For a given such module, all its bases are equivalent under
    left-multiplication by a unimodular matrix, that is, a square matrix which
    has determinant in `\Bold{K}\setminus\{0\}`.

    There are bases which are called reduced or minimal: their rows have the
    minimal degree possible among all bases of this module; here the degree of
    a row is the maximum of the degrees of the entries of the row. An
    equivalent condition is that the leading matrix of this basis has full rank
    (see :meth:`leading_matrix`, :meth:`reduced_form`, :meth:`is_reduced`).
    There is a unique minimal basis, called the Popov basis of the module,
    which satisfies some additional normalization condition (see
    :meth:`popov_form`, :meth:`is_popov`).

    These notions can be extended via a more general degree measure, involving
    a tuple of integers which is called shift and acts as column degree shifts
    in the definition of row degree. Precisely, for given `s_1,\ldots,s_n \in
    \ZZ` and a row vector `[p_1 \; \cdots \; p_n] \in \Bold{K}[x]^{1 \times
    n}`, its shifted row degree is the maximum of `\deg(p_j) + s_j` for `1 \leq
    j \leq n` (see :meth:`row_degrees`). Then, reduced bases and Popov bases
    are defined similarly, with respect to this notion of degree.

    Another important canonical basis is the Hermite basis, which is an upper
    triangular matrix satisfying a normalization condition similar to that for
    the Popov basis. In fact, if `d` is the largest degree appearing in the
    Hermite basis, then the Hermite basis coincide with the shifted Popov basis
    with the shifts `((n-1)d,\ldots,2d,d,0)`.
    """

    def _check_shift_dimension(self, shifts, row_wise=True):
        r"""
        Raises an exception if the ``shifts`` argument does not have the right
        length.

        For an `m \times n` polynomial matrix, if working row-wise then
        ``shifts`` should have `n` entries; if working column-wise, it should
        have `m` entries.

        INPUT:

        - ``shifts`` -- list of integers, or ``None``

        - ``row_wise`` -- boolean (default: ``True``); if ``True``
          then shifts apply to the columns of the matrix and otherwise to its
          rows (see the class description for more details)

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [[3*x+1, 0, 1], [x^3+3, 0, 0]])
            sage: M._check_shift_dimension(shifts=[1,3,2])

            sage: M._check_shift_dimension(shifts=[1,3,2], row_wise=False)
            Traceback (most recent call last):
            ...
            ValueError: shifts length should be the row dimension
        """
        if shifts is not None and (not row_wise) and len(shifts) != self.nrows():
            raise ValueError('shifts length should be the row dimension')
        if shifts is not None and (row_wise and len(shifts) != self.ncols()):
            raise ValueError('shifts length should be the column dimension')

    def degree(self):
        r"""
        Return the degree of this matrix.

        For a given polynomial matrix, its degree is the maximum of the degrees
        of all its entries. If the matrix is nonzero, this is a nonnegative
        integer; here, the degree of the zero matrix is -1.

        OUTPUT: integer

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [[3*x+1, 0, 1], [x^3+3, 0, 0]])
            sage: M.degree()
            3

        The zero matrix has degree ``-1``::

            sage: M = Matrix(pR, 2, 3)
            sage: M.degree()
            -1

        For an empty matrix, the degree is not defined::

            sage: M = Matrix(pR, 3, 0)
            sage: M.degree()
            Traceback (most recent call last):
            ...
            ValueError: empty matrix does not have a degree
        """
        if self.nrows() == 0 or self.ncols() == 0:
            raise ValueError('empty matrix does not have a degree')
        return max(self[i, j].degree()
                   for i in range(self.nrows()) for j in range(self.ncols()))

    def degree_matrix(self, shifts=None, row_wise=True):
        r"""
        Return the matrix of the (shifted) degrees in this matrix.

        For a given polynomial matrix `M = (M_{i,j})_{i,j}`, its degree matrix
        is the matrix `(\deg(M_{i,j}))_{i,j}` formed by the degrees of its
        entries. Here, the degree of the zero polynomial is `-1`.

        For given shifts `s_1,\ldots,s_m \in \ZZ`, the shifted degree
        matrix of `M` is either `(\deg(M_{i,j})+s_j)_{i,j}` if working
        row-wise, or `(\deg(M_{i,j})+s_i)_{i,j}` if working column-wise. In the
        former case, `m` has to be the number of columns of `M`; in the latter
        case, the number of its rows. Here, if `M_{i,j}=0` then the
        corresponding entry in the shifted degree matrix is
        `\min(s_1,\ldots,s_m)-1`. For more on shifts and working row-wise
        versus column-wise, see the class documentation.

        INPUT:

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        - ``row_wise`` -- boolean (default: ``True``); if ``True``
          then shifts apply to the columns of the matrix and otherwise to its
          rows (see the class description for more details)

        OUTPUT: integer matrix

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [[3*x+1, 0, 1], [x^3+3, 0, 0]])
            sage: M.degree_matrix()
            [ 1 -1  0]
            [ 3 -1 -1]

            sage: M.degree_matrix(shifts=[0,1,2])
            [ 1 -1  2]
            [ 3 -1 -1]

        The zero entries in the polynomial matrix can be identified in the
        (shifted) degree matrix as the entries equal to ``min(shifts)-1``::

            sage: M.degree_matrix(shifts=[-2,1,2])
            [-1 -3  2]
            [ 1 -3 -3]

        Using ``row_wise=False``, the function supports shifts applied to the
        rows of the matrix (which, in terms of modules, means that we are
        working column-wise, see the class documentation)::

            sage: M.degree_matrix(shifts=[-1,2], row_wise=False)
            [ 0 -2 -1]
            [ 5 -2 -2]
        """
        self._check_shift_dimension(shifts,row_wise)
        if shifts is None:
            return self.apply_map(lambda x: x.degree())
        from sage.matrix.constructor import matrix
        zero_degree = min(shifts) - 1
        if row_wise:
            return matrix( ZZ, [[ self[i,j].degree() + shifts[j]
                if self[i,j] != 0 else zero_degree
                for j in range(self.ncols()) ] for i in range(self.nrows())] )
        else:
            return matrix( ZZ, [[ self[i,j].degree() + shifts[i]
                if self[i,j] != 0 else zero_degree
                for j in range(self.ncols()) ] for i in range(self.nrows())] )

    def constant_matrix(self):
        r"""
        Return the constant coefficient of this matrix seen as a polynomial
        with matrix coefficients; this is also this matrix evaluated at zero.

        OUTPUT: a matrix over the base field

        EXAMPLES::

            sage: pR.<x> = GF(7)[]

            sage: M = Matrix([
            ....:    [  x^3+5*x^2+5*x+1,       5,       6*x+4,         0],
            ....:    [      6*x^2+3*x+1,       1,           2,         0],
            ....:    [2*x^3+4*x^2+6*x+4, 5*x + 1, 2*x^2+5*x+5, x^2+5*x+6]
            ....:     ])
            sage: M.constant_matrix()
            [1 5 4 0]
            [1 1 2 0]
            [4 1 5 6]
        """
        from sage.matrix.constructor import matrix
        return matrix([[self[i,j].constant_coefficient()
            for j in range(self.ncols())] for i in range(self.nrows())])

    def is_constant(self):
        r"""
        Return whether this polynomial matrix is constant,
        that is, all its entries are constant.

        OUTPUT: boolean

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix([
            ....:    [  x^3+5*x^2+5*x+1,       5,       6*x+4,         0],
            ....:    [      6*x^2+3*x+1,       1,           2,         0],
            ....:    [2*x^3+4*x^2+6*x+4, 5*x + 1, 2*x^2+5*x+5, x^2+5*x+6]
            ....:     ])
            sage: M.is_constant()
            False
            sage: M = Matrix(pR, [[1,5,2], [3,1,5]]); M.is_constant()
            True
            sage: M = Matrix.zero(pR, 3, 5); M.is_constant()
            True

        .. SEEALSO::

            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.is_constant`
        """
        return all([self[i,j].is_constant()
            for j in range(self.ncols()) for i in range(self.nrows())])

    def coefficient_matrix(self,d,row_wise=True):
        r"""
        Return the constant matrix which is obtained from this matrix by taking
        the coefficient of its entries with degree specified by `d`.

        - if `d` is an integer, this selects the coefficient of `d` for all
          entries;
        - if `d` is a list `(d_1,\ldots,d_m)` and ``row_wise`` is ``True``,
          this selects the coefficient of degree `d_i` for all entries of the
          `i`-th row for each `i`;
        - if `d` is a list `(d_1,\ldots,d_n)` and ``row_wise`` is ``False``,
          this selects the coefficient of degree `d_i` for all entries of the
          `j`-th column for each `j`.

        INPUT:

        - ``d`` -- list of integers, or an integer,

        - ``row_wise`` -- boolean (default: ``True``); if ``True``
          (resp. ``False``) then `d` should be a list of length equal to the
          row (resp. column) dimension of this matrix

        OUTPUT: a matrix over the base field

        EXAMPLES::

            sage: pR.<x> = GF(7)[]

            sage: M = Matrix([
            ....:    [  x^3+5*x^2+5*x+1,       5,       6*x+4,         0],
            ....:    [      6*x^2+3*x+1,       1,           2,         0],
            ....:    [2*x^3+4*x^2+6*x+4, 5*x + 1, 2*x^2+5*x+5, x^2+5*x+6]
            ....:     ])
            sage: M.coefficient_matrix(2)
            [5 0 0 0]
            [6 0 0 0]
            [4 0 2 1]
            sage: M.coefficient_matrix(0) == M.constant_matrix()
            True

        Row-wise and column-wise coefficient extraction are available::

            sage: M.coefficient_matrix([3,2,1])
            [1 0 0 0]
            [6 0 0 0]
            [6 5 5 5]

            sage: M.coefficient_matrix([2,0,1,3], row_wise=False)
            [5 5 6 0]
            [6 1 0 0]
            [4 1 5 0]

        Negative degrees give zero coefficients::

            sage: M.coefficient_matrix([-1,0,1,3], row_wise=False)
            [0 5 6 0]
            [0 1 0 0]
            [0 1 5 0]

        Length of list of degrees is checked::

            sage: M.coefficient_matrix([2,1,1,2])
            Traceback (most recent call last):
            ...
            ValueError: length of input degree list should be the row
            dimension of the input matrix

            sage: M.coefficient_matrix([3,2,1], row_wise=False)
            Traceback (most recent call last):
            ...
            ValueError: length of input degree list should be the column
            dimension of the input matrix
        """
        m = self.nrows()
        n = self.ncols()

        # if d is an integer, make it a uniform list
        if not isinstance(d,list):
            d = [d]*m if row_wise else [d]*n

        # raise an error if d does not have the right length
        if row_wise and len(d) != m:
            raise ValueError("length of input degree list should be the "
                             "row dimension of the input matrix")
        elif (not row_wise) and len(d) != n:
            raise ValueError("length of input degree list should be the "
                             "column dimension of the input matrix")

        from sage.matrix.constructor import matrix
        return matrix(self.base_ring().base_ring(), m, n,
                [[self[i,j][d[i]] if row_wise else self[i,j][d[j]]
            for j in range(n)] for i in range(m)])

    def truncate(self, d, row_wise=True):
        r"""
        Return the matrix which is obtained from this matrix after truncating
        all its entries according to precisions specified by `d`.

        - if `d` is an integer, the truncation is at precision `d` for all
          entries;
        - if `d` is a list `(d_1,\ldots,d_m)` and ``row_wise`` is ``True``, all
          entries of the `i`-th row are truncated at precision `d_i` for each
          `i`;
        - if `d` is a list `(d_1,\ldots,d_n)` and ``row_wise`` is ``False``,
          all entries of the `j`-th column are truncated at precision `d_j` for
          each `j`.

        Here the convention for univariate polynomials is to take zero
        for the truncation for a negative `d`.

        INPUT:

        - ``d`` -- list of integers, or an integer,

        - ``row_wise`` -- boolean (default: ``True``); if ``True``
          (resp. ``False``) then `d` should be a list of length equal to the
          row (resp. column) dimension of this matrix

        OUTPUT: a polynomial matrix

        EXAMPLES::

            sage: pR.<x> = GF(7)[]

            sage: M = Matrix([
            ....:    [  x^3+5*x^2+5*x+1,       5,       6*x+4,         0],
            ....:    [      6*x^2+3*x+1,       1,           2,         0],
            ....:    [2*x^3+4*x^2+6*x+4, 5*x + 1, 2*x^2+5*x+5, x^2+5*x+6]
            ....:     ])
            sage: M.truncate(2)
            [5*x + 1       5 6*x + 4       0]
            [3*x + 1       1       2       0]
            [6*x + 4 5*x + 1 5*x + 5 5*x + 6]
            sage: M.truncate(1) == M.constant_matrix()
            True

        Row-wise and column-wise truncation are available::

            sage: M.truncate([3,2,1])
            [5*x^2 + 5*x + 1               5         6*x + 4               0]
            [        3*x + 1               1               2               0]
            [              4               1               5               6]

            sage: M.truncate([2,1,1,2], row_wise=False)
            [5*x + 1       5       4       0]
            [3*x + 1       1       2       0]
            [6*x + 4       1       5 5*x + 6]

        Length of list of truncation orders is checked::

            sage: M.truncate([2,1,1,2])
            Traceback (most recent call last):
            ...
            ValueError: length of input precision list should be the row
            dimension of the input matrix

            sage: M.truncate([3,2,1], row_wise=False)
            Traceback (most recent call last):
            ...
            ValueError: length of input precision list should be the column
            dimension of the input matrix

        .. SEEALSO::

            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.truncate` .
        """
        m = self.nrows()
        n = self.ncols()
        from sage.matrix.constructor import matrix

        # if d is an integer, make it a uniform list
        if not isinstance(d,list):
            d = [d]*m if row_wise else [d]*n

        # raise an error if d does not have the right length
        if row_wise and len(d) != m:
            raise ValueError("length of input precision list should be the "
                             "row dimension of the input matrix")
        elif (not row_wise) and len(d) != n:
            raise ValueError("length of input precision list should be the "
                             "column dimension of the input matrix")

        return matrix(self.base_ring(), m, n, [[self[i,j].truncate(d[i])
            if row_wise else self[i,j].truncate(d[j])
            for j in range(n)] for i in range(m)])

    def shift(self, d, row_wise=True):
        r"""
        Return the matrix which is obtained from this matrix after shifting
        all its entries as specified by `d`.

        - if `d` is an integer, the shift is by `d` for all entries;
        - if `d` is a list `(d_1,\ldots,d_m)` and ``row_wise`` is ``True``, all
          entries of the `i`-th row are shifted by `d_i` for each `i`;
        - if `d` is a list `(d_1,\ldots,d_n)` and ``row_wise`` is ``False``,
          all entries of the `j`-th column are shifted by `d_j` for each `j`.

        Shifting by `d` means multiplying by the variable to the power `d`; if
        `d` is negative then terms of negative degree after shifting are
        discarded.

        INPUT:

        - ``d`` -- list of integers, or an integer,

        - ``row_wise`` -- boolean (default: ``True``); if ``True``
          (resp. ``False``) then `d` should be a list of length equal to the
          row (resp. column) dimension of this matrix

        OUTPUT: a polynomial matrix

        EXAMPLES::

            sage: pR.<x> = GF(7)[]

            sage: M = Matrix([
            ....:    [  x^3+5*x^2+5*x+1,       5,       6*x+4,         0],
            ....:    [      6*x^2+3*x+1,       1,           2,         0],
            ....:    [2*x^3+4*x^2+6*x+4, 5*x + 1, 2*x^2+5*x+5, x^2+5*x+6]
            ....:     ])
            sage: M.shift(-2)
            [  x + 5       0       0       0]
            [      6       0       0       0]
            [2*x + 4       0       2       1]

        Row-wise and column-wise shifting are available::

            sage: M.shift([-1,2,-2])
            [      x^2 + 5*x + 5                   0                   6                   0]
            [6*x^4 + 3*x^3 + x^2                 x^2               2*x^2                   0]
            [            2*x + 4                   0                   2                   1]

            sage: M.shift([-1,1,0,0], row_wise=False)
            [  x^2 + 5*x + 5             5*x         6*x + 4               0]
            [        6*x + 3               x               2               0]
            [2*x^2 + 4*x + 6       5*x^2 + x 2*x^2 + 5*x + 5   x^2 + 5*x + 6]

            sage: M.shift([-d for d in M.row_degrees()]) == M.leading_matrix()
            True

        Length of input shift degree list is checked::

            sage: M.shift([1,3,1,4])
            Traceback (most recent call last):
            ...
            ValueError: length of input shift list should be the row
            dimension of the input matrix

            sage: M.shift([5,2,-1], row_wise=False)
            Traceback (most recent call last):
            ...
            ValueError: length of input shift list should be the column
            dimension of the input matrix

        .. SEEALSO::

            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.shift` .
        """
        m = self.nrows()
        n = self.ncols()
        from sage.matrix.constructor import matrix

        # if d is an integer, make it a uniform list
        if not isinstance(d,list):
            d = [d]*m if row_wise else [d]*n

        # raise an error if d does not have the right length
        if row_wise and len(d) != m:
            raise ValueError("length of input shift list should be the "
                             "row dimension of the input matrix")
        elif (not row_wise) and len(d) != n:
            raise ValueError("length of input shift list should be the "
                             "column dimension of the input matrix")

        return matrix(self.base_ring(), m, n, [[self[i,j].shift(d[i])
            if row_wise else self[i,j].shift(d[j])
            for j in range(n)] for i in range(m)])

    def reverse(self, degree=None, row_wise=True, entry_wise=False):
        r"""
        Return the matrix which is obtained from this matrix after reversing
        all its entries with respect to the degree as specified by ``degree``.

        Reversing a polynomial with respect to an integer `d` follows the
        convention for univariate polynomials, in particular it uses truncation
        or zero-padding as necessary if `d` differs from the degree of this
        polynomial.

        If ``entry_wise`` is ``True``: the input ``degree`` and ``row_wise``
        are ignored, and all entries of the matrix are reversed with respect to
        their respective degrees.

        If ``entry_wise`` is ``False`` (the default):

        - if ``degree`` is an integer, all entries are reversed with respect to
          it;
        - if ``degree`` is not provided, then all entries are reversed with
          respect to the degree of the whole matrix;
        - if ``degree`` is a list `(d_1,\ldots,d_m)` and ``row_wise`` is
          ``True``, all entries of the `i`-th row are reversed with respect to
          `d_i` for each `i`;
        - if ``degree`` is a list `(d_1,\ldots,d_n)` and ``row_wise`` is
          ``False``, all entries of the `j`-th column are reversed with respect
          to `d_j` for each `j`.

        INPUT:

        - ``degree`` -- (default: ``None``) a list of nonnegative
          integers, or a nonnegative integer

        - ``row_wise`` -- boolean (default: ``True``); if ``True``
          (resp. ``False``) then ``degree`` should be a list of length equal to
          the row (resp. column) dimension of this matrix

        - ``entry_wise`` -- boolean (default: ``False``); if ``True``
          then the input ``degree`` and ``row_wise`` are ignored

        OUTPUT: a polynomial matrix

        EXAMPLES::

            sage: pR.<x> = GF(7)[]

            sage: M = Matrix([
            ....:    [  x^3+5*x^2+5*x+1,       5,       6*x+4,         0],
            ....:    [      6*x^2+3*x+1,       1,           2,         0],
            ....:    [2*x^3+4*x^2+6*x+4, 5*x + 1, 2*x^2+5*x+5, x^2+5*x+6]
            ....:     ])
            sage: M.reverse()
            [  x^3 + 5*x^2 + 5*x + 1                   5*x^3           4*x^3 +
            6*x^2                       0]
            [      x^3 + 3*x^2 + 6*x                     x^3
            2*x^3                       0]
            [4*x^3 + 6*x^2 + 4*x + 2             x^3 + 5*x^2     5*x^3 + 5*x^2
            + 2*x       6*x^3 + 5*x^2 + x]

            sage: M.reverse(1)
            [  x + 5     5*x 4*x + 6       0]
            [  x + 3       x     2*x       0]
            [4*x + 6   x + 5 5*x + 5 6*x + 5]

            sage: M.reverse(0) == M.constant_matrix()
            True

        Entry-wise reversing with respect to each entry's degree::

            sage: M.reverse(entry_wise=True)
            [  x^3 + 5*x^2 + 5*x + 1                       5
            4*x + 6                       0]
            [          x^2 + 3*x + 6                       1
            2                       0]
            [4*x^3 + 6*x^2 + 4*x + 2                   x + 5         5*x^2 +
            5*x + 2         6*x^2 + 5*x + 1]

        Row-wise and column-wise degree reversing are available::

            sage: M.reverse([2,3,1])
            [    x^2 + 5*x + 5             5*x^2       4*x^2 + 6*x
            0]
            [x^3 + 3*x^2 + 6*x               x^3             2*x^3
            0]
            [          4*x + 6             x + 5           5*x + 5
            6*x + 5]

            sage: M.reverse(M.column_degrees(), row_wise=False)
            [  x^3 + 5*x^2 + 5*x + 1                     5*x             4*x^2
            + 6*x                       0]
            [      x^3 + 3*x^2 + 6*x                       x
            2*x^2                       0]
            [4*x^3 + 6*x^2 + 4*x + 2                   x + 5         5*x^2 +
            5*x + 2         6*x^2 + 5*x + 1]

        Wrong length or negativity of input degree raise errors:

            sage: M.reverse([1,3,1,4])
            Traceback (most recent call last):
            ...
            ValueError: length of input degree list should be the row
            dimension of the input matrix

            sage: M.reverse([5,2,1], row_wise=False)
            Traceback (most recent call last):
            ...
            ValueError: length of input degree list should be the column
            dimension of the input matrix

            sage: M.reverse([2,3,-1])
            Traceback (most recent call last):
            ...
            ValueError: degree argument must be a nonnegative integer, got -1

        .. SEEALSO::

            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.reverse` .
        """
        m = self.nrows()
        n = self.ncols()
        from sage.matrix.constructor import matrix

        # if entry_wise, just return the matrix with all entries reversed
        if entry_wise:
            return matrix(self.base_ring(), m, n, [[self[i,j].reverse()
                for j in range(n)] for i in range(m)])

        # if degree is None, make it the matrix degree
        if degree is None:
            degree = self.degree()
        # if degree is an integer, make it a uniform list
        if not isinstance(degree,list):
            degree = [degree]*m if row_wise else [degree]*n

        # raise an error if degree does not have the right length
        if row_wise and len(degree) != m:
            raise ValueError("length of input degree list should be the "
                             "row dimension of the input matrix")
        elif (not row_wise) and len(degree) != n:
            raise ValueError("length of input degree list should be the "
                             "column dimension of the input matrix")

        return matrix(self.base_ring(), m, n, [[self[i,j].reverse(degree[i])
            if row_wise else self[i,j].reverse(degree[j])
            for j in range(n)] for i in range(m)])

    def inverse_series_trunc(self, d):
        r"""
        Return a matrix polynomial approximation of precision ``d`` of the
        inverse series of this matrix polynomial.

        Here matrix polynomial means that ``self`` is seen as a univariate
        polynomial with matrix coefficients, meaning that this method has the
        same output as if one: 1) converts this matrix to a univariate
        polynomial with matrix coefficients, 2) calls

        :meth:`sage.rings.polynomial.polynomial_element.Polynomial.inverse_series_trunc`

        on that univariate polynomial, and 3) converts back to a matrix of
        polynomials.

        Raises a :exc:`ZeroDivisionError` if the constant matrix of ``self`` is
        not invertible (i.e. has zero determinant); raises an
        :exc:`ArithmeticError` if ``self`` is nonsquare; and raises a
        :exc:`ValueError` if the precision ``d`` is not positive.

        INPUT:

        - ``d`` -- positive integer

        OUTPUT: the unique polynomial matrix `B` of degree less than `d` such
        that `AB` and `BA` are the identity matrix modulo `x^d`, where `A` is
        ``self``.

        ALGORITHM: This uses Newton iteration, performing about `\log(d)`
        polynomial matrix multiplications in size `m \times m` and in degree
        less than `2d`, where `m` is the row dimension of ``self``.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: A = Matrix(pR, 3, 3,
            ....:            [[4*x+5,           5*x^2 + x + 1, 4*x^2 + 4],
            ....:             [6*x^2 + 6*x + 6, 4*x^2 + 5*x,   4*x^2 + x + 3],
            ....:             [3*x^2 + 2,       4*x + 1,       x^2 + 3*x]])
            sage: B = A.inverse_series_trunc(4); B
            [    x^3 + 5*x^2 + x + 4   x^3 + 5*x^2 + 6*x + 4         6*x^2 + 5*x + 3]
            [        4*x^2 + 5*x + 6     6*x^3 + x^2 + x + 6       3*x^3 + 2*x^2 + 2]
            [5*x^3 + 5*x^2 + 6*x + 6 4*x^3 + 2*x^2 + 6*x + 4   6*x^3 + x^2 + 6*x + 1]
            sage: (B*A).truncate(4) == 1
            True

            sage: A.inverse_series_trunc(0)
            Traceback (most recent call last):
            ...
            ValueError: the precision must be positive

            sage: A[:2,:].inverse_series_trunc(4)
            Traceback (most recent call last):
            ...
            ArithmeticError: the input matrix must be square

            sage: A[0,:] = A[0,:] - A[0,:](0) + A[1,:](0) + A[2,:](0)
            sage: A.inverse_series_trunc(4)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: the constant matrix term self(0) must be invertible

        .. SEEALSO::

            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.inverse_series_trunc` .

        .. TODO::

            in the current state of polynomial matrix multiplication (July
            2021), it would be highly beneficial to use conversions and rely on
            polynomials with matrix coefficients when the matrix size is
            "large" and the degree "small", see
            :issue:`31472#comment:5`.
        """
        if d <= 0:
            raise ValueError("the precision must be positive")
        try:
            inv_trunc = self.constant_matrix().inverse()
        except ZeroDivisionError:
            raise ZeroDivisionError("the constant matrix term self(0)"
                                    " must be invertible")
        except ArithmeticError:
            raise ArithmeticError("the input matrix must be square")

        # in comments below, A=self, B=inv_trunc
        # at this point, B = A^{-1} mod x
        from sage.misc.misc import newton_method_sizes
        for next_prec in newton_method_sizes(d)[1:]:
            # Newton iteration: B = (2I - B*A)*B mod x^next_prec
            mul_trunc = - (inv_trunc*self).truncate(next_prec)
            for i in range(self.nrows()):
                mul_trunc[i,i] += 2
            inv_trunc = (mul_trunc * inv_trunc).truncate(next_prec)
        return inv_trunc

    def solve_left_series_trunc(self, B, d):
        r"""
        Try to find a solution `X` to the equation `X A = B`, at precision
        ``d``.

        If ``self`` is a matrix `A`, then this function returns a vector or
        matrix `X` such that `X A = B \bmod x^d`. If `B` is a vector then `X`
        is a vector, and if `B` is a matrix then `X` is a matrix.

        Raises :exc:`ValueError` if ``d`` is not strictly positive, or if there is
        a dimension mismatch between `A` and `B`, or if there is no solution to
        the given matrix equation at the specified precision.

        INPUT:

        - ``B`` -- a polynomial matrix or polynomial vector

        - ``d`` -- positive integer

        OUTPUT:

        A solution to the matrix equation, returned as polynomial matrix of
        degree less than ``d`` if ``B`` is a polynomial matrix; a polynomial
        vector of degree less than ``d`` if `B` is a polynomial vector.

        ALGORITHM:

        If `A` is square with invertible constant term, then the unique
        solution is found by calling :meth:`inverse_series_trunc` and using the
        Dixon & Moenck-Carter iteration. Otherwise, a left minimal approximant
        basis of a matrix formed by `A` and `B` is computed, for an appropriate
        shift which ensures that this basis reveals either a solution `X` or
        the fact that no such solution exists.

        EXAMPLES::


            sage: pR.<x> = GF(7)[]
            sage: A = Matrix(pR, 3, 3,
            ....:            [[4*x+5,           5*x^2 + x + 1, 4*x^2 + 4],
            ....:             [6*x^2 + 6*x + 6, 4*x^2 + 5*x,   4*x^2 + x + 3],
            ....:             [3*x^2 + 2,       4*x + 1,       x^2 + 3*x]])
            sage: A.is_square() and A.constant_matrix().is_invertible()
            True
            sage: B = vector([2*x^2 + 6*x + 6, 0, x + 6])
            sage: X = A.solve_left_series_trunc(B, 4); X
            (3*x^3 + 3*x^2 + 2*x + 4, 4*x^3 + x^2 + 2*x + 6, 6*x^3 + x + 3)
            sage: B == X*A % x**4
            True

            sage: B = Matrix(pR, 2, 3,
            ....:            [[3*x, x^2 + x + 2, x^2 + 2*x + 3],
            ....:             [  0,   6*x^2 + 1,             1]])
            sage: A.solve_left_series_trunc(B, 3)
            [6*x^2 + 2*x + 2         4*x + 3     2*x^2 + 3*x]
            [3*x^2 + 4*x + 5       4*x^2 + 3   x^2 + 6*x + 3]
            sage: X = A.solve_left_series_trunc(B, 37); B == X*A % x**37
            True

        Dimensions of input are checked::

            sage: A.solve_left_series_trunc(B[:,:2], 3)
            Traceback (most recent call last):
            ...
            ValueError: number of columns of self must equal number of columns of right-hand side

        Raises an exception when no solution::

            sage: A[2:,:].solve_left_series_trunc(B, 4)
            Traceback (most recent call last):
            ...
            ValueError: matrix equation has no solutions

            sage: Ax = x*A; C = vector(pR, [1,1,1])
            sage: Ax.solve_left_series_trunc(C, 5)
            Traceback (most recent call last):
            ...
            ValueError: matrix equation has no solutions

        Supports rectangular and rank-deficient cases::

            sage: A[:,:2].solve_left_series_trunc(B[:,:2], 4)
            [5*x^2 + 2*x + 5         5*x + 5         2*x + 4]
            [5*x^3 + 2*x + 1 2*x^2 + 2*x + 5           4*x^2]

            sage: V = Matrix([[3*x^2 + 4*x + 1, 4*x]])
            sage: A[:2,:].solve_left_series_trunc(V*A[:2,:], 4) == V
            True

            sage: A[1,:] = (x+1) * A[0,:]; A[2,:] = (x+5) * A[0,:]
            sage: B = (3*x^3+x^2+2)*A[0,:]
            sage: A.solve_left_series_trunc(B, 6)
            [4*x^2 + 6*x + 2       3*x^2 + x               0]

        .. SEEALSO::

            :meth:`solve_right_series_trunc` .
        """
        from sage.structure.element import Vector
        if isinstance(B, Vector):
            if self.ncols() != B.degree():
                raise ValueError("number of columns of self must equal "
                                 "degree of right-hand side")
        else:
            if self.ncols() != B.ncols():
                raise ValueError("number of columns of self must equal "
                                 "number of columns of right-hand side")

        if d <= 0:
            raise ValueError("the precision must be positive")
        try:
            # case where self is square, with invertible constant term
            precA = 1+self.degree()
            if isinstance(B, Vector):
                BB = B.row()
                X = B.row().parent().zero().__copy__()
            else:
                BB = B.__copy__()
                X = B.parent().zero().__copy__()
            inv_self = self.inverse_series_trunc(precA)
            for k in range(0,(d/precA).ceil()):
                # compute XX = BB * invA mod x^precA
                XX = (BB * inv_self).truncate(precA)
                # compute BB = (BB - XX*A) x^(-precA)
                BB = (BB - XX*self).shift(-precA)
                # update X = X + x^(k*precA) * XX
                X = X + XX.shift(k*precA)
            return X.truncate(d)[0] if isinstance(B, Vector) else X.truncate(d)
        except (ZeroDivisionError,ArithmeticError):
            # general case (possibly no solution)
            m = self.nrows()
            from sage.matrix.constructor import matrix
            if isinstance(B, Vector):
                F = matrix.block([[self],[-B.row()]])
                s = [0]*m + [d]
            else:
                F = matrix.block([[self],[-B]])
                s = [0]*m + [d]*(B.nrows())
            # Warning: the next call works, with the current implementation
            # (Beckermann-Labahn style) of minimal approximant basis, but this
            # implementation is modified then one might have to require Popov
            # basis with the option "normal_form=True"
            P = F.minimal_approximant_basis(d,s)
            if P[m:,m:] != 1:
                raise ValueError("matrix equation has no solutions")
            else:
                return P[m][:m] if isinstance(B, Vector) else P[m:,:m]

    def solve_right_series_trunc(self, B, d):
        r"""
        Try to find a solution `X` to the equation `A X = B`, at precision
        ``d``.

        If ``self`` is a matrix `A`, then this function returns a vector or
        matrix `X` such that `A X = B \bmod x^d`. If `B` is a vector then `X`
        is a vector, and if `B` is a matrix then `X` is a matrix.

        Raises :exc:`ValueError` if ``d`` is not strictly positive, or if there is
        a dimension mismatch between `A` and `B`, or if there is no solution to
        the given matrix equation at the specified precision.

        INPUT:

        - ``B`` -- a polynomial matrix or polynomial vector

        - ``d`` -- positive integer

        OUTPUT:

        A solution to the matrix equation, returned as polynomial matrix of
        degree less than ``d`` if ``B`` is a polynomial matrix; a polynomial
        vector of degree less than ``d`` if `B` is a polynomial vector.

        ALGORITHM:

        If `A` is square with invertible constant term, then the unique
        solution is found by calling :meth:`inverse_series_trunc` and using the
        Dixon & Moenck-Carter iteration. Otherwise, a right minimal approximant
        basis of a matrix formed by `A` and `B` is computed, for an appropriate
        shift which ensures that this basis reveals either a solution `X` or
        the fact that no such solution exists.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: A = Matrix(pR, 3, 3,
            ....:     [[4*x+5,           5*x^2 + x + 1, 4*x^2 + 4],
            ....:      [6*x^2 + 6*x + 6, 4*x^2 + 5*x,   4*x^2 + x + 3],
            ....:      [3*x^2 + 2,       4*x + 1,       x^2 + 3*x]])
            sage: A.is_square() and A.constant_matrix().is_invertible()
            True
            sage: B = vector([2*x^2 + 6*x + 6, 0, x + 6])
            sage: X = A.solve_right_series_trunc(B, 4); X
            (2*x^3 + x^2, 5*x^3 + x^2 + 5*x + 6, 4*x^3 + 6*x^2 + 4*x)
            sage: B == A*X % x**4
            True
            sage: B = Matrix(pR, 3, 2,
            ....:            [[5*x^2 + 6*x + 3, 4*x^2 + 6*x + 4],
            ....:             [  x^2 + 4*x + 2,         5*x + 2],
            ....:             [        5*x + 3,               0]])
            sage: A.solve_right_series_trunc(B, 3)
            [  3*x^2 + x + 1 5*x^2 + 4*x + 3]
            [6*x^2 + 3*x + 1         4*x + 1]
            [      6*x^2 + 1   2*x^2 + x + 4]
            sage: X = A.solve_right_series_trunc(B, 37); B == A*X % x**37
            True

        Dimensions of input are checked::

            sage: A.solve_right_series_trunc(B[:2,:], 3)
            Traceback (most recent call last):
            ...
            ValueError: number of rows of self must equal number of rows of right-hand side

        Raises an exception when no solution::

            sage: A[:,2:].solve_right_series_trunc(B, 4)
            Traceback (most recent call last):
            ...
            ValueError: matrix equation has no solutions

            sage: Ax = x*A; C = vector(pR, [1,1,1])
            sage: Ax.solve_right_series_trunc(C, 5)
            Traceback (most recent call last):
            ...
            ValueError: matrix equation has no solutions

        Supports rectangular and rank-deficient cases::

            sage: A[:2,:].solve_right_series_trunc(B[:2,:],4)
            [    5*x^2 + 4*x           x + 4]
            [  x^2 + 3*x + 5 3*x^2 + 4*x + 4]
            [        5*x + 3         3*x + 2]

            sage: V = Matrix([[2*x^2 + 5*x + 1], [3*x^2+4]])
            sage: A[:,:2].solve_right_series_trunc(A[:,:2]*V, 4) == V
            True

            sage: A[:,1] = (x+1) * A[:,0]; A[:,2] = (x+5) * A[:,0]
            sage: B = (3*x^3+x^2+2)*A[:,0]
            sage: A.solve_right_series_trunc(B, 6)
            [4*x^2 + 6*x + 2]
            [      3*x^2 + x]
            [              0]

        .. SEEALSO::

            :meth:`solve_left_series_trunc` .
        """
        from sage.structure.element import Vector
        if isinstance(B, Vector):
            try:
                return self.transpose().solve_left_series_trunc(B, d)
            except ValueError as e:
                raise ValueError(str(e).replace('column', 'row'))
        else:
            try:
                return self.transpose().solve_left_series_trunc(B.transpose(), d).transpose()
            except ValueError as e:
                raise ValueError(str(e).replace('column', 'row'))

    def row_degrees(self, shifts=None):
        r"""
        Return the (shifted) row degrees of this matrix.

        For a given polynomial matrix `M = (M_{i,j})_{i,j}` with `m` rows and
        `n` columns, its row degrees is the tuple `(d_1,\ldots,d_m)` where `d_i
        = \max_j(\deg(M_{i,j}))` for `1\leq i \leq m`. Thus, `d_i=-1` if
        the `i`-th row of `M` is zero, and `d_i \geq 0` otherwise.

        For given shifts `s_1,\ldots,s_n \in \ZZ`, the shifted row degrees of
        `M` is `(d_1,\ldots,d_m)` where `d_i = \max_j(\deg(M_{i,j})+s_j)`.
        Here, if the `i`-th row of `M` is zero then `d_i
        =\min(s_1,\ldots,s_n)-1`; otherwise, `d_i` is larger than this value.

        INPUT:

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        OUTPUT: list of integers

        REFERENCES:

        - [Wol1974]_ (Section 2.5, without shifts), and [VBB1992]_ (Section 3).

        - Up to changes of signs, shifted row degrees coincide with the notion
          of *defect* commonly used in the rational approximation literature
          (see for example [Bec1992]_ ).

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [[3*x+1, 0, 1], [x^3+3, 0, 0]])
            sage: M.row_degrees()
            [1, 3]

            sage: M.row_degrees(shifts=[0,1,2])
            [2, 3]

        A zero row in a polynomial matrix can be identified in the (shifted)
        row degrees as the entries equal to ``min(shifts)-1``::

            sage: M = Matrix(pR, [[3*x+1, 0, 1], [x^3+3, 0, 0], [0, 0, 0]])
            sage: M.row_degrees()
            [1, 3, -1]

            sage: M.row_degrees(shifts=[-2,1,2])
            [2, 1, -3]

        The row degrees of an empty matrix (`0\times n` or `m\times 0`) is
        not defined::

            sage: M = Matrix(pR, 0, 3)
            sage: M.row_degrees()
            Traceback (most recent call last):
            ...
            ValueError: empty matrix does not have row degrees

            sage: M = Matrix(pR, 3, 0)
            sage: M.row_degrees()
            Traceback (most recent call last):
            ...
            ValueError: empty matrix does not have row degrees
        """
        self._check_shift_dimension(shifts,row_wise=True)
        if self.ncols() == 0 or self.nrows() == 0:
            raise ValueError('empty matrix does not have row degrees')
        if shifts is None:
            return [ max([ self[i,j].degree() for j in range(self.ncols()) ])
                    for i in range(self.nrows()) ]
        zero_degree = min(shifts) - 1
        return [ max([ self[i,j].degree() + shifts[j]
            if self[i,j] != 0 else zero_degree
            for j in range(self.ncols()) ]) for i in range(self.nrows()) ]

    def column_degrees(self, shifts=None):
        r"""
        Return the (shifted) column degrees of this matrix.

        For a given polynomial matrix `M = (M_{i,j})_{i,j}` with `m` rows and
        `n` columns, its column degrees is the tuple `(d_1,\ldots,d_n)` where
        `d_j = \max_i(\deg(M_{i,j}))` for `1\leq j \leq n`. Thus, `d_j=-1` if
        the `j`-th column of `M` is zero, and `d_j \geq 0` otherwise.

        For given shifts `s_1,\ldots,s_m \in \ZZ`, the shifted column degrees of
        `M` is `(d_1,\ldots,d_n)` where `d_j = \max_i(\deg(M_{i,j})+s_i)`.
        Here, if the `j`-th column of `M` is zero then `d_j =
        \min(s_1,\ldots,s_m)-1`; otherwise `d_j` is larger than this value.

        INPUT:

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        OUTPUT: list of integers

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [[3*x+1, 0, 1], [x^3+3, 0, 0]])
            sage: M.column_degrees()
            [3, -1, 0]

            sage: M.column_degrees(shifts=[0,2])
            [5, -1, 0]

        A zero column in a polynomial matrix can be identified in the (shifted)
        column degrees as the entries equal to ``min(shifts)-1``::

            sage: M.column_degrees(shifts=[-2,1])
            [4, -3, -2]

        The column degrees of an empty matrix (`0\times n` or `m\times 0`) is
        not defined::

            sage: M = Matrix(pR, 0, 3)
            sage: M.column_degrees()
            Traceback (most recent call last):
            ...
            ValueError: empty matrix does not have column degrees

            sage: M = Matrix(pR, 3, 0)
            sage: M.column_degrees()
            Traceback (most recent call last):
            ...
            ValueError: empty matrix does not have column degrees

        .. SEEALSO::

            The documentation of :meth:`row_degrees`.
        """
        self._check_shift_dimension(shifts,row_wise=False)
        if self.ncols() == 0 or self.nrows() == 0:
            raise ValueError('empty matrix does not have column degrees')
        if shifts is None:
            return [ max([ self[i,j].degree() for i in range(self.nrows()) ])
                    for j in range(self.ncols()) ]
        zero_degree = min(shifts) - 1
        return [ max([ self[i,j].degree() + shifts[i]
            if self[i,j] != 0 else zero_degree
            for i in range(self.nrows()) ]) for j in range(self.ncols()) ]

    def leading_matrix(self, shifts=None, row_wise=True):
        r"""
        Return the (shifted) leading matrix of this matrix.

        Let `M` be a univariate polynomial matrix in `\Bold{K}[x]^{m \times
        n}`. Working row-wise and without shifts, its leading matrix is the
        matrix in `\Bold{K}^{m \times n}` formed by the leading coefficients of
        the entries of `M` which reach the degree of the corresponding row.

        More precisely, if working row-wise, let `s_1,\ldots,s_n \in \ZZ`
        be a shift, and let `(d_1,\ldots,d_m)` denote the shifted row degrees of
        `M`. Then, the shifted leading matrix of `M` is the matrix in
        `\Bold{K}^{m \times n}` whose entry `i,j` is the coefficient of degree
        `d_i-s_j` of the entry `i,j` of `M`.

        If working column-wise, let `s_1,\ldots,s_m \in \ZZ` be a shift,
        and let `(d_1,\ldots,d_n)` denote the shifted column degrees of `M`.
        Then, the shifted leading matrix of `M` is the matrix in `\Bold{K}^{m
        \times n}` whose entry `i,j` is the coefficient of degree `d_j-s_i` of
        the entry `i,j` of `M`.

        INPUT:

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        - ``row_wise`` -- boolean (default: ``True``); ``True`` if
          working row-wise (see the class description)

        OUTPUT: a matrix over the base field

        REFERENCES:

        [Wol1974]_ (Section 2.5, without shifts) and [VBB1992]_ (Section 3).

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [[3*x+1, 0, 1], [x^3+3, 0, 0]])
            sage: M.leading_matrix()
            [3 0 0]
            [1 0 0]

            sage: M.leading_matrix().base_ring()
            Finite Field of size 7

            sage: M.leading_matrix(shifts=[0,1,2])
            [0 0 1]
            [1 0 0]

            sage: M.leading_matrix(row_wise=False)
            [0 0 1]
            [1 0 0]

            sage: M.leading_matrix(shifts=[-2,1], row_wise=False)
            [0 0 1]
            [1 0 0]

            sage: M.leading_matrix(shifts=[2,0], row_wise=False)
            [3 0 1]
            [1 0 0]
        """
        self._check_shift_dimension(shifts,row_wise)
        from sage.matrix.constructor import matrix
        if row_wise:
            row_degrees = self.row_degrees(shifts)
            if shifts is None:
                return matrix([ [ self[i,j].leading_coefficient()
                    if self[i,j].degree() == row_degrees[i] else 0
                    for j in range(self.ncols()) ]
                    for i in range(self.nrows()) ])
            return matrix([ [ self[i,j].leading_coefficient()
                if self[i,j].degree() + shifts[j] == row_degrees[i] else 0
                for j in range(self.ncols()) ]
                for i in range(self.nrows()) ])
        else:
            column_degrees = self.column_degrees(shifts)
            if shifts is None:
                return matrix([ [ self[i,j].leading_coefficient()
                    if self[i,j].degree() == column_degrees[j] else 0
                    for j in range(self.ncols()) ]
                    for i in range(self.nrows()) ])
            return matrix([ [ self[i,j].leading_coefficient()
                if self[i,j].degree() + shifts[i] == column_degrees[j] else 0
                for j in range(self.ncols()) ]
                for i in range(self.nrows()) ])

    def _is_empty_popov(self, row_wise=True, include_zero_vectors=True):
        r"""
        Assuming that this matrix is empty, that is, of dimensions `0\times n`
        or `m\times 0`, return a boolean indicating if it is in shifted Popov
        form. If zero vectors are allowed in shifted reduced forms, this always
        returns true. Otherwise, by convention and if working row-wise, for
        `n\geq 0` the `0\times n` matrix is in shifted Popov form for all
        shifts, and for `m>0` the `m \times 0` matrix is not in shifted Popov
        form for any shift. The convention is similar if working column-wise.

        The behaviour of this method for non-empty matrices is not defined.

        INPUT:

        - ``row_wise`` -- boolean (default: ``True``); ``True`` if
          one considers the row-wise shifted Popov form

        - ``include_zero_vectors`` -- boolean (default: ``True``);
          ``False`` if one does not allow zero rows in row reduced forms (resp.
          zero columns in column reduced forms)

        OUTPUT: boolean

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, 0, 0)
            sage: M._is_empty_popov()
            True
            sage: M._is_empty_popov(include_zero_vectors=False)
            True

            sage: M = Matrix(pR, 0, 3)
            sage: M._is_empty_popov(include_zero_vectors=False)
            True
            sage: M._is_empty_popov(row_wise=False)
            True
            sage: M._is_empty_popov(row_wise=False,include_zero_vectors=False)
            False

        .. SEEALSO::

            :meth:`is_popov` .
        """
        if include_zero_vectors:
            return True
        else:
            # assume we work row-wise, self is in shifted Popov form iff self.nrows()==0:
            # --> if self.nrows()==0, then self is in shifted Popov form
            # --> if self.nrows()>0, then self.ncols()==0 and thus self is not
            # in shifted Popov form
            return self.nrows() == 0 if row_wise else self.ncols() == 0

    def is_reduced(self,
            shifts=None,
            row_wise=True,
            include_zero_vectors=True):
        r"""
        Return a boolean indicating whether this matrix is in (shifted) reduced
        form.

        An `m \times n` univariate polynomial matrix `M` is said to be in
        shifted row reduced form if it has `k` nonzero rows with `k \leq n` and
        its shifted leading matrix has rank `k`. Equivalently, when considering
        all the matrices obtained by left-multiplying `M` by a unimodular
        matrix, then the shifted row degrees of `M` -- once sorted in
        nondecreasing order -- is lexicographically minimal.

        Similarly, `M` is said to be in shifted column reduced form if it has
        `k` nonzero columns with `k \leq m` and its shifted leading matrix has
        rank `k`.

        Sometimes, one forbids `M` to have zero rows (resp. columns) in the
        above definitions; an optional parameter allows one to adopt this more
        restrictive setting.

        INPUT:

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        - ``row_wise`` -- boolean (default: ``True``); ``True`` if
          working row-wise (see the class description)

        - ``include_zero_vectors`` -- boolean (default: ``True``);
          ``False`` if one does not allow zero rows in row reduced forms (resp.
          zero columns in column reduced forms)

        OUTPUT: boolean

        REFERENCES:

        [Wol1974]_ (Section 2.5, without shifts) and [VBB1992]_ (Section 3).

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [[3*x+1, 0, 1], [x^3+3, 0, 0]])
            sage: M.is_reduced()
            False

            sage: M.is_reduced(shifts=[0,1,2])
            True

            sage: M.is_reduced(shifts=[2,0], row_wise=False)
            True

            sage: M.is_reduced(shifts=[2,0], row_wise=False,
            ....:              include_zero_vectors=False)
            False

            sage: M = Matrix(pR, [[3*x+1, 0, 1], [x^3+3, 0, 0], [0, 1, 0]])
            sage: M.is_reduced(shifts=[2,0,0], row_wise=False)
            True

        .. SEEALSO::

            :meth:`leading_matrix` ,
            :meth:`reduced_form` .
        """
        self._check_shift_dimension(shifts,row_wise)
        if self.ncols() == 0 or self.nrows() == 0:
            return self._is_empty_popov(row_wise,include_zero_vectors)
        if include_zero_vectors:
            number_generators =                                           \
                [self[i,:] != 0 for i in range(self.nrows())].count(True) \
                if row_wise else                                          \
                [self[:,j] != 0 for j in range(self.ncols())].count(True)
        else:
            number_generators = self.nrows() if row_wise else self.ncols()
        return number_generators == \
            self.leading_matrix(shifts, row_wise).rank()

    def leading_positions(self,
            shifts=None,
            row_wise=True,
            return_degree=False):
        r"""
        Return the (shifted) leading positions (also known as the pivot
        indices), and optionally the (shifted) pivot degrees of this matrix.

        If working row-wise, for a given shift `s_1,\ldots,s_n \in
        \ZZ`, taken as `(0,\ldots,0)` by default, and a row vector of
        univariate polynomials `[p_1,\ldots,p_n]`, the leading position of
        this vector is the index `j` of the rightmost nonzero entry `p_j` such
        that `\deg(p_j) + s_j` is equal to the shifted row degree of the vector.
        Then the pivot degree of the vector is the degree `\deg(p_j)`.

        For the zero row, both the leading positions and degree are `-1`.  For
        a `m \times n` polynomial matrix, the leading positions and pivot
        degrees are the two lists containing the leading positions and the
        pivot degrees of its rows.

        The definition is similar if working column-wise (instead of rightmost
        nonzero entry, we choose the bottommost nonzero entry).

        INPUT:

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        - ``row_wise`` -- boolean (default: ``True``); ``True`` if
          working row-wise (see the class description)

        - ``return_degree`` -- boolean (default: ``False``); ``True``
          implies that the pivot degrees are returned

        OUTPUT: list of integers if ``return_degree=False``; a pair of lists
        of integers otherwise

        REFERENCES:

        [Kai1980]_ (Section 6.7.2, without shifts).

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [[3*x+1, 0, 1], [x^3+3, 0, 0]])
            sage: M.leading_positions()
            [0, 0]

            sage: M.leading_positions(return_degree=True)
            ([0, 0], [1, 3])

            sage: M.leading_positions(shifts=[0,5,2], return_degree=True)
            ([2, 0], [0, 3])

            sage: M.leading_positions(row_wise=False, return_degree=True)
            ([1, -1, 0], [3, -1, 0])

            sage: M.leading_positions(shifts=[1,2], row_wise=False,
            ....:                     return_degree=True)
            ([1, -1, 0], [3, -1, 0])

        In case several entries in the row (resp. column) reach the shifted row
        (resp. column) degree, the leading position is chosen as the rightmost
        (resp. bottommost) such entry::

            sage: M.leading_positions(shifts=[0,5,1], return_degree=True)
            ([2, 0], [0, 3])

            sage: M.leading_positions(shifts=[2,0], row_wise=False,
            ....:                     return_degree=True)
            ([1, -1, 0], [3, -1, 0])

        The leading positions and pivot degrees of an empty matrix (`0\times n`
        or `m\times 0`) is not defined::

            sage: M = Matrix(pR, 0, 3)
            sage: M.leading_positions()
            Traceback (most recent call last):
            ...
            ValueError: empty matrix does not have leading positions

            sage: M.leading_positions(row_wise=False)
            Traceback (most recent call last):
            ...
            ValueError: empty matrix does not have leading positions

            sage: M = Matrix(pR, 3, 0)
            sage: M.leading_positions(row_wise=False)
            Traceback (most recent call last):
            ...
            ValueError: empty matrix does not have leading positions
        """
        self._check_shift_dimension(shifts,row_wise)

        if self.ncols() == 0 or self.nrows() == 0:
            raise ValueError('empty matrix does not have leading positions')

        if row_wise:
            row_degrees = self.row_degrees(shifts)
            if shifts is None:
                pivot_index = [ (-1 if row_degrees[i] == -1 else
                    max( [ j for j in range(self.ncols()) if
                    (self[i,j].degree() == row_degrees[i]) ] ))
                    for i in range(self.nrows()) ]
            else:
                zero_degree=min(shifts) - 1
                pivot_index = [ (-1 if row_degrees[i] == zero_degree else
                    max( [ j for j in range(self.ncols()) if
                    (self[i,j] != 0 and
                    self[i,j].degree() + shifts[j] == row_degrees[i]) ] ))
                    for i in range(self.nrows()) ]
            pivot_degree = [ (-1 if pivot_index[i] == -1 else
                self[i,pivot_index[i]].degree())
                for i in range(self.nrows()) ]
            return (pivot_index,pivot_degree) if return_degree else pivot_index

        # now in the column-wise case
        column_degrees = self.column_degrees(shifts)
        if shifts is None:
            pivot_index = [ (-1 if column_degrees[j] == -1 else
                max( [ i for i in range(self.nrows()) if
                (self[i,j].degree() == column_degrees[j]) ] ))
                for j in range(self.ncols()) ]
        else:
            zero_degree=min(shifts) - 1
            pivot_index = [ (-1 if column_degrees[j] == zero_degree else
                max( [ i for i in range(self.nrows()) if
                (self[i,j] != 0 and
                self[i,j].degree() + shifts[i] == column_degrees[j]) ] ))
                for j in range(self.ncols()) ]
        pivot_degree = [ (-1 if pivot_index[j] == -1 else
            self[pivot_index[j],j].degree())
            for j in range(self.ncols()) ]
        return (pivot_index,pivot_degree) if return_degree else pivot_index

    def is_weak_popov(self,
            shifts=None,
            row_wise=True,
            ordered=False,
            include_zero_vectors=True):
        r"""
        Return a boolean indicating whether this matrix is in (shifted)
        (ordered) weak Popov form.

        If working row-wise (resp. column-wise), a polynomial matrix is said to
        be in weak Popov form if the leading positions of its nonzero rows
        (resp. columns) are pairwise distinct. For the ordered weak Popov form,
        these positions must be strictly increasing, except for the possibly
        repeated -1 entries which are at the end. For the shifted variants, see
        the class description for an introduction to shifts.

        INPUT:

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        - ``row_wise`` -- boolean (default: ``True``); ``True`` if
          working row-wise (see the class description)

        - ``ordered`` -- boolean (default: ``False``); ``True`` if
          checking for an ordered weak Popov form

        - ``include_zero_vectors`` -- boolean (default: ``True``);
          ``False`` if one does not allow zero rows (resp. zero columns) in
          (ordered) weak Popov forms

        OUTPUT: boolean

        REFERENCES:

        [Kai1980]_ (Section 6.7.2, square case without shifts), [MS2003]_
        (without shifts), [BLV1999]_ .

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix([ [x^3+3*x^2+6*x+6, 3*x^2+3*x+6, 4*x^2+x+3],
            ....:              [5,               1,           0        ],
            ....:              [2*x^2+2,         2*x+5,       x^2+4*x+6] ])
            sage: M.is_weak_popov()
            True

        One can check whether the leading positions, in addition to being
        pairwise distinct, are actually in increasing order::

            sage: M.is_weak_popov(ordered=True)
            True

            sage: N = M.with_swapped_rows(1, 2)
            sage: N.is_weak_popov()
            True
            sage: N.is_weak_popov(ordered=True)
            False

        Shifts and orientation (row-wise or column-wise) are supported::

            sage: M.is_weak_popov(shifts=[2,3,1])
            False

            sage: M.is_weak_popov(shifts=[0,2,0], row_wise=False,
            ....:                 ordered=True)
            True

        Rectangular matrices are supported::

            sage: M = Matrix([
            ....:    [  x^3+5*x^2+5*x+1,       5,       6*x+4,         0],
            ....:    [      6*x^2+3*x+1,       1,           2,         0],
            ....:    [2*x^3+4*x^2+6*x+4, 5*x + 1, 2*x^2+5*x+5, x^2+5*x+6]
            ....:     ])
            sage: M.is_weak_popov(shifts=[0,2,1,3])
            True

            sage: M.is_weak_popov(shifts=[0,2,1,3], ordered=True)
            True

        Zero rows (resp. columns) can be forbidden::

            sage: M = Matrix([
            ....:   [      6*x+4,       0,             5*x+1, 0],
            ....:   [          2, 5*x + 1,       6*x^2+3*x+1, 0],
            ....:   [2*x^2+5*x+5,       1, 2*x^3+4*x^2+6*x+4, 0]
            ....:   ])
            sage: M.is_weak_popov(shifts=[2,1,0], row_wise=False,
            ....:                 ordered=True)
            True

            sage: M.is_weak_popov(shifts=[2,1,0], row_wise=False,
            ....:                 include_zero_vectors=False)
            False

        .. SEEALSO::

            :meth:`weak_popov_form` .
        """
        self._check_shift_dimension(shifts,row_wise)
        if self.ncols() == 0 or self.nrows() == 0:
            return self._is_empty_popov(row_wise)
        leading_positions = self.leading_positions(shifts, row_wise)
        # here, because of the below sorting and of the convention that zero
        # rows (resp. columns) are at the bottom (resp. right) of the matrix in
        # the row-wise case (resp. column-wise case), it will be convenient to
        # have leading position ncols (resp. nrows) for these zero vectors
        pos_zero_vec = self.ncols() if row_wise else self.nrows()
        leading_positions = [pos if pos >= 0 else pos_zero_vec + 1
                             for pos in leading_positions]
        # leading positions should not have duplicates, which is equivalent to:
        # once sorted, it doesn't contain a pair of equal successive entries
        if not ordered:
            leading_positions.sort()
        # check that there is no zero vector, if it is forbidden
        # (in the ordered case, even though we have not sorted, this test of
        # the last leading position is sufficient: in any case, if there is a
        # zero vector followed by a nonzero one, the testing loop below will
        # return False)
        if leading_positions[-1] > pos_zero_vec and not include_zero_vectors:
            return False
        # it remains to test whether leading_positions is strictly increasing
        # until it reaches a zero vector (note this also checks that there is
        # no zero vector encountered after a nonzero one)
        for index,next_leading_position in enumerate(leading_positions[1:]):
            if next_leading_position <= pos_zero_vec and \
                    next_leading_position <= leading_positions[index]:
                return False
        return True

    def is_popov(self,
            shifts=None,
            row_wise=True,
            up_to_permutation=False,
            include_zero_vectors=True):
        r"""
        Return a boolean indicating whether this matrix is in (shifted) Popov
        form.

        If working row-wise (resp. column-wise), a polynomial matrix is said to
        be in Popov form if it has no zero row above a nonzero row (resp. no
        zero column to the left of a nonzero column), the leading positions of
        its nonzero rows (resp. columns) are strictly increasing, and for each
        row (resp. column) the pivot entry is monic and has degree strictly
        larger than the other entries in its column (resp. row).

        Since other conventions concerning the ordering of the rows (resp.
        columns) are sometimes useful, an optional argument allows one to test
        whether the matrix is in Popov form up to row (resp. column)
        permutation. For example, there is an alternative definition which
        replaces "leading positions strictly increasing" by "row (resp. column)
        degree nondecreasing, and for rows (resp. columns) of same degree,
        leading positions increasing".

        INPUT:

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        - ``row_wise`` -- boolean (default: ``True``); ``True`` if
          working row-wise (see the class description)

        - ``up_to_permutation`` -- (option, default: ``False``) boolean,
          ``True`` if testing Popov form up to row permutation (if working
          row-wise).

        - ``include_zero_vectors`` -- boolean (default: ``True``);
          ``False`` if one does not allow zero rows (resp. zero columns) in
          Popov forms

        OUTPUT: boolean

        REFERENCES:

        For the square case, without shifts: [Pop1972]_ and [Kai1980]_ (Section
        6.7.2). For the general case: [BLV2006]_ .

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [[x^4+6*x^3+4*x+4, 3*x+6,     3  ],
            ....:                 [x^2+6*x+6,       x^2+5*x+5, 2  ],
            ....:                 [3*x,             6*x+5,     x+5]])
            sage: M.is_popov()
            True

            sage: M.is_popov(shifts=[0,1,2])
            True

            sage: M[:,:2].is_popov()
            False

            sage: M[:2,:].is_popov(shifts=[0,1,2])
            True

            sage: M = Matrix(pR, [[x^4+3*x^3+x^2+2*x+6, x^3+5*x^2+5*x+1],
            ....:                 [6*x+1,               x^2+4*x+1      ],
            ....:                 [6,                   6              ]])
            sage: M.is_popov(row_wise=False)
            False

            sage: M.is_popov(shifts=[0,2,3], row_wise=False)
            True

        One can forbid zero rows (or columns if not working row-wise)::

            sage: N = Matrix(pR, [[x^4+3*x^3+x^2+2*x+6, 6*x+1     ],
            ....:                 [5*x^2+5*x+1,         x^2+4*x+1 ],
            ....:                 [0,                   0         ]])

            sage: N.is_popov()
            True

            sage: N.is_popov(include_zero_vectors=False)
            False

        One can verify Popov form up to row permutation (or column permutation
        if not working row-wise)::

            sage: M.swap_columns(0, 1)
            sage: M.is_popov(shifts=[0,2,3], row_wise=False)
            False

            sage: M.is_popov(shifts=[0,2,3], row_wise=False,
            ....:            up_to_permutation=True)
            True

            sage: N.swap_rows(0, 2)

            sage: N.is_popov()
            False

            sage: N.is_popov(up_to_permutation=True)
            True
        """
        # the matrix should be in weak Popov form (ordered except if
        # up_to_permutation==True)
        if self.ncols() == 0 or self.nrows() == 0:
            return self._is_empty_popov(row_wise)
        if not self.is_weak_popov(shifts,
                                  row_wise,
                                  not up_to_permutation,
                                  include_zero_vectors):
            return False

        # pivot entries should be monic, and pivot degrees must be the greatest
        # in their column (or in their row if column-wise)
        leading_positions,pivot_degree = self.leading_positions(shifts, row_wise,
                return_degree=True)
        for i,index in enumerate(leading_positions):
            if index >= 0:
                if row_wise:
                    if not self[i,index].is_monic():
                        return False
                    for k in range(self.nrows()):
                        if k == i:
                            continue
                        if self[k,index].degree() >= pivot_degree[i]:
                            return False
                else: # now column-wise
                    if not self[index,i].is_monic():
                        return False
                    for k in range(self.ncols()):
                        if k == i:
                            continue
                        if self[index,k].degree() >= pivot_degree[i]:
                            return False
        return True

    def is_hermite(self,
            row_wise=True,
            lower_echelon=False,
            include_zero_vectors=True):
        r"""
        Return a boolean indicating whether this matrix is in Hermite form.

        If working row-wise, a polynomial matrix is said to be in Hermite form
        if it is in row echelon form with all pivot entries monic and such that
        all entries above a pivot have degree less than this pivot. Being in
        row echelon form means that all zero rows are gathered at the bottom of
        the matrix, and in each nonzero row the pivot (leftmost nonzero entry)
        is strictly to the right of the pivot of the row just above this row.

        Note that, for any integer `d` strictly greater than all degrees
        appearing in the Hermite form, then the Hermite form coincides with the
        shifted Popov form with the shifts `((n-1)d,\ldots,2d,d,0)`, where `n`
        is the column dimension.

        If working column-wise, a polynomial matrix is said to be in Hermite
        form if it is in column echelon form with all pivot entries monic and
        such that all entries to the left of a pivot have degree less than this
        pivot. Being in column echelon form means that all zero columns are
        gathered at the right-hand side of the matrix, and in each nonzero
        column the pivot (topmost nonzero entry) is strictly below the pivot of
        the column just to the left of this row.

        Optional arguments provide support of alternative definitions,
        concerning the choice of upper or lower echelon forms and concerning
        whether zero rows (resp. columns) are allowed.

        INPUT:

        - ``row_wise`` -- boolean (default: ``True``); ``True`` if
          working row-wise (see the class description)

        - ``lower_echelon`` -- boolean (default: ``False``);
          ``False`` if working with upper triangular Hermite forms, ``True`` if
          working with lower triangular Hermite forms.

        - ``include_zero_vectors`` -- boolean (default: ``True``);
          ``False`` if one does not allow zero rows (resp. zero columns) in
          Hermite forms

        OUTPUT: boolean

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [[x^4+6*x^3+4*x+4, 3*x+6,     3  ],
            ....:                 [0,               x^2+5*x+5, 2  ],
            ....:                 [0,               0,         x+5]])

            sage: M.is_hermite()
            True
            sage: M.is_hermite(row_wise=False)
            True
            sage: M.is_hermite(row_wise=False, lower_echelon=True)
            False

            sage: N = Matrix(pR, [[x+5, 0,               0        ],
            ....:                 [2,   x^4+6*x^3+4*x+4, 0        ],
            ....:                 [3,   3*x^3+6,         x^2+5*x+5]])
            sage: N.is_hermite()
            False
            sage: N.is_hermite(lower_echelon=True)
            True
            sage: N.is_hermite(row_wise=False)
            False
            sage: N.is_hermite(row_wise=False, lower_echelon=True)
            False

        Rectangular matrices with zero rows are supported. Zero rows (resp.
        columns) can be forbidden, and otherwise they should be at the bottom
        (resp. the right-hand side) of the matrix::

            sage: N[:,1:].is_hermite(lower_echelon=True)
            False
            sage: N[[1,2,0],1:].is_hermite(lower_echelon=True)
            True
            sage: N[:2,:].is_hermite(row_wise=False, lower_echelon=True)
            True
            sage: N[:2,:].is_hermite(row_wise=False,
            ....:                    lower_echelon=True,
            ....:                    include_zero_vectors=False)
            False

        .. SEEALSO::

            :meth:`hermite_form` .
        """
        # shift for lower echelon
        shift = [j*(self.degree() + 1) for j in range(self.ncols())] \
                if row_wise else \
                [(self.nrows() - j)*(self.degree() + 1) for j in range(self.nrows())]
        # if upper echelon, reverse shift
        if not lower_echelon:
            shift.reverse()
        return self.is_popov(shifts=shift,
                row_wise=row_wise,
                include_zero_vectors=include_zero_vectors)

    def weak_popov_form(self,
            transformation=False,
            shifts=None,
            row_wise=True,
            ordered=False,
            include_zero_vectors=True):
        r"""
        Return a (shifted) (ordered) weak Popov form of this matrix.

        See :meth:`is_weak_popov` for a definition of weak Popov forms. If the
        input matrix is `A`, a weak Popov form of `A` is any matrix `P` in weak
        Popov form and such that `UA = P` for some unimodular matrix `U`. The
        latter matrix is called the transformation, and the first optional
        argument allows one to specify whether to return this transformation.

        Sometimes, one forbids weak Popov forms to have zero rows (resp.
        columns) in the above definitions; an optional parameter allows one to
        adopt this more restrictive setting. If zero rows (resp. columns) are
        allowed, the convention here is to place them as the bottommost rows
        (resp. the rightmost columns) of the output weak Popov form.

        Note that, if asking for the transformation and discarding zero vectors
        (i.e. ``transformation=True`` and ``include_zero_vectors=False``), then
        the returned transformation is still the complete unimodular matrix,
        including its bottommost rows (resp. rightmost columns) which
        correspond to zero rows (resp. columns) of the complete weak Popov
        form. In fact, this bottom part of the transformation yields a basis of
        the left (resp. right) kernel of the input matrix.

        INPUT:

        - ``transformation`` -- (default: ``False``) if this
          is ``True``, the transformation matrix `U` will be returned as well

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        - ``row_wise`` -- boolean (default: ``True``); ``True`` if
          working row-wise (see the class description)

        - ``ordered`` -- boolean (default: ``False``); ``True`` if
          seeking an ordered weak Popov form

        - ``include_zero_vectors`` -- boolean (default: ``True``);
          ``False`` if zero rows (resp. zero columns) should be discarded from
          the (ordered) weak Popov forms

        OUTPUT:

        A polynomial matrix which is a weak Popov form of ``self`` if
        ``transformation`` is ``False``; otherwise two polynomial matrices
        which are a weak Popov form of ``self`` and the corresponding
        unimodular transformation.

        ALGORITHM:

        This method implements the Mulders-Storjohann algorithm of [MS2003]_,
        straightforwardly extended to the case of shifted forms.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [
            ....:    [      6*x+4,       5*x^3+5*x,       6*x^2+2*x+2],
            ....:    [4*x^2+5*x+2, x^4+5*x^2+2*x+4, 4*x^3+6*x^2+6*x+5]])
            sage: P, U = M.weak_popov_form(transformation=True)
            sage: P
            [              4             x^2   6*x^2 + x + 2]
            [              2 4*x^2 + 2*x + 4               5]
            sage: U
            [2*x^2 + 1       4*x]
            [      4*x         1]
            sage: P.is_weak_popov() and U.is_invertible() and U*M == P
            True

        Demonstrating the ``ordered`` option::

            sage: P.leading_positions()
            [2, 1]
            sage: PP = M.weak_popov_form(ordered=True); PP
            [              2 4*x^2 + 2*x + 4               5]
            [              4             x^2   6*x^2 + x + 2]
            sage: PP.leading_positions()
            [1, 2]

        Demonstrating shifts::

            sage: P = M.weak_popov_form(shifts=[0,2,4]); P
            [            6*x^2 + 6*x + 4 5*x^4 + 4*x^3 + 5*x^2 + 5*x                     2*x + 2]
            [                          2             4*x^2 + 2*x + 4                           5]
            sage: P == M.weak_popov_form(shifts=[-10,-8,-6])
            True

        Column-wise form is the row-wise form of the transpose::

            sage: M.weak_popov_form() == M.T.weak_popov_form(row_wise=False).T
            True

        Zero vectors can be discarded::

            sage: M.weak_popov_form(row_wise=False)
            [x + 4     6     0]
            [    5     1     0]

            sage: # needs sage.combinat
            sage: P, U = M.weak_popov_form(transformation=True,
            ....:                          row_wise=False,
            ....:                          include_zero_vectors=False)
            sage: P
            [x + 4     6]
            [    5     1]
            sage: U
            [                5*x + 2         5*x^2 + 4*x + 4 3*x^3 + 3*x^2 + 2*x + 4]
            [                      1                       1                 2*x + 1]
            [                5*x + 5                       2                       6]
            sage: M*U[:,:2] == P and (M*U[:,2]).is_zero()
            True

        .. SEEALSO::

            :meth:`is_weak_popov` ,
            :meth:`reduced_form` ,
            :meth:`popov_form` ,
            :meth:`hermite_form` .
        """
        # if column-wise, call the algorithm on transpose
        if not row_wise:
            W = self.T.weak_popov_form(transformation,
                        shifts,
                        True,
                        ordered,
                        include_zero_vectors)
            return (W[0].T,W[1].T) if transformation else W.T
        # --> now, below, we are working row-wise
        # row dimension:
        m = self.nrows()
        # make shift nonnegative, required by main call _weak_popov_form
        self._check_shift_dimension(shifts,row_wise=True)
        if shifts is None:
            nonnegative_shifts = None
        else:
            min_shifts = min(shifts)
            nonnegative_shifts = [s-min_shifts for s in shifts]
        # call main procedure to compute weak Popov and transformation
        M = self.__copy__()
        U = M._weak_popov_form(transformation=transformation,
                shifts=nonnegative_shifts)
        # move zero rows to the bottom of the matrix
        from sage.combinat.permutation import Permutation
        zero_rows = []
        nonzero_rows = []
        for i in range(m):
            # note the "i+1" due to the format of permutation used below
            if M[i].is_zero():
                zero_rows.append(i+1)
            else:
                nonzero_rows.append(i+1)
        M.permute_rows(Permutation(nonzero_rows + zero_rows))
        if transformation:
            U.permute_rows(Permutation(nonzero_rows + zero_rows))
        # order other rows by increasing leading positions
        if ordered:
            lpos = M.leading_positions(nonnegative_shifts,row_wise=True)
            # find permutation that sorts leading_positions in increasing order
            # --> force max value to zero rows so that they remain bottom rows
            if include_zero_vectors: # otherwise, zero rows already removed
                for i in range(m):
                    if lpos[i] == -1:
                        lpos[i] = m
            sorted_lpos = sorted([(lpos[i],i+1) for i in range(m)])
            row_permutation = Permutation([elt[1] for elt in sorted_lpos])
            # apply permutation to weak Popov form and the transformation
            M.permute_rows(row_permutation)
            if transformation:
                U.permute_rows(row_permutation)
        # remove zero rows if asked to
        if not include_zero_vectors:
            M = M.delete_rows(range(m-len(zero_rows),m))
        # set immutable and return
        M.set_immutable()
        if transformation:
            U.set_immutable()
        return (M,U) if transformation else M

    def _weak_popov_form(self, transformation=False, shifts=None):
        """
        Transform this matrix in-place into weak Popov form.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(2^4, 'a')
            sage: PF.<x> = F[]
            sage: A = matrix(PF,[[1,  a*x^17 + 1 ],
            ....:                [0,  a*x^11 + a^2*x^7 + 1 ]])
            sage: M = A.__copy__()
            sage: U = M._weak_popov_form(transformation=True)
            sage: U * A == M
            True
            sage: M.is_weak_popov()
            True
            sage: U.is_invertible()
            True

            sage: PF.<x> = QQ[]
            sage: A = matrix(PF,3,[x,   x^2, x^3,
            ....:                  x^2, x^1, 0,
            ....:                  x^3, x^3, x^3])
            sage: A.weak_popov_form()
            [        x       x^2       x^3]
            [      x^2         x         0]
            [  x^3 - x x^3 - x^2         0]
            sage: M = A.__copy__()
            sage: U = M._weak_popov_form(transformation=True, shifts=[16,8,0])
            sage: M
            [               x              x^2              x^3]
            [               0         -x^2 + x       -x^4 + x^3]
            [               0                0 -x^5 + x^4 + x^3]
            sage: U * A == M
            True
        """
        cdef Py_ssize_t i, j
        cdef Py_ssize_t c, d, best, bestp

        cdef Py_ssize_t m = self.nrows()
        cdef Py_ssize_t n = self.ncols()

        cdef Matrix M = self
        cdef Matrix U

        cdef list to_row, conflicts

        R = self.base_ring()
        one = R.one()

        if transformation:
            from sage.matrix.constructor import identity_matrix
            U = identity_matrix(R, m)

        # initialise to_row and conflicts list
        to_row = [[] for i in range(n)]
        conflicts = []
        for i in range(m):
            bestp = -1
            best = -1
            for c in range(n):
                d = M.get_unsafe(i,c).degree()

                if shifts and d >= 0:
                    d += shifts[c]

                if d >= best:
                    bestp = c
                    best = d

            if best >= 0:
                to_row[bestp].append((i,best))
                if len(to_row[bestp]) > 1:
                    conflicts.append(bestp)

        # while there is a conflict, do a simple transformation
        while conflicts:
            c = conflicts.pop()
            row = to_row[c]
            i,ideg = row.pop()
            j,jdeg = row.pop()

            if jdeg > ideg:
                i,j = j,i
                ideg,jdeg = jdeg,ideg

            coeff = - M.get_unsafe(i,c).lc() / M.get_unsafe(j,c).lc()
            s = coeff * one.shift(ideg - jdeg)

            M.add_multiple_of_row_c(i, j, s, 0)
            if transformation:
                U.add_multiple_of_row_c(i, j, s, 0)

            row.append((j,jdeg))

            bestp = -1
            best = -1
            for c in range(n):
                d = M.get_unsafe(i,c).degree()

                if shifts and d >= 0:
                    d += shifts[c]

                if d >= best:
                    bestp = c
                    best = d

            if best >= 0:
                to_row[bestp].append((i,best))
                if len(to_row[bestp]) > 1:
                    conflicts.append(bestp)

        if transformation:
            return U

    def popov_form(self,
            transformation=False,
            shifts=None,
            row_wise=True,
            include_zero_vectors=True):
        r"""
        Return the (shifted) Popov form of this matrix.

        See :meth:`is_popov` for a definition of Popov forms. If the input
        matrix is `A`, the (shifted) Popov form of `A` is the unique matrix `P`
        in (shifted) Popov form and such that `UA = P` for some unimodular
        matrix `U`. The latter matrix is called the transformation, and the
        first optional argument allows one to specify whether to return this
        transformation. We refer to the description of :meth:`weak_popov_form`
        for an explanation of the option ``include_zero_vectors`` .

        INPUT:

        - ``transformation`` -- (default: ``False``) if this
          is ``True``, the transformation matrix `U` will be returned as well

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        - ``row_wise`` -- boolean (default: ``True``); ``True`` if
          working row-wise (see the class description)

        - ``include_zero_vectors`` -- boolean (default: ``True``);
          ``False`` if zero rows (resp. zero columns) should be discarded from
          the Popov forms

        OUTPUT:

        - A polynomial matrix which is the Popov form of ``self`` if
          ``transformation`` is ``False``; otherwise two polynomial matrices
          which are the Popov form of ``self`` and the corresponding unimodular
          transformation.

        ALGORITHM:

        This method implements the Mulders-Storjohann algorithm of [MS2003]_
        for transforming a weak Popov form into Popov form, straightforwardly
        extended to the case of shifted forms.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [
            ....:     [      6*x+4,       5*x^3+5*x,       6*x^2+2*x+2],
            ....:     [4*x^2+5*x+2, x^4+5*x^2+2*x+4, 4*x^3+6*x^2+6*x+5]])

            sage: # needs sage.combinat
            sage: P, U = M.popov_form(transformation=True)
            sage: P
            [            4 x^2 + 4*x + 1             3]
            [            0       4*x + 1 x^2 + 6*x + 1]
            sage: U
            [            x             2]
            [5*x^2 + x + 6       3*x + 2]
            sage: P.is_popov() and U.is_invertible() and U*M == P
            True

        Demonstrating shifts and specific case of Hermite form::

            sage: # needs sage.combinat
            sage: P = M.popov_form(shifts=[0,2,4]); P
            [              4*x^2 + 3*x + 4 x^4 + 3*x^3 + 5*x^2 + 5*x + 5                             0]
            [                            6               5*x^2 + 6*x + 5                             1]
            sage: P.is_popov(shifts=[0,2,4])
            True
            sage: P == M.popov_form(shifts=[-6,-4,-2])
            True
            sage: dd = sum(M.row_degrees()) + 1
            sage: M.popov_form(shifts=[2*dd,dd,0]) == M.hermite_form()
            True

        Column-wise form is the row-wise form of the transpose::

            sage: M.popov_form() == M.T.popov_form(row_wise=False).T
            True

        Zero vectors can be discarded::

            sage: M.popov_form(row_wise=False)
            [x + 2     6     0]
            [    0     1     0]

            sage: # needs sage.combinat
            sage: P, U = M.popov_form(transformation=True,
            ....:                     row_wise=False,
            ....:                     include_zero_vectors=False)
            sage: P
            [x + 2     6]
            [    0     1]
            sage: U
            [        3*x^2 + 6*x + 3         5*x^2 + 4*x + 4 3*x^3 + 3*x^2 + 2*x + 4]
            [                      3                       1                 2*x + 1]
            [                5*x + 2                       2                       6]
            sage: M*U[:,:2] == P and (M*U[:,2]).is_zero()
            True

        .. SEEALSO::

            :meth:`is_popov` ,
            :meth:`reduced_form` ,
            :meth:`weak_popov_form` ,
            :meth:`hermite_form` .
        """
        # if column-wise, call the algorithm on transpose
        if not row_wise:
            P = self.T.popov_form(transformation,
                        shifts,
                        True,
                        include_zero_vectors)
            return (P[0].T,P[1].T) if transformation else P.T
        # --> now, below, we are working row-wise
        # row dimension:
        nrows_zero = self.nrows()

        # compute row-wise weak Popov form:
        # -> non-ordered since we will soon order rows otherwise anyway
        # -> without zero rows, we will re-insert them later if asked to
        WP = self.weak_popov_form(transformation,shifts,True,False,False)
        if transformation:
            P,UU = WP[0].__copy__(),WP[1]
        else:
            P = WP.__copy__()
        m = P.nrows()
        # for now, only consider rows of transformation corresponding to
        # nonzero rows, other rows will be reinserted later
        if transformation:
            U = UU[:m].__copy__()

        # compute leading positions and shifted row degrees
        lpos,rdeg = P.leading_positions(shifts,True,True)
        if shifts is not None:
            rdeg = [rdeg[i] + shifts[lpos[i]] for i in range(m)]

        # 1/ transform P into ascending order (as defined in
        # [Mulders&Storjohann, 2003, p394], recall here P has no zero rows)
        # -> sort the (degree,pivot) couples by lex order,
        #        keeping track of the performed permutation
        # -> and permute P,U,lpos,rdeg accordingly
        from sage.combinat.permutation import Permutation
        sorted_rdeg_lpos = sorted([(rdeg[i],lpos[i],i+1) for i in range(m)])
        rdeg = [elt[0] for elt in sorted_rdeg_lpos]
        lpos = [elt[1] for elt in sorted_rdeg_lpos]
        row_permutation = Permutation([elt[2] for elt in sorted_rdeg_lpos])
        P.permute_rows(row_permutation)
        if transformation:
            U.permute_rows(row_permutation)

        # 2/ ensure all pivots are monic
        for i in range(m):
            inv_lc = 1/P[i,lpos[i]].leading_coefficient()
            P.rescale_row(i,inv_lc)
            if transformation:
                U.rescale_row(i,inv_lc)

        # 3/ reduce degrees as much as possible, row by row
        # (this works because of the above ascending order)
        for i in range(1,m):
            # use rows k=0...i-1 to reduce degrees of row i in column lpos[k]
            delta = 0
            while delta >= 0:
                # see [Mulders&Storjohann, 2003, Algo. PopovForm, p396]
                delta = -1
                j = -1
                for k in range(i):
                    if P[i,lpos[k]].degree() - P[k,lpos[k]].degree() > delta:
                        delta = P[i,lpos[k]].degree() - P[k,lpos[k]].degree()
                        j = k
                if delta>=0:
                    # recall the leading coefficient of P[j,lpos[j]] is 1
                    c = - P[i,lpos[j]].leading_coefficient()
                    shifted_row_Pj = c * P[j,:].shift(delta)
                    P[i,:] = P[i,:] + shifted_row_Pj
                    if transformation:
                        shifted_row_Uj = c * U[j,:].shift(delta)
                        U[i,:] = U[i,:] + shifted_row_Uj

        # 4/ transform so as to have increasing leading positions
        sorted_lpos = sorted([(lpos[i],i+1) for i in range(m)])
        row_permutation = Permutation([elt[1] for elt in sorted_lpos])
        P.permute_rows(row_permutation)
        if transformation:
            U.permute_rows(row_permutation)

        # reinsert zero rows: in U in all cases, in P if asked to
        if transformation:
            U = U.stack(UU[m:,:])
        if include_zero_vectors:
            from sage.matrix.constructor import matrix
            P = P.stack(matrix(self.base_ring(),nrows_zero-m,self.ncols()))
        # return
        return (P,U) if transformation else P

    def reduced_form(self,
            transformation=None,
            shifts=None,
            row_wise=True,
            include_zero_vectors=True):
        r"""
        Return a row reduced form of this matrix (resp. a column reduced form
        if the optional parameter ``row_wise`` is set to ``False``).

        An `m \times n` univariate polynomial matrix `M` is said to be in
        (shifted) row reduced form if it has `k` nonzero rows with `k \leq n`
        and its (shifted) leading matrix has rank `k`. See :meth:`is_reduced`
        for more information.

        Currently, the implementation of this method is a direct call to
        :meth:`weak_popov_form`.

        INPUT:

        - ``transformation`` -- (default: ``False``) if this
          is ``True``, the transformation matrix `U` will be returned as well:
          this is a unimodular matrix over `\Bold{K}[x]` such that ``self``
          equals `UR`, where `R` is the output matrix.

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        - ``row_wise`` -- boolean (default: ``True``); ``True`` if
          working row-wise (see the class description)

        - ``include_zero_vectors`` -- boolean (default: ``True``);
          ``False`` if one does not allow zero rows in row reduced forms (resp.
          zero columns in column reduced forms)

        OUTPUT:

        - A polynomial matrix `R` which is a reduced form of ``self`` if
          ``transformation=False``; otherwise two polynomial matrices `R, U`
          such that `UA = R` and `R` is reduced and `U` is unimodular where `A`
          is ``self``.

        EXAMPLES::

            sage: pR.<x> = GF(3)[]
            sage: A = matrix(pR, 3, [x,   x^2, x^3,
            ....:                    x^2, x^1, 0,
            ....:                    x^3, x^3, x^3])
            sage: R = A.reduced_form(); R
            [        x           x^2       x^3]
            [      x^2             x         0]
            [  x^3 + 2*x x^3 + 2*x^2         0]
            sage: R.is_reduced()
            True
            sage: R2 = A.reduced_form(shifts=[6,3,0]); R2
            [                x               x^2               x^3]
            [                0         2*x^2 + x       2*x^4 + x^3]
            [                0                 0 2*x^5 + x^4 + x^3]
            sage: R2.is_reduced(shifts=[6,3,0])
            True
            sage: R2.is_reduced()
            False

        If the matrix is an `n \times 1` matrix with at least one nonzero entry,
        `R` has a single nonzero entry and that entry is a scalar multiple of
        the greatest-common-divisor of the entries of the matrix::

            sage: A = matrix([[x*(x-1)*(x+1)], [x*(x-2)*(x+2)], [x]])
            sage: R = A.reduced_form()
            sage: R
            [x]
            [0]
            [0]

        A zero matrix is already reduced::

            sage: A = matrix(pR, 2, [0,0,0,0])
            sage: A.reduced_form()
            [0 0]
            [0 0]

        In the following example, the original matrix is already reduced, but
        the output is a different matrix: currently this method is an alias for
        :meth:`weak_popov_form`, which is a stronger reduced form::

            sage: R.<x> = QQ['x']
            sage: A = matrix([[x,x,x],[0,0,x]]); A
            [x x x]
            [0 0 x]
            sage: A.is_reduced()
            True
            sage: W = A.reduced_form(); W
            [ x  x  x]
            [-x -x  0]
            sage: W.is_weak_popov()
            True

        The last example shows the usage of the transformation parameter::

            sage: # needs sage.rings.finite_rings
            sage: Fq.<a> = GF(2^3)
            sage: pR.<x> = Fq[]
            sage: A = matrix(pR, [[x^2+a,  x^4+a],
            ....:                 [  x^3,  a*x^4]])
            sage: W, U = A.reduced_form(transformation=True)
            sage: W, U
            (
            [          x^2 + a           x^4 + a]  [1 0]
            [x^3 + a*x^2 + a^2               a^2], [a 1]
            )
            sage: W.is_reduced()
            True
            sage: U*W == A
            True
            sage: U.is_invertible()
            True

        .. SEEALSO::

            :meth:`is_reduced` ,
            :meth:`weak_popov_form` ,
            :meth:`popov_form` ,
            :meth:`hermite_form` .
        """
        return self.weak_popov_form(transformation,
                shifts,
                row_wise,
                False,
                include_zero_vectors)

    def hermite_form(self, include_zero_rows=True, transformation=False):
        """
        Return the Hermite form of this matrix.

        See :meth:`is_hermite` for a definition of Hermite forms. If the input
        is a matrix `A`, then its Hermite form is the unique matrix `H` in Hermite
        form such that `UA = H` for some unimodular matrix `U`.

        INPUT:

        - ``include_zero_rows`` -- boolean (default: ``True``); if ``False``,
          the zero rows in the output matrix are deleted

        - ``transformation`` -- boolean (default: ``False``); if ``True``,
          return the transformation matrix

        OUTPUT:

        - the Hermite normal form `H` of this matrix `A` .

        - (optional) transformation matrix `U` such that `UA = H` .

        EXAMPLES::

            sage: M.<x> = GF(7)[]
            sage: A = matrix(M, 2, 3, [x, 1, 2*x, x, 1+x, 2])
            sage: A.hermite_form()
            [      x       1     2*x]
            [      0       x 5*x + 2]
            sage: A.hermite_form(transformation=True)
            (
            [      x       1     2*x]  [1 0]
            [      0       x 5*x + 2], [6 1]
            )
            sage: A = matrix(M, 2, 3, [x, 1, 2*x, 2*x, 2, 4*x])
            sage: A.hermite_form(transformation=True, include_zero_rows=False)
            ([  x   1 2*x], [0 4])
            sage: H, U = A.hermite_form(transformation=True,
            ....:                       include_zero_rows=True); H, U
            (
            [  x   1 2*x]  [0 4]
            [  0   0   0], [5 1]
            )
            sage: U * A == H
            True
            sage: H, U = A.hermite_form(transformation=True,
            ....:                       include_zero_rows=False)
            sage: U * A
            [  x   1 2*x]
            sage: U * A == H
            True

        .. SEEALSO::

            :meth:`is_hermite` ,
            :meth:`popov_form` .
        """
        A = self.__copy__()
        U = A._hermite_form_euclidean(transformation=transformation,
                                      normalization=lambda p: ~p.lc())
        if not include_zero_rows:
            i = A.nrows() - 1
            while i >= 0 and A.row(i) == 0:
                i -= 1
            A = A[:i+1]
            if transformation:
                U = U[:i+1]

        A.set_immutable()
        if transformation:
            U.set_immutable()

        return (A, U) if transformation else A

    def left_quo_rem(self, B):
        r"""
        Return, if it exists, the quotient and remainder `(Q,R)` such that
        ``self`` is `BQ+R`, where `R` has row degrees less than those of `B`
        entrywise.

        This method directly calls :meth:`right_quo_rem` on transposed
        matrices, and transposes the result. See :meth:`right_quo_rem` for a
        complete documentation and more examples.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: A = Matrix(pR, 3, 2,
            ....:            [[      3*x^3 + 3*x,         2*x^3 + 4],
            ....:             [  3*x^3 + 6*x + 5, 6*x^3 + 5*x^2 + 1],
            ....:             [  2*x^3 + 2*x + 6,   3*x^2 + 2*x + 2]])
            sage: B = Matrix(pR, 3, 3,
            ....:            [[              3,       x + 3,               6],
            ....:             [3*x^3 + 3*x + 1, 4*x^2 + 3*x,   6*x^3 + x + 4],
            ....:             [  4*x^2 + x + 4, 3*x^2 + 4*x, 3*x^2 + 3*x + 2]])
            sage: Q, R = A.left_quo_rem(B); Q, R
            (
            [2*x^2 + 4*x + 6 6*x^2 + 4*x + 1]  [              3               1]
            [    3*x^2 + 5*x   2*x^2 + x + 5]  [              6 5*x^2 + 2*x + 3]
            [    6*x^2 + 3*x 4*x^2 + 6*x + 1], [        2*x + 3         6*x + 3]
            )
            sage: rdegR = R.row_degrees(); rdegB = B.row_degrees()
            sage: A == B*Q+R and all(rdegR[i] < rdegB[i] for i in range(3))
            True

            sage: A[:2,:].left_quo_rem(B)
            Traceback (most recent call last):
            ...
            ValueError: row dimension of self should be the row dimension of
            the input matrix

        Rectangular or rank-deficient matrices are supported but there may be
        no quotient and remainder (unless the matrix has full row rank, see
        :meth:`right_quo_rem`)::

            sage: Q, R = A[:2,:].left_quo_rem(B[:2,:]); Q, R
            (
            [      3*x + 3       2*x + 1]
            [  3*x^2 + 5*x 2*x^2 + x + 5]  [            5             0]
            [            0             0], [4*x^2 + x + 2     4*x^2 + x]
            )
            sage: rdegR = R.row_degrees(); rdegB = B[:2,:].row_degrees()
            sage: A[:2,:] == B[:2,:]*Q+R
            True
            sage: all([rdegR[i] < rdegB[i] for i in range(len(rdegR))])
            True

            sage: A.left_quo_rem(B[:,:2])
            Traceback (most recent call last):
            ...
            ValueError: division of these matrices does not admit a remainder
            with the required degree property

        .. SEEALSO::

            :meth:`right_quo_rem` ,
            :meth:`reduce` .
        """
        if self.nrows() != B.nrows():
            raise ValueError("row dimension of self should be the"
                             " row dimension of the input matrix")
        (Q,R) = self.T.right_quo_rem(B.T)
        return (Q.T,R.T)

    def right_quo_rem(self, B):
        r"""
        Return, if it exists, the quotient and remainder `(Q,R)` such that
        ``self`` is `QB+R`, where `R` has column degrees less than those of `B`
        entrywise.

        If ``self`` is a `k \times m` polynomial matrix (written `A` below),
        and the input `B` is an `m \times m` polynomial matrix in column
        reduced form, then `(Q,R)` is unique. Both `Q` and `R` have dimensions
        `k \times m`. In this case, this method implements Newton iteration of
        a reversed polynomial matrix `B`, generalizing to matrices the fast
        division of univariate polynomials.

        If `A` is a `k \times n` polynomial matrix, and the input `B` is an `m
        \times n` polynomial matrix such that `B` has full column rank, or more
        generally such that the matrix equation `A = XB` has a rational
        solution, then there exists such a `(Q,R)` but it may not be unique;
        the algorithm returns one such quotient and remainder. Here `Q` is `k
        \times m` and `R` is `k \times n`. In this case, this method follows
        the folklore approach based on solving the matrix equation `A = XB` and
        splitting `X` into its polynomial part and proper fraction part.

        Finally, if the matrix equation `A = XB` has no rational solution, this
        method computes the normal form `R` and quotient `Q` of the rows of `A`
        with respect to the row space of `B` (see :meth:`reduce`). Doing this
        for a well-chosen shift ensures that either `R` does not have column
        degrees less than those of `B`, and then there is no valid quotient and
        remainder, or it does satisfy this degree constraint, and then this `R`
        can be returned as a remainder along with the quotient `Q`.

        A :exc:`ValueError` is raised if the dimensions of ``self`` and `B` are
        not conformal, or if there exists no quotient and remainder.

        EXAMPLES:

        Case where `B` is a square, column reduced matrix::

            sage: pR.<x> = GF(7)[]
            sage: A = Matrix(pR, 2, 3,
            ....:     [[3*x^3 + 3*x, 3*x^3 + 6*x + 5,   2*x^3 + 2*x + 6],
            ....:      [2*x^3 + 4,   6*x^3 + 5*x^2 + 1, 3*x^2 + 2*x + 2]])
            sage: B = Matrix(pR, 3, 3,
            ....:     [[4*x^2 + 3*x + 3, 3*x^2 + 3*x + 1,   4*x^2 + x + 4],
            ....:      [6*x^2 + 2*x + 3,     4*x^2 + 3*x,     3*x^2 + 4*x],
            ....:      [5*x^2 + 3*x + 6,   6*x^2 + x + 4, 3*x^2 + 3*x + 2]])
            sage: B.is_reduced(row_wise=False)
            True
            sage: Q, R = A.right_quo_rem(B); Q, R
            (
            [    4*x   x + 2 6*x + 1]  [  x + 2 6*x + 1 5*x + 4]
            [4*x + 3   x + 6 3*x + 4], [4*x + 2 2*x + 3 4*x + 3]
            )
            sage: A == Q*B+R and R.degree() < 2
            True
            sage: A[:,:2].right_quo_rem(B)
            Traceback (most recent call last):
            ...
            ValueError: column dimension of self should be the column dimension
            of the input matrix

            sage: B = Matrix(pR, 3, 3,
            ....:     [[3,     3*x^3 + 3*x + 1, 4*x^2 + x + 4],
            ....:      [x + 3, 4*x^2 + 3*x,     3*x^2 + 4*x],
            ....:      [6,     6*x^3 + x + 4,   3*x^2 + 3*x + 2]])
            sage: B.is_reduced(row_wise=False)
            True
            sage: Q, R = A.right_quo_rem(B); Q, R
            (
            [2*x^2 + 4*x + 6     3*x^2 + 5*x     6*x^2 + 3*x]
            [6*x^2 + 4*x + 1   2*x^2 + x + 5 4*x^2 + 6*x + 1],
            <BLANKLINE>
            [              3               6         2*x + 3]
            [              1 5*x^2 + 2*x + 3         6*x + 3]
            )
            sage: cdegR = R.column_degrees(); cdegB = B.column_degrees()
            sage: A == Q*B+R and all([cdegR[i] < cdegB[i] for i in range(3)])
            True

        With a nonsingular but also non-reduced matrix, there exists a
        solution, but it might not be unique::

            sage: B = Matrix(pR, 3, 3,
            ....:     [[              5,               0, 2*x + 6],
            ....:      [            4*x, 3*x^2 + 4*x + 5,   x + 1],
            ....:      [3*x^2 + 5*x + 2, 6*x^3 + 4*x + 6,       3]])
            sage: B.det() != 0 and (not B.is_reduced(row_wise=False))
            True
            sage: Q, R = A.right_quo_rem(B); Q, R
            (
            [    6*x^2 + 3*x 4*x^2 + 3*x + 1         5*x + 1]
            [  x^2 + 5*x + 5 5*x^2 + 3*x + 5           x + 2],
            <BLANKLINE>
            [      4*x + 5 x^2 + 2*x + 1             2]
            [      6*x + 3     5*x^2 + 6             3]
            )
            sage: cdegR = R.column_degrees(); cdegB = B.column_degrees()
            sage: A == Q*B+R and all(cdegR[i] < cdegB[i] for i in range(3))
            True

            sage: Q2 = Matrix(pR, 2, 3,
            ....:      [[6*x^2 + 3*x + 1, 4*x^2 + 3*x + 6, 5*x + 1],
            ....:       [  x^2 + 5*x + 3, 5*x^2 + 3*x + 2,   x + 2]])
            sage: R2 = Matrix(pR, 2, 3,
            ....:      [[    5*x, 3*x + 4, 5],
            ....:       [4*x + 6,     5*x, 4]])
            sage: A == Q2*B + R2
            True

        The same remark holds more generally for full column rank matrices:
        there exists a solution, but it might not be unique. However, for all
        other cases (rank-deficient matrix `B` or matrix `B` having strictly
        fewer rows than columns) there may be no solution::

            sage: C = B.stack(B[1,:] + B[2,:])  # 4 x 3, full column rank
            sage: Q, R = A.right_quo_rem(C); Q, R
            (
            [    6*x^2 + 3*x 4*x^2 + 3*x + 1         5*x + 1               0]
            [  x^2 + 5*x + 5 5*x^2 + 3*x + 5           x + 2               0],
            <BLANKLINE>
            [      4*x + 5 x^2 + 2*x + 1             2]
            [      6*x + 3     5*x^2 + 6             3]
            )

            sage: A.right_quo_rem(B[:2,:])  # matrix 2 x 3, full row rank               # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: division of these matrices does not admit a remainder
            with the required degree property
            sage: D = copy(B); D[2,:] = B[0,:]+B[1,:]  # square, singular
            sage: A.right_quo_rem(D)
            Traceback (most recent call last):
            ...
            ValueError: division of these matrices does not admit a remainder
            with the required degree property

        In the latter case (rank-deficient or strictly fewer rows than columns,
        with no solution to `A = XB`), there might still be a quotient and
        remainder, in which case this method will find it via normal form
        computation::

            sage: B = Matrix(pR, 1, 2, [[x, x]])
            sage: A = Matrix(pR, 1, 2, [[x, x+2]])
            sage: A.right_quo_rem(B)
            ([1], [0 2])
            sage: A == 1*B + Matrix([[0,2]])
            True

        .. SEEALSO::

            :meth:`left_quo_rem` ,
            :meth:`reduce` .
        """
        if self.ncols() != B.ncols():
            raise ValueError("column dimension of self should be the"
                             " column dimension of the input matrix")
        if B.is_square() and \
           B.is_reduced(row_wise=False,include_zero_vectors=False):
            # case of B column reduced (without zero columns):
            # direct matrix version of univariate polynomial quo_rem
            return self._right_quo_rem_reduced(B)
        try:
            # more generally, method via solving A = XB over rationals
            # (always possible if B has full column rank, otherwise might still
            # work if this matrix equation has a solution "by luck")
            return self._right_quo_rem_solve(B)
        except ValueError:
            # more generally, any B
            # (compute normal form w.r.t a well-chosen shift and check degrees
            # are as expected)
            s = [-d for d in B.column_degrees()]
            (Q,R) = self.reduce(B,shifts=s,return_quotient=True)
            cdeg = R.column_degrees()
            if all([cdeg[i] + s[i] < 0 for i in range(B.ncols())]):
                return (Q, R)
            raise ValueError("division of these matrices does not admit a "
                             "remainder with the required degree property")

    def _right_quo_rem_reduced(self, B):
        r"""
        If ``self`` is a `k x m` polynomial matrix (written `A` below), and the
        input `B` is an `m x m` polynomial matrix in column reduced form, this
        computes the unique couple `(Q,R)` of `k x m` polynomial matrices such
        that `A = QB + R`, with the column degrees of `R` entrywise less than
        those of `B`.

        The fact that `B` is in column reduced form is required, and not
        checked. Reference: we follow the folklore algorithm which generalizes
        the fast division of univariate polynomials; precisely we implement
        [Neiger-Vu, ISSAC 2017, Algorithm 1] .

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: A = Matrix(pR, 2, 3,
            ....:     [[3*x^3 + 3*x, 3*x^3 + 6*x + 5,   2*x^3 + 2*x + 6],
            ....:      [2*x^3 + 4,   6*x^3 + 5*x^2 + 1, 3*x^2 + 2*x + 2]])
            sage: B = Matrix(pR, 3, 3,
            ....:     [[4*x^2 + 3*x + 3, 3*x^2 + 3*x + 1,   4*x^2 + x + 4],
            ....:      [6*x^2 + 2*x + 3,     4*x^2 + 3*x,     3*x^2 + 4*x],
            ....:      [5*x^2 + 3*x + 6,   6*x^2 + x + 4, 3*x^2 + 3*x + 2]])
            sage: B.is_reduced(row_wise=False)
            True
            sage: Q, R = A._right_quo_rem_reduced(B); Q, R
            (
            [    4*x   x + 2 6*x + 1]  [  x + 2 6*x + 1 5*x + 4]
            [4*x + 3   x + 6 3*x + 4], [4*x + 2 2*x + 3 4*x + 3]
            )
            sage: A == Q*B+R and R.degree() < 2
            True

            sage: B = Matrix(pR, 3, 3,
            ....:     [[4*x + 3*x + 3, 3*x^3 + 3*x + 1,   4*x^2 + x + 4],
            ....:      [6*x + 2*x + 3,     4*x^2 + 3*x,     3*x^2 + 4*x],
            ....:      [6,             6*x^3 + x + 4,   3*x^2 + 3*x + 2]])
            sage: B.is_reduced(row_wise=False)
            True
            sage: Q, R = A._right_quo_rem_reduced(B); Q, R
            (
            [2*x^2 + 4*x + 6     3*x^2 + 5*x     6*x^2 + 3*x]
            [6*x^2 + 4*x + 1   2*x^2 + x + 5 4*x^2 + 6*x + 1],
            <BLANKLINE>
            [              3               6         2*x + 3]
            [              1 5*x^2 + 2*x + 3         6*x + 3]
            )
            sage: cdegR = R.column_degrees(); cdegB = B.column_degrees()
            sage: A == Q*B+R and all(cdegR[i] < cdegB[i] for i in range(3))
            True
        """
        # Step 0: find parameter d  (delta in above reference)
        cdegA = self.column_degrees() # zero columns of A --> entries -1 in cdegA
        cdeg = B.column_degrees()  # all nonnegative since column reduced
        d = max([cdegA[i]-cdeg[i]+1 for i in range(B.nrows())])
        if d<=0: # A already reduced modulo B, quotient is zero
            return (self.parent().zero().__copy__(), self)
        # Step 1: reverse input matrices
        # Brev = B(1/x) diag(x^(cdeg[i]))
        # Arev = A(1/x) diag(x^(d+cdeg[i]-1))
        Brev = B.reverse(degree=cdeg, row_wise=False)
        Arev = self.reverse(degree=[d+c-1 for c in cdeg], row_wise=False)
        # Step 2: compute quotient
        # compute Qrev = Arev Brev^{-1} mod x^d
        # then quotient is the reverse Q = x^(d-1) Qrev(1/x)
        Q = Brev.solve_left_series_trunc(Arev, d).reverse(degree=d-1)
        # Step 3: deduce remainder and return
        R = self - Q*B
        return Q,R

    def _right_quo_rem_solve(self, B):
        r"""
        If ``self`` is a `k x n` polynomial matrix (written `A` below), and the
        input `B` is an `m x n` polynomial matrix with full column rank, this
        computes a couple `(Q,R)` polynomial matrices, of sizes `k x m` and `k
        x n` respectively, such that `A = QB + R`, with the column degrees of
        `R` entrywise less than those of `B`.

        Algorithm: this method follows the folklore approach based on solving
        the matrix equation `A = XB` and separating `X` into its polynomial
        part and proper fraction part.

        The fact that `B` is in column reduced form is not required, and not
        checked. If `B` does not have full column rank, two situations can
        arise. Either the above matrix equation has a solution `X`, which
        implies the existence of a quotient and remainder as described above,
        and such a quotient and remainder is returned by the method. Or this
        matrix equation has no solution and this method fails: this raises
        :exc:`ValueError`; however this is not a proof that there is no valid
        division with remainder (see the last example below).

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: A = Matrix(pR, 2, 3,
            ....:     [[3*x^3 + 3*x, 3*x^3 + 6*x + 5,   2*x^3 + 2*x + 6],
            ....:      [2*x^3 + 4,   6*x^3 + 5*x^2 + 1, 3*x^2 + 2*x + 2]])
            sage: B = Matrix(pR, 3, 3,
            ....:     [[4*x + 3*x + 3, 3*x^3 + 3*x + 1,   4*x^2 + x + 4],
            ....:      [6*x + 2*x + 3,     4*x^2 + 3*x,     3*x^2 + 4*x],
            ....:      [6,             6*x^3 + x + 4,   3*x^2 + 3*x + 2]])
            sage: Q, R = A._right_quo_rem_solve(B); Q, R
            (
            [2*x^2 + 4*x + 6     3*x^2 + 5*x     6*x^2 + 3*x]
            [6*x^2 + 4*x + 1   2*x^2 + x + 5 4*x^2 + 6*x + 1],
            <BLANKLINE>
            [              3               6         2*x + 3]
            [              1 5*x^2 + 2*x + 3         6*x + 3]
            )
            sage: B.is_reduced(row_wise=False)
            True
            sage: cdegR = R.column_degrees(); cdegB = B.column_degrees()
            sage: A == Q*B+R and all([cdegR[i] < cdegB[i] for i in range(3)])
            True

        With a nonsingular but also non-reduced matrix, there exists a solution
        and one is found by this method, but it might not be unique::

            sage: B = Matrix(pR, 3, 3,
            ....:     [[              5,               0, 2*x + 6],
            ....:      [            4*x, 3*x^2 + 4*x + 5,   x + 1],
            ....:      [3*x^2 + 5*x + 2, 6*x^3 + 4*x + 6,       3]])
            sage: B.det() != 0 and not B.is_reduced(row_wise=False)
            True
            sage: Q, R = A._right_quo_rem_solve(B); Q, R
            (
            [    6*x^2 + 3*x 4*x^2 + 3*x + 1         5*x + 1]
            [  x^2 + 5*x + 5 5*x^2 + 3*x + 5           x + 2],
            <BLANKLINE>
            [      4*x + 5 x^2 + 2*x + 1             2]
            [      6*x + 3     5*x^2 + 6             3]
            )
            sage: cdegR = R.column_degrees(); cdegB = B.column_degrees()
            sage: A == Q*B+R and all(cdegR[i] < cdegB[i] for i in range(3))
            True

            sage: Q2 = Matrix(pR, 2, 3,
            ....:      [[6*x^2 + 3*x + 1, 4*x^2 + 3*x + 6, 5*x + 1],
            ....:       [  x^2 + 5*x + 3, 5*x^2 + 3*x + 2,   x + 2]])
            sage: R2 = Matrix(pR, 2, 3,
            ....:      [[    5*x, 3*x + 4, 5],
            ....:       [4*x + 6,     5*x, 4]])
            sage: A == Q2*B + R2
            True

        The same remark holds more generally for full column rank matrices:
        there exists a solution and this method will find one. However for all
        other cases (rank-deficient or strictly fewer rows than columns) there
        might be no solution::

            sage: C = B.stack(B[1,:] + B[2,:])  # 4 x 3, full column rank
            sage: Q, R = A._right_quo_rem_solve(C); Q, R
            (
            [    6*x^2 + 3*x 4*x^2 + 3*x + 1         5*x + 1               0]
            [  x^2 + 5*x + 5 5*x^2 + 3*x + 5           x + 2               0],
            <BLANKLINE>
            [      4*x + 5 x^2 + 2*x + 1             2]
            [      6*x + 3     5*x^2 + 6             3]
            )

            sage: A._right_quo_rem_solve(B[:2,:])  # 2 x 3, full row rank
            Traceback (most recent call last):
            ...
            ValueError: dividing via system solving yields no solution
            sage: D = copy(B); D[2,:] = B[0,:]+B[1,:]  # square, singular
            sage: A._right_quo_rem_solve(D)
            Traceback (most recent call last):
            ...
            ValueError: dividing via system solving yields no solution

        In the latter case (rank-deficient or strictly fewer rows than
        columns), even when there is a solution, this method might not find
        it::

            sage: B = Matrix(pR, 1, 2, [[x, x]])
            sage: A = Matrix(pR, 1, 2, [[x, x+2]])
            sage: A == 1*B + Matrix([[0,2]])    # a valid quo_rem
            True
            sage: A._right_quo_rem_solve(B)
            Traceback (most recent call last):
            ...
            ValueError: dividing via system solving yields no solution
        """
        k = self.nrows()
        m = B.nrows()
        # find rational Q such that QB = A
        try:
            X = B.solve_left(self)
        except ValueError:
            raise ValueError("dividing via system solving yields no solution")
        from sage.arith.functions import lcm
        f = lcm([X[i,j].denom() for j in range(m) for i in range(k)])
        # Write X = Q + R/r with Q and R having polynomial entries;
        # keep only the Q part
        for i in range(k):
            for j in range(m):
                X[i,j] = ((X[i,j].numer() * f) // X[i,j].denom()) // f
        Q = X.change_ring(self.base_ring())
        R = self - Q*B
        return (Q,R)

    def reduce(self, B, shifts=None, row_wise=True, return_quotient=False):
        r"""
        Reduce ``self``, i.e. compute its normal form, modulo the row space of
        `B` with respect to ``shifts``.

        If ``self`` is a `k \times n` polynomial matrix (written `A` below),
        and the input `B` is an `m \times n` polynomial matrix, this computes
        the normal form `R` of `A` with respect the row space of `B` and the
        monomial order defined by ``shifts`` (written `s` below). This means
        that the `i` th row of `R` is equal to the `i` th row of `A` up to
        addition of an element in the row space of `B`, and if `J =
        (j_1,\ldots,j_r)` are the `s`-leading positions of the `s`-Popov form
        `P` of `A`, then the submatrix `R_{*,J}` (submatrix of `R` formed by
        its columns in `J`) has column degrees smaller entrywise than the
        column degrees of `P_{*,J}`.

        If the option ``row_wise`` is set to ``False``, the same operation is
        performed, but with everything considered column-wise: column space of
        `B`, `i` th column of `R` and `A`, column-wise `s`-leading positions
        and `s`-Popov form, and submatrices `R_{J,*}` and `P_{J,*}`.

        The operation above can be seen as a matrix generalization of division
        with remainder for univariate polynomials. If the option
        ``return_quotient`` is set to ``True``, this method returns both the
        normal form `R` and a quotient matrix `Q` such that `A = QB + R` (or `A
        = BQ + R` if ``row_wise`` is ``False``). Whereas the remainder is
        unique, this quotient is not unique unless `B` has a trivial left
        kernel i.e. has full row rank (or right kernel, full column rank if
        ``row_wise`` is ``False``).

        This method checks whether `B` is in `s`-Popov form, and if not,
        computes the `s`-Popov form `P` of `B`, which can take some time.
        Therefore, if `P` is already known or is to be re-used, this method
        should be called directly with `P`, yielding the same normal form `R`
        since `P` and `B` have the same row space (or column space, if
        ``row_wise`` is ``False``).

        A :exc:`ValueError` is raised if the dimensions of the shifts and/or of
        the matrices are not conformal.

        INPUT:

        - ``B`` -- polynomial matrix

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        - ``row_wise`` -- boolean (default: ``True``); ``True`` if
          working row-wise (see the class description)

        - ``return_quotient`` -- (default: ``False``) if this
          is ``True``, the quotient will be returned as well

        OUTPUT: a polynomial matrix if ``return_quotient=False``, two
        polynomial matrices otherwise.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: B = Matrix(pR, [
            ....:     [      6*x+4,       5*x^3+5*x,       6*x^2+2*x+2],
            ....:     [4*x^2+5*x+2, x^4+5*x^2+2*x+4, 4*x^3+6*x^2+6*x+5]])
            sage: A = Matrix(pR, 1, 3, [
            ....:     [3*x^4+3*x^3+4*x^2+5*x+1, x^4+x^3+5*x^2+4*x+4, 4*x^4+2*x^3+x]])

            sage: Q, R = A.reduce(B,return_quotient=True); R
            [3*x^4 + 3*x^3 + 4*x + 3                 2*x + 2                 2*x + 6]
            sage: A == Q*B + R
            True
            sage: P = B.popov_form(); P.leading_positions(return_degree=True)
            ([1, 2], [2, 2])
            sage: R.degree_matrix()
            [4 1 1]
            sage: A.reduce(P) == R
            True
            sage: A.reduce(P[:,:2])
            Traceback (most recent call last):
            ...
            ValueError: column dimension of self should be the column
            dimension of the input matrix

        Demonstrating shifts::

            sage: Qs, Rs = A.reduce(B, shifts=[0,2,4], return_quotient=True); Rs
            [3*x^4 + 3*x^3 + 6*x + 2             4*x^3 + 5*x                       0]
            sage: A == Qs*B + Rs
            True
            sage: Ps = B.popov_form(shifts=[0,2,4])
            sage: Ps.leading_positions(shifts=[0,2,4], return_degree=True)
            ([1, 2], [4, 0])
            sage: Rs.degree_matrix()
            [ 4  3 -1]
            sage: A.reduce(Ps, shifts=[0,2,4]) == Rs
            True

        If ``return_quotient`` is ``False``, only the normal form is returned::

            sage: R == A.reduce(B) and Rs == A.reduce(B, shifts=[0,2,4])
            True

        Demonstrating column-wise normal forms, with a matrix `A` which has
        several columns, and a matrix `B` which does not have full column rank
        (its column-wise Popov form has a zero column)::

            sage: A = Matrix(pR, 2, 2,
            ....:     [[5*x^3 + 2*x^2 + 4*x + 1,           x^3 + 4*x + 4],
            ....:      [2*x^3 + 5*x^2 + 2*x + 4,         2*x^3 + 3*x + 2]])
            sage: (Q,R) = A.reduce(B,row_wise=False, return_quotient=True); R
            [0 3]
            [0 0]
            sage: A == B*Q + R
            True
            sage: P = B.popov_form(row_wise=False); P
            [x + 2     6     0]
            [    0     1     0]
            sage: P.leading_positions(row_wise=False, return_degree=True)
            ([0, 1, -1], [1, 0, -1])
            sage: R.degree_matrix()
            [-1  0]
            [-1 -1]

        .. SEEALSO::

            :meth:`left_quo_rem` ,
            :meth:`right_quo_rem` .
        """
        if row_wise and self.ncols() != B.ncols():
            raise ValueError("column dimension of self should be the"
                             " column dimension of the input matrix")
        if not row_wise and (self.nrows() != B.nrows()):
            raise ValueError("row dimension of self should be the"
                             " row dimension of the input matrix")
        # note: is_popov calls B._check_shift_dimension(shifts,row_wise)
        # --> no need to check here again
        if B.is_popov(shifts,row_wise,False,False):
            lpos = B.leading_positions(shifts=shifts,row_wise=row_wise)
            set_lpos = set(lpos) # to make things faster for huge matrices..
            if row_wise:
                non_lpos = [j for j in range(B.ncols()) if j not in set_lpos]
                # A_{*,J} = Q B_{*,J} + R0
                Q,R0 = self[:,lpos]._right_quo_rem_reduced(B[:,lpos])
                # other columns are given by A_{*,not J} - Q B_{*, not J}
                R = self.parent().zero().__copy__()
                R[:,lpos] = R0
                R[:,non_lpos] = self[:,non_lpos] - Q * B[:,non_lpos]
                return (Q,R) if return_quotient else R
            else:
                non_lpos = [i for i in range(B.nrows()) if i not in set_lpos]
                # A_{I,*} = B_{I,*} Q + R0
                Q,R0 = self[lpos,:].T._right_quo_rem_reduced(B[lpos,:].T)
                Q = Q.T
                R0 = R0.T
                # other columns are given by A_{not I,*} - B_{not I,*} Q
                R = self.parent().zero().__copy__()
                R[lpos,:] = R0
                R[non_lpos,:] = self[non_lpos,:] - B[non_lpos,:] * Q
                return (Q,R) if return_quotient else R
        elif return_quotient:
            P,U = B.popov_form(True,shifts,row_wise,False)
            Q,R = self.reduce(P,shifts,row_wise,True)
            # row-wise: UB = P and A = QP + R ==> A = QUB + R
            # --> careful: the last rows of U may correspond to zero rows of P,
            # which have been discarded... so not exactly UB = P
            return (Q*U[:P.nrows(),:], R) if row_wise \
                    else (U[:,:P.ncols()] * Q, R)
        else:
            P = B.popov_form(False,shifts,row_wise,False)
            return self.reduce(P,shifts,row_wise,False)

    def is_minimal_approximant_basis(self,
            pmat,
            order,
            shifts=None,
            row_wise=True,
            normal_form=False):
        r"""
        Return whether this matrix is an approximant basis in
        ``shifts``-ordered weak Popov form for the polynomial matrix ``pmat``
        at order ``order``.

        If ``normal_form`` is ``True``, then the polynomial matrix must
        furthermore be in ``shifts``-Popov form. An error is raised if the
        input dimensions are not sound. If a single integer is provided for
        ``order``, then it is interpreted as a list of repeated integers with
        this value. (See :meth:`minimal_approximant_basis` for definitions and
        more details.)

        INPUT:

        - ``pmat`` -- a polynomial matrix

        - ``order`` -- list of positive integers, or a positive integer

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        - ``row_wise`` -- boolean (default: ``True``); if ``True``
          then the basis considered row-wise and operates on the left of
          ``pmat``. Otherwise it is column-wise and operates on the right of
          ``pmat``.

        - ``normal_form`` -- boolean (default: ``False``); if
          ``True`` then checks for a basis in ``shifts``-Popov form

        OUTPUT: boolean

        ALGORITHM:

        Verification that the matrix is formed by approximants is done via a
        truncated matrix product; verification that the matrix is square,
        nonsingular and in shifted weak Popov form is done via
        :meth:`is_weak_popov`; verification that the matrix generates the
        module of approximants is done via the characterization in Theorem 2.1
        of [GN2018]_ .

        EXAMPLES::

            sage: pR.<x> = GF(97)[]

        We consider the following example from [Arne Storjohann, Notes on
        computing minimal approximant bases, 2006]::

            sage: order = 8; shifts = [1,1,0,0,0]
            sage: pmat = Matrix(pR, 5, 1, [
            ....:     pR([35,  0, 41, 87,  3, 42, 22, 90]),
            ....:     pR([80, 15, 62, 87, 14, 93, 24,  0]),
            ....:     pR([42, 57, 90, 87, 22, 80, 71, 53]),
            ....:     pR([37, 72, 74,  6,  5, 75, 23, 47]),
            ....:     pR([36, 10, 74,  1, 29, 44, 87, 74])])
            sage: appbas = Matrix(pR, [
            ....:     [x+47,   57, 58*x+44,     9*x+23,      93*x+76],
            ....:     [  15, x+18, 52*x+23,     15*x+58,     93*x+88],
            ....:     [  17,   86, x^2+77*x+16, 76*x+29,     90*x+78],
            ....:     [  44,   36, 3*x+42,      x^2+50*x+26, 85*x+44],
            ....:     [   2,   22, 54*x+94,     73*x+24,     x^2+2*x+25]])
            sage: appbas.is_minimal_approximant_basis(                                  # needs sage.libs.pari
            ....:     pmat, order, shifts, row_wise=True, normal_form=True)
            True

        The matrix `x^8 \mathrm{Id}_5` is square, nonsingular, in Popov form,
        and its rows are approximants for ``pmat`` at order 8. However, it is
        not an approximant basis since its rows generate a module strictly
        contained in the set of approximants for ``pmat`` at order 8::

            sage: M = x^8 * Matrix.identity(pR, 5)
            sage: M.is_minimal_approximant_basis(pmat, 8)                               # needs sage.libs.pari
            False

        Since ``pmat`` is a single column, with nonzero constant coefficient,
        its column-wise approximant bases at order 8 are all `1\times 1`
        matrices `[c x^8]` for some nonzero field element `c`::

            sage: M = Matrix(pR, [x^8])
            sage: M.is_minimal_approximant_basis(
            ....:     pmat, 8, row_wise=False, normal_form=True)
            True

        Exceptions are raised if input dimensions are not sound::

            sage: appbas.is_minimal_approximant_basis(pmat, [8,8], shifts)
            Traceback (most recent call last):
            ...
            ValueError: order length should be the column dimension
                        of the input matrix

            sage: appbas.is_minimal_approximant_basis(
            ....:     pmat, order, shifts, row_wise=False)
            Traceback (most recent call last):
            ...
            ValueError: shifts length should be the column dimension
                        of the input matrix

            sage: Matrix(pR, [x^8]).is_minimal_approximant_basis(pmat, 8)
            Traceback (most recent call last):
            ...
            ValueError: column dimension should be the row dimension of the
            input matrix

        .. SEEALSO::

            :meth:`minimal_approximant_basis` .
        """
        m = pmat.nrows()
        n = pmat.ncols()

        # set default shifts / check shifts dimension
        if shifts is None:
            shifts = [0] * m if row_wise else [0] * n
        elif row_wise and len(shifts) != m:
            raise ValueError('shifts length should be the row dimension of'
                             ' the input matrix')
        elif (not row_wise) and len(shifts) != n:
            raise ValueError('shifts length should be the column dimension'
                             ' of the input matrix')

        # set default order / check order dimension
        if not isinstance(order,list):
            order = [order]*n if row_wise else [order]*m

        if row_wise and len(order) != n:
            raise ValueError("order length should be the column dimension"
                             " of the input matrix")
        elif (not row_wise) and len(order) != m:
            raise ValueError("order length should be the row dimension of"
                             " the input matrix")

        # raise an error if self does not have the right dimension
        if row_wise and self.ncols() != m:
            raise ValueError("column dimension should be the row dimension"
                             " of the input matrix")
        elif (not row_wise) and self.nrows() != n:
            raise ValueError("row dimension should be the column dimension"
                             " of the input matrix")

        # check square
        if not self.is_square():
            return False
        # check nonsingular and shifts-(ordered weak) Popov form
        if normal_form and (not self.is_popov(shifts, row_wise, False, False)):
            return False
        if (not normal_form) and (not self.is_weak_popov(shifts, row_wise, True, False)):
            return False

        # check that self is a basis of the set of approximants
        if row_wise:
            # check that self * pmat is 0 bmod x^order
            # and compute certificate matrix ``cert_mat`` which is
            # the constant term of (self * pmat) * x^(-order)
            residual = self * pmat
            if not residual.truncate(order,row_wise=False).is_zero():
                return False
            cert_mat = residual.coefficient_matrix(order,row_wise=False)

            # check that self generates the set of approximants
            # 1/ determinant of self should be a monomial c*x^d,
            # with d the sum of pivot degrees
            d = sum([self[i,i].degree() for i in range(m)])
            X = self.base_ring().gen()
            if self.determinant() != (self(1).determinant() * X**d):
                return False
            # 2/ the m x (m+n) constant matrix [self(0) | cert_mat] should have
            # full rank, that is, rank m
            from sage.matrix.constructor import block_matrix
            if block_matrix([[self.constant_matrix(), cert_mat]]).rank() < m:
                return False

        else:
            # check that pmat * self is 0 bmod x^order
            # and compute certificate matrix ``cert_mat`` which is
            # the constant term of x^(-order) * (pmat * self)
            residual = pmat * self
            if not residual.truncate(order).is_zero():
                return False
            cert_mat = residual.coefficient_matrix(order)

            # check that self generates the set of approximants
            # 1/ determinant of self should be a monomial c*x^d,
            # with d the sum of pivot degrees
            d = sum([self[i,i].degree() for i in range(n)])
            X = self.base_ring().gen()
            if self.determinant() != (self(1).determinant() * X**d):
                return False
            # 2/ the (m+n) x n constant matrix [self(0).T | cert_mat.T].T
            # should have full rank, that is, rank n
            from sage.matrix.constructor import block_matrix
            if block_matrix([[self.constant_matrix()], [cert_mat]]).rank() < n:
                return False

        return True

    def minimal_approximant_basis(self,
            order,
            shifts=None,
            row_wise=True,
            normal_form=False):
        r"""
        Return an approximant basis in ``shifts``-ordered weak Popov form for
        this polynomial matrix at order ``order``.

        Assuming we work row-wise, if `F` is an `m \times n` polynomial matrix
        and `(d_0,\ldots,d_{n-1})` are positive integers, then an approximant
        basis for `F` at order `(d_0,\ldots,d_{n-1})` is a polynomial matrix
        whose rows form a basis of the module of approximants for `F` at order
        `(d_0,\ldots,d_{n-1})`. The latter approximants are the polynomial
        vectors `p` of size `m` such that the column `j` of `p F` has valuation
        at least `d_j`, for all `0 \le j \le n-1`.

        If ``normal_form`` is ``True``, then the output basis `P` is
        furthermore in ``shifts``-Popov form. By default, `P` is considered
        row-wise, that is, its rows are left-approximants for ``self``; if
        ``row_wise`` is ``False`` then its columns are right-approximants for
        ``self``. It is guaranteed that the degree of the output basis is at
        most the maximum of the entries of ``order``, independently of
        ``shifts``.

        An error is raised if the input dimensions are not sound: if working
        row-wise (resp. column-wise), the length of ``order`` must be the
        number of columns (resp. rows) of ``self``, while the length of
        ``shifts`` must be the number of rows (resp. columns) of ``self``.

        If a single integer is provided for ``order``, then it is converted
        into a list of repeated integers with this value.

        INPUT:

        - ``order`` -- list of positive integers, or a positive integer

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        - ``row_wise`` -- boolean (default: ``True``); if ``True``
          then the output basis is considered row-wise and operates on the left
          of ``self``. Otherwise it is column-wise and operates on the right
          of ``self``.

        - ``normal_form`` -- boolean (default: ``False``); if
          ``True`` then the output basis is in ``shifts``-Popov form

        OUTPUT: a polynomial matrix

        ALGORITHM:

        The implementation is inspired from the iterative algorithms described
        in [VBB1992]_ and [BL1994]_ ; for obtaining the normal form, it relies
        directly on Lemmas 3.3 and 4.1 in [JNSV2016]_ .

        EXAMPLES::

            sage: pR.<x> = GF(7)[]

            sage: order = [4, 3]; shifts = [-1, 2, 0]
            sage: F = Matrix(pR, [[5*x^3 + 4*x^2 + 4*x + 6, 5*x^2 + 4*x + 1],
            ....:                 [        2*x^2 + 2*x + 3, 6*x^2 + 6*x + 3],
            ....:                 [4*x^3         +   x + 1, 4*x^2 + 2*x + 3]])
            sage: P = F.minimal_approximant_basis(order, shifts)
            sage: P.is_minimal_approximant_basis(F, order, shifts)
            True

        By default, the computed basis is not required to be in normal form
        (and will not be except in rare special cases)::

            sage: P.is_minimal_approximant_basis(F, order, shifts,
            ....:                                normal_form=True)
            False
            sage: P = F.minimal_approximant_basis(order, shifts,
            ....:                                 normal_form=True)
            sage: P.is_minimal_approximant_basis(F, order, shifts,
            ....:                                normal_form=True)
            True

        If shifts are not specified, they are chosen as uniform `[0,\ldots,0]`
        by default. Besides, if the orders are all the same, one can rather
        give a single integer::

            sage: (F.minimal_approximant_basis(3) ==
            ....:  F.minimal_approximant_basis([3,3], shifts=None))
            True

        One can work column-wise by specifying ``row_wise=False``::

            sage: P = F.minimal_approximant_basis([5,2,2], [0,1],
            ....:                                 row_wise=False)
            sage: P.is_minimal_approximant_basis(F, [5,2,2], shifts=[0,1],
            ....:                                row_wise=False)
            True
            sage: (F.minimal_approximant_basis(3, row_wise=True) ==
            ....:  F.transpose().minimal_approximant_basis(
            ....:      3, row_wise=False).transpose())
            True

        Errors are raised if the input dimensions are not sound::

            sage: P = F.minimal_approximant_basis([4], shifts)
            Traceback (most recent call last):
            ...
            ValueError: order length should be the column dimension

            sage: P = F.minimal_approximant_basis(order, [0,0,0,0])
            Traceback (most recent call last):
            ...
            ValueError: shifts length should be the row dimension

        An error is raised if order does not contain only positive integers::

            sage: P = F.minimal_approximant_basis([1,0], shifts)
            Traceback (most recent call last):
            ...
            ValueError: order should consist of positive integers
        """
        m = self.nrows()
        n = self.ncols()

        # set default shifts / check shifts dimension
        if shifts is None:
            shifts = [0] * m if row_wise else [0] * n
        elif row_wise and len(shifts) != m:
            raise ValueError('shifts length should be the row dimension')
        elif (not row_wise) and len(shifts) != n:
            raise ValueError('shifts length should be the column dimension')

        # set default order / check order dimension
        if not isinstance(order,list):
            order = [order]*n if row_wise else [order]*m

        if row_wise and len(order) != n:
            raise ValueError("order length should be the column dimension")
        elif (not row_wise) and len(order) != m:
            raise ValueError("order length should be the row dimension")

        for o in order:
            if o < 1:
                raise ValueError("order should consist of positive integers")

        # compute approximant basis
        # if required, normalize it into shifted Popov form
        if row_wise:
            P,rdeg = self._approximant_basis_iterative(order, shifts)
            if normal_form:
                # compute the list "- pivot degree"
                # (since weak Popov, pivot degree is rdeg-shifts entrywise)
                # Note: -deg(P[i,i]) = shifts[i] - rdeg[i]
                degree_shifts = [shifts[i] - rdeg[i] for i in range(m)]
                # compute approximant basis with that list as shifts
                P,rdeg = self._approximant_basis_iterative(order,
                        degree_shifts)
                # left-multiply by inverse of leading matrix
                lmat = P.leading_matrix(shifts=degree_shifts)
                P = lmat.inverse() * P
        else:
            P,rdeg = self.transpose()._approximant_basis_iterative(order,
                    shifts)
            if normal_form:
                # compute the list "- pivot degree"
                # (since weak Popov, pivot degree is rdeg-shifts entrywise)
                degree_shifts = [shifts[i] - rdeg[i] for i in range(n)]
                # compute approximant basis with that list as shifts
                P, rdeg = self.transpose()._approximant_basis_iterative(
                    order, degree_shifts)
                P = P.transpose()
                # right-multiply by inverse of leading matrix
                lmat = P.leading_matrix(shifts=degree_shifts, row_wise=False)
                P = P * lmat.inverse()
            else:
                P = P.transpose()

        return P

    def _approximant_basis_iterative(self, order, shifts):
        r"""
        Return a ``shifts``-ordered weak Popov approximant basis for this
        polynomial matrix at order ``order``
        (see :meth:`minimal_approximant_basis` for definitions).

        The output basis is considered row-wise, that is, its rows are
        left-approximants for the columns of ``self``. It is guaranteed that
        the degree of the output basis is at most the maximum of the entries of
        ``order``, independently of ``shifts``.

        The input dimensions are supposed to be sound: the length of ``order``
        must be the number of columns of ``self``, while the length of
        ``shifts`` must be the number of rows of ``self``.

        INPUT:

        - ``order`` -- list of positive integers

        - ``shifts`` -- list of integers

        OUTPUT:

        - a polynomial matrix (the approximant basis ``P``).

        - a list of integers (the shifts-row degrees of ``P``).

        ALGORITHM:

        This is inspired from the iterative algorithms described in [VBB1992]_
        and [BL1994]_ .

        EXAMPLES::

            sage: pR.<x> = GF(7)[]

        This method supports any number of columns or rows, as well as
        arbitrary shifts and orders::

            sage: order = [4, 1, 2]; shifts = [-3, 4]
            sage: pmat = Matrix(pR, [[5*x^3 + 4*x^2 + 4*x + 6, 5, 4],
            ....:                    [2*x^3 + 2*x^2 + 2*x + 3, 6, 6*x + 3]])
            sage: appbas, rdeg = pmat._approximant_basis_iterative(order,
            ....:                                                  shifts)
            sage: appbas.is_minimal_approximant_basis(pmat, order, shifts)
            True

        The returned list is the shifted row degrees of ``appbas``::

            sage: rdeg == appbas.row_degrees(shifts)
            True

        Approximant bases for the zero matrix are all constant unimodular
        matrices; in fact, this algorithm returns the identity::

            sage: pmat = Matrix(pR, 3, 2)
            sage: appbas,rdeg = pmat._approximant_basis_iterative([2,5],
            ....:                                                 [5,0,-4])
            sage: rdeg == [5,0,-4] and appbas == Matrix.identity(pR, 3)
            True
        """
        # Define parameters and perform some sanity checks
        m = self.nrows()
        n = self.ncols()
        polynomial_ring = self.base_ring()
        X = polynomial_ring.gen()

        # 'rest_order': the orders that remains to be dealt with
        # 'rest_index': indices of orders that remains to be dealt with
        rest_order = list(order)
        rest_index = list(range(n))

        # initialization of the residuals (= input self)
        # and of the approximant basis (= identity matrix)
        from sage.matrix.constructor import identity_matrix
        appbas = identity_matrix(polynomial_ring, m)
        residuals = self.__copy__()

        # throughout the algorithm, 'rdeg' will be the shifts-row degrees of
        # 'appbas'
        # --> initially, 'rdeg' is the shift-row degree of the identity matrix
        rdeg = list(shifts)

        while rest_order:
            # invariant:
            #   * appbas is a shifts-ordered weak Popov approximant basis for
            #   (self,doneorder)
            #   where doneorder = the already processed order, that is, the
            #   tuple order-rest_order (entrywise subtraction)
            #   * rdeg is the shifts-row degree of appbas
            #   * residuals is the submatrix of columns (appbas * self)[:,j]
            #   for all j such that rest_order[j] > 0

            # choice for the next coefficient to be dealt with: first of the
            # largest entries in order (--> process 'self' degree-wise, and
            # left to right)
            # Note: one may also consider the first one in order (--> process
            # 'self' columnwise, from left column to right column, set j=0
            # instead of the below), but it seems to often be (barely) slower
            max_rest_order = max(rest_order)
            for ind, value in enumerate(rest_order):
                if value == max_rest_order:
                    j = ind
                    break
            d = order[rest_index[j]] - rest_order[j]

            # coefficient = the coefficient of degree d of the column j of the
            # residual matrix
            # --> this is very likely nonzero and we want to make it zero, so
            # that this column becomes zero mod X^{d+1}
            coefficient = [residuals[i, j][d] for i in range(m)]

            # Lambda: collect rows [i] with nonzero coefficient[i]
            # pi: index of the first row with smallest shift, among those in
            # Lambda
            Lambda = []
            pi = -1
            for i in range(m):
                if coefficient[i] != 0:
                    Lambda.append(i)
                    if pi < 0 or rdeg[i] < rdeg[pi]:
                        pi = i
            if Lambda: # otherwise, nothing to do
                # update all rows in Lambda--{pi}
                Lambda.remove(pi)
                for row in Lambda:
                    scalar = -coefficient[row]/coefficient[pi]
                    appbas.add_multiple_of_row(row, pi, scalar)
                    residuals.add_multiple_of_row(row, pi, scalar)
                # update row pi
                rdeg[pi] += 1
                appbas.rescale_row(pi, X)
                residuals.rescale_row(pi, X)

            # Decrement rest_order[j], unless there is no more work to do in
            # this column, i.e. if rest_order[j] was 1:
            # in this case remove the column j of
            # residual,rest_order,rest_index
            if rest_order[j] == 1:
                residuals = residuals.delete_columns([j])
                rest_order.pop(j)
                rest_index.pop(j)
            else:
                rest_order[j] -= 1
        return appbas, rdeg

    def is_minimal_kernel_basis(self,
            pmat,
            shifts=None,
            row_wise=True,
            normal_form=False):
        r"""
        Return whether this matrix is a left kernel basis in
        ``shifts``-ordered weak Popov form for the polynomial matrix ``pmat``.

        If ``normal_form`` is ``True``, then the kernel basis must furthermore
        be in ``shifts``-Popov form. An error is raised if the input dimensions
        are not sound.

        INPUT:

        - ``pmat`` -- a polynomial matrix

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        - ``row_wise`` -- boolean (default: ``True``); if ``True``
          then the basis is considered row-wise and operates on the left of
          ``pmat``. Otherwise it is column-wise and operates on the right of
          ``pmat``.

        - ``normal_form`` -- boolean (default: ``False``); if
          ``True`` then checks for a basis in ``shifts``-Popov form

        OUTPUT: boolean

        ALGORITHM:

        Verification that the matrix has full rank and is in shifted weak
        Popov form is done via :meth:`is_weak_popov`; verification that the
        matrix is a left kernel basis is done by checking that the rank is
        correct, that the product is indeed zero, and that the matrix is
        saturated, i.e. it has unimodular column bases (see Lemma 6.10 of
        https://arxiv.org/pdf/1807.01272.pdf for details).

        EXAMPLES::

            sage: pR.<x> = GF(97)[]
            sage: pmat = Matrix(pR, [[1], [x], [x**2]])

            sage: kerbas = Matrix(pR, [[x,-1,0], [0,x,-1]])
            sage: kerbas.is_minimal_kernel_basis(pmat)
            True

        A matrix in Popov form which has the right rank, all rows in the
        kernel, but does not generate the kernel::

            sage: kerbas = Matrix(pR, [[x**2,0,-1], [0,x,-1]])
            sage: kerbas.is_minimal_kernel_basis(pmat)
            False

        Shifts and right kernel bases are supported (with ``row_wise``), and one can test whether the kernel basis is normalized in shifted-Popov form (with ``normal_form``)::

            sage: kerbas = Matrix(pR, [[-x,-x**2], [1,0], [0,1]])
            sage: kerbas.is_minimal_kernel_basis(
            ....:     pmat.transpose(), row_wise=False,
            ....:     normal_form=True, shifts=[0,1,2])
            True
        """
        m = pmat.nrows()
        n = pmat.ncols()

        # set default shifts / check shifts dimension
        if shifts is None:
            shifts = [0] * m if row_wise else [0] * n
        elif row_wise and len(shifts) != m:
            raise ValueError('shifts length should be the row dimension of'
                             ' the input matrix')
        elif (not row_wise) and len(shifts) != n:
            raise ValueError('shifts length should be the column dimension'
                             ' of the input matrix')

        # raise an error if self does not have the right dimension
        if row_wise and self.ncols() != m:
            raise ValueError("column dimension should be the row dimension"
                             " of the input matrix")
        elif (not row_wise) and self.nrows() != n:
            raise ValueError("row dimension should be the column dimension"
                             " of the input matrix")

        # check full rank and shifts-(ordered weak) Popov form
        if normal_form and (not self.is_popov(shifts, row_wise, False, False)):
            return False
        if (not normal_form) and (not self.is_weak_popov(shifts, row_wise, True, False)):
            return False

        # check that self consists of kernel vectors
        if row_wise and self * pmat != 0:
            return False
        if (not row_wise) and pmat * self != 0:
            return False

        # check self.rank() is right (the above weak Popov test ensures self
        # has full row rank if row wise, and full column rank if column wise)
        rk = pmat.rank()
        if row_wise and self.nrows()!=m-rk:
            return False
        if (not row_wise) and self.ncols()!=n-rk:
            return False

        # final check: self is row saturated (assuming row wise),
        # since self has full rank this is equivalent to the fact that its
        # columns generate the identity matrix of size m; in particular the
        # column Popov and column Hermite forms of self are the identity
        if row_wise:
            hnf = self.transpose().hermite_form(include_zero_rows=False)
            # TODO benchmark to see if should be replaced by popov_form
        else:
            hnf = self.hermite_form(include_zero_rows=False)
            # TODO benchmark to see if should be replaced by popov_form
        if hnf != 1:
            return False

        return True

    def minimal_kernel_basis(self,
            shifts=None,
            row_wise=True,
            normal_form=False):
        r"""
        Return a left kernel basis in ``shifts``-ordered weak Popov form for
        this polynomial matrix.

        Assuming we work row-wise, if `F` is an `m \times n` polynomial matrix,
        then a left kernel basis for `F` is a polynomial matrix whose rows form
        a basis of the left kernel of `F`, which is the module of polynomial
        vectors `p` of size `m` such that `p F` is zero.

        If ``normal_form`` is ``True``, then the output basis `P` is
        furthermore in ``shifts``-Popov form. By default, `P` is considered
        row-wise, that is, its rows are left kernel vectors for ``self``; if
        ``row_wise`` is ``False`` then its columns are right kernel vectors for
        ``self``.

        An error is raised if the input dimensions are not sound: if working
        row-wise (resp. column-wise), the length of ``shifts`` must be the
        number of rows (resp. columns) of ``self``.

        INPUT:

        - ``shifts`` -- (default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``

        - ``row_wise`` -- boolean (default: ``True``); if ``True``
          then the output basis considered row-wise and operates on the left
          of ``self``. Otherwise it is column-wise and operates on the right
          of ``self``.

        - ``normal_form`` -- boolean (default: ``False``); if
          ``True`` then the output basis is in ``shifts``-Popov form

        OUTPUT: a polynomial matrix

        ALGORITHM: uses minimal approximant basis computation at a
        sufficiently large order so that the approximant basis contains
        a kernel basis as a submatrix.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: pmat = Matrix([[(x+1)*(x+3)], [(x+1)*(x+3)+1]])
            sage: pmat.minimal_kernel_basis()
            [6*x^2 + 3*x + 3   x^2 + 4*x + 3]

            sage: pmat = Matrix([[(x+1)*(x+3)], [(x+1)*(x+4)]])
            sage: pmat.minimal_kernel_basis()
            [6*x + 3   x + 3]

            sage: pmat.minimal_kernel_basis(row_wise=False)
            []

            sage: pmat = Matrix(pR, [[1, x, x**2]])
            sage: pmat.minimal_kernel_basis(row_wise=False, normal_form=True)
            [x 0]
            [6 x]
            [0 6]

            sage: pmat.minimal_kernel_basis(row_wise=False, normal_form=True,
            ....:                           shifts=[0,1,2])
            [  6*x 6*x^2]
            [    1     0]
            [    0     1]

        Some particular cases (matrix is zero, dimension is zero, column is zero)::

            sage: Matrix(pR, 2, 1).minimal_kernel_basis()
            [1 0]
            [0 1]

            sage: Matrix(pR, 2, 0).minimal_kernel_basis()
            [1 0]
            [0 1]

            sage: Matrix(pR, 0, 2).minimal_kernel_basis()
            []

            sage: Matrix(pR, 3, 2, [[1,0],[1,0],[1,0]]).minimal_kernel_basis()
            [6 1 0]
            [6 0 1]

            sage: Matrix(pR, 3, 2, [[x,0],[1,0],[x+1,0]]).minimal_kernel_basis()
            [6 x 0]
            [6 6 1]

        TESTS:

        We check that PR #37208 is fixed::

            sage: Matrix(pR, 2, 0).minimal_kernel_basis().is_sparse()
            False
        """
        from sage.matrix.constructor import matrix

        m = self.nrows()
        n = self.ncols()

        # set default shifts / check shifts dimension
        if shifts is None:
            shifts = [0] * m if row_wise else [0] * n
        elif row_wise and len(shifts) != m:
            raise ValueError('shifts length should be the row dimension')
        elif (not row_wise) and len(shifts) != n:
            raise ValueError('shifts length should be the column dimension')

        # compute kernel basis
        if row_wise:
            if m <= n and self.constant_matrix().rank() == m:
                # early exit: kernel is empty; note: this covers the case m==0
                return matrix(self.base_ring(), 0, m)

            if n == 0: # early exit: kernel is identity
                return matrix.identity(self.base_ring(), m)

            d = self.degree() # well defined since m > 0 and n > 0
            if d == -1: # matrix is zero: kernel is identity
                return matrix.identity(self.base_ring(), m)

            # degree bounds on the kernel basis
            degree_bound = min(m,n)*d+max(shifts)
            degree_bounds = [degree_bound - shifts[i] for i in range(m)]

            # orders for approximation
            orders = self.column_degrees(degree_bounds)
            for i in range(n): orders[i] = orders[i]+1

            # note: minimal_approximant_basis requires orders[i] > 0
            # -> if d>0, then degree_bounds > 0 entry-wise and this tuple
            # `orders` already has all entries strictly positive
            # -> if d==0, then `orders[i]` is zero exactly when the column i
            # of self is zero; we may as well take orders[i] == 1 for such
            # columns which do not influence the left kernel
            if d == 0:
                for i in range(n):
                    if orders[i] == 0:
                        orders[i] = 1

            # compute approximant basis and retrieve kernel rows
            P = self.minimal_approximant_basis(orders,shifts,True,normal_form)
            row_indices = []
            for i in range(m):
                if P[i, i].degree() + shifts[i] <= degree_bound:
                    row_indices.append(i)
            return P[row_indices,:]

        else:
            if n <= m and self.constant_matrix().rank() == n:
                # early exit: kernel is empty; this covers the case n==0
                return matrix(self.base_ring(), n, 0)

            if m == 0: # early exit: kernel is identity
                return matrix.identity(self.base_ring(), n)

            d = self.degree() # well defined since m > 0 and n > 0
            if d == -1: # matrix is zero
                return matrix.identity(self.base_ring(), n)

            # degree bounds on the kernel basis
            degree_bound = min(m,n)*d+max(shifts)
            degree_bounds = [degree_bound - shifts[i] for i in range(n)]

            # orders for approximation
            orders = self.row_degrees(degree_bounds)
            for i in range(m): orders[i] = orders[i]+1

            # note: minimal_approximant_basis requires orders[i] > 0
            # -> if d>0, then degree_bounds > 0 entry-wise and this tuple
            # `orders` already has all entries strictly positive
            # -> if d==0, then `orders[i]` is zero exactly when the row i
            # of self is zero; we may as well take orders[i] == 1 for such
            # rows which do not influence the right kernel
            if d == 0:
                for i in range(m):
                    if orders[i] == 0:
                        orders[i] = 1

            # compute approximant basis and retrieve kernel columns
            P = self.minimal_approximant_basis(orders,shifts,False,normal_form)
            column_indices = []
            for j in range(n):
                if P[j, j].degree() + shifts[j] <= degree_bound:
                    column_indices.append(j)
            return P[:,column_indices]

    def _basis_completion_via_reversed_approx(self):
        r"""
        Return a Smith form-preserving nonsingular completion of a row basis of
        this matrix. For a more detailed description, see
        :meth:`basis_completion`, in the row-wise case.

        EXAMPLES:

        Three polynomials whose GCD is `1` can be completed into a unimodular
        matrix::

            sage: ring.<x> = GF(7)[]
            sage: mat = matrix([[x*(x-1)*(x-2), (x-2)*(x-3)*(x-4), (x-4)*(x-5)*(x-6)]])
            sage: mat
            [      x^3 + 4*x^2 + 2*x   x^3 + 5*x^2 + 5*x + 4   x^3 + 6*x^2 + 4*x + 6]
            sage: rcomp = mat._basis_completion_via_reversed_approx(); rcomp
            [        5*x^2 + 4*x + 1             5*x^2 + 2*x                   5*x^2]
            [          2*x^3 + 4*x^2 2*x^3 + 6*x^2 + 2*x + 1       2*x^3 + x^2 + 3*x]
            sage: basis = mat.stack(rcomp); basis
            [      x^3 + 4*x^2 + 2*x   x^3 + 5*x^2 + 5*x + 4   x^3 + 6*x^2 + 4*x + 6]
            [        5*x^2 + 4*x + 1             5*x^2 + 2*x                   5*x^2]
            [          2*x^3 + 4*x^2 2*x^3 + 6*x^2 + 2*x + 1       2*x^3 + x^2 + 3*x]
            sage: basis.determinant()
            6

        The following matrix has rank `2` and trivial Smith form. It can be
        completed row-wise into a `3 \times 3` unimodular matrix::

            sage: mat = matrix(ring, 2, 3, \
                    [[x^2 + 5*x + 5,   3*x^2 + x + 3, 4*x^2 + 5*x + 4], \
                     [5*x^2 + 4*x,   3*x^2 + 4*x + 5, 5*x^2 + 5*x + 3]])
            sage: rcomp = mat._basis_completion_via_reversed_approx(); rcomp
            [  2*x^2 + 1 4*x^2 + 3*x 2*x^2 + 3*x]
            sage: mat.stack(rcomp).determinant()
            3

        The following matrix has rank 1 and its nonzero Smith factor is `x+3`.
        A row-wise completion has a single nonzero row, whereas a column-wise
        completion has two columns; in both cases, the Smith form is preserved::

            sage: mat = matrix(ring, 3, 2, \
                    [[    x^3 + x^2 + 5*x + 5,         2*x^3 + 2*x + 4], \
                     [  3*x^3 + 2*x^2 + x + 3,   6*x^3 + 5*x^2 + x + 1], \
                     [2*x^3 + 5*x^2 + 3*x + 4, 4*x^3 + 6*x^2 + 5*x + 6]])
            sage: mat.smith_form(transformation=False)
            [x + 3     0]
            [    0     0]
            [    0     0]
            sage: rcomp = mat._basis_completion_via_reversed_approx(); rcomp
            [x + 1   2*x]
            sage: ccomp = mat.transpose()._basis_completion_via_reversed_approx().transpose(); ccomp
            [3*x + 1 4*x + 4]
            [    2*x 5*x + 1]
            [    6*x       x]
            sage: rcomp.stack(mat).smith_form(transformation=False)
            [    1     0]
            [    0 x + 3]
            [    0     0]
            [    0     0]
            sage: ccomp.augment(mat).smith_form(transformation=False)
            [    1     0     0     0]
            [    0     1     0     0]
            [    0     0 x + 3     0]

        TESTS:

        Corner cases are handled correctly::

            sage: matrix(ring, 0, 0)._basis_completion_via_reversed_approx()
            []
            sage: matrix(ring, 0, 2)._basis_completion_via_reversed_approx()
            [1 0]
            [0 1]
            sage: matrix(ring, 2, 0)._basis_completion_via_reversed_approx()
            []
        """
        from sage.matrix.constructor import matrix  # for identity

        ring = self.base_ring()
        m = self.nrows()
        n = self.ncols()

        # corner cases: after this, m>0 and n>0
        if m == 0:
            return matrix.identity(ring, n)
        if n == 0:
            return matrix(ring, 0, 0)

        # find column degrees (zero columns have degree -1)
        cdeg = self.column_degrees()

        # list zero and nonzero columns
        zcols = []
        nonzcols = []
        for j in range(n):
            if cdeg[j] < 0:
                zcols.append(j)
            else:
                nonzcols.append(j)

        if len(nonzcols) == 0:
            return matrix.identity(ring, n)

        # restrict to nonzero columns, and reverse entries column-wise
        mat = self.matrix_from_columns(nonzcols)
        cdeg = [cdeg[j] for j in nonzcols]  # cdeg >= 0 entrywise
        mat_rev = mat.reverse(row_wise=False, degree=cdeg)

        # compute shifted-minimal right kernel basis
        kernel_basis = mat_rev.minimal_kernel_basis(row_wise=False, shifts=cdeg)
        # if kernel_basis has zero columns, then mat is full column rank, there is
        # nothing to complete, just return the trivial completion
        if kernel_basis.ncols() == 0:
            return matrix.identity(ring, n).matrix_from_rows(zcols)

        # compute shifted-minimal left approximant basis
        # with shifts -cdeg and approximation orders cdeg(kernel_basis)+1
        mcdeg = [-c for c in cdeg]
        orders = [d+1 for d in kernel_basis.column_degrees(shifts=cdeg)]
        approx_basis = kernel_basis.minimal_approximant_basis(order=orders, shifts=mcdeg)

        # idea: by choice of parameters (shifts and orders), this approximant
        # basis has a subset of rows which form a saturated basis of the
        # rational vector space generated by mat_rev; this is precisely all
        # rows `row` of approx_basis that satisfy row * kernel_basis == 0
        # --> we can actually detect them without computing the product,
        # these rows are those which have "small" -cdeg-shifted row degree

        # select the appropriate rows from approx basis
        rdeg = approx_basis.row_degrees(shifts=mcdeg)
        completion_indices = [i for i in range(len(rdeg)) if rdeg[i] > 0]
        completion_rev = approx_basis.matrix_from_rows(completion_indices)
        rdeg = [rdeg[i] for i in completion_indices]

        # now we know the completion len(completion_indices) rows,
        # plus the trivial rows related to zero columns we have discarded
        completion = matrix(ring, len(zcols)+len(completion_indices), n)

        # fill trivial part of result matrix:
        # to each zero column of input mat corresponds an identity row
        # -> the following loop has the same effect as the next commented line:
        #   completion[0:len(zcols),:] = matrix.identity(ring, n).matrix_from_rows(zcols)
        for k in range(len(zcols)):
            completion[k,zcols[k]] = 1

        # fill nontrivial part of result matrix: this is
        # diag(x**rdeg) * completion_rev(1/x) * diag(x**cdeg)
        #
        # warning: cannot use reverse w.r.t rdeg or cdeg here,
        # deg(completion_rev) may a priori be greater than max(cdeg) and
        # max(rdeg), so there would be negative degree coefficients that would
        # be discarded
        # -> the following code has the same effect as the next 3 commented lines:
        # left_diag = matrix.diagonal([ring.monomial(d) for d in rdeg])
        # right_diag = matrix.diagonal([ring.monomial(d) for d in cdeg])
        # completion[len(zcols):,nonzcols] = left_diag * completion_rev(var**(-1)) * right_diag
        # where var = ring.gen()
        for i in range(len(completion_indices)):
            for j in range(len(nonzcols)):
                completion[i + len(zcols), nonzcols[j]] = completion_rev[i,j].reverse(rdeg[i]+cdeg[j])

        return completion

    def basis_completion(self, row_wise=True, algorithm='approximant'):
        r"""
        Return a Smith form-preserving nonsingular completion of a basis of
        this matrix: row-wise completing a row basis if ``row_wise`` is True;
        column-wise completing a column basis if it is False.

        For a more detailed description, consider the row-wise case (the
        column-wise case is the same up to matrix transposition). Let `A` be
        the input matrix, `m \times n` over univariate polynomials
        `\Bold{K}[x]`, for some field `\Bold{K}`, and let `r` be the rank of
        `A`, which is unknown a priori. This computes a matrix `C` of
        dimensions `(n-r) \times n` such that stacking both matrices one above
        the other, say `[[A],[C]]`, gives a matrix of maximal rank `n` and with
        the same nontrivial Smith factors as `A`. In particular, `C` has full
        row rank, and the rank of the input matrix may be recovered from the
        number of rows of `C`.

        As a consequence, if `B` is a basis of the module generated by the rows
        of `A` (for example `B = A` if `A` has full row rank), then `[[B],[C]]`
        is nonsingular, and its determinant is the product of the nonzero Smith
        factors of `A` up to multiplication by a nonzero element of `\Bold{K}`.

        In particular, for `A` with full row rank: if the rows `A` can be
        completed into a basis of `\Bold{K}[x]^{n}` (or equivalently, `A` has
        unimodular column bases, or also, if the rows of `A` generate all
        polynomial vectors in the rational row space of `A`), then `C` provides
        such a completion. In this case, `[[A],[C]]` is unimodular: it is
        invertible over `\Bold{K}[x]`, and `det([[A],[C]])` is a nonzero
        element of the base field `\Bold{K}`.

        INPUT:

        - ``row_wise`` -- boolean (default: ``True``); if ``True`` then
          compute a row-wise completion, else compute a column-wise completion

        - ``algorithm`` -- (default: ``'approximant'``) selects the
          approach for computing the completion; currently supported:
          ``'approximant'`` and ``'smith'``

        OUTPUT: a matrix over the same base ring as the input matrix, which
        forms a completion as defined above

        ALGORITHM:

        - ``'approximant'`` -- the approximant-based algorithm follows the ideas in
          [ZL2014]_ , based on polynomial reversals combined with the
          computation of a minimal kernel basis and a minimal approximant
          basis.

        - ``'smith'`` -- the Smith form-based algorithm computes the Smith form of
          this matrix along with corresponding unimodular transformations, from
          which a completion is readily obtained.

        EXAMPLES:

        Three polynomials whose GCD is `1` can be completed into a unimodular
        matrix::

            sage: ring.<x> = GF(7)[]
            sage: mat = matrix([[x*(x-1)*(x-2), (x-2)*(x-3)*(x-4), (x-4)*(x-5)*(x-6)]])
            sage: mat
            [      x^3 + 4*x^2 + 2*x   x^3 + 5*x^2 + 5*x + 4   x^3 + 6*x^2 + 4*x + 6]
            sage: rcomp = mat.basis_completion(); rcomp
            [        5*x^2 + 4*x + 1             5*x^2 + 2*x                   5*x^2]
            [          2*x^3 + 4*x^2 2*x^3 + 6*x^2 + 2*x + 1       2*x^3 + x^2 + 3*x]
            sage: basis = mat.stack(rcomp); basis
            [      x^3 + 4*x^2 + 2*x   x^3 + 5*x^2 + 5*x + 4   x^3 + 6*x^2 + 4*x + 6]
            [        5*x^2 + 4*x + 1             5*x^2 + 2*x                   5*x^2]
            [          2*x^3 + 4*x^2 2*x^3 + 6*x^2 + 2*x + 1       2*x^3 + x^2 + 3*x]
            sage: basis.determinant()
            6

        The following matrix has rank `2` and trivial Smith form. It can be
        completed row-wise into a `3 \times 3` unimodular matrix (column-wise,
        there is nothing to complete)::

            sage: mat = matrix(ring, 2, 3, \
                    [[x^2 + 5*x + 5,   3*x^2 + x + 3, 4*x^2 + 5*x + 4], \
                     [5*x^2 + 4*x,   3*x^2 + 4*x + 5, 5*x^2 + 5*x + 3]])
            sage: rcomp = mat.basis_completion(); rcomp
            [  2*x^2 + 1 4*x^2 + 3*x 2*x^2 + 3*x]
            sage: mat.stack(rcomp).determinant()
            3
            sage: mat.basis_completion(row_wise=False)
            []

        The following matrix has rank 1 and its nonzero Smith factor is `x+3`.
        A row-wise completion has a single nonzero row, whereas a column-wise
        completion has two columns; in both cases, the Smith form is preserved::

            sage: mat = matrix(ring, 3, 2, \
                    [[    x^3 + x^2 + 5*x + 5,         2*x^3 + 2*x + 4], \
                     [  3*x^3 + 2*x^2 + x + 3,   6*x^3 + 5*x^2 + x + 1], \
                     [2*x^3 + 5*x^2 + 3*x + 4, 4*x^3 + 6*x^2 + 5*x + 6]])
            sage: mat.smith_form(transformation=False)
            [x + 3     0]
            [    0     0]
            [    0     0]
            sage: rcomp = mat.basis_completion(); rcomp
            [x + 1   2*x]
            sage: ccomp = mat.basis_completion(row_wise=False); ccomp
            [3*x + 1 4*x + 4]
            [    2*x 5*x + 1]
            [    6*x       x]
            sage: rcomp.stack(mat).smith_form(transformation=False)
            [    1     0]
            [    0 x + 3]
            [    0     0]
            [    0     0]
            sage: ccomp.augment(mat).smith_form(transformation=False)
            [    1     0     0     0]
            [    0     1     0     0]
            [    0     0 x + 3     0]

        Here are a few more examples, similar to the above but over fields
        other than ``GF(7)``::

            sage: ring.<x> = QQ[]
            sage: mat = matrix([[x*(x-1)*(x-2), (x-2)*(x-3)*(x-4), (x-4)*(x-5)*(x-6)]])
            sage: mat
            [        x^3 - 3*x^2 + 2*x   x^3 - 9*x^2 + 26*x - 24 x^3 - 15*x^2 + 74*x - 120]
            sage: rcomp = mat.basis_completion(algorithm='smith'); rcomp
            [        -1/12*x - 1/12         -1/12*x + 5/12                      0]
            [                  1/12                   1/12 1/24*x^2 - 13/24*x + 2]
            sage: mat.stack(rcomp).determinant()
            1

            sage: mat = matrix([[x*(x-1), x*(x-2)], \
                                [x*(x-2), x*(x-3)], \
                                [(x-1)*(x-2), (x-1)*(x-3)]])
            sage: mat.smith_form(transformation=False)
            [1 0]
            [0 x]
            [0 0]
            sage: ccomp = mat.basis_completion(row_wise=False, algorithm='smith')
            sage: ccomp
            [1/2*x - 1/2]
            [  1/2*x - 1]
            [1/2*x - 3/2]
            sage: ccomp.augment(mat).smith_form(transformation=False)
            [    1     0     0]
            [    0     1     0]
            [    0     0 1/2*x]

            sage: field.<a> = NumberField(x**2 - 2)
            sage: ring.<y> = field[]
            sage: mat = matrix([[3*a*y - 1, (-8*a - 1)*y - 2*a + 1]])
            sage: rcomp = mat.basis_completion(algorithm='smith'); rcomp
            [ 39/119*a - 30/119 -99/119*a + 67/119]
            sage: mat.stack(rcomp).determinant()
            1

        TESTS:

        Corner cases are handled correctly::

            sage: matrix(ring, 0, 0).basis_completion()
            []
            sage: matrix(ring, 0, 0).basis_completion(row_wise=False)
            []
            sage: matrix(ring, 0, 0).basis_completion(algorithm='smith')
            []
            sage: matrix(ring, 0, 0).basis_completion(row_wise=False,algorithm='smith')
            []
            sage: matrix(ring, 0, 2).basis_completion()
            [1 0]
            [0 1]
            sage: matrix(ring, 0, 2).basis_completion(row_wise=False)
            []
            sage: matrix(ring, 0, 2).basis_completion(algorithm='smith')
            [1 0]
            [0 1]
            sage: matrix(ring, 0, 2).basis_completion(row_wise=False,algorithm='smith')
            []
            sage: matrix(ring, 2, 0).basis_completion()
            []
            sage: matrix(ring, 2, 0).basis_completion(row_wise=False)
            [1 0]
            [0 1]
            sage: matrix(ring, 2, 0).basis_completion(algorithm='smith')
            []
            sage: matrix(ring, 2, 0).basis_completion(row_wise=False,algorithm='smith')
            [1 0]
            [0 1]
        """
        if algorithm == "approximant":
            if row_wise:
                return self._basis_completion_via_reversed_approx()
            else:
                Ctrsp = self.transpose()._basis_completion_via_reversed_approx()
                return Ctrsp.transpose()

        elif algorithm == "smith":
            ring = self.base_ring()
            m = self.nrows()
            n = self.ncols()

            # Smith form, S == U * self * V
            S, U, V = self.smith_form(transformation=True)

            # S is m x n with rk = rank(self) first diagonal entries nonzero, find rk
            rk = 0
            while rk < min(m, n) and not S[rk,rk].is_zero():
                rk += 1

            # case: matrix is zero (including if m == 0 || n == 0)
            from sage.matrix.constructor import matrix
            if rk == 0:
                if row_wise:
                    return matrix.identity(ring, n)
                else:
                    return matrix.identity(ring, m)

            # now, matrix is nonzero (and nonempty)
            if row_wise:
                C = matrix.identity(ring, n)[rk:,:]
                VV = V.inverse_of_unit()
                return C * VV
            else:
                C = matrix.identity(ring, m)[:,rk:]
                UU = U.inverse_of_unit()
                return UU * C

        else:
            raise ValueError("algorithm must be one of \"approximant\" or \"smith\".")

    def _is_basis_completion(self, mat, row_wise=True):
        r"""
        Return whether this matrix is a basis completion for ``mat``. For the
        definition of basis completion, including its orientation row-wise or
        column-wise, see the documentation of :meth:`basis_completion`.

        INPUT:

        - ``row_wise`` -- boolean (default: ``True``); if ``True`` then
          check for row-wise completion, else check for column-wise completion

        OUTPUT: boolean indicating whether this matrix is a completion of ``mat``

        EXAMPLES:

        Using the same examples as :meth:`basis_completion`::

            sage: ring.<x> = GF(7)[]
            sage: mat1 = matrix([[x*(x-1)*(x-2), (x-2)*(x-3)*(x-4), (x-4)*(x-5)*(x-6)]])
            sage: rcomp1 = matrix(ring, 2, 3, \
                    [[5*x^2 + 4*x + 1, 5*x^2 + 2*x, 5*x^2], \
                     [2*x^3 + 4*x^2, 2*x^3 + 6*x^2 + 2*x + 1, 2*x^3 + x^2 + 3*x]])
            sage: rcomp1._is_basis_completion(mat1)
            True

            sage: mat2 = matrix(ring, 2, 3, \
                    [[x^2 + 5*x + 5,   3*x^2 + x + 3, 4*x^2 + 5*x + 4], \
                     [5*x^2 + 4*x,   3*x^2 + 4*x + 5, 5*x^2 + 5*x + 3]])
            sage: rcomp2 = matrix(ring, 1, 3, [[2*x^2 + 1, 4*x^2 + 3*x, 2*x^2 + 3*x]])
            sage: rcomp2._is_basis_completion(mat2)
            True
            sage: ccomp2 = matrix(ring, 2, 0)
            sage: ccomp2._is_basis_completion(mat2, row_wise=False)
            True

            sage: mat3 = matrix(ring, 3, 2, \
                    [[    x^3 + x^2 + 5*x + 5,         2*x^3 + 2*x + 4], \
                     [  3*x^3 + 2*x^2 + x + 3,   6*x^3 + 5*x^2 + x + 1], \
                     [2*x^3 + 5*x^2 + 3*x + 4, 4*x^3 + 6*x^2 + 5*x + 6]])
            sage: rcomp3 = matrix(ring, 1, 2, [[x + 1, 2*x]])
            sage: rcomp3._is_basis_completion(mat3)
            True
            sage: ccomp3 = matrix(ring, 3, 2, \
                                    [[3*x + 1, 4*x + 4], \
                                     [    2*x, 5*x + 1], \
                                     [    6*x,       x]])
            sage: ccomp3._is_basis_completion(mat3, row_wise=False)
            True

        A row-wise completion is generally not a column-wise completion (most
        often, matrix dimensions are not even compatible), one exception being
        the completions of square zero matrices::

            sage: rcomp2._is_basis_completion(mat2, row_wise=False)
            Traceback (most recent call last):
            ...
            TypeError: number of rows must be the same, 2 != 1
            sage: ccomp2._is_basis_completion(mat2, row_wise=True)
            Traceback (most recent call last):
            ...
            TypeError: number of columns must be the same, not 3 and 0
            sage: rcomp3._is_basis_completion(mat3, row_wise=False)
            Traceback (most recent call last):
            ...
            TypeError: number of rows must be the same, 3 != 1
            sage: ccomp3._is_basis_completion(mat3, row_wise=True)
            False

            sage: zero_mat = matrix(ring, 2, 2)
            sage: comp = zero_mat.basis_completion(); comp
            [1 0]
            [0 1]
            sage: comp._is_basis_completion(zero_mat, row_wise=True)
            True
            sage: comp._is_basis_completion(zero_mat, row_wise=False)
            True

        Completions that do not preserve the Smith factors are not valid,
        even when the sought rank is reached::

            sage: (x * rcomp2)._is_basis_completion(mat2)
            False
            sage: ((x+2) * rcomp3)._is_basis_completion(mat3)
            False
            sage: (ccomp3 * matrix.diagonal([x,1]))._is_basis_completion(mat3, row_wise=False)
            False

        Preserving Smith factors without reaching full rank is not a valid
        completion::

            sage: mat = matrix(ring, [[1,0,0]])
            sage: matrix(ring, [[0,1,0]])._is_basis_completion(mat)
            False
            sage: matrix(ring, [[0,1,0], [0,0,1]])._is_basis_completion(mat)
            True

        TESTS:

        Corner cases are handled correctly::

            sage: empty_mat = matrix(ring, 0, 0)
            sage: empty_rows = matrix(ring, 2, 0)
            sage: empty_columns = matrix(ring, 0, 2)
            sage: id22 = matrix.identity(ring, 2)

            sage: empty_mat._is_basis_completion(empty_mat)
            True
            sage: empty_mat._is_basis_completion(empty_mat, row_wise=False)
            True
            sage: empty_mat._is_basis_completion(empty_rows)
            True
            sage: empty_mat._is_basis_completion(empty_columns, row_wise=False)
            True
            sage: empty_mat._is_basis_completion(empty_columns)
            Traceback (most recent call last):
            ...
            TypeError: number of columns must be the same, not 2 and 0

            sage: empty_columns._is_basis_completion(id22)
            True
            sage: empty_columns._is_basis_completion(id22, row_wise=False)
            Traceback (most recent call last):
            ...
            TypeError: number of rows must be the same, 2 != 0
            sage: empty_columns._is_basis_completion(empty_columns)
            False
            sage: empty_rows._is_basis_completion(id22, row_wise=False)
            True
            sage: empty_rows._is_basis_completion(empty_rows, row_wise=True)
            False
        """

        # compute completed matrix
        if row_wise:
            cmat = mat.stack(self)
        else:
            cmat = mat.augment(self)

        # compute (nonzero) invariant factors for mat, from largest to nonzero smallest
        snf_mat = mat.smith_form(transformation=False)
        mat_factors = [snf_mat[i,i] for i in range(min(mat.nrows(),mat.ncols())-1, -1, -1)
                                            if not snf_mat[i,i].is_zero()]
        rk = len(mat_factors)

        # check completion has right completing dimension
        if row_wise and self.nrows() != (mat.ncols() - rk):
            return False
        if (not row_wise) and self.ncols() != (mat.nrows() - rk):
            return False

        # compute invariant factors for completed matrix, from largest to smallest
        snf_cmat = cmat.smith_form(transformation=False)
        cmat_factors = [snf_cmat[i,i] for i in range(min(cmat.nrows(),cmat.ncols())-1, -1, -1)]

        # check first rk = rank(mat) factors match
        # note: largest factors may not be monic, compare their monic counterparts
        if rk > 0 and mat_factors[0].monic() != cmat_factors[0].monic():
            return False
        # other factors are monic
        if mat_factors[1:] != cmat_factors[1:rk]:
            return False

        # check remaining factors for cmat are 1
        # this guarantees Smith-form preserving property, and also
        # the fact that cmat has the required rank
        # (this rank being mat.ncols() if row-wise; mat.nrows() if column-wise)
        if not all(f.is_one() for f in cmat_factors[rk:]):
            return False

        return True
