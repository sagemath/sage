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
"""
# ****************************************************************************
#       Copyright (C) 2016 Kwankyu Lee <ekwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.matrix.matrix2 cimport Matrix
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ


cdef class Matrix_polynomial_dense(Matrix_generic_dense):
    r"""
    Dense matrix over a univariate polynomial ring over a field.

    For a field $\Bold{K}$, we consider matrices over the univariate
    polynomial ring $\Bold{K}[x]$.
    
    They are often used to represent bases of some $\Bold{K}[x]$-modules. In
    this context, there are two possible representations which are both
    commonly used in the literature.

    - Working column-wise: each column of the matrix is a vector in the basis;
      then, a $\\Bold{K}[x]$-submodule of $\\Bold{K}[x]^{m}$ of rank $n$ is
      represented by an $m \\times n$ matrix, whose columns span the module
      (via $\\Bold{K}[x]$-linear combinations). This matrix has full rank,
      and $n \\leq m$.

    - Working row-wise: each row of the matrix is a vector in the basis; then,
      a $\\Bold{K}[x]$-submodule of $\\Bold{K}[x]^{n}$ of rank $m$ is
      represented by an $m \\times n$ matrix, whose rows span the module (via
      $\\Bold{K}[x]$-linear combinations). This matrix has full rank, and $m
      \\leq n$.

    For the rest of this class description, we assume that one is working
    row-wise. For a given such module, all its bases are equivalent under
    left-multiplication by a unimodular matrix, that is, a matrix which has
    determinant in $\Bold{K}\setminus\{0\}$.

    There are bases which are called reduced or minimal: their rows have the
    minimal degree possible among all bases of this module. The degree of a row
    is the maximum of the degrees of the entries of the row. An equivalent
    condition is that the leading matrix of this basis has full rank (see the
    description of :meth:`leading_matrix`). There is a unique minimal basis,
    called the Popov basis of the module, which satisfies some additional
    normalization condition (see the description of :meth:`row_degrees`).

    These notions can be extended via a more general degree measure, involving
    a tuple of integers which is called shift and acts as column degree shifts
    in the definition of row degree. Precisely, for given $s_1,\ldots,s_n \in
    \ZZ$ and a row vector $[p_1 \; \cdots \; p_n] \in \Bold{K}[x]^{1 \times
    n}$, its shifted row degree is the maximum of $\deg(p_j) + s_j$ for $1 \leq
    j \leq n$. Then, minimal bases and Popov bases are defined similarly, with
    respect to this notion of degree.

    Another important canonical basis is the Hermite basis, which is a lower
    triangular matrix satisfying a normalization condition similar to that for
    the Popov basis. In fact, if $d$ is the largest degree appearing in the
    Hermite basis, then the Hermite basis coincide with the shifted Popov basis
    with the shifts $(0,d,2d,\ldots,(n-1)d)$.
    """

    def _check_shift_dimension(self, shifts, row_wise=True):
        r"""
        Raises an exception if the ``shifts`` argument does not have the right
        length.

        For an $m \times n$ polynomial matrix, if working row-wise then
        ``shifts`` should have $n$ entries; if working column-wise, it should
        have $m$ entries.

        INPUT:

        - ``shifts`` -- list of integers, or ``None``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, if ``True``
          then shifts apply to the columns of the matrix and otherwise to its
          rows (see the class description for more details).

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix( pR, [[3*x+1, 0, 1], [x^3+3, 0, 0]])
            sage: M._check_shift_dimension(shifts=[1,3,2])

            sage: M._check_shift_dimension(shifts=[1,3,2], row_wise=False)
            Traceback (most recent call last):
            ...
            ValueError: shifts length should be the row dimension
        """
        if shifts != None and (not row_wise) and len(shifts) != self.nrows():
            raise ValueError('shifts length should be the row dimension')
        if shifts != None and (row_wise and len(shifts) != self.ncols()):
            raise ValueError('shifts length should be the column dimension')

    def degree(self):
        r"""
        Return the degree of this matrix.

        For a given polynomial matrix, its degree is the maximum of the degrees
        of all its entries. If the matrix is nonzero, this is a nonnegative
        integer; here, the degree of the zero matrix is -1.

        OUTPUT: an integer.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix( pR, [[3*x+1, 0, 1], [x^3+3, 0, 0]])
            sage: M.degree()
            3

        The zero matrix has degree ``-1``::

            sage: M = Matrix( pR, 2, 3 )
            sage: M.degree()
            -1

        For an empty matrix, the degree is not defined::

            sage: M = Matrix( pR, 3, 0 )
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

        For a given polynomial matrix $M = (M_{i,j})_{i,j}$, its degree matrix
        is the matrix $(\deg(M_{i,j}))_{i,j}$ formed by the degrees of its
        entries. Here, the degree of the zero polynomial is $-1$.

        For given shifts $s_1,\ldots,s_m \in \ZZ$, the shifted degree
        matrix of $M$ is either $(\deg(M_{i,j})+s_j)_{i,j}$ if working
        row-wise, or $(\deg(M_{i,j})+s_i)_{i,j}$ if working column-wise. In the
        former case, $m$ has to be the number of columns of $M$; in the latter
        case, the number of its rows. Here, if $M_{i,j}=0$ then the
        corresponding entry in the shifted degree matrix is
        $\min(s_1,\ldots,s_m)-1$. For more on shifts and working row-wise
        versus column-wise, see the class documentation.

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, if ``True``
          then shifts apply to the columns of the matrix and otherwise to its
          rows (see the class description for more details).

        OUTPUT: an integer matrix.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix( pR, [[3*x+1, 0, 1], [x^3+3, 0, 0]])
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
        from sage.matrix.constructor import Matrix
        zero_degree = min(shifts) - 1
        if row_wise: 
            return Matrix( ZZ, [[ self[i,j].degree() + shifts[j]
                if self[i,j] != 0 else zero_degree
                for j in range(self.ncols()) ] for i in range(self.nrows())] )
        else:
            return Matrix( ZZ, [[ self[i,j].degree() + shifts[i]
                if self[i,j] != 0 else zero_degree
                for j in range(self.ncols()) ] for i in range(self.nrows())] )

    def constant_matrix(self):
        r"""
        Return the constant coefficient of this matrix seen as a polynomial
        with matrix coefficients; this is also this matrix evaluated at zero.

        OUTPUT: a matrix over the base field.

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
        from sage.matrix.constructor import Matrix
        return Matrix([[self[i,j].constant_coefficient()
            for j in range(self.ncols())] for i in range(self.nrows())])

    def is_constant(self):
        r"""
        Return ``True`` if and only if this polynomial matrix is constant,
        that is, all its entries are constant.

        OUTPUT: a boolean.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]

            sage: M = Matrix([
            ....:    [  x^3+5*x^2+5*x+1,       5,       6*x+4,         0],
            ....:    [      6*x^2+3*x+1,       1,           2,         0],
            ....:    [2*x^3+4*x^2+6*x+4, 5*x + 1, 2*x^2+5*x+5, x^2+5*x+6]
            ....:     ])
            sage: M.is_constant()
            False
            sage: M = Matrix(pR,[[1,5,2],[3,1,5]]); M.is_constant()
            True
            sage: M = Matrix.zero(pR,3,5); M.is_constant()
            True

        .. SEEALSO::

            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.is_constant` .
        """
        return all([self[i,j].is_constant()
            for j in range(self.ncols()) for i in range(self.nrows())])

    def coefficient_matrix(self,d,row_wise=True):
        r"""
        Return the constant matrix which is obtained from this matrix by taking
        the coefficient of its entries with degree specified by `d`.

        - if `d` is an integer, this selects the coefficient of `d` for all
          entries;
        - if `d` is a list $(d_1,\ldots,d_m)$ and ``row_wise`` is ``True``,
          this selects the coefficient of degree $d_i$ for all entries of the
          $i$th row for each $i$;
        - if `d` is a list $(d_1,\ldots,d_n)$ and ``row_wise`` is ``False``,
          this selects the coefficient of degree $d_i$ for all entries of the
          $j$th column for each $j$.

        INPUT:

        - ``d`` -- a list of integers, or an integer,

        - ``row_wise`` -- (optional, default: ``True``) boolean, if ``True``
          (resp. ``False``) then `d` should be a list of length equal to the
          row (resp. column) dimension of this matrix.

        OUTPUT: a matrix over the base field.

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
            raise ValueError("length of input degree list should be the " \
                                      + "row dimension of the input matrix")
        elif (not row_wise) and len(d) != n:
            raise ValueError("length of input degree list should be the " \
                                      + "column dimension of the input matrix")

        from sage.matrix.constructor import Matrix
        return Matrix(self.base_ring().base_ring(), m, n,
                [[self[i,j][d[i]] if row_wise else self[i,j][d[j]]
            for j in range(n)] for i in range(m)])

    def truncate(self, d, row_wise=True):
        r"""
        Return the matrix which is obtained from this matrix after truncating
        all its entries according to precisions specified by `d`.

        - if `d` is an integer, the truncation is at precision `d` for all
          entries;
        - if `d` is a list $(d_1,\ldots,d_m)$ and ``row_wise`` is ``True``, all
          entries of the $i$th row are truncated at precision $d_i$ for each
          $i$;
        - if `d` is a list $(d_1,\ldots,d_n)$ and ``row_wise`` is ``False``,
          all entries of the $j$th column are truncated at precision $d_j$ for
          each $j$.

        Here the convention for univariate polynomials is to take zero
        for the truncation for a negative `d`.

        INPUT:

        - ``d`` -- a list of integers, or an integer,

        - ``row_wise`` -- (optional, default: ``True``) boolean, if ``True``
          (resp. ``False``) then `d` should be a list of length equal to the
          row (resp. column) dimension of this matrix.

        OUTPUT: a polynomial matrix.

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
        from sage.matrix.constructor import Matrix

        # if d is an integer, make it a uniform list
        if not isinstance(d,list):
            d = [d]*m if row_wise else [d]*n

        # raise an error if d does not have the right length
        if row_wise and len(d) != m:
            raise ValueError("length of input precision list should be the " \
                                      + "row dimension of the input matrix")
        elif (not row_wise) and len(d) != n:
            raise ValueError("length of input precision list should be the " \
                                      + "column dimension of the input matrix")

        return Matrix(self.base_ring(), m, n, [[self[i,j].truncate(d[i])
            if row_wise else self[i,j].truncate(d[j])
            for j in range(n)] for i in range(m)])

    def shift(self, d, row_wise=True):
        r"""
        Return the matrix which is obtained from this matrix after shifting
        all its entries as specified by `d`.

        - if `d` is an integer, the shift is by `d` for all entries;
        - if `d` is a list $(d_1,\ldots,d_m)$ and ``row_wise`` is ``True``, all
          entries of the $i$th row are shifted by $d_i$ for each $i$;
        - if `d` is a list $(d_1,\ldots,d_n)$ and ``row_wise`` is ``False``,
          all entries of the $j$th column are shifted by $d_j$ for each $j$.

        Shifting by `d` means multiplying by the variable to the power `d`; if
        `d` is negative then terms of negative degree after shifting are
        discarded.

        INPUT:

        - ``d`` -- a list of integers, or an integer,

        - ``row_wise`` -- (optional, default: ``True``) boolean, if ``True``
          (resp. ``False``) then `d` should be a list of length equal to the
          row (resp. column) dimension of this matrix.

        OUTPUT: a polynomial matrix.

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
        from sage.matrix.constructor import Matrix

        # if d is an integer, make it a uniform list
        if not isinstance(d,list):
            d = [d]*m if row_wise else [d]*n

        # raise an error if d does not have the right length
        if row_wise and len(d) != m:
            raise ValueError("length of input shift list should be the " \
                                      + "row dimension of the input matrix")
        elif (not row_wise) and len(d) != n:
            raise ValueError("length of input shift list should be the " \
                                      + "column dimension of the input matrix")

        return Matrix(self.base_ring(), m, n, [[self[i,j].shift(d[i])
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
        - if ``degree`` is a list $(d_1,\ldots,d_m)$ and ``row_wise`` is
          ``True``, all entries of the $i$th row are reversed with respect to
          $d_i$ for each $i$;
        - if ``degree`` is a list $(d_1,\ldots,d_n)$ and ``row_wise`` is
          ``False``, all entries of the $j$th column are reversed with respect
          to $d_j$ for each $j$.

        INPUT:

        - ``degree`` -- (optional, default: ``None``) a list of nonnegative
          integers, or a nonnegative integer,

        - ``row_wise`` -- (optional, default: ``True``) boolean, if ``True``
          (resp. ``False``) then ``degree`` should be a list of length equal to
          the row (resp. column) dimension of this matrix.

        - ``entry_wise`` -- (optional, default: ``False``) boolean, if ``True``
          then the input ``degree`` and ``row_wise`` are ignored.

        OUTPUT: a polynomial matrix.

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

            sage: M.reverse(M.column_degrees(),row_wise=False)
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
            OverflowError: can't convert negative value to unsigned long

        .. SEEALSO::

            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.reverse` .
        """
        m = self.nrows()
        n = self.ncols()
        from sage.matrix.constructor import Matrix

        # if entry_wise, just return the matrix with all entries reversed
        if entry_wise:
            return Matrix(self.base_ring(), m, n, [[self[i,j].reverse()
                for j in range(n)] for i in range(m)])

        # if degree is None, make it the matrix degree
        if degree==None:
            degree = self.degree()
        # if degree is an integer, make it a uniform list
        if not isinstance(degree,list):
            degree = [degree]*m if row_wise else [degree]*n

        # raise an error if degree does not have the right length
        if row_wise and len(degree) != m:
            raise ValueError("length of input degree list should be the " \
                                      + "row dimension of the input matrix")
        elif (not row_wise) and len(degree) != n:
            raise ValueError("length of input degree list should be the " \
                                      + "column dimension of the input matrix")

        return Matrix(self.base_ring(), m, n, [[self[i,j].reverse(degree[i])
            if row_wise else self[i,j].reverse(degree[j])
            for j in range(n)] for i in range(m)])

    def row_degrees(self, shifts=None):
        r"""
        Return the (shifted) row degrees of this matrix.

        For a given polynomial matrix $M = (M_{i,j})_{i,j}$ with $m$ rows and
        $n$ columns, its row degrees is the tuple $(d_1,\ldots,d_m)$ where $d_i
        = \max_j(\deg(M_{i,j}))$ for $1\leq i \leq m$. Thus, $d_i=-1$ if
        the $i$-th row of $M$ is zero, and $d_i \geq 0$ otherwise.

        For given shifts $s_1,\ldots,s_n \in \ZZ$, the shifted row degrees of
        $M$ is $(d_1,\ldots,d_m)$ where $d_i = \max_j(\deg(M_{i,j})+s_j)$.
        Here, if the $i$-th row of $M$ is zero then $d_i
        =\min(s_1,\ldots,s_n)-1$; otherwise, $d_i$ is larger than this value.

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        OUTPUT: a list of integers.

        REFERENCES:

        - [Wol1974]_ (Section 2.5, without shifts), and [VBB1992]_ (Section 3).

        - Up to changes of signs, shifted row degrees coincide with the notion
          of *defect* commonly used in the rational approximation literature
          (see for example [Bec1992]_ ).

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0] ])
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

        The row degrees of an empty matrix ($0\times n$ or $m\times 0$) is
        not defined::
            
            sage: M = Matrix( pR, 0, 3 )
            sage: M.row_degrees()
            Traceback (most recent call last):
            ...
            ValueError: empty matrix does not have row degrees

            sage: M = Matrix( pR, 3, 0 )
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

        For a given polynomial matrix $M = (M_{i,j})_{i,j}$ with $m$ rows and
        $n$ columns, its column degrees is the tuple $(d_1,\ldots,d_n)$ where
        $d_j = \max_i(\deg(M_{i,j}))$ for $1\leq j \leq n$. Thus, $d_j=-1$ if
        the $j$-th column of $M$ is zero, and $d_j \geq 0$ otherwise.

        For given shifts $s_1,\ldots,s_m \in \ZZ$, the shifted column degrees of
        $M$ is $(d_1,\ldots,d_n)$ where $d_j = \max_i(\deg(M_{i,j})+s_i)$.
        Here, if the $j$-th column of $M$ is zero then $d_j =
        \min(s_1,\ldots,s_m)-1$; otherwise $d_j$ is larger than this value.

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        OUTPUT: a list of integers.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0] ])
            sage: M.column_degrees()
            [3, -1, 0]

            sage: M.column_degrees(shifts=[0,2])
            [5, -1, 0]

        A zero column in a polynomial matrix can be identified in the (shifted)
        column degrees as the entries equal to ``min(shifts)-1``::

            sage: M.column_degrees(shifts=[-2,1])
            [4, -3, -2]

        The column degrees of an empty matrix ($0\times n$ or $m\times 0$) is
        not defined::

            sage: M = Matrix( pR, 0, 3 )
            sage: M.column_degrees()
            Traceback (most recent call last):
            ...
            ValueError: empty matrix does not have column degrees

            sage: M = Matrix( pR, 3, 0 )
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

        Let $M$ be a univariate polynomial matrix in $\Bold{K}[x]^{m \times
        n}$. Working row-wise and without shifts, its leading matrix is the
        matrix in $\Bold{K}^{m \times n}$ formed by the leading coefficients of
        the entries of $M$ which reach the degree of the corresponding row.
  
        More precisely, if working row-wise, let $s_1,\ldots,s_n \in \ZZ$
        be a shift, and let $(d_1,\ldots,d_m)$ denote the shifted row degrees of
        $M$. Then, the shifted leading matrix of $M$ is the matrix in
        $\Bold{K}^{m \times n}$ whose entry $i,j$ is the coefficient of degree
        $d_i-s_j$ of the entry $i,j$ of $M$.

        If working column-wise, let $s_1,\ldots,s_m \in \ZZ$ be a shift,
        and let $(d_1,\ldots,d_n)$ denote the shifted column degrees of $M$.
        Then, the shifted leading matrix of $M$ is the matrix in $\Bold{K}^{m
        \times n}$ whose entry $i,j$ is the coefficient of degree $d_j-s_i$ of
        the entry $i,j$ of $M$.

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, ``True`` if
          working row-wise (see the class description).

        OUTPUT: a matrix over the base field.

        REFERENCES:
        
        [Wol1974]_ (Section 2.5, without shifts) and [VBB1992]_ (Section 3).

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0] ])
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
        from sage.matrix.constructor import Matrix
        if row_wise:
            row_degrees = self.row_degrees(shifts)
            if shifts is None:
                return Matrix([ [ self[i,j].leading_coefficient()
                    if self[i,j].degree() == row_degrees[i] else 0
                    for j in range(self.ncols()) ]
                    for i in range(self.nrows()) ])
            return Matrix([ [ self[i,j].leading_coefficient()
                if self[i,j].degree() + shifts[j] == row_degrees[i] else 0
                for j in range(self.ncols()) ]
                for i in range(self.nrows()) ])
        else:
            column_degrees = self.column_degrees(shifts)
            if shifts is None:
                return Matrix([ [ self[i,j].leading_coefficient()
                    if self[i,j].degree() == column_degrees[j] else 0
                    for j in range(self.ncols()) ]
                    for i in range(self.nrows()) ])
            return Matrix([ [ self[i,j].leading_coefficient()
                if self[i,j].degree() + shifts[i] == column_degrees[j] else 0
                for j in range(self.ncols()) ]
                for i in range(self.nrows()) ])

    def _is_empty_popov(self, row_wise=True, include_zero_vectors=True):
        r"""
        Assuming that this matrix is empty, that is, of dimensions $0\times n$
        or $m\times 0$, return a boolean indicating if it is in shifted Popov
        form. If zero vectors are allowed in shifted reduced forms, this always
        returns true. Otherwise, by convention and if working row-wise, for
        $n\geq 0$ the $0\times n$ matrix is in shifted Popov form for all
        shifts, and for $m>0$ the $m \times 0$ matrix is not in shifted Popov
        form for any shift. The convention is similar if working column-wise.

        The behaviour of this method for non-empty matrices is not defined.

        INPUT:

        - ``row_wise`` -- (optional, default: ``True``) boolean, ``True`` if
          one considers the row-wise shifted Popov form.

        - ``include_zero_vectors`` -- (optional, default: ``True``) boolean,
          ``False`` if one does not allow zero rows in row reduced forms (resp.
          zero columns in column reduced forms).

        OUTPUT: a boolean.

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

        An $m \times n$ univariate polynomial matrix $M$ is said to be in
        shifted row reduced form if it has $k$ nonzero rows with $k \leq n$ and
        its shifted leading matrix has rank $k$. Equivalently, when considering
        all the matrices obtained by left-multiplying $M$ by a unimodular
        matrix, then the shifted row degrees of $M$ -- once sorted in
        nondecreasing order -- is lexicographically minimal.

        Similarly, $M$ is said to be in shifted column reduced form if it has
        $k$ nonzero columns with $k \leq m$ and its shifted leading matrix has
        rank $k$.

        Sometimes, one forbids $M$ to have zero rows (resp. columns) in the
        above definitions; an optional parameter allows one to adopt this more
        restrictive setting.

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, ``True`` if
          working row-wise (see the class description).

        - ``include_zero_vectors`` -- (optional, default: ``True``) boolean,
          ``False`` if one does not allow zero rows in row reduced forms (resp.
          zero columns in column reduced forms).

        OUTPUT: a boolean value.

        REFERENCES:
        
        [Wol1974]_ (Section 2.5, without shifts) and [VBB1992]_ (Section 3).

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0] ])
            sage: M.is_reduced()
            False

            sage: M.is_reduced(shifts=[0,1,2])
            True

            sage: M.is_reduced(shifts=[2,0], row_wise=False)
            True

            sage: M.is_reduced(shifts=[2,0], row_wise=False,
            ....:                           include_zero_vectors=False)
            False

            sage: M = Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0], [0, 1, 0] ])
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

        If working row-wise, for a given shift $s_1,\ldots,s_n \in
        \ZZ$, taken as $(0,\ldots,0)$ by default, and a row vector of
        univariate polynomials $[p_1,\ldots,p_n]$, the leading position of
        this vector is the index $j$ of the rightmost nonzero entry $p_j$ such
        that $\deg(p_j) + s_j$ is equal to the shifted row degree of the vector.
        Then the pivot degree of the vector is the degree $\deg(p_j)$.
        
        For the zero row, both the leading positions and degree are $-1$.  For
        a $m \times n$ polynomial matrix, the leading positions and pivot
        degrees are the two lists containing the leading positions and the
        pivot degrees of its rows.

        The definition is similar if working column-wise (instead of rightmost
        nonzero entry, we choose the bottommost nonzero entry).

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, ``True`` if
          working row-wise (see the class description).

        - ``return_degree`` -- (optional, default: ``False``) boolean, ``True``
          implies that the pivot degrees are returned.

        OUTPUT: a list of integers if ``return_degree=False``; a pair of lists
        of integers otherwise.

        REFERENCES:
        
        [Kai1980]_ (Section 6.7.2, without shifts).

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0] ])
            sage: M.leading_positions()
            [0, 0]

            sage: M.leading_positions(return_degree=True)
            ([0, 0], [1, 3])

            sage: M.leading_positions(shifts=[0,5,2], return_degree=True)
            ([2, 0], [0, 3])

            sage: M.leading_positions(row_wise=False, return_degree=True)
            ([1, -1, 0], [3, -1, 0])

            sage: M.leading_positions(shifts=[1,2], row_wise=False,
            ....:   return_degree=True)
            ([1, -1, 0], [3, -1, 0])

        In case several entries in the row (resp. column) reach the shifted row
        (resp. column) degree, the leading position is chosen as the rightmost
        (resp. bottommost) such entry::

            sage: M.leading_positions(shifts=[0,5,1],return_degree=True)
            ([2, 0], [0, 3])

            sage: M.leading_positions(shifts=[2,0], row_wise=False,return_degree=True)
            ([1, -1, 0], [3, -1, 0])

        The leading positions and pivot degrees of an empty matrix ($0\times n$
        or $m\times 0$) is not defined::

            sage: M = Matrix( pR, 0, 3 )
            sage: M.leading_positions()
            Traceback (most recent call last):
            ...
            ValueError: empty matrix does not have leading positions

            sage: M.leading_positions(row_wise=False)
            Traceback (most recent call last):
            ...
            ValueError: empty matrix does not have leading positions

            sage: M = Matrix( pR, 3, 0 )
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
        (resp. columns) are pairwise distinct (for the ordered weak Popov form,
        these positions must be strictly increasing, except for the possibly
        repeated -1 entries which are at the end).

        Sometimes, one forbids $M$ to have zero rows (resp. columns) in the
        above definitions; an optional parameter allows one to adopt this more
        restrictive setting.

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, ``True`` if
          working row-wise (see the class description).

        - ``ordered`` -- (optional, default: ``False``) boolean, ``True`` if
          checking for an ordered weak Popov form.

        - ``include_zero_vectors`` -- (optional, default: ``True``) boolean,
          ``False`` if one does not allow zero rows (resp. zero columns) in
          (ordered) weak Popov forms.

        OUTPUT: a boolean.

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

            sage: N = M.with_swapped_rows(1,2)
            sage: N.is_weak_popov()
            True
            sage: N.is_weak_popov(ordered=True)
            False

        Shifts and orientation (row-wise or column-wise) are supported::

            sage: M.is_weak_popov(shifts=[2,3,1])
            False

            sage: M.is_weak_popov(shifts=[0,2,0],row_wise=False,ordered=True)
            True

        Rectangular matrices are supported::

            sage: M = Matrix([
            ....:    [  x^3+5*x^2+5*x+1,       5,       6*x+4,         0],
            ....:    [      6*x^2+3*x+1,       1,           2,         0],
            ....:    [2*x^3+4*x^2+6*x+4, 5*x + 1, 2*x^2+5*x+5, x^2+5*x+6]
            ....:     ])
            sage: M.is_weak_popov(shifts=[0,2,1,3])
            True

            sage: M.is_weak_popov(shifts=[0,2,1,3],ordered=True)
            True

        Zero rows (resp. columns) can be forbidden::

            sage: M = Matrix([
            ....:   [      6*x+4,       0,             5*x+1, 0],
            ....:   [          2, 5*x + 1,       6*x^2+3*x+1, 0],
            ....:   [2*x^2+5*x+5,       1, 2*x^3+4*x^2+6*x+4, 0] 
            ....:   ])
            sage: M.is_weak_popov(shifts=[2,1,0], row_wise=False, ordered=True)
            True

            sage: M.is_weak_popov(shifts=[2,1,0], row_wise=False,
            ....:    include_zero_vectors=False)
            False

        .. SEEALSO::

            :meth:`weak_popov_form` .
        """
        self._check_shift_dimension(shifts,row_wise)
        if self.ncols() == 0 or self.nrows() == 0:
            return self._is_empty_popov(row_wise)
        leading_positions = self.leading_positions(shifts, row_wise)
        # here, it will be convenient to have leading position
        # larger than ncols for zero/empty rows
        leading_positions = [pos if pos>=0 else self.ncols() + 1 for pos in leading_positions]
        # leading positions should not have duplicates, which is equivalent to:
        # once sorted, it doesn't contain a pair of equal successive entries
        if not ordered:
            leading_positions.sort()
        # check that there is no zero vector, if it is forbidden
        if leading_positions[-1] > self.ncols() and not include_zero_vectors:
            return False
        # now leading_positions is nondecreasing: it remains to test whether
        # it is strictly increasing (at least until the zero rows part)
        for index,next_leading_position in enumerate(leading_positions[1:]):
            if next_leading_position <= self.ncols() and \
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

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, ``True`` if
          working row-wise (see the class description).

        - ``up_to_permutation`` -- (option, default: ``False``) boolean,
          ``True`` if testing Popov form up to row permutation (if working
          row-wise).

        - ``include_zero_vectors`` -- (optional, default: ``True``) boolean,
          ``False`` if one does not allow zero rows (resp. zero columns) in
          Popov forms.

        OUTPUT: a boolean.

        REFERENCES:
        
        For the square case, without shifts: [Pop1972]_ and [Kai1980]_ (Section
        6.7.2). For the general case: [BLV2006]_ .

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [ [x^4+6*x^3+4*x+4, 3*x+6,     3  ],
            ....:                  [x^2+6*x+6,       x^2+5*x+5, 2  ],
            ....:                  [3*x,             6*x+5,     x+5] ])
            sage: M.is_popov()
            True

            sage: M.is_popov(shifts=[0,1,2])
            True

            sage: M[:,:2].is_popov()
            False

            sage: M[:2,:].is_popov(shifts=[0,1,2])
            True

            sage: M = Matrix(pR, [ [x^4+3*x^3+x^2+2*x+6, x^3+5*x^2+5*x+1],
            ....:                  [6*x+1,               x^2+4*x+1      ],
            ....:                  [6,                   6              ] ])
            sage: M.is_popov(row_wise=False)
            False

            sage: M.is_popov(shifts=[0,2,3], row_wise=False)
            True

        One can forbid zero rows (or columns if not working row-wise)::

            sage: N = Matrix(pR, [ [x^4+3*x^3+x^2+2*x+6, 6*x+1     ],
            ....:                  [5*x^2+5*x+1,         x^2+4*x+1 ],
            ....:                  [0,                   0         ] ])

            sage: N.is_popov()
            True

            sage: N.is_popov(include_zero_vectors=False)
            False

        One can verify Popov form up to row permutation (or column permutation
        if not working row-wise)::

            sage: M.swap_columns(0,1)
            sage: M.is_popov(shifts=[0,2,3], row_wise=False)
            False

            sage: M.is_popov(shifts=[0,2,3], row_wise=False,
            ....:   up_to_permutation=True)
            True

            sage: N.swap_rows(0,2)

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

        Note that, for any integer $d$ strictly greater than all degrees
        appearing in the Hermite form, then the Hermite form coincides with the
        shifted Popov form with the shifts $(0,d,2d,\ldots,(n-1)d)$, where $n$
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

        - ``row_wise`` -- (optional, default: ``True``) boolean, ``True`` if
          working row-wise (see the class description).

        - ``lower_echelon`` -- (optional, default: ``False``) boolean,
          ``False`` if working with upper triangular Hermite forms, ``True`` if
          working with lower triangular Hermite forms.

        - ``include_zero_vectors`` -- (optional, default: ``True``) boolean,
          ``False`` if one does not allow zero rows (resp. zero columns) in
          Hermite forms.

        OUTPUT: a boolean.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [ [x^4+6*x^3+4*x+4, 3*x+6,     3  ],
            ....:                  [0,               x^2+5*x+5, 2  ],
            ....:                  [0,               0,         x+5] ])

            sage: M.is_hermite()
            True
            sage: M.is_hermite(row_wise=False)
            True
            sage: M.is_hermite(row_wise=False, lower_echelon=True)
            False

            sage: N = Matrix(pR, [ [x+5, 0,               0        ],
            ....:                  [2,   x^4+6*x^3+4*x+4, 0        ],
            ....:                  [3,   3*x^3+6,         x^2+5*x+5] ])
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

    def weak_popov_form(self, transformation=False, shifts=None):
        r"""
        Return a (row-wise) weak Popov form of this matrix.

        A polynomial matrix is said to be in (row-wise) weak Popov form if the
        (shifted) leading positions of its nonzero rows are pairwise distinct.
        The leading position of a row is the right-most position whose entry has
        the maximal degree in the row, see :meth:`leading_positions`. See the
        class description for an introduction to shifts.

        The weak Popov form is non-canonical, so an input matrix have many weak
        Popov forms (for any given shifts).

        INPUT:

        - ``transformation`` -- (optional, default: ``False``). If this
          is ``True``, the transformation matrix `U` will be returned as well:
          this is a unimodular matrix over `\Bold{K}[x]` such that ``self``
          equals `UW`, where `W` is the output matrix.

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        OUTPUT:

        - A polynomial matrix `W` which is a weak Popov form of ``self`` if
          ``transformation=False``; otherwise two polynomial matrices `W, U`
          such that `UA = W` and `W` is in weak Popov form and `U` is unimodular
          where `A` is ``self``.

        ALGORITHM:

        This method implements the Mulders-Storjohann algorithm of [MS2003]_.

        EXAMPLES::

            sage: PF.<x> = GF(11)[]
            sage: A = matrix(PF,[[1,  3*x^9 + x^2 + 1 ],
            ....:                [0,         x^11 + x ]])
            sage: W, U = A.weak_popov_form(transformation=True); W
            [              8*x^7 + 3*x^5 + 1    9*x^6 + 3*x^5 + 2*x^4 + x^2 + 1]
            [                          7*x^2                  7*x^4 + 7*x^2 + x]
            sage: U * A == W
            True
            sage: W.is_weak_popov()
            True
            sage: U.is_invertible()
            True

        Demonstrating shifts::

            sage: A.weak_popov_form(shifts=[2, 0])
            [              8*x^7 + 1 8*x^7 + 9*x^6 + x^2 + 1]
            [                  7*x^2       7*x^4 + 7*x^2 + x]
            sage: A.weak_popov_form(shifts=[10, 0]) == A
            True

        A zero matrix will return itself::

            sage: Z = matrix(PF,2,2)
            sage: Z.weak_popov_form()
            [0 0]
            [0 0]

        .. SEEALSO::

            :meth:`is_weak_popov` ,
            :meth:`reduced_form` ,
            :meth:`hermite_form` .
        """
        self._check_shift_dimension(shifts,row_wise=True)
        M = self.__copy__()
        U = M._weak_popov_form(transformation=transformation, shifts=shifts)
        M.set_immutable()
        if transformation:
            U.set_immutable()
        return (M,U) if transformation else M

    def _weak_popov_form(self, transformation=False, shifts=None):
        """
        Transform this matrix in-place into weak Popov form.

        EXAMPLES::

            sage: F.<a> = GF(2^4,'a')
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

        if shifts and len(shifts) != M.ncols():
            raise ValueError("the number of shifts must equal the number of columns")

        # initialise to_row and conflicts list
        to_row = [[] for i in range(n)]
        conflicts = []
        for i in range(m):
            bestp = -1
            best = -1
            for c in range(n):
                d = M.get_unsafe(i,c).degree()

                if shifts and d >= 0 :
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

    def reduced_form(self, transformation=None, shifts=None, row_wise=True):
        r"""
        Return a row reduced form of this matrix (resp. a column reduced form
        if the optional parameter `row_wise` is set to `False`).

        An $m \times n$ univariate polynomial matrix $M$ is said to be in
        (shifted) row reduced form if it has $k$ nonzero rows with $k \leq n$
        and its (shifted) leading matrix has rank $k$. See :meth:`is_reduced`
        for more information.

        A row reduced form is non-canonical and a given matrix has many row
        reduced forms; this method returns just one.

        INPUT:

        - ``transformation`` -- (optional, default: ``False``). If this
          is ``True``, the transformation matrix `U` will be returned as well:
          this is a unimodular matrix over `\Bold{K}[x]` such that ``self``
          equals `UR`, where `R` is the output matrix.

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, ``True`` if
          working row-wise (see the class description).

        OUTPUT:

        - A polynomial matrix `R` which is a reduced form of ``self`` if
          ``transformation=False``; otherwise two polynomial matrices `R, U`
          such that `UA = R` and `R` is reduced and `U` is unimodular where `A`
          is ``self``.

        EXAMPLES::

            sage: pR.<x> = GF(3)[]
            sage: A = matrix(pR,3,[x,   x^2, x^3,
            ....:                  x^2, x^1, 0,
            ....:                  x^3, x^3, x^3])
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

        If the matrix is an `n \times 1` matrix with at least one non-zero entry,
        `R` has a single non-zero entry and that entry is a scalar multiple of
        the greatest-common-divisor of the entries of the matrix::

            sage: A = matrix([[x*(x-1)*(x+1)],[x*(x-2)*(x+2)],[x]])
            sage: R = A.reduced_form()
            sage: R
            [0]
            [0]
            [x]

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

            sage: Fq.<a> = GF(2^3)
            sage: pR.<x> = Fq[]
            sage: A = matrix(pR, [[x^2+a,  x^4+a],
            ....:                  [  x^3,  a*x^4]])
            sage: W,U = A.reduced_form(transformation=True)
            sage: W,U
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
            :meth:`weak_popov_form` .
        """
        self._check_shift_dimension(shifts,row_wise)
        if not row_wise:
            return self.T.reduced_form(transformation, shifts, row_wise=True).T
        return self.weak_popov_form(transformation, shifts)

    def hermite_form(self, include_zero_rows=True, transformation=False):
        """
        Return the Hermite form of this matrix.

        The Hermite form is also normalized, i.e., the pivot polynomials
        are monic.

        INPUT:

        - ``include_zero_rows`` -- boolean (default: ``True``); if ``False``,
          the zero rows in the output matrix are deleted

        - ``transformation`` -- boolean (default: ``False``); if ``True``,
          return the transformation matrix

        OUTPUT:

        - the Hermite normal form `H` of this matrix `A`

        - (optional) transformation matrix `U` such that `UA = H`
 
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
            sage: H, U = A.hermite_form(transformation=True, include_zero_rows=True); H, U
            (
            [  x   1 2*x]  [0 4]
            [  0   0   0], [5 1]
            )
            sage: U * A == H
            True
            sage: H, U = A.hermite_form(transformation=True, include_zero_rows=False)
            sage: U * A
            [  x   1 2*x]
            sage: U * A == H
            True

        .. SEEALSO::
        
            :meth:`is_hermite` .
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

    def is_minimal_approximant_basis(self,
            pmat,
            order,
            shifts=None,
            row_wise=True,
            normal_form=False):
        r"""
        Return ``True`` if and only if this matrix is an approximant basis in
        ``shifts``-ordered weak Popov form for the polynomial matrix ``pmat``
        at order ``order``.
        
        If ``normal_form`` is ``True``, then the polynomial matrix must
        furthermore be in ``shifts``-Popov form. An error is raised if the
        input dimensions are not sound. If a single integer is provided for
        ``order``, then it is interpreted as a list of repeated integers with
        this value. (See :meth:`minimal_approximant_basis` for definitions and 
        more details.)

        INPUT:

        - ``pmat`` -- a polynomial matrix.

        - ``order`` -- a list of positive integers, or a positive integer.

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, if ``True``
          then the basis considered row-wise and operates on the left of
          ``pmat``; otherwise it is column-wise and operates on the right of
          ``pmat``.

        - ``normal_form`` -- (optional, default: ``False``) boolean, if
          ``True`` then checks for a basis in ``shifts``-Popov form.

        OUTPUT: a boolean.

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
            sage: pmat = Matrix(pR, 5, 1, [ \
                    pR([35,  0, 41, 87,  3, 42, 22, 90]), \
                    pR([80, 15, 62, 87, 14, 93, 24,  0]), \
                    pR([42, 57, 90, 87, 22, 80, 71, 53]), \
                    pR([37, 72, 74,  6,  5, 75, 23, 47]), \
                    pR([36, 10, 74,  1, 29, 44, 87, 74]) ])
            sage: appbas = Matrix(pR, [ \
                   [x+47,   57, 58*x+44,     9*x+23,      93*x+76], \
                   [  15, x+18, 52*x+23,     15*x+58,     93*x+88], \
                   [  17,   86, x^2+77*x+16, 76*x+29,     90*x+78], \
                   [  44,   36, 3*x+42,      x^2+50*x+26, 85*x+44], \
                   [   2,   22, 54*x+94,     73*x+24,     x^2+2*x+25] ])
            sage: appbas.is_minimal_approximant_basis(pmat,\
                    order, shifts, row_wise=True, normal_form=True)
            True

        The matrix `x^8 \mathrm{Id}_5` is square, nonsingular, in Popov form,
        and its rows are approximants for ``pmat`` at order 8. However, it is
        not an approximant basis since its rows generate a module strictly
        contained in the set of approximants for ``pmat`` at order 8::

            sage: (x^8*Matrix.identity(pR, 5)).is_minimal_approximant_basis(\
                                                                    pmat, 8)
            False

        Since ``pmat`` is a single column, with nonzero constant coefficient,
        its column-wise approximant bases at order 8 are all `1\times 1`
        matrices `[c x^8]` for some nonzero field element `c`::

            sage: Matrix(pR, [x^8]).is_minimal_approximant_basis(pmat, \
                    8, row_wise=False, normal_form=True)
            True

        Exceptions are raised if input dimensions are not sound::

            sage: appbas.is_minimal_approximant_basis(pmat, [8,8], shifts)
            Traceback (most recent call last):
            ...
            ValueError: order length should be the column dimension 
                        of the input matrix

            sage: appbas.is_minimal_approximant_basis(pmat, \
                    order, shifts, row_wise=False)
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
            raise ValueError('shifts length should be the row dimension of' \
                                                      + ' the input matrix')
        elif (not row_wise) and len(shifts) != n:
            raise ValueError('shifts length should be the column dimension' \
                                                   + ' of the input matrix')

        # set default order / check order dimension
        if not isinstance(order,list):
            order = [order]*n if row_wise else [order]*m

        if row_wise and len(order) != n:
            raise ValueError("order length should be the column dimension" \
                                                  + " of the input matrix")
        elif (not row_wise) and len(order) != m:
            raise ValueError("order length should be the row dimension of" \
                                                     + " the input matrix")

        # raise an error if self does not have the right dimension
        if row_wise and self.ncols() != m:
            raise ValueError("column dimension should be the row dimension" \
                                                    + " of the input matrix")
        elif (not row_wise) and self.nrows() != n:
            raise ValueError("row dimension should be the column dimension" \
                                                    + " of the input matrix")

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

        - ``order`` -- a list of positive integers, or a positive integer.

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, if ``True``
          then the output basis is considered row-wise and operates on the left
          of ``self``; otherwise it is column-wise and operates on the right
          of ``self``.

        - ``normal_form`` -- (optional, default: ``False``) boolean, if
          ``True`` then the output basis is in ``shifts``-Popov form.

        OUTPUT: a polynomial matrix.

        ALGORITHM:

        The implementation is inspired from the iterative algorithms described
        in [VBB1992]_ and [BL1994]_ ; for obtaining the normal form, it relies
        directly on Lemmas 3.3 and 4.1 in [JNSV2016]_ .

        EXAMPLES::

            sage: pR.<x> = GF(7)[]

            sage: order = [4, 3]; shifts = [-1, 2, 0]
            sage: F = Matrix(pR, [[5*x^3 + 4*x^2 + 4*x + 6, 5*x^2 + 4*x + 1], \
                                  [        2*x^2 + 2*x + 3, 6*x^2 + 6*x + 3], \
                                  [4*x^3         +   x + 1, 4*x^2 + 2*x + 3] ])
            sage: P = F.minimal_approximant_basis(order, shifts)
            sage: P.is_minimal_approximant_basis(F, order, shifts)
            True

        By default, the computed basis is not required to be in normal form
        (and will not be except in rare special cases)::

            sage: P.is_minimal_approximant_basis(F, order, shifts, \
                                                    normal_form=True)
            False
            sage: P = F.minimal_approximant_basis(order, shifts, normal_form=True)
            sage: P.is_minimal_approximant_basis(F, order, shifts, \
                                                    normal_form=True)
            True

        If shifts are not specified, they are chosen as uniform `[0,\ldots,0]`
        by default. Besides, if the orders are all the same, one can rather
        give a single integer::

            sage: F.minimal_approximant_basis(3) == \
                    F.minimal_approximant_basis([3,3], shifts=None)
            True

        One can work column-wise by specifying ``row_wise=False``::

            sage: P = F.minimal_approximant_basis([5,2,2], [0,1], row_wise=False)
            sage: P.is_minimal_approximant_basis(F, [5,2,2], \
                                shifts=[0,1], row_wise=False)
            True
            sage: F.minimal_approximant_basis(3, row_wise=True) == \
                F.transpose().minimal_approximant_basis(3, row_wise=False).transpose()
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
                P,rdeg = self.transpose()._approximant_basis_iterative( \
                                                order, degree_shifts)
                P = P.transpose()
                # right-multiply by inverse of leading matrix
                lmat = P.leading_matrix(shifts=degree_shifts,row_wise=False)
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

        - ``order`` -- a list of positive integers.

        - ``shifts`` -- a list of integers.

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
            sage: pmat = Matrix(pR, [[5*x^3 + 4*x^2 + 4*x + 6, 5, 4], \
                    [2*x^3 + 2*x^2 + 2*x + 3, 6, 6*x + 3]])
            sage: appbas,rdeg = pmat._approximant_basis_iterative(order, \
                                                                    shifts)
            sage: appbas.is_minimal_approximant_basis(pmat, order, shifts)
            True

        The returned list is the shifted row degrees of ``appbas``::

            sage: rdeg == appbas.row_degrees(shifts)
            True

        Approximant bases for the zero matrix are all constant unimodular
        matrices; in fact, this algorithm returns the identity::

            sage: pmat = Matrix(pR, 3, 2)
            sage: appbas,rdeg = pmat._approximant_basis_iterative([2,5], \
                                                                [5,0,-4])
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
        Return ``True`` if and only if this matrix is a left kernel basis in
        ``shifts``-ordered weak Popov form for the polynomial matrix ``pmat``.

        If ``normal_form`` is ``True``, then the kernel basis must furthermore
        be in ``shifts``-Popov form. An error is raised if the input dimensions
        are not sound.

        INPUT:

        - ``pmat`` -- a polynomial matrix.

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, if ``True``
          then the basis is considered row-wise and operates on the left of
          ``pmat``; otherwise it is column-wise and operates on the right of
          ``pmat``.

        - ``normal_form`` -- (optional, default: ``False``) boolean, if
          ``True`` then checks for a basis in ``shifts``-Popov form.

        OUTPUT: a boolean.

        ALGORITHM:

        Verification that the matrix has full rank and is in shifted weak
        Popov form is done via :meth:`is_weak_popov`; verification that the
        matrix is a left kernel basis is done by checking that the rank is
        correct, that the product is indeed zero, and that the matrix is
        saturated, i.e. it has unimodular column bases (see Lemma 6.10 of
        https://arxiv.org/pdf/1807.01272.pdf for details).

        EXAMPLES::

            sage: pR.<x> = GF(97)[]
            sage: pmat = Matrix(pR, [[1],[x],[x**2]])

            sage: kerbas = Matrix(pR, [[x,-1,0],[0,x,-1]])
            sage: kerbas.is_minimal_kernel_basis(pmat)
            True

        A matrix in Popov form which has the right rank, all rows in the
        kernel, but does not generate the kernel::

            sage: kerbas = Matrix(pR, [[x**2,0,-1],[0,x,-1]])
            sage: kerbas.is_minimal_kernel_basis(pmat)
            False

        Shifts and right kernel bases are supported (with ``row_wise``), and one can test whether the kernel basis is normalized in shifted-Popov form (with ``normal_form``)::

            sage: kerbas = Matrix(pR, [[-x,-x**2],[1,0],[0,1]])
            sage: kerbas.is_minimal_kernel_basis(pmat.transpose(),row_wise=False,normal_form=True,shifts=[0,1,2])
            True
        """
        m = pmat.nrows()
        n = pmat.ncols()

        # set default shifts / check shifts dimension
        if shifts is None:
            shifts = [0] * m if row_wise else [0] * n
        elif row_wise and len(shifts) != m:
            raise ValueError('shifts length should be the row dimension of' \
                                                      + ' the input matrix')
        elif (not row_wise) and len(shifts) != n:
            raise ValueError('shifts length should be the column dimension' \
                                                   + ' of the input matrix')

        # raise an error if self does not have the right dimension
        if row_wise and self.ncols() != m:
            raise ValueError("column dimension should be the row dimension" \
                                                    + " of the input matrix")
        elif (not row_wise) and self.nrows() != n:
            raise ValueError("row dimension should be the column dimension" \
                                                    + " of the input matrix")

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
            # TODO replace by popov_form (likely more efficient) once it is written
        else:
            hnf = self.hermite_form(include_zero_rows=False)
            # TODO replace by popov_form (likely more efficient) once it is written
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

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, if ``True``
          then the output basis considered row-wise and operates on the left
          of ``self``; otherwise it is column-wise and operates on the right
          of ``self``.

        - ``normal_form`` -- (optional, default: ``False``) boolean, if
          ``True`` then the output basis is in ``shifts``-Popov form.

        OUTPUT: a polynomial matrix.

        ALGORITHM: uses minimal approximant basis computation at a
        sufficiently large order so that the approximant basis contains
        a kernel basis as a submatrix.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: pmat = Matrix([[(x+1)*(x+3)],[(x+1)*(x+3)+1]])
            sage: pmat.minimal_kernel_basis()
            [6*x^2 + 3*x + 3   x^2 + 4*x + 3]

            sage: pmat = Matrix([[(x+1)*(x+3)],[(x+1)*(x+4)]])
            sage: pmat.minimal_kernel_basis()
            [6*x + 3   x + 3]

            sage: pmat.minimal_kernel_basis(row_wise=False)
            []

            sage: pmat = Matrix(pR, [[1,x,x**2]])
            sage: pmat.minimal_kernel_basis(row_wise=False,normal_form=True)
            [x 0]
            [6 x]
            [0 6]

            sage: pmat.minimal_kernel_basis(row_wise=False,normal_form=True,shifts=[0,1,2])
            [  6*x 6*x^2]
            [    1     0]
            [    0     1]
        """
        m = self.nrows()
        n = self.ncols()
        d = self.degree()

        # set default shifts / check shifts dimension
        if shifts is None:
            shifts = [0] * m if row_wise else [0] * n
        elif row_wise and len(shifts) != m:
            raise ValueError('shifts length should be the row dimension')
        elif (not row_wise) and len(shifts) != n:
            raise ValueError('shifts length should be the column dimension')

        # compute kernel basis
        if row_wise:
            if d is -1: # matrix is zero
                from sage.matrix.constructor import Matrix
                return Matrix.identity(self.base_ring(), m, m)

            if m <= n and self.constant_matrix().rank() == m:
                # early exit: kernel is empty
                from sage.matrix.constructor import Matrix
                return Matrix(self.base_ring(), 0, m)

            # degree bounds on the kernel basis
            degree_bound = min(m,n)*d+max(shifts)
            degree_bounds = [degree_bound - shifts[i] for i in range(m)]

            # orders for approximation
            orders = self.column_degrees(degree_bounds)
            for i in range(n): orders[i] = orders[i]+1

            # compute approximant basis and retrieve kernel rows
            P = self.minimal_approximant_basis(orders,shifts,True,normal_form)
            row_indices = []
            for i in range(m):
                if P[i,i].degree() + shifts[i] <= degree_bound:
                    row_indices.append(i)
            return P[row_indices,:]

        else:
            if d is -1: # matrix is zero
                from sage.matrix.constructor import Matrix
                return Matrix.identity(self.base_ring(), n, n)

            if n <= m and self.constant_matrix().rank() == n:
                # early exit: kernel is empty
                from sage.matrix.constructor import Matrix
                return Matrix(self.base_ring(), n, 0)

            # degree bounds on the kernel basis
            degree_bound = min(m,n)*d+max(shifts)
            degree_bounds = [degree_bound - shifts[i] for i in range(n)]

            # orders for approximation
            orders = self.row_degrees(degree_bounds)
            for i in range(m): orders[i] = orders[i]+1

            # compute approximant basis and retrieve kernel columns
            P = self.minimal_approximant_basis(orders,shifts,False,normal_form)
            column_indices = []
            for j in range(n):
                if P[j,j].degree() + shifts[j] <= degree_bound:
                    column_indices.append(j)
            return P[:,column_indices]
