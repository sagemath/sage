"""
Dense matrices over univariate polynomials over fields

The implementation inherits from Matrix_generic_dense but some algorithms
are optimized for polynomial matrices.

AUTHORS:

- Kwankyu Lee (2016-12-15): initial version with code moved from other files.

- Johan Rosenkilde (2017-02-07): added weak_popov_form()

"""
#*****************************************************************************
#       Copyright (C) 2016 Kwankyu Lee <ekwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.matrix.matrix2 cimport Matrix

cdef class Matrix_polynomial_dense(Matrix_generic_dense):
    """
    Dense matrix over a univariate polynomial ring over a field.
    """

    def is_weak_popov(self):
        r"""
        Return ``True`` if the matrix is in weak Popov form.

        OUTPUT:

        A matrix over a polynomial ring is in weak Popov form if all
        leading positions are different [MS2003]_. A leading position
        is the position `i` in a row with the highest degree; in case of tie,
        the maximal `i` is used (i.e. furthest to the right).

        EXAMPLES:

        A matrix with the same leading position in two rows is not in weak
        Popov form::

            sage: PF = PolynomialRing(GF(2^12,'a'),'x')
            sage: A = matrix(PF,3,[x,   x^2, x^3,\
                                   x^2, x^2, x^2,\
                                   x^3, x^2, x    ])
            sage: A.is_weak_popov()
            False

        If a matrix has different leading positions, it is in weak Popov
        form::

            sage: B = matrix(PF,3,[1,    1,  x^3,\
                                   x^2,  1,  1,\
                                   1,x^  2,  1  ])
            sage: B.is_weak_popov()
            True

        Weak Popov form is not restricted to square matrices::

            sage: PF = PolynomialRing(GF(7),'x')
            sage: D = matrix(PF,2,4,[x^2+1, 1, 2, x,\
                                     3*x+2, 0, 0, 0 ])
            sage: D.is_weak_popov()
            False

        Even a matrix with more rows than columns can still be in weak Popov
        form::

            sage: E = matrix(PF,4,2,[4*x^3+x, x^2+5*x+2,\
                                     0,       0,\
                                     4,       x,\
                                     0,       0         ])
            sage: E.is_weak_popov()
            True

        A matrix with fewer columns than non-zero rows is never in weak
        Popov form::

            sage: F = matrix(PF,3,2,[x^2,   x,\
                                     x^3+2, x,\
                                     4,     5])
            sage: F.is_weak_popov()
            False

        TESTS:

        Verify tie breaking by selecting right-most index::

            sage: F = matrix(PF,2,2,[x^2, x^2,\
                                     x,   5   ])
            sage: F.is_weak_popov()
            True

        .. SEEALSO::

            - :meth:`weak_popov_form <sage.matrix.matrix_polynomial_dense.weak_popov_form>`

        AUTHOR:

        - David Moedinger (2014-07-30)
        """
        t = set()
        for r in range(self.nrows()):
            max = -1
            for c in range(self.ncols()):
                if self[r, c].degree() >= max:
                    max = self[r, c].degree()
                    p = c
            if not max == -1:
                if p in t:
                    return False
                t.add(p)
        return True

    def weak_popov_form(self, transformation=False, shifts=None):
        r"""
        Return a weak Popov form of the matrix.

        A matrix is in weak Popov form if the leading positions of the nonzero
        rows are all different. The leading position of a row is the right-most
        position whose entry has the maximal degree in the row.

        The weak Popov form is non-canonical, so an input matrix have many weak
        Popov forms.

        INPUT:

        - ``transformation`` -- boolean (default: ``False``) If ``True``, the
          transformation matrix is returned together with the weak Popov form.

        - ``shifts`` -- (default: ``None``) A tuple or list of integers
          `s_1, \ldots, s_n`, where `n` is the number of columns of the matrix.
          If given, a "shifted weak Popov form" is computed, i.e. such that the
          matrix `A\,\mathrm{diag}(x^{s_1}, \ldots, x^{s_n})` is in weak Popov
          form, where `\mathrm{diag}` denotes a diagonal matrix.

        ALGORITHM:

        This method implements the Mulders-Storjohann algorithm of [MS2003]_.

        EXAMPLES::

            sage: F.<a> = GF(2^4,'a')
            sage: PF.<x> = F[]
            sage: A = matrix(PF,[[1,  a*x^17 + 1 ],\
                                 [0,  a*x^11 + a^2*x^7 + 1 ]])
            sage: M, U = A.weak_popov_form(transformation=True)
            sage: U * A == M
            True
            sage: M.is_weak_popov()
            True
            sage: U.is_invertible()
            True

        A zero matrix will return itself::

            sage: Z = matrix(PF,5,3)
            sage: Z.weak_popov_form()
            [0 0 0]
            [0 0 0]
            [0 0 0]
            [0 0 0]
            [0 0 0]

        Shifted weak popov form is computed if ``shifts`` is given::

            sage: PF.<x> = QQ[]
            sage: A = matrix(PF,3,[x,   x^2, x^3,\
                                   x^2, x^1, 0,\
                                   x^3, x^3, x^3])
            sage: A.weak_popov_form()
            [        x       x^2       x^3]
            [      x^2         x         0]
            [  x^3 - x x^3 - x^2         0]
            sage: H,U = A.weak_popov_form(transformation=True, shifts=[16,8,0])
            sage: H
            [               x              x^2              x^3]
            [               0         -x^2 + x       -x^4 + x^3]
            [               0                0 -x^5 + x^4 + x^3]
            sage: U * A == H
            True

        .. SEEALSO::

            :meth:`is_weak_popov <sage.matrix.matrix_polynomial_dense.is_weak_popov>`
        """
        M = self.__copy__()
        U = M._weak_popov_form(transformation=transformation, shifts=shifts)
        M.set_immutable()
        if transformation:
            U.set_immutable()
        return (M,U) if transformation else M

    def _weak_popov_form(self, transformation=False, shifts=None):
        """
        Transform the matrix in place into weak Popov form.

        EXAMPLES::

            sage: F.<a> = GF(2^4,'a')
            sage: PF.<x> = F[]
            sage: A = matrix(PF,[[1,  a*x^17 + 1 ],\
                                 [0,  a*x^11 + a^2*x^7 + 1 ]])
            sage: M = A.__copy__()
            sage: U = M._weak_popov_form(transformation=True)
            sage: U * A == M
            True
            sage: M.is_weak_popov()
            True
            sage: U.is_invertible()
            True

            sage: PF.<x> = QQ[]
            sage: A = matrix(PF,3,[x,   x^2, x^3,\
                                   x^2, x^1, 0,\
                                   x^3, x^3, x^3])
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

    def row_reduced_form(self, transformation=None, shifts=None):
        r"""
        Return a row reduced form of this matrix.

        A matrix `M` is row reduced if the (row-wise) leading term matrix has
        the same rank as `M`. The (row-wise) leading term matrix of a polynomial
        matrix `M` is the matrix over `k` whose `(i,j)`'th entry is the
        `x^{d_i}` coefficient of `M[i,j]`, where `d_i` is the greatest degree
        among polynomials in the `i`'th row of `M_0`.

        A row reduced form is non-canonical so a given matrix has many row
        reduced forms.

        INPUT:

        - ``transformation`` -- (default: ``False``). If this ``True``, the
          transformation matrix `U` will be returned as well: this is an
          invertible matrix over `k[x]` such that ``self`` equals `UW`, where
          `W` is the output matrix.

        - ``shifts`` -- (default: ``None``) A tuple or list of integers
          `s_1, \ldots, s_n`, where `n` is the number of columns of the matrix.
          If given, a "shifted row reduced form" is computed, i.e. such that the
          matrix `A\,\mathrm{diag}(x^{s_1}, \ldots, x^{s_n})` is row reduced, where
          `\mathrm{diag}` denotes a diagonal matrix.

        OUTPUT:

        - `W` -- a row reduced form of this matrix.

        EXAMPLES::

            sage: R.<t> = GF(3)['t']
            sage: M = matrix([[(t-1)^2],[(t-1)]])
            sage: M.row_reduced_form()
            [    0]
            [t + 2]

        If the matrix is an `n \times 1` matrix with at least one non-zero entry,
        `W` has a single non-zero entry and that entry is a scalar multiple of
        the greatest-common-divisor of the entries of the matrix::

            sage: M1 = matrix([[t*(t-1)*(t+1)],[t*(t-2)*(t+2)],[t]])
            sage: output1 = M1.row_reduced_form()
            sage: output1
            [0]
            [0]
            [t]

        The following is the first half of example 5 in [Hes2002]_ *except* that we
        have transposed the matrix; [Hes2002]_ uses column operations and we use row::

            sage: R.<t> = QQ['t']
            sage: M = matrix([[t^3 - t,t^2 - 2],[0,t]]).transpose()
            sage: M.row_reduced_form()
            [      t    -t^2]
            [t^2 - 2       t]

        The next example demonstrates what happens when the matrix is a zero matrix::

            sage: R.<t> = GF(5)['t']
            sage: M = matrix(R, 2, [0,0,0,0])
            sage: M.row_reduced_form()
            [0 0]
            [0 0]

        In the following example, the original matrix is already row reduced, but
        the output is a different matrix. This is because currently this method
        simply computes a weak Popov form, which is always also a row reduced matrix
        (see :meth:`weak_popov_form`). This behavior is likely to change when a faster
        algorithm designed specifically for row reduced form is implemented in Sage::

            sage: R.<t> = QQ['t']
            sage: M = matrix([[t,t,t],[0,0,t]]); M
            [t t t]
            [0 0 t]
            sage: M.row_reduced_form()
            [ t  t  t]
            [-t -t  0]

        The last example shows the usage of the transformation parameter::

            sage: Fq.<a> = GF(2^3)
            sage: Fx.<x> = Fq[]
            sage: A = matrix(Fx,[[x^2+a,x^4+a],[x^3,a*x^4]])
            sage: W,U = A.row_reduced_form(transformation=True);
            sage: W,U
            (
            [          x^2 + a           x^4 + a]  [1 0]
            [x^3 + a*x^2 + a^2               a^2], [a 1]
            )
            sage: U*W == A
            True
            sage: U.is_invertible()
            True

        """
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
