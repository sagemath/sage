# sage_setup: distribution = sagemath-modules
"""
Dense matrices over multivariate polynomials over fields.

AUTHOR:

- Enrique Artal (2023-??): initial version
"""

# *****************************************************************************
#       Copyright (C) 2023 Enrique Artal <artal@unizar.es>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from sage.matrix.constructor import identity_matrix


cdef class Matrix_laurent_mpolynomial_dense(Matrix_generic_dense):
    """
    Dense matrix over a Laurent multivariate polynomial ring over a field.
    """
    def laurent_matrix_reduction(self):
        """
        From a matrix ``self`` of Laurent polynomials, apply elementary operations
        to obtain a matrix ``P`` of polynomials such that the variables do not divide
        any column and any row.

        OUTPUT:

        Three matrices ``L``, ``P``, ``R`` such that ``self`` equals ``L P R``,
        where ``L`` and ``R`` are diagonal with monomial entries.

        EXAMPLES:

            sage: R.<x, y> = LaurentPolynomialRing(QQ)
            sage: L = [1/3*x^-1*y - 6*x^-2*y^2 - 1/2*x^-2*y, 1/5*x + 1/2*y + 1/6]
            sage: L += [1/2 - 5*x^-1*y - 2*x^-1, -1/3*y^-2 - 4*x^-1*y^-1 + 11*x^-1*y^-2]
            sage: A = matrix(R, 2, L)
            sage: lf, P, rg = A.laurent_matrix_reduction()
            sage: lf
            [     x^-2         0]
            [        0 x^-1*y^-2]
            sage: P
            [            1/3*x - 6*y - 1/2 1/5*x^3 + 1/2*x^2*y + 1/6*x^2]
            [        1/2*x*y - 5*y^2 - 2*y             -1/3*x - 4*y + 11]
            sage: rg
            [y 0]
            [0 1]
        """
        R = self.base_ring()
        n_rows, n_cols = self.dimensions()
        mat_l = identity_matrix(R, n_rows)
        mat_r = identity_matrix(R, n_cols)
        res = self.__copy__()
        for j, rw in enumerate(res.rows()):
            for t in R.gens():
                n = min(mon.degree(t) for a in rw for _, mon in a)
                res.rescale_row(j, t ** -n)
                mat_l.rescale_col(j, t ** n)
        for j, cl in enumerate(res.columns()):
            for t in R.gens():
                n = min(mon.degree(t) for a in cl for _, mon in a)
                res.rescale_col(j, t ** -n)
                mat_r.rescale_row(j, t ** n)
        res = res.change_ring(R.polynomial_ring())
        return mat_l, res, mat_r

    def _fitting_ideal(self, i):
        r"""
        Return the `i`-th Fitting ideal of the matrix.

        This is the ideal generated
        by the `n - i` minors, where `n` is the number of columns.

        INPUT:

        - ``i`` -- integer

        OUTPUT: an ideal on the base ring

        EXAMPLES::

            sage: R.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: M = matrix(R, [[2*x^-1-z, 0, y-z^-2, 0], [0, z - y^-1, z - x, 0],[z - y, x^-2 - y, 0, z]])
            sage: M
            [-z + 2*x^-1           0    y - z^-2           0]
            [          0    z - y^-1      -x + z           0]
            [     -y + z   -y + x^-2           0           z]
            sage: M.fitting_ideal(0)
            Ideal (0) of Multivariate Laurent Polynomial Ring in x, y, z over Rational Field
            sage: M.fitting_ideal(1) == M._fitting_ideal(1)
            True
            sage: M.fitting_ideal(1).groebner_basis()
            (x^4 - 2*x^3*y - x*z^3 - 4*x^2*y + 8*x*y^2 + 4*x*y*z + 2*z^2 - 8*y,
            x*y*z^2 - x*z - 2*y*z + 2,
            x^2*z - x*z^2 - 2*x + 2*z,
            y^2*z + 1/4*x^2 - 1/2*x*y - 1/4*x*z - y + 1/2)
            sage: M.fitting_ideal(2).groebner_basis()
            (1,)
            sage: M.fitting_ideal(3).groebner_basis()
            (1,)
            sage: M.fitting_ideal(4).groebner_basis()
            (1,)
            sage: [R.ideal(M.minors(i)) == M._fitting_ideal(4 - i) for i in range(5)]
            [True, True, True, True, True]
        """
        R = self.base_ring()
        S = R.polynomial_ring()
        A = self.laurent_matrix_reduction()[1].change_ring(S)
        J = A._fitting_ideal(i)
        return J.change_ring(R)
