from sage.matrix.constructor import identity_matrix
from sage.rings.polynomial.laurent_polynomial_ring_base import LaurentPolynomialRing_generic


def laurent_matrix_reduction(A):
    """
    From a matrix `A` of Laurent polynomials, apply elementary operations
    to obtain a matrix `P` of polynomials such that the variables do not divide
    no column and no row.

    OUTPUT: Three matrices `L`, `P`, `R` such that `A = L P R`, where `L` and
            `R` are diagonal with monomial entries.

    EXAMPLES:

        sage: from sage.rings.polynomial.laurent_reduction import laurent_matrix_reduction
        sage: R.<x, y> = LaurentPolynomialRing(QQ)
        sage: L = [1/3*x^-1*y - 6*x^-2*y^2 - 1/2*x^-2*y, 1/5*x + 1/2*y + 1/6]
        sage: L += [1/2 - 5*x^-1*y - 2*x^-1, -1/3*y^-2 - 4*x^-1*y^-1 + 11*x^-1*y^-2]
        sage: A = matrix(R, 2, L)
        sage: lf, P, rg = laurent_matrix_reduction(A)
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
    R = A.base_ring()
    n_rows, n_cols = A.dimensions()
    mat_l = identity_matrix(R, n_rows)
    mat_r = identity_matrix(R, n_cols)
    if not isinstance(R, LaurentPolynomialRing_generic):
        return mat_l, A, mat_r
    res = A
    for j, rw in enumerate(A.rows()):
        for t in R.gens():
            n = min(mon.degree(t) for a in rw for cf, mon in a)
            res.rescale_row(j, t ** -n)
            mat_l.rescale_col(j, t ** n)
    for j, cl in enumerate(A.columns()):
        for t in R.gens():
            n = min(mon.degree(t) for a in cl for cf, mon in a)
            res.rescale_col(j, t ** -n)
            mat_r.rescale_row(j, t ** n)
    res = res.change_ring(R.polynomial_ring())
    return mat_l, res, mat_r
