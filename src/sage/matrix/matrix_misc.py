# sage_setup: distribution = sagemath-modules
"""
Miscellaneous matrix functions
"""

# ****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.fields import Fields
_Fields = Fields()


def row_iterator(A):
    for i in range(A.nrows()):
        yield A.row(i)


def prm_mul(p1, p2, mask_free, prec):
    """
    Return the product of ``p1`` and ``p2``, putting free variables in
    ``mask_free`` to `1`.

    This function is mainly use as a subroutine of
    :func:`permanental_minor_polynomial`.

    INPUT:

    - `p1,p2` -- polynomials as dictionaries

    - ``mask_free`` -- integer mask that give the list of free variables
      (the `i`-th variable is free if the `i`-th bit of ``mask_free`` is `1`)

    - ``prec`` -- if ``prec`` is not ``None``, truncate the product at precision ``prec``

    EXAMPLES::

        sage: from sage.matrix.matrix_misc import prm_mul
        sage: t = polygen(ZZ, 't')
        sage: p1 = {0: 1, 1: t, 4: t}
        sage: p2 = {0: 1, 1: t, 2: t}
        sage: prm_mul(p1, p2, 1, None)
        {0: 2*t + 1, 2: t^2 + t, 4: t^2 + t, 6: t^2}
    """
    p = {}
    if not p2:
        return p
    for exp1, v1 in p1.items():
        if v1.is_zero():
            continue
        for exp2, v2 in p2.items():
            if exp1 & exp2:
                continue
            v = v1 * v2
            if prec is not None:
                v._unsafe_mutate(prec, 0)
            exp = exp1 | exp2
            exp = exp ^ (exp & mask_free)
            if exp not in p:
                p[exp] = v
            else:
                p[exp] += v
    return p


def permanental_minor_polynomial(A, permanent_only=False, var='t', prec=None):
    r"""
    Return the polynomial of the sums of permanental minors of ``A``.

    INPUT:

    - ``A`` -- a matrix

    - ``permanent_only`` -- if ``True``, return only the permanent of `A`

    - ``var`` -- name of the polynomial variable

    - ``prec`` -- if prec is not None, truncate the polynomial at precision `prec`


    The polynomial of the sums of permanental minors is

    .. MATH::

        \sum_{i=0}^{min(nrows, ncols)} p_i(A) x^i

    where `p_i(A)` is the `i`-th permanental minor of `A` (that can also be
    obtained through the method
    :meth:`~sage.matrix.matrix2.Matrix.permanental_minor` via
    ``A.permanental_minor(i)``).

    The algorithm implemented by that function has been developed by P. Butera
    and M. Pernici, see [BP2015]_. Its complexity is `O(2^n m^2 n)` where `m` and
    `n` are the number of rows and columns of `A`.  Moreover, if `A` is a banded
    matrix with width `w`, that is `A_{ij}=0` for `|i - j| > w` and `w < n/2`,
    then the complexity of the algorithm is `O(4^w (w+1) n^2)`.

    INPUT:

    - ``A`` -- matrix

    - ``permanent_only`` -- boolean (default: ``False``); if ``True``, only the
      permanent is computed (might be faster)

    - ``var`` -- a variable name

    EXAMPLES::

        sage: from sage.matrix.matrix_misc import permanental_minor_polynomial
        sage: m = matrix([[1,1],[1,2]])
        sage: permanental_minor_polynomial(m)
        3*t^2 + 5*t + 1
        sage: permanental_minor_polynomial(m, permanent_only=True)
        3
        sage: permanental_minor_polynomial(m, prec=2)
        5*t + 1

    ::

        sage: M = MatrixSpace(ZZ,4,4)
        sage: A = M([1,0,1,0,1,0,1,0,1,0,10,10,1,0,1,1])
        sage: permanental_minor_polynomial(A)
        84*t^3 + 114*t^2 + 28*t + 1
        sage: [A.permanental_minor(i) for i in range(5)]
        [1, 28, 114, 84, 0]

    An example over `\QQ`::

        sage: M = MatrixSpace(QQ,2,2)
        sage: A = M([1/5,2/7,3/2,4/5])
        sage: permanental_minor_polynomial(A, True)
        103/175

    An example with polynomial coefficients::

        sage: R.<a> = PolynomialRing(ZZ)
        sage: A = MatrixSpace(R,2)([[a,1], [a,a+1]])
        sage: permanental_minor_polynomial(A, True)
        a^2 + 2*a

    A usage of the ``var`` argument::

        sage: m = matrix(ZZ,4,[0,1,2,3,1,2,3,0,2,3,0,1,3,0,1,2])
        sage: permanental_minor_polynomial(m, var='x')
        164*x^4 + 384*x^3 + 172*x^2 + 24*x + 1

    ALGORITHM:

        The permanent `perm(A)` of a `n \times n` matrix `A` is the coefficient
        of the `x_1 x_2 \ldots x_n` monomial in

        .. MATH::

            \prod_{i=1}^n \left( \sum_{j=1}^n A_{ij} x_j \right)

        Evaluating this product one can neglect `x_i^2`, that is `x_i`
        can be considered to be nilpotent of order `2`.

        To formalize this procedure, consider the algebra
        `R = K[\eta_1, \eta_2, \ldots, \eta_n]` where the `\eta_i` are
        commuting, nilpotent of order `2` (i.e. `\eta_i^2 = 0`).
        Formally it is the quotient ring of the polynomial
        ring in `\eta_1, \eta_2, \ldots, \eta_n` quotiented by the ideal
        generated by the `\eta_i^2`.

        We will mostly consider the ring `R[t]` of polynomials over `R`. We
        denote a generic element of `R[t]` by `p(\eta_1, \ldots, \eta_n)` or
        `p(\eta_{i_1}, \ldots, \eta_{i_k})` if we want to emphasize that some
        monomials in the `\eta_i` are missing.

        Introduce an "integration" operation `\langle p \rangle` over `R` and
        `R[t]` consisting in the sum of the coefficients of the non-vanishing
        monomials in `\eta_i` (i.e. the result of setting all variables `\eta_i`
        to `1`). Let us emphasize that this is *not* a morphism of algebras as
        `\langle \eta_1 \rangle^2 = 1` while `\langle \eta_1^2 \rangle = 0`!

        Let us consider an example of computation.
        Let `p_1 = 1 + t \eta_1 + t \eta_2` and
        `p_2 = 1 + t \eta_1 + t \eta_3`. Then

        .. MATH::

            p_1 p_2 = 1 + 2t \eta_1 +
                    t (\eta_2 + \eta_3) +
                    t^2 (\eta_1 \eta_2 + \eta_1 \eta_3 + \eta_2 \eta_3)

        and

        .. MATH::

            \langle p_1 p_2 \rangle = 1 + 4t + 3t^2

        In this formalism, the permanent is just

        .. MATH::

            perm(A) = \langle \prod_{i=1}^n \sum_{j=1}^n A_{ij} \eta_j \rangle

        A useful property of `\langle . \rangle` which makes this algorithm
        efficient for band matrices is the following: let
        `p_1(\eta_1, \ldots, \eta_n)` and `p_2(\eta_j, \ldots, \eta_n)` be
        polynomials in `R[t]` where `j \ge 1`.  Then one has

        .. MATH::

            \langle p_1(\eta_1, \ldots, \eta_n) p_2 \rangle =
            \langle p_1(1, \ldots, 1, \eta_j, \ldots, \eta_n) p_2 \rangle

        where `\eta_1,..,\eta_{j-1}` are replaced by `1` in `p_1`. Informally,
        we can "integrate" these variables *before* performing the product. More
        generally, if a monomial `\eta_i` is missing in one of the terms of a
        product of two terms, then it can be integrated in the other term.

        Now let us consider an `m \times n` matrix with `m \leq n`. The *sum of
        permanental `k`-minors of `A`* is

        .. MATH::

            perm(A, k) = \sum_{r,c} perm(A_{r,c})

        where the sum is over the `k`-subsets `r` of rows and `k`-subsets `c` of
        columns and `A_{r,c}` is the submatrix obtained from `A` by keeping only
        the rows `r` and columns `c`. Of course
        `perm(A, \min(m,n)) = perm(A)` and note that `perm(A,1)` is just the sum
        of all entries of the matrix.

        The generating function of these sums of permanental minors is

        .. MATH::

            g(t) = \left\langle
            \prod_{i=1}^m \left(1 + t \sum_{j=1}^n A_{ij} \eta_j\right)
            \right\rangle

        In fact the `t^k` coefficient of `g(t)` corresponds to choosing
        `k` rows of `A`;  `\eta_i` is associated to the `i`-th column;
        nilpotency avoids having twice the same column in a product of `A`'s.

        For more details, see the article [BP2015]_.

        From a technical point of view, the product in
        `K[\eta_1, \ldots, \eta_n][t]` is implemented as a subroutine in
        :func:`prm_mul`. The indices of the rows and columns actually start at
        `0`, so the variables are  `\eta_0, \ldots, \eta_{n-1}`. Polynomials are
        represented in dictionary form: to a variable `\eta_i` is associated
        the key `2^i` (or in Python ``1 << i``). The keys associated to products
        are obtained by considering the development in base `2`: to the monomial
        `\eta_{i_1} \ldots \eta_{i_k}` is associated the key
        `2^{i_1} + \ldots + 2^{i_k}`. So the product `\eta_1 \eta_2` corresponds
        to the key `6 = (110)_2` while `\eta_0 \eta_3` has key `9 = (1001)_2`.
        In particular all operations on monomials are implemented via bitwise
        operations on the keys.
    """
    if permanent_only:
        prec = None
    elif prec is not None:
        prec = int(prec)
        if prec == 0:
            raise ValueError('the argument `prec` must be a positive integer')

    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

    K = PolynomialRing(A.base_ring(), var)
    nrows = A.nrows()
    ncols = A.ncols()
    A = A.rows()
    p = {0: K.one()}
    t = K.gen()
    vars_to_do = list(range(ncols))
    for i in range(nrows):
        # build the polynomial p1 = 1 + t sum A_{ij} eta_j
        if permanent_only:
            p1 = {}
        else:
            p1 = {0: K.one()}
        a = A[i]   # the i-th row of A
        for j in range(len(a)):
            if a[j]:
                p1[1 << j] = a[j] * t

        # make the product with the preceding polynomials, taking care of
        # variables that can be integrated
        mask_free = 0
        j = 0
        while j < len(vars_to_do):
            jj = vars_to_do[j]
            if all(A[k][jj] == 0 for k in range(i+1, nrows)):
                mask_free += 1 << jj
                vars_to_do.remove(jj)
            else:
                j += 1
        p = prm_mul(p, p1, mask_free, prec)

    if not p:
        return K.zero()

    if len(p) != 1 or 0 not in p:
        raise RuntimeError("Something is wrong! Certainly a problem in the"
                           " algorithm... please contact sage-devel@googlegroups.com")

    p = p[0]
    return p[min(nrows,ncols)] if permanent_only else p
