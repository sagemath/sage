# sage.doctest: needs sage.libs.pari       (for charpoly, minimal_polynomial in __init__)
r"""
`J`-ideals of matrices

Let `B` be an `n\times n`-matrix over a principal ideal domain `D`.

For an ideal `J`, the `J`-ideal of `B` is defined to be
`N_J(B) = \{ f\in D[X] \mid f(B) \in M_n(J) \}`.

For a prime element `p` of `D` and `t\ge 0`, a `(p^t)`-minimal polynomial of `B`
is a monic polynomial `f\in N_{(p^t)}(B)` of minimal degree.

This module computes these minimal polynomials.

Let `p` be a prime element of `D`. Then there is a finite set `\mathcal{S}_p` of
positive integers and monic polynomials `\nu_{ps}` for `s\in\mathcal{S}_p` such
that for `t\ge 1`,

.. MATH::

   N_{(p^t)}(B) = \mu_BD[X] + p^tD[X]
       + \sum_{\substack{s\in\mathcal{S}_p \\ s \le b(t) }}
           p^{\max\{0,t-s\}}\nu_{ps}D[X]

holds where `b(t) = \min\{r\in \mathcal{S}_p \mid r \ge s\}`. The
degree of `\nu_{ps}` is strictly increasing in `s\in \mathcal{S}_p` and
`\nu_{ps}` is a `(p^s)`-minimal polynomial. If `t\le \max\mathcal{S}_p`,
then the summand `\mu_BD[X]` can be omitted.

All computations are done by the class
:class:`ComputeMinimalPolynomials` where various intermediate results
are cached. It provides the following methods:

* :meth:`~ComputeMinimalPolynomials.p_minimal_polynomials`
  computes `\mathcal{S}_p` and the monic polynomials `\nu_{ps}`.

* :meth:`~ComputeMinimalPolynomials.null_ideal` determines `N_{(p^t)}(B)`.

* :meth:`~ComputeMinimalPolynomials.prime_candidates` determines all primes `p`
  where `\mathcal{S}_p` might be non-empty.

* :meth:`~ComputeMinimalPolynomials.integer_valued_polynomials_generators`
  determines the generators of the ring `\{f \in K[X] \mid f(B) \in M_n(D)\}`
  of integer valued polynomials on `B`.

EXAMPLES::

    sage: from sage.matrix.compute_J_ideal import ComputeMinimalPolynomials
    sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
    sage: C = ComputeMinimalPolynomials(B)
    sage: C.prime_candidates()
    [2, 3, 5]
    sage: for t in range(4):
    ....:     print(C.null_ideal(2^t))
    Principal ideal (1) of
        Univariate Polynomial Ring in x over Integer Ring
    Ideal (2, x^2 + x) of
        Univariate Polynomial Ring in x over Integer Ring
    Ideal (4, x^2 + 3*x + 2) of
        Univariate Polynomial Ring in x over Integer Ring
    Ideal (8, x^3 + x^2 - 12*x - 20, 2*x^2 + 6*x + 4) of
        Univariate Polynomial Ring in x over Integer Ring
    sage: C.p_minimal_polynomials(2)
    {2: x^2 + 3*x + 2}
    sage: C.integer_valued_polynomials_generators()
    (x^3 + x^2 - 12*x - 20, [1, 1/4*x^2 + 3/4*x + 1/2])

The last output means that

.. MATH::

   \{f \in \QQ[X] \mid f(B) \in M_3(\ZZ)\} =
       (x^3 + x^2 - 12x - 20)\QQ[X] + \ZZ[X]
       + \frac{1}{4}(x^2 + 3x + 2) \ZZ[X].

.. TODO::

   Test code over PIDs other than ZZ.

   This requires implementation of
   :meth:`~sage.matrix.matrix_integer_dense.Matrix_integer_dense.frobenius`
   over more general domains than ZZ.

   Additionally, :func:`lifting` requires modification or a bug
   needs fixing, see
   `AskSage Question 35555 <https://ask.sagemath.org/question/35555/lifting-a-matrix-from-mathbbqyy-1/>`_.

REFERENCES:

   [Ris2016]_, [HR2016]_

AUTHORS:

- Clemens Heuberger (2016)
- Roswitha Rissner (2016)

ACKNOWLEDGEMENT:

- Clemens Heuberger is supported by the Austrian Science Fund (FWF):
  P 24644-N26.

- Roswitha Rissner is supported by the Austrian Science Fund (FWF):
  P 27816-N26.


Classes and Methods
===================
"""

# *****************************************************************************
# Copyright (C) 2016 Clemens Heuberger <clemens.heuberger@aau.at>
#               2016 Roswitha Rissner <roswitha.rissner@tugraz.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# https://www.gnu.org/licenses/
# *****************************************************************************

from sage.matrix.constructor import matrix
from sage.structure.sage_object import SageObject


def lifting(p, t, A, G):
    r"""
    Compute generators of `\{f \in D[X]^d \mid Af \equiv 0 \pmod{p^{t}}\}` given
    generators of `\{f\in D[X]^d \mid Af \equiv 0\pmod{p^{t-1}}\}`.

    INPUT:

    - ``p`` -- a prime element of some principal ideal domain `D`

    - ``t`` -- nonnegative integer

    - ``A`` -- a `c\times d` matrix over `D[X]`

    - ``G`` -- a matrix over `D[X]`. The columns of
      `\begin{pmatrix}p^{t-1}I& G\end{pmatrix}` are generators
      of `\{ f\in D[X]^d \mid Af \equiv 0\pmod{p^{t-1}}\}`;
      can be set to ``None`` if ``t`` is zero

    OUTPUT:

    A matrix `F` over `D[X]` such that the columns of
    `\begin{pmatrix}p^tI&F&pG\end{pmatrix}` are generators of
    `\{ f\in D[X]^d \mid Af \equiv 0\pmod{p^t}\}`.

    EXAMPLES::

        sage: from sage.matrix.compute_J_ideal import lifting
        sage: X = polygen(ZZ, 'X')
        sage: A = matrix([[1, X], [2*X, X^2]])
        sage: G0 = lifting(5, 0, A, None)
        sage: G1 = lifting(5, 1, A, G0); G1
        []
        sage: (A*G1 % 5).is_zero()
        True
        sage: A = matrix([[1, X, X^2], [2*X, X^2, 3*X^3]])
        sage: G0 = lifting(5, 0, A, None)
        sage: G1 = lifting(5, 1, A, G0); G1
        [3*X^2]
        [    X]
        [    1]
        sage: (A*G1 % 5).is_zero()
        True
        sage: G2 = lifting(5, 2, A, G1); G2
        [15*X^2 23*X^2]
        [   5*X      X]
        [     5      1]
        sage: (A*G2 % 25).is_zero()
        True
        sage: lifting(5, 10, A, G1)
        Traceback (most recent call last):
        ...
        ValueError: A*G not zero mod 5^9

    ALGORITHM:

    [HR2016]_, Algorithm 1.

    TESTS::

        sage: A = matrix([[1, X], [X, X^2]])
        sage: G0 = lifting(5, 0, A, None)
        sage: G1 = lifting(5, 1, A, G0); G1
        Traceback (most recent call last):
        ...
        ValueError: [  1   X|]
        [  X X^2|] does not have full rank.
    """
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

    DX = A.parent().base()
    (X,) = DX.variable_names()
    D = DX.base_ring()
    d = A.ncols()
    c = A.nrows()

    if t == 0:
        return matrix(DX, d, 0)

    if not (A*G % p**(t-1)).is_zero():
        raise ValueError("A*G not zero mod %s^%s" % (p, t-1))

    R = A*G/p**(t-1)
    R.change_ring(DX)

    AR = matrix.block([[A, R]])
    Fp = D.quotient(p*D)
    FpX = PolynomialRing(Fp, name=X)

    ARb = AR.change_ring(FpX)
    (Db, Sb, Tb) = ARb.smith_form()
    #assert Sb * ARb * Tb == Db
    #assert all(i == j or Db[i, j].is_zero()
    #           for i in range(Db.nrows())
    #           for j in range(Db.ncols()))

    r = Db.rank()
    if r != c:
        raise ValueError("{} does not have full rank.".format(ARb))

    T = Tb.change_ring(DX)

    F1 = matrix.block([[p**(t-1) * matrix.identity(d), G]])*T
    F = F1.matrix_from_columns(range(r, F1.ncols()))
    assert (A*F % (p**t)).is_zero(), "A*F=%s" % (A*F)

    return F


def p_part(f, p):
    r"""
    Compute the `p`-part of a polynomial.

    INPUT:

    - ``f`` -- a polynomial over `D`

    - ``p`` -- a prime in `D`

    OUTPUT:

    A polynomial `g` such that `\deg g \le \deg f` and
    all nonzero coefficients of `f - p g` are not divisible by `p`.

    EXAMPLES::

        sage: from sage.matrix.compute_J_ideal import p_part
        sage: X = polygen(ZZ, 'X')
        sage: f = X^3 + 5*X + 25
        sage: g = p_part(f, 5); g
        X + 5
        sage: f - 5*g
        X^3

    TESTS:

    Return value is supposed to be a polynomial, see :issue:`22402`

        sage: g = p_part(X+1, 2)
        sage: g.parent()
        Univariate Polynomial Ring in X over Integer Ring
    """
    DX = f.parent()
    (X,) = DX.gens()
    return DX(sum(c//p * X**i for i, c in enumerate(f.list())
               if c % p == 0))


class ComputeMinimalPolynomials(SageObject):
    r"""
    Create an object for computing `(p^t)`-minimal polynomials and `J`-ideals.

    For an ideal `J` and a square matrix `B` over a principal ideal
    domain `D`, the `J`-ideal of `B` is defined to be
    `N_J(B) = \{ f\in D[X] \mid f(B) \in M_n(J) \}`.

    For a prime element `p` of `D` and `t\ge 0`, a `(p^t)`-minimal
    polynomial of `B` is a monic polynomial `f\in N_{(p^t)}(B)` of
    minimal degree.

    The characteristic polynomial of `B` is denoted by `\chi_B`; `n`
    is the size of `B`.

    INPUT:

    - ``B`` -- a square matrix over a principal ideal domain `D`

    OUTPUT:

    An object which allows to call :meth:`p_minimal_polynomials`,
    :meth:`null_ideal` and :meth:`integer_valued_polynomials_generators`.

    EXAMPLES::

        sage: from sage.matrix.compute_J_ideal import ComputeMinimalPolynomials
        sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
        sage: C = ComputeMinimalPolynomials(B)
        sage: C.prime_candidates()
        [2, 3, 5]
        sage: for t in range(4):
        ....:     print(C.null_ideal(2^t))
        Principal ideal (1) of
            Univariate Polynomial Ring in x over Integer Ring
        Ideal (2, x^2 + x) of
            Univariate Polynomial Ring in x over Integer Ring
        Ideal (4, x^2 + 3*x + 2) of
            Univariate Polynomial Ring in x over Integer Ring
        Ideal (8, x^3 + x^2 - 12*x - 20, 2*x^2 + 6*x + 4) of
            Univariate Polynomial Ring in x over Integer Ring
        sage: C.p_minimal_polynomials(2)
        {2: x^2 + 3*x + 2}
        sage: C.integer_valued_polynomials_generators()
        (x^3 + x^2 - 12*x - 20, [1, 1/4*x^2 + 3/4*x + 1/2])
    """
    def __init__(self, B):
        r"""
        Initialize the ComputeMinimalPolynomials class.

        INPUT:

        - ``B`` -- a square matrix

        TESTS::

            sage: from sage.matrix.compute_J_ideal import ComputeMinimalPolynomials
            sage: ComputeMinimalPolynomials(matrix([[1, 2]]))
            Traceback (most recent call last):
            ...
            TypeError: square matrix required
        """
        from sage.rings.polynomial.polynomial_ring import polygen

        super().__init__()
        if not B.is_square():
            raise TypeError("square matrix required")

        self._B = B
        self._D = B.base_ring()
        X = polygen(self._D)
        adjugate = (X - B).adjugate()
        d = B.nrows()**2
        b = matrix(d, 1, adjugate.list())
        self.chi_B = B.charpoly(X)
        self.mu_B = B.minimal_polynomial()
        self._A = matrix.block([[b , -self.chi_B*matrix.identity(d)]])
        self._DX = X.parent()
        self._cache = {}

    def find_monic_replacements(self, p, t, pt_generators, prev_nu):
        r"""
        Replace possibly non-monic generators of `N_{(p^t)}(B)` by monic
        generators.

        INPUT:

        - ``p`` -- a prime element of `D`

        - ``t`` -- nonnegative integer

        - ``pt_generators`` -- list `(g_1, \ldots, g_s)` of polynomials in
          `D[X]` such that `N_{(p^t)}(B) = (g_1, \ldots, g_s) + pN_{(p^{t-1})}(B)`

        - ``prev_nu`` -- a `(p^{t-1})`-minimal polynomial of `B`

        OUTPUT:

        A list `(h_1, \ldots, h_r)` of monic polynomials such that
        `N_{(p^t)}(B) = (h_1, \ldots, h_r) + pN_{(p^{t-1})}(B)`.

        EXAMPLES::

            sage: from sage.matrix.compute_J_ideal import ComputeMinimalPolynomials
            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: C = ComputeMinimalPolynomials(B)
            sage: x = polygen(ZZ, 'x')
            sage: nu_1 = x^2 + x
            sage: generators_4 = [2*x^2 + 2*x, x^2 + 3*x + 2]
            sage: C.find_monic_replacements(2, 2, generators_4, nu_1)
            [x^2 + 3*x + 2]

        TESTS::

            sage: C.find_monic_replacements(2, 3, generators_4, nu_1)
            Traceback (most recent call last):
            ...
            ValueError: [2*x^2 + 2*x, x^2 + 3*x + 2] not in N_{(2^3)}(B)
            sage: C.find_monic_replacements(2, 2, generators_4, x^2)
            Traceback (most recent call last):
            ...
            ValueError: x^2 not in N_{(2^1)}(B)

        ALGORITHM:

        [HR2016]_, Algorithms 2 and 3.
        """
        from sage.arith.misc import xgcd

        if not all((g(self._B) % p**t).is_zero()
                   for g in pt_generators):
            raise ValueError("%s not in N_{(%s^%s)}(B)" %
                             (pt_generators, p, t))

        if not (prev_nu(self._B) % p**(t-1)).is_zero():
            raise ValueError("%s not in N_{(%s^%s)}(B)" % (prev_nu, p, t-1))

        (X,) = self._DX.gens()

        replacements = []
        for f in pt_generators:
            g = f
            p_prt = p_part(g, p)

            while g != p*p_prt:
                r = p_prt.quo_rem(prev_nu)[1]
                g1 = g - p*p_prt
                d, u, v = xgcd(g1.leading_coefficient(), p)
                h = u*(p*r + g1) + v*p*prev_nu*X**(g1.degree()-prev_nu.degree())
                replacements.append(h % p**t)
                #reduce coefficients mod p^t to keep coefficients small
                g = g.quo_rem(h)[1]
                p_prt = p_part(g, p)

        replacements = list(set(replacements))
        assert all(g.is_monic() for g in replacements),\
            "Something went wrong in find_monic_replacements"

        return replacements

    def current_nu(self, p, t, pt_generators, prev_nu):
        r"""
        Compute `(p^t)`-minimal polynomial of `B`.

        INPUT:

        - ``p`` -- a prime element of `D`

        - ``t`` -- positive integer

        - ``pt_generators`` -- list `(g_1, \ldots, g_s)` of polynomials in
          `D[X]` such that `N_{(p^t)}(B) = (g_1, \ldots, g_s) + pN_{(p^{t-1})}(B)`

        - ``prev_nu`` -- a `(p^{t-1})`-minimal polynomial of `B`

        OUTPUT:

        A `(p^t)`-minimal polynomial of `B`.

        EXAMPLES::

            sage: from sage.matrix.compute_J_ideal import ComputeMinimalPolynomials
            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: C = ComputeMinimalPolynomials(B)
            sage: x = polygen(ZZ, 'x')
            sage: nu_1 = x^2 + x
            sage: generators_4 = [2*x^2 + 2*x, x^2 + 3*x + 2]
            sage: C.current_nu(2, 2, generators_4, nu_1)
            x^2 + 3*x + 2

        ALGORITHM:

        [HR2016]_, Algorithm 4.

        TESTS::

            sage: C.current_nu(2, 3, generators_4, nu_1)
            Traceback (most recent call last):
            ...
            ValueError: [2*x^2 + 2*x, x^2 + 3*x + 2] not in N_{(2^3)}(B)
            sage: C.current_nu(2, 2, generators_4, x^2)
            Traceback (most recent call last):
            ...
            ValueError: x^2 not in N_{(2^1)}(B)
        """
        import heapq

        from sage.misc.verbose import verbose

        if not all((g(self._B) % p**t).is_zero()
                   for g in pt_generators):
            raise ValueError("%s not in N_{(%s^%s)}(B)" %
                             (pt_generators, p, t))

        if not (prev_nu(self._B) % p**(t-1)).is_zero():
            raise ValueError("%s not in N_{(%s^%s)}(B)" % (prev_nu, p, t-1))

        generators = self.find_monic_replacements(p, t, pt_generators, prev_nu)

        verbose("------------------------------------------")
        verbose(pt_generators)
        verbose("Generators with (p^t)-generating property:")
        verbose(generators)

        heap = list((f.degree(), f) for f in generators)
        heapq.heapify(heap)

        # find poly of minimal degree
        deg_g, g = heapq.heappop(heap)

        # find nu
        while heap:
            deg_f, f = heapq.heappop(heap)
            #take first element in generators not equal g
            r = (f.quo_rem(g)[1]) % p**t
            if r != 0:
                for h in self.find_monic_replacements(p, t, [r], prev_nu):
                    heapq.heappush(heap, (h.degree(), h))
                if heap and heap[0][0] < deg_g:
                    deg_g, g = heapq.heappushpop(heap, (deg_g, g))

            verbose([g] + [h for (deg_h, h) in heap])

        return g

    def mccoy_column(self, p, t, nu):
        r"""
        Compute matrix for McCoy's criterion.

        INPUT:

        - ``p`` -- a prime element in `D`

        - ``t`` -- positive integer

        - ``nu`` -- a `(p^t)`-minimal polynomial of `B`

        OUTPUT:

        An `(n^2 + 1) \times 1` matrix `g` with first entry ``nu`` such that
        `\begin{pmatrix}b& -\chi_B I\end{pmatrix}g \equiv 0\pmod{p^t}` where `b`
        consists of the entries of `\operatorname{adj}(X-B)`.

        EXAMPLES::

            sage: from sage.matrix.compute_J_ideal import ComputeMinimalPolynomials
            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: C = ComputeMinimalPolynomials(B)
            sage: x = polygen(ZZ, 'x')
            sage: nu_4 = x^2 + 3*x + 2
            sage: g = C.mccoy_column(2, 2, nu_4)
            sage: b = matrix(9, 1, (x - B).adjugate().list())
            sage: M = matrix.block([[b, -B.charpoly(x)*matrix.identity(9)]])
            sage: (M*g % 4).is_zero()
            True

        ALGORITHM:

        [HR2016]_, Algorithm 5.

        TESTS::

            sage: nu_2 = x^2 + x
            sage: C.mccoy_column(2, 2, nu_2)
            Traceback (most recent call last):
            ...
            ValueError: x^2 + x not in (2^2)-ideal
        """
        if not (nu(self._B) % p**t).is_zero():
            raise ValueError(
                "%s not in (%s^%s)-ideal" % (nu, p, t))

        column = matrix(self._DX, self._A.ncols(), 1,
                        [nu] + [(nu*b).quo_rem(self.chi_B)[0]
                                  for b in self._A[:, 0].list()])

        assert (self._A * column % p**t).is_zero(),\
                                 "McCoy column incorrect"

        return column

    def p_minimal_polynomials(self, p, s_max=None):
        r"""
        Compute `(p^s)`-minimal polynomials `\nu_s` of `B`.

        Compute a finite subset `\mathcal{S}` of the positive
        integers and `(p^s)`-minimal polynomials
        `\nu_s` for `s \in \mathcal{S}`.

        For `0 < t \le \max \mathcal{S}`, a `(p^t)`-minimal polynomial is
        given by `\nu_s` where
        `s = \min\{ r \in \mathcal{S} \mid r\ge t \}`.
        For `t > \max \mathcal{S}`, the minimal polynomial of `B` is
        also a `(p^t)`-minimal polynomial.

        INPUT:

        - ``p`` -- a prime in `D`

        - ``s_max`` -- positive integer (default: ``None``); if set, only
          `(p^s)`-minimal polynomials for ``s <= s_max`` are computed
          (see below for details)

        OUTPUT:

        A dictionary. Keys are the finite set `\mathcal{S}`, the values
        are the associated `(p^s)`-minimal polynomials `\nu_s`,
        `s \in \mathcal{S}`.

        Setting ``s_max`` only affects the output if ``s_max`` is at
        most `\max\mathcal{S}` where `\mathcal{S}` denotes the full
        set. In that case, only those `\nu_s` with ``s <= s_max`` are
        returned where ``s_max`` is always included even if it is not
        included in the full set `\mathcal{S}`.

        EXAMPLES::

            sage: from sage.matrix.compute_J_ideal import ComputeMinimalPolynomials
            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: C = ComputeMinimalPolynomials(B)
            sage: C.p_minimal_polynomials(2)
            {2: x^2 + 3*x + 2}
            sage: set_verbose(1)
            sage: C = ComputeMinimalPolynomials(B)
            sage: C.p_minimal_polynomials(2)
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            ------------------------------------------
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            p = 2, t = 1:
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            Result of lifting:
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            F =
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            [x^2 + x]
            [      x]
            [      0]
            [      1]
            [      1]
            [  x + 1]
            [      1]
            [      0]
            [      0]
            [  x + 1]
            verbose 1 (...: compute_J_ideal.py, current_nu)
            ------------------------------------------
            verbose 1 (...: compute_J_ideal.py, current_nu)
            (x^2 + x)
            verbose 1 (...: compute_J_ideal.py, current_nu)
            Generators with (p^t)-generating property:
            verbose 1 (...: compute_J_ideal.py, current_nu)
            [x^2 + x]
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            nu = x^2 + x
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            corresponding columns for G
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            [x^2 + x]
            [  x + 2]
            [      0]
            [      1]
            [      1]
            [  x - 1]
            [     -1]
            [     10]
            [      0]
            [  x + 1]
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            ------------------------------------------
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            p = 2, t = 2:
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            Result of lifting:
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            F =
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            [  2*x^2 + 2*x x^2 + 3*x + 2]
            [          2*x         x + 4]
            [            0             0]
            [            2             1]
            [            2             1]
            [      2*x + 2         x + 1]
            [            2            -1]
            [            0            10]
            [            0             0]
            [      2*x + 2         x + 3]
            verbose 1 (...: compute_J_ideal.py, current_nu)
            ------------------------------------------
            verbose 1 (...: compute_J_ideal.py, current_nu)
            (2*x^2 + 2*x, x^2 + 3*x + 2)
            verbose 1 (...: compute_J_ideal.py, current_nu)
            Generators with (p^t)-generating property:
            verbose 1 (...: compute_J_ideal.py, current_nu)
            [x^2 + 3*x + 2]
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            nu = x^2 + 3*x + 2
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            corresponding columns for G
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            [x^2 + 3*x + 2]
            [        x + 4]
            [            0]
            [            1]
            [            1]
            [        x + 1]
            [           -1]
            [           10]
            [            0]
            [        x + 3]
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            ------------------------------------------
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            p = 2, t = 3:
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            Result of lifting:
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            F =
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            [x^3 + 7*x^2 + 6*x x^3 + 3*x^2 + 2*x]
            [        x^2 + 8*x         x^2 + 4*x]
            [                0                 0]
            [                x             x + 4]
            [            x + 4                 x]
            [    x^2 + 5*x + 4           x^2 + x]
            [           -x + 4                -x]
            [             10*x              10*x]
            [                0                 0]
            [        x^2 + 7*x     x^2 + 3*x + 4]
            verbose 1 (...: compute_J_ideal.py, current_nu)
            ------------------------------------------
            verbose 1 (...: compute_J_ideal.py, current_nu)
            (x^3 + 7*x^2 + 6*x, x^3 + 3*x^2 + 2*x)
            verbose 1 (...: compute_J_ideal.py, current_nu)
            Generators with (p^t)-generating property:
            verbose 1 (...: compute_J_ideal.py, current_nu)
            ...
            verbose 1 (...: compute_J_ideal.py, current_nu)
            [x^3 + 3*x^2 + 2*x]
            verbose 1 (...: compute_J_ideal.py, p_minimal_polynomials)
            nu = x^3 + 3*x^2 + 2*x
            {2: x^2 + 3*x + 2}
            sage: set_verbose(0)
            sage: C.p_minimal_polynomials(2, s_max=1)
            {1: x^2 + x}
            sage: C.p_minimal_polynomials(2, s_max=2)
            {2: x^2 + 3*x + 2}
            sage: C.p_minimal_polynomials(2, s_max=3)
            {2: x^2 + 3*x + 2}

        ALGORITHM:

        [HR2016]_, Algorithm 5.
        """
        from sage.misc.verbose import verbose
        from sage.rings.infinity import Infinity

        deg_mu = self.mu_B.degree()
        if s_max is None:
            s_max = Infinity

        if p in self._cache:
            (t, G, p_min_polys) = self._cache[p]
            if t < Infinity:
                nu = G[0][0]
        else:
            t = 0
            p_min_polys = {}
            nu = self._DX(1)
            d = self._A.ncols()
            G = matrix(self._DX, d, 0)

        while t < s_max:
            deg_prev_nu = nu.degree()
            t += 1
            verbose("------------------------------------------")
            verbose("p = %s, t = %s:" % (p, t))

            verbose("Result of lifting:")
            verbose("F =")
            F = lifting(p, t, self._A, G)
            verbose(F)

            nu = self.current_nu(p, t, F[0], nu)

            verbose("nu = %s" % nu)
            if nu.degree() >= deg_mu:
                t = Infinity
                break

            if nu.degree() == deg_prev_nu:
                G = G.delete_columns([G.ncols() - 1])
                del p_min_polys[t-1]

            column = self.mccoy_column(p, t, nu)
            verbose("corresponding columns for G")
            verbose(column)

            G = matrix.block([[p * G, column]])
            p_min_polys[t] = nu

        self._cache[p] = (t, G, p_min_polys)

        if s_max < t:
            result = {r: polynomial
                      for r, polynomial in p_min_polys.items() if r < s_max}
            next_t_candidates = list(r for r in p_min_polys if r >= s_max)
            if next_t_candidates:
                next_t = min(next_t_candidates)
                result.update({s_max: p_min_polys[next_t] % p**s_max})

            return result

        return p_min_polys

    def null_ideal(self, b=0):
        r"""
        Return the `(b)`-ideal `N_{(b)}(B)=\{f\in D[X] \mid f(B)\in M_n(bD)\}`.

        INPUT:

        - ``b`` -- an element of `D` (default: 0)

        OUTPUT: an ideal in `D[X]`

        EXAMPLES::

            sage: from sage.matrix.compute_J_ideal import ComputeMinimalPolynomials
            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: C = ComputeMinimalPolynomials(B)
            sage: C.null_ideal()
            Principal ideal (x^3 + x^2 - 12*x - 20) of
                Univariate Polynomial Ring in x over Integer Ring
            sage: C.null_ideal(2)
            Ideal (2, x^2 + x) of
                Univariate Polynomial Ring in x over Integer Ring
            sage: C.null_ideal(4)
            Ideal (4, x^2 + 3*x + 2) of
                Univariate Polynomial Ring in x over Integer Ring
            sage: C.null_ideal(8)
            Ideal (8, x^3 + x^2 - 12*x - 20, 2*x^2 + 6*x + 4) of
                Univariate Polynomial Ring in x over Integer Ring
            sage: C.null_ideal(3)
            Ideal (3, x^3 + x^2 - 12*x - 20) of
                Univariate Polynomial Ring in x over Integer Ring
            sage: C.null_ideal(6)
            Ideal (6, 2*x^3 + 2*x^2 - 24*x - 40, 3*x^2 + 3*x) of
                Univariate Polynomial Ring in x over Integer Ring
        """
        from sage.arith.misc import factor

        mu_B_coefficients = []
        generators = []

        if b == 0:
            mu_B_coefficients = [1]
        else:
            for (p, t) in factor(b):
                cofactor = b // p**t
                p_polynomials = self.p_minimal_polynomials(p, t)
                generators += [cofactor*p**(t-s)*nu
                               for s, nu in p_polynomials.items()]
                if not p_polynomials or max(p_polynomials) < t:
                    mu_B_coefficients.append(cofactor)

            assert all((g(self._B) % b).is_zero() for g in generators), \
                "Polynomials not in %s-ideal" % (b,)

        if mu_B_coefficients:
            (mu_B_coefficient,) = self._D.ideal(mu_B_coefficients).gens()
            generators = [mu_B_coefficient * self.mu_B] + generators

        if b != 0:
            generators = [self._DX(b)] + generators

        return self._DX.ideal(generators)

    def prime_candidates(self):
        r"""
        Determine those primes `p` where `\mu_B` might not be a
        `(p)`-minimal polynomial.

        OUTPUT: list of primes

        EXAMPLES::

            sage: from sage.matrix.compute_J_ideal import ComputeMinimalPolynomials
            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: C = ComputeMinimalPolynomials(B)
            sage: C.prime_candidates()
            [2, 3, 5]
            sage: C.p_minimal_polynomials(2)
            {2: x^2 + 3*x + 2}
            sage: C.p_minimal_polynomials(3)
            {}
            sage: C.p_minimal_polynomials(5)
            {}

        This means that `3` and `5` were candidates, but actually, `\mu_B` turns
        out to be a `(3)`-minimal polynomial and a `(5)`-minimal polynomial.
        """
        from sage.arith.misc import factor

        _, T = self._B.frobenius_form(2)

        return [p for p, _ in factor(T.det())]

    def integer_valued_polynomials_generators(self):
        r"""
        Determine the generators of the ring of integer valued polynomials on `B`.

        OUTPUT:

        A pair ``(mu_B, P)`` where ``P`` is a list of polynomials in `K[X]`
        such that

        .. MATH::

           \{f \in K[X] \mid f(B) \in M_n(D)\} = \mu_B K[X]
               + \sum_{g\in P} g D[X]

        where `K` denotes the fraction field of `D`.

        EXAMPLES::

            sage: from sage.matrix.compute_J_ideal import ComputeMinimalPolynomials
            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: C = ComputeMinimalPolynomials(B)
            sage: C.integer_valued_polynomials_generators()
            (x^3 + x^2 - 12*x - 20, [1, 1/4*x^2 + 3/4*x + 1/2])
        """
        return (self.mu_B, [self._DX(1)] +
                [nu / p**s
                 for p in self.prime_candidates()
                 for s, nu in self.p_minimal_polynomials(p).items()])
