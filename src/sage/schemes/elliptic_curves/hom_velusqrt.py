r"""
√élu Algorithm for Elliptic-Curve Isogenies

The √élu algorithm computes isogenies of elliptic curves in time
`\tilde O(\sqrt\ell)` rather than naïvely `O(\ell)`, where `\ell`
is the degree.

The core idea is to reindex the points in the kernel subgroup in
a baby-step-giant-step manner, then use fast resultant computations
to evaluate "elliptic polynomials"
(see :class:`FastEllipticPolynomial`)
in essentially square-root time.

Based on experiments with Sage version 9.7,
the isogeny degree where
:class:`EllipticCurveHom_velusqrt`
begins to outperform
:class:`~sage.schemes.elliptic_curves.ell_curve_isogeny.EllipticCurveIsogeny`
can be as low as `\approx 100`,
but is typically closer to `\approx 1000`,
depending on the exact situation.

REFERENCES: [BDLS2020]_

EXAMPLES::

    sage: from sage.schemes.elliptic_curves.hom_velusqrt import EllipticCurveHom_velusqrt
    sage: E = EllipticCurve(GF(6666679), [5,5])
    sage: K = E(9970, 1003793, 1)
    sage: K.order()
    10009
    sage: phi = EllipticCurveHom_velusqrt(E, K)
    sage: phi
    Elliptic-curve isogeny (using √élu) of degree 10009:
      From: Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 6666679
      To:   Elliptic Curve defined by y^2 = x^3 + 227975*x + 3596133 over Finite Field of size 6666679
    sage: phi.codomain()
    Elliptic Curve defined by y^2 = x^3 + 227975*x + 3596133 over Finite Field of size 6666679

Note that the isogeny is usually not identical to the one computed by
:class:`~sage.schemes.elliptic_curves.ell_curve_isogeny.EllipticCurveIsogeny`::

    sage: psi = EllipticCurveIsogeny(E, K)
    sage: psi
    Isogeny of degree 10009
        from Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 6666679
        to Elliptic Curve defined by y^2 = x^3 + 5344836*x + 3950273 over Finite Field of size 6666679

However, they are certainly separable isogenies with the same kernel
and must therefore be equal *up to post-isomorphism*::

    sage: isos = psi.codomain().isomorphisms(phi.codomain())
    sage: sum(iso * psi == phi for iso in isos)
    1

Just like
:class:`~sage.schemes.elliptic_curves.ell_curve_isogeny.EllipticCurveIsogeny`,
the constructor supports a ``model`` keyword argument::

    sage: E = EllipticCurve(GF(6666679), [1,1])
    sage: K = E(9091, 517864)
    sage: phi = EllipticCurveHom_velusqrt(E, K, model='montgomery')
    sage: phi
    Elliptic-curve isogeny (using √élu) of degree 2999:
      From: Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 6666679
      To:   Elliptic Curve defined by y^2 = x^3 + 1559358*x^2 + x over Finite Field of size 6666679

Internally, :class:`EllipticCurveHom_velusqrt` works on short
Weierstraß curves, but it performs the conversion automatically::

    sage: E = EllipticCurve(GF(101), [1,2,3,4,5])
    sage: K = E(1, 2)
    sage: K.order()
    37
    sage: EllipticCurveHom_velusqrt(E, K)
    Elliptic-curve isogeny (using √élu) of degree 37:
      From: Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Finite Field of size 101
      To:   Elliptic Curve defined by y^2 = x^3 + 66*x + 86 over Finite Field of size 101

However, this does imply not all elliptic curves are supported.
Curves without a short Weierstraß model exist in characteristics
`2` and `3`::

    sage: F.<t> = GF(3^3)
    sage: E = EllipticCurve(F, [1,1,1,1,1])
    sage: P = E(t^2+2, 1)
    sage: P.order()
    19
    sage: EllipticCurveHom_velusqrt(E, P)
    Traceback (most recent call last):
    ...
    NotImplementedError: only implemented for curves having a short Weierstrass model

Furthermore, the implementation is restricted to finite fields,
since this appears to be the most relevant application for the
√élu algorithm::

    sage: E = EllipticCurve('26b1')
    sage: P = E(1,0)
    sage: P.order()
    7
    sage: EllipticCurveHom_velusqrt(E, P)
    Traceback (most recent call last):
    ...
    NotImplementedError: only implemented for elliptic curves over finite fields

.. NOTE::

    Currently :class:`EllipticCurveHom_velusqrt` does not implement
    all methods of :class:`EllipticCurveHom`. This will hopefully
    change in the future.

AUTHORS:

- Lorenz Panny (2022)
"""

# ****************************************************************************
#       Copyright (C) 2022 Lorenz Panny
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sequence import Sequence
from sage.structure.all import coercion_model as cm

from sage.misc.misc_c import prod

from sage.structure.richcmp import op_EQ

from sage.rings.integer import Integer

from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
from sage.schemes.elliptic_curves.hom import EllipticCurveHom, compare_via_evaluation


#TODO: This is general. It should be elsewhere.
class ProductTree:
    r"""
    A simple product tree.

    INPUT:

    - ``leaves`` -- a sequence of elements in a common ring

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.hom_velusqrt import ProductTree
        sage: R.<x> = GF(101)[]
        sage: vs = [x - i for i in range(1,10)]
        sage: tree = ProductTree(vs)
        sage: tree.value()
        x^9 + 56*x^8 + 62*x^7 + 44*x^6 + 47*x^5 + 42*x^4 + 15*x^3 + 11*x^2 + 12*x + 13
        sage: tree.remainders(x^7 + x + 1)
        [3, 30, 70, 27, 58, 72, 98, 98, 23]
        sage: tree.remainders(x^100)
        [1, 1, 1, 1, 1, 1, 1, 1, 1]

    ::

        sage: vs = prime_range(100)
        sage: tree = ProductTree(vs)
        sage: tree.value().factor()
        2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41 * 43 * 47 * 53 * 59 * 61 * 67 * 71 * 73 * 79 * 83 * 89 * 97
        sage: tree.remainders(3599)
        [1, 2, 4, 1, 2, 11, 12, 8, 11, 3, 3, 10, 32, 30, 27, 48, 0, 0, 48, 49, 22, 44, 30, 39, 10]

    We can access the individual layers of the tree::

        sage: tree.layers
        [(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97),
         (6, 35, 143, 323, 667, 1147, 1763, 2491, 3599, 4757, 5767, 7387, 97),
         (210, 46189, 765049, 4391633, 17120443, 42600829, 97),
         (9699690, 3359814435017, 729345064647247, 97),
         (32589158477190044730, 70746471270782959),
         (2305567963945518424753102147331756070,)]

    .. NOTE::

        Use this class if you need the :meth:`remainders` method.
        To compute just the product, :func:`prod` is likely faster.
    """
    def __init__(self, leaves):
        r"""
        Initialize a product tree having the given ring elements
        as its leaves.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import ProductTree
            sage: vs = prime_range(100)
            sage: tree = ProductTree(vs)
        """
        V = tuple(leaves)
        self.layers = [V]
        while len(V) > 1:
            V = tuple(prod(V[i:i+2]) for i in range(0,len(V),2))
            self.layers.append(V)

    def __len__(self):
        r"""
        Return the number of leaves of this product tree.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import ProductTree
            sage: R.<x> = GF(101)[]
            sage: vs = [x - i for i in range(1,10)]
            sage: tree = ProductTree(vs)
            sage: len(tree)
            9
            sage: len(tree) == len(vs)
            True
            sage: len(tree.remainders(x^2))
            9
        """
        return len(self.layers[0])

    def value(self):
        r"""
        Return the value represented by this product tree
        (i.e., the product of all leaves).

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import ProductTree
            sage: R.<x> = GF(101)[]
            sage: vs = [x - i for i in range(1,10)]
            sage: tree = ProductTree(vs)
            sage: tree.value()
            x^9 + 56*x^8 + 62*x^7 + 44*x^6 + 47*x^5 + 42*x^4 + 15*x^3 + 11*x^2 + 12*x + 13
            sage: tree.value() == prod(vs)
            True
        """
        assert len(self.layers[-1]) == 1
        return self.layers[-1][0]

    def remainders(self, x):
        r"""
        Given a value `x`, return a list of all remainders of `x`
        modulo the leaves of this product tree.

        The base ring must support the ``%`` operator for this
        method to work.

        INPUT:

        - ``x`` -- an element of the base ring of this product tree

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import ProductTree
            sage: vs = prime_range(100)
            sage: tree = ProductTree(vs)
            sage: n = 1085749272377676749812331719267
            sage: tree.remainders(n)
            [1, 1, 2, 1, 9, 1, 7, 15, 8, 20, 15, 6, 27, 11, 2, 6, 0, 25, 49, 5, 51, 4, 19, 74, 13]
            sage: [n % v for v in vs]
            [1, 1, 2, 1, 9, 1, 7, 15, 8, 20, 15, 6, 27, 11, 2, 6, 0, 25, 49, 5, 51, 4, 19, 74, 13]
        """
        X = [x]
        for V in reversed(self.layers):
            X = [X[i//2] % V[i] for i in range(len(V))]
        return X

#TODO: This is general. It should be elsewhere.
def prod_with_derivative(pairs):
    r"""
    Given a list of pairs `(f, \partial f)` of ring elements, return
    the pair `(\prod f, \partial \prod f)`, assuming `\partial` is an
    operator obeying the standard product rule.

    This function is entirely algebraic, hence still works when the
    elements `f` and `\partial f` are all passed through some ring
    homomorphism first. (See the polynomial-evaluation example below.)

    INPUT:

    - ``pairs`` -- a sequence of tuples `(f, \partial f)` of elements
      of a common ring

    ALGORITHM:

    This function wraps the given pairs in a thin helper class that
    automatically applies the product rule whenever multiplication
    is invoked, then calls :func:`prod` on the wrapped pairs.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.hom_velusqrt import prod_with_derivative
        sage: R.<x> = ZZ[]
        sage: fs = [x^2 + 2*x + 3, 4*x + 5, 6*x^7 + 8*x + 9]
        sage: prod(fs)
        24*x^10 + 78*x^9 + 132*x^8 + 90*x^7 + 32*x^4 + 140*x^3 + 293*x^2 + 318*x + 135
        sage: prod(fs).derivative()
        240*x^9 + 702*x^8 + 1056*x^7 + 630*x^6 + 128*x^3 + 420*x^2 + 586*x + 318
        sage: F, dF = prod_with_derivative((f, f.derivative()) for f in fs)
        sage: F
        24*x^10 + 78*x^9 + 132*x^8 + 90*x^7 + 32*x^4 + 140*x^3 + 293*x^2 + 318*x + 135
        sage: dF
        240*x^9 + 702*x^8 + 1056*x^7 + 630*x^6 + 128*x^3 + 420*x^2 + 586*x + 318

    The main reason for this function to exist is that it allows us to
    *evaluate* the derivative of a product of polynomials at a point
    `\alpha` without ever fully expanding the product *as a polynomial*::

        sage: alpha = 42
        sage: F(alpha)
        442943981574522759
        sage: dF(alpha)
        104645261461514994
        sage: us = [f(alpha) for f in fs]
        sage: vs = [f.derivative()(alpha) for f in fs]
        sage: prod_with_derivative(zip(us, vs))
        (442943981574522759, 104645261461514994)
    """
    class _aux:
        def __init__(self, f, df):
            self.f, self.df = f, df

        def __mul__(self, other):
            return _aux(self.f * other.f, self.df * other.f + self.f * other.df)

        def __iter__(self):
            yield self.f
            yield self.df

    return tuple(prod(_aux(*tup) for tup in pairs))


def _choose_IJK(n):
    r"""
    Helper function to choose an "index system" for the set
    `\{1,3,5,7,...,n-2\}` where `n \geq 5` is an odd integer.

    INPUT:

    - ``n`` -- odd :class:`~sage.rings.integer.Integer` `\geq 5`

    REFERENCES: [BDLS2020]_, Examples 4.7 and 4.12

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.hom_velusqrt import _choose_IJK
        sage: IJK = _choose_IJK(101); IJK
        (range(10, 91, 20), range(1, 10, 2), range(101, 101, 2))
        sage: I,J,K = IJK
        sage: sorted([i + s*j for i in iter(I) for j in iter(J) for s in (+1,-1)] + list(iter(K)))
        [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51,
            53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99]

    TESTS::

        sage: for n in range(5,1000,2):
        ....:     I,J,K = _choose_IJK(ZZ(n))
        ....:     assert sorted([i + s*j for i in iter(I) for j in iter(J) for s in (+1,-1)] + list(iter(K))) == sorted(range(1,n,2))
    """
    if n % 2 != 1 or n < 5:
        raise ValueError('n must be odd and >= 5')
    b = (n-1).isqrt() // 2
    c = (n-1) // (4*b)
    I = range(2*b, 2*b*(2*c-1)+1, 4*b)
    J = range(1, 2*b, 2)
    K = range(4*b*c+1, n, 2)
    return I, J, K

def _points_range(rr, P, Q=None):
    r"""
    Return an iterator yielding all points `Q + [i]P` where `i` runs
    through the :class:`range` object ``rr``.

    INPUT:

    - ``rr`` -- :class:`range` object defining a sequence `S \subseteq \ZZ`
    - ``P`` -- element of an additive abelian group
    - ``Q`` -- element of the same group, or ``None``

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.hom_velusqrt import _points_range
        sage: E = EllipticCurve(GF(1123), [4,5])
        sage: P = E(1, 75)
        sage: 2*P
        (1038 : 498 : 1)
        sage: 5*P
        (236 : 598 : 1)
        sage: 8*P
        (717 : 530 : 1)
        sage: list(_points_range(range(2,10,3), P))
        [(1038 : 498 : 1), (236 : 598 : 1), (717 : 530 : 1)]
        sage: Q = E(7, 202)
        sage: Q + 2*P
        (65 : 717 : 1)
        sage: Q + 5*P
        (1119 : 788 : 1)
        sage: Q + 8*P
        (949 : 315 : 1)
        sage: list(_points_range(range(2,10,3), P, Q))
        [(65 : 717 : 1), (1119 : 788 : 1), (949 : 315 : 1)]
    """
    if not rr:
        return
    a,b,s = rr.start, rr.stop, rr.step
    R = a*P if Q is None else Q + a*P
    yield R
    sP = s*P
    for _ in range(a+s, b, s):
        yield (R := R + sP)

class FastEllipticPolynomial:
    r"""
    A class to represent and evaluate an *elliptic polynomial*,
    and optionally its derivative, in essentially square-root time.

    The elliptic polynomials computed by this class are of the form

    .. MATH::

        h_S(Z) = \prod_{i\in S} (Z - x(Q + [i]P))

    where `P` is a point of odd order `n \geq 5` and `Q` is either ``None``,
    in which case it is assumed to be `\infty`, or an arbitrary point which is
    not a multiple of `P`.

    The index set `S` is chosen as follows:

    - If `Q` is given, then `S = \{0,1,2,3,...,n-1\}`.

    - If `Q` is omitted, then `S = \{1,3,5,...,n-2\}`. Note that in this case,
      `h_{\{1,2,3,...,n-1\}}` can be computed as `h_S^2` since `n` is odd.

    INPUT:

    - ``E`` -- an elliptic curve in short Weierstraß form
    - ``n`` -- an odd integer `\geq 5`
    - ``P`` -- a point on `E`
    - ``Q`` -- a point on `E`, or ``None``

    ALGORITHM: [BDLS2020]_, Algorithm 2

    .. NOTE::

        Currently only implemented for short Weierstraß curves.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.hom_velusqrt import FastEllipticPolynomial
        sage: E = EllipticCurve(GF(71), [5,5])
        sage: P = E(4, 35)
        sage: hP = FastEllipticPolynomial(E, P.order(), P); hP
        Fast elliptic polynomial prod(Z - x(i*P) for i in range(1,n,2)) with n = 19, P = (4 : 35 : 1)
        sage: hP(7)
        19
        sage: prod(7 - (i*P).xy()[0] for i in range(1,P.order(),2))
        19

    Passing `Q` changes the index set::

        sage: Q = E(0, 17)
        sage: hPQ = FastEllipticPolynomial(E, P.order(), P, Q)
        sage: hPQ(7)
        58
        sage: prod(7 - (Q+i*P).xy()[0] for i in range(P.order()))
        58

    The call syntax has an optional keyword argument ``derivative``, which
    will make the function return the pair `(h_S(\alpha), h_S'(\alpha))`
    instead of just `h_S(\alpha)`::

        sage: hP(7, derivative=True)
        (19, 15)
        sage: R.<Z> = E.base_field()[]
        sage: HP = prod(Z - (i*P).xy()[0] for i in range(1,P.order(),2))
        sage: HP
        Z^9 + 16*Z^8 + 57*Z^7 + 6*Z^6 + 45*Z^5 + 31*Z^4 + 46*Z^3 + 10*Z^2 + 28*Z + 41
        sage: HP(7)
        19
        sage: HP.derivative()(7)
        15

    ::

        sage: hPQ(7, derivative=True)
        (58, 62)
        sage: R.<Z> = E.base_field()[]
        sage: HPQ = prod(Z - (Q+i*P).xy()[0] for i in range(P.order()))
        sage: HPQ
        Z^19 + 53*Z^18 + 67*Z^17 + 39*Z^16 + 56*Z^15 + 32*Z^14 + 44*Z^13 + 6*Z^12 + 27*Z^11 + 29*Z^10 + 38*Z^9 + 48*Z^8 + 38*Z^7 + 43*Z^6 + 21*Z^5 + 25*Z^4 + 33*Z^3 + 49*Z^2 + 60*Z
        sage: HPQ(7)
        58
        sage: HPQ.derivative()(7)
        62

    The input can be an element of any algebra over the base ring::

        sage: R.<T> = GF(71)[]
        sage: S.<t> = R.quotient(T^2)
        sage: hP(7 + t)
        15*t + 19
    """
    def __init__(self, E, n, P, Q=None):
        r"""
        Initialize this elliptic polynomial and precompute some
        input-independent data required for evaluation.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import FastEllipticPolynomial
            sage: E = EllipticCurve(GF(71), [5,5])
            sage: P = E(0, 17)
            sage: FastEllipticPolynomial(E, P.order(), P)
            Fast elliptic polynomial prod(Z - x(i*P) for i in range(1,n,2)) with n = 57, P = (0 : 17 : 1)
        """
        if any(E.a_invariants()[:-2]):
            raise NotImplementedError('only implemented for short Weierstrass curves')

        n = Integer(n)

        if Q is None:
            IJK = _choose_IJK(n)        # [1,3,5,7,...,n-4,n-2]
        else:
            IJK = _choose_IJK(2*n+1)    # [1,3,5,7,...,2n-1] = [0,1,2,3,...,n-2,n-1]

        self.base = E.base_ring()
        R, Z = self.base['x'].objgen()

        # Cassels, Lectures on Elliptic Curves, p.132
        A,B = E.a_invariants()[-2:]
        Fs = lambda X,Y: (
                (X - Y)**2,
                -2 * (X*Y + A) * (X + Y) - 4*B,
                (X*Y - A)**2 - 4*B*(X+Y),
            )

        I, J, K = IJK
        xI = (R.xy()[0] for R in _points_range(I, P, Q))
        xJ = [R.xy()[0] for R in _points_range(J, P   )]
        xK = (R.xy()[0] for R in _points_range(K, P, Q))

        self.hItree = ProductTree(Z - xi for xi in xI)

        self.EJparts = [Fs(Z,xj) for xj in xJ]

        DJ = prod(F0j for F0j,_,_ in self.EJparts)
        self.DeltaIJ = self._hI_resultant(DJ)

        self.hK = R(prod(Z - xk for xk in xK))
        self.dhK = self.hK.derivative()

        if Q is None:
            self._repr = f"Fast elliptic polynomial prod(Z - x(i*P) for i in range(1,n,2)) with {n = }, {P = }"
        else:
            self._repr = f"Fast elliptic polynomial prod(Z - x(Q+i*P) for i in range(n)) with {n = }, {P = }, {Q = }"

    def __call__(self, alpha, *, derivative=False):
        r"""
        Evaluate this elliptic polynomial at a point `\alpha`,
        and if ``derivative`` is set to ``True`` also return
        the evaluation of the derivative at `\alpha`.

        INPUT:

        - ``alpha`` -- an element of any algebra over the base ring
        - ``derivative`` -- boolean (default: ``False``)

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import FastEllipticPolynomial
            sage: E = EllipticCurve(GF(71), [5,5])
            sage: P = E(4, 35)
            sage: hP = FastEllipticPolynomial(E, P.order(), P); hP
            Fast elliptic polynomial prod(Z - x(i*P) for i in range(1,n,2)) with n = 19, P = (4 : 35 : 1)
            sage: hP(7)
            19
            sage: hP(7, derivative=True)
            (19, 15)
        """
        base = cm.common_parent(self.base, alpha)

        EJparts = [tuple(F.base_extend(base) for F in part) for part in self.EJparts]

        EJfacs = [(F0j * alpha + F1j) * alpha + F2j for F0j,F1j,F2j in EJparts]
        if not derivative:
            EJ = prod(EJfacs)
        else:
            dEJfacs = [2 * F0j * alpha + F1j for F0j,F1j,_ in EJparts]
            EJ, dEJ = prod_with_derivative(zip(EJfacs, dEJfacs))

        EJrems = self.hItree.remainders(EJ)
        R = self._hI_resultant(EJ, EJrems)
        hK = self.hK(alpha)
        res = hK * R / self.DeltaIJ

        if not derivative:
            return res

        dEJrems = self.hItree.remainders(dEJ)
        cnt = EJrems.count(0)
        if cnt == 0:
            dR = sum(R // EJrem * dEJrem for EJrem, dEJrem in zip(EJrems, dEJrems))
        elif cnt == 1:
            dR = prod(EJrem or dEJrem for EJrem, dEJrem in zip(EJrems, dEJrems))
        else:
            dR = 0
        dhK = self.dhK(alpha)
        dres = (dhK * R + hK * dR) / self.DeltaIJ

        return res, dres

    def _hI_resultant(self, poly, rems=None):
        r"""
        Internal helper function to evaluate a resultant with `h_I` quickly,
        using the product tree constructed in :meth:`__init__`.

        INPUT:

        - ``poly`` -- an element of the base ring of this product tree,
          which must be a polynomial ring supporting ``%``
        - ``rems`` -- result of ``self.hItree.remainders(poly)``, or ``None``

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import FastEllipticPolynomial
            sage: E = EllipticCurve(GF(71), [5,5])
            sage: P = E(4, 35)
            sage: hP = FastEllipticPolynomial(E, P.order(), P)
            sage: f = GF(71)['x']([5,4,3,2,1])
            sage: hP._hI_resultant(f)
            66
            sage: prod(f(r) for fi in hP.hItree.layers[0]
            ....:            for r in fi.roots(multiplicities=False))
            66

        ::

            sage: Q = E(0, 17)
            sage: hPQ = FastEllipticPolynomial(E, P.order(), P, Q)
            sage: f = GF(71)['x']([9,8,7,6,5,4,3,2,1])
            sage: hPQ._hI_resultant(f)
            36
            sage: prod(f(r) for fi in hPQ.hItree.layers[0]
            ....:            for r in fi.roots(multiplicities=False))
            36
        """
        if rems is None:
            rems = self.hItree.remainders(poly)
        r = prod(rems)
        s = -1 if len(self.hItree)%2 == 1 == poly.degree() else 1
        assert r.is_constant()
        return s * r[0]

    def __repr__(self):
        r"""
        Return a string representation of this elliptic polynomial.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import FastEllipticPolynomial
            sage: E = EllipticCurve(GF(71), [5,5])
            sage: P = E(4, 35)
            sage: FastEllipticPolynomial(E, P.order(), P)
            Fast elliptic polynomial prod(Z - x(i*P) for i in range(1,n,2)) with n = 19, P = (4 : 35 : 1)
            sage: Q = E(0, 17)
            sage: FastEllipticPolynomial(E, P.order(), P, Q)
            Fast elliptic polynomial prod(Z - x(Q+i*P) for i in range(n)) with n = 19, P = (4 : 35 : 1), Q = (0 : 17 : 1)
        """
        return self._repr


def _point_outside_subgroup(P):
    r"""
    Simple helper function to return a point on an elliptic
    curve `E` that is not a multiple of a given point `P`.
    The base field is extended if (and only if) necessary.

    INPUT:

    - ``P`` -- a point on an elliptic curve over a finite field

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.hom_velusqrt import _point_outside_subgroup
        sage: E = EllipticCurve(GF(71), [5,5])
        sage: P = E(4, 35)
        sage: Q = _point_outside_subgroup(P); Q     # random
        (14 : 11 : 1)
        sage: Q.curve()(P).discrete_log(Q)
        Traceback (most recent call last):
        ...
        ValueError: ECDLog problem has no solution (...)

    An example where `P` generates `E(\mathbb F_q)`::

        sage: E.<P> = EllipticCurve(GF(71), [5,5])
        sage: P.order() == E.cardinality()
        True
        sage: Q = _point_outside_subgroup(P); Q     # random
        (35*z2 + 7 : 24*z2 + 7 : 1)
        sage: Q.curve()(P).discrete_log(Q)
        Traceback (most recent call last):
        ...
        ValueError: ECDLog problem has no solution (...)

    An example where the group is non-cyclic:

        sage: E.<P,_> = EllipticCurve(GF(71^2), [0,1])
        sage: E.abelian_group()
        Additive abelian group isomorphic to Z/72 + Z/72 embedded in Abelian group of points on Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field in z2 of size 71^2
        sage: P = E.random_point()
        sage: Q = _point_outside_subgroup(P); Q     # random
        (18*z2 + 46 : 58*z2 + 61 : 1)
        sage: Q in E
        True
        sage: P.discrete_log(Q)
        Traceback (most recent call last):
        ...
        ValueError: ECDLog problem has no solution (...)

    .. NOTE::

        The field extension is only needed when `P` generates the
        entire rational subgroup of `E`. But in that case, the
        isogeny defined by `P` is simply `\pi-1` (where `\pi` is
        Frobenius). Thus, once `\pi-1` can be represented in Sage,
        we may just return that in
        :meth:`~sage.schemes.elliptic_curves.ell_field.EllipticCurve_field.isogeny`
        rather than insisting on using √élu.
    """
    E = P.curve()
    n = P.order()
    if n == E.order():
        d = 2 + (n == 7 and E.base_field().cardinality() == 3)
        F = E.base_field().extension(d)
        E = E.base_extend(F)
        P = E(P)
#    assert E.cardinality() > n
    for _ in range(1000):
        Q = E.random_point()
        if n*Q or not P.weil_pairing(Q,n).is_one():
            return Q
    else:
        raise NotImplementedError('could not find a point outside the kernel')

class EllipticCurveHom_velusqrt(EllipticCurveHom):
    r"""
    This class implements separable odd-degree isogenies of elliptic
    curves over finite fields using the √élu algorithm.

    The complexity is `\tilde O(\sqrt{\ell})` base-field operations,
    where `\ell` is the degree.

    REFERENCES: [BDLS2020]_

    INPUT:

    - ``E`` -- an elliptic curve over a finite field
    - ``P`` -- a point on `E` of odd order `\geq 9`
    - ``codomain`` -- codomain elliptic curve (optional)
    - ``model`` -- string (optional); input to
      :meth:`~sage.schemes.elliptic_curves.ell_field.compute_model`
    - ``Q`` -- a point on `E` outside `\langle P\rangle`, or ``None``

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.hom_velusqrt import EllipticCurveHom_velusqrt
        sage: F.<t> = GF(10009^3)
        sage: E = EllipticCurve(F, [t,t])
        sage: K = E(2154*t^2 + 5711*t + 2899, 7340*t^2 + 4653*t + 6935)
        sage: phi = EllipticCurveHom_velusqrt(E, K); phi
        Elliptic-curve isogeny (using √élu) of degree 601:
          From: Elliptic Curve defined by y^2 = x^3 + t*x + t over Finite Field in t of size 10009^3
          To:   Elliptic Curve defined by y^2 = x^3 + (263*t^2+3173*t+4759)*x + (3898*t^2+6111*t+9443) over Finite Field in t of size 10009^3
        sage: phi(K)
        (0 : 1 : 0)
        sage: P = E(2, 3163*t^2 + 7293*t + 5999)
        sage: phi(P)
        (6085*t^2 + 855*t + 8720 : 8078*t^2 + 9889*t + 6030 : 1)
        sage: Q = E(6, 5575*t^2 + 6607*t + 9991)
        sage: phi(Q)
        (626*t^2 + 9749*t + 1291 : 5931*t^2 + 8549*t + 3111 : 1)
        sage: phi(P + Q)
        (983*t^2 + 4894*t + 4072 : 5047*t^2 + 9325*t + 336 : 1)
        sage: phi(P) + phi(Q)
        (983*t^2 + 4894*t + 4072 : 5047*t^2 + 9325*t + 336 : 1)

    TESTS:

    Check on a random example that the isogeny is a well-defined
    group homomorphism with the correct kernel::

        sage: from sage.schemes.elliptic_curves.hom_velusqrt import _random_example_for_testing
        sage: E, K = _random_example_for_testing()
        sage: phi = EllipticCurveHom_velusqrt(E, K)
        sage: not phi(K)
        True
        sage: not phi(randrange(2^99) * K)
        True
        sage: P = E.random_point()
        sage: phi(P) in phi.codomain()
        True
        sage: Q = E.random_point()
        sage: phi(Q) in phi.codomain()
        True
        sage: phi(P + Q) == phi(P) + phi(Q)
        True

    Check that the isogeny preserves the field of definition::

        sage: Sequence(K).universe() == phi.domain().base_field()
        True
        sage: phi.codomain().base_field() == phi.domain().base_field()
        True

    Check that the isogeny affects the Weil pairing in the correct way::

        sage: m = lcm(P.order(), Q.order())
        sage: e1 = P.weil_pairing(Q, m)
        sage: e2 = phi(P).weil_pairing(phi(Q), m)
        sage: e2 == e1^phi.degree()
        True

    Check that the isogeny matches (up to isomorphism) the one from
    :class:`~sage.schemes.elliptic_curves.ell_curve_isogeny.EllipticCurveIsogeny`::

        sage: psi = EllipticCurveIsogeny(E, K)
        sage: check = lambda iso: all(iso(psi(Q)) == phi(Q) for Q in E.gens())
        sage: any(map(check, psi.codomain().isomorphisms(phi.codomain())))
        True

    .. SEEALSO::

        :class:`~sage.schemes.elliptic_curves.ell_curve_isogeny.EllipticCurveIsogeny`
    """
    def __init__(self, E, P, *, codomain=None, model=None, Q=None):
        r"""
        Initialize this √élu isogeny from a kernel point of odd order.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import EllipticCurveHom_velusqrt
            sage: E = EllipticCurve(GF(71), [5,5])
            sage: P = E(-2, 22)
            sage: EllipticCurveHom_velusqrt(E, P)
            Elliptic-curve isogeny (using √élu) of degree 19:
              From: Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 71
              To:   Elliptic Curve defined by y^2 = x^3 + 13*x + 11 over Finite Field of size 71

        ::

            sage: E.<P> = EllipticCurve(GF(419), [1,0])
            sage: K = 4*P
            sage: EllipticCurveHom_velusqrt(E, K)
            Elliptic-curve isogeny (using √élu) of degree 105:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 419
              To:   Elliptic Curve defined by y^2 = x^3 + 301*x + 86 over Finite Field of size 419
            sage: E2 = EllipticCurve(GF(419), [0,6,0,385,42])
            sage: EllipticCurveHom_velusqrt(E, K, codomain=E2)
            Elliptic-curve isogeny (using √élu) of degree 105:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 419
              To:   Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 385*x + 42 over Finite Field of size 419
            sage: EllipticCurveHom_velusqrt(E, K, model="montgomery")
            Elliptic-curve isogeny (using √élu) of degree 105:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 419
              To:   Elliptic Curve defined by y^2 = x^3 + 6*x^2 + x over Finite Field of size 419

        Note that the implementation in fact also works in almost all
        cases when the degree is `5` or `7`. The reason we restrict to
        degrees `\geq 9` is that (only!) when trying to compute a
        `7`-isogeny from a rational point on an elliptic curve defined
        over `\GF{3}`, the point `Q` required in the formulas has to be
        defined over a cubic extension rather than an at most quadratic
        extension, which can result in the constructed isogeny being
        irrational. See :trac:`34467`. The assertion in the following
        example currently fails if the minimum degree is lowered::

            sage: E = EllipticCurve(GF(3), [2,1])
            sage: P, = E.gens()
            sage: P.order()
            7
            sage: psi = E.isogeny(P)
            sage: phi = E.isogeny(P, algorithm='velusqrt')              # not tested
            sage: phi._Q.base_ring()                                    # not tested
            Finite Field in z3 of size 3^3
            sage: assert phi.codomain().is_isomorphic(psi.codomain())   # not tested
        """
        if not isinstance(E, EllipticCurve_finite_field):
            raise NotImplementedError('only implemented for elliptic curves over finite fields')

        if codomain is not None and model is not None:
            raise ValueError('cannot specify a codomain curve and model name simultaneously')

        try:
            self._raw_domain = E.short_weierstrass_model()
        except ValueError:
            raise NotImplementedError('only implemented for curves having a short Weierstrass model')
        self._pre_iso = E.isomorphism_to(self._raw_domain)

        try:
            P = E(P)
        except TypeError:
            raise ValueError('given kernel point P does not lie on E')
        self._P = self._pre_iso(P)

        self._degree = self._P.order()
        if self._degree % 2 != 1 or self._degree < 9:
            raise NotImplementedError('only implemented for odd degrees >= 9')

        if Q is not None:
            self._Q = E(Q)
            EE = E
        else:
            self._Q = _point_outside_subgroup(self._P)  # may extend base field
            EE = self._Q.curve()
            self._P = EE(self._P)

        self._base_ring = EE.base_ring()

        self._h0 = FastEllipticPolynomial(EE, self._degree, self._P)
        self._h1 = FastEllipticPolynomial(EE, self._degree, self._P, self._Q)

        self._domain = E
        self._compute_codomain(model=model)

        if codomain is not None:
            self._post_iso = self._codomain.isomorphism_to(codomain) * self._post_iso
            self._codomain = codomain

        super().__init__(self._domain, self._codomain)

    def _raw_eval(self, x, y=None):
        r"""
        Evaluate the "inner" √élu isogeny (i.e., without applying
        pre- and post-isomorphism) at either just an `x`-coordinate
        or a pair `(x,y)` of coordinates.

        If the given point lies in the kernel, the empty tuple
        ``()`` is returned.

        No checking of the input coordinates is performed.

        ALGORITHM:

        - [Ren2018]_, Theorem 1
        - :class:`FastEllipticPolynomial`

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import EllipticCurveHom_velusqrt
            sage: E = EllipticCurve(GF(65537), [1,1])
            sage: P = E(2112, 803)
            sage: phi = EllipticCurveHom_velusqrt(E, P, Q=(32924,0))
            sage: phi._raw_domain is E
            True
            sage: phi._raw_codomain
            Elliptic Curve defined by y^2 = x^3 + ... over Finite Field of size 65537
            sage: Q = E(42, 15860)
            sage: phi._raw_eval(Q.xy()[0])
            11958
            sage: phi._raw_eval(*Q.xy())
            (11958, 42770)
            sage: phi._raw_codomain.defining_polynomial()(*phi._raw_eval(*Q.xy()), 1)
            0

        No checking is performed::

            sage: E.defining_polynomial()(123, 456, 1)
            50907
            sage: phi._raw_eval(123, 456)
            (3805, 29941)
        """
        if y is None:
            h0 = self._h0(x)
            h1 = self._h1(x)
        else:
            h0, h0d = self._h0(x, derivative=True)
            h1, h1d = self._h1(x, derivative=True)

#        assert h0 == prod(x - (        i*self._P).xy()[0] for i in range(1,self._P.order(),2))
#        assert h1 == prod(x - (self._Q+i*self._P).xy()[0] for i in range(  self._P.order()  ))

        if not h0:
            return ()

        xx = h1 / h0**2

        if y is None:
            return xx

#        assert h0d == sum(prod(x - (        i*self._P).xy()[0] for i in range(1,self._P.order(),2) if i!=j) for j in range(1,self._P.order(),2))
#        assert h1d == sum(prod(x - (self._Q+i*self._P).xy()[0] for i in range(  self._P.order()  ) if i!=j) for j in range(  self._P.order()  ))

        yy = y * (h1d - 2 * h1 / h0 * h0d) / h0**2

        return xx, yy

    def _compute_codomain(self, model=None):
        r"""
        Helper method to compute the codomain of this √élu isogeny
        once the data for :meth:`_raw_eval` has been initialized.

        Called by the constructor.

        INPUT:

        - ``model`` -- string (optional); input to
          :meth:`~sage.schemes.elliptic_curves.ell_field.compute_model`

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import EllipticCurveHom_velusqrt
            sage: E = EllipticCurve(GF(71), [0,5,0,1,0])
            sage: P = E(4, 19)
            sage: phi = EllipticCurveHom_velusqrt(E, P)
            sage: phi._raw_codomain
            Elliptic Curve defined by y^2 = x^3 + ... over Finite Field of size 71
            sage: phi._codomain
            Elliptic Curve defined by y^2 = x^3 + 8*x + 34 over Finite Field of size 71
            sage: phi.codomain()
            Elliptic Curve defined by y^2 = x^3 + 8*x + 34 over Finite Field of size 71

        Passing a ``model`` parameter is supported::

            sage: phi._compute_codomain('montgomery')
            sage: phi
            Elliptic-curve isogeny (using √élu) of degree 19:
              From: Elliptic Curve defined by y^2 = x^3 + 5*x^2 + x over Finite Field of size 71
              To:   Elliptic Curve defined by y^2 = x^3 + 40*x^2 + x over Finite Field of size 71

        TESTS::

            sage: F.<t> = GF(5^2)
            sage: E = EllipticCurve([3*t, 2*t+4, 3*t+2, t+4, 3*t])
            sage: K = E(3*t, 2)
            sage: EllipticCurveHom_velusqrt(E, K)   # indirect doctest
            Elliptic-curve isogeny (using √élu) of degree 19:
              From: Elliptic Curve defined by y^2 + 3*t*x*y + (3*t+2)*y = x^3 + (2*t+4)*x^2 + (t+4)*x + 3*t over Finite Field in t of size 5^2
              To:   Elliptic Curve defined by y^2 = x^3 + (4*t+3)*x + 2 over Finite Field in t of size 5^2
        """
        R, Z = self._base_ring['Z'].objgen()
        poly = self._raw_domain.two_division_polynomial().monic()(Z)

        f = 1
        for g,_ in poly.factor():
            if g.degree() == 1:
                f *= Z - self._raw_eval(-g[0])
            else:
                K, X0 = self._base_ring.extension(g,'T').objgen()
                imX0 = self._raw_eval(X0)
                try:
                    imX0 = imX0.polynomial()    # K is a FiniteField
                except AttributeError:
                    imX0 = imX0.lift()          # K is a PolynomialQuotientRing
                V = R['V'].gen()
                f *= (Z - imX0(V)).resultant(g(V))

        a6,a4,a2,_ = f.monic().list()

        self._raw_codomain = EllipticCurve(self._domain.base_ring(), [0,a2,0,a4,a6])

        if model is None:
            model = 'short_weierstrass'

        from sage.schemes.elliptic_curves.ell_field import compute_model
        self._codomain = compute_model(self._raw_codomain, model)
        self._post_iso = self._raw_codomain.isomorphism_to(self._codomain)

    def _eval(self, P):
        r"""
        Evaluate this √élu isogeny at a point.

        INPUT:

        - ``P`` -- point on the domain, defined over any algebra over the base field

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import EllipticCurveHom_velusqrt
            sage: E = EllipticCurve(GF(71), [0,5,0,1,0])
            sage: K = E(4, 19)
            sage: phi = EllipticCurveHom_velusqrt(E, K, model='montgomery')
            sage: phi
            Elliptic-curve isogeny (using √élu) of degree 19:
              From: Elliptic Curve defined by y^2 = x^3 + 5*x^2 + x over Finite Field of size 71
              To:   Elliptic Curve defined by y^2 = x^3 + 40*x^2 + x over Finite Field of size 71
            sage: phi(K)
            (0 : 1 : 0)
            sage: phi(5*K)
            (0 : 1 : 0)
            sage: phi(E(0))
            (0 : 1 : 0)
            sage: phi(E(0,0))
            (0 : 0 : 1)
            sage: phi(E(7,13))
            (70 : 31 : 1)

        TESTS::

            sage: P,Q = (E.random_point() for _ in 'PQ')
            sage: assert phi(P) in phi.codomain()
            sage: assert phi(Q) in phi.codomain()
            sage: assert phi(P + Q) == phi(P) + phi(Q)

        Randomized testing::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import EllipticCurveHom_velusqrt
            sage: from sage.schemes.elliptic_curves.hom_velusqrt import _random_example_for_testing
            sage: E, K = _random_example_for_testing()
            sage: phi = EllipticCurveHom_velusqrt(E, K)
            sage: phi.degree() == K.order()
            True
            sage: P = E.random_point()
            sage: phi(P) in phi.codomain()
            True
            sage: Q = E.random_point()
            sage: phi(Q) in phi.codomain()
            True
            sage: phi(P + Q) == phi(P) + phi(Q)
            True
        """
        if self._domain.defining_polynomial()(*P):
            raise ValueError(f'{P} not on {self._domain}')

        k = Sequence(P).universe()

        if not P:
            return self._codomain(0).change_ring(k)

        P = self._pre_iso._eval(P)

        xy = self._raw_eval(*P.xy())

        if xy == ():
            return self._codomain(0).change_ring(k)

        return self._post_iso._eval(Sequence(xy, k) + [1])

    _call_ = _eval

    def _repr_(self):
        r"""
        Return basic information about this √élu isogeny as a string.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import EllipticCurveHom_velusqrt
            sage: E.<P> = EllipticCurve(GF(71), [5,5])
            sage: phi = EllipticCurveHom_velusqrt(E, P)
            sage: phi   # indirect doctest
            Elliptic-curve isogeny (using √élu) of degree 57:
              From: Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 71
              To:   Elliptic Curve defined by y^2 = x^3 + 19*x + 45 over Finite Field of size 71
        """
        return f'Elliptic-curve isogeny (using √élu) of degree {self._degree}:' \
                f'\n  From: {self._domain}' \
                f'\n  To:   {self._codomain}'

    @staticmethod
    def _comparison_impl(left, right, op):
        r"""
        Compare a √élu isogeny to another elliptic-curve morphism.

        Called by :meth:`EllipticCurveHom._richcmp_`.

        INPUT:

        - ``left, right`` -- :class:`~sage.schemes.elliptic_curves.hom.EllipticCurveHom` objects

        ALGORITHM:

        :func:`~sage.schemes.elliptic_curves.hom.compare_via_evaluation`

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import EllipticCurveHom_velusqrt
            sage: E = EllipticCurve(GF(101), [5,5,5,5,5])
            sage: phi = EllipticCurveHom_velusqrt(E, E.lift_x(11)); phi
            Elliptic-curve isogeny (using √élu) of degree 59:
              From: Elliptic Curve defined by y^2 + 5*x*y + 5*y = x^3 + 5*x^2 + 5*x + 5 over Finite Field of size 101
              To:   Elliptic Curve defined by y^2 = x^3 + 15*x + 25 over Finite Field of size 101
            sage: psi = EllipticCurveHom_velusqrt(E, E.lift_x(-1)); psi
            Elliptic-curve isogeny (using √élu) of degree 59:
              From: Elliptic Curve defined by y^2 + 5*x*y + 5*y = x^3 + 5*x^2 + 5*x + 5 over Finite Field of size 101
              To:   Elliptic Curve defined by y^2 = x^3 + 15*x + 25 over Finite Field of size 101
            sage: phi == psi
            True
        """
        if op != op_EQ:
            return NotImplemented
        return compare_via_evaluation(left, right)


def _random_example_for_testing():
    r"""
    Function to generate somewhat random valid √élu inputs
    for testing purposes.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.hom_velusqrt import _random_example_for_testing
        sage: E, K = _random_example_for_testing()
        sage: E                 # random
        Elliptic Curve defined by y^2 + (t^3+6*t^2)*x*y + (t^3+3*t^2+2*t+2)*y = x^3 + (6*t^3+2*t^2+t)*x^2 + (3*t^3+2*t^2+6*t+1)*x + (t^3+2*t^2+2) over Finite Field in t of size 7^4
        sage: E.short_weierstrass_model()
        Elliptic Curve defined by y^2 = x^3 + ... over Finite Field ...
        sage: K                 # random
        (3*t^3 + 4*t^2 + 4*t + 3 : 6*t^3 + 5*t^2 + 5*t : 1)
        sage: K.order()         # random
        101
        sage: K in E
        True
        sage: K.order() % 2
        1
        sage: 5 <= K.order()
        True
    """
    from sage.all import prime_range, choice, randrange, GF, gcd
    while True:
        p = choice(prime_range(2, 100))
        e = randrange(1,5)
        F,t = GF((p,e),'t').objgen()
        try:
            E = EllipticCurve([F.random_element() for _ in range(5)])
        except ArithmeticError:
            continue
        try:
            E.short_weierstrass_model()
        except ValueError:
            continue
        if E.cardinality() < 9:
            continue
        A = E.abelian_group()
        ds = max(A.invariants()).prime_to_m_part(2).divisors()
        ds = [d for d in ds if 9 <= d < 1000]
        if ds:
            deg = choice(ds)
            break
    G = A.torsion_subgroup(deg)
    while True:
        v = [randrange(deg) for _ in range(G.ngens())]
        if gcd([deg] + v) == 1:
            break
    K = G(v).element()
    assert K.order() == deg
    return E, K

