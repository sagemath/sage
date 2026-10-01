r"""
Analytic certificates for the Birch and Swinnerton-Dyer formula

This is an internal module for :mod:`sage.schemes.elliptic_curves.BSD`.
In particular, the rational reconstruction below does *not* assume that
the analytic order of Sha is an integer, or that BSD is true.

All caches belong to one computation.  A numerical value returned by
``Sha.an()`` is never entered in a cache of certified values.

AUTHORS:

- Barinder Banwait (2026): analytic certification for the BSD rewrite

EXAMPLES::

    sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
    sage: data = _AnalyticData(EllipticCurve('11a3'), rigorous=True)
    sage: data.rank, data.sha(EllipticCurve('11a3'))
    (0, 1)
    sage: data = _AnalyticData(EllipticCurve('37a1'), rigorous=True)
    sage: data.rank, data.sha(EllipticCurve('37a1'))
    (1, 1)
"""

# Distributed under the terms of the GNU General Public License (GPL).

from sage.rings.real_arb import RealBallField

from sage.arith.misc import kronecker
from sage.libs.pari import pari
from sage.rings.complex_arb import ComplexBallField
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ


def _unique_integer(ball):
    """
    Return the unique integer in ``ball``, or ``None``.

    No rounding of an approximate midpoint is used.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.bsd_certificates import _unique_integer
        sage: _unique_integer(RBF(7, rad=1/4))
        7
        sage: _unique_integer(RBF(7, rad=2)) is None
        True
        sage: _unique_integer(RBF(7/2, rad=1/8)) is None
        True
    """
    if not ball.is_finite():
        return None
    lo = ball.lower().ceil()
    hi = ball.upper().floor()
    return ZZ(lo) if lo == hi else None


def _derivative_ball(E, prec):
    r"""
    Enclose `L'(E,1)` for a curve with root number -1.

    Use the series in :meth:`~sage.schemes.elliptic_curves.lseries_ell.Lseries_ell.deriv_at1`,
    evaluating `E_1(x)=\Gamma(0,x)` in Arb.  Its geometric tail bound is
    valid for the chosen truncation `k\geq\sqrt N`.  Ball arithmetic
    accounts for roundoff, including that in the tail bound.

    Indeed `|a_n|\leq d(n)\sqrt n\leq2n`, and
    `E_1(x)\leq e^{-x}/x`.  With `a=2\pi/\sqrt N` and `n>k`, each
    omitted term has absolute value at most `4e^{-an}/(an)\leq2e^{-an}`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.bsd_certificates import _derivative_ball
        sage: L = _derivative_ball(EllipticCurve('37a1'), 64)
        sage: L > 0
        True
        sage: abs(L - RBF('0.3059997738340523')) < RBF('1e-15')
        True
        sage: _derivative_ball(EllipticCurve('11a1'), 64)
        Traceback (most recent call last):
        ...
        ValueError: the derivative certificate requires root number -1
    """
    if E.root_number() != -1:
        raise ValueError("the derivative certificate requires root number -1")
    R = RealBallField(prec + 32)
    a = 2 * R.pi() / R(E.conductor()).sqrt()
    k = max(R(E.conductor()).sqrt().upper().ceil(),
            ((prec + 16) * R(2).log() / a).upper().ceil())
    an = E.anlist(k)
    value = 2 * sum((R(an[n]) / n * R(0).gamma(a * n)
                     for n in range(1, k + 1)), R.zero())
    tail = 2 * (-a * (k + 1)).exp() / (1 - (-a).exp())
    return value.add_error(tail)


def _real_period_ball(E, prec):
    r"""
    Enclose the total real Neron period of a minimal curve.

    The algebraic inputs to the AGM are the exact inputs also used by
    :mod:`sage.schemes.elliptic_curves.period_lattice`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.bsd_certificates import _real_period_ball
        sage: for label in ['11a1', '37a1']:
        ....:     E = EllipticCurve(label)
        ....:     omega = _real_period_ball(E, 100)
        ....:     approx = E.period_lattice().real_period(prec=100) * E.real_components()
        ....:     assert abs(omega - omega.parent()(approx)) < RBF('1e-27')
    """
    R = RealBallField(prec)
    lattice = E.period_lattice()
    if E.discriminant() > 0:
        a, b, _ = (R(x) for x in lattice._abc)
        omega = R.pi() / a.agm(b)
    else:
        a = ComplexBallField(prec)(lattice._abc[0])
        omega = R.pi() / abs(a).agm(abs(a.real()))
    return E.real_components() * omega


def _finite_height_terms(P):
    r"""
    Return rational coefficients of logarithms in the finite local height.

    The curve must be globally minimal.  This is Silverman's local-height
    algorithm [Sil1988]_, Section 5, with Sage's doubled normalization.
    Factor only the discriminant; the denominator contributes one logarithm.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.bsd_certificates import _finite_height_terms
        sage: _finite_height_terms(EllipticCurve('37a1')(0, 0))
        []
        sage: E = EllipticCurve('123a1')
        sage: _finite_height_terms(E.gens()[0])
        [(3, -4/5)]
    """
    E = P.curve()
    a1, a2, a3, a4, _ = E.a_invariants()
    b2, b4, b6, b8 = E.b_invariants()
    x, y = P.xy()
    den = x.denominator()
    terms = [(den, QQ.one())] if den != 1 else []
    for p, n in E.discriminant().factor():
        if den % p == 0:
            continue
        A = (3*x**2 + 2*a2*x + a4 - a1*y).valuation(p)
        B = (2*y + a1*x + a3).valuation(p)
        C = (3*x**4 + b2*x**3 + 3*b4*x**2 + 3*b6*x + b8).valuation(p)
        if A <= 0 or B <= 0:
            r = QQ.zero()
        elif E.c4().valuation(p) == 0:
            m = min(B, QQ(n)/2)
            r = -m*(n-m)/n
        elif C >= 3*B:
            r = -QQ(2)*B/3
        else:
            r = -QQ(C)/4
        if r:
            terms.append((p, r))
    return terms


def _height_ball(P, prec):
    r"""
    Enclose the canonical height of a nonzero point on a minimal curve.

    Let `F=(F_0,F_1)` be the homogeneous quartic duplication map on
    x-coordinates.  On each of the charts `(t,1)` and `(1,t)`, `|t|\leq1`,
    a polynomial Bezout identity gives an exact positive lower bound for
    `\|F\|_\infty`; coefficient sums give an upper bound.  By homogeneity
    these give bounds `L,U` on the unit 1-norm sphere.  Set
    `B=\max(|\log L|,|\log U|)`.

    Normalize the projective coordinates after each duplication.  The
    archimedean local height is the logarithm of the initial 1-norm plus
    `\sum_{j\geq0}4^{-j-1}\log\|F(X_j,Z_j)\|_1`.  After `n` terms the
    absolute tail is at most `B/(3\cdot4^n)`.  Add the exact finite local
    contributions from :func:`_finite_height_terms`.

    This also explains why changing from max-norm to 1-norm does not change
    the limit.  No numerical height-error allowance is assumed.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.bsd_certificates import _height_ball
        sage: P = EllipticCurve('37a1')(0, 0)
        sage: H = _height_ball(P, 80)
        sage: H > 0
        True
        sage: (_height_ball(3*P, 80) - 9*H).contains_zero()
        True
        sage: abs(H - H.parent()(P.height(100))) < RBF('1e-23')
        True
    """
    E = P.curve()
    if not E.is_minimal() or P.is_zero():
        raise ValueError("a nonzero point on a minimal curve is required")
    b2, b4, b6, b8 = E.b_invariants()
    t = PolynomialRing(QQ, 't').gen()
    f = t**4 - b4*t**2 - 2*b6*t - b8
    g = 4*t**3 + b2*t**2 + 2*b4*t + b6
    fr = 1 - b4*t**2 - 2*b6*t**3 - b8*t**4
    gr = 4*t + b2*t**2 + 2*b4*t**3 + b6*t**4
    lower = []
    for ff, gg in ((f, g), (fr, gr)):
        d, a, b = ff.xgcd(gg)
        if d.degree() != 0 or not d:
            raise ArithmeticError("singular duplication map")
        lower.append(abs(d[0]) / (sum(map(abs, a)) + sum(map(abs, b))))
    L = min(lower) / 16
    U = 2 * max(sum(map(abs, f)), sum(map(abs, g)))
    R = RealBallField(prec + 32)
    B = max(abs(R(L).log()).upper(), abs(R(U).log()).upper())
    n = max(1, ((R(B).log()/R(2).log() + prec + 8)/2).upper().ceil())
    workprec = 2*prec + 2*n + 64
    while True:
        R = RealBallField(workprec)
        scale = abs(P[0]) + 1
        X, Z = R(P[0]/scale), R(1/scale)
        h = R(scale).log()
        for j in range(n):
            FX = X**4 - b4*X**2*Z**2 - 2*b6*X*Z**3 - b8*Z**4
            FZ = 4*X**3*Z + b2*X**2*Z**2 + 2*b4*X*Z**3 + b6*Z**4
            s = abs(FX) + abs(FZ)
            if not s > 0:
                break
            h += s.log() / ZZ(4)**(j+1)
            X, Z = FX/s, FZ/s
        else:
            h += sum((r*R(p).log() for p, r in _finite_height_terms(P)), R.zero())
            h = h.add_error(R(B)/(3*ZZ(4)**n))
            if h.is_finite() and R(h.rad()) < R(2)**(-prec):
                return h
        workprec *= 2


class _AnalyticData:
    r"""
    Per-call analytic inputs, with separate rigorous and assumed paths.

    The rank-one lattice follows from [GZ1986]_, Chapter V, Theorem 2.1.
    Write `A=L'(E,1)t^2/(\Omega_E R C)`, where the real period is total
    and `R` uses Sage's canonical height.  For an odd fundamental Heegner
    discriminant `D<-4`, let `s_D=\sqrt{|D|}L(E^D,1)/\Omega_E^-`.
    This is an exact quadratic-character sum of minus modular symbols.

    Gross--Zagier gives `4\hat h_{\QQ}(y_D)/R=A/\lambda`, where
    `\lambda=t^2/(2c^2s_D C)`.  Here the height over the quadratic field
    is twice the normalized height, and
    `\operatorname{area}(E)=\Omega_E\Omega_E^-/2`.
    Since the twist has rank zero, `2y_D` belongs to `E(\QQ)` modulo
    torsion.  Consequently `A/\lambda` is an integer square, without
    assuming the BSD formula.  We use integrality for reconstruction
    and squareness as an additional consistency check.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
        sage: E = EllipticCurve('91b1')
        sage: a = _AnalyticData(E, True)
        sage: a.sha(E), a.certificates[E.a_invariants()]['certified']
        (1, True)
        sage: b = _AnalyticData(E, False)
        sage: b.sha(E), b.certificates[E.a_invariants()]['certified']
        (1, False)
        sage: _AnalyticData(EllipticCurve('389a1'), True).rank is None
        True
    """
    def __init__(self, E, rigorous):
        self.E = E.minimal_model()
        self.rigorous = bool(rigorous)
        self._symbols = {}
        self._sha = {}
        self._gens = {}
        self._derivatives = {}
        self._heegner = {}
        self.certificates = {}
        self.rank_certificate = {'certified': False}
        if not self.rigorous:
            self.rank = ZZ(self.E.analytic_rank())
            self.rank_certificate['method'] = 'assumed analytic_rank()'
        elif self.E.root_number() == 1:
            value = self.symbol(self.E, 1, 0)
            self.rank = ZZ.zero() if value else None
            self.rank_certificate = {'certified': bool(value),
                                     'method': 'exact modular symbol', 'value': value}
        else:
            self.rank = None
            for prec in (64, 128, 256, 512):
                value = self.derivative(prec)
                if not value.contains_zero():
                    self.rank = ZZ.one()
                    self.rank_certificate = {'certified': True,
                                             'method': 'root number and derivative ball',
                                             'value': value}
                    break

    def symbol(self, E, sign, cusp):
        """
        Evaluate an exactly normalized modular symbol from infinity to ``cusp``.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
            sage: E = EllipticCurve('11a3')
            sage: _AnalyticData(E, True).symbol(E, 1, 0)
            1/25
        """
        key = (E.a_invariants(), sign)
        if key not in self._symbols:
            self._symbols[key] = pari.msfromell(E.pari_curve(), sign)
        M, symbol = self._symbols[key]
        return QQ(pari.mseval(M, symbol, [infinity, cusp]))

    def derivative(self, prec):
        """
        Cache derivative enclosures shared by the isogeny class.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
            sage: a = _AnalyticData(EllipticCurve('37a1'), True)
            sage: a.derivative(64) > 0
            True
        """
        if prec not in self._derivatives:
            self._derivatives[prec] = _derivative_ball(self.E, prec)
        return self._derivatives[prec]

    def gens(self, E):
        """
        Obtain a saturated basis, independently of cached numerical regulators.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
            sage: E = EllipticCurve('37a1')
            sage: len(_AnalyticData(E, True).gens(E))
            1
        """
        key = E.a_invariants()
        if key not in self._gens:
            if self.rank == 0:
                self._gens[key] = []
            else:
                # Analytic rank one already supplies the upper bound by
                # Gross--Zagier and Kolyvagin (conditionally on that rank
                # in the assumed mode).  Use gens only to *find* a point;
                # explicitly prove saturation below, even for cached points.
                points = E.gens(proof=False)
                if not points:
                    raise NotImplementedError('no non-torsion point was found')
                if len(points) != 1 or points[0].has_finite_order():
                    raise ArithmeticError("rank-one certificate and Mordell-Weil basis disagree")
                self._gens[key] = E.saturation(points)[0]
        return self._gens[key]

    def sha(self, E):
        """
        Return the exact (or explicitly assumed) analytic order of Sha.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
            sage: E = EllipticCurve('681b1')
            sage: _AnalyticData(E, True).sha(E)
            9
        """
        key = E.a_invariants()
        if key in self._sha:
            return self._sha[key]
        certificate = {'certified': self.rigorous}
        t, C = E.torsion_order(), E.tamagawa_product()
        if not self.rigorous:
            value = QQ(E.sha().an())
            certificate['method'] = 'assumed sha().an()'
        elif self.rank == 0:
            ratio = self.symbol(E, 1, 0) / E.real_components()
            value = ratio * t**2 / C
            certificate.update(method='exact modular symbols', L_ratio=ratio)
        elif self.rank == 1:
            heegner = self.heegner_data(E)
            lattice = heegner['lattice']
            P = self.gens(E)[0]
            prec = 64
            while True:
                omega = _real_period_ball(E, prec + 32)
                height = _height_ball(P, prec + 16)
                enclosure = self.derivative(prec) * t**2 / (omega * height * C)
                index_square = _unique_integer(enclosure/lattice)
                if index_square is not None:
                    if index_square <= 0 or not index_square.is_square():
                        raise ArithmeticError("Gross-Zagier reconstruction is not a positive square")
                    value = index_square * lattice
                    break
                prec *= 2
            certificate.update(heegner)
            certificate.update(method='Gross-Zagier rational reconstruction',
                               index_square=index_square,
                               enclosure=enclosure, generator=P)
        else:
            raise ValueError("analytic rank zero or one is required")
        if value <= 0:
            raise ArithmeticError("the analytic order of Sha is not positive")
        self._sha[key] = value
        self.certificates[key] = certificate
        return value

    def heegner_data(self, E):
        r"""
        Compute an exact Gross--Zagier lattice on the original curve.

        This calculation also supplies the Heegner point used in the
        Kolyvagin--Jetchev bound.  It is needed in conditional mode too:
        an assumed analytic order alone is not the index over a specified
        quadratic field.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
            sage: E = EllipticCurve('37a1')
            sage: a = _AnalyticData(E, True)
            sage: a.heegner_data(E)
            {'discriminant': -7, 'lattice': 1/4, 'manin_constant': 1, 'twist_sum': 2}

        The unconditional Manin algorithm also handles nonoptimal curves::

            sage: E = EllipticCurve('11a3')
            sage: E.pari_curve().ellmaninconstant()
            5
        """
        if self.rank != 1:
            raise ValueError('a Heegner lattice requires analytic rank one')
        key = E.a_invariants()
        if key not in self._heegner:
            # PARI 2.17 has only the unconditional modular-symbol algorithm;
            # later versions offer a non-default, conditional table lookup.
            c = ZZ(E.pari_curve().ellmaninconstant())
            D = ZZ(-7)
            while True:
                if D.is_fundamental_discriminant() and all(
                        kronecker(D, q) == 1 for q in E.conductor().prime_divisors()):
                    s = sum(kronecker(D, a) * self.symbol(E, -1, QQ(a)/abs(D))
                            for a in range(1, abs(D)))
                    if s:
                        break
                D -= 4
            lattice = QQ(E.torsion_order()**2)/(2*c**2*s*E.tamagawa_product())
            self._heegner[key] = {'discriminant': D, 'manin_constant': c,
                                  'twist_sum': s, 'lattice': lattice}
        return self._heegner[key]


def _padic_sha_bound(E, p, analytic):
    r"""
    Compute the Iwasawa upper bound with explicit analytic inputs.

    The caller checks the hypotheses of [SW2013]_ and [Wuthrich2014]_.
    This is the rank-at-most-one specialization of ``Sha.an_padic``, with
    certified generators and exactly normalized PARI modular symbols.
    No quadratic-twist optimization or cached ``an_padic`` value is used.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData, _padic_sha_bound
        sage: E = EllipticCurve('11a3')
        sage: _padic_sha_bound(E, 5, _AnalyticData(E, True))
        0
        sage: E = EllipticCurve('123a1')
        sage: _padic_sha_bound(E, 41, _AnalyticData(E, True))  # long time
        0
    """
    from sage.modules.free_module_element import vector

    from sage.misc.cachefunc import cached_method
    from sage.rings.padics.factory import Qp
    from sage.schemes.elliptic_curves.padic_lseries import (
        pAdicLseriesOrdinary,
        pAdicLseriesSupersingular,
    )

    if analytic.rank == 0:
        # Exactly the untwisted rank-zero shortcut in Sha.an_padic.
        return analytic.sha(E).valuation(p)

    base = pAdicLseriesOrdinary if E.is_ordinary(p) else pAdicLseriesSupersingular

    class ExactSymbols(base):
        # Use the existing series/precision algorithms with a private set of
        # symbols.  Initializing these attributes directly avoids constructing
        # and then discarding an uncertified L_ratio normalization in the base
        # constructor.  No curve-level or public p-adic L-series cache is changed.
        def __init__(self):
            self._E = E
            self._p = ZZ(p)
            self._implementation = 'pari'
            self._normalize = 'L_ratio'
            self._modular_symbol = lambda r: analytic.symbol(E, 1, r)
            self._negative_modular_symbol = lambda r: analytic.symbol(E, -1, r)

        @cached_method
        def _c_bound(self, sign=1):
            # Every path is an integral combination of PARI's path generators.
            # Their values therefore bound all symbol denominators, without
            # a conjectural optimal-curve identification or period comparison.
            analytic.symbol(E, sign, 0)
            M, symbol = analytic._symbols[(E.a_invariants(), sign)]
            values = pari.mseval(M, symbol)
            return max([ZZ.zero()] + [QQ(v).denominator().valuation(p) for v in values])

    lp = ExactSymbols()
    P = analytic.gens(E)[0]
    r = 1
    for prec in (20, 40, 80):
        if E.has_multiplicative_reduction(p):
            reg = E.tate_curve(p).padic_height(prec=prec+4)(P)
        elif E.is_ordinary(p):
            reg = E.padic_height(p, prec=prec)(P)
        else:
            h = lp.Dp_valued_height(prec=prec)(P)
            phi = lp.frobenius(prec+2)
            a, c = phi[0, 0], phi[1, 0]
            reg = vector([h[0] - a/c*h[1], h[1]/c])
        if reg != 0:
            break
    else:
        raise NotImplementedError('p-adic regulator is zero to the computed precision')
    lg = Qp(p, prec)(1+p).log()
    factor = QQ(E.tamagawa_product()) / E.torsion_order()**2
    if E.is_ordinary(p):
        if E.has_good_reduction(p):
            bsdp = factor * reg * (1 - 1/lp.alpha(prec=prec))**2 / lg
        else:
            r += 1
            bsdp = factor * reg * E.tate_curve(p).L_invariant(prec=prec) / lg**r
        if bsdp == 0:
            raise NotImplementedError('p-adic leading-term denominator is not certified nonzero')
        n = max(2, bsdp.valuation())
        while lp._prec_bounds(n, r+1)[r] <= bsdp.valuation():
            n += 1
        while True:
            leading = lp.series(n, prec=r+1)[r]
            if leading != 0:
                quotient = leading / bsdp
                break
            n += 1
    else:
        bsdp = factor * reg / lg
        n = max(3, min(x.valuation() for x in bsdp if x != 0) + 2)
        while True:
            series = lp.Dp_valued_series(n, prec=2)
            quotients = [series[i][1]/bsdp[i] for i in (0, 1)
                         if bsdp[i] != 0 and series[i][1] != 0]
            if quotients:
                if len(quotients) == 2 and quotients[0] - quotients[1] != 0:
                    raise ArithmeticError('inconsistent supersingular p-adic BSD quotients')
                quotient = max(quotients, key=lambda x: x.precision_relative())
                break
            n += 1
    if quotient.precision_relative() <= 0:
        raise NotImplementedError('p-adic bound has insufficient relative precision')
    bound = quotient.valuation()
    if bound < 0:
        raise ArithmeticError('negative Iwasawa upper bound for the order of Sha')
    return ZZ(bound)
