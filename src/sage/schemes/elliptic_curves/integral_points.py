r"""
Integral points on elliptic curves over `\QQ`

The elliptic logarithm method [ST1994b]_ and David's lower bound
[Dav1995]_, Theorem 2.1, give a finite coefficient bound. Arb balls
enclose the analytic quantities, exact Gram--Schmidt lengths certify
lattice reduction, and interval enumeration has rational fallbacks.

Write `h(a/b)=\log\max(|a|,|b|)` for coprime integers `a,b`, and
`H(P)=\lim_{n\to\infty}h(x(2^nP))/4^n`. This is Sage's canonical
height, twice Silverman's `\widehat h_S(P)`. The corrected comparison
[Sil1990]_, Theorem 1.1, used by ``E.silverman_height_bound()`` is

.. MATH::

    |H(P)-h(x(P))|&=|2\widehat h_S(P)-h(x(P))|\\
       &\leq2\max\bigl(h(j)/24+\mu(E)+0.961,\mu(E)+1.07\bigr),\\
    \mu(E)&=\frac{\log|\Delta|}{12}
       +\frac{\log\max(1,|j|)}{12}
       +\frac{\log\max(1,|b_2/12|)}{2}+\frac{\log 2^*}{2},

where `2^*=2` if `b_2\ne0` and `2^*=1` otherwise; `0.961` is
Bremner's correction of Silverman's `0.973`. We evaluate an enlargement in
balls, using `3` in place of `1.922` and `2.14`, and `2` in place of
`2^*`, to obtain an integer upper bound. Exact doubling and
polarization enclose the height matrix; interval LDL decomposition
certifies `H(\sum m_iP_i+T)\geq\lambda M^2`, where
`M=\max|m_i|` and `\lambda>0`. Both height-comparison directions are
needed here and below. Sage's ``CPS_height_bound()`` only bounds `h(x)-H`,
so it cannot replace Silverman's two-sided bound.

Let `B_h` bound that absolute difference. For `s=x+b_2/12`, let `R_0`
bound the moduli of the roots of `f(s)=s^3-c_4s/48-c_6/864`, and put
`Z=\max(2R_0,|b_2/12|,1)`. A real point has `s\geq-R_0`, so
`s\leq Z` implies `|x|\leq Z+|b_2/12|`. Integral points in this
region satisfy `H(P)\leq\log\max(1,Z+|b_2/12|)+B_h`.
For `s>Z`, each root is at distance at least `s/2` and
`f(t)\geq t^3/8` for `t\geq s`. The principal signed logarithm has
`|u(P)|\leq\int_s^\infty\sqrt{2}t^{-3/2}dt
=2\sqrt{2}/\sqrt{s}`. Also `0<x\leq2s` and
`H(P)\leq\log x+B_h`. Thus
`|2u(P)|\leq K\exp(-H(P)/2)`, with
`K=8\exp(B_h/2)`, and hence at most `K\exp(-\lambda M^2/2)`.

The points `2P_i` lie on the identity component, and
`2u(P)=\sum m_i u(2P_i)+b\omega`, where `\omega` is the least positive
real period and the denominator of `b` divides the rational torsion
exponent `t`. Choose `|u(2P_i)|\leq\omega/2`;
then `|b|\leq1+rM/2`. For `M\geq1`, the integer coefficients of
`t\,2u(P)` have logarithmic heights at most `\log(t(r+2)M)`.

Here is the application of [Dav1995]_, pp. 9--10, Theorem 2.1.
David's `k` counts elliptic logarithms, not the extra additive
coordinate. We take `k=r+1`, `D=1`, `\beta_0=0`,
`u_i=u(2P_i)` for `1\leq i\leq r`, and `u_k=\omega`, a pole of
`\wp`, so `\gamma_k=O`.
The curve has `g_2=c_4/12`, `g_3=c_6/216`; its height is
`h_D=\max(1,h(1,g_2,g_3),h(j))`. David takes `\tau=\omega_2/\omega_1`
in the fundamental domain, so
`\mathcal A=|\omega_1|^2\operatorname{Im}\tau
=\operatorname{Im}(\overline{\omega_1}\omega_2)` is the covolume of
the period lattice, independent of the oriented basis; here it is
`\omega` times the least positive imaginary part of a period, as
computed in ``_Logarithms.enclosures``.
David's conditions (1)--(4), p. 10, are respectively

.. MATH::

    D\log B_D&\geq\log V_1,\qquad V_1\geq\cdots\geq V_k,\\
    \log B_D&\geq\max(eh_D,h(\beta_i):0\leq i\leq k),\\
    \log V_i&\geq\max(\widehat h(\gamma_i),h_D,
                           3\pi|u_i|^2/(D\mathcal A)),\\
    e\leq E_D&\leq\min_i
           \frac{e\sqrt{D\log V_i}\sqrt{\mathcal A}}{\sqrt{3\pi}|u_i|}.

If `L(v)\neq0`, the theorem gives

.. MATH::

    \log|L(v)|>-C_D D^{2k+2}(\log B_D+\log(DE_D))
       (\log E_D)^{-2k-1}
       (\log\log B_D+h_D+\log(DE_D))^{k+1}\prod_i\log V_i,

where David's constant is
`C_D=2.9\,10^{6k+6}4^{2k^2}(k+1)^{2k^2+9k+12.3}`.
Here `L(v)=t\,2u(P)`, which is nonzero because `P\neq O`.
Set `E_D=e`, `\log V_i=A`, and
`\log B_D=\beta=\log M+A+3h+\lceil\log(t(r+2))\rceil`, with
`h\geq h_D`. The code chooses `A\geq\max(1,h,3\pi\omega^2/
\mathcal A,8H(P_i))`. David's projective canonical height of `2P_i`
is `3H(2P_i)/2=6H(P_i)`; using `8H(P_i)` only enlarges this
requirement. The period has height zero. These choices satisfy (1)--(4):
in (4) the numerator includes `e`, so `A\geq3\pi|u_i|^2/\mathcal A`
suffices. In (2), `3h\geq eh_D` and the coefficient heights are
bounded above as stated.

Since `\log\beta\leq\beta`, David's right-hand magnitude is at most
`C_D F`, where `F=A^k(\beta+1)(\beta+h+1)^{k+1}`. Round David's
constant up to `C=3\,10^{6k+6}4^{2k^2}(k+1)^{2k^2+9k+13}` and
divide the form by `t` to obtain `\lambda M^2/2<CF+\log K+\log t`.
The stopping test checks the strict reverse inequality at a power of ten.
For `s=\log M\geq(k+1)\log10`, the logarithmic derivative of `F` is
at most `(k+2)/\beta<2`; adding the positive constant `\log K+\log t`
cannot increase it. The left side grows as `\exp(2s)`, so the test
stays true for every larger `M`. Combining with the bounded region gives
the initial coefficient bound.

For reduction, round normalized logarithms at an integer scale,
enclosing the rounding errors. Every nonzero lattice vector has
squared length at least the minimum exact squared Gram--Schmidt
length. Use the resulting lower bound for the linear form only after
certifying positivity, increasing scale and precision as needed.
In lattice dimensions at most eight, [FP1985]_ enumeration computes
the exact shortest-vector length. The final search intersects the
coefficient box with a certified height ellipsoid and includes all
torsion translates and both signs. These bounds apply to any independent
set of rational points: the list is complete in their integer span
together with the torsion subgroup. No full-rank or saturation
hypothesis is used. Automatic generators come from ``E.gens()``.

AUTHORS:

- Manvir Jaswal (2026): certified bounds and enumeration

EXAMPLES::

    sage: from sage.schemes.elliptic_curves.integral_points import integral_points
    sage: [P[0] for P in integral_points(EllipticCurve([1,-1,1,-11,16]))]
    [-1, 2, 20829482]
    sage: E = EllipticCurve('389a1')
    sage: [P[0] for P in integral_points(E, [E(0, 0)])]
    [0, 3]
"""

# *****************************************************************************
#       Copyright (C) 2026 Manvir Jaswal
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.
# *****************************************************************************

from math import isqrt

from sage.matrix.constructor import matrix
from sage.rings.complex_arb import ComplexBallField
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.qqbar import QQbar
from sage.rings.rational_field import QQ
from sage.rings.real_arb import RealBallField
from sage.rings.real_double import RDF


def _upper(x):
    """
    Return a rational upper endpoint of a finite real ball.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import _upper
        sage: _upper(RealBallField(53)('nan'))
        Traceback (most recent call last):
        ...
        ArithmeticError: a finite enclosure could not be computed
    """
    if not x.is_finite():
        raise ArithmeticError("a finite enclosure could not be computed")
    return x.upper().exact_rational()


def _ceil_upper(x):
    """
    Round a real ball upwards, including its error.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import _ceil_upper
        sage: _ceil_upper(RealBallField(53)(2).add_error(1/4))
        3
        sage: _ceil_upper(RealBallField(53)(-3/2))
        -1
    """
    return _upper(x).ceil()


def _ceil_sqrt_upper(x):
    r"""
    Bound the square root of the nonnegative part of an enclosure.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import _ceil_sqrt_upper
        sage: _ceil_sqrt_upper(RealBallField(53)(0).add_error(1/4))
        1
        sage: _ceil_sqrt_upper(RealBallField(53)(-1))
        0
    """
    # A ball straddling zero has no real ball square root; use its upper endpoint.
    return _ceil_upper(x.parent()(max(QQ(0), _upper(x))).sqrt())


def _validate_mw_base(E, mw_base):
    """
    Check that supplied generators are nontorsion points on this curve.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import _validate_mw_base
        sage: E = EllipticCurve('37a1')
        sage: _validate_mw_base(E, (E(0, 0),))
        [(0 : 0 : 1)]
    """
    try:
        len(mw_base)
    except TypeError:
        raise TypeError("mw_base must be a list")
    from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
    points = list(mw_base)
    if not all(isinstance(P, EllipticCurvePoint_field) and P.curve() is E
               for P in points):
        raise ValueError("points are not on the correct curve")
    if any(not P or P.has_finite_order() for P in points):
        raise ValueError("mw_base contains a torsion point")
    return points


def _height_difference(E, R):
    r"""
    Evaluate an enlargement of ``E.silverman_height_bound()`` in balls.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import _height_difference
        sage: E = EllipticCurve('37a1')
        sage: B = _height_difference(E, RealBallField(128)); B
        7
        sage: B >= E.silverman_height_bound()
        True
    """
    j = E.j_invariant()
    b = E.b2()/12
    h_j = R(max(abs(j.numerator()), j.denominator())).log()
    # Silverman's bound, enlarged to avoid its decimal constants.
    value = (h_j/12 + R(abs(E.discriminant())).log()/6
             + R(max(QQ(1), abs(j))).log()/6
             + R(max(QQ(1), abs(b))).log() + R(2).log() + 3)
    return _ceil_upper(value)


def _positive_definite(H):
    """
    Verify positive definiteness by interval LDL decomposition.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import _positive_definite
        sage: R = RealBallField(53)
        sage: _positive_definite([[R(2), R(1)], [R(1), R(2)]])
        True
        sage: _positive_definite([[R(1), R(2)], [R(2), R(1)]])
        False
    """
    n = len(H)
    L = [[0 for j in range(n)] for i in range(n)]
    diagonal = []
    for j in range(n):
        d = H[j][j] - sum(L[j][k]**2 * diagonal[k] for k in range(j))
        if not d > 0:
            return False
        diagonal.append(d)
        for i in range(j+1, n):
            L[i][j] = (H[i][j] - sum(L[i][k]*L[j][k]*diagonal[k]
                                     for k in range(j))) / d
    return True


def _height_matrix_bound(E, basis, B, R):
    r"""
    Enclose canonical heights and certify a positive lower eigenvalue.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import _height_difference, _height_matrix_bound
        sage: E = EllipticCurve('389a1'); basis = [E(-1, 1), E(0, 0)]; R = RealBallField(128)
        sage: lower, _, H = _height_matrix_bound(E, basis, _height_difference(E, R), R)
        sage: M = E.height_pairing_matrix(basis)
        sage: all(H[i][j].contains_exact(M[i, j]) for i in range(2) for j in range(2))
        True
        sage: 0 < lower < min(M.change_ring(RDF).eigenvalues())
        True
    """
    r = len(basis)
    points = list(basis)
    pairs = [(i, j) for i in range(r) for j in range(i)]
    points.extend(basis[i] + basis[j] for i, j in pairs)
    multiples = [8*P for P in points]
    for n in range(3, 11):
        scale = ZZ(4)**n
        heights = []
        for P in multiples:
            if not P:
                heights.append(R(0))
            else:
                x = P[0]
                naive = R(max(abs(x.numerator()), x.denominator())).log()
                heights.append((naive/scale).add_error(R(B)/scale))
        H = [[R(0) for j in range(r)] for i in range(r)]
        for i in range(r):
            H[i][i] = heights[i]
        for h, (i, j) in zip(heights[r:], pairs):
            H[i][j] = H[j][i] = (h - heights[i] - heights[j])/2
        if r == 1:
            lower = H[0][0].lower().exact_rational()
            if lower > 0 and lower > 3*H[0][0].mid()/4:
                return lower, heights[:r], H
        elif _positive_definite(H):
            mid = matrix(RDF, [[float(x.mid()) for x in row] for row in H])
            estimate = min(x.real() for x in mid.eigenvalues())
            uncertainty = max(sum(_upper(abs(x-R(x.mid()))) for x in row)
                              for row in H)
            if uncertainty > QQ(estimate)/32:
                multiples = [2*P for P in multiples]
                continue
            candidate = QQ(estimate)*QQ(7)/8
            for _ in range(16):
                if candidate > 0:
                    shifted = [[H[i][j] - (candidate if i == j else 0)
                                for j in range(r)] for i in range(r)]
                    if _positive_definite(shifted):
                        return candidate, heights[:r], H
                candidate /= 2
        multiples = [2*P for P in multiples]
    raise ArithmeticError("a positive height bound could not be certified")


class _Logarithms:
    """Enclose real logarithms and periods using Carlson's integral."""

    def __init__(self, E, basis):
        """
        Find the roots of `f(s)=s^3-c_4s/48-c_6/864`, sorting the real ones.

        TESTS::

            sage: from sage.schemes.elliptic_curves.integral_points import _Logarithms
            sage: [len(_Logarithms(EllipticCurve(l), []).real) for l in ['37a1', '11a1']]
            [3, 1]
        """
        self.E = E
        self.basis = basis
        f = PolynomialRing(QQ, 'x')([-E.c6()/864, -E.c4()/48, 0, 1])
        self.roots = f.roots(QQbar, multiplicities=False)
        self.real = sorted(z for z in self.roots if z.imag() == 0)
        self.root = self.real[-1]

    def enclosures(self, precision):
        r"""
        Return enclosures of `u(2P)/\omega` for each basis point `P`, of the
        real period `\omega` and of the covolume of the period lattice, and a
        rational upper bound for the moduli of the roots of `f`.

        TESTS:

        Compare with the period lattice for each sign of the discriminant::

            sage: from sage.schemes.elliptic_curves.integral_points import _Logarithms
            sage: for E, P in [(EllipticCurve('37a1'), (0, 0)), (EllipticCurve([0,0,0,0,-2]), (3, 5))]:
            ....:     P = E(P); L = E.period_lattice()
            ....:     logs, w, area, _ = _Logarithms(E, [P]).enclosures(100)
            ....:     t = (2*P).elliptic_logarithm().real()/w.mid() - logs[0].mid()
            ....:     assert abs(w.mid() - L.real_period()) < 1e-12
            ....:     assert abs(area.mid() - L.complex_area()) < 1e-12
            ....:     assert abs(t - t.round()) < 1e-12
        """
        E = self.E
        C = ComplexBallField(precision)
        R = RealBallField(precision)
        e1 = C(self.root)
        e2, e3 = [C(z) for z in self.roots if z != self.root]
        w = (2*C(0).elliptic_rf(e1-e2, e1-e3)).real()
        if E.discriminant() > 0:
            a, b, c = self.real
            v = (2*C(0).elliptic_rf(C(b-a), C(c-a))).real()
        else:
            v = C(0).elliptic_rf(e2-e1, e3-e1).real()
        if not (w > 0 and v > 0):
            raise ArithmeticError("period enclosures are not positive")
        logarithms = []
        for P in self.basis:
            Q = 2*P
            s = Q[0] + E.b2()/12
            if not R(s) - e1.real() > 0:
                raise ArithmeticError("the real logarithm is not separated from a branch point")
            u = C(C(s)-e1).elliptic_rf(C(s)-e2, C(s)-e3).real()
            # The derivative of wp is negative on (0, w/2).
            if Q[1] + (E.a1()*Q[0] + E.a3())/2 > 0:
                u = -u
            logarithms.append(u/w)
        root_bound = max(_upper(abs(C(z))) for z in self.roots)
        return logarithms, w, w*v, root_bound


def _shortest_vector_squared(basis, orthogonal=None):
    r"""
    Return the exact squared minimum of an integer row lattice of rank at most eight.

    Use Fincke--Pohst enumeration [FP1985]_ with rational Gram--Schmidt
    coefficients and outward integer bounds; no floating-point test prunes a vector.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import _shortest_vector_squared
        sage: _shortest_vector_squared(matrix(ZZ, [[2, 1], [1, 2]]))
        2
        sage: _shortest_vector_squared(matrix(ZZ, [[3]]))
        9
        sage: _shortest_vector_squared(identity_matrix(ZZ, 8))
        1
    """
    dimension = basis.nrows()
    if not 1 <= dimension <= 8:
        raise ValueError("the lattice rank must be between one and eight")
    if orthogonal is None:
        orthogonal = basis.gram_schmidt()[0]
    stars = orthogonal.rows()
    lengths = [row.dot_product(row) for row in stars]
    if any(d <= 0 for d in lengths):
        raise ArithmeticError("the lattice basis is not independent")
    rows = basis.rows()
    mu = [[rows[j].dot_product(stars[i])/lengths[i] for i in range(j)]
          for j in range(dimension)]
    best = min(row.dot_product(row) for row in rows)
    coefficients = [ZZ(0)] * dimension

    def visit(i, partial, nonzero):
        nonlocal best
        if i < 0:
            if nonzero:
                best = min(best, partial)
            return
        center = sum(mu[j][i]*coefficients[j] for j in range(i+1, dimension))
        radius = ZZ(isqrt(int(((best-partial)/lengths[i]).ceil())) + 1)
        first, last = (-center-radius).ceil(), (-center+radius).floor()
        if not nonzero:
            first = max(first, ZZ(0))
        for n in range(first, last+1):
            total = partial + lengths[i]*(n+center)**2
            if total < best:
                coefficients[i] = ZZ(n)
                visit(i-1, total, nonzero or n != 0)

    visit(dimension-1, QQ(0), False)
    return ZZ(best)


def _coefficient_bound(E, basis, verbose=False):
    r"""
    Return a certified bound for free coefficients of integral points.

    The points in ``basis`` must be independent modulo torsion; they need
    not generate the full Mordell--Weil group. Positive definiteness of
    their height matrix is certified or an error is raised.
    The subgroup includes all torsion translates. Return the coefficient
    bound, the enclosed height matrix and a rational height upper bound.
    The curve must be an integral model.

    TESTS:

    The point with x-coordinate 20829482 is a seventh multiple::

        sage: from sage.schemes.elliptic_curves.integral_points import _coefficient_bound
        sage: E = EllipticCurve([1,-1,1,-11,16])
        sage: P = E(2,-2)
        sage: (7*P)[0]
        20829482
        sage: _coefficient_bound(E, [P])[0] >= 7
        True
    """
    r = len(basis)
    if not r:
        return ZZ(0), [], QQ(0)
    R = RealBallField(128)
    B = _height_difference(E, R)
    eigenvalue, heights, height_matrix = _height_matrix_bound(E, basis, B, R)
    geometry = _Logarithms(E, basis)
    precision = 128
    while True:
        try:
            logs, w, area, root_bound = geometry.enclosures(precision)
            break
        except ArithmeticError:
            precision *= 2
            if precision > 65536:
                raise ArithmeticError("real logarithms could not be enclosed")
    R = w.parent()
    # Bounded region and K: see the module docstring.
    shift = abs(E.b2()/12)
    z = max(2*root_bound, shift, QQ(1))
    small_height = R(max(QQ(1), z+shift)).log()+B
    small_bound = _ceil_upper((small_height/eigenvalue).sqrt())
    height_limit = None
    log_upper = (R(4)*R(2).sqrt()).log() + (R(B)+R(2).log())/2
    # David's conditions (1)-(4): see the module docstring.
    g2, g3 = E.c4()/12, E.c6()/216
    d = g2.denominator().lcm(g3.denominator())
    h = max(ZZ(1), _ceil_upper(R(max(d, abs(d*g2), abs(d*g3))).log()),
            _ceil_upper(R(max(abs(E.j_invariant().numerator()),
                             E.j_invariant().denominator())).log()))
    A = max(ZZ(1), h, _ceil_upper(3*R.pi()*w**2/area),
            *(_ceil_upper(8*height) for height in heights))
    torsion_exponent = E.torsion_subgroup().exponent()
    k = r+1
    constant = (ZZ(3)*10**(6*k+6)*4**(2*k*k)
                * (k+1)**(2*k*k+9*k+13))
    log_t = _upper(R(torsion_exponent).log())
    offset = A + 3*h + _ceil_upper(R(torsion_exponent*(r+2)).log())
    power = ZZ(k+1)
    # Persistence of this stopping test is proved in the module docstring.
    while True:
        s = power*R(10).log()
        beta = s + offset
        rhs = constant*R(A)**k*(beta+1)*(beta+h+1)**(k+1)+log_upper+log_t
        if R(eigenvalue)*R(10)**(2*power)/2 > rhs:
            break
        power *= 2
    N = max(ZZ(10)**power, small_bound)
    if verbose:
        print("Height difference bound:", B)
        print("Initial coefficient bound:", N)

    for _ in range(32):
        if N <= small_bound:
            break
        # Start near the coefficient scale when the exact minimum is
        # available. A failed positivity test only enlarges the scale.
        scale = ((torsion_exponent*N)**(r+1) if r+1 <= 8
                 else (10*torsion_exponent*(r+1)*N)**(r+1))
        for _ in range(32):
            needed = scale.nbits()+64
            if needed > precision:
                precision = needed
                logs, w, area, root_bound = geometry.enclosures(precision)
            R = w.parent()
            rounded = [ZZ((scale*x).mid().round()) for x in logs]
            errors = [_upper(abs(scale*x-a)) for x, a in zip(logs, rounded)]
            lattice = matrix(ZZ, r+1, r+1, 1)
            for i, a in enumerate(rounded):
                lattice[i, r] = a
            lattice[r, r] = scale
            reduced = lattice.LLL()
            gs = reduced.gram_schmidt()[0]
            shortest_squared = (_shortest_vector_squared(reduced, gs) if r+1 <= 8
                                else min(row.dot_product(row) for row in gs.rows()))
            gap = QQ(shortest_squared)/torsion_exponent**2 - r*N**2
            if gap > 0:
                lower = (R(gap).sqrt()-N*sum(errors))/scale
                if lower > 0:
                    # Re-enclose the analytic upper constant at this precision.
                    log_K = ((R(4)*R(2).sqrt()).log()
                             + (R(B)+R(2).log())/2-w.log())
                    analytic_height = 2*(log_K-lower.log())
                    limit = max(_upper(small_height), _upper(analytic_height))
                    height_limit = limit if height_limit is None else min(height_limit, limit)
                    new_bound = max(small_bound, _ceil_sqrt_upper(analytic_height/eigenvalue))
                    break
            scale *= 16
        else:
            raise ArithmeticError("lattice reduction could not certify a bound")
        if new_bound >= N:
            break
        N = ZZ(new_bound)
    else:
        raise ArithmeticError("coefficient reduction did not terminate")
    if verbose:
        print("Coefficient bound:", N)
    if height_limit is None:
        height_limit = sum(_upper(abs(x))*N**2 for row in height_matrix for x in row)
    grid = ZZ(2)**64
    return ZZ(N), height_matrix, QQ((height_limit*grid).ceil())/grid


def _height_ellipsoid(height_matrix, height_limit):
    r"""
    Return an exact LDL ellipsoid containing the certified height ellipsoid.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import _height_ellipsoid
        sage: R = RealBallField(128)
        sage: _height_ellipsoid([[R(2), R(1)], [R(1), R(2)]], 10)
        ([3/2, 2], [[], [1/2]], 10)
        sage: _height_ellipsoid([[R(1), R(2)], [R(2), R(1)]], 10) is None
        True
    """
    r = len(height_matrix)
    grid = ZZ(2)**64
    mid = [[QQ((x.mid().exact_rational()*grid).round())/grid for x in row]
           for row in height_matrix]
    error = max(sum(_upper(abs(x-x.parent()(m))) for x, m in zip(row, centers))
                for row, centers in zip(height_matrix, mid))
    error = QQ((error*grid).ceil())/grid
    # A symmetric error matrix has operator norm at most its maximum
    # absolute row sum. Thus mid-error*I is below the true height matrix.
    A = [[mid[i][j] - (error if i == j else 0) for j in range(r)] for i in range(r)]
    L = [[QQ(0) for j in range(r)] for i in range(r)]
    diagonal = []
    for j in range(r):
        d = A[j][j] - sum(L[j][k]**2*diagonal[k] for k in range(j))
        if d <= 0:
            return None
        diagonal.append(d)
        for i in range(j+1, r):
            L[i][j] = (A[i][j] - sum(L[i][k]*L[j][k]*diagonal[k] for k in range(j)))/d
    mu = [[L[r-1-k][r-1-j] for k in range(j)] for j in range(r)]
    return list(reversed(diagonal)), mu, QQ(height_limit)


def _bounded_integral_x_coordinates(E, basis, N, height_data=None, precision=128):
    r"""
    Enumerate integral x-coordinates with coefficients bounded by ``N``.

    With ``height_data=(matrix, bound)``, intersect the box with its certified
    height ellipsoid. Include all torsion translates, with exact fallbacks.
    The curve must be an integral model, where integral x implies integral y.
    Use ``precision`` bits for the real enclosures.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import _bounded_integral_x_coordinates
        sage: E = EllipticCurve('389a1')
        sage: sorted(_bounded_integral_x_coordinates(E, [E(-1, 1), E(0, 0)], 4))
        [-2, -1, 0, 1, 3, 4, 6, 39, 133, 188]
        sage: R = RealBallField(128)
        sage: H = [[R(1), R(0)], [R(0), R(1)]]
        sage: sorted(_bounded_integral_x_coordinates(E, [E(0, 0), E(-1, 1)], 4, (H, 1)))
        [-1, 0]
        sage: E = EllipticCurve('11a1')
        sage: sorted(_bounded_integral_x_coordinates(E, [], 0))
        [5, 16]
        sage: _bounded_integral_x_coordinates(E, [], -1)
        Traceback (most recent call last):
        ...
        ValueError: the coefficient bound must be nonnegative

    Low-precision enclosures fall back to exact group operations::

        sage: E = EllipticCurve([0,0,0,-25,0]); P = E(-4,6)
        sage: xs = _bounded_integral_x_coordinates(E, [P,P], 2)
        sage: sorted(xs)
        [-5, -4, 0, 5, 45]
        sage: _bounded_integral_x_coordinates(E, [P,P], 2, precision=2) == xs
        True
    """
    N = ZZ(N)
    if N < 0:
        raise ValueError("the coefficient bound must be nonnegative")
    torsion = E.torsion_points()
    zero = E(0)
    xs = set()

    def use(P):
        if P and P[0].is_integral():
            xs.add(ZZ(P[0]))

    r = len(basis)
    ellipsoid = _height_ellipsoid(*height_data) if height_data is not None and r else None
    if ellipsoid is not None:
        diagonal, mu, height_limit = ellipsoid
        basis = list(reversed(basis))
        if r == 1:
            N = min(N, ZZ(isqrt(int((height_limit/diagonal[0]).floor()))))
    if r <= 1:
        P = zero
        for n in range(N+1 if r else 1):
            for T in torsion:
                use(P+T)
            if r:
                P += basis[0]
        return xs

    R = RealBallField(precision)
    a1, a2, a3 = map(R, E.a_invariants()[:3])
    has_a1, has_a2 = bool(E.a1()), bool(E.a2())

    def enclose(P):
        return (R(P[0]), R(P[1])) if P else None

    def add(P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        denominator = x2-x1
        if denominator.contains_zero():
            raise ArithmeticError("ambiguous ball addition")
        slope = (y2-y1)/denominator
        x3 = slope**2-x1-x2
        if has_a1:
            x3 += a1*slope
        if has_a2:
            x3 -= a2
        y3 = -y1-a1*x3-a3+slope*(x1-x3)
        return x3, y3

    multiples = []
    for P in basis:
        row = {0: zero}
        Q = zero
        for n in range(1, N+1):
            Q += P
            row[n], row[-n] = Q, -Q
        multiples.append({n: (Q, enclose(Q)) for n, Q in row.items()})
    # Absorb torsion into the final summand. At the last level only its
    # x-coordinate is needed, so no final y-coordinate enclosure is formed.
    last_points = {}
    for n, (Q, ball) in multiples[-1].items():
        translates = [Q+T if T else Q for T in torsion]
        last_points[n] = [(P, enclose(P)) for P in translates]

    def choices(j, negative, coefficients, partial):
        last = N if negative else ZZ(0)
        if ellipsoid is None:
            for n in range(-N, last+1):
                yield n, partial
            return
        center = sum(a*n for a, n in zip(mu[j], coefficients))
        radius = ZZ(isqrt(int(((height_limit-partial)/diagonal[j]).ceil())) + 1)
        first = max(-N, (-center-radius).ceil())
        last = min(last, (-center+radius).floor())
        for n in range(first, last+1):
            total = partial + diagonal[j]*(n+center)**2
            if total <= height_limit:
                yield n, total

    def visit(j, enclosure, terms, negative, coefficients, partial):
        if j < r-1:
            for n, total in choices(j, negative, coefficients, partial):
                Q, ball = multiples[j][n]
                new_terms = terms + (Q,)
                try:
                    new_ball = add(enclosure, ball)
                except ArithmeticError:
                    new_ball = enclose(sum(new_terms, zero))
                visit(j+1, new_ball, new_terms, negative or n < 0,
                      coefficients + (n,), total)
            return
        exact = None
        for n, total in choices(j, negative, coefficients, partial):
            for Q, ball in last_points[n]:
                if enclosure is None:
                    if ball is None:
                        continue
                    x = ball[0]
                elif ball is None:
                    x = enclosure[0]
                else:
                    x1, y1 = enclosure
                    x2, y2 = ball
                    denominator = x2-x1
                    if denominator.contains_zero():
                        x = None
                    else:
                        slope = (y2-y1)/denominator
                        x = slope**2-x1-x2
                        if has_a1:
                            x += a1*slope
                        if has_a2:
                            x -= a2
                # An ambiguous or non-finite enclosure must be checked
                # exactly. Cache the exact prefix across final summands.
                if x is None or not x.is_finite() or x.contains_integer():
                    if exact is None:
                        exact = sum(terms, zero)
                    use(exact+Q)

    visit(0, None, (), False, (), QQ(0))
    return xs


def integral_points(E, mw_base='auto', both_signs=False, verbose=False):
    r"""
    Compute integral points in the span of generators and rational torsion.

    See :meth:`~sage.schemes.elliptic_curves.ell_rational_field.EllipticCurve_rational_field.integral_points`
    for the arguments and completeness condition.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import integral_points
        sage: E = EllipticCurve([0,0,0,-16,0])
        sage: [P[0] for P in integral_points(E)]
        [-4, 0, 4]
        sage: points = integral_points(E, both_signs=True)
        sage: all(-P in points for P in points)
        True

    Curves constructed directly from the class are supported::

        sage: from sage.schemes.elliptic_curves.ell_rational_field import EllipticCurve_rational_field
        sage: E = EllipticCurve_rational_field([0,1,1,-2,0])
        sage: [P[0] for P in integral_points(E)]
        [-2, -1, 0, 1, 3, 4, 6, 39, 133, 188]
    """
    if not E.is_integral():
        raise ValueError("integral_points() can only be called on an integral model")
    if mw_base == 'auto':
        try:
            basis = E.gens()
        except RuntimeError:
            raise RuntimeError("Unable to compute Mordell-Weil basis of {}, hence unable to compute integral points.".format(E))
    else:
        basis = _validate_mw_base(E, mw_base)
    if basis:
        M = E.height_pairing_matrix(basis)
        R = RealBallField(M.base_ring().precision())
        # A positive height matrix is needed at the working precision.
        if not _positive_definite([[R(x) for x in row] for row in M]):
            raise ArithmeticError("the given points are dependent, or the precision is insufficient")
        basis, _ = E.lll_reduce(basis, M)
    if verbose:
        print("Generators:", basis)
    N, height_matrix, height_limit = _coefficient_bound(E, basis, verbose)
    xs = _bounded_integral_x_coordinates(E, basis, N, (height_matrix, height_limit))
    if both_signs:
        result = [P for x in xs for P in E.lift_x(x, all=True)]
    else:
        result = [E.lift_x(x) for x in xs]
    if verbose:
        print("Found %s integral points." % len(result))
    return sorted(result)
