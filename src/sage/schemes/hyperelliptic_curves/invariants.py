# -*- coding: utf-8 -*-
r"""
Compute invariants of quintics and sextics via 'Ueberschiebung'

.. TODO::

    * Implement invariants in small positive characteristic.

    * Cardona-Quer and additional invariants for classifying automorphism groups.

AUTHOR:

- Nick Alexander

"""
from sage.rings.integer_ring import ZZ
from sage.rings.all import PolynomialRing


def diffxy(f, x, xtimes, y, ytimes):
    r"""
    Differentiate a polynomial ``f``, ``xtimes`` with respect to ``x``, and
    ```ytimes`` with respect to ``y``.

    EXAMPLES::

        sage: R.<u, v> = QQ[]
        sage: sage.schemes.hyperelliptic_curves.invariants.diffxy(u^2*v^3, u, 0, v, 0)
        u^2*v^3
        sage: sage.schemes.hyperelliptic_curves.invariants.diffxy(u^2*v^3, u, 2, v, 1)
        6*v^2
        sage: sage.schemes.hyperelliptic_curves.invariants.diffxy(u^2*v^3, u, 2, v, 2)
        12*v
        sage: sage.schemes.hyperelliptic_curves.invariants.diffxy(u^2*v^3 + u^4*v^4, u, 2, v, 2)
        144*u^2*v^2 + 12*v
    """
    h = f
    for i in range(xtimes):
        h = h.derivative(x)
    for j in range(ytimes):
        h = h.derivative(y)
    return h


def differential_operator(f, g, k):
    r"""
    Return the differential operator `(f g)_k` symbolically in the polynomial ring in ``dfdx, dfdy, dgdx, dgdy``.

    This is defined by Mestre on p 315 [Mes1991]_:

    .. MATH::

        (f g)_k = \frac{(m - k)! (n - k)!}{m! n!} \left(
        \frac{\partial f}{\partial x} \frac{\partial g}{\partial y} -
        \frac{\partial f}{\partial y} \frac{\partial g}{\partial x} \right)^k .

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves.invariants import differential_operator
        sage: R.<x, y> = QQ[]
        sage: differential_operator(x, y, 0)
        1
        sage: differential_operator(x, y, 1)
        -dfdy*dgdx + dfdx*dgdy
        sage: differential_operator(x*y, x*y, 2)
        1/4*dfdy^2*dgdx^2 - 1/2*dfdx*dfdy*dgdx*dgdy + 1/4*dfdx^2*dgdy^2
        sage: differential_operator(x^2*y, x*y^2, 2)
        1/36*dfdy^2*dgdx^2 - 1/18*dfdx*dfdy*dgdx*dgdy + 1/36*dfdx^2*dgdy^2
        sage: differential_operator(x^2*y, x*y^2, 4)
        1/576*dfdy^4*dgdx^4 - 1/144*dfdx*dfdy^3*dgdx^3*dgdy + 1/96*dfdx^2*dfdy^2*dgdx^2*dgdy^2 - 1/144*dfdx^3*dfdy*dgdx*dgdy^3 + 1/576*dfdx^4*dgdy^4
    """
    (x, y) = f.parent().gens()
    n = max(ZZ(f.degree()), ZZ(k))
    m = max(ZZ(g.degree()), ZZ(k))
    R, (fx, fy, gx, gy) = PolynomialRing(f.base_ring(), 4, 'dfdx,dfdy,dgdx,dgdy').objgens()
    const = (m - k).factorial() * (n - k).factorial() / (m.factorial() * n.factorial())
    U = f.base_ring()(const) * (fx*gy - fy*gx)**k
    return U


def diffsymb(U, f, g):
    r"""
    Given a differential operator ``U`` in ``dfdx, dfdy, dgdx, dgdy``,
    represented symbolically by ``U``, apply it to ``f, g``.

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves.invariants import diffsymb
        sage: R.<x, y> = QQ[]
        sage: S.<dfdx, dfdy, dgdx, dgdy> = QQ[]
        sage: [ diffsymb(dd, x^2, y*0 + 1) for dd in S.gens() ]
        [2*x, 0, 0, 0]
        sage: [ diffsymb(dd, x*0 + 1, y^2) for dd in S.gens() ]
        [0, 0, 0, 2*y]
        sage: [ diffsymb(dd, x^2, y^2) for dd in S.gens() ]
        [2*x*y^2, 0, 0, 2*x^2*y]

        sage: diffsymb(dfdx + dfdy*dgdy, y*x^2, y^3)
        2*x*y^4 + 3*x^2*y^2
    """
    (x, y) = f.parent().gens()
    R, (fx, fy, gx, gy) = PolynomialRing(f.base_ring(), 4, 'dfdx,dfdy,dgdx,dgdy').objgens()
    res = 0
    for coeff, mon in list(U):
        mon = R(mon)
        a = diffxy(f, x, mon.degree(fx), y, mon.degree(fy))
        b = diffxy(g, x, mon.degree(gx), y, mon.degree(gy))
        temp = coeff * a * b
        res = res + temp
    return res


def Ueberschiebung(f, g, k):
    r"""
    Return the differential operator `(f g)_k`.

    This is defined by Mestre on page 315 [Mes1991]_:

    .. MATH::

        (f g)_k = \frac{(m - k)! (n - k)!}{m! n!} \left(
        \frac{\partial f}{\partial x} \frac{\partial g}{\partial y} -
        \frac{\partial f}{\partial y} \frac{\partial g}{\partial x} \right)^k .

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves.invariants import Ueberschiebung as ub
        sage: R.<x, y> = QQ[]
        sage: ub(x, y, 0)
        x*y
        sage: ub(x^5 + 1, x^5 + 1, 1)
        0
        sage: ub(x^5 + 5*x + 1, x^5 + 5*x + 1, 0)
        x^10 + 10*x^6 + 2*x^5 + 25*x^2 + 10*x + 1
    """
    U = differential_operator(f, g, k)
    # U is the (f g)_k = ... of Mestre, p315, symbolically
    return diffsymb(U, f, g)


def ubs(f):
    r"""
    Given a sextic form `f`, return a dictionary of the invariants of Mestre, p 317 [Mes1991]_.

    `f` may be homogeneous in two variables or inhomogeneous in one.

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves.invariants import ubs
        sage: x = QQ['x'].0
        sage: ubs(x^6 + 1)
        {'A': 2,
         'B': 2/3,
         'C': -2/9,
         'D': 0,
         'Delta': -2/3*x^2*h^2,
         'f': x^6 + h^6,
         'i': 2*x^2*h^2,
         'y1': 0,
         'y2': 0,
         'y3': 0}

        sage: R.<u, v> = QQ[]
        sage: ubs(u^6 + v^6)
        {'A': 2,
         'B': 2/3,
         'C': -2/9,
         'D': 0,
         'Delta': -2/3*u^2*v^2,
         'f': u^6 + v^6,
         'i': 2*u^2*v^2,
         'y1': 0,
         'y2': 0,
         'y3': 0}

        sage: R.<t> = GF(31)[]
        sage: ubs(t^6 + 2*t^5 + t^2 + 3*t + 1)
        {'A': 0,
         'B': -12,
         'C': -15,
         'D': -15,
         'Delta': -10*t^4 + 12*t^3*h + 7*t^2*h^2 - 5*t*h^3 + 2*h^4,
         'f': t^6 + 2*t^5*h + t^2*h^4 + 3*t*h^5 + h^6,
         'i': -4*t^4 + 10*t^3*h + 2*t^2*h^2 - 9*t*h^3 - 7*h^4,
         'y1': 4*t^2 - 10*t*h - 13*h^2,
         'y2': 6*t^2 - 4*t*h + 2*h^2,
         'y3': 4*t^2 - 4*t*h - 9*h^2}
    """
    ub = Ueberschiebung
    if f.parent().ngens() == 1:
        f = PolynomialRing(f.parent().base_ring(), 1, f.parent().variable_name())(f)
        x1, x2 = f.homogenize().parent().gens()
        f = sum([ f[i]*x1**i*x2**(6-i) for i in range(7) ])
    U = {}
    U['f'] = f
    U['i'] = ub(f, f, 4)
    U['Delta'] = ub(U['i'], U['i'], 2)
    U['y1'] = ub(f, U['i'], 4)
    U['y2'] = ub(U['i'], U['y1'], 2)
    U['y3'] = ub(U['i'], U['y2'], 2)
    U['A'] = ub(f, f, 6)
    U['B'] = ub(U['i'], U['i'], 4)
    U['C'] = ub(U['i'], U['Delta'], 4)
    U['D'] = ub(U['y3'], U['y1'], 2)
    return U


def clebsch_to_igusa(A, B, C, D):
    r"""
    Convert Clebsch invariants `A, B, C, D` to Igusa invariants `I_2, I_4, I_6, I_{10}`.

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves.invariants import clebsch_to_igusa, igusa_to_clebsch
        sage: clebsch_to_igusa(2, 3, 4, 5)
        (-240, 17370, 231120, -103098906)
        sage: igusa_to_clebsch(*clebsch_to_igusa(2, 3, 4, 5))
        (2, 3, 4, 5)

        sage: Cs = tuple(map(GF(31), (2, 3, 4, 5))); Cs
        (2, 3, 4, 5)
        sage: clebsch_to_igusa(*Cs)
        (8, 10, 15, 26)
        sage: igusa_to_clebsch(*clebsch_to_igusa(*Cs))
        (2, 3, 4, 5)
    """
    I2 = -120*A
    I4 = -720*A**2 + 6750*B
    I6 = 8640*A**3 - 108000*A*B + 202500*C
    I10 = -62208*A**5 + 972000*A**3*B + 1620000*A**2*C - 3037500*A*B**2 - 6075000*B*C - 4556250*D
    return (I2, I4, I6, I10)


def igusa_to_clebsch(I2, I4, I6, I10):
    r"""
    Convert Igusa invariants `I_2, I_4, I_6, I_{10}` to Clebsch invariants `A, B, C, D`.

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves.invariants import clebsch_to_igusa, igusa_to_clebsch
        sage: igusa_to_clebsch(-2400, 173700, 23112000, -10309890600)
        (20, 342/5, 2512/5, 43381012/1125)
        sage: clebsch_to_igusa(*igusa_to_clebsch(-2400, 173700, 23112000, -10309890600))
        (-2400, 173700, 23112000, -10309890600)

        sage: Is = tuple(map(GF(31), (-2400, 173700, 23112000, -10309890600))); Is
        (18, 7, 12, 27)
        sage: igusa_to_clebsch(*Is)
        (20, 25, 25, 12)
        sage: clebsch_to_igusa(*igusa_to_clebsch(*Is))
        (18, 7, 12, 27)
    """
    A = -(+ I2) / 120
    B = -(- I2**2 - 20*I4)/135000
    C = -(+ I2**3 + 80*I2*I4 - 600*I6)/121500000
    D = -(+ 9*I2**5 + 700*I2**3*I4 - 3600*I2**2*I6 - 12400*I2*I4**2 + 48000*I4*I6 + 10800000*I10) / 49207500000000
    return (A, B, C, D)


def clebsch_invariants(f):
    r"""
    Given a sextic form `f`, return the Clebsch invariants `(A, B, C, D)` of
    Mestre, p 317, [Mes1991]_.

    `f` may be homogeneous in two variables or inhomogeneous in one.

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves.invariants import clebsch_invariants
        sage: R.<x, y> = QQ[]
        sage: clebsch_invariants(x^6 + y^6)
        (2, 2/3, -2/9, 0)
        sage: R.<x> = QQ[]
        sage: clebsch_invariants(x^6 + x^5 + x^4 + x^2 + 2)
        (62/15, 15434/5625, -236951/140625, 229930748/791015625)

        sage: magma(x^6 + 1).ClebschInvariants() # optional - magma
        [ 2, 2/3, -2/9, 0 ]
        sage: magma(x^6 + x^5 + x^4 + x^2 + 2).ClebschInvariants() # optional - magma
        [ 62/15, 15434/5625, -236951/140625, 229930748/791015625 ]
    """
    R = f.parent().base_ring()
    if R.characteristic() in [2, 3, 5]:
        raise NotImplementedError("Invariants of binary sextics/genus 2 hyperelliptic "
                                  "curves not implemented in characteristics 2, 3, and 5")

    U = ubs(f)
    L = U['A'], U['B'], U['C'], U['D']
    assert all(t.is_constant() for t in L)
    return tuple([ t.constant_coefficient() for t in L ])


def igusa_clebsch_invariants(f):
    r"""
    Given a sextic form `f`, return the Igusa-Clebsch invariants `I_2, I_4,
    I_6, I_{10}` of Igusa and Clebsch [IJ1960]_.

    `f` may be homogeneous in two variables or inhomogeneous in one.

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves.invariants import igusa_clebsch_invariants
        sage: R.<x, y> = QQ[]
        sage: igusa_clebsch_invariants(x^6 + y^6)
        (-240, 1620, -119880, -46656)
        sage: R.<x> = QQ[]
        sage: igusa_clebsch_invariants(x^6 + x^5 + x^4 + x^2 + 2)
        (-496, 6220, -955932, -1111784)

        sage: magma(x^6 + 1).IgusaClebschInvariants() # optional - magma
        [ -240, 1620, -119880, -46656 ]
        sage: magma(x^6 + x^5 + x^4 + x^2 + 2).IgusaClebschInvariants() # optional - magma
        [ -496, 6220, -955932, -1111784 ]

    TESTS:

    Let's check a symbolic example::

        sage: R.<a, b, c, d, e> = QQ[]
        sage: S.<x> = R[]
        sage: igusa_clebsch_invariants(x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e)[0]
        6*b^2 - 16*a*c + 40*d

        sage: absolute_igusa_invariants_wamelen(GF(5)['x'](x^6 - 2*x))
        Traceback (most recent call last):
        ...
        NotImplementedError: Invariants of binary sextics/genus 2 hyperelliptic curves not implemented in characteristics 2, 3, and 5
    """
    return clebsch_to_igusa(*clebsch_invariants(f))


def absolute_igusa_invariants_wamelen(f):
    r"""
    Given a sextic form `f`, return the three absolute Igusa invariants used by van Wamelen [Wam1999]_.

    `f` may be homogeneous in two variables or inhomogeneous in one.

    REFERENCES:

    - [Wam1999]_

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves.invariants import absolute_igusa_invariants_wamelen
        sage: R.<x> = QQ[]
        sage: absolute_igusa_invariants_wamelen(x^5 - 1)
        (0, 0, 0)

    The following example can be checked against van Wamelen's paper::

        sage: i1, i2, i3 = absolute_igusa_invariants_wamelen(-x^5 + 3*x^4 + 2*x^3 - 6*x^2 - 3*x + 1)
        sage: list(map(factor, (i1, i2, i3)))
        [2^7 * 3^15, 2^5 * 3^11 * 5, 2^4 * 3^9 * 31]

    TESTS::

        sage: absolute_igusa_invariants_wamelen(GF(3)['x'](x^5 - 2*x))
        Traceback (most recent call last):
        ...
        NotImplementedError: Invariants of binary sextics/genus 2 hyperelliptic curves not implemented in characteristics 2, 3, and 5
    """
    I2, I4, I6, I10 = igusa_clebsch_invariants(f)
    i1 = I2**5/I10
    i2 = I2**3*I4/I10
    i3 = I2**2*I6/I10
    return (i1, i2, i3)


def absolute_igusa_invariants_kohel(f):
    r"""
    Given a sextic form `f`, return the three absolute Igusa invariants used by Kohel [KohECHIDNA]_.

    `f` may be homogeneous in two variables or inhomogeneous in one.

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves.invariants import absolute_igusa_invariants_kohel
        sage: R.<x> = QQ[]
        sage: absolute_igusa_invariants_kohel(x^5 - 1)
        (0, 0, 0)
        sage: absolute_igusa_invariants_kohel(x^5 - x)
        (100, -20000, -2000)

    The following example can be checked against Kohel's database [KohECHIDNA]_ ::

        sage: i1, i2, i3 = absolute_igusa_invariants_kohel(-x^5 + 3*x^4 + 2*x^3 - 6*x^2 - 3*x + 1)
        sage: list(map(factor, (i1, i2, i3)))
        [2^2 * 3^5 * 5 * 31, 2^5 * 3^11 * 5, 2^4 * 3^9 * 31]
        sage: list(map(factor, (150660, 28343520, 9762768)))
        [2^2 * 3^5 * 5 * 31, 2^5 * 3^11 * 5, 2^4 * 3^9 * 31]

    TESTS::

        sage: absolute_igusa_invariants_kohel(GF(2)['x'](x^5 - x))
        Traceback (most recent call last):
        ...
        NotImplementedError: Invariants of binary sextics/genus 2 hyperelliptic curves not implemented in characteristics 2, 3, and 5
    """
    I2, I4, I6, I10 = igusa_clebsch_invariants(f)
    i1 = I4*I6/I10
    i2 = I2**3*I4/I10
    i3 = I2**2*I6/I10
    return (i1, i2, i3)
