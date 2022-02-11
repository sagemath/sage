# -*- coding: utf-8 -*-
r"""
Weierstrass `\wp`-function for elliptic curves

The Weierstrass `\wp` function associated to an elliptic curve over a
field `k` is a Laurent series of the form

.. MATH::

    \wp(z) = \frac{1}{z^2} +  c_2 \cdot z^2 + c_4 \cdot z^4 + \cdots.

If the field is contained in `\mathbb{C}`, then this is the series
expansion of the map from `\mathbb{C}` to `E(\mathbb{C})` whose kernel
is the period lattice of `E`.

Over other fields, like finite fields, this still makes sense as a
formal power series with coefficients in `k` - at least its first `p-2`
coefficients where `p` is the characteristic of `k`. It can be defined
via the formal group as `x+c` in the variable `z=\log_E(t)` for a
constant `c` such that the constant term `c_0` in `\wp(z)` is zero.

EXAMPLES::

    sage: E = EllipticCurve([0,1])
    sage: E.weierstrass_p()
    z^-2 - 1/7*z^4 + 1/637*z^10 - 1/84721*z^16 + O(z^20)

REFERENCES:

- [BMSS2006]_

AUTHORS:

- Dan Shumov 04/09: original implementation

- Chris Wuthrich 11/09: major restructuring

- Jeroen Demeyer (2014-03-06): code clean up, fix characteristic bound
  for quadratic algorithm (see :trac:`15855`)

"""

#*****************************************************************************
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.power_series_ring import PowerSeriesRing

def weierstrass_p(E, prec=20, algorithm=None):
    r"""
    Computes the Weierstrass `\wp`-function on an elliptic curve.

    INPUT:

    - ``E`` -- an elliptic curve

    - ``prec`` -- precision

    - ``algorithm`` -- string (default:``None``) an algorithm identifier
      indicating the ``pari``, ``fast`` or ``quadratic`` algorithm.
      If the algorithm is ``None``, then this function determines the
      best algorithm to use.

    OUTPUT:

    a Laurent series in one variable `z` with coefficients in the base
    field `k` of `E`.

    EXAMPLES::

        sage: E = EllipticCurve('11a1')
        sage: E.weierstrass_p(prec=10)
        z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8 + O(z^10)
        sage: E.weierstrass_p(prec=8)
        z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + O(z^8)
        sage: Esh = E.short_weierstrass_model()
        sage: Esh.weierstrass_p(prec=8)
        z^-2 + 13392/5*z^2 + 1080432/7*z^4 + 59781888/25*z^6 + O(z^8)

        sage: E.weierstrass_p(prec=8, algorithm='pari')
        z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + O(z^8)
        sage: E.weierstrass_p(prec=8, algorithm='quadratic')
        z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + O(z^8)

        sage: k = GF(11)
        sage: E = EllipticCurve(k, [1,1])
        sage: E.weierstrass_p(prec=6, algorithm='fast')
        z^-2 + 2*z^2 + 3*z^4 + O(z^6)
        sage: E.weierstrass_p(prec=7, algorithm='fast')
        Traceback (most recent call last):
        ...
        ValueError: for computing the Weierstrass p-function via the fast algorithm, the characteristic (11) of the underlying field must be greater than prec + 4 = 11
        sage: E.weierstrass_p(prec=8)
        z^-2 + 2*z^2 + 3*z^4 + 5*z^6 + O(z^8)
        sage: E.weierstrass_p(prec=8, algorithm='quadratic')
        z^-2 + 2*z^2 + 3*z^4 + 5*z^6 + O(z^8)
        sage: E.weierstrass_p(prec=8, algorithm='pari')
        z^-2 + 2*z^2 + 3*z^4 + 5*z^6 + O(z^8)
        sage: E.weierstrass_p(prec=9)
        Traceback (most recent call last):
        ...
        NotImplementedError: currently no algorithms for computing the Weierstrass p-function for that characteristic / precision pair is implemented. Lower the precision below char(k) - 2
        sage: E.weierstrass_p(prec=9, algorithm="quadratic")
        Traceback (most recent call last):
        ...
        ValueError: for computing the Weierstrass p-function via the quadratic algorithm, the characteristic (11) of the underlying field must be greater than prec + 2 = 11
        sage: E.weierstrass_p(prec=9, algorithm='pari')
        Traceback (most recent call last):
        ...
        ValueError: for computing the Weierstrass p-function via pari, the characteristic (11) of the underlying field must be greater than prec + 2 = 11

    TESTS::

        sage: E.weierstrass_p(prec=4, algorithm='foo')
        Traceback (most recent call last):
        ...
        ValueError: unknown algorithm for computing the Weierstrass p-function
    """
    k = E.base_ring()
    p = k.characteristic()

    # if the algorithm is not set, try to determine algorithm from input
    if algorithm is None:
        if p == 0 or p > prec + 4:
            algorithm = "fast"
        elif p > prec + 2:
            algorithm = "pari"
        else:
            raise NotImplementedError("currently no algorithms for computing the Weierstrass p-function for that characteristic / precision pair is implemented. Lower the precision below char(k) - 2")

    if algorithm == "pari":
        if 0 < p <= prec + 2:
            raise ValueError("for computing the Weierstrass p-function via pari, the characteristic (%s) of the underlying field must be greater than prec + 2 = %s"%(p,prec+2))
        return compute_wp_pari(E, prec)

    # quadratic and fast algorithms require short Weierstrass model
    Esh = E.short_weierstrass_model()

    A = Esh.a4()
    B = Esh.a6()

    if algorithm == "quadratic":
        if 0 < p <= prec + 2:
            raise ValueError("for computing the Weierstrass p-function via the quadratic algorithm, the characteristic (%s) of the underlying field must be greater than prec + 2 = %s"%(p,prec+2))
        wp = compute_wp_quadratic(k, A, B, prec)
    elif algorithm == "fast":
        if 0 < p <= prec + 4:
            raise ValueError("for computing the Weierstrass p-function via the fast algorithm, the characteristic (%s) of the underlying field must be greater than prec + 4 = %s"%(p,prec+4))
        wp = compute_wp_fast(k, A, B, prec)
    else:
        raise ValueError("unknown algorithm for computing the Weierstrass p-function")

    R = wp.parent()
    z = R.gen()
    u = E.isomorphism_to(Esh).u
    return wp(z*u) * u**2

def compute_wp_pari(E,prec):
    r"""
    Computes the Weierstrass `\wp`-function with the ``ellwp`` function
    from PARI.

    EXAMPLES::

        sage: E = EllipticCurve([0,1])
        sage: from sage.schemes.elliptic_curves.ell_wp import compute_wp_pari
        sage: compute_wp_pari(E, prec=20)
        z^-2 - 1/7*z^4 + 1/637*z^10 - 1/84721*z^16 + O(z^20)
        sage: compute_wp_pari(E, prec=30)
        z^-2 - 1/7*z^4 + 1/637*z^10 - 1/84721*z^16 + 3/38548055*z^22 - 4/8364927935*z^28 + O(z^30)
    """
    ep = E.__pari__()
    wpp = ep.ellwp(n=prec)
    k = E.base_ring()
    R = LaurentSeriesRing(k,'z')
    z = R.gen()
    wp = z**(-2)
    for i in range(prec):
        wp += k(wpp[i]) * z**i
    wp = wp.add_bigoh(prec)
    return wp


def compute_wp_quadratic(k, A, B, prec):
    r"""
    Compute the truncated Weierstrass function of an elliptic curve
    defined by short Weierstrass model: `y^2 = x^3 + Ax + B`. Uses an
    algorithm that is of complexity `O(prec^2)`.

    Let p be the characteristic of the underlying field. Then we must
    have either p = 0, or p > prec + 2.

    INPUT:

     - ``k`` - the field of definition of the curve
     - ``A`` - and
     - ``B`` - the coefficients of the elliptic curve
     - ``prec`` - the precision to which we compute the series.

    OUTPUT:

    A Laurent series approximating the Weierstrass `\wp`-function to precision ``prec``.

    ALGORITHM:

    This function uses the algorithm described in section 3.2 of [BMSS2006]_.

    EXAMPLES::

        sage: E = EllipticCurve([7,0])
        sage: E.weierstrass_p(prec=10, algorithm='quadratic')
        z^-2 - 7/5*z^2 + 49/75*z^6 + O(z^10)

        sage: E = EllipticCurve(GF(103),[1,2])
        sage: E.weierstrass_p(algorithm='quadratic')
        z^-2 + 41*z^2 + 88*z^4 + 11*z^6 + 57*z^8 + 55*z^10 + 73*z^12 + 11*z^14 + 17*z^16 + 50*z^18 + O(z^20)

        sage: from sage.schemes.elliptic_curves.ell_wp import compute_wp_quadratic
        sage: compute_wp_quadratic(E.base_ring(), E.a4(), E.a6(), prec=10)
        z^-2 + 41*z^2 + 88*z^4 + 11*z^6 + 57*z^8 + O(z^10)

    """
    m = (prec + 1)//2
    c = [0 for j in range(m)]
    c[0] = -A/5
    c[1] = -B/7

    # first Z represent z^2
    R = LaurentSeriesRing(k,'z')
    Z = R.gen()
    pe = Z**-1 + c[0]*Z + c[1]*Z**2

    for i in range(3, m):
        t = 0
        for j in range(1, i - 1):
            t += c[j-1]*c[i-2-j]
        ci = (3*t)/((i-2)*(2*i+3))
        pe += ci * Z**i
        c[i-1] = ci

    return pe(Z**2).add_bigoh(prec)

def compute_wp_fast(k, A, B, m):
    r"""
    Computes the Weierstrass function of an elliptic curve defined by short Weierstrass model: `y^2 = x^3 + Ax + B`. It does this with as fast as polynomial of degree `m` can be multiplied together in the base ring, i.e. `O(M(n))` in the notation of [BMSS2006]_.

    Let `p` be the characteristic of the underlying field: Then we must have either `p=0`, or `p > m + 3`.

    INPUT:

     - ``k`` - the base field of the curve
     - ``A`` - and
     - ``B`` - as the coefficients of the short Weierstrass model `y^2 = x^3 +Ax +B`, and
     - ``m`` - the precision to which the function is computed to.

    OUTPUT:

    the Weierstrass `\wp` function as a Laurent series to precision `m`.

    ALGORITHM:

    This function uses the algorithm described in section 3.3 of
    [BMSS2006]_.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_wp import compute_wp_fast
        sage: compute_wp_fast(QQ, 1, 8, 7)
        z^-2 - 1/5*z^2 - 8/7*z^4 + 1/75*z^6 + O(z^7)

        sage: k = GF(37)
        sage: compute_wp_fast(k, k(1), k(8), 5)
        z^-2 + 22*z^2 + 20*z^4 + O(z^5)

    """
    R = PowerSeriesRing(k,'z',default_prec=m+5)
    z = R.gen()
    s = 2
    f1 = z.add_bigoh(m+3)
    n = 2*m + 4

    # solve the nonlinear differential equation
    while (s < n):
        f1pr = f1.derivative()
        next_s = 2*s - 1

        a = 2*f1pr
        b = -(6*B*(f1**5) + 4*A*(f1**3))
        c = B*(f1**6) + A*f1**4 + 1 - (f1pr**2)

        # we should really be computing only mod z^next_s here.
        # but we loose only a factor 2
        f2 = solve_linear_differential_system(a, b, c, 0)
        # sometimes we get to 0 quicker than s reaches n
        if f2 == 0:
            break
        f1 = f1 + f2
        s = next_s

    R = f1
    Q = R**2
    pe = 1/Q

    return pe


def solve_linear_differential_system(a, b, c, alpha):
    r"""
    Solves a system of linear differential equations: `af' + bf = c` and `f'(0) = \alpha`
    where `a`, `b`, and `c` are power series in one variable and `\alpha` is a constant in the coefficient ring.

    ALGORITHM:

    due to Brent and Kung '78.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_wp import solve_linear_differential_system
        sage: k = GF(17)
        sage: R.<x> = PowerSeriesRing(k)
        sage: a = 1+x+O(x^7); b = x+O(x^7); c = 1+x^3+O(x^7); alpha = k(3)
        sage: f = solve_linear_differential_system(a,b,c,alpha)
        sage: f
        3 + x + 15*x^2 + x^3 + 10*x^5 + 3*x^6 + 13*x^7 + O(x^8)
        sage: a*f.derivative()+b*f - c
        O(x^7)
        sage: f(0) == alpha
        True

    """
    a_recip = 1/a
    B =  b * a_recip
    C =  c * a_recip
    int_B = B.integral()
    J = int_B.exp()
    J_recip = 1/J
    CJ = C * J
    int_CJ = CJ.integral()
    f =  J_recip * (alpha + int_CJ)

    return f
