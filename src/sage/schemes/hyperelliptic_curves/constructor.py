r"""
Constructor for hyperelliptic curves using the smooth model

Uses weighted projective space `\mathbb{P}(1 : g + 1 : 1)`.

TODO:

- We do not yet support the construction of hyperelliptic curves over rings

AUTHORS:

- David Kohel (2006): initial version
- Anna Somoza (2019-04): dynamic class creation
- Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2024): adapt to smooth model
"""

# ****************************************************************************
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#                     2019 Anna Somoza <anna.somoza.henares@gmail.com>
#                     2025 Sabrina Kunzweiler <sabrina.kunzweiler@math.u-bordeaux.fr>
#                     2025 Gareth Ma <grhkm21@gmail.com>
#                     2025 Giacomo Pope <giacomopope@gmail.com>
#
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.finite_fields import FiniteFields
from sage.functions.other import ceil
from sage.rings.abc import pAdicField
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import RationalField
from sage.schemes.hyperelliptic_curves.hyperelliptic_finite_field import (
    HyperellipticCurve_finite_field,
)
from sage.schemes.hyperelliptic_curves.hyperelliptic_g2 import (
    HyperellipticCurve_g2,
    HyperellipticCurve_g2_finite_field,
    HyperellipticCurve_g2_padic_field,
    HyperellipticCurve_g2_rational_field,
)
from sage.schemes.hyperelliptic_curves.hyperelliptic_generic import (
    HyperellipticCurve_generic,
)
from sage.schemes.hyperelliptic_curves.hyperelliptic_padic_field import (
    HyperellipticCurve_padic_field,
)
from sage.schemes.hyperelliptic_curves.hyperelliptic_rational_field import (
    HyperellipticCurve_rational_field,
)


def HyperellipticCurve(
    f, h=0, names = ['x', 'y'], check_squarefree: bool = True, distinguished_point=None
):
    r"""
    Constructor function for creating a hyperelliptic curve with
    smooth model with polynomials `f`, `h`.

    In Sage, a hyperelliptic curve of genus `g` is always
    specified by an (affine) equation in Weierstrass form

    .. MATH::
        y^2 + h(x) y = f(x),

    for some polynomials `h` and `f`. This defines a smooth
    model in weighted projective space `\mathbb{P}(1 : g + 1 : 1)`

    .. MATH::
        Y^2 + H(X,Z) Y = F(X,Z),

    where `H` is the degree `g + 1` homogenization of `h`,
    and `F` is the degree `2 g + 2` homogenization of `f`.

    INPUT:

    -  ``f`` -- polynomial

    -  ``h`` (default: ``0``) -- polynomial

    - ``names`` -- (default: ``['x', 'y']``) names for the coordinate functions

    -  ``check_squarefree`` (default: ``True``) -- test if
       the input defines a hyperelliptic curve

    - ``distinguished_point`` (default: ``None``) -- set a custom
       point to be the distinguished point for arithmetic. One of
       points at infinity if one so exists.

    EXAMPLES:

    We create a hyperelliptic curve over the rationals::

        sage: R.<x> = PolynomialRing(QQ)
        sage: H = HyperellipticCurve(x^8+1, x^4+1); H
        Hyperelliptic Curve over Rational Field defined by y^2 + (x^4 + 1)*y = x^8 + 1

    This hyperelliptic curve has no points at infinity, i.e. `H` is inert::

        sage: H.points_at_infinity()
        []
        sage: H.is_inert()
        True

    We can extend the base field to obtain a hyperelliptic curve with two points at infinity::

        sage: K.<alpha> = QQ.extension(x^2+x-1)
        sage: HK = H.change_ring(K)
        sage: HK.points_at_infinity()
        [(1 : alpha : 0), (1 : -alpha - 1 : 0)]
        sage: HK.is_split()
        True

    The construction of hyperelliptic curves is supported over different fields. The correct class is chosen automatically::

        sage: F = FiniteField(13)
        sage: S.<x> = PolynomialRing(F)
        sage: HF = HyperellipticCurve(x^5 + x^4 + x^3 + x^2 + x + 1, x^3 + x); HF
        Hyperelliptic Curve over Finite Field of size 13 defined by y^2 + (x^3 + x)*y = x^5 + x^4 + x^3 + x^2 + x + 1
        sage: type(HF)
        <class 'sage.schemes.hyperelliptic_curves.hyperelliptic_g2.HyperellipticCurve_g2_finite_field_with_category'>
        sage: Q5 = Qp(5,10)
        sage: T.<x> = Q5[]
        sage: H5 = HyperellipticCurve(x^7 + 1); H5
        Hyperelliptic Curve over 5-adic Field with capped relative precision 10 defined by (1 + O(5^10))*y^2 = (1 + O(5^10))*x^7 + 1 + O(5^10)
        sage: type(H5)
        <class 'sage.schemes.hyperelliptic_curves.hyperelliptic_padic_field.HyperellipticCurve_padic_field_with_category'>

    The input polynomials need not be monic::

        sage: R.<x> = QQ[]
        sage: HyperellipticCurve(3*x^5+1)
        Hyperelliptic Curve over Rational Field defined by y^2 = 3*x^5 + 1

    The polynomials `f` and `h` need to define a smooth curve of genus at
    least one. In particular polynomials defining elliptic curves are
    allowed as input::

        sage: E = HyperellipticCurve(x^3+1)
        sage: E.genus()
        1
        sage: HyperellipticCurve(x)
        Traceback (most recent call last):
        ...
        ValueError: arguments f = x and h = 0 must define a curve of genus at least one.

    The following polynomials define a singular curve and are
    not allowed as input::

        sage: C = HyperellipticCurve(x^6 + 2*x - 1, 2*x - 2)
        Traceback (most recent call last):
        ...
        ValueError: singularity in the provided affine patch

    The constructor accepts a distinguished point as an argument, allowing
    the user to pick which point at infinity is selected in the case there
    is more than one for the given curve model::

        sage: C = HyperellipticCurve(x^5 + x + 2, distinguished_point=(1, -2))
        sage: C.distinguished_point()
        (1 : -2 : 1)

    We can change the names of the variables in the output::

        sage: k.<a> = GF(9); R.<x> = k[]                                                # needs sage.rings.finite_rings
        sage: HyperellipticCurve(x^3 + x - 1, x + a, names=['u','v'])                   # needs sage.rings.finite_rings
        Hyperelliptic Curve over Finite Field in a of size 3^2
         defined by v^2 + (u + a)*v = u^3 + u + 2

    """

    # ---------------------------
    # Internal Helper functions
    # ---------------------------

    def genus(f, h) -> Integer:
        r"""
        Helper function to compute the genus of a hyperelliptic curve
        defined by `y^2 + h(x)y = f(x)`.
        """
        # Some classes still have issues with degrees returning `int`
        # rather than Sage Integer types, so convert the output.
        df = f.degree()
        dh_2 = 2 * h.degree()
        if dh_2 < df:
            return Integer((df - 1) // 2)
        return Integer((dh_2 - 1) // 2)

    def check_no_affine_singularities(f, h):
        r"""
        Helper function which determines whether there are any
        affine singularities in the curve `y^2 + h(x)y = f(x)`.
        """
        if f.base_ring().characteristic() == 2:
            if h.is_zero():
                return False
            if h.is_constant():
                return True
            return h.gcd(f.derivative() ** 2 - f * h.derivative() ** 2).is_one()

        if h.is_zero():
            return f.gcd(f.derivative()).is_one()

        g1 = h**2 + 4 * f
        g2 = 2 * f.derivative() + h * h.derivative()
        return g1.gcd(g2).is_one()

    def homogenize_defining_polynomial(f, h):
        r"""
        Compute the (homogenised weighted projective) defining polynomial of
        the hyperelliptic curve.
        """
        X, Y, Z = PolynomialRing(f.base_ring(), names="X, Y, Z").gens()

        # Some classes still have issues with degrees returning `int`
        d = Integer(max(h.degree(), ceil(f.degree() / 2)))
        F = sum(f[i] * X**i * Z ** (2 * d - i) for i in range(2 * d + 1))

        if h.is_zero():
            return Y**2 - F
        H = sum(h[i] * X**i * Z ** (d - i) for i in range(d + 1))
        return Y**2 + H * Y - F

    # -------------------------------------------
    # Typechecking and projective model creation
    # -------------------------------------------

    # Check the polynomials are of the right type
    # D is the discriminant; use this for the type check
    # rather than f and h, one of which might be constant.
    D = h**2 + 4 * f
    if not isinstance(D, Polynomial):
        raise TypeError(f"arguments f = {f} and h = {h} must be polynomials")

    # Store the hyperelliptic polynomials as the correct type
    polynomial_ring = D.parent()
    base_ring = D.base_ring()
    f = polynomial_ring(f)
    h = polynomial_ring(h)

    # Ensure that there are no affine singular points
    if check_squarefree and not check_no_affine_singularities(f, h):
        raise ValueError("singularity in the provided affine patch")

    # Compute the genus of the curve from f, h
    curve_genus = genus(f, h)
    if curve_genus == 0:
        raise ValueError(
            f"arguments f = {f} and h = {h} must define a curve of genus at least one."
        )

    # Compute the smooth model for the hyperelliptic curve
    # using a weighted projective space
    defining_polynomial = homogenize_defining_polynomial(f, h)

    # -----------------------
    # Class selection
    # -----------------------

    # Special class for finite fields
    if base_ring in FiniteFields():
        if curve_genus == 2:
            cls = HyperellipticCurve_g2_finite_field
        else:
            cls = HyperellipticCurve_finite_field
    # Special class for pAdic fields
    elif isinstance(base_ring, pAdicField):
        if curve_genus == 2:
            cls = HyperellipticCurve_g2_padic_field
        else:
            cls = HyperellipticCurve_padic_field
    # Special class for rational fields
    elif isinstance(base_ring, RationalField):
        if curve_genus == 2:
            cls = HyperellipticCurve_g2_rational_field
        else:
            cls = HyperellipticCurve_rational_field
    # Default class for all other fields
    elif curve_genus == 2:
        cls = HyperellipticCurve_g2
    else:
        cls = HyperellipticCurve_generic

    H = cls(defining_polynomial, f, h, curve_genus, names=names)
    if distinguished_point:
        H.set_distinguished_point(H(distinguished_point))
    return H
