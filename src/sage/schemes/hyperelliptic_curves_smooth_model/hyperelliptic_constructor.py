"""
Constructor for Hyperelliptic Curves using the smooth model

Adapted from /hyperelliptic/constructor.py

AUTHORS:

- David Kohel (2006): initial version
- HyperellipticCurveSmooth-Team (TODO)
"""

# ****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#                2019 Anna Somoza <anna.somoza.henares@gmail.com>
#                2024 Hyperelliptic-Team (TODO)
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.finite_fields import FiniteFields
from sage.rings.abc import pAdicField
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.rational_field import is_RationalField
from sage.schemes.hyperelliptic_curves_smooth_model.hyperelliptic_finite_field import (
    HyperellipticCurveSmoothModel_finite_field,
)
from sage.schemes.hyperelliptic_curves_smooth_model.hyperelliptic_g2 import (
    HyperellipticCurveSmoothModel_g2,
    HyperellipticCurveSmoothModel_g2_finite_field,
    HyperellipticCurveSmoothModel_g2_padic_field,
    HyperellipticCurveSmoothModel_g2_rational_field,
)
from sage.schemes.hyperelliptic_curves_smooth_model.hyperelliptic_generic import (
    HyperellipticCurveSmoothModel_generic,
)
from sage.schemes.hyperelliptic_curves_smooth_model.hyperelliptic_padic_field import (
    HyperellipticCurveSmoothModel_padic_field,
)
from sage.schemes.hyperelliptic_curves_smooth_model.hyperelliptic_rational_field import (
    HyperellipticCurveSmoothModel_rational_field,
)
from sage.schemes.toric.library import toric_varieties

"""
TODO:

- We currently cannot support the construction of curves over rings

"""


def HyperellipticCurveSmoothModel(f, h=0, check_squarefree=True):
    r"""
    Returns the hyperelliptic curve `y^2 + h y = f`, represented using
    the smooth model over weighted projective space. This is built from
    univariate polynomials `h` and `f`. If `h` is not given, then it
    defaults to 0.

    INPUT:

    -  ``f`` -- univariate polynomial

    -  ``h`` -- optional univariate polynomial

    -  ``check_squarefree`` (default: ``True``) -- test if
       the input defines a hyperelliptic curve when f is
       homogenized to degree `2g+2` and h to degree
       `g+1` for some g.

    .. WARNING::

        Unlike the original HyperellipticCurve which uses the projective
        model, this implementation currently only supports Hyperelliptic
        curves defined over fields (not generic rings).

    .. NOTE::

        The words "hyperelliptic curve" are normally only used for curves of
        genus at least two, but this class allows more general smooth double
        covers of the projective line (conics and elliptic curves), even though
        the class is not meant for those and some outputs may be incorrect.

    EXAMPLES:

    Basic examples::

        sage: R.<x> = QQ[]
        sage: HyperellipticCurveSmoothModel(x^5 + x + 1)
        Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x + 1
        sage: HyperellipticCurveSmoothModel(x^19 + x + 1, x - 2)
        Hyperelliptic Curve over Rational Field defined by y^2 + (x - 2)*y = x^19 + x + 1

        sage: k.<a> = GF(9); R.<x> = k[]                                                # needs sage.rings.finite_rings
        sage: HyperellipticCurveSmoothModel(x^3 + x - 1, x+a)                           # needs sage.rings.finite_rings
        Hyperelliptic Curve over Finite Field in a of size 3^2
         defined by y^2 + (x + a)*y = x^3 + x + 2

    Characteristic two::

        sage: # needs sage.rings.finite_rings
        sage: P.<x> = GF(8, 'a')[]
        sage: HyperellipticCurveSmoothModel(x^7 + 1, x)
        Hyperelliptic Curve over Finite Field in a of size 2^3
         defined by y^2 + x*y = x^7 + 1
        sage: HyperellipticCurveSmoothModel(x^8 + x^7 + 1, x^4 + 1)
        Hyperelliptic Curve over Finite Field in a of size 2^3
         defined by y^2 + (x^4 + 1)*y = x^8 + x^7 + 1
        sage: HyperellipticCurveSmoothModel(x^8 + 1, x)
        Traceback (most recent call last):
        ...
        ValueError: not a hyperelliptic curve: highly singular at infinity
        sage: HyperellipticCurveSmoothModel(x^8 + x^7 + 1, x^4)
        Traceback (most recent call last):
        ...
        ValueError: not a hyperelliptic curve: singularity in the provided affine patch

        sage: F.<t> = PowerSeriesRing(FiniteField(2))
        sage: P.<x> = PolynomialRing(FractionField(F))
        sage: HyperellipticCurveSmoothModel(x^5 + t, x)
        Hyperelliptic Curve over Laurent Series Ring in t over Finite Field of size 2
         defined by y^2 + x*y = x^5 + t

    We can change the names of the variables in the output::

        sage: k.<a> = GF(9); R.<x> = k[]                                                # needs sage.rings.finite_rings
        sage: HyperellipticCurveSmoothModel(x^3 + x - 1, x + a, names=['X','Y'])                   # needs sage.rings.finite_rings
        Hyperelliptic Curve over Finite Field in a of size 3^2
         defined by Y^2 + (X + a)*Y = X^3 + X + 2

    This class also allows curves of genus zero or one, which are strictly
    speaking not hyperelliptic::

        sage: P.<x> = QQ[]
        sage: HyperellipticCurveSmoothModel(x^2 + 1)
        Hyperelliptic Curve over Rational Field defined by y^2 = x^2 + 1
        sage: HyperellipticCurveSmoothModel(x^4 - 1)
        Hyperelliptic Curve over Rational Field defined by y^2 = x^4 - 1
        sage: HyperellipticCurveSmoothModel(x^3 + 2*x + 2)
        Hyperelliptic Curve over Rational Field defined by y^2 = x^3 + 2*x + 2

    Double roots::

        sage: P.<x> = GF(7)[]
        sage: HyperellipticCurve((x^3-x+2)^2*(x^6-1))
        Traceback (most recent call last):
        ...
        ValueError: not a hyperelliptic curve: singularity in the provided affine patch

        sage: HyperellipticCurve((x^3-x+2)^2*(x^6-1), check_squarefree=False)
        Hyperelliptic Curve over Finite Field of size 7 defined by
         y^2 = x^12 + 5*x^10 + 4*x^9 + x^8 + 3*x^7 + 3*x^6 + 2*x^4 + 3*x^3 + 6*x^2 + 4*x + 3

    The input for a (smooth) hyperelliptic curve of genus `g` should not
    contain polynomials of degree greater than `2g+2`. In the following
    example, the hyperelliptic curve has genus 2 and there exists a model
    `y^2 = F` of degree 6, so the model `y^2 + yh = f` of degree 200 is not
    allowed.::

        sage: P.<x> = QQ[]
        sage: h = x^100
        sage: F = x^6 + 1
        sage: f = F - h^2/4
        sage: HyperellipticCurve(f, h)
        Traceback (most recent call last):
        ...
        ValueError: not a hyperelliptic curve: highly singular at infinity

        sage: HyperellipticCurve(F)
        Hyperelliptic Curve over Rational Field defined by y^2 = x^6 + 1

    An example with a singularity over an inseparable extension of the
    base field::

        sage: F.<t> = GF(5)[]
        sage: P.<x> = F[]
        sage: HyperellipticCurve(x^5 + t)
        Traceback (most recent call last):
        ...
        ValueError: not a hyperelliptic curve: singularity in the provided affine patch

    Input with integer coefficients creates objects with the integers
    as base ring, but only checks smoothness over `\QQ`, not over Spec(`\ZZ`).
    In other words, it is checked that the discriminant is non-zero, but it is
    not checked whether the discriminant is a unit in `\ZZ^*`.::

        sage: P.<x> = ZZ[]
        sage: HyperellipticCurve(3*x^7 + 6*x + 6)
        Hyperelliptic Curve over Integer Ring defined by y^2 = 3*x^7 + 6*x + 6

    TESTS:

    Check that `f` can be a constant (see :issue:`15516`)::

        sage: R.<u> = PolynomialRing(Rationals())
        sage: HyperellipticCurve(-12, u^4 + 7)
        Hyperelliptic Curve over Rational Field defined by y^2 + (x^4 + 7)*y = -12

    Check that two curves with the same class name have the same class type::

        sage: # needs sage.rings.finite_rings
        sage: R.<t> = PolynomialRing(GF(next_prime(10^9)))
        sage: C = HyperellipticCurve(t^5 + t + 1)
        sage: C2 = HyperellipticCurve(t^5 + 3*t + 1)
        sage: type(C2) == type(C)
        True

    Check that the inheritance is correct::

        sage: # needs sage.rings.finite_rings
        sage: R.<t> = PolynomialRing(GF(next_prime(10^9)))
        sage: C = HyperellipticCurve(t^5 + t + 1)
        sage: type(C).mro()
        [<class 'sage.schemes.hyperelliptic_curves.constructor.HyperellipticCurve_g2_FiniteField_with_category'>,
         <class 'sage.schemes.hyperelliptic_curves.constructor.HyperellipticCurve_g2_FiniteField'>,
         <class 'sage.schemes.hyperelliptic_curves.hyperelliptic_g2.HyperellipticCurve_g2'>,
         <class 'sage.schemes.hyperelliptic_curves.hyperelliptic_finite_field.HyperellipticCurve_finite_field'>,
         <class 'sage.schemes.hyperelliptic_curves.hyperelliptic_generic.HyperellipticCurve_generic'>,
        ...]
    """

    # -----------------
    # Helper functions
    # -----------------

    def __genus(f, h):
        # Some classes still have issues with degrees returning `int`
        df = Integer(f.degree())
        dh_2 = 2 * Integer(h.degree())
        if dh_2 < df:
            return (df - 1) // 2
        return (dh_2 - 1) // 2

    def __check_no_affine_singularities(f, h):
        if f.base_ring().characteristic() == 2:
            if h.is_zero():
                return False
            elif h.is_constant():
                return True
            return h.gcd(f.derivative() ** 2 - f * h.derivative() ** 2).is_one()

        if h.is_zero():
            return f.gcd(f.derivative()).is_one()

        g1 = h**2 + 4 * f
        g2 = 2 * f.derivative() + h * h.derivative()
        return g1.gcd(g2).is_one()

    def __projective_model(f, h, genus):
        """
        Compute the weighted projective model (1 : g + 1 : 1)

        EXAMPLES::

        TODO
        """
        T = toric_varieties.WP(
            [1, genus + 1, 1], base_ring=f.base_ring(), names="X, Y, Z"
        )
        (X, Y, Z) = T.gens()

        # Some classes still have issues with degrees returning `int`
        d = max(Integer(h.degree()), (Integer(f.degree()) / 2).ceil())
        F = sum(f[i] * X**i * Z ** (2 * d - i) for i in range(2 * d + 1))

        if h.is_zero():
            G = Y**2 - F
        else:
            H = sum(h[i] * X**i * Z ** (d - i) for i in range(d + 1))
            G = Y**2 + H * Y - F

        return T.subscheme(G)

    # -------------------------------------------
    # Typechecking and projective model creation
    # -------------------------------------------

    # Check the polynomials are of the right type
    F = h**2 + 4 * f
    if not isinstance(F, Polynomial):
        raise TypeError(f"arguments {f = } and {h = } must be polynomials")

    # Store the hyperelliptic polynomials as the correct type
    polynomial_ring = F.parent()
    base_ring = F.base_ring()
    f = polynomial_ring(f)
    h = polynomial_ring(h)

    # Ensure that there are no affine singular points
    if check_squarefree and not __check_no_affine_singularities(f, h):
        raise ValueError("singularity in the provided affine patch")

    # Compute the genus of the curve from f, h
    genus = __genus(f, h)

    # Compute the smooth model for the hyperelliptic curve
    # using a weighted projective space (via Toric Variety)
    projective_model = __projective_model(f, h, genus)

    # -----------------------
    # Class selection
    # -----------------------

    # Special class for finite fields
    if base_ring in FiniteFields():
        if genus == 2:
            cls = HyperellipticCurveSmoothModel_g2_finite_field
        else:
            cls = HyperellipticCurveSmoothModel_finite_field
    # Special class for padic fields
    elif isinstance(base_ring, pAdicField):
        if genus == 2:
            cls = HyperellipticCurveSmoothModel_g2_padic_field
        else:
            cls = HyperellipticCurveSmoothModel_padic_field
    # Special class for rational fields
    elif is_RationalField(base_ring):
        if genus == 2:
            cls = HyperellipticCurveSmoothModel_g2_rational_field
        else:
            cls = HyperellipticCurveSmoothModel_rational_field
    # Default class for all other fields
    elif genus == 2:
        cls = HyperellipticCurveSmoothModel_g2
    else:
        cls = HyperellipticCurveSmoothModel_generic

    return cls(projective_model, f, h, genus)
