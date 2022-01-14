# -*- coding: utf-8 -*-
r"""
Isogenies

An isogeny `\varphi: E_1\to E_2` between two elliptic curves `E_1` and
`E_2` is a morphism of curves that sends the origin of `E_1` to the
origin of `E_2`. Such a morphism is automatically a morphism of group
schemes and the kernel is a finite subgroup scheme of `E_1`.  Such a
subscheme can either be given by a list of generators, which have to
be torsion points, or by a polynomial in the coordinate `x` of the
Weierstrass equation of `E_1`.

The usual way to create and work with isogenies is illustrated with
the following example::

    sage: k = GF(11)
    sage: E = EllipticCurve(k,[1,1])
    sage: Q = E(6,5)
    sage: phi = E.isogeny(Q)
    sage: phi
    Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 11 to Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 11
    sage: P = E(4,5)
    sage: phi(P)
    (10 : 0 : 1)
    sage: phi.codomain()
    Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 11
    sage: phi.rational_maps()
    ((x^7 + 4*x^6 - 3*x^5 - 2*x^4 - 3*x^3 + 3*x^2 + x - 2)/(x^6 + 4*x^5 - 4*x^4 - 5*x^3 + 5*x^2), (x^9*y - 5*x^8*y - x^7*y + x^5*y - x^4*y - 5*x^3*y - 5*x^2*y - 2*x*y - 5*y)/(x^9 - 5*x^8 + 4*x^6 - 3*x^4 + 2*x^3))

The functions directly accessible from an elliptic curve ``E`` over a
field are ``isogeny`` and ``isogeny_codomain``.

The most useful functions that apply to isogenies are

- ``codomain``
- ``degree``
- ``domain``
- ``dual``
- ``rational_maps``
- ``kernel_polynomial``

.. WARNING::

   Only cyclic, separable isogenies are implemented (except for [2]). Some
   algorithms may need the isogeny to be normalized.

AUTHORS:

- Daniel Shumow <shumow@gmail.com>: 2009-04-19: initial version

- Chris Wuthrich : 7/09: changes: add check of input, not the full list is needed.
  10/09: eliminating some bugs.

- John Cremona 2014-08-08: tidying of code and docstrings, systematic
  use of univariate vs. bivariate polynomials and rational functions.

"""

# ****************************************************************************
#       Copyright (C) 2009 Daniel Shumow <shumow@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy

from sage.structure.sequence import Sequence

from sage.schemes.elliptic_curves.hom import EllipticCurveHom

from sage.rings.all import PolynomialRing, Integer, LaurentSeriesRing
from sage.rings.polynomial.polynomial_element import is_Polynomial
from sage.schemes.elliptic_curves.all import EllipticCurve
from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve

from sage.rings.number_field.number_field_base import is_NumberField

from sage.schemes.elliptic_curves.weierstrass_morphism \
        import WeierstrassIsomorphism, isomorphisms, baseWI

from sage.sets.set import Set

#
# Private function for parsing input to determine the type of
# algorithm
#
def isogeny_determine_algorithm(E, kernel):
    r"""
    Helper function that allows the various isogeny functions to infer
    the algorithm type from the parameters passed in.

    INPUT:

    - ``E`` (elliptic curve) -- an elliptic curve

    - ``kernel`` -- either a list of points on ``E``, or a univariate
      polynomial or list of coefficients of a univariate polynomial.

    OUTPUT:

    (string) either 'velu' or 'kohel'

    If ``kernel`` is a list of points on the EllipticCurve `E`, then
    we will try to use Velu's algorithm.

    If ``kernel`` is a list of coefficients or a univariate
    polynomial, we will try to use the Kohel's algorithms.

    EXAMPLES:

    This helper function will be implicitly called by the following examples::

        sage: R.<x> = GF(5)[]
        sage: E = EllipticCurve(GF(5), [0,0,0,1,0])

    We can construct the same isogeny from a kernel polynomial::

        sage: phi = EllipticCurveIsogeny(E, x+3)

    or from a list of coefficients of a kernel polynomial::

        sage: phi == EllipticCurveIsogeny(E, [3,1])
        True

    or from a rational point which generates the kernel::

        sage: phi == EllipticCurveIsogeny(E,  E((2,0)) )
        True

    In the first two cases, Kohel's algorithm will be used, while in
    the third case it is Velu::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogeny_determine_algorithm
        sage: isogeny_determine_algorithm(E, x+3)
        'kohel'
        sage: isogeny_determine_algorithm(E, [3, 1])
        'kohel'
        sage: isogeny_determine_algorithm(E, E((2,0)))
        'velu'
    """
    kernel_is_list = isinstance(kernel, list)

    if not kernel_is_list and kernel in E:
        kernel = [kernel]
        kernel_is_list = True

    if (is_Polynomial(kernel) or ( kernel_is_list) and (kernel[0] in E.base_ring()) ):
        algorithm = "kohel"
    elif (kernel_is_list) and (kernel[0] in E):
        # note that if kernel[0] is on an extension of E this
        # condition will be false
        algorithm = "velu"
    else:
        raise ValueError("Invalid Parameters to EllipticCurveIsogeny constructor.")
    return algorithm

def isogeny_codomain_from_kernel(E, kernel, degree=None):
    r"""
    Compute the isogeny codomain given a kernel.

    INPUT:

    - ``E`` - The domain elliptic curve.

    - ``kernel`` - Either a list of points in the kernel of the isogeny, or a
                   kernel polynomial (specified as a either a univariate
                   polynomial or a coefficient list.)

    - ``degree`` - an integer, (default:``None``)  optionally specified degree
                   of the kernel.

    OUTPUT:

    (elliptic curve) the codomain of the separable normalized isogeny
    from this kernel

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogeny_codomain_from_kernel
        sage: E = EllipticCurve(GF(7), [1,0,1,0,1])
        sage: R.<x> = GF(7)[]
        sage: isogeny_codomain_from_kernel(E, [4,1], degree=3)
        Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + 6 over Finite Field of size 7
        sage: EllipticCurveIsogeny(E, [4,1]).codomain() == isogeny_codomain_from_kernel(E, [4,1], degree=3)
        True
        sage: isogeny_codomain_from_kernel(E, x^3 + x^2 + 4*x + 3)
        Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + 6 over Finite Field of size 7
        sage: isogeny_codomain_from_kernel(E, x^3 + 2*x^2 + 4*x + 3)
        Elliptic Curve defined by y^2 + x*y + y = x^3 + 5*x + 2 over Finite Field of size 7

        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: kernel_list = [E((15,10)), E((10,3)),E((6,5))]
        sage: isogeny_codomain_from_kernel(E, kernel_list)
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 3*x + 15 over Finite Field of size 19

    """

    algorithm = isogeny_determine_algorithm(E, kernel)

    if ("velu"==algorithm):
        # if we are using Velu's formula, just instantiate the isogeny
        # and return the codomain
        return EllipticCurveIsogeny(E, kernel).codomain()
    elif ("kohel"==algorithm):
        return compute_codomain_kohel(E, kernel, degree)

def compute_codomain_formula(E, v, w):
    r"""
    Compute the codomain curve given parameters `v` and `w` (as in
    Velu / Kohel / etc formulas).

    INPUT:

    - ``E`` -- an elliptic curve

    - ``v``, ``w`` -- elements of the base field of ``E``

    OUTPUT:

    The elliptic curve with invariants
    `[a_1,a_2,a_3,a_4-5v,a_6-(a_1^2+4a_2)v-7w]` where
    `E=[a_1,a_2,a_3,a_4,a_6]`.

    EXAMPLES:

    This formula is used by every Isogeny instantiation::

        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: phi = EllipticCurveIsogeny(E, E((1,2)) )
        sage: phi.codomain()
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 9*x + 13 over Finite Field of size 19
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_codomain_formula
        sage: v = phi._EllipticCurveIsogeny__v
        sage: w = phi._EllipticCurveIsogeny__w
        sage: compute_codomain_formula(E, v, w) == phi.codomain()
        True
    """
    a1,a2,a3,a4,a6 = E.ainvs()

    A4 = a4 - 5*v
    A6 = a6 - (a1**2 + 4*a2)*v - 7*w

    return EllipticCurve([a1, a2, a3, A4, A6])

def compute_vw_kohel_even_deg1(x0, y0, a1, a2, a4):
    r"""
    Compute Velu's (v,w) using Kohel's formulas for isogenies of
    degree exactly divisible by 2.

    INPUT:

    - ``x0``, ``y0`` -- coordinates of a 2-torsion point on an elliptic curve E

    - ``a1``, ``a2``, ``a4`` -- invariants of E

    OUTPUT:

    (tuple) Velu's isogeny parameters (v,w).

    EXAMPLES:

    This function will be implicitly called by the following example::

        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: phi = EllipticCurveIsogeny(E, [9,1]); phi
        Isogeny of degree 2 from Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Finite Field of size 19 to Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 9*x + 8 over Finite Field of size 19
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_vw_kohel_even_deg1
        sage: a1,a2,a3,a4,a6 = E.ainvs()
        sage: x0 = -9
        sage: y0 = -(a1*x0 + a3)/2
        sage: compute_vw_kohel_even_deg1(x0, y0, a1, a2, a4)
        (18, 9)
    """
    v = (3*x0**2 + 2*a2*x0 + a4 - a1*y0)
    w = x0*v

    return (v,w)

def compute_vw_kohel_even_deg3(b2,b4,s1,s2,s3):
    r"""
    Compute Velu's (v,w) using Kohel's formulas for isogenies of
    degree divisible by 4.

    INPUT:

    - ``b2``, ``b4`` -- invariants of an elliptic curve E

    - ``s1``, ``s2``, ``s3`` -- signed coefficients of the 2-division
      polynomial of E

    OUTPUT:

    (tuple) Velu's isogeny parameters (v,w).

    EXAMPLES:

    This function will be implicitly called by the following example::

        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: R.<x> = GF(19)[]
        sage: phi = EllipticCurveIsogeny(E, x^3 + 7*x^2 + 15*x + 12); phi
        Isogeny of degree 4 from Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Finite Field of size 19 to Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 3*x + 15 over Finite Field of size 19
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_vw_kohel_even_deg3
        sage: (b2,b4) = (E.b2(), E.b4())
        sage: (s1, s2, s3) = (-7, 15, -12)
        sage: compute_vw_kohel_even_deg3(b2, b4, s1, s2, s3)
        (4, 7)
    """
    temp1 = (s1**2 - 2*s2)
    v = 3*temp1 + b2*s1/2 + 3*b4/2
    w = 3*(s1**3 - 3*s1*s2 + 3*s3) + b2*temp1/2 + b4*s1/2

    return (v,w)


def compute_vw_kohel_odd(b2,b4,b6,s1,s2,s3,n):
    r"""
    Compute Velu's (v,w) using Kohel's formulas for isogenies of odd
    degree.

    INPUT:

    - ``b2``, ``b4``, ``b6`` -- invariants of an elliptic curve E

    - ``s1``, ``s2``, ``s3`` -- signed coefficients of lowest powers
      of x in the kernel polynomial.

    - ``n`` (int) -- the degree

    OUTPUT:

    (tuple) Velu's isogeny parameters (v,w).

    EXAMPLES:

    This function will be implicitly called by the following example::

        sage: E = EllipticCurve(GF(19), [18,17,16,15,14])
        sage: R.<x> = GF(19)[]
        sage: phi = EllipticCurveIsogeny(E, x^3 + 14*x^2 + 3*x + 11); phi
        Isogeny of degree 7 from Elliptic Curve defined by y^2 + 18*x*y + 16*y = x^3 + 17*x^2 + 15*x + 14 over Finite Field of size 19 to Elliptic Curve defined by y^2 + 18*x*y + 16*y = x^3 + 17*x^2 + 18*x + 18 over Finite Field of size 19
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_vw_kohel_odd
        sage: (b2,b4,b6) = (E.b2(), E.b4(), E.b6())
        sage: (s1,s2,s3) = (-14,3,-11)
        sage: compute_vw_kohel_odd(b2,b4,b6,s1,s2,s3,3)
        (7, 1)

    """
    v = 6*(s1**2 - 2*s2) + b2*s1 + n*b4
    w = 10*(s1**3 - 3*s1*s2 + 3*s3) + 2*b2*(s1**2 - 2*s2) + 3*b4*s1 + n*b6

    return (v,w)


def compute_codomain_kohel(E, kernel, degree):
    r"""
    Compute the codomain from the kernel polynomial using Kohel's
    formulas.

    INPUT:

    - ``E`` -- an elliptic curve

    - ``kernel`` (polynomial or list) -- the kernel polynomial, or a
      list of its coefficients

    - ``degree`` (int) -- degree of the isogeny

    OUTPUT:

    (elliptic curve) -- the codomain elliptic curve ``E``/``kernel``

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_codomain_kohel
        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: phi = EllipticCurveIsogeny(E, [9,1])
        sage: phi.codomain() == isogeny_codomain_from_kernel(E, [9,1])
        True
        sage: compute_codomain_kohel(E, [9,1], 2)
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 9*x + 8 over Finite Field of size 19
        sage: R.<x> = GF(19)[]
        sage: E = EllipticCurve(GF(19), [18,17,16,15,14])
        sage: phi = EllipticCurveIsogeny(E, x^3 + 14*x^2 + 3*x + 11)
        sage: phi.codomain() == isogeny_codomain_from_kernel(E, x^3 + 14*x^2 + 3*x + 11)
        True
        sage: compute_codomain_kohel(E, x^3 + 14*x^2 + 3*x + 11, 7)
        Elliptic Curve defined by y^2 + 18*x*y + 16*y = x^3 + 17*x^2 + 18*x + 18 over Finite Field of size 19
        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: phi = EllipticCurveIsogeny(E, x^3 + 7*x^2 + 15*x + 12)
        sage: isogeny_codomain_from_kernel(E, x^3 + 7*x^2 + 15*x + 12) == phi.codomain()
        True
        sage: compute_codomain_kohel(E, x^3 + 7*x^2 + 15*x + 12,4)
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 3*x + 15 over Finite Field of size 19

    .. NOTE::

       This function uses the formulas of Section 2.4 of [Koh1996]_.
    """
    # First set up the polynomial ring

    base_field = E.base_ring()
    poly_ring = PolynomialRing(base_field,'x')

    if (is_Polynomial(kernel)):
        psi = poly_ring(kernel)
        kernel_list = psi.list()
    elif isinstance(kernel, list) and (kernel[0] in base_field):
        kernel_list = kernel
        psi = poly_ring(kernel_list)
    else:
        raise ValueError("Invalid input to compute_codomain_kohel")

    # next determine the even / odd part of the isogeny
    psi_2tor = two_torsion_part(E, psi)

    if (0 != psi_2tor.degree()): # even degree case

        psi_quo = psi//psi_2tor

        if (0 != psi_quo.degree()):
            raise ArithmeticError("For basic Kohel's algorithm, if the kernel degree is even then the kernel must be contained in the two torsion.")

        n = psi_2tor.degree()

        if (1 == n): # degree divisible exactly by 2

            a1,a2,a3,a4,a6 = E.ainvs()

            x0 = -psi_2tor.constant_coefficient()

            # determine y0
            if (2 == base_field.characteristic()):
                y0 = (x0**3 + a2*x0**2 + a4*x0 + a6).sqrt()
            else:
                y0 = -(a1*x0 + a3)/2

            # now (x0,y0) is the 2-torsion point in the kernel

            (v,w) = compute_vw_kohel_even_deg1(x0,y0,a1,a2,a4)

        elif (3 == n): # psi_2tor is the full 2-division polynomial

            b2 = E.b2()
            b4 = E.b4()

            s = psi_2tor.list()
            s1 = -s[n-1]
            s2 = s[n-2]
            s3 = -s[n-3]

            (v,w) = compute_vw_kohel_even_deg3(b2,b4,s1,s2,s3)

    else:  # odd degree case

        n = psi.degree()

        b2 = E.b2()
        b4 = E.b4()
        b6 = E.b6()

        s1 = 0
        s2 = 0
        s3 = 0

        if 1 <= n:
            s1 = -kernel_list[n-1]

        if 2 <= n:
            s2 = kernel_list[n-2]

        if 3 <= n:
            s3 = -kernel_list[n-3]

        # initializing these allows us to calculate E2.
        (v,w) = compute_vw_kohel_odd(b2,b4,b6,s1,s2,s3,n)

    return compute_codomain_formula(E, v, w)


def two_torsion_part(E, psi):
    r"""
    Return the greatest common divisor of ``psi`` and the 2 torsion
    polynomial of `E`.

    INPUT:

    - ``E`` -- an elliptic curve

    - ``psi`` -- a univariate polynomial over the base field of ``E``

    OUTPUT:

    (polynomial) the gcd of psi and the 2-torsion polynomial of ``E``.

    EXAMPLES:

    Every function that computes the kernel polynomial via Kohel's
    formulas will call this function::

        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: R.<x> = GF(19)[]
        sage: phi = EllipticCurveIsogeny(E, x + 13)
        sage: isogeny_codomain_from_kernel(E, x + 13) == phi.codomain()
        True
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import two_torsion_part
        sage: two_torsion_part(E, x+13)
        x + 13

    """
    x = psi.parent().gen() # NB psi is univariate but could be constant
    psi_2 = E.two_division_polynomial(x)
    return psi.gcd(psi_2)

class EllipticCurveIsogeny(EllipticCurveHom):
    r"""
    Class Implementing Isogenies of Elliptic Curves

    This class implements cyclic, separable, normalized isogenies of
    elliptic curves.

    Several different algorithms for computing isogenies are
    available.  These include:

    - Velu's Formulas: Velu's original formulas for computing
      isogenies.  This algorithm is selected by giving as the
      ``kernel`` parameter a list of points which generate a finite
      subgroup.

    - Kohel's Formulas: Kohel's original formulas for computing
      isogenies.  This algorithm is selected by giving as the
      ``kernel`` parameter a monic polynomial (or a coefficient list
      (little endian)) which will define the kernel of the isogeny.

    INPUT:

    - ``E`` -- an elliptic curve, the domain of the isogeny to
      initialize.

    - ``kernel`` -- a kernel, either a point in ``E``, a list of
      points in ``E``, a monic kernel polynomial, or ``None``.  If
      initializing from a domain/codomain, this must be set to None.

    - ``codomain`` -- an elliptic curve (default:``None``).  If
      ``kernel`` is ``None``, then this must be the codomain of a
      cyclic, separable, normalized isogeny, furthermore, ``degree``
      must be the degree of the isogeny from ``E`` to ``codomain``. If
      ``kernel`` is not ``None``, then this must be isomorphic to the
      codomain of the cyclic normalized separable isogeny defined by
      ``kernel``, in this case, the isogeny is post composed with an
      isomorphism so that this parameter is the codomain.

    - ``degree`` -- an integer (default:``None``).  If ``kernel`` is
      ``None``, then this is the degree of the isogeny from ``E`` to
      ``codomain``.  If ``kernel`` is not ``None``, then this is used
      to determine whether or not to skip a gcd of the kernel
      polynomial with the two torsion polynomial of ``E``.

    - ``model`` -- a string (default:``None``).  Only supported
      variable is ``minimal``, in which case if ``E`` is a curve over
      the rationals or over a number field, then the codomain is a
      global minimum model where this exists.

    - ``check`` (default: ``True``) checks if the input is valid to
      define an isogeny

    EXAMPLES:

    A simple example of creating an isogeny of a field of small
    characteristic::

        sage: E = EllipticCurve(GF(7), [0,0,0,1,0])
        sage: phi = EllipticCurveIsogeny(E, E((0,0)) ); phi
        Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 3*x over Finite Field of size 7
        sage: phi.degree() == 2
        True
        sage: phi.kernel_polynomial()
        x
        sage: phi.rational_maps()
        ((x^2 + 1)/x, (x^2*y - y)/x^2)
        sage: phi == loads(dumps(phi))  # known bug
        True

    A more complicated example of a characteristic 2 field::

        sage: E = EllipticCurve(GF(2^4,'alpha'), [0,0,1,0,1])
        sage: P = E((1,1))
        sage: phi_v = EllipticCurveIsogeny(E, P); phi_v
        Isogeny of degree 3 from Elliptic Curve defined by y^2 + y = x^3 + 1 over Finite Field in alpha of size 2^4 to Elliptic Curve defined by y^2 + y = x^3 over Finite Field in alpha of size 2^4
        sage: phi_ker_poly = phi_v.kernel_polynomial()
        sage: phi_ker_poly
        x + 1
        sage: ker_poly_list = phi_ker_poly.list()
        sage: phi_k = EllipticCurveIsogeny(E, ker_poly_list)
        sage: phi_k == phi_v
        True
        sage: phi_k.rational_maps()
        ((x^3 + x + 1)/(x^2 + 1), (x^3*y + x^2*y + x*y + x + y)/(x^3 + x^2 + x + 1))
        sage: phi_v.rational_maps()
        ((x^3 + x + 1)/(x^2 + 1), (x^3*y + x^2*y + x*y + x + y)/(x^3 + x^2 + x + 1))
        sage: phi_k.degree() == phi_v.degree() == 3
        True
        sage: phi_k.is_separable()
        True
        sage: phi_v(E(0))
        (0 : 1 : 0)
        sage: alpha = E.base_field().gen()
        sage: Q = E((0, alpha*(alpha + 1)))
        sage: phi_v(Q)
        (1 : alpha^2 + alpha : 1)
        sage: phi_v(P) == phi_k(P)
        True
        sage: phi_k(P) == phi_v.codomain()(0)
        True

    We can create an isogeny that has kernel equal to the full 2
    torsion::

        sage: E = EllipticCurve(GF(3), [0,0,0,1,1])
        sage: ker_list = E.division_polynomial(2).list()
        sage: phi = EllipticCurveIsogeny(E, ker_list); phi
        Isogeny of degree 4 from Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 3 to Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 3
        sage: phi(E(0))
        (0 : 1 : 0)
        sage: phi(E((0,1)))
        (1 : 0 : 1)
        sage: phi(E((0,2)))
        (1 : 0 : 1)
        sage: phi(E((1,0)))
        (0 : 1 : 0)
        sage: phi.degree()
        4

    We can also create trivial isogenies with the trivial kernel::

        sage: E = EllipticCurve(GF(17), [11, 11, 4, 12, 10])
        sage: phi_v = EllipticCurveIsogeny(E, E(0))
        sage: phi_v.degree()
        1
        sage: phi_v.rational_maps()
        (x, y)
        sage: E == phi_v.codomain()
        True
        sage: P = E.random_point()
        sage: phi_v(P) == P
        True

        sage: E = EllipticCurve(GF(31), [23, 1, 22, 7, 18])
        sage: phi_k = EllipticCurveIsogeny(E, [1]); phi_k
        Isogeny of degree 1 from Elliptic Curve defined by y^2 + 23*x*y + 22*y = x^3 + x^2 + 7*x + 18 over Finite Field of size 31 to Elliptic Curve defined by y^2 + 23*x*y + 22*y = x^3 + x^2 + 7*x + 18 over Finite Field of size 31
        sage: phi_k.degree()
        1
        sage: phi_k.rational_maps()
        (x, y)
        sage: phi_k.codomain() == E
        True
        sage: phi_k.kernel_polynomial()
        1
        sage: P = E.random_point(); P == phi_k(P)
        True

    Velu and Kohel also work in characteristic 0::

        sage: E = EllipticCurve(QQ, [0,0,0,3,4])
        sage: P_list = E.torsion_points()
        sage: phi = EllipticCurveIsogeny(E, P_list); phi
        Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 3*x + 4 over Rational Field to Elliptic Curve defined by y^2 = x^3 - 27*x + 46 over Rational Field
        sage: P = E((0,2))
        sage: phi(P)
        (6 : -10 : 1)
        sage: phi_ker_poly = phi.kernel_polynomial()
        sage: phi_ker_poly
        x + 1
        sage: ker_poly_list = phi_ker_poly.list()
        sage: phi_k = EllipticCurveIsogeny(E, ker_poly_list); phi_k
        Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 3*x + 4 over Rational Field to Elliptic Curve defined by y^2 = x^3 - 27*x + 46 over Rational Field
        sage: phi_k(P) == phi(P)
        True
        sage: phi_k == phi
        True
        sage: phi_k.degree()
        2
        sage: phi_k.is_separable()
        True

    A more complicated example over the rationals (of odd degree)::

        sage: E = EllipticCurve('11a1')
        sage: P_list = E.torsion_points()
        sage: phi_v = EllipticCurveIsogeny(E, P_list); phi_v
        Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
        sage: P = E((16,-61))
        sage: phi_v(P)
        (0 : 1 : 0)
        sage: ker_poly = phi_v.kernel_polynomial(); ker_poly
        x^2 - 21*x + 80
        sage: ker_poly_list = ker_poly.list()
        sage: phi_k = EllipticCurveIsogeny(E, ker_poly_list); phi_k
        Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
        sage: phi_k == phi_v
        True
        sage: phi_v(P) == phi_k(P)
        True
        sage: phi_k.is_separable()
        True

    We can also do this same example over the number field defined by
    the irreducible two torsion polynomial of `E`::

        sage: E = EllipticCurve('11a1')
        sage: P_list = E.torsion_points()
        sage: K.<alpha> = NumberField(x^3 - 2* x^2 - 40*x - 158)
        sage: EK = E.change_ring(K)
        sage: P_list = [EK(P) for P in P_list]
        sage: phi_v = EllipticCurveIsogeny(EK, P_list); phi_v
        Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in alpha with defining polynomial x^3 - 2*x^2 - 40*x - 158 to Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-7820)*x + (-263580) over Number Field in alpha with defining polynomial x^3 - 2*x^2 - 40*x - 158
        sage: P = EK((alpha/2,-1/2))
        sage: phi_v(P)
        (122/121*alpha^2 + 1633/242*alpha - 3920/121 : -1/2 : 1)
        sage: ker_poly = phi_v.kernel_polynomial()
        sage: ker_poly
        x^2 - 21*x + 80
        sage: ker_poly_list = ker_poly.list()
        sage: phi_k = EllipticCurveIsogeny(EK, ker_poly_list)
        sage: phi_k
        Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in alpha with defining polynomial x^3 - 2*x^2 - 40*x - 158 to Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-7820)*x + (-263580) over Number Field in alpha with defining polynomial x^3 - 2*x^2 - 40*x - 158
        sage: phi_v == phi_k
        True
        sage: phi_k(P) == phi_v(P)
        True
        sage: phi_k == phi_v
        True
        sage: phi_k.degree()
        5
        sage: phi_v.is_separable()
        True

    The following example shows how to specify an isogeny from domain
    and codomain::

        sage: E = EllipticCurve('11a1')
        sage: R.<x> = QQ[]
        sage: f = x^2 - 21*x + 80
        sage: phi = E.isogeny(f)
        sage: E2 = phi.codomain()
        sage: phi_s = EllipticCurveIsogeny(E, None, E2, 5)
        sage: phi_s
        Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
        sage: phi_s == phi
        True
        sage: phi_s.rational_maps() == phi.rational_maps()
        True

    However only cyclic normalized isogenies can be constructed this
    way. So it won't find the isogeny [3]::

        sage: E.isogeny(None, codomain=E,degree=9)
        Traceback (most recent call last):
        ...
        ValueError: The two curves are not linked by a cyclic normalized isogeny of degree 9

    Also the presumed isogeny between the domain and codomain must be
    normalized::

        sage: E2.isogeny(None,codomain=E,degree=5)
        Traceback (most recent call last):
        ...
        ValueError: The two curves are not linked by a cyclic normalized isogeny of degree 5
        sage: phihat = phi.dual(); phihat
        Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
        sage: phihat.is_normalized()
        False

    Here an example of a construction of a endomorphisms with cyclic
    kernel on a CM-curve::

        sage: K.<i> = NumberField(x^2+1)
        sage: E = EllipticCurve(K, [1,0])
        sage: RK.<X> = K[]
        sage: f = X^2 - 2/5*i + 1/5
        sage: phi= E.isogeny(f)
        sage: isom = phi.codomain().isomorphism_to(E)
        sage: phi = isom * phi
        sage: phi.codomain() == phi.domain()
        True
        sage: phi.rational_maps()
        (((4/25*i + 3/25)*x^5 + (4/5*i - 2/5)*x^3 - x)/(x^4 + (-4/5*i + 2/5)*x^2 + (-4/25*i - 3/25)), ((11/125*i + 2/125)*x^6*y + (-23/125*i + 64/125)*x^4*y + (141/125*i + 162/125)*x^2*y + (3/25*i - 4/25)*y)/(x^6 + (-6/5*i + 3/5)*x^4 + (-12/25*i - 9/25)*x^2 + (2/125*i - 11/125)))

    Domain and codomain tests (see :trac:`12880`)::

        sage: E = EllipticCurve(QQ, [0,0,0,1,0])
        sage: phi = EllipticCurveIsogeny(E,  E(0,0))
        sage: phi.domain() == E
        True
        sage: phi.codomain()
        Elliptic Curve defined by y^2 = x^3 - 4*x over Rational Field

        sage: E = EllipticCurve(GF(31), [1,0,0,1,2])
        sage: phi = EllipticCurveIsogeny(E, [17, 1])
        sage: phi.domain()
        Elliptic Curve defined by y^2 + x*y = x^3 + x + 2 over Finite Field of size 31
        sage: phi.codomain()
        Elliptic Curve defined by y^2 + x*y = x^3 + 24*x + 6 over Finite Field of size 31

    Composition tests (see :trac:`16245`)::

        sage: E = EllipticCurve(j=GF(7)(0))
        sage: phi = E.isogeny([E(0), E((0,1)), E((0,-1))]); phi
        Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7
        sage: phi2 = phi * phi; phi2
        Composite map:
          From: Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7
          To:   Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7
          Defn:   Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7
                then
                  Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7

    Examples over relative number fields used not to work (see :trac:`16779`)::

        sage: pol26 = hilbert_class_polynomial(-4*26)
        sage: pol = NumberField(pol26,'a').optimized_representation()[0].polynomial()
        sage: K.<a> = NumberField(pol)
        sage: j = pol26.roots(K)[0][0]
        sage: E = EllipticCurve(j=j)
        sage: L.<b> = K.extension(x^2+26)
        sage: EL = E.change_ring(L)
        sage: iso2 = EL.isogenies_prime_degree(2); len(iso2)
        1
        sage: iso3 = EL.isogenies_prime_degree(3); len(iso3)
        2

    Examples over function fields used not to work (see :trac:`11327`)::

        sage: F.<t> = FunctionField(QQ)
        sage: E = EllipticCurve([0,0,0,-t^2,0])
        sage: isogs = E.isogenies_prime_degree(2)
        sage: isogs[0]
        Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + (-t^2)*x over Rational function field in t over Rational Field to Elliptic Curve defined by y^2 = x^3 + 4*t^2*x over Rational function field in t over Rational Field
        sage: isogs[0].rational_maps()
        ((x^2 - t^2)/x, (x^2*y + t^2*y)/x^2)
        sage: duals = [phi.dual() for phi in isogs]
        sage: duals[0]
        Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 4*t^2*x over Rational function field in t over Rational Field to Elliptic Curve defined by y^2 = x^3 + (-t^2)*x over Rational function field in t over Rational Field
        sage: duals[0].rational_maps()
        ((1/4*x^2 + t^2)/x, (1/8*x^2*y + (-1/2*t^2)*y)/x^2)
        sage: duals[0]
        Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 4*t^2*x over Rational function field in t over Rational Field to Elliptic Curve defined by y^2 = x^3 + (-t^2)*x over Rational function field in t over Rational Field
    """

    ####################
    # member variables
    ####################

    __check = None
    #
    # variables common to all algorithms
    #
    __E1 = None  # domain curve
    __E2 = None  # codomain curve

    __degree = None

    __algorithm = None

    __check = None
    #
    # pre isomorphism
    #
    __intermediate_domain = None
    __pre_isomorphism = None
    __prei_x_coord_ratl_map = None
    __prei_y_coord_ratl_map = None

    #
    # post isomorphism
    #

    __intermediate_codomain = None
    __post_isomorphism = None
    __posti_x_coord_ratl_map = None
    __posti_y_coord_ratl_map = None

    #
    # algebraic structs
    #
    __base_field = None
    __poly_ring = None # univariate in x over __base_field
    __mpoly_ring = None # bivariate in x, y over __base_field

    #
    # Rational Maps
    #
    __rational_maps_initialized = False
    __X_coord_rational_map = None
    __Y_coord_rational_map = None

    #
    # The dual
    #
    __dual = None

    #
    # Kernel Data
    #

    __kernel_list = None  # list of elements in the kernel

    __kernel_polynomial_list = None # polynomial stored as a little endian list of coefficients

    __kernel_polynomial = None # polynomial with roots at x values for x-coordinate of points in the kernel

    __inner_kernel_polynomial = None # the inner kernel polynomial (ignoring preisomorphism)

    __n = None


    #
    # member variables common to Velu's formula
    #

    # we keep track of the 2 torsion and non2torsion separately
    __kernel_2tor = None
    __kernel_non2tor = None

    # variables used in Velu's formula (as well as Kohel's variant)
    __v = None
    __w = None


    #
    # member variables specific to Kohel's algorithm.
    #

    __psi = None # psi polynomial
    __phi = None # phi polynomial
    __omega = None # omega polynomial


    #
    # Python Special Functions
    #

    def __init__(self, E, kernel, codomain=None, degree=None, model=None, check=True):
        r"""
        Constructor for EllipticCurveIsogeny class.

        EXAMPLES::

            sage: E = EllipticCurve(GF(2), [0,0,1,0,1])
            sage: phi = EllipticCurveIsogeny(E, [1,1]); phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 + y = x^3 + 1 over Finite Field of size 2 to Elliptic Curve defined by y^2 + y = x^3 over Finite Field of size 2

            sage: E = EllipticCurve(GF(31), [0,0,0,1,0])
            sage: P = E((2,17))
            sage: phi = EllipticCurveIsogeny(E, P); phi
            Isogeny of degree 8 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 31 to Elliptic Curve defined by y^2 = x^3 + 10*x + 28 over Finite Field of size 31

            sage: E = EllipticCurve('17a1')
            sage: phi = EllipticCurveIsogeny(E, [41/3, -55, -1, -1, 1]); phi
            Isogeny of degree 9 from Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 - x - 14 over Rational Field to Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 - 56*x - 10124 over Rational Field

            sage: E = EllipticCurve('37a1')
            sage: triv = EllipticCurveIsogeny(E, E(0)); triv
            Isogeny of degree 1 from Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: triv.rational_maps()
            (x, y)

            sage: E = EllipticCurve('49a3')
            sage: R.<X> = QQ[]
            sage: EllipticCurveIsogeny(E,X^3-13*X^2-58*X+503,check=False)
            Isogeny of degree 7 from Elliptic Curve defined by y^2 + x*y = x^3 - x^2 - 107*x + 552 over Rational Field to Elliptic Curve defined by y^2 + x*y = x^3 - x^2 - 5252*x - 178837 over Rational Field

        """
        if not is_EllipticCurve(E):
            raise ValueError("E parameter must be an EllipticCurve.")

        if not isinstance(kernel, list) and kernel in E:
            # a single point was given, we put it in a list
            # the first condition assures that [1,1] is treated as x+1
            kernel = [kernel]

        # if the kernel is None and the codomain isn't
        # calculate the kernel polynomial
        pre_isom = None
        post_isom = None

        self.__check = check

        if (kernel is None) and (codomain is not None):

            if (degree is None):
                raise ValueError("If specifying isogeny by domain and codomain, degree parameter must be set.")

            # save the codomain: really used now (trac #7096)
            old_codomain = codomain

            (pre_isom, post_isom, E, codomain, kernel) = compute_sequence_of_maps(E, codomain, degree)

        self.__init_algebraic_structs(E)

        algorithm = isogeny_determine_algorithm(E, kernel)

        self.__algorithm = algorithm

        if ("velu"==algorithm):
            self.__init_from_kernel_list(kernel)
        elif ("kohel"==algorithm):
            self.__init_from_kernel_polynomial(kernel)

        self.__compute_E2()

        self.__setup_post_isomorphism(codomain, model)

        if (pre_isom is not None):
            self._set_pre_isomorphism(pre_isom)

        if (post_isom is not None):
            self.__set_post_isomorphism(old_codomain, post_isom)   #(trac #7096)

        # Inheritance house keeping

        self.__perform_inheritance_housekeeping()

    def _eval(self, P):
        r"""
        Less strict evaluation method for internal use.

        In particular, this can be used to evaluate ``self`` at a
        point defined over an extension field.

        INPUT: a sequence of 3 coordinates defining a point on ``self``

        OUTPUT: the result of evaluating ``self'' at the given point

        EXAMPLES::

            sage: E = EllipticCurve([1,0]); E
            Elliptic Curve defined by y^2 = x^3 + x over Rational Field
            sage: phi = E.isogeny(E(0,0))
            sage: P = E.change_ring(QQbar).lift_x(QQbar.random_element())
            sage: phi._eval(P).curve()
            Elliptic Curve defined by y^2 = x^3 + (-4)*x over Algebraic Field

        ::

            sage: E = EllipticCurve(j=Mod(0,419))
            sage: K = next(filter(bool, E(0).division_points(5)))
            sage: psi = E.isogeny(K)
            sage: Ps = E.change_ring(GF(419**2))(0).division_points(5)
            sage: {psi._eval(P).curve() for P in Ps}
            {Elliptic Curve defined by y^2 = x^3 + 140*x + 214 over Finite Field in z2 of size 419^2}
        """
        if self._domain.defining_polynomial()(*P):
            raise ValueError(f'{P} not on {self._domain}')

        if not P:
            k = Sequence(tuple(P)).universe()
            return self._codomain(0).change_ring(k)

        Q = P.xy()
        pre_iso = self.get_pre_isomorphism()
        if pre_iso is not None:
            Q = baseWI.__call__(pre_iso, Q)

        if self.__algorithm == 'velu':
            compute = self.__compute_via_velu
        elif self.__algorithm == 'kohel':
            compute = self.__compute_via_kohel
        else:
            raise NotImplementedError

        try:
            Q = compute(*Q)
        except ZeroDivisionError:
            Q = (0,1,0)

        post_iso = self.get_post_isomorphism()
        if post_iso is not None:
            Q = baseWI.__call__(post_iso, Q)

        k = Sequence(tuple(P) + tuple(Q)).universe()
        return self._codomain.base_extend(k).point(Q)


    def _call_(self, P):
        r"""
        Function that implements the call-ability of elliptic curve
        isogenies.

        EXAMPLES::

            sage: E = EllipticCurve(GF(17), [1, 9, 5, 4, 3])
            sage: phi = EllipticCurveIsogeny(E, [6,13,1])
            sage: phi(E((1,0)))
            (15 : 13 : 1)

            sage: E = EllipticCurve(GF(23), [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E, E((0,0)))
            sage: phi(E((1,5)))
            (2 : 0 : 1)

            sage: E = EllipticCurve(QQ, [0,0,0,3,0])
            sage: P = E((1,2))
            sage: phi = EllipticCurveIsogeny(E, [0,1])
            sage: phi(P)
            (4 : -4 : 1)
            sage: phi(-P)
            (4 : 4 : 1)

            sage: E = EllipticCurve(GF(17), [0,-1,0,-3,-1])
            sage: Q = E((16,0))
            sage: tau = E.isogeny([Q],E)
            sage: tau(Q)
            (0 : 1 : 0)

        TESTS:

        Tests for :trac:`10888`::

            sage: K.<th> = NumberField(x^2+3)
            sage: E = EllipticCurve(K,[7,0])
            sage: phi = E.isogeny(E(0,0))
            sage: P = E(-3,4*th)
            sage: phi(P)
            (-16/3 : 8/9*th : 1)
            sage: Q = phi(P)
            sage: phihat = phi.dual()
            sage: phihat(Q)
            (-1/48 : 127/576*th : 1)

        Call a composed isogeny (added for :trac:`16238`)::

            sage: E = EllipticCurve(j=GF(7)(0))
            sage: phi = E.isogeny([E(0), E((0,1)), E((0,-1))])
            sage: phi(E.points()[0])
            (0 : 1 : 0)
            sage: phi2 = phi * phi
            sage: phi2(E.points()[0])
            (0 : 1 : 0)

        Coercion works fine with :meth:`_call_` (added for :trac:`16238`)::

            sage: K.<th> = NumberField(x^2+3)
            sage: E = EllipticCurve(K,[7,0])
            sage: E2=EllipticCurve(K,[5,0])
            sage: phi=E.isogeny(E(0))
            sage: phi(E2(0))
            (0 : 1 : 0)
            sage: E2(20,90)
            (20 : 90 : 1)
            sage: phi(E2(20,90))
            Traceback (most recent call last):
            ...
            TypeError: (20 : 90 : 1) fails to convert into the map's domain Elliptic Curve defined by y^2 = x^3 + 7*x over Number Field in th with defining polynomial x^2 + 3, but a `pushforward` method is not properly implemented

        """
        if P.is_zero():
            return self.__E2(0)

        (xP, yP) = P.xy()
        # if there is a pre isomorphism, apply it
        if self.__pre_isomorphism is not None:
            temp_xP = self.__prei_x_coord_ratl_map(xP)
            temp_yP = self.__prei_y_coord_ratl_map(xP, yP)
            (xP, yP) = (temp_xP, temp_yP)

        if "velu" == self.__algorithm:
            outP = self.__compute_via_velu_numeric(xP, yP)
        elif "kohel" == self.__algorithm:
            outP = self.__compute_via_kohel_numeric(xP,yP)

        # the intermediate functions return the point at infinity
        # if the input point is in the kernel
        if (outP == self.__intermediate_codomain(0)):
            return self.__E2(0)

        # if there is a post isomorphism, apply it
        if (self.__post_isomorphism is not None):
            tempX = self.__posti_x_coord_ratl_map(outP[0])
            tempY = self.__posti_y_coord_ratl_map(outP[0], outP[1])
            outP = (tempX, tempY)

        return self.__E2(outP)

    def __getitem__(self, i):
        r"""
        Return one of the rational map components.

        .. NOTE::

           Both components are returned as elements of the function
           field `F(x,y)` in two variables over the base field `F`,
           though the first only involves `x`.  To obtain the
           `x`-coordinate function as a rational function in `F(x)`,
           use :meth:`x_rational_map`.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,2,0,1,-1])
            sage: phi = EllipticCurveIsogeny(E, [1])
            sage: phi[0]
            x
            sage: phi[1]
            y

            sage: E = EllipticCurve(GF(17), [0,0,0,3,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi[0]
            (x^2 + 3)/x
            sage: phi[1]
            (x^2*y - 3*y)/x^2
        """
        return self.rational_maps()[i]

    def __iter__(self):
        r"""
        Return an iterator through the rational map components.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,2,0,1,-1])
            sage: phi = EllipticCurveIsogeny(E, [1])
            sage: for c in phi: print(c)
            x
            y

            sage: E = EllipticCurve(GF(17), [0,0,0,3,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: for c in phi: print(c)
            (x^2 + 3)/x
            (x^2*y - 3*y)/x^2
        """
        return iter(self.rational_maps())

    def __neg__(self):
        r"""
        Function to implement unary negation (-) operator on
        isogenies. Return a copy of this isogeny that has been
        negated.

        EXAMPLES:

        The following examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(17)(0))
            sage: phi = EllipticCurveIsogeny(E,  E((-1,0)) )
            sage: negphi = -phi
            sage: phi(E((0,1))) + negphi(E((0,1))) == 0
            True

            sage: E = EllipticCurve(j=GF(19)(1728))
            sage: R.<x> = GF(19)[]
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: negphi = -phi
            sage: phi(E((3,7))) + negphi(E((3,12))) == phi(2*E((3,7)))
            True
            sage: negphi(E((18,6)))
            (17 : 0 : 1)

            sage: R.<x> = QQ[]
            sage: E = EllipticCurve('17a1')
            sage: R.<x> = QQ[]
            sage: f = x - 11/4
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: negphi = -phi
            sage: phi.rational_maps()[0] == negphi.rational_maps()[0]
            True
            sage: P = E((7,13))
            sage: phi(P) + negphi(P) == 0
            True

        """
        output = copy(self)
        E2 = output.__E2
        iso = WeierstrassIsomorphism(E2, (-1,0,-E2.a1(),-E2.a3()))
        output._set_post_isomorphism(iso)
        return output

    #
    # Sage Special Functions
    #

    def _repr_(self):
        r"""
        Special sage specific function that implement the
        functionality to display the isogeny self as a string.

        EXAMPLES::

            sage: E = EllipticCurve(GF(31), [1,0,1,1,0])
            sage: phi = EllipticCurveIsogeny(E, E((0,0)) )
            sage: phi._repr_()
            'Isogeny of degree 17 from Elliptic Curve defined by y^2 + x*y + y = x^3 + x over Finite Field of size 31 to Elliptic Curve defined by y^2 + x*y + y = x^3 + 14*x + 9 over Finite Field of size 31'

            sage: E = EllipticCurve(QQ, [1,0,0,1,9])
            sage: phi = EllipticCurveIsogeny(E, [2,1])
            sage: phi._repr_()
            'Isogeny of degree 2 from Elliptic Curve defined by y^2 + x*y = x^3 + x + 9 over Rational Field to Elliptic Curve defined by y^2 + x*y = x^3 - 59*x + 165 over Rational Field'

        """
        return 'Isogeny of degree %r from %r to %r' % (
                self.__degree, self.__E1, self.__E2)

    def _latex_(self):
        r"""
        Special sage specific function that implements functionality
        to display an isogeny object as a latex string.

        This function returns a latex string representing the isogeny
        self as the `x` and `y` coordinate rational functions.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,0,0,1,-1])
            sage: phi = EllipticCurveIsogeny(E, E(0))
            sage: phi._latex_()
            '\\left( x , y \\right)'

            sage: E = EllipticCurve(GF(17), [0,0,0,1,-1])
            sage: R.<X> = GF(17)[]
            sage: phi = EllipticCurveIsogeny(E, X+11)
            sage: phi._latex_()
            '\\left( \\frac{x^{2} + 11 x + 7}{x + 11} , \\frac{x^{2} y + 5 x y + 12 y}{x^{2} + 5 x + 2} \\right)'


        """
        ratl_maps = self.rational_maps()
        return '\\left( %s , %s \\right)' % (ratl_maps[0]._latex_(), ratl_maps[1]._latex_())


    ###########################
    # Private Common Functions
    ###########################

    # delete the hash value
    def __clear_cached_values(self):
        r"""
        A private function to clear the hash if the codomain has been
        modified by a pre or post isomorphism.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: E = EllipticCurve(QQ, [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: old_ratl_maps = phi.rational_maps()
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: phi.set_post_isomorphism(WeierstrassIsomorphism(phi.codomain(), (-1,0,0,0)))
            ...
            sage: old_ratl_maps == phi.rational_maps()
            False
            sage: old_ratl_maps[1] == -phi.rational_maps()[1]
            True

            sage: F = GF(127); R.<x> = F[]
            sage: E = EllipticCurve(j=F(1728))
            sage: f = x^5 + 43*x^4 + 97*x^3 + 81*x^2 + 42*x + 82
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: old_ratl_maps = phi.rational_maps()
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: phi.set_post_isomorphism(WeierstrassIsomorphism(phi.codomain(), (-13,13,-13,13)))
            ...
            sage: old_ratl_maps == phi.rational_maps()
            False
            sage: phi._EllipticCurveIsogeny__clear_cached_values()

        """
        self.__rational_maps_initialized = False
        self.__X_coord_rational_map = None
        self.__Y_coord_rational_map = None
        self.__dual = None

    # performs the inheritance house keeping
    def __perform_inheritance_housekeeping(self):
        r"""
        Internal helper function, sets values on the super classes of
        this class.

        EXAMPLES:

        The following examples will implicitly exercise this
        function::

            sage: E = EllipticCurve(GF(43), [2,3,5,7,11])
            sage: R.<x> = GF(43)[]; f = x + 42
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi._EllipticCurveIsogeny__perform_inheritance_housekeeping()
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E2 = phi.codomain()
            sage: post_isom = WeierstrassIsomorphism(E2, (41, 37, 31, 29))
            sage: phi.set_post_isomorphism(post_isom)
            ...
            sage: E1pr = WeierstrassIsomorphism(E, (-1, 2, -3, 4)).codomain()
            sage: pre_isom = E1pr.isomorphism_to(E)
            sage: phi.set_pre_isomorphism(pre_isom)
            ...

        """
        # one of the superclasses uses these fields
        self._domain = self.__E1
        self._codomain = self.__E2

        # sets up the parent
        EllipticCurveHom.__init__(self, self._domain, self._codomain)

    def __init_algebraic_structs(self, E):
        r"""
        An internal function for EllipticCurveIsogeny objects that
        sets up the member variables necessary for algebra.

        EXAMPLES::

            sage: E = EllipticCurve(j=GF(17)(0))
            sage: phi = EllipticCurveIsogeny(E,  E((-1,0)))

        The constructor calls this function itself, so the fields it
        sets are already defined::

            sage: phi._EllipticCurveIsogeny__E1
            Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 17
            sage: phi._EllipticCurveIsogeny__base_field
            Finite Field of size 17
            sage: phi._EllipticCurveIsogeny__poly_ring
            Univariate Polynomial Ring in x over Finite Field of size 17
            sage: phi._EllipticCurveIsogeny__mpoly_ring
            Multivariate Polynomial Ring in x, y over Finite Field of size 17
            sage: phi._EllipticCurveIsogeny__intermediate_domain
            Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 17

        Now, calling the initialization function does nothing more::

            sage: phi._EllipticCurveIsogeny__init_algebraic_structs(E)
            sage: phi._EllipticCurveIsogeny__E1
            Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 17
            sage: phi._EllipticCurveIsogeny__base_field
            Finite Field of size 17
            sage: phi._EllipticCurveIsogeny__poly_ring
            Univariate Polynomial Ring in x over Finite Field of size 17
            sage: phi._EllipticCurveIsogeny__mpoly_ring
            Multivariate Polynomial Ring in x, y over Finite Field of size 17
            sage: phi._EllipticCurveIsogeny__intermediate_domain
            Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 17

            sage: E = EllipticCurve(QQ, [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi._EllipticCurveIsogeny__init_algebraic_structs(E)
            sage: phi._EllipticCurveIsogeny__E1
            Elliptic Curve defined by y^2 = x^3 + x over Rational Field
            sage: phi._EllipticCurveIsogeny__base_field
            Rational Field
            sage: phi._EllipticCurveIsogeny__poly_ring
            Univariate Polynomial Ring in x over Rational Field
            sage: phi._EllipticCurveIsogeny__mpoly_ring
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: phi._EllipticCurveIsogeny__intermediate_domain
            Elliptic Curve defined by y^2 = x^3 + x over Rational Field

            sage: F = GF(19); R.<x> = F[]
            sage: E = EllipticCurve(j=GF(19)(0))
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: phi._EllipticCurveIsogeny__init_algebraic_structs(E)
            sage: phi._EllipticCurveIsogeny__E1
            Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 19
            sage: phi._EllipticCurveIsogeny__base_field
            Finite Field of size 19
            sage: phi._EllipticCurveIsogeny__poly_ring
            Univariate Polynomial Ring in x over Finite Field of size 19
            sage: phi._EllipticCurveIsogeny__mpoly_ring
            Multivariate Polynomial Ring in x, y over Finite Field of size 19
            sage: phi._EllipticCurveIsogeny__intermediate_domain
            Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 19
        """
        self.__E1 = E
        self.__base_field = E.base_ring()
        self.__poly_ring = PolynomialRing(self.__base_field, ['x'])
        self.__mpoly_ring = PolynomialRing(self.__base_field, ['x','y'])
        from sage.rings.all import FractionField
        self.__xfield = FractionField(self.__poly_ring)
        self.__xyfield = FractionField(self.__mpoly_ring)
        self.__intermediate_domain = E

    def __compute_E2(self):
        r"""
        Private function that computes and sets the isogeny codomain.

        EXAMPLES:

        These examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(7)(1728))
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.codomain()
            Elliptic Curve defined by y^2 = x^3 + 3*x over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__compute_E2()

            sage: R.<x> = GF(7)[]
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: phi.codomain()
            Elliptic Curve defined by y^2 = x^3 + 3*x over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__compute_E2()

        """

        if ("velu" == self.__algorithm):
            E2 = self.__compute_E2_via_velu()
        elif ("kohel" == self.__algorithm):
            E2 = self.__compute_E2_via_kohel()

        self.__E2 = E2
        self.__intermediate_codomain = E2

    # initializes the rational maps fields
    def __initialize_rational_maps(self, precomputed_maps=None):
        r"""
        Private function that computes and initializes the rational
        maps.

        INPUT:

        - ``precomputed_maps`` (default None) -- tuple (X,Y) of
          rational functions in x,y

        EXAMPLES:

        The following examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(7)(1728))
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi._EllipticCurveIsogeny__initialize_rational_maps()
            sage: phi.rational_maps()
            ((x^2 + 1)/x, (x^2*y - y)/x^2)

            sage: R.<x> = GF(7)[]
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: phi.rational_maps()
            ((x^2 + 1)/x, (x^2*y - y)/x^2)
            sage: phi._EllipticCurveIsogeny__initialize_rational_maps()

            sage: E = EllipticCurve([1,2,3,4,5])
            sage: Eshort = E.short_weierstrass_model()
            sage: phi = E.isogeny(E(0), Eshort)
            sage: phiX, phiY = phi.rational_maps()
            sage: phiX(1,2), phiY(1,2)
            (63, 864)
        """
        if self.__rational_maps_initialized:
            return

        if precomputed_maps is None:
            if ("velu"==self.__algorithm):
                (X_map, Y_map) = self.__initialize_rational_maps_via_velu()

            if ("kohel"==self.__algorithm):
                (X_map, Y_map) = self.__initialize_rational_maps_via_kohel()
        else:
            X_map, Y_map = precomputed_maps
            # cannot coerce directly in xfield for some reason
            X_map = self.__poly_ring(X_map.numerator())/self.__poly_ring(X_map.denominator())

        if self.__prei_x_coord_ratl_map is not None:
            prei_X_map = self.__prei_x_coord_ratl_map
            prei_Y_map = self.__prei_y_coord_ratl_map
            X_map = X_map(prei_X_map)
            Y_map = Y_map([prei_X_map, prei_Y_map])

        if self.__posti_x_coord_ratl_map is not None:
            # Do not reverse the order here!
            Y_map = self.__posti_y_coord_ratl_map([X_map, Y_map])
            X_map = self.__posti_x_coord_ratl_map(X_map)

        self.__X_coord_rational_map = self.__xfield(X_map)
        self.__Y_coord_rational_map = self.__xyfield(Y_map)
        self.__rational_maps_initialized = True


    def __init_kernel_polynomial(self):
        r"""
        Private function that initializes the kernel polynomial (if
        the algorithm does not take it as a parameter).

        EXAMPLES:

        The following examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(7)(1728))
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.kernel_polynomial()
            x
            sage: phi._EllipticCurveIsogeny__init_kernel_polynomial()
            [0, 1]
        """

        if (self.__kernel_polynomial_list is not None):
            return self.__kernel_polynomial_list

        if ("velu" == self.__algorithm):
            ker_poly_list = self.__init_kernel_polynomial_velu()
        else:
            raise ValueError("The kernel polynomial should already be defined!")

        return ker_poly_list

    def __set_pre_isomorphism(self, domain, isomorphism):
        r"""
        Private function to set the pre isomorphism and domain (and
        keep track of the domain of the isogeny).

        EXAMPLES::

            sage: E = EllipticCurve(GF(43), [2,3,5,7,11])
            sage: R.<x> = GF(43)[]; f = x + 42
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi._EllipticCurveIsogeny__perform_inheritance_housekeeping()
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E1pr = WeierstrassIsomorphism(E, (-1, 2, -3, 4)).codomain()
            sage: pre_isom = E1pr.isomorphism_to(E)
            sage: phi.set_pre_isomorphism(pre_isom)
            ...
            sage: phi._EllipticCurveIsogeny__set_pre_isomorphism(E, WeierstrassIsomorphism(E, (-1, 3, -3, 4)))
            sage: E == phi.domain()
            True

        """

        self.__E1 = domain

        # set the isomorphism
        self.__pre_isomorphism = isomorphism

        # calculate the isomorphism as a rational map.

        u, r, s, t = [self.__base_field(c) for c in isomorphism.tuple()]
        uinv = 1/u
        uinv2 = uinv**2
        uinv3 = uinv*uinv2

        x = self.__poly_ring.gen()
        y = self.__xyfield.gen(1) # not mpoly_ring.gen(1) else we end
                                  # up in K(x)[y] and trouble ensues

        self.__prei_x_coord_ratl_map = (x - r) * uinv2
        self.__prei_y_coord_ratl_map = (y - s*(x-r) - t) * uinv3

        if (self.__kernel_polynomial is not None):
            ker_poly = self.__kernel_polynomial
            ker_poly = ker_poly(self.__prei_x_coord_ratl_map)
            self.__kernel_polynomial = ker_poly.monic()

        self.__perform_inheritance_housekeeping()


    def __set_post_isomorphism(self, codomain, isomorphism):
        r"""
        Private function to set the post isomorphism and codomain (and
        keep track of the codomain of the isogeny).

        EXAMPLES:

        The following examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(7)(1728))
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E2 = phi.codomain()
            sage: isom = WeierstrassIsomorphism(E2, (-1,2,-3,4))
            sage: phi.set_post_isomorphism(isom)
            ...
            sage: phi._EllipticCurveIsogeny__set_post_isomorphism(E2, WeierstrassIsomorphism(phi.codomain(), (1,-2,3,-4)))
            sage: E2 == phi.codomain()
            True

        """

        # set the codomains
        self.__E2 = codomain

        # set the isomorphism
        self.__post_isomorphism = isomorphism

        # calculate the isomorphism as a rational map.

        u, r, s, t = [self.__base_field(c) for c in isomorphism.tuple()]
        uinv = 1/u
        uinv2 = uinv**2
        uinv3 = uinv*uinv2

        x = self.__poly_ring.gen()
        y = self.__xyfield.gen(1)

        self.__posti_x_coord_ratl_map = (x - r) * uinv2
        self.__posti_y_coord_ratl_map = (y - s*(x-r) - t) * uinv3

        self.__perform_inheritance_housekeeping()

    def __setup_post_isomorphism(self, codomain, model):
        r"""
        Private function to set up the post isomorphism given the
        codomain.

        EXAMPLES:

        The following examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(7)(1728))
            sage: E2 = EllipticCurve(GF(7), [0,0,0,5,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)), E2); phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 5*x over Finite Field of size 7
            sage: E3 = EllipticCurve(GF(7), [0,0,0,6,0])
            sage: phi._EllipticCurveIsogeny__setup_post_isomorphism(E3, None)
            sage: phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7

            sage: R.<x> = QQ[]
            sage: E = EllipticCurve(j=1728)
            sage: f = x^3 - x
            sage: phi = EllipticCurveIsogeny(E, f, model='minimal'); phi
            Isogeny of degree 4 from Elliptic Curve defined by y^2 = x^3 - x over Rational Field to Elliptic Curve defined by y^2 = x^3 - x over Rational Field

            sage: phi = EllipticCurveIsogeny(E, f, model=None)
            sage: phi._EllipticCurveIsogeny__setup_post_isomorphism(None, 'minimal')
            sage: phi
            Isogeny of degree 4 from Elliptic Curve defined by y^2 = x^3 - x over Rational Field to Elliptic Curve defined by y^2 = x^3 - x over Rational Field

        """
        # TODO: add checks to make sure that codomain and model
        # parameters are consistent with the algorithm used.

        post_isom = None
        newE2 = None

        oldE2 = self.__E2

        if (model is not None):

            if (codomain is not None):
                raise ValueError("Cannot specify a codomain and model flag simultaneously.")

            if ('minimal' == model):

                if (not is_NumberField(oldE2.base_field())):
                    raise ValueError("specifying minimal for model flag only valid with curves over number fields.")

                newE2 = oldE2.global_minimal_model(semi_global=True)
                post_isom = oldE2.isomorphism_to(newE2)

            else:
                raise ValueError("Unknown value of model flag.")

        elif (codomain is not None):
            if (not is_EllipticCurve(codomain)):
                raise ValueError("Codomain parameter must be an elliptic curve.")

            if (not oldE2.is_isomorphic(codomain)):
                raise ValueError("Codomain parameter must be isomorphic to computed codomain isogeny")

            newE2 = codomain
            post_isom = oldE2.isomorphism_to(newE2)

        if (post_isom is not None):
            self.__set_post_isomorphism(newE2, post_isom)

        return

    ###########################
    # Velu's Formula Functions
    ###########################

    #
    # Setup function for Velu's formula
    #

    def __init_from_kernel_list(self, kernel_gens):
        r"""
        Private function that initializes the isogeny from a list of
        points which generate the kernel (For Velu's formulas.)

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0))); phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 4*x over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__init_from_kernel_list([E(0), E((0,0))])

        The following example demonstrates the necessity of avoiding any calls
        to P.order(), since such calls involve factoring the group order which
        could take a long time. ::

            sage: p = 12 * next_prime(2^180) * next_prime(2^194) - 1
            sage: F = FiniteField(p, proof=False)
            sage: E = EllipticCurve([F(1), F(0)])
            sage: P = E(0).division_points(3)[1]
            sage: EllipticCurveIsogeny(E, P)
            Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 461742260113997803268895001173557974278278194575766957660028841364655249961609425998827452443620996655395008156411 to Elliptic Curve defined by y^2 = x^3 + 80816485163488178037199320944019099858815874115367810482828676054000067654558381377552245721755005198633191074893*x + 301497584865165444049833326660609767433467459033532853758006118022998267706948164646650354324860226263546558337993 over Finite Field of size 461742260113997803268895001173557974278278194575766957660028841364655249961609425998827452443620996655395008156411

        """
        if self.__check :
            for P in kernel_gens:
                if not P.has_finite_order():
                    raise ValueError("The points in the kernel must be of finite order.")

        # Compute a list of points in the subgroup generated by the
        # points in kernel_gens.  This is very naive: when finite
        # subgroups are implemented better, this could be simplified,
        # but it won't speed things up too much.

        kernel_set = Set([self.__E1(0)])
        from sage.misc.all import flatten

        def all_multiples(itr, terminal):
            mult_list = [terminal]
            R = terminal + itr
            while R != terminal:
                mult_list.append(R)
                R = R + itr
            return mult_list
        for P in kernel_gens:
            kernel_set += Set(flatten([all_multiples(P,Q)
                                       for Q in kernel_set]))
        self.__kernel_list = kernel_set.list()
        self.__kernel_2tor = {}
        self.__kernel_non2tor = {}
        self.__degree = Integer(len(kernel_set))
        self.__sort_kernel_list()

    #
    # Precompute the values in Velu's Formula.
    #
    def __sort_kernel_list(self):
        r"""
        Private function that sorts the list of points in the kernel
        (For Velu's formulas). Sorts out the 2 torsion points, and
        puts them in a dictionary.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P); phi
            Isogeny of degree 4 from Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__kernel_2tor = {}
            sage: phi._EllipticCurveIsogeny__kernel_non2tor = {}
            sage: phi._EllipticCurveIsogeny__sort_kernel_list()

        """

        a1,a2,a3,a4,a6 = self.__E1.ainvs()

        v = 0
        w = 0

        for Q in self.__kernel_list:

            if Q.is_zero():
                continue

            (xQ,yQ) = Q.xy()

            gxQ = 3*xQ**2 + 2*a2*xQ + a4 - a1*yQ
            gyQ = -2*yQ - a1*xQ - a3

            uQ = gyQ**2

            # sort torsion points:
            if (2*yQ == -a1*xQ - a3): # Q is 2-torsion
                vQ = gxQ
                self.__kernel_2tor[xQ] = (xQ,yQ,gxQ,gyQ,vQ,uQ)
                v = v + vQ
                w = w + (uQ + xQ*vQ)
            elif xQ not in self.__kernel_non2tor: # Q is not a 2-torsion
                vQ = 2*gxQ - a1*gyQ
                self.__kernel_non2tor[xQ] = (xQ,yQ,gxQ,gyQ,vQ,uQ)
                v = v + vQ
                w = w + (uQ + xQ*vQ)

        self.__v = v
        self.__w = w

    #
    # Velu's formula computing the codomain curve
    #
    def __compute_E2_via_velu(self):
        r"""
        Private function that computes the codomain via Velu's
        formulas.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P)
            sage: phi.codomain()
            Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__compute_E2_via_velu()
            Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 7

        """
        v = self.__v
        w = self.__w

        return compute_codomain_formula(self.__E1, v,w)


    def __velu_sum_helper(self, Qvalues, a1, a3, x, y):
        r"""
        Private function for Velu's formulas, helper function to help
        perform the summation.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P)
            sage: Q = E((0,0)); phi(Q)
            (0 : 0 : 1)
            sage: phi.rational_maps()
            ((x^4 - 2*x^3 + x^2 - 3*x)/(x^3 - 2*x^2 + 3*x - 2), (x^5*y - 2*x^3*y - x^2*y - 2*x*y + 2*y)/(x^5 + 3*x^3 + 3*x^2 + x - 1))

            sage: F = GF(7)
            sage: E = EllipticCurve(F, [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)) )
            sage: Qvals = phi._EllipticCurveIsogeny__kernel_2tor[0]
            sage: phi._EllipticCurveIsogeny__velu_sum_helper(Qvals, 0, 0, F(5), F(5))
            (3, 3)
            sage: R.<x,y> = GF(7)[]
            sage: phi._EllipticCurveIsogeny__velu_sum_helper(Qvals, 0, 0, x, y)
            (1/x, y/x^2)

        """
        xQ = Qvalues[0]
        yQ = Qvalues[1]
        gxQ = Qvalues[2]
        gyQ = Qvalues[3]
        vQ = Qvalues[4]
        uQ = Qvalues[5]

        t1 = x - xQ
        inv_t1 = t1**-1
        inv_t1_2 = inv_t1**2
        inv_t1_3 = inv_t1_2*inv_t1

        tX = vQ*inv_t1 + uQ*(inv_t1_2)

        tY0 = uQ*(2*y + a1*x + a3)
        tY1 = vQ*(a1*t1 + y - yQ)
        tY2 = a1*uQ - gxQ*gyQ

        # Without this explicit coercion, tY ends up in K(x)[y]
        # instead of K(x,y), and trouble ensues!
        from sage.rings.all import FractionField
        F = FractionField(y.parent())
        tY =  ( tY0*F(inv_t1_3) + (tY1 + tY2)*F(inv_t1_2) )

        return (tX, tY)


    def __compute_via_velu_numeric(self, xP, yP):
        r"""
        Private function that sorts the list of points in the kernel
        (for Velu's formulas). Sorts out the 2 torsion points, and
        puts them in a dictionary.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: F = GF(7)
            sage: E = EllipticCurve(F, [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P)
            sage: Q = E((0,0)); phi(Q)
            (0 : 0 : 1)
            sage: Q = E((-1,0)); phi(Q)
            (0 : 0 : 1)
            sage: phi._EllipticCurveIsogeny__compute_via_velu_numeric(F(0), F(0))
            (0, 0)
            sage: phi._EllipticCurveIsogeny__compute_via_velu_numeric(F(-1), F(0))
            (0, 0)

        """
        # first check if the point is in the kernel
        if xP in self.__kernel_2tor or xP in self.__kernel_non2tor:
            return self.__intermediate_codomain(0)

        outP = self.__compute_via_velu(xP,yP)

        return outP


    def __compute_via_velu(self, xP, yP):
        r"""
        Private function for Velu's formulas, to perform the summation.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: F = GF(7)
            sage: E = EllipticCurve(F, [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P)
            sage: Q = E((0,0)); phi(Q)
            (0 : 0 : 1)
            sage: phi.rational_maps()
            ((x^4 - 2*x^3 + x^2 - 3*x)/(x^3 - 2*x^2 + 3*x - 2), (x^5*y - 2*x^3*y - x^2*y - 2*x*y + 2*y)/(x^5 + 3*x^3 + 3*x^2 + x - 1))
            sage: phi._EllipticCurveIsogeny__compute_via_velu(F(0), F(0))
            (0, 0)
            sage: R.<x,y> = GF(7)[]
            sage: phi._EllipticCurveIsogeny__compute_via_velu(x, y)
            ((x^4 - 2*x^3 + x^2 - 3*x)/(x^3 - 2*x^2 + 3*x - 2),
             (x^5*y - 2*x^3*y - x^2*y - 2*x*y + 2*y)/(x^5 + 3*x^3 + 3*x^2 + x - 1))
        """
        ker_2tor = self.__kernel_2tor
        ker_non2tor = self.__kernel_non2tor

        X = 0
        Y = 0

        a1 = self.__E1.a1()
        a3 = self.__E1.a3()

        # next iterate over the 2torsion points of the kernel
        for Qvalues in ker_2tor.values():
            (tX, tY) = self.__velu_sum_helper(Qvalues, a1, a3, xP, yP)
            X = X + tX
            Y = Y + tY

        for Qvalues in ker_non2tor.values():
            (tX, tY) = self.__velu_sum_helper(Qvalues, a1, a3, xP, yP)
            X = X + tX
            Y = Y + tY

        X = xP + X
        Y = yP - Y

        return (X,Y)


    def __initialize_rational_maps_via_velu(self):
        r"""
        Private function for Velu's formulas, helper function to
        initialize the rational maps.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P)
            sage: phi.rational_maps()
            ((x^4 - 2*x^3 + x^2 - 3*x)/(x^3 - 2*x^2 + 3*x - 2), (x^5*y - 2*x^3*y - x^2*y - 2*x*y + 2*y)/(x^5 + 3*x^3 + 3*x^2 + x - 1))
            sage: phi._EllipticCurveIsogeny__initialize_rational_maps_via_velu()
            ((x^4 + 5*x^3 + x^2 + 4*x)/(x^3 + 5*x^2 + 3*x + 5), (x^5*y - 2*x^3*y - x^2*y - 2*x*y + 2*y)/(x^5 + 3*x^3 + 3*x^2 + x - 1))
        """
        x = self.__poly_ring.gen()
        y = self.__mpoly_ring.gen(1)

        return self.__compute_via_velu(x,y)


    def __init_kernel_polynomial_velu(self):
        r"""
        Private function for Velu's formulas, helper function to
        initialize the rational maps.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P)
            sage: phi.kernel_polynomial()
            x^2 + 2*x + 4
            sage: phi._EllipticCurveIsogeny__init_kernel_polynomial_velu()
            [4, 2, 1]
        """
        poly_ring = self.__poly_ring
        x = poly_ring.gen()

        invX = 0

        if (self.__pre_isomorphism is not None):
            pre_isom = self.__pre_isomorphism
            u = pre_isom.u
            r = pre_isom.r
            invX = (u**2)*x + r
        else:
            invX = x

        psi = poly_ring(1)

        for Qvalues in self.__kernel_2tor.values():
            xQ = invX(x=Qvalues[0])
            psi = psi*(x - xQ)

        for Qvalues in self.__kernel_non2tor.values():
            xQ = invX(x=Qvalues[0])
            psi = psi*(x - xQ)

        ker_poly_list = psi.list()

        self.__kernel_polynomial_list = ker_poly_list
        self.__kernel_polynomial = psi

        return ker_poly_list



    ###################################
    # Kohel's Variant of Velu's Formula
    ###################################

    def __init_from_kernel_polynomial(self, kernel_polynomial):
        r"""
        Private function that initializes the isogeny from a kernel
        polynomial.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: phi = EllipticCurveIsogeny(E, x);phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 4*x over Finite Field of size 7

            sage: phi._EllipticCurveIsogeny__init_from_kernel_polynomial(x)

            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3); phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 1 over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 4*x + 2 over Finite Field of size 7

            sage: phi._EllipticCurveIsogeny__init_from_kernel_polynomial(x+6)

        """
        poly_ring = self.__poly_ring
        E = self.__E1

        # Convert to a univariate polynomial, even if it had a
        # bivariate parent, or was given as a list:
        self.__kernel_polynomial = psi = poly_ring(kernel_polynomial)

        if psi.leading_coefficient() != 1:
            raise ValueError("The kernel polynomial must be monic.")

        self.__kernel_polynomial_list = psi.list()

        #
        # Determine if kernel polynomial is entirely a two torsion
        #
        psi_G = two_torsion_part(E, psi).monic()

        if (0 != psi_G.degree()): # even degree case

            psi_quo = psi//psi_G

            if (0 != psi_quo.degree()):
                raise NotImplementedError("For basic Kohel's algorithm, if the kernel degree is even then the kernel must be contained in the two torsion.")

            (phi, omega, v, w, n, d) = self.__init_even_kernel_polynomial(E, psi_G)

        else: # odd degree case

            (phi, omega, v, w, n, d) = self.__init_odd_kernel_polynomial(E, psi)


        #
        # Set up the necessary instance variables
        #

        self.__kernel_polynomial = psi
        self.__inner_kernel_polynomial = psi

        self.__degree = Integer(d)  # degree of the isogeny

        # As a rational map, the isogeny maps (x,y) to (X,Y), where
        # X=phi(x)/psi(x)^2 and Y=omega(x,y)/psi(x)^3.  Both phi and
        # psi are univariate polynomials in x, while omega is a
        # bivariate polynomial in x, y.  The names are compatible so
        # that univariate polynomials automatically coerce into the
        # bivariate polynomial ring.

        self.__psi = psi
        self.__phi = phi
        self.__omega = omega

        self.__v = v
        self.__w = w

    def __init_even_kernel_polynomial(self, E, psi_G):
        r"""
        Return the isogeny parameters for the 2-part of an isogeny.

        INPUT:

        - ``E`` -- an elliptic curve

        - ``psi_G`` -- a univariate polynomial over the base field of
          ``E`` of degree 1 or 3 dividing its 2-division polynomial

        OUTPUT:

        (phi, omega, v, w, n, d) where:

        - ``phi`` is a univariate polynomial, the numerator of the
          `X`-coordinate of the isogeny;

        - ``omega`` is a bivariate polynomial, the numerator of the
          `Y`-coordinate of the isogeny;

        - ``v``, ``w`` are the Velu parameters of the isogeny;

        - ``n`` is the degree of ``psi``;

        - ``d`` is the degree of the isogeny.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [-1,0])
            sage: phi = EllipticCurveIsogeny(E, x); phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 4*x over Finite Field of size 7

            sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import two_torsion_part
            sage: psig = two_torsion_part(E,x)
            sage: phi._EllipticCurveIsogeny__init_even_kernel_polynomial(E,psig)
            (x^3 + 6*x, x^3*y + x*y, 6, 0, 1, 2)

            sage: F = GF(2^4, 'alpha'); R.<x> = F[]
            sage: E = EllipticCurve(F, [1,1,0,1,0])
            sage: phi = EllipticCurveIsogeny(E, x); phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + x over Finite Field in alpha of size 2^4 to Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + 1 over Finite Field in alpha of size 2^4

            sage: psig = two_torsion_part(E,x)
            sage: phi._EllipticCurveIsogeny__init_even_kernel_polynomial(E,psig)
            (x^3 + x, x^3*y + x^2 + x*y, 1, 0, 1, 2)

            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: R.<x> = GF(7)[]
            sage: f = x^3 + 6*x^2 + 1
            sage: phi = EllipticCurveIsogeny(E, f); phi
            Isogeny of degree 4 from Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 1 over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 2*x + 5 over Finite Field of size 7
            sage: psig = two_torsion_part(E,f)
            sage: psig = two_torsion_part(E,f)
            sage: phi._EllipticCurveIsogeny__init_even_kernel_polynomial(E,psig)
            (x^7 + 5*x^6 + 2*x^5 + 6*x^4 + 3*x^3 + 5*x^2 + 6*x + 3, x^9*y - 3*x^8*y + 2*x^7*y - 3*x^3*y + 2*x^2*y + x*y - y, 1, 6, 3, 4)
        """
        #check if the polynomial really divides the two_torsion_polynomial
        if  self.__check and E.division_polynomial(2, x=self.__poly_ring.gen()) % psi_G  != 0 :
            raise ValueError("The polynomial {} does not define a finite subgroup of {}.".format(psi_G,E))

        n = psi_G.degree() # 1 or 3
        d = n+1            # 2 or 4

        base_field = self.__base_field
        char = base_field.characteristic()

        a1,a2,a3,a4,a6 = E.ainvs()
        b2,b4,_,_ = E.b_invariants()
        x = self.__poly_ring.gen()
        y = self.__mpoly_ring.gen(1)

        if (1 == n):
            x0 = -psi_G.constant_coefficient()

            # determine y0
            if (2 == char):
                y0 = (x0**3 + a2*x0**2 + a4*x0 + a6).sqrt()
            else:
                y0 = -(a1*x0 + a3)/2

            (v,w) = compute_vw_kohel_even_deg1(x0,y0,a1,a2,a4)

            phi = (x*psi_G + v)*psi_G
            omega = (y*psi_G**2 - v*(a1*psi_G + (y - y0)))*psi_G

        elif (3 == n):
            s = psi_G.list()
            s1 = -s[n-1]
            s2 = s[n-2]
            s3 = -s[n-3]

            psi_G_pr = psi_G.derivative()
            psi_G_prpr = psi_G_pr.derivative()

            phi = (psi_G_pr**2) + (-2*psi_G_prpr + (4*x - s1))*psi_G
            phi_pr = phi.derivative(x)

            psi_2 = 2*y + a1*x + a3

            omega = (psi_2*(phi_pr*psi_G - phi*psi_G_pr) - (a1*phi + a3*psi_G)*psi_G)/2

            phi = phi*psi_G
            omega = omega*psi_G

            (v,w) = compute_vw_kohel_even_deg3(b2,b4,s1,s2,s3)

        else:
            raise ValueError("input polynomial must be of degree 1 or 3, not %d" % n)

        return (phi, omega, v, w, n, d)


    def __init_odd_kernel_polynomial(self, E, psi):
        r"""
        Return the isogeny parameters for a cyclic isogeny of odd degree.

        INPUT:

        - ``E`` -- an elliptic curve

        - ``psi`` -- a univariate polynomial over the base field of
          ``E``, assumed to be a kernel polynomial

        OUTPUT:

        (phi, omega, v, w, n, d) where:

        - ``phi`` is a univariate polynomial, the numerator of the
          `X`-coordinate of the isogeny;

        - ``omega`` is a bivariate polynomial, the numerator of the
          `Y`-coordinate of the isogeny;

        - ``v``, ``w`` are the Velu parameters of the isogeny;

        - ``n`` is the degree of ``psi``;

        - ``d`` is the degree of the isogeny.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3); phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 1 over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 4*x + 2 over Finite Field of size 7

            sage: R.<x> = GF(7)[]
            sage: phi._EllipticCurveIsogeny__init_odd_kernel_polynomial(E, x+6)
            (x^3 + 5*x^2 + 3*x + 2, x^3*y - 3*x^2*y + x*y, 2, 6, 1, 3)

            sage: F = GF(2^4, 'alpha'); R.<x> = F[]
            sage: alpha = F.gen()
            sage: E = EllipticCurve(F, [1,1,F.gen(),F.gen()^2+1,1])
            sage: f = x + alpha^2 + 1
            sage: phi = EllipticCurveIsogeny(E, f); phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 + x*y + alpha*y = x^3 + x^2 + (alpha^2+1)*x + 1 over Finite Field in alpha of size 2^4 to Elliptic Curve defined by y^2 + x*y + alpha*y = x^3 + x^2 + alpha*x + alpha^3 over Finite Field in alpha of size 2^4

            sage: R.<x> = F[]
            sage: f = x + alpha^2 + 1
            sage: phi._EllipticCurveIsogeny__init_odd_kernel_polynomial(E, f)
            (x^3 + (alpha^2 + 1)*x + alpha^3 + alpha^2 + alpha, x^3*y + (alpha^2 + 1)*x^2*y + (alpha^2 + alpha + 1)*x^2 + (alpha^2 + 1)*x*y + (alpha^2 + alpha)*x + alpha*y + alpha, alpha^2 + alpha + 1, alpha^3 + alpha^2 + alpha, 1, 3)

            sage: E = EllipticCurve(j=-262537412640768000)
            sage: f = (E.isogenies_prime_degree()[0]).kernel_polynomial()
            sage: f.degree()
            81
            sage: E.isogeny(kernel=f, check=False)
            Isogeny of degree 163 from Elliptic Curve defined by y^2 + y = x^3 - 2174420*x + 1234136692 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - 57772164980*x - 5344733777551611 over Rational Field
        """
        n = psi.degree()
        d = 2*n + 1

        # check if the polynomial really divides the torsion polynomial :
        if self.__check:
            from .isogeny_small_degree import is_kernel_polynomial
            if not is_kernel_polynomial(E, d, psi):
                raise ValueError("The polynomial {} does not define a finite subgroup of {}.".format(psi,E))

        b2, b4, b6, _ = E.b_invariants()

        psi_coeffs = psi.list()

        s1 = 0
        s2 = 0
        s3 = 0

        if 1 <= n:
            s1 = -psi_coeffs[n-1]

        if 2 <= n:
            s2 = psi_coeffs[n-2]

        if 3 <= n:
            s3 = -psi_coeffs[n-3]

        # initializing these allows us to calculate E2.
        (v,w) = compute_vw_kohel_odd(b2,b4,b6,s1,s2,s3,n)

        # initialize the polynomial temporary variables

        psi_pr = psi.derivative()
        psi_prpr = psi_pr.derivative()

        x = self.__poly_ring.gen()

        phi = (4*x**3 + b2*x**2 + 2*b4*x + b6)*(psi_pr**2 - psi_prpr*psi) - \
                (6*x**2 + b2*x + b4)*psi_pr*psi + (d*x - 2*s1)*psi**2

        phi_pr = phi.derivative(x)

        if (2 != self.__base_field.characteristic()):
            omega = self.__compute_omega_fast(E, psi, psi_pr, phi, phi_pr)
        else:
            omega = self.__compute_omega_general(E, psi, psi_pr, phi, phi_pr)

        return (phi, omega, v, w, n, d)


    #
    # This is the fast omega computation that works when characteristic is not 2
    #
    def __compute_omega_fast(self, E, psi, psi_pr, phi, phi_pr):
        r"""
        Return omega from phi, psi and their derivatives, used when
        the characteristic field is not 2.

        INPUT:

        - ``E`` -- an elliptic curve.

        - ``psi, psi_pr, phi, phi_pr`` -- univariate polynomials over
          the base field of ``E``, where ``psi`` is the kernel
          polynomial and ``phi`` the numerator of the `X`-coordinate
          of the isogeny, together with their derivatives.

        OUTPUT:

        - ``omega`` -- a bivariate polynomial giving the numerator of
          the `Y`-coordinate of the isogeny.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3); phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 1 over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 4*x + 2 over Finite Field of size 7

            sage: R.<x,y> = GF(7)[]
            sage: psi = phi._EllipticCurveIsogeny__psi
            sage: psi_pr = psi.derivative()
            sage: fi = phi._EllipticCurveIsogeny__phi
            sage: fi_pr = fi.derivative()
            sage: phi._EllipticCurveIsogeny__compute_omega_fast(E, psi, psi_pr, fi, fi_pr)
            x^3*y - 3*x^2*y + x*y

        """

        a1 = E.a1()
        a3 = E.a3()

        x, y = self.__mpoly_ring.gens()

        psi_2 = 2*y + a1*x + a3

        # note, the formula below is correct
        # the formula in Kohel's thesis has some typos
        # notably the first plus sign should be a minus
        # as it is here below.

        return phi_pr*psi*psi_2/2 - phi*psi_pr*psi_2 - (a1*phi + a3*psi**2)*psi/2

    def __compute_omega_general(self, E, psi, psi_pr, phi, phi_pr):
        r"""
        Return omega from phi, psi and their derivatives, in any
        characteristic.

        INPUT:

        - ``E`` -- an elliptic curve.

        - ``psi, psi_pr, phi, phi_pr`` -- univariate polynomials over
          the base field of ``E``, where ``psi`` is the kernel
          polynomial and ``phi`` the numerator of the `X`-coordinate
          of the isogeny, together with their derivatives.

        OUTPUT:

        - ``omega`` -- a bivariate polynomial giving the numerator of
          the `Y`-coordinate of the isogeny.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: F = GF(2^4, 'alpha'); R.<x> = F[]
            sage: alpha = F.gen()
            sage: E = EllipticCurve(F, [1,1,F.gen(),F.gen()^2+1,1])
            sage: f = x + alpha^2 + 1
            sage: phi = EllipticCurveIsogeny(E, f); phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 + x*y + alpha*y = x^3 + x^2 + (alpha^2+1)*x + 1 over Finite Field in alpha of size 2^4 to Elliptic Curve defined by y^2 + x*y + alpha*y = x^3 + x^2 + alpha*x + alpha^3 over Finite Field in alpha of size 2^4

            sage: R.<x,y> = F[]
            sage: psi = phi._EllipticCurveIsogeny__psi
            sage: psi_pr = psi.derivative()
            sage: fi = phi._EllipticCurveIsogeny__phi
            sage: fi_pr = fi.derivative()
            sage: phi._EllipticCurveIsogeny__compute_omega_general(E, psi, psi_pr, fi, fi_pr)
            x^3*y + (alpha^2 + 1)*x^2*y + (alpha^2 + alpha + 1)*x^2 + (alpha^2 + 1)*x*y + (alpha^2 + alpha)*x + alpha*y + alpha

        A bug fixed in :trac:`7907`::

            sage: F = GF(128,'a')
            sage: a = F.gen()
            sage: E = EllipticCurve([1,0,0,0,(a**6+a**4+a**2+a)])
            sage: x = polygen(F)
            sage: ker =  (x^6 + (a^6 + a^5 + a^4 + a^3 + a^2 + a)*x^5 + (a^6 + a^5 + a^2 + 1)*x^4 + (a^6 + a^5 + a^4 + a^3 + a^2 + 1)*x^3 + (a^6 + a^3 + a)*x^2 + (a^4 + a^3 + 1)*x + a^5 + a^4 + a)
            sage: E.isogeny(ker)
            Isogeny of degree 13 from Elliptic Curve defined by y^2 + x*y = x^3 + (a^6+a^4+a^2+a) over Finite Field in a of size 2^7 to Elliptic Curve defined by y^2 + x*y = x^3 + (a^6+a^5+a^4+a^3+a^2+a)*x + (a^5+a^3) over Finite Field in a of size 2^7


        """
        a1, a2, a3, a4, a6 = E.ainvs()
        b2, b4, _, _ = E.b_invariants()

        n = psi.degree()
        d = 2 * n + 1

        x, y = self.__mpoly_ring.gens()

        psi_2 = 2 * y + a1 * x + a3

        psi_coeffs = psi.list()

        if (0 < n):
            s1 = -psi_coeffs[n - 1]
        else:
            s1 = 0

        psi_prpr = 0
        cur_x_pow = 1

        # Note: we now get the "derivatives" of psi
        # these are not actually the derivatives
        # furthermore, the formulas in Kohel's
        # thesis are wrong, the correct formulas
        # are coded below

        from sage.arith.all import binomial

        for j in range(n - 1):
            psi_prpr += binomial(j+2, 2) * psi_coeffs[(j+2)] * cur_x_pow
            cur_x_pow = x * cur_x_pow

        psi_prprpr = 0
        cur_x_pow = 1

        for j in range(n - 2):
            psi_prprpr += (3 * binomial(j+3,3)) * psi_coeffs[(j+3)] * cur_x_pow
            cur_x_pow = x * cur_x_pow

        omega = phi_pr*psi*y - phi*psi_pr*psi_2 + \
                ((a1*x + a3)*(psi_2**2)*(psi_prpr*psi_pr-psi_prprpr*psi) + \
                (a1*psi_2**2 - 3*(a1*x + a3)*(6*x**2 + b2*x + b4))*psi_prpr*psi + \
                (a1*x**3 + 3*a3*x**2 + (2*a2*a3 - a1*a4)*x + (a3*a4 - 2*a1*a6))*psi_pr**2 + \
                (-(3*a1*x**2 + 6*a3*x + (-a1*a4 + 2*a2*a3)) + \
                (a1*x + a3)*(d*x - 2*s1) )*psi_pr*psi + (a1*s1 + a3*n)*psi**2)*psi

        return omega

    def __compute_via_kohel_numeric(self, xP, yP):
        r"""
        Private function that computes the image of a point under this
        isogeny, using Kohel's formulas.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3)
            sage: P = E((0,1)); phi(P)
            (2 : 0 : 1)
            sage: P = E((1,1)); phi(P)
            (0 : 1 : 0)
            sage: phi._EllipticCurveIsogeny__compute_via_kohel_numeric(0, 1)
            (2, 0)
            sage: phi._EllipticCurveIsogeny__compute_via_kohel_numeric(1, 1)
            (0 : 1 : 0)

        """
        # first check if this point is in the kernel:
        if 0 == self.__inner_kernel_polynomial(x=xP):
            return self.__intermediate_codomain(0)

        return self.__compute_via_kohel(xP, yP)

    def __compute_via_kohel(self, xP, yP):
        r"""
        Private function that applies Kohel's formulas.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3)
            sage: P = E((0,1)); phi(P)
            (2 : 0 : 1)
            sage: phi.rational_maps()
            ((x^3 - 2*x^2 + 3*x + 2)/(x^2 - 2*x + 1), (x^3*y - 3*x^2*y + x*y)/(x^3 - 3*x^2 + 3*x - 1))
            sage: phi._EllipticCurveIsogeny__compute_via_kohel(0,1)
            (2, 0)
            sage: R.<x,y> = GF(7)[]
            sage: phi._EllipticCurveIsogeny__compute_via_kohel(x,y)
            ((x^3 - 2*x^2 + 3*x + 2)/(x^2 - 2*x + 1), (x^3*y - 3*x^2*y + x*y)/(x^3 - 3*x^2 + 3*x - 1))

        """
        a = self.__phi(xP)
        b = self.__omega(xP, yP)
        c = self.__psi(xP)
        return (a/c**2, b/c**3)

    def __initialize_rational_maps_via_kohel(self):
        r"""
        Private function that computes and initializes the rational
        maps of this isogeny.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3)
            sage: phi.rational_maps()
            ((x^3 - 2*x^2 + 3*x + 2)/(x^2 - 2*x + 1), (x^3*y - 3*x^2*y + x*y)/(x^3 - 3*x^2 + 3*x - 1))
            sage: phi._EllipticCurveIsogeny__initialize_rational_maps_via_kohel()
            ((x^3 + 5*x^2 + 3*x + 2)/(x^2 + 5*x + 1), (x^3*y - 3*x^2*y + x*y)/(x^3 - 3*x^2 + 3*x - 1))


        """
        x = self.__poly_ring.gen()
        y = self.__xyfield.gen(1)
        return self.__compute_via_kohel(x,y)

    #
    # Kohel's formula computing the codomain curve
    #
    def __compute_E2_via_kohel(self):
        r"""
        Private function that computes and initializes the codomain of
        the isogeny (via Kohel's.)

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3)
            sage: phi.codomain()
            Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 4*x + 2 over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__compute_E2_via_kohel()
            Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 4*x + 2 over Finite Field of size 7

        """

        v = self.__v
        w = self.__w

        return compute_codomain_formula(self.__E1, v,w)

    #
    # public isogeny methods
    #

    def degree(self):
        r"""
        Return the degree of this isogeny.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.degree()
            2
            sage: phi = EllipticCurveIsogeny(E, [0,1,0,1])
            sage: phi.degree()
            4

            sage: E = EllipticCurve(GF(31), [1,0,0,1,2])
            sage: phi = EllipticCurveIsogeny(E, [17, 1])
            sage: phi.degree()
            3

        """
        return self.__degree

    def rational_maps(self):
        r"""
        Return the pair of rational maps defining this isogeny.

        .. NOTE::

           Both components are returned as elements of the function
           field `F(x,y)` in two variables over the base field `F`,
           though the first only involves `x`.  To obtain the
           `x`-coordinate function as a rational function in `F(x)`,
           use :meth:`x_rational_map`.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,2,0,1,-1])
            sage: phi = EllipticCurveIsogeny(E, [1])
            sage: phi.rational_maps()
            (x, y)

            sage: E = EllipticCurve(GF(17), [0,0,0,3,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.rational_maps()
            ((x^2 + 3)/x, (x^2*y - 3*y)/x^2)
        """
        if (not self.__rational_maps_initialized):
            self.__initialize_rational_maps()
        return (self.__xyfield(self.__X_coord_rational_map),
                self.__Y_coord_rational_map)

    def x_rational_map(self):
        r"""
        Return the rational map giving the `x`-coordinate of this isogeny.

        .. NOTE::

           This function returns the `x`-coordinate component of the
           isogeny as a rational function in `F(x)`, where `F` is the
           base field.  To obtain both coordinate functions as
           elements of the function field `F(x,y)` in two variables,
           use :meth:`rational_maps`.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,2,0,1,-1])
            sage: phi = EllipticCurveIsogeny(E, [1])
            sage: phi.x_rational_map()
            x

            sage: E = EllipticCurve(GF(17), [0,0,0,3,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.x_rational_map()
            (x^2 + 3)/x
        """
        if (not self.__rational_maps_initialized):
            self.__initialize_rational_maps()
        return self.__X_coord_rational_map

    def kernel_polynomial(self):
        r"""
        Return the kernel polynomial of this isogeny.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,0,0,2,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.kernel_polynomial()
            x

            sage: E = EllipticCurve('11a1')
            sage: phi = EllipticCurveIsogeny(E, E.torsion_points())
            sage: phi.kernel_polynomial()
            x^2 - 21*x + 80

            sage: E = EllipticCurve(GF(17), [1,-1,1,-1,1])
            sage: phi = EllipticCurveIsogeny(E, [1])
            sage: phi.kernel_polynomial()
            1

            sage: E = EllipticCurve(GF(31), [0,0,0,3,0])
            sage: phi = EllipticCurveIsogeny(E, [0,3,0,1])
            sage: phi.kernel_polynomial()
            x^3 + 3*x
        """
        if self.__kernel_polynomial is None:
            self.__init_kernel_polynomial()

        return self.__kernel_polynomial


    def set_pre_isomorphism(self, preWI):
        r"""
        Modify this isogeny by precomposing with a Weierstrass isomorphism.

        .. WARNING::

            Isogenies will be immutable in a future release of Sage.
            This method is deprecated in favor of using the ``*`` operator
            to compose elliptic-curve morphisms.

        EXAMPLES::

            sage: E = EllipticCurve(GF(31), [1,1,0,1,-1])
            sage: R.<x> = GF(31)[]
            sage: f = x^3 + 9*x^2 + x + 30
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: Epr = E.short_weierstrass_model()
            sage: isom = Epr.isomorphism_to(E)
            sage: phi.set_pre_isomorphism(isom)
            ...
            sage: phi.rational_maps()
            ((-6*x^4 - 3*x^3 + 12*x^2 + 10*x - 1)/(x^3 + x - 12), (3*x^7 + x^6*y - 14*x^6 - 3*x^5 + 5*x^4*y + 7*x^4 + 8*x^3*y - 8*x^3 - 5*x^2*y + 5*x^2 - 14*x*y + 14*x - 6*y - 6)/(x^6 + 2*x^4 + 7*x^3 + x^2 + 7*x - 11))
            sage: phi(Epr((0,22)))
            (13 : 21 : 1)
            sage: phi(Epr((3,7)))
            (14 : 17 : 1)

            sage: E = EllipticCurve(GF(29), [0,0,0,1,0])
            sage: R.<x> = GF(29)[]
            sage: f = x^2 + 5
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi
            Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 29 to Elliptic Curve defined by y^2 = x^3 + 20*x over Finite Field of size 29
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: inv_isom = WeierstrassIsomorphism(E, (1,-2,5,10))
            sage: Epr = inv_isom.codomain()
            sage: isom = Epr.isomorphism_to(E)
            sage: phi.set_pre_isomorphism(isom)
            ...
            sage: phi
            Isogeny of degree 5 from Elliptic Curve defined by y^2 + 10*x*y + 20*y = x^3 + 27*x^2 + 6 over Finite Field of size 29 to Elliptic Curve defined by y^2 = x^3 + 20*x over Finite Field of size 29
            sage: phi(Epr((12,1)))
            (26 : 0 : 1)
            sage: phi(Epr((2,9)))
            (0 : 0 : 1)
            sage: phi(Epr((21,12)))
            (3 : 0 : 1)
            sage: phi.rational_maps()[0]
            (x^5 - 10*x^4 - 6*x^3 - 7*x^2 - x + 3)/(x^4 - 8*x^3 + 5*x^2 - 14*x - 6)

            sage: E = EllipticCurve('11a1')
            sage: R.<x> = QQ[]
            sage: f = x^2 - 21*x + 80
            sage: phi = EllipticCurveIsogeny(E, f); phi
            Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: Epr = E.short_weierstrass_model()
            sage: isom = Epr.isomorphism_to(E)
            sage: phi.set_pre_isomorphism(isom)
            ...
            sage: phi
            Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 - 13392*x - 1080432 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
            sage: phi(Epr((168,1188)))
            (0 : 1 : 0)
        """
        from sage.misc.superseded import deprecation
        deprecation(32388, 'Elliptic-curve isogenies will be immutable in a future release of Sage.'
                          ' Use phi*psi instead of phi.set_pre_isomorphism(psi) to obtain the composite isogeny.')
        return self._set_pre_isomorphism(preWI)

    def _set_pre_isomorphism(self, preWI):
        """
        Modify this isogeny by pre-composing with a
        :class:`sage.schemes.elliptic_curves.weierstrass_morphism.WeierstrassIsomorphism`.

        For internal use only.

        TESTS:

        These tests were copied from :meth:`set_pre_isomorphism`::

            sage: E = EllipticCurve(GF(31), [1,1,0,1,-1])
            sage: R.<x> = GF(31)[]
            sage: f = x^3 + 9*x^2 + x + 30
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: Epr = E.short_weierstrass_model()
            sage: isom = Epr.isomorphism_to(E)
            sage: phi._set_pre_isomorphism(isom)
            sage: phi.rational_maps()
            ((-6*x^4 - 3*x^3 + 12*x^2 + 10*x - 1)/(x^3 + x - 12), (3*x^7 + x^6*y - 14*x^6 - 3*x^5 + 5*x^4*y + 7*x^4 + 8*x^3*y - 8*x^3 - 5*x^2*y + 5*x^2 - 14*x*y + 14*x - 6*y - 6)/(x^6 + 2*x^4 + 7*x^3 + x^2 + 7*x - 11))
            sage: phi(Epr((0,22)))
            (13 : 21 : 1)
            sage: phi(Epr((3,7)))
            (14 : 17 : 1)

            sage: E = EllipticCurve(GF(29), [0,0,0,1,0])
            sage: R.<x> = GF(29)[]
            sage: f = x^2 + 5
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi
            Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 29 to Elliptic Curve defined by y^2 = x^3 + 20*x over Finite Field of size 29
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: inv_isom = WeierstrassIsomorphism(E, (1,-2,5,10))
            sage: Epr = inv_isom.codomain()
            sage: isom = Epr.isomorphism_to(E)
            sage: phi._set_pre_isomorphism(isom)
            sage: phi
            Isogeny of degree 5 from Elliptic Curve defined by y^2 + 10*x*y + 20*y = x^3 + 27*x^2 + 6 over Finite Field of size 29 to Elliptic Curve defined by y^2 = x^3 + 20*x over Finite Field of size 29
            sage: phi(Epr((12,1)))
            (26 : 0 : 1)
            sage: phi(Epr((2,9)))
            (0 : 0 : 1)
            sage: phi(Epr((21,12)))
            (3 : 0 : 1)
            sage: phi.rational_maps()[0]
            (x^5 - 10*x^4 - 6*x^3 - 7*x^2 - x + 3)/(x^4 - 8*x^3 + 5*x^2 - 14*x - 6)

            sage: E = EllipticCurve('11a1')
            sage: R.<x> = QQ[]
            sage: f = x^2 - 21*x + 80
            sage: phi = EllipticCurveIsogeny(E, f); phi
            Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: Epr = E.short_weierstrass_model()
            sage: isom = Epr.isomorphism_to(E)
            sage: phi._set_pre_isomorphism(isom)
            sage: phi
            Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 - 13392*x - 1080432 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
            sage: phi(Epr((168,1188)))
            (0 : 1 : 0)
        """
        WIdom = preWI.domain()
        WIcod = preWI.codomain()

        if not isinstance(preWI, WeierstrassIsomorphism):
            raise ValueError("Invalid parameter: isomorphism must be of type Weierstrass isomorphism.")

        if (self.__E1 != WIcod):
            raise ValueError("Invalid parameter: isomorphism must have codomain curve equal to this isogenies' domain.")

        if (self.__pre_isomorphism is None):
            isom = preWI
            domain = WIdom
        else:
            isom = self.__pre_isomorphism*preWI
            domain = WIdom

        self.__clear_cached_values()

        self.__set_pre_isomorphism(domain, isom)

        return


    def set_post_isomorphism(self, postWI):
        r"""
        Modify this isogeny by postcomposing with a Weierstrass isomorphism.

        .. WARNING::

            Isogenies will be immutable in a future release of Sage.
            This method is deprecated in favor of using the ``*`` operator
            to compose elliptic-curve morphisms.

        EXAMPLES::

            sage: E = EllipticCurve(j=GF(31)(0))
            sage: R.<x> = GF(31)[]
            sage: phi = EllipticCurveIsogeny(E, x+18)
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: phi.set_post_isomorphism(WeierstrassIsomorphism(phi.codomain(), (6,8,10,12)))
            ...
            sage: phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 31 to Elliptic Curve defined by y^2 + 24*x*y + 7*y = x^3 + 22*x^2 + 16*x + 20 over Finite Field of size 31

            sage: E = EllipticCurve(j=GF(47)(0))
            sage: f = E.torsion_polynomial(3)/3
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: E2 = phi.codomain()
            sage: post_isom = E2.isomorphism_to(E)
            sage: phi.set_post_isomorphism(post_isom)
            ...
            sage: phi.rational_maps() == E.multiplication_by_m(3)
            False
            sage: phi.switch_sign()
            ...
            sage: phi.rational_maps() == E.multiplication_by_m(3)
            True

        Example over a number field::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^2 + 2)
            sage: E = EllipticCurve(j=K(1728))
            sage: ker_list = E.torsion_points()
            sage: phi = EllipticCurveIsogeny(E, ker_list)
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: post_isom = WeierstrassIsomorphism(phi.codomain(), (a,2,3,5))
            sage: phi
            Isogeny of degree 4 from Elliptic Curve defined by y^2 = x^3 + x over Number Field in a with defining polynomial x^2 + 2 to Elliptic Curve defined by y^2 = x^3 + (-44)*x + 112 over Number Field in a with defining polynomial x^2 + 2

        """
        from sage.misc.superseded import deprecation
        deprecation(32388, 'Elliptic-curve isogenies will be immutable in a future release of Sage.'
                          ' Use psi*phi instead of phi.set_post_isomorphism(psi) to obtain the composite isogeny.')
        return self._set_post_isomorphism(postWI)

    def _set_post_isomorphism(self, postWI):
        """
        Modify this isogeny by post-composing with a
        :class:`sage.schemes.elliptic_curves.weierstrass_morphism.WeierstrassIsomorphism`.

        For internal use only.

        TESTS:

        These tests were copied from :meth:`set_post_isomorphism`::

            sage: E = EllipticCurve(j=GF(31)(0))
            sage: R.<x> = GF(31)[]
            sage: phi = EllipticCurveIsogeny(E, x+18)
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: phi._set_post_isomorphism(WeierstrassIsomorphism(phi.codomain(), (6,8,10,12)))
            sage: phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 31 to Elliptic Curve defined by y^2 + 24*x*y + 7*y = x^3 + 22*x^2 + 16*x + 20 over Finite Field of size 31

            sage: E = EllipticCurve(j=GF(47)(0))
            sage: f = E.torsion_polynomial(3)/3
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: E2 = phi.codomain()
            sage: post_isom = E2.isomorphism_to(E)
            sage: phi._set_post_isomorphism(post_isom)
            sage: phi.rational_maps() == E.multiplication_by_m(3)
            False
            sage: phi = -phi
            sage: phi.rational_maps() == E.multiplication_by_m(3)
            True

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^2 + 2)
            sage: E = EllipticCurve(j=K(1728))
            sage: ker_list = E.torsion_points()
            sage: phi = EllipticCurveIsogeny(E, ker_list)
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: post_isom = WeierstrassIsomorphism(phi.codomain(), (a,2,3,5))
            sage: phi
            Isogeny of degree 4 from Elliptic Curve defined by y^2 = x^3 + x over Number Field in a with defining polynomial x^2 + 2 to Elliptic Curve defined by y^2 = x^3 + (-44)*x + 112 over Number Field in a with defining polynomial x^2 + 2
        """
        WIdom = postWI.domain()
        WIcod = postWI.codomain()

        if not isinstance(postWI, WeierstrassIsomorphism):
            raise ValueError("Invalid parameter: isomorphism must be of type Weierstrass isomorphism.")

        if (self.__E2 != WIdom):
            raise ValueError("Invalid parameter: isomorphism must have domain curve equal to this isogenies' codomain.")

        if (self.__post_isomorphism is None):
            isom = postWI
            codomain = WIcod
        else:
            isom = postWI*self.__post_isomorphism
            codomain = WIcod

        self.__clear_cached_values()

        self.__set_post_isomorphism(codomain, isom)

        return


    def get_pre_isomorphism(self):
        r"""
        Return the pre-isomorphism of this isogeny, or ``None``.

        EXAMPLES::

            sage: E = EllipticCurve(GF(31), [1,1,0,1,-1])
            sage: R.<x> = GF(31)[]
            sage: f = x^3 + 9*x^2 + x + 30
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi.get_post_isomorphism()
            sage: Epr = E.short_weierstrass_model()
            sage: isom = Epr.isomorphism_to(E)
            sage: phi.set_pre_isomorphism(isom)
            ...
            sage: isom == phi.get_pre_isomorphism()
            True

            sage: E = EllipticCurve(GF(83), [1,0,1,1,0])
            sage: R.<x> = GF(83)[]; f = x+24
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: E2 = phi.codomain()
            sage: phi2 = EllipticCurveIsogeny(E, None, E2, 2)
            sage: phi2.get_pre_isomorphism()
            Elliptic-curve morphism:
              From: Elliptic Curve defined by y^2 + x*y + y = x^3 + x over Finite Field of size 83
              To:   Elliptic Curve defined by y^2 = x^3 + 62*x + 74 over Finite Field of size 83
              Via:  (u,r,s,t) = (1, 76, 41, 3)
        """
        return self.__pre_isomorphism

    def get_post_isomorphism(self):
        r"""
        Return the post-isomorphism of this isogeny, or ``None``.

        EXAMPLES::

            sage: E = EllipticCurve(j=GF(31)(0))
            sage: R.<x> = GF(31)[]
            sage: phi = EllipticCurveIsogeny(E, x+18)
            sage: phi.get_post_isomorphism()
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (6,8,10,12))
            sage: phi.set_post_isomorphism(isom)
            ...
            sage: isom == phi.get_post_isomorphism()
            True

            sage: E = EllipticCurve(GF(83), [1,0,1,1,0])
            sage: R.<x> = GF(83)[]; f = x+24
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: E2 = phi.codomain()
            sage: phi2 = EllipticCurveIsogeny(E, None, E2, 2)
            sage: phi2.get_post_isomorphism()
            Elliptic-curve morphism:
            From: Elliptic Curve defined by y^2 = x^3 + 65*x + 69 over Finite Field of size 83
            To:   Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + 16 over Finite Field of size 83
            Via:  (u,r,s,t) = (1, 7, 42, 42)
        """
        return self.__post_isomorphism


    def switch_sign(self):
        r"""
        Compose this isogeny with `[-1]` (negation).

        .. WARNING::

            Isogenies will be immutable in a future version of Sage.
            This method is deprecated in favor of using the unary ``-``
            operator to negate elliptic-curve morphisms.

        EXAMPLES::

            sage: E = EllipticCurve(GF(23), [0,0,0,1,0])
            sage: f = E.torsion_polynomial(3)/3
            sage: phi = EllipticCurveIsogeny(E, f, E)
            sage: phi.rational_maps() == E.multiplication_by_m(3)
            False
            sage: phi.switch_sign()
            ...
            sage: phi.rational_maps() == E.multiplication_by_m(3)
            True

            sage: E = EllipticCurve(GF(17), [-2, 3, -5, 7, -11])
            sage: R.<x> = GF(17)[]
            sage: f = x+6
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 + 15*x*y + 12*y = x^3 + 3*x^2 + 7*x + 6 over Finite Field of size 17 to Elliptic Curve defined by y^2 + 15*x*y + 12*y = x^3 + 3*x^2 + 4*x + 8 over Finite Field of size 17
            sage: phi.rational_maps()
            ((x^2 + 6*x + 4)/(x + 6), (x^2*y - 5*x*y + 8*x - 2*y)/(x^2 - 5*x + 2))
            sage: phi.switch_sign()
            ...
            sage: phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 + 15*x*y + 12*y = x^3 + 3*x^2 + 7*x + 6 over Finite Field of size 17 to Elliptic Curve defined by y^2 + 15*x*y + 12*y = x^3 + 3*x^2 + 4*x + 8 over Finite Field of size 17
            sage: phi.rational_maps()
            ((x^2 + 6*x + 4)/(x + 6),
             (2*x^3 - x^2*y - 5*x^2 + 5*x*y - 4*x + 2*y + 7)/(x^2 - 5*x + 2))

            sage: E = EllipticCurve('11a1')
            sage: R.<x> = QQ[]
            sage: f = x^2 - 21*x + 80
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: (xmap1, ymap1) = phi.rational_maps()
            sage: phi.switch_sign()
            ...
            sage: (xmap2, ymap2) = phi.rational_maps()
            sage: xmap1 == xmap2
            True
            sage: ymap1 == -ymap2 - E.a1()*xmap2 - E.a3()
            True

            sage: K.<a> = NumberField(x^2 + 1)
            sage: E = EllipticCurve(K, [0,0,0,1,0])
            sage: R.<x> = K[]
            sage: phi = EllipticCurveIsogeny(E, x-a)
            sage: phi.rational_maps()
            ((x^2 + (-a)*x - 2)/(x + (-a)), (x^2*y + (-2*a)*x*y + y)/(x^2 + (-2*a)*x - 1))
            sage: phi.switch_sign()
            ...
            sage: phi.rational_maps()
            ((x^2 + (-a)*x - 2)/(x + (-a)), (-x^2*y + (2*a)*x*y - y)/(x^2 + (-2*a)*x - 1))

        """
        from sage.misc.superseded import deprecation
        deprecation(32388, 'Elliptic-curve isogenies will be immutable in a future release of Sage.'
                          ' Use -phi instead of phi.switch_sign() to obtain the negated isogeny.')
        self._set_post_isomorphism(WeierstrassIsomorphism(self.__E2, (-1,0,-self.__E2.a1(),-self.__E2.a3())))

    def dual(self):
        r"""
        Return the isogeny dual to this isogeny.

        .. NOTE::

           If `\varphi\colon E \to E_2` is the given isogeny and `n`
           is its degree, then the dual is by definition the unique
           isogeny `\hat\varphi\colon E_2\to E` such that the
           compositions `\hat\varphi\circ\varphi` and
           `\varphi\circ\hat\varphi` are the multiplication-by-`n`
           maps on `E` and `E_2`, respectively.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: R.<x> = QQ[]
            sage: f = x^2 - 21*x + 80
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi_hat = phi.dual()
            sage: phi_hat.domain() == phi.codomain()
            True
            sage: phi_hat.codomain() == phi.domain()
            True
            sage: (X, Y) = phi.rational_maps()
            sage: (Xhat, Yhat) = phi_hat.rational_maps()
            sage: Xm = Xhat.subs(x=X, y=Y)
            sage: Ym = Yhat.subs(x=X, y=Y)
            sage: (Xm, Ym) == E.multiplication_by_m(5)
            True

            sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
            sage: R.<x> = GF(37)[]
            sage: f = x^3 + x^2 + 28*x + 33
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi_hat = phi.dual()
            sage: phi_hat.codomain() == phi.domain()
            True
            sage: phi_hat.domain() == phi.codomain()
            True
            sage: (X, Y) = phi.rational_maps()
            sage: (Xhat, Yhat) = phi_hat.rational_maps()
            sage: Xm = Xhat.subs(x=X, y=Y)
            sage: Ym = Yhat.subs(x=X, y=Y)
            sage: (Xm, Ym) == E.multiplication_by_m(7)
            True

            sage: E = EllipticCurve(GF(31), [0,0,0,1,8])
            sage: R.<x> = GF(31)[]
            sage: f = x^2 + 17*x + 29
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi_hat = phi.dual()
            sage: phi_hat.codomain() == phi.domain()
            True
            sage: phi_hat.domain() == phi.codomain()
            True
            sage: (X, Y) = phi.rational_maps()
            sage: (Xhat, Yhat) = phi_hat.rational_maps()
            sage: Xm = Xhat.subs(x=X, y=Y)
            sage: Ym = Yhat.subs(x=X, y=Y)
            sage: (Xm, Ym) == E.multiplication_by_m(5)
            True

        Test for :trac:`23928`::

            sage: E = EllipticCurve(j=GF(431**2)(4))
            sage: phi = E.isogeny(E.lift_x(0))
            sage: phi.dual()
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 427*x over Finite Field in z2 of size 431^2 to Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 431^2

        Test (for :trac:`7096`)::

            sage: E = EllipticCurve('11a1')
            sage: phi = E.isogeny(E(5,5))
            sage: phi.dual().dual() == phi
            True

            sage: k = GF(103)
            sage: E = EllipticCurve(k,[11,11])
            sage: phi = E.isogeny(E(4,4))
            sage: phi
            Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + 11*x + 11 over Finite Field of size 103 to Elliptic Curve defined by y^2 = x^3 + 25*x + 80 over Finite Field of size 103
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: phi = WeierstrassIsomorphism(phi.codomain(),(5,0,1,2)) * phi
            sage: phi.dual().dual() == phi
            True

            sage: E = EllipticCurve(GF(103),[1,0,0,1,-1])
            sage: phi = E.isogeny(E(60,85))
            sage: phi.dual()
            Isogeny of degree 7 from Elliptic Curve defined by y^2 + x*y = x^3 + 84*x + 34 over Finite Field of size 103 to Elliptic Curve defined by y^2 + x*y = x^3 + x + 102 over Finite Field of size 103

        Check that :trac:`17293` is fixed::

            sage: k.<s> = QuadraticField(2)
            sage: E = EllipticCurve(k, [-3*s*(4 + 5*s), 2*s*(2 + 14*s + 11*s^2)])
            sage: phi = E.isogenies_prime_degree(3)[0]
            sage: (-phi).dual() == -(phi.dual())
            True
            sage: phi._EllipticCurveIsogeny__clear_cached_values()  # forget the dual
            sage: -(phi.dual()) == (-phi).dual()
            True

        """
        if (self.__base_field.characteristic() in [2,3]):
            raise NotImplementedError("Computation of dual isogenies not yet implemented in characteristics 2 and 3")

        if (self.__dual is not None):
            return self.__dual

        # trac 7096
        (E1, E2pr, pre_isom, post_isom) = compute_intermediate_curves(self.codomain(), self.domain())

        F = self.__base_field
        d = self.__degree

        # trac 7096
        if F(d) == 0:
            raise NotImplementedError("The dual isogeny is not separable: only separable isogenies are currently implemented")

        # trac 7096
        # this should take care of the case when the isogeny is not normalized.
        u = self.formal()[1]
        isom = WeierstrassIsomorphism(E2pr, (u/F(d), 0, 0, 0))

        E2 = isom.codomain()

        pre_isom = self.__E2.isomorphism_to(E1)
        post_isom = E2.isomorphism_to(self.__E1)

        phi_hat = EllipticCurveIsogeny(E1, None, E2, d)

        phi_hat._set_pre_isomorphism(pre_isom)
        phi_hat._set_post_isomorphism(post_isom)
        phi_hat.__perform_inheritance_housekeeping()

        assert phi_hat.codomain() == self.domain()

        # trac 7096 : this adjusts a posteriori the automorphism on
        # the codomain of the dual isogeny.  we used _a_ Weierstrass
        # isomorphism to get to the original curve, but we may have to
        # change it by an automorphism.  We impose the condition that
        # the composition has the degree as a leading coefficient in
        # the formal expansion.

        phi_sc = self.formal()[1]
        phihat_sc = phi_hat.formal()[1]

        sc = phi_sc * phihat_sc/F(d)

        if sc == 0:
            raise RuntimeError("Bug in computing dual isogeny: sc = 0")

        if sc != 1:
            auts = self.__E1.automorphisms()
            aut = [a for a in auts if a.u == sc]
            if len(aut) != 1:
                raise ValueError("There is a bug in dual().")
            phi_hat._set_post_isomorphism(aut[0])

        self.__dual = phi_hat

        return phi_hat


    @staticmethod
    def _composition_impl(left, right):
        r"""
        Return the composition of an ``EllipticCurveIsogeny``
        with another elliptic-curve morphism.

        Called by :meth:`EllipticCurveHom._composition_`.

        EXAMPLES::

            sage: E = EllipticCurve(GF(127), [5,2])
            sage: phi = E.isogeny(E.lift_x(47)); E2 = phi.codomain()
            sage: iso1 = E.change_weierstrass_model(1,1,1,1).isomorphism_to(E)
            sage: iso2 = E2.isomorphism_to(E2.change_weierstrass_model(39,0,0,0))
            sage: phi * iso1            # indirect doctest
            Isogeny of degree 11 from Elliptic Curve defined by y^2 + 2*x*y + 2*y = x^3 + 2*x^2 + 6*x + 7 over Finite Field of size 127 to Elliptic Curve defined by y^2 = x^3 + 37*x + 85 over Finite Field of size 127
            sage: iso2 * phi            # indirect doctest
            Isogeny of degree 11 from Elliptic Curve defined by y^2 = x^3 + 5*x + 2 over Finite Field of size 127 to Elliptic Curve defined by y^2 = x^3 + 117*x + 58 over Finite Field of size 127
            sage: iso2 * phi * iso1     # indirect doctest
            Isogeny of degree 11 from Elliptic Curve defined by y^2 + 2*x*y + 2*y = x^3 + 2*x^2 + 6*x + 7 over Finite Field of size 127 to Elliptic Curve defined by y^2 = x^3 + 117*x + 58 over Finite Field of size 127

        TESTS:

        We should return ``NotImplemented`` when passed a combination of
        elliptic-curve morphism types that we don't handle here::

            sage: phi._composition_impl(iso1, iso1**-1)
            NotImplemented
        """
        if isinstance(left, WeierstrassIsomorphism) and isinstance(right, EllipticCurveIsogeny):
            result = copy(right)
            result._set_post_isomorphism(left)
            return result

        if isinstance(left, EllipticCurveIsogeny) and isinstance(right, WeierstrassIsomorphism):
            assert isinstance(left, EllipticCurveIsogeny)
            result = copy(left)
            result._set_pre_isomorphism(right)
            return result

        return NotImplemented

    def n(self):
        r"""
        Numerical Approximation inherited from Map (through morphism),
        nonsensical for isogenies.

        EXAMPLES::

            sage: E = EllipticCurve(j=GF(7)(0))
            sage: phi = EllipticCurveIsogeny(E, [ E((0,1)), E((0,-1))])
            sage: phi.n()
            Traceback (most recent call last):
            ...
            NotImplementedError: Numerical approximations do not make sense for Elliptic Curve Isogenies
        """
        raise NotImplementedError("Numerical approximations do not make sense for Elliptic Curve Isogenies")

def compute_isogeny_starks(E1, E2, ell):
    r"""
    Return the kernel polynomials of an isogeny of degree ``ell``
    between ``E1`` and ``E2``.

    INPUT:

    - ``E1``  - an elliptic curve in short Weierstrass form.
    - ``E2``  - an elliptic curve in short Weierstrass form.
    - ``ell`` - the degree of the isogeny from E1 to E2.

    OUTPUT:

    polynomial over the field of definition of ``E1``, ``E2``, that is
    the kernel polynomial of the isogeny from ``E1`` to ``E2``.

    .. NOTE::

       There must be a degree ``ell``, separable, normalized cyclic
       isogeny from ``E1`` to ``E2``, or an error will be raised.

    ALGORITHM:

    This function uses Starks Algorithm as presented in section 6.2 of
    [BMSS2006]_.

    .. NOTE::

       As published in [BMSS2006]_, the algorithm is incorrect, and a
       correct version (with slightly different notation) can be found
       in [Mo2009]_.  The algorithm originates in [Sta1973]_.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_starks, compute_sequence_of_maps

        sage: E = EllipticCurve(GF(97), [1,0,1,1,0])
        sage: R.<x> = GF(97)[]; f = x^5 + 27*x^4 + 61*x^3 + 58*x^2 + 28*x + 21
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: (isom1, isom2, E1pr, E2pr, ker_poly) = compute_sequence_of_maps(E, E2, 11)
        sage: compute_isogeny_starks(E1pr, E2pr, 11)
        x^10 + 37*x^9 + 53*x^8 + 66*x^7 + 66*x^6 + 17*x^5 + 57*x^4 + 6*x^3 + 89*x^2 + 53*x + 8

        sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
        sage: R.<x> = GF(37)[]
        sage: f = (x + 14) * (x + 30)
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_isogeny_starks(E, E2, 5)
        x^4 + 14*x^3 + x^2 + 34*x + 21
        sage: f**2
        x^4 + 14*x^3 + x^2 + 34*x + 21

        sage: E = EllipticCurve(QQ, [0,0,0,1,0])
        sage: R.<x> = QQ[]
        sage: f = x
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_isogeny_starks(E, E2, 2)
        x
    """
    K = E1.base_field()
    R = PolynomialRing(K, 'x')
    x = R.gen()

    wp1 = E1.weierstrass_p(prec=4*ell+4)  #BMSS2006 claim 2*ell is enough, but it is not M09
    wp2 = E2.weierstrass_p(prec=4*ell+4)

    # viewed them as power series in Z = z^2
    S = LaurentSeriesRing(K, 'Z')
    Z = S.gen()
    pe1 = 1/Z
    pe2 = 1/Z
    for i in range(2 * ell + 1):
        pe1 += wp1[2*i] * Z**i
        pe2 += wp2[2*i] * Z**i
    pe1 = pe1.add_bigoh(2*ell+2)
    pe2 = pe2.add_bigoh(2*ell+2)

    n = 1
    q = [R(1), R(0)]
    T = pe2

    while ( q[n].degree() < (ell-1) ):
        n += 1
        a_n = 0
        r = -T.valuation()
        while (0 <= r):
            t_r = T[-r]
            a_n = a_n + t_r * x**r
            T = T - t_r*pe1**r
            r = -T.valuation()

        q_n = a_n*q[n-1] + q[n-2]
        q.append(q_n)

        if (n == ell+1 or T == 0):
            if (T == 0 or T.valuation()<2):
                raise ValueError("The two curves are not linked by a cyclic normalized isogeny of degree %s" % ell)
            break

        T = 1/T

    qn = q[n]
    qn = (1/qn.leading_coefficient())*qn

    return qn

def split_kernel_polynomial(poly):
    r"""
    Internal helper function for ``compute_isogeny_kernel_polynomial``.

    INPUT:

    - ``poly`` -- a nonzero univariate polynomial.

    OUTPUT:

    The maximum separable divisor of ``poly``.  If the input is a full
    kernel polynomial where the roots which are `x`-coordinates of
    points of order greater than 2 have multiplicity 1, the output
    will be a polynomial with the same roots, all of multiplicity 1.

    EXAMPLES:

    The following example implicitly exercises this function::

        sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
        sage: R.<x> = GF(37)[]
        sage: f = (x + 10) * (x + 12) * (x + 16)
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_starks
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import split_kernel_polynomial
        sage: ker_poly = compute_isogeny_starks(E, E2, 7); ker_poly
        x^6 + 2*x^5 + 20*x^4 + 11*x^3 + 36*x^2 + 35*x + 16
        sage: ker_poly.factor()
        (x + 10)^2 * (x + 12)^2 * (x + 16)^2
        sage: poly = split_kernel_polynomial(ker_poly); poly
        x^3 + x^2 + 28*x + 33
        sage: poly.factor()
        (x + 10) * (x + 12) * (x + 16)
    """
    from sage.misc.misc_c import prod
    return prod([p for p,e in poly.squarefree_decomposition()])

def compute_isogeny_kernel_polynomial(E1, E2, ell, algorithm="starks"):
    r"""
    Return the kernel polynomial of an isogeny of degree ``ell``
    between ``E1`` and ``E2``.

    INPUT:

    - ``E1``        - an elliptic curve in short Weierstrass form.

    - ``E2``        - an elliptic curve in short Weierstrass form.

    - ``ell``       - the degree of the isogeny from ``E1`` to ``E2``.

    - ``algorithm`` - currently only ``starks`` (default) is implemented.

    OUTPUT:

    polynomial over the field of definition of ``E1``, ``E2``, that is
    the kernel polynomial of the isogeny from ``E1`` to ``E2``.


    .. NOTE::

       If there is no degree ``ell``, cyclic, separable, normalized
       isogeny from ``E1`` to ``E2`` then an error will be raised.


    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_kernel_polynomial

        sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
        sage: R.<x> = GF(37)[]
        sage: f = (x + 14) * (x + 30)
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_isogeny_kernel_polynomial(E, E2, 5)
        x^2 + 7*x + 13
        sage: f
        x^2 + 7*x + 13

        sage: R.<x> = QQ[]
        sage: K.<i> = NumberField(x^2 + 1)
        sage: E = EllipticCurve(K, [0,0,0,1,0])
        sage: E2 = EllipticCurve(K, [0,0,0,16,0])
        sage: compute_isogeny_kernel_polynomial(E, E2, 4)
        x^3 + x
    """
    return split_kernel_polynomial(compute_isogeny_starks(E1, E2, ell))

def compute_intermediate_curves(E1, E2):
    r"""
    Return intermediate curves and isomorphisms.

    .. NOTE::

       This is used so we can compute `\wp` functions from the short
       Weierstrass model more easily.

    .. WARNING::

       The base field must be of characteristic not equal to 2,3.

    INPUT:

    - ``E1`` - an elliptic curve
    - ``E2`` - an elliptic curve

    OUTPUT:

    tuple (``pre_isomorphism``, ``post_isomorphism``,
    ``intermediate_domain``, ``intermediate_codomain``):

    - ``intermediate_domain``: a short Weierstrass model isomorphic to
      ``E1``

    - ``intermediate_codomain``: a short Weierstrass model isomorphic
      to ``E2``

    - ``pre_isomorphism``: normalized isomorphism from ``E1`` to
      intermediate_domain

    - ``post_isomorphism``: normalized isomorphism from
      intermediate_codomain to ``E2``

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_intermediate_curves
        sage: E = EllipticCurve(GF(83), [1,0,1,1,0])
        sage: R.<x> = GF(83)[]; f = x+24
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_intermediate_curves(E, E2)
        (Elliptic Curve defined by y^2 = x^3 + 62*x + 74 over Finite Field of size 83,
         Elliptic Curve defined by y^2 = x^3 + 65*x + 69 over Finite Field of size 83,
         Elliptic-curve morphism:
          From: Elliptic Curve defined by y^2 + x*y + y = x^3 + x over Finite Field of size 83
          To:   Elliptic Curve defined by y^2 = x^3 + 62*x + 74 over Finite Field of size 83
          Via:  (u,r,s,t) = (1, 76, 41, 3),
         Elliptic-curve morphism:
          From: Elliptic Curve defined by y^2 = x^3 + 65*x + 69 over Finite Field of size 83
          To:   Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + 16 over Finite Field of size 83
          Via:  (u,r,s,t) = (1, 7, 42, 42))

        sage: R.<x> = QQ[]
        sage: K.<i> = NumberField(x^2 + 1)
        sage: E = EllipticCurve(K, [0,0,0,1,0])
        sage: E2 = EllipticCurve(K, [0,0,0,16,0])
        sage: compute_intermediate_curves(E, E2)
        (Elliptic Curve defined by y^2 = x^3 + x over Number Field in i with defining polynomial x^2 + 1,
         Elliptic Curve defined by y^2 = x^3 + 16*x over Number Field in i with defining polynomial x^2 + 1,
         Elliptic-curve endomorphism of Elliptic Curve defined by y^2 = x^3 + x over Number Field in i with defining polynomial x^2 + 1
          Via:  (u,r,s,t) = (1, 0, 0, 0),
         Elliptic-curve endomorphism of Elliptic Curve defined by y^2 = x^3 + 16*x over Number Field in i with defining polynomial x^2 + 1
          Via:  (u,r,s,t) = (1, 0, 0, 0))

    """
    if (E1.base_ring().characteristic() in [2,3]):
        raise NotImplementedError("compute_intermediate_curves is only defined for characteristics not 2 or 3")

    # We cannot just use
    # E1w = E1.short_weierstrass_model()
    # E2w = E2.short_weierstrass_model()
    # as the resulting isomorphisms would not be normalised (u=1)

    c4, c6 = E1.c_invariants()
    E1w = EllipticCurve([0,0,0,-c4/48, -c6/864])
    c4, c6 = E2.c_invariants()
    E2w = EllipticCurve([0,0,0,-c4/48, -c6/864])

    # We cannot even just use pre_iso = E1.isomorphism_to(E1w) since
    # it may have u=-1; similarly for E2

    urst = [w for w in isomorphisms(E1,E1w) if w[0]==1][0]
    pre_iso = WeierstrassIsomorphism(E1,urst,E1w)
    urst = [w for w in isomorphisms(E2w,E2) if w[0]==1][0]
    post_iso = WeierstrassIsomorphism(E2w,urst,E2)
    return (E1w, E2w, pre_iso, post_iso)

def compute_sequence_of_maps(E1, E2, ell):
    r"""
    Return intermediate curves, isomorphisms and kernel polynomial.

    INPUT:

    - ``E1``, ``E2`` -- elliptic curves.

    - ``ell`` -- a prime such that there is a degree ``ell`` separable
      normalized isogeny from ``E1`` to ``E2``.

    OUTPUT:

    (pre_isom, post_isom, E1pr, E2pr, ker_poly) where:

    - ``E1pr`` is an elliptic curve in short Weierstrass form
      isomorphic to ``E1``;

    - ``E2pr`` is an elliptic curve in short Weierstrass form
      isomorphic to ``E2``;

    - ``pre_isom`` is a normalised isomorphism from ``E1`` to
      ``E1pr``;

    - ``post_isom`` is a normalised isomorphism from ``E2pr`` to
      ``E2``;

    - ``ker_poly`` is the kernel polynomial of an ``ell``-isogeny from
      ``E1pr`` to ``E2pr``.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_sequence_of_maps
        sage: E = EllipticCurve('11a1')
        sage: R.<x> = QQ[]; f = x^2 - 21*x + 80
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_sequence_of_maps(E, E2, 5)
        (Elliptic-curve morphism:
          From: Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
          To:   Elliptic Curve defined by y^2 = x^3 - 31/3*x - 2501/108 over Rational Field
          Via:  (u,r,s,t) = (1, 1/3, 0, -1/2),
         Elliptic-curve morphism:
          From: Elliptic Curve defined by y^2 = x^3 - 23461/3*x - 28748141/108 over Rational Field
          To:   Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
          Via:  (u,r,s,t) = (1, -1/3, 0, 1/2),
         Elliptic Curve defined by y^2 = x^3 - 31/3*x - 2501/108 over Rational Field,
         Elliptic Curve defined by y^2 = x^3 - 23461/3*x - 28748141/108 over Rational Field,
         x^2 - 61/3*x + 658/9)

        sage: K.<i> = NumberField(x^2 + 1)
        sage: E = EllipticCurve(K, [0,0,0,1,0])
        sage: E2 = EllipticCurve(K, [0,0,0,16,0])
        sage: compute_sequence_of_maps(E, E2, 4)
        (Elliptic-curve endomorphism of Elliptic Curve defined by y^2 = x^3 + x over Number Field in i with defining polynomial x^2 + 1
          Via:  (u,r,s,t) = (1, 0, 0, 0),
         Elliptic-curve endomorphism of Elliptic Curve defined by y^2 = x^3 + 16*x over Number Field in i with defining polynomial x^2 + 1
          Via:  (u,r,s,t) = (1, 0, 0, 0),
         Elliptic Curve defined by y^2 = x^3 + x over Number Field in i with defining polynomial x^2 + 1,
         Elliptic Curve defined by y^2 = x^3 + 16*x over Number Field in i with defining polynomial x^2 + 1,
         x^3 + x)

        sage: E = EllipticCurve(GF(97), [1,0,1,1,0])
        sage: R.<x> = GF(97)[]; f = x^5 + 27*x^4 + 61*x^3 + 58*x^2 + 28*x + 21
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_sequence_of_maps(E, E2, 11)
        (Elliptic-curve morphism:
          From: Elliptic Curve defined by y^2 + x*y + y = x^3 + x over Finite Field of size 97
          To:   Elliptic Curve defined by y^2 = x^3 + 52*x + 31 over Finite Field of size 97
          Via:  (u,r,s,t) = (1, 8, 48, 44),
         Elliptic-curve morphism:
          From: Elliptic Curve defined by y^2 = x^3 + 41*x + 66 over Finite Field of size 97
          To:   Elliptic Curve defined by y^2 + x*y + y = x^3 + 87*x + 26 over Finite Field of size 97
          Via:  (u,r,s,t) = (1, 89, 49, 49),
         Elliptic Curve defined by y^2 = x^3 + 52*x + 31 over Finite Field of size 97,
         Elliptic Curve defined by y^2 = x^3 + 41*x + 66 over Finite Field of size 97,
         x^5 + 67*x^4 + 13*x^3 + 35*x^2 + 77*x + 69)
    """
    (E1pr, E2pr, pre_isom, post_isom) = compute_intermediate_curves(E1, E2)

    ker_poly = compute_isogeny_kernel_polynomial(E1pr, E2pr, ell)

    return (pre_isom, post_isom, E1pr, E2pr, ker_poly)


# Utility function for manipulating isogeny degree matrices


def fill_isogeny_matrix(M):
    """
    Return a filled isogeny matrix giving all degrees from one giving only prime degrees.

    INPUT:

    - ``M`` -- a square symmetric matrix whose off-diagonal `i`, `j`
      entry is either a prime `l` (if the `i`'th and `j`'th curves
      have an `l`-isogeny between them), otherwise is 0.

    OUTPUT:

    (matrix) a square matrix with entries `1` on the diagonal, and in
    general the `i`, `j` entry is `d>0` if `d` is the minimal degree
    of an isogeny from the `i`'th to the `j`'th curve,

    EXAMPLES::

        sage: M = Matrix([[0, 2, 3, 3, 0, 0], [2, 0, 0, 0, 3, 3], [3, 0, 0, 0, 2, 0], [3, 0, 0, 0, 0, 2], [0, 3, 2, 0, 0, 0], [0, 3, 0, 2, 0, 0]]); M
        [0 2 3 3 0 0]
        [2 0 0 0 3 3]
        [3 0 0 0 2 0]
        [3 0 0 0 0 2]
        [0 3 2 0 0 0]
        [0 3 0 2 0 0]
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import fill_isogeny_matrix
        sage: fill_isogeny_matrix(M)
        [ 1  2  3  3  6  6]
        [ 2  1  6  6  3  3]
        [ 3  6  1  9  2 18]
        [ 3  6  9  1 18  2]
        [ 6  3  2 18  1  9]
        [ 6  3 18  2  9  1]
    """
    from sage.matrix.all import Matrix
    from sage.rings.infinity import Infinity

    n = M.nrows()
    M0 = copy(M)
    for i in range(n):
        M0[i,i]=1

    def fix(d):
        if d == 0:
            return Infinity
        return d

    def fix2(d):
        if d == Infinity:
            return 0
        return d

    def pr(M1, M2):
        return Matrix([[fix2(min([fix(M1[i,k]*M2[k,j]) for k in range(n)])) for i in range(n)] for j in range(n)])

    M1 = M0
    M2 = pr(M0,M1)
    while M1!=M2:
        M1 = M2
        M2 = pr(M0,M1)

    return M1

def unfill_isogeny_matrix(M):
    """
    Reverses the action of ``fill_isogeny_matrix``.

    INPUT:

    - ``M`` -- a square symmetric matrix of integers.

    OUTPUT:

    (matrix) a square symmetric matrix obtained from ``M`` by
    replacing non-prime entries with `0`.

    EXAMPLES::

        sage: M = Matrix([[0, 2, 3, 3, 0, 0], [2, 0, 0, 0, 3, 3], [3, 0, 0, 0, 2, 0], [3, 0, 0, 0, 0, 2], [0, 3, 2, 0, 0, 0], [0, 3, 0, 2, 0, 0]]); M
        [0 2 3 3 0 0]
        [2 0 0 0 3 3]
        [3 0 0 0 2 0]
        [3 0 0 0 0 2]
        [0 3 2 0 0 0]
        [0 3 0 2 0 0]
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import fill_isogeny_matrix, unfill_isogeny_matrix
        sage: M1 = fill_isogeny_matrix(M); M1
        [ 1  2  3  3  6  6]
        [ 2  1  6  6  3  3]
        [ 3  6  1  9  2 18]
        [ 3  6  9  1 18  2]
        [ 6  3  2 18  1  9]
        [ 6  3 18  2  9  1]
        sage: unfill_isogeny_matrix(M1)
        [0 2 3 3 0 0]
        [2 0 0 0 3 3]
        [3 0 0 0 2 0]
        [3 0 0 0 0 2]
        [0 3 2 0 0 0]
        [0 3 0 2 0 0]
        sage: unfill_isogeny_matrix(M1) == M
        True
    """
    n = M.nrows()
    M1 = copy(M)
    zero = Integer(0)
    for i in range(n):
        M1[i, i] = zero
        for j in range(i):
            if not M1[i, j].is_prime():
                M1[i, j] = zero
                M1[j, i] = zero
    return M1
