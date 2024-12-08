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
    sage: E = EllipticCurve(k, [1,1])
    sage: Q = E(6,5)
    sage: phi = E.isogeny(Q)
    sage: phi
    Isogeny of degree 7
     from Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 11
       to Elliptic Curve defined by y^2 = x^3 + 7*x + 8
    over Finite Field of size 11
    sage: P = E(4,5)
    sage: phi(P)
    (10 : 0 : 1)
    sage: phi.codomain()
    Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 11
    sage: phi.rational_maps()
    ((x^7 + 4*x^6 - 3*x^5 - 2*x^4
       - 3*x^3 + 3*x^2 + x - 2)/(x^6 + 4*x^5 - 4*x^4 - 5*x^3 + 5*x^2),
     (x^9*y - 5*x^8*y - x^7*y + x^5*y - x^4*y
       - 5*x^3*y - 5*x^2*y - 2*x*y - 5*y)/(x^9 - 5*x^8 + 4*x^6 - 3*x^4 + 2*x^3))

The methods directly accessible from an elliptic curve ``E`` over a
field are
:meth:`~sage.schemes.elliptic_curves.ell_field.EllipticCurve_field.isogeny`
and
:meth:`~sage.schemes.elliptic_curves.ell_field.EllipticCurve_field.isogeny_codomain`.

The most useful methods that apply to isogenies are:

- ``.domain()``
- ``.codomain()``
- :meth:`~EllipticCurveHom.degree`
- :meth:`~EllipticCurveIsogeny.dual`
- :meth:`~EllipticCurveIsogeny.rational_maps`
- :meth:`~EllipticCurveIsogeny.kernel_polynomial`

.. WARNING::

    This class only implements separable isogenies. When using Kohel's
    algorithm, only cyclic isogenies can be computed (except for `[2]`).

    Working with other kinds of isogenies may be possible using other
    child classes of :class:`EllipticCurveHom`.

    Some algorithms may need the isogeny to be normalized.

AUTHORS:

- Daniel Shumow <shumow@gmail.com>: 2009-04-19: initial version

- Chris Wuthrich: 7/09: add check of input, not the full list is needed.
  10/09: eliminating some bugs.

- John Cremona 2014-08-08: tidying of code and docstrings, systematic
  use of univariate vs. bivariate polynomials and rational functions.

- Lorenz Panny (2022-04): major cleanup of code and documentation

- Lorenz Panny (2022): inseparable duals

- Rémy Oudompheng (2023): implementation of the BMSS algorithm
"""

# ****************************************************************************
#       Copyright (C) 2009 Daniel Shumow <shumow@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy, deepcopy

from sage.misc.cachefunc import cached_method
from sage.structure.sequence import Sequence

from sage.schemes.elliptic_curves.hom import EllipticCurveHom

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer import Integer
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.fraction_field import FractionField

from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic

from sage.schemes.elliptic_curves.weierstrass_morphism \
        import WeierstrassIsomorphism, _isomorphisms, baseWI, negation_morphism

#
# Private function for parsing input to determine the type of
# algorithm
#


def _isogeny_determine_algorithm(E, kernel):
    r"""
    Helper function to infer the algorithm to be used from the
    parameters passed to the various isogeny functions.

    If ``kernel`` is a list of points on the elliptic curve `E`,
    we will try to use Vélu's algorithm.

    If ``kernel`` is a list of coefficients or a univariate
    polynomial, we will try to use the Kohel's algorithms.

    INPUT:

    - ``E`` -- domain elliptic curve

    - ``kernel`` -- either a list of points on ``E``, or a univariate
      polynomial or list of coefficients of a univariate polynomial

    OUTPUT: string; either ``'velu'`` or ``'kohel'``.

    EXAMPLES:

    This helper function will be implicitly called by the following examples::

        sage: R.<x> = GF(5)[]
        sage: E = EllipticCurve(GF(5), [0,0,0,1,0])  # indirect doctest

    We can construct the same isogeny from a kernel polynomial::

        sage: phi = EllipticCurveIsogeny(E, x + 3)  # indirect doctest

    or from a list of coefficients of a kernel polynomial::

        sage: phi == EllipticCurveIsogeny(E, [3,1])  # indirect doctest
        True

    or from a rational point which generates the kernel::

        sage: phi == EllipticCurveIsogeny(E, E((2,0)))  # indirect doctest
        True

    In the first two cases, Kohel's algorithm will be used, while in
    the third case it is Vélu::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import _isogeny_determine_algorithm
        sage: _isogeny_determine_algorithm(E, x + 3)
        'kohel'
        sage: _isogeny_determine_algorithm(E, [3, 1])
        'kohel'
        sage: _isogeny_determine_algorithm(E, E((2,0)))
        'velu'
    """
    kernel_is_list = isinstance(kernel, list)

    if not kernel_is_list and kernel in E:
        kernel = [kernel]
        kernel_is_list = True

    if isinstance(kernel, Polynomial) or (kernel_is_list and kernel[0] in E.base_ring()):
        return "kohel"

    if kernel_is_list and kernel[0] in E:
        # note that if kernel[0] is on an extension of E this
        # condition will be false
        return "velu"

    raise ValueError("invalid parameters to EllipticCurveIsogeny constructor")


def isogeny_codomain_from_kernel(E, kernel):
    r"""
    Compute the isogeny codomain given a kernel.

    INPUT:

    - ``E`` -- domain elliptic curve

    - ``kernel`` -- either a list of points in the kernel of the isogeny,
                    or a kernel polynomial (specified as either a
                    univariate polynomial or a coefficient list)

    OUTPUT: elliptic curve) The codomain of the separable normalized isogeny
    defined by this kernel.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogeny_codomain_from_kernel
        sage: E = EllipticCurve(GF(7), [1,0,1,0,1])
        sage: R.<x> = GF(7)[]
        sage: isogeny_codomain_from_kernel(E, [4,1])
        Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + 6
         over Finite Field of size 7
        sage: (EllipticCurveIsogeny(E, [4,1]).codomain()
        ....:  == isogeny_codomain_from_kernel(E, [4,1]))
        True
        sage: isogeny_codomain_from_kernel(E, x^3 + x^2 + 4*x + 3)
        Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + 6
         over Finite Field of size 7
        sage: isogeny_codomain_from_kernel(E, x^3 + 2*x^2 + 4*x + 3)
        Elliptic Curve defined by y^2 + x*y + y = x^3 + 5*x + 2
         over Finite Field of size 7

        sage: # needs sage.rings.finite_rings
        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: kernel_list = [E((15,10)), E((10,3)), E((6,5))]
        sage: isogeny_codomain_from_kernel(E, kernel_list)
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 3*x + 15
         over Finite Field of size 19
    """
    algorithm = _isogeny_determine_algorithm(E, kernel)

    if algorithm == 'velu':
        # if we are using Velu's formula, just instantiate the isogeny
        # and return the codomain
        return EllipticCurveIsogeny(E, kernel).codomain()

    if algorithm == 'kohel':
        return compute_codomain_kohel(E, kernel)

    raise NotImplementedError


def compute_codomain_formula(E, v, w):
    r"""
    Compute the codomain curve given parameters `v` and `w` (as in
    Vélu/Kohel/etc. formulas).

    INPUT:

    - ``E`` -- an elliptic curve

    - ``v``, ``w`` -- elements of the base field of ``E``

    OUTPUT:

    The elliptic curve with invariants
    `[a_1,a_2,a_3,a_4-5v,a_6-(a_1^2+4a_2)v-7w]` where
    `E = [a_1,a_2,a_3,a_4,a_6]`.

    EXAMPLES:

    This formula is used by every invocation of the
    :class:`EllipticCurveIsogeny` constructor::

        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: phi = EllipticCurveIsogeny(E, E((1,2)) )
        sage: phi.codomain()
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 9*x + 13
         over Finite Field of size 19
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_codomain_formula
        sage: v = phi._EllipticCurveIsogeny__v
        sage: w = phi._EllipticCurveIsogeny__w
        sage: compute_codomain_formula(E, v, w) == phi.codomain()
        True
    """
    a1, a2, a3, a4, a6 = E.a_invariants()

    A4 = a4 - 5*v
    A6 = a6 - (a1**2 + 4*a2)*v - 7*w

    return EllipticCurve([a1, a2, a3, A4, A6])


def compute_vw_kohel_even_deg1(x0, y0, a1, a2, a4):
    r"""
    Compute Vélu's `(v,w)` using Kohel's formulas for isogenies of
    degree exactly divisible by `2`.

    INPUT:

    - ``x0``, ``y0`` -- coordinates of a 2-torsion point on an elliptic curve `E`

    - ``a1``, ``a2``, ``a4`` -- invariants of `E`

    OUTPUT: Vélu's isogeny parameters `(v,w)`.

    EXAMPLES:

    This function will be implicitly called by the following example::

        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: phi = EllipticCurveIsogeny(E, [9,1]); phi
        Isogeny of degree 2
         from Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5
              over Finite Field of size 19
           to Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 9*x + 8
              over Finite Field of size 19
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_vw_kohel_even_deg1
        sage: a1,a2,a3,a4,a6 = E.a_invariants()
        sage: x0 = -9
        sage: y0 = -(a1*x0 + a3)/2
        sage: compute_vw_kohel_even_deg1(x0, y0, a1, a2, a4)
        (18, 9)
    """
    v = 3*x0**2 + 2*a2*x0 + a4 - a1*y0
    w = x0 * v
    return v, w


def compute_vw_kohel_even_deg3(b2, b4, s1, s2, s3):
    r"""
    Compute Vélu's `(v,w)` using Kohel's formulas for isogenies of
    degree divisible by `4`.

    INPUT:

    - ``b2``, ``b4`` -- invariants of an elliptic curve `E`

    - ``s1``, ``s2``, ``s3`` -- signed coefficients of the 2-division
      polynomial of `E`

    OUTPUT: Vélu's isogeny parameters `(v,w)`.

    EXAMPLES:

    This function will be implicitly called by the following example::

        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: R.<x> = GF(19)[]
        sage: phi = EllipticCurveIsogeny(E, x^3 + 7*x^2 + 15*x + 12); phi
        Isogeny of degree 4
         from Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5
              over Finite Field of size 19
           to Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 3*x + 15
              over Finite Field of size 19
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_vw_kohel_even_deg3
        sage: b2,b4 = E.b2(), E.b4()
        sage: s1, s2, s3 = -7, 15, -12
        sage: compute_vw_kohel_even_deg3(b2, b4, s1, s2, s3)
        (4, 7)
    """
    temp1 = s1**2 - 2*s2
    v = 3*temp1 + (b2*s1 + 3*b4)/2
    w = 3*(s1**3 - 3*s1*s2 + 3*s3) + (b2*temp1 + b4*s1)/2
    return v, w


def compute_vw_kohel_odd(b2, b4, b6, s1, s2, s3, n):
    r"""
    Compute Vélu's `(v,w)` using Kohel's formulas for isogenies of odd
    degree.

    INPUT:

    - ``b2``, ``b4``, ``b6`` -- invariants of an elliptic curve `E`

    - ``s1``, ``s2``, ``s3`` -- signed coefficients of lowest powers
      of `x` in the kernel polynomial

    - ``n`` -- integer; the degree

    OUTPUT: Vélu's isogeny parameters `(v,w)`

    EXAMPLES:

    This function will be implicitly called by the following example::

        sage: E = EllipticCurve(GF(19), [18,17,16,15,14])
        sage: R.<x> = GF(19)[]
        sage: phi = EllipticCurveIsogeny(E, x^3 + 14*x^2 + 3*x + 11); phi
        Isogeny of degree 7
         from Elliptic Curve defined by y^2 + 18*x*y + 16*y = x^3 + 17*x^2 + 15*x + 14
              over Finite Field of size 19
           to Elliptic Curve defined by y^2 + 18*x*y + 16*y = x^3 + 17*x^2 + 18*x + 18
              over Finite Field of size 19
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_vw_kohel_odd
        sage: b2,b4,b6 = E.b2(), E.b4(), E.b6()
        sage: s1,s2,s3 = -14,3,-11
        sage: compute_vw_kohel_odd(b2,b4,b6,s1,s2,s3,3)
        (7, 1)
    """
    v = 6*(s1**2 - 2*s2) + b2*s1 + n*b4
    w = 10*(s1**3 - 3*s1*s2 + 3*s3) + 2*b2*(s1**2 - 2*s2) + 3*b4*s1 + n*b6
    return v, w


def compute_codomain_kohel(E, kernel):
    r"""
    Compute the codomain from the kernel polynomial using Kohel's
    formulas.

    INPUT:

    - ``E`` -- domain elliptic curve

    - ``kernel`` -- polynomial or list; the kernel polynomial, or a
      list of its coefficients

    OUTPUT: elliptic curve; the codomain elliptic curve of the isogeny
    defined by ``kernel``

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_codomain_kohel
        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: phi = EllipticCurveIsogeny(E, [9,1])
        sage: phi.codomain() == isogeny_codomain_from_kernel(E, [9,1])
        True
        sage: compute_codomain_kohel(E, [9,1])
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 9*x + 8
         over Finite Field of size 19
        sage: R.<x> = GF(19)[]
        sage: E = EllipticCurve(GF(19), [18,17,16,15,14])
        sage: phi = EllipticCurveIsogeny(E, x^3 + 14*x^2 + 3*x + 11)
        sage: phi.codomain() == isogeny_codomain_from_kernel(E, x^3 + 14*x^2 + 3*x + 11)
        True
        sage: compute_codomain_kohel(E, x^3 + 14*x^2 + 3*x + 11)
        Elliptic Curve defined by y^2 + 18*x*y + 16*y = x^3 + 17*x^2 + 18*x + 18
         over Finite Field of size 19
        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: phi = EllipticCurveIsogeny(E, x^3 + 7*x^2 + 15*x + 12)
        sage: isogeny_codomain_from_kernel(E, x^3 + 7*x^2 + 15*x + 12) == phi.codomain()
        True
        sage: compute_codomain_kohel(E, x^3 + 7*x^2 + 15*x + 12)
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 3*x + 15
         over Finite Field of size 19

    ALGORITHM:

    This function uses the formulas of Section 2.4 of [Koh1996]_.
    """
    # First set up the polynomial ring
    base_field = E.base_ring()
    poly_ring = PolynomialRing(base_field,'x')

    try:
        psi = poly_ring(kernel)
    except TypeError:
        raise ValueError("invalid input to compute_codomain_kohel")

    # next determine the even / odd part of the isogeny
    psi_2tor = two_torsion_part(E, psi)

    if psi_2tor.degree() != 0: # even degree case

        psi_quo = psi//psi_2tor

        if psi_quo.degree() != 0:
            raise NotImplementedError("Kohel's algorithm currently only supports cyclic isogenies (except for [2])")

        n = psi_2tor.degree()

        if n == 1: # degree divisible exactly by 2

            a1, a2, a3, a4, a6 = E.a_invariants()

            x0 = -psi_2tor.constant_coefficient()

            # determine y0
            if base_field.characteristic() == 2:
                y0 = (x0**3 + a2*x0**2 + a4*x0 + a6).sqrt()
            else:
                y0 = -(a1*x0 + a3)/2

            # now (x0,y0) is the 2-torsion point in the kernel

            v, w = compute_vw_kohel_even_deg1(x0, y0, a1, a2, a4)

        elif n == 3:  # psi_2tor is the full 2-division polynomial

            b2, b4, _, _ = E.b_invariants()

            s1 = -psi_2tor[n - 1]
            s2 = psi_2tor[n - 2]
            s3 = -psi_2tor[n - 3]

            v, w = compute_vw_kohel_even_deg3(b2, b4, s1, s2, s3)

    else:  # odd degree case

        n = psi.degree()

        b2, b4, b6, _ = E.b_invariants()

        s1 = -psi[n - 1] if n >= 1 else 0
        s2 = psi[n - 2] if n >= 2 else 0
        s3 = -psi[n - 3] if n >= 3 else 0

        # initializing these allows us to calculate E2.
        v, w = compute_vw_kohel_odd(b2, b4, b6, s1, s2, s3, n)

    return compute_codomain_formula(E, v, w)


def two_torsion_part(E, psi):
    r"""
    Return the greatest common divisor of ``psi`` and the 2-torsion
    polynomial of `E`.

    INPUT:

    - ``E`` -- an elliptic curve

    - ``psi`` -- a univariate polynomial over the base field of ``E``

    OUTPUT: polynomial; the `\gcd` of ``psi`` and the 2-torsion polynomial of ``E``

    EXAMPLES:

    Every function that computes the kernel polynomial via Kohel's
    formulas will call this function::

        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: R.<x> = GF(19)[]
        sage: phi = EllipticCurveIsogeny(E, x + 13)
        sage: isogeny_codomain_from_kernel(E, x + 13) == phi.codomain()
        True
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import two_torsion_part
        sage: two_torsion_part(E, x + 13)
        x + 13
    """
    x = psi.parent().gen() # NB psi is univariate but could be constant
    psi_2 = E.two_division_polynomial(x)
    return psi.gcd(psi_2)


class EllipticCurveIsogeny(EllipticCurveHom):
    r"""
    This class implements separable isogenies of elliptic curves.

    Several different algorithms for computing isogenies are
    available.  These include:

    - Vélu's Formulas: Vélu's original formulas for computing
      isogenies.  This algorithm is selected by giving as the
      ``kernel`` parameter a single point, or a list of points,
      generating a finite subgroup.

    - Kohel's Formulas: Kohel's original formulas for computing
      isogenies.  This algorithm is selected by giving as the
      ``kernel`` parameter a monic polynomial (or a coefficient list)
      which will define the kernel of the isogeny.
      Kohel's algorithm is currently only implemented for cyclic
      isogenies, with the exception of `[2]`.

    INPUT:

    - ``E`` -- an elliptic curve; the domain of the isogeny to initialize

    - ``kernel`` -- a kernel; either a point on ``E``, a list of
      points on ``E``, a monic kernel polynomial, or ``None``.
      If initializing from a domain/codomain, this must be ``None``.

    - ``codomain`` -- an elliptic curve (default: ``None``)

      - If ``kernel`` is ``None``, then ``degree`` must be given as well
        and the given ``codomain`` must be the codomain of a cyclic,
        separable, normalized isogeny of the given degree.

      - If ``kernel`` is not ``None``, then this must be isomorphic to
        the codomain of the separable isogeny defined by ``kernel``; in
        this case, the isogeny is post-composed with an isomorphism so
        that the codomain equals the given curve.

    - ``degree`` -- integer (default: ``None``)

      - If ``kernel`` is ``None``, then this is the degree of the isogeny
        from ``E`` to ``codomain``.

      - If ``kernel`` is not ``None``, then this is used to determine
        whether or not to skip a `\gcd` of the given kernel polynomial
        with the two-torsion polynomial of ``E``.

    - ``model`` -- string (default: ``None``); supported values
      (cf. :func:`~sage.schemes.elliptic_curves.ell_field.compute_model`):

      - ``'minimal'`` -- if ``E`` is a curve over the rationals or
        over a number field, then the codomain is a global minimal
        model where this exists.

      - ``'short_weierstrass'`` -- the codomain is a short Weierstrass curve,
        assuming one exists.

      - ``'montgomery'`` -- the codomain is an (untwisted) Montgomery
        curve, assuming one exists over this field.

    - ``check`` -- boolean (default: ``True``); check whether the input is valid.
      Setting this to ``False`` can lead to significant speedups.

    EXAMPLES:

    A simple example of creating an isogeny of a field of small
    characteristic::

        sage: E = EllipticCurve(GF(7), [0,0,0,1,0])
        sage: phi = EllipticCurveIsogeny(E, E((0,0)) ); phi
        Isogeny of degree 2
         from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7
           to Elliptic Curve defined by y^2 = x^3 + 3*x over Finite Field of size 7
        sage: phi.degree() == 2
        True
        sage: phi.kernel_polynomial()
        x
        sage: phi.rational_maps()
        ((x^2 + 1)/x, (x^2*y - y)/x^2)
        sage: phi == loads(dumps(phi))          # known bug
        True

    A more complicated example of a characteristic-2 field::

        sage: # needs sage.rings.finite_rings
        sage: E = EllipticCurve(GF(2^4,'alpha'), [0,0,1,0,1])
        sage: P = E((1,1))
        sage: phi_v = EllipticCurveIsogeny(E, P); phi_v
        Isogeny of degree 3
         from Elliptic Curve defined by y^2 + y = x^3 + 1
              over Finite Field in alpha of size 2^4
           to Elliptic Curve defined by y^2 + y = x^3
              over Finite Field in alpha of size 2^4
        sage: phi_ker_poly = phi_v.kernel_polynomial()
        sage: phi_ker_poly
        x + 1
        sage: phi_k = EllipticCurveIsogeny(E, phi_ker_poly)
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

    We can create an isogeny whose kernel equals the full 2-torsion::

        sage: # needs sage.rings.finite_rings
        sage: E = EllipticCurve(GF((3,2)), [0,0,0,1,1])
        sage: ker_poly = E.division_polynomial(2)
        sage: phi = EllipticCurveIsogeny(E, ker_poly); phi
        Isogeny of degree 4
         from Elliptic Curve defined by y^2 = x^3 + x + 1
              over Finite Field in z2 of size 3^2
           to Elliptic Curve defined by y^2 = x^3 + x + 1
              over Finite Field in z2 of size 3^2
        sage: P1,P2,P3 = filter(bool, E(0).division_points(2))
        sage: phi(P1)
        (0 : 1 : 0)
        sage: phi(P2)
        (0 : 1 : 0)
        sage: phi(P3)
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
        Isogeny of degree 1
         from Elliptic Curve defined by y^2 + 23*x*y + 22*y = x^3 + x^2 + 7*x + 18
              over Finite Field of size 31
           to Elliptic Curve defined by y^2 + 23*x*y + 22*y = x^3 + x^2 + 7*x + 18
              over Finite Field of size 31
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

    Vélu and Kohel also work in characteristic `0`::

        sage: E = EllipticCurve(QQ, [0,0,0,3,4])
        sage: P_list = E.torsion_points()
        sage: phi = EllipticCurveIsogeny(E, P_list); phi
        Isogeny of degree 2
         from Elliptic Curve defined by y^2 = x^3 + 3*x + 4 over Rational Field
           to Elliptic Curve defined by y^2 = x^3 - 27*x + 46 over Rational Field
        sage: P = E((0,2))
        sage: phi(P)
        (6 : -10 : 1)
        sage: phi_ker_poly = phi.kernel_polynomial()
        sage: phi_ker_poly
        x + 1
        sage: phi_k = EllipticCurveIsogeny(E, phi_ker_poly); phi_k
        Isogeny of degree 2
         from Elliptic Curve defined by y^2 = x^3 + 3*x + 4 over Rational Field
           to Elliptic Curve defined by y^2 = x^3 - 27*x + 46 over Rational Field
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
        Isogeny of degree 5
         from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
           to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
        sage: P = E((16,-61))
        sage: phi_v(P)
        (0 : 1 : 0)
        sage: ker_poly = phi_v.kernel_polynomial(); ker_poly
        x^2 - 21*x + 80
        sage: phi_k = EllipticCurveIsogeny(E, ker_poly); phi_k
        Isogeny of degree 5
         from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
           to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
        sage: phi_k == phi_v
        True
        sage: phi_v(P) == phi_k(P)
        True
        sage: phi_k.is_separable()
        True

    We can also do this same example over the number field defined by
    the irreducible two-torsion polynomial of `E`::

        sage: # needs sage.rings.number_field
        sage: E = EllipticCurve('11a1')
        sage: P_list = E.torsion_points()
        sage: x = polygen(ZZ, 'x')
        sage: K.<alpha> = NumberField(x^3 - 2* x^2 - 40*x - 158)
        sage: EK = E.change_ring(K)
        sage: P_list = [EK(P) for P in P_list]
        sage: phi_v = EllipticCurveIsogeny(EK, P_list); phi_v
        Isogeny of degree 5
         from Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20)
              over Number Field in alpha with defining polynomial x^3 - 2*x^2 - 40*x - 158
           to Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-7820)*x + (-263580)
              over Number Field in alpha with defining polynomial x^3 - 2*x^2 - 40*x - 158
        sage: P = EK((alpha/2,-1/2))
        sage: phi_v(P)
        (122/121*alpha^2 + 1633/242*alpha - 3920/121 : -1/2 : 1)
        sage: ker_poly = phi_v.kernel_polynomial()
        sage: ker_poly
        x^2 - 21*x + 80
        sage: phi_k = EllipticCurveIsogeny(EK, ker_poly); phi_k
        Isogeny of degree 5
         from Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20)
              over Number Field in alpha with defining polynomial x^3 - 2*x^2 - 40*x - 158
           to Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-7820)*x + (-263580)
              over Number Field in alpha with defining polynomial x^3 - 2*x^2 - 40*x - 158
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
        sage: phi_s = EllipticCurveIsogeny(E, None, E2, 5); phi_s
        Isogeny of degree 5
         from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
           to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
        sage: phi_s == phi
        True
        sage: phi_s.rational_maps() == phi.rational_maps()
        True

    However, only cyclic normalized isogenies can be constructed this way.
    The non-cyclic multiplication-by-`3` isogeny won't be found::

        sage: E.isogeny(None, codomain=E, degree=9)
        Traceback (most recent call last):
        ...
        ValueError: the two curves are not linked by a cyclic normalized isogeny of degree 9

    Non-normalized isogeny also won't be found::

        sage: E2.isogeny(None, codomain=E, degree=5)
        Traceback (most recent call last):
        ...
        ValueError: the two curves are not linked by a cyclic normalized isogeny of degree 5
        sage: phihat = phi.dual(); phihat
        Isogeny of degree 5
         from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580
              over Rational Field
           to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
        sage: phihat.is_normalized()
        False

    Here an example of a construction of a endomorphisms with cyclic
    kernel on a CM-curve::

        sage: # needs sage.rings.number_field
        sage: K.<i> = NumberField(x^2 + 1)
        sage: E = EllipticCurve(K, [1,0])
        sage: RK.<X> = K[]
        sage: f = X^2 - 2/5*i + 1/5
        sage: phi= E.isogeny(f)
        sage: isom = phi.codomain().isomorphism_to(E)
        sage: phi = isom * phi
        sage: phi.codomain() == phi.domain()
        True
        sage: phi.rational_maps()
        (((4/25*i + 3/25)*x^5 + (4/5*i - 2/5)*x^3 - x)/(x^4 + (-4/5*i + 2/5)*x^2 + (-4/25*i - 3/25)),
         ((11/125*i + 2/125)*x^6*y + (-23/125*i + 64/125)*x^4*y + (141/125*i + 162/125)*x^2*y + (3/25*i - 4/25)*y)/(x^6 + (-6/5*i + 3/5)*x^4 + (-12/25*i - 9/25)*x^2 + (2/125*i - 11/125)))

    TESTS:

    Domain and codomain tests (see :issue:`12880`)::

        sage: E = EllipticCurve(QQ, [0,0,0,1,0])
        sage: phi = EllipticCurveIsogeny(E, E(0,0))
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

    Composition tests (see :issue:`16245`, cf. :issue:`34410`)::

        sage: E = EllipticCurve(j=GF(7)(0))
        sage: phi = E.isogeny([E(0), E((0,1)), E((0,-1))]); phi
        Composite morphism of degree 3:
          From: Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7
          To:   Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7
        sage: phi2 = phi * phi; phi2
        Composite morphism of degree 9 = 3^2:
          From: Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7
          To:   Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7

    Examples over relative number fields used not to work (see :issue:`16779`)::

        sage: # long time, needs sage.rings.number_field
        sage: pol26 = hilbert_class_polynomial(-4*26)
        sage: F = NumberField(pol26,'a')
        sage: pol = F.optimized_representation()[0].polynomial()
        sage: K.<a> = NumberField(pol)
        sage: j = pol26.roots(K)[0][0]
        sage: E = EllipticCurve(j=j)
        sage: L.<b> = K.extension(x^2 + 26)
        sage: EL = E.change_ring(L)
        sage: iso2 = EL.isogenies_prime_degree(2); len(iso2)
        1
        sage: iso3 = EL.isogenies_prime_degree(3); len(iso3)
        2

    Examples over function fields used not to work (see :issue:`11327`)::

        sage: F.<t> = FunctionField(QQ)
        sage: E = EllipticCurve([0,0,0,-t^2,0])
        sage: isogs = E.isogenies_prime_degree(2)
        sage: isogs[0]
        Isogeny of degree 2
         from Elliptic Curve defined by y^2 = x^3 + (-t^2)*x
              over Rational function field in t over Rational Field
           to Elliptic Curve defined by y^2 = x^3 + 4*t^2*x
              over Rational function field in t over Rational Field
        sage: isogs[0].rational_maps()
        ((x^2 - t^2)/x, (x^2*y + t^2*y)/x^2)
        sage: duals = [phi.dual() for phi in isogs]
        sage: duals[0]
        Isogeny of degree 2
         from Elliptic Curve defined by y^2 = x^3 + 4*t^2*x
              over Rational function field in t over Rational Field
           to Elliptic Curve defined by y^2 = x^3 + (-t^2)*x
              over Rational function field in t over Rational Field
        sage: duals[0].rational_maps()
        ((1/4*x^2 + t^2)/x, (1/8*x^2*y + (-1/2*t^2)*y)/x^2)
        sage: duals[0]
        Isogeny of degree 2
         from Elliptic Curve defined by y^2 = x^3 + 4*t^2*x
              over Rational function field in t over Rational Field
           to Elliptic Curve defined by y^2 = x^3 + (-t^2)*x
              over Rational function field in t over Rational Field
    """

    ####################
    # member variables
    ####################

    #
    # variables common to all algorithms
    #
    _domain = None
    _codomain = None

    _degree = None

    __algorithm = None

    __check = None

    #
    # pre-isomorphism
    #
    __pre_isomorphism = None
    __prei_ratl_maps = None

    #
    # post-isomorphism
    #
    __post_isomorphism = None
    __posti_ratl_maps = None

    #
    # algebraic structs
    #
    __base_field = None
    __poly_ring = None # univariate in x over __base_field
    __mpoly_ring = None # __base_field[x][y], internal use only

    #
    # Rational Maps
    #
    __ratl_maps = None

    #
    # Kernel Data
    #

    __kernel_list = None  # list of elements in the kernel

    __kernel_polynomial = None # polynomial with roots at x values for x-coordinate of points in the kernel

    __inner_kernel_polynomial = None # the inner kernel polynomial (ignoring preisomorphism)

    #
    # member variables common to Velu's formula
    #

    # a full set of representatives of the kernel subgroup modulo negation
    __kernel_mod_sign = None

    # variables used in Velu's formula (as well as Kohel's variant)
    __v = None
    __w = None

    #
    # member variables specific to Kohel's algorithm.
    #
    __psi = None # psi polynomial
    __phi = None # phi polynomial
    __omega = None # omega polynomial, an element of k[x][y]

    # Other members: __xfield, __xyfield

    #
    # Python Special Functions
    #

    def __init__(self, E, kernel, codomain=None, degree=None, model=None, check=True):
        r"""
        Constructor for ``EllipticCurveIsogeny`` class.

        EXAMPLES::

            sage: E = EllipticCurve(GF(2), [0,0,1,0,1])
            sage: phi = EllipticCurveIsogeny(E, [1,1]); phi
            Isogeny of degree 3
             from Elliptic Curve defined by y^2 + y = x^3 + 1 over Finite Field of size 2
               to Elliptic Curve defined by y^2 + y = x^3 over Finite Field of size 2

            sage: E = EllipticCurve(GF(31), [0,0,0,1,0])
            sage: P = E((2,17))
            sage: phi = EllipticCurveIsogeny(E, P); phi
            Isogeny of degree 8
             from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 31
               to Elliptic Curve defined by y^2 = x^3 + 10*x + 28 over Finite Field of size 31

            sage: E = EllipticCurve('17a1')
            sage: phi = EllipticCurveIsogeny(E, [41/3, -55, -1, -1, 1]); phi
            Isogeny of degree 9
             from Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 - x - 14
                  over Rational Field
               to Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 - 56*x - 10124
                  over Rational Field

            sage: E = EllipticCurve('37a1')
            sage: triv = EllipticCurveIsogeny(E, E(0)); triv
            Isogeny of degree 1
             from Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
               to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: triv.rational_maps()
            (x, y)

            sage: E = EllipticCurve('49a3')
            sage: R.<X> = QQ[]
            sage: EllipticCurveIsogeny(E, X^3 - 13*X^2 - 58*X + 503, check=False)
            Isogeny of degree 7
             from Elliptic Curve defined by y^2 + x*y = x^3 - x^2 - 107*x + 552
                  over Rational Field
               to Elliptic Curve defined by y^2 + x*y = x^3 - x^2 - 5252*x - 178837
                  over Rational Field
        """
        if not isinstance(E, EllipticCurve_generic):
            raise ValueError("given E is not an elliptic curve")

        if not isinstance(kernel, list) and kernel in E:
            # a single point was given, we put it in a list
            # the first condition assures that [1,1] is treated as x+1
            kernel = [kernel]

        # if the kernel is None and the codomain isn't
        # calculate the kernel polynomial
        pre_isom = None
        post_isom = None

        self.__check = check

        if kernel is None and codomain is not None:

            if degree is None:
                raise ValueError("degree must be given when specifying isogeny by domain and codomain")

            # save the codomain: really used now (trac #7096)
            old_codomain = codomain

            pre_isom, post_isom, E, codomain, kernel = compute_sequence_of_maps(E, codomain, degree)

        self.__init_algebraic_structs(E)

        algorithm = _isogeny_determine_algorithm(E, kernel)

        self.__algorithm = algorithm

        if algorithm == 'velu':
            self.__init_from_kernel_gens(kernel)
        elif algorithm == 'kohel':
            self.__init_from_kernel_polynomial(kernel)
        else:
            raise NotImplementedError

        self.__compute_codomain()

        self.__setup_post_isomorphism(codomain, model)

        if pre_isom is not None:
            self.__set_pre_isomorphism(pre_isom.domain(), pre_isom)

        if post_isom is not None:
            self.__set_post_isomorphism(old_codomain, post_isom)   #(trac #7096)

        # Inheritance house keeping
        self.__perform_inheritance_housekeeping()

    def _compose_with_isomorphism(self, pre_isomorphism=None, post_isomorphism=None):
        """
        Returns a modified isogeny by pre-composing or post-composing with a
        :class:`sage.schemes.elliptic_curves.weierstrass_morphism.WeierstrassIsomorphism`.

        Old tests from ``__clear_cached_values``::

            sage: R.<x> = QQ[]
            sage: E = EllipticCurve(QQ, [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: old_ratl_maps = phi.rational_maps()
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import negation_morphism
            sage: phi = phi._compose_with_isomorphism(post_isomorphism=negation_morphism(phi.codomain()))
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
            sage: phi = phi._compose_with_isomorphism(post_isomorphism=WeierstrassIsomorphism(phi.codomain(),
            ....:                                                  (-13,13,-13,13)))
            sage: old_ratl_maps == phi.rational_maps()
            False

        Coverage tests::

            sage: E = EllipticCurve(GF(7), [0, 1])
            sage: phi = EllipticCurveIsogeny(E, E.0)
            sage: E2 = phi.codomain()
            sage: phi._compose_with_isomorphism(pre_isomorphism=E2.isomorphism_to(E2))
            Traceback (most recent call last):
            ...
            ValueError: invalid parameter: pre-isomorphism must have codomain curve equal to this isogenies' domain
            sage: phi._compose_with_isomorphism(post_isomorphism=E2.isomorphism_to(E2))
            Isogeny of degree 6 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 6*x + 1 over Finite Field of size 7
            sage: phi._compose_with_isomorphism(pre_isomorphism=E.isomorphism_to(E))
            Isogeny of degree 6 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 6*x + 1 over Finite Field of size 7
            sage: phi._compose_with_isomorphism(post_isomorphism=E.isomorphism_to(E))
            Traceback (most recent call last):
            ...
            ValueError: invalid parameter: post-isomorphism must have domain curve equal to this isogenies' codomain
            sage: phi._compose_with_isomorphism(pre_isomorphism=E2.isogeny(E2.0))
            Traceback (most recent call last):
            ...
            ValueError: invalid parameter: pre-isomorphism must be a WeierstrassIsomorphism
            sage: phi._compose_with_isomorphism(post_isomorphism=E2.isogeny(E2.0))
            Traceback (most recent call last):
            ...
            ValueError: invalid parameter: post-isomorphism must be a WeierstrassIsomorphism
        """
        final_pre_isomorphism = self.__pre_isomorphism
        final_post_isomorphism = self.__post_isomorphism
        if pre_isomorphism is not None:
            if not isinstance(pre_isomorphism, WeierstrassIsomorphism):
                raise ValueError("invalid parameter: pre-isomorphism must be a WeierstrassIsomorphism")
            if self._domain != pre_isomorphism.codomain():
                raise ValueError("invalid parameter: pre-isomorphism must have codomain curve equal to this isogenies' domain")
            final_pre_isomorphism = final_pre_isomorphism * pre_isomorphism if final_pre_isomorphism is not None else pre_isomorphism

        if post_isomorphism is not None:
            if not isinstance(post_isomorphism, (WeierstrassIsomorphism, type(None))):
                raise ValueError("invalid parameter: post-isomorphism must be a WeierstrassIsomorphism")
            if self._codomain != post_isomorphism.domain():
                raise ValueError("invalid parameter: post-isomorphism must have domain curve equal to this isogenies' codomain")
            final_post_isomorphism = post_isomorphism * final_post_isomorphism if final_post_isomorphism is not None else post_isomorphism

        result = EllipticCurveIsogeny.__new__(EllipticCurveIsogeny)
        # we're bypassing calling __init__, make sure we do an equivalent amount of work

        for attr in ["_degree", "__algorithm", "__check", "__prei_ratl_maps",
                     "__posti_ratl_maps", "__base_field", "__poly_ring",
                     "__mpoly_ring", "__kernel_list", "__kernel_polynomial",
                     "__inner_kernel_polynomial", "__kernel_mod_sign", "__v",
                     "__w", "__psi", "__phi", "__omega", "__xfield",
                     "__xyfield"]:
            if attr.startswith("__"):
                attr = "_EllipticCurveIsogeny" + attr
            setattr(result, attr, getattr(self, attr))

        if final_pre_isomorphism is not None:
            result.__set_pre_isomorphism(
                    self._domain if final_pre_isomorphism is None else final_pre_isomorphism.domain(),
                    final_pre_isomorphism)
        else:
            result._domain = self._domain

        if final_post_isomorphism is not None:
            result.__set_post_isomorphism(
                    self._codomain if final_post_isomorphism is None else final_post_isomorphism.codomain(),
                    final_post_isomorphism)
        else:
            result._codomain = self._codomain

        result.__perform_inheritance_housekeeping()
        return result

    def _eval(self, P):
        r"""
        Less strict evaluation method for internal use.

        In particular, this can be used to evaluate ``self`` at a
        point defined over an extension field.

        INPUT:

        - ``P`` -- a sequence of 3 coordinates defining a point on ``self``

        OUTPUT: the result of evaluating ``self`` at the given point

        EXAMPLES::

            sage: E = EllipticCurve([1,0]); E
            Elliptic Curve defined by y^2 = x^3 + x over Rational Field
            sage: phi = E.isogeny(E(0,0))
            sage: P = E.change_ring(QQbar).lift_x(QQbar.random_element())               # needs sage.rings.number_field
            sage: phi._eval(P).curve()                                                  # needs sage.rings.number_field
            Elliptic Curve defined by y^2 = x^3 + (-4)*x over Algebraic Field

        ::

            sage: E = EllipticCurve(j=Mod(0,419))
            sage: K = next(filter(bool, E(0).division_points(5)))
            sage: psi = E.isogeny(K)
            sage: Ps = E.change_ring(GF(419**2))(0).division_points(5)                  # needs sage.rings.number_field
            sage: {psi._eval(P).curve() for P in Ps}                                    # needs sage.rings.number_field
            {Elliptic Curve defined by y^2 = x^3 + 140*x + 214 over Finite Field in z2 of size 419^2}
        """
        if self._domain.defining_polynomial()(*P):
            raise ValueError(f"{P} not on {self._domain}")

        k = Sequence(P).universe()

        if not P:
            return self._codomain(0).change_ring(k)

        Q = P.xy()

        if self.__pre_isomorphism is not None:
            Q = baseWI.__call__(self.__pre_isomorphism, Q)

        if self.__algorithm == "velu":
            compute = self.__compute_via_velu
        elif self.__algorithm == "kohel":
            compute = self.__compute_via_kohel
        else:
            raise NotImplementedError

        try:
            Q = compute(*Q)
        except ZeroDivisionError:
            Q = (0, 1, 0)

        if self.__post_isomorphism is not None:
            Q = baseWI.__call__(self.__post_isomorphism, Q)

        return self._codomain.base_extend(k).point(Q)

    def _call_(self, P):
        r"""
        Implement evaluation of elliptic-curve isogenies using the
        function-call syntax.

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
            sage: tau = E.isogeny([Q], E)
            sage: tau(Q)
            (0 : 1 : 0)

        TESTS:

        Tests for :issue:`10888`::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<th> = NumberField(x^2 + 3)
            sage: E = EllipticCurve(K, [7,0])
            sage: phi = E.isogeny(E(0,0))
            sage: P = E(-3,4*th)
            sage: phi(P)
            (-16/3 : 8/9*th : 1)
            sage: Q = phi(P)
            sage: phihat = phi.dual()
            sage: phihat(Q)
            (-1/48 : 127/576*th : 1)

        Call a composed isogeny (added for :issue:`16238`)::

            sage: E = EllipticCurve(j=GF(7)(0))
            sage: phi = E.isogeny([E(0), E((0,1)), E((0,-1))])
            sage: phi(E.points()[0])
            (0 : 1 : 0)
            sage: phi2 = phi * phi
            sage: phi2(E.points()[0])
            (0 : 1 : 0)

        Coercion works fine with :meth:`_call_` (added for :issue:`16238`)::

            sage: # needs sage.rings.number_field
            sage: K.<th> = NumberField(x^2 + 3)
            sage: E = EllipticCurve(K, [7,0])
            sage: E2 = EllipticCurve(K, [5,0])
            sage: phi = E.isogeny(E(0))
            sage: phi(E2(0))
            (0 : 1 : 0)
            sage: E2(20,90)
            (20 : 90 : 1)
            sage: phi(E2(20,90))
            Traceback (most recent call last):
            ...
            TypeError: (20 : 90 : 1) fails to convert into the map's domain
            Elliptic Curve defined by y^2 = x^3 + 7*x over
            Number Field in th with defining polynomial x^2 + 3,
            but a `pushforward` method is not properly implemented

        Check that copying the order over works::

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve(GF(431), [1,0])
            sage: P, = E.gens()
            sage: Q = 2^99*P; Q.order()
            27
            sage: phi = E.isogeny(3^99*P)
            sage: phi(Q)._order
            27

        Test for :issue:`35983`::

            sage: E = EllipticCurve([1,0,0,-1,0])
            sage: P = E([1,0])
            sage: P.order()
            +Infinity
            sage: phi = E.isogenies_prime_degree(2)[0]
            sage: Q = phi(P); Q
            (0 : 1 : 1)
            sage: Q.order()
            +Infinity
        """
        if P.is_zero():
            return self._codomain(0)

        xP, yP = P.xy()

        # if there is a pre-isomorphism, apply it
        if self.__pre_isomorphism is not None:
            yP = self.__prei_ratl_maps[1](xP, yP)
            xP = self.__prei_ratl_maps[0](xP)

        if self.__algorithm == 'velu':
            outP = self.__compute_via_velu_numeric(xP, yP)
        elif self.__algorithm == 'kohel':
            outP = self.__compute_via_kohel_numeric(xP, yP)
        else:
            raise NotImplementedError

        # the intermediate functions return the empty tuple ()
        # if the input point is in the kernel
        if outP == ():
            return self._codomain(0)

        xP, yP = outP

        # if there is a post-isomorphism, apply it
        if self.__post_isomorphism is not None:
            yP = self.__posti_ratl_maps[1](xP, yP)
            xP = self.__posti_ratl_maps[0](xP)

        Q = self._codomain(xP, yP)
        if hasattr(P, '_order'):
            if P.has_infinite_order() or P._order.gcd(self._degree).is_one():
                Q._order = P._order
            # TODO: For non-coprime degree, the order of the point
            # may get reduced by a divisor of the degree when passing
            # through the isogeny. We could run something along the
            # lines of order_from_multiple() to determine the new
            # order, but this probably shouldn't happen by default
            # as it'll be detrimental to performance in some cases.
        return Q

    def __getitem__(self, i):
        r"""
        Return one of the rational-map components.

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
            sage: phi = EllipticCurveIsogeny(E, E((0,0)))
            sage: phi[0]
            (x^2 + 3)/x
            sage: phi[1]
            (x^2*y - 3*y)/x^2
        """
        return self.rational_maps()[i]

    def __iter__(self):
        r"""
        Return an iterator through the rational-map components.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,2,0,1,-1])
            sage: phi = EllipticCurveIsogeny(E, [1])
            sage: for c in phi: print(c)
            x
            y

            sage: E = EllipticCurve(GF(17), [0,0,0,3,0])
            sage: phi = EllipticCurveIsogeny(E, E((0,0)))
            sage: for c in phi: print(c)
            (x^2 + 3)/x
            (x^2*y - 3*y)/x^2
        """
        return iter(self.rational_maps())

    def __neg__(self):
        r"""
        Return a copy of the isogeny that has been negated.

        This implements the unary `-` operator.

        EXAMPLES:

        The following examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(17)(0))
            sage: phi = EllipticCurveIsogeny(E, E((-1,0)) )
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

            sage: E = EllipticCurve(GF(23), [0,0,0,1,0])
            sage: f = E.torsion_polynomial(3)/3
            sage: phi = EllipticCurveIsogeny(E, f, E)
            sage: phi.rational_maps() == E.multiplication_by_m(3)
            False
            sage: negphi = -phi
            sage: negphi.rational_maps() == E.multiplication_by_m(3)
            True

            sage: E = EllipticCurve(GF(17), [-2, 3, -5, 7, -11])
            sage: R.<x> = GF(17)[]
            sage: f = x+6
            sage: phi = EllipticCurveIsogeny(E, f); phi
            Isogeny of degree 2
             from Elliptic Curve defined by y^2 + 15*x*y + 12*y = x^3 + 3*x^2 + 7*x + 6 over Finite Field of size 17
               to Elliptic Curve defined by y^2 + 15*x*y + 12*y = x^3 + 3*x^2 + 4*x + 8 over Finite Field of size 17
            sage: phi.rational_maps()
            ((x^2 + 6*x + 4)/(x + 6), (x^2*y - 5*x*y + 8*x - 2*y)/(x^2 - 5*x + 2))
            sage: negphi = -phi
            sage: negphi
            Isogeny of degree 2
             from Elliptic Curve defined by y^2 + 15*x*y + 12*y = x^3 + 3*x^2 + 7*x + 6 over Finite Field of size 17
               to Elliptic Curve defined by y^2 + 15*x*y + 12*y = x^3 + 3*x^2 + 4*x + 8 over Finite Field of size 17
            sage: negphi.rational_maps()
            ((x^2 + 6*x + 4)/(x + 6),
             (2*x^3 - x^2*y - 5*x^2 + 5*x*y - 4*x + 2*y + 7)/(x^2 - 5*x + 2))

            sage: E = EllipticCurve('11a1')
            sage: R.<x> = QQ[]
            sage: f = x^2 - 21*x + 80
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: (xmap1, ymap1) = phi.rational_maps()
            sage: negphi = -phi
            sage: (xmap2, ymap2) = negphi.rational_maps()
            sage: xmap1 == xmap2
            True
            sage: ymap1 == -ymap2 - E.a1()*xmap2 - E.a3()
            True

            sage: # needs sage.rings.number_field
            sage: K.<a> = NumberField(x^2 + 1)
            sage: E = EllipticCurve(K, [0,0,0,1,0])
            sage: R.<x> = K[]
            sage: phi = EllipticCurveIsogeny(E, x - a)
            sage: phi.rational_maps()
            ((x^2 + (-a)*x - 2)/(x + (-a)), (x^2*y + (-2*a)*x*y + y)/(x^2 + (-2*a)*x - 1))
            sage: negphi = -phi
            sage: negphi.rational_maps()
            ((x^2 + (-a)*x - 2)/(x + (-a)), (-x^2*y + (2*a)*x*y - y)/(x^2 + (-2*a)*x - 1))
        """
        return self._compose_with_isomorphism(post_isomorphism=negation_morphism(self._codomain))

    #
    # Sage Special Functions
    #

    def _repr_(self):
        r"""
        Return basic information about the isogeny as a string.

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
        return f'Isogeny of degree {self._degree} from {self._domain} to {self._codomain}'

    def _latex_(self):
        r"""
        Return the rational maps of the isogeny as a LaTeX string.

        This function returns a latex string representing the isogeny
        ``self`` as the `x` and `y` coordinate rational functions.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,0,0,1,-1])
            sage: phi = EllipticCurveIsogeny(E, E(0))
            sage: phi._latex_()
            '\\left( x , y \\right)'

            sage: E = EllipticCurve(GF(17), [0,0,0,1,-1])
            sage: R.<X> = GF(17)[]
            sage: phi = EllipticCurveIsogeny(E, X + 11)
            sage: phi._latex_()
            '\\left( \\frac{x^{2} + 11 x + 7}{x + 11} , \\frac{x^{2} y + 5 x y + 12 y}{x^{2} + 5 x + 2} \\right)'
        """
        fx,fy = self.rational_maps()
        return fr'\left( {fx._latex_()} , {fy._latex_()} \right)'

    ###########################
    # Private Common Functions
    ###########################

    def __perform_inheritance_housekeeping(self):
        r"""
        Internal helper function, sets values on the super classes of
        this class.

        This should only be called during :meth:`__init__`,
        and should not be called more than once.

        EXAMPLES:

        The following examples will implicitly exercise this
        function::

            sage: E = EllipticCurve(GF(43), [2,3,5,7,11])
            sage: R.<x> = GF(43)[]; f = x + 42
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E2 = phi.codomain()
            sage: post_isom = WeierstrassIsomorphism(E2, (41, 37, 31, 29))
            sage: phi = phi._compose_with_isomorphism(post_isomorphism=post_isom)
            sage: E1pr = WeierstrassIsomorphism(E, (-1, 2, -3, 4)).codomain()
            sage: pre_isom = E1pr.isomorphism_to(E)
            sage: phi = phi._compose_with_isomorphism(pre_isomorphism=pre_isom)
        """
        EllipticCurveHom.__init__(self, self._domain, self._codomain)

    def __init_algebraic_structs(self, E):
        r"""
        An internal function for EllipticCurveIsogeny objects that
        sets up the member variables necessary for algebra.

        EXAMPLES::

            sage: E = EllipticCurve(j=GF(17)(0))
            sage: phi = EllipticCurveIsogeny(E, E((-1,0)))

        The constructor calls this function itself, so the fields it
        sets are already defined::

            sage: phi._domain
            Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 17
            sage: phi._EllipticCurveIsogeny__base_field
            Finite Field of size 17
            sage: phi._EllipticCurveIsogeny__poly_ring
            Univariate Polynomial Ring in x over Finite Field of size 17
            sage: phi._EllipticCurveIsogeny__mpoly_ring
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Finite Field of size 17

        Now, calling the initialization function does nothing more::

            sage: phi._EllipticCurveIsogeny__init_algebraic_structs(E)
            sage: phi._domain
            Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 17
            sage: phi._EllipticCurveIsogeny__base_field
            Finite Field of size 17
            sage: phi._EllipticCurveIsogeny__poly_ring
            Univariate Polynomial Ring in x over Finite Field of size 17
            sage: phi._EllipticCurveIsogeny__mpoly_ring
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Finite Field of size 17

            sage: E = EllipticCurve(QQ, [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E, E((0,0)))
            sage: phi._EllipticCurveIsogeny__init_algebraic_structs(E)
            sage: phi._domain
            Elliptic Curve defined by y^2 = x^3 + x over Rational Field
            sage: phi._EllipticCurveIsogeny__base_field
            Rational Field
            sage: phi._EllipticCurveIsogeny__poly_ring
            Univariate Polynomial Ring in x over Rational Field
            sage: phi._EllipticCurveIsogeny__mpoly_ring
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: F = GF(19); R.<x> = F[]
            sage: E = EllipticCurve(j=GF(19)(0))
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: phi._EllipticCurveIsogeny__init_algebraic_structs(E)
            sage: phi._domain
            Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 19
            sage: phi._EllipticCurveIsogeny__base_field
            Finite Field of size 19
            sage: phi._EllipticCurveIsogeny__poly_ring
            Univariate Polynomial Ring in x over Finite Field of size 19
            sage: phi._EllipticCurveIsogeny__mpoly_ring
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Finite Field of size 19
        """
        self._domain = E
        self.__base_field = E.base_ring()
        self.__poly_ring = PolynomialRing(self.__base_field, ['x'])
        self.__mpoly_ring = PolynomialRing(self.__poly_ring, ['y'])
        # The fraction fields are implicitly part of the public API, being the parents
        # of the rational maps.
        self.__xfield = FractionField(self.__poly_ring)
        self.__xyfield = FractionField(PolynomialRing(self.__base_field, ['x', 'y']))

    def __compute_codomain(self):
        r"""
        Private function that computes and sets the isogeny codomain.

        EXAMPLES:

        These examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(7)(1728))
            sage: phi = EllipticCurveIsogeny(E, E((0,0)))
            sage: phi.codomain()
            Elliptic Curve defined by y^2 = x^3 + 3*x over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__compute_codomain()

            sage: R.<x> = GF(7)[]
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: phi.codomain()
            Elliptic Curve defined by y^2 = x^3 + 3*x over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__compute_codomain()
        """
        if self.__algorithm == 'velu':
            self._codomain = self.__compute_codomain_via_velu()
        elif self.__algorithm == 'kohel':
            self._codomain = self.__compute_codomain_via_kohel()
        else:
            raise NotImplementedError

    def __initialize_rational_maps(self, precomputed_maps=None):
        r"""
        Private function that computes and initializes the rational
        maps.

        INPUT:

        - ``precomputed_maps`` -- (default: ``None``) tuple `(X,Y)`
          of rational functions in `x,y`

        EXAMPLES:

        The following examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(7)(1728))
            sage: phi = EllipticCurveIsogeny(E, E((0,0)))
            sage: phi.rational_maps()  # implicit doctest
            ((x^2 + 1)/x, (x^2*y - y)/x^2)

            sage: R.<x> = GF(7)[]
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: phi.rational_maps()  # implicit doctest
            ((x^2 + 1)/x, (x^2*y - y)/x^2)

            sage: E = EllipticCurve([1,2,3,4,5])
            sage: Eshort = E.short_weierstrass_model()
            sage: phi = E.isogeny(E(0), Eshort)
            sage: phiX, phiY = phi.rational_maps()  # implicit doctest
            sage: phiX(1,2), phiY(1,2)
            (63, 864)
        """
        if self.__ratl_maps is not None:
            return

        if precomputed_maps is None:
            if self.__algorithm == 'velu':
                X_map, Y_map = self.__initialize_rational_maps_via_velu()
            elif self.__algorithm == 'kohel':
                X_map, Y_map = self.__initialize_rational_maps_via_kohel()
            else:
                raise NotImplementedError

        else:
            X_map, Y_map = precomputed_maps
            # cannot coerce directly in xfield for some reason
            X_map = self.__poly_ring(X_map.numerator()) \
                    / self.__poly_ring(X_map.denominator())

        if self.__prei_ratl_maps is not None:
            prei_X_map, prei_Y_map = self.__prei_ratl_maps
            X_map = X_map(prei_X_map)
            Y_map = Y_map([prei_X_map, prei_Y_map])

        if self.__posti_ratl_maps is not None:
            posti_X_map, posti_Y_map = self.__posti_ratl_maps
            # Do not reverse the order here!
            Y_map = posti_Y_map([X_map, Y_map])
            X_map = posti_X_map(X_map)

        self.__ratl_maps = self.__xfield(X_map), self.__xyfield(Y_map)

    def __init_kernel_polynomial(self):
        r"""
        Private function that initializes the kernel polynomial (if
        the algorithm does not take it as a parameter).

        EXAMPLES:

        The following examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(7)(1728))
            sage: phi = EllipticCurveIsogeny(E, E((0,0)))
            sage: phi.kernel_polynomial()  # implicit doctest
            x
        """
        if self.__kernel_polynomial is None:
            if self.__algorithm == 'velu':
                self.__init_kernel_polynomial_velu()
            else:
                assert False, "the kernel polynomial should already be defined!"

    def __set_pre_isomorphism(self, domain, isomorphism):
        r"""
        Private function to set the pre-isomorphism and domain (and
        keep track of the domain of the isogeny).

        This should only be called during initialization,
        :meth:`__perform_inheritance_housekeeping` should be called after.

        TESTS::

            sage: E = EllipticCurve(GF(31), [1,1,0,1,-1])
            sage: R.<x> = GF(31)[]
            sage: f = x^3 + 9*x^2 + x + 30
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: Epr = E.short_weierstrass_model()
            sage: isom = Epr.isomorphism_to(E)
            sage: phi = phi._compose_with_isomorphism(isom)
            sage: phi.rational_maps()
            ((-6*x^4 - 3*x^3 + 12*x^2 + 10*x - 1)/(x^3 + x - 12),
             (3*x^7 + x^6*y - 14*x^6 - 3*x^5 + 5*x^4*y + 7*x^4 + 8*x^3*y - 8*x^3 - 5*x^2*y + 5*x^2 - 14*x*y + 14*x - 6*y - 6)/(x^6 + 2*x^4 + 7*x^3 + x^2 + 7*x - 11))
            sage: phi(Epr((0,22)))
            (13 : 21 : 1)
            sage: phi(Epr((3,7)))
            (14 : 17 : 1)

            sage: E = EllipticCurve(GF(29), [0,0,0,1,0])
            sage: R.<x> = GF(29)[]
            sage: f = x^2 + 5
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi
            Isogeny of degree 5
             from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 29
               to Elliptic Curve defined by y^2 = x^3 + 20*x over Finite Field of size 29
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: inv_isom = WeierstrassIsomorphism(E, (1,-2,5,10))
            sage: Epr = inv_isom.codomain()
            sage: isom = Epr.isomorphism_to(E)
            sage: phi = phi._compose_with_isomorphism(isom)
            sage: phi
            Isogeny of degree 5
             from Elliptic Curve defined by y^2 + 10*x*y + 20*y = x^3 + 27*x^2 + 6 over Finite Field of size 29
               to Elliptic Curve defined by y^2 = x^3 + 20*x over Finite Field of size 29
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
            Isogeny of degree 5
             from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
               to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: Epr = E.short_weierstrass_model()
            sage: isom = Epr.isomorphism_to(E)
            sage: phi = phi._compose_with_isomorphism(isom)
            sage: phi
            Isogeny of degree 5
             from Elliptic Curve defined by y^2 = x^3 - 13392*x - 1080432 over Rational Field
               to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
            sage: phi(Epr((168,1188)))
            (0 : 1 : 0)
        """
        self._domain = domain
        self.__pre_isomorphism = isomorphism

        # calculate the isomorphism as a rational map.

        u, r, s, t = (self.__base_field(c) for c in isomorphism.tuple())
        uinv = 1/u
        uinv2 = uinv**2
        uinv3 = uinv*uinv2

        x = self.__poly_ring.gen()
        y = self.__xyfield.gen(1) # not mpoly_ring.gen(1) else we end
                                  # up in K(x)[y] and trouble ensues

        self.__prei_ratl_maps = (x - r) * uinv2, (y - s*(x-r) - t) * uinv3

        if self.__kernel_polynomial is not None:
            ker_poly = self.__kernel_polynomial
            ker_poly = ker_poly(self.__prei_ratl_maps[0])
            self.__kernel_polynomial = ker_poly.monic()

    def __set_post_isomorphism(self, codomain, isomorphism):
        r"""
        Private function to set the post-isomorphism and codomain (and
        keep track of the codomain of the isogeny).

        This should only be called during initialization,
        :meth:`__perform_inheritance_housekeeping` should be called after.

        TESTS::

            sage: E = EllipticCurve(j=GF(31)(0))
            sage: R.<x> = GF(31)[]
            sage: phi = EllipticCurveIsogeny(E, x+18)
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: phi = phi._compose_with_isomorphism(post_isomorphism=WeierstrassIsomorphism(phi.codomain(), (6,8,10,12)))
            sage: phi
            Isogeny of degree 3
             from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 31
               to Elliptic Curve defined by y^2 + 24*x*y + 7*y = x^3 + 22*x^2 + 16*x + 20 over Finite Field of size 31

            sage: E = EllipticCurve(j=GF(47)(0))
            sage: f = E.torsion_polynomial(3)/3
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: E2 = phi.codomain()
            sage: post_isom = E2.isomorphism_to(E)
            sage: phi = phi._compose_with_isomorphism(post_isomorphism=post_isom)
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
            Isogeny of degree 4
             from Elliptic Curve defined by y^2 = x^3 + x over Number Field in a with defining polynomial x^2 + 2
               to Elliptic Curve defined by y^2 = x^3 + (-44)*x + 112 over Number Field in a with defining polynomial x^2 + 2
        """
        self._codomain = codomain
        self.__post_isomorphism = isomorphism

        # calculate the isomorphism as a rational map.

        u, r, s, t = (self.__base_field(c) for c in isomorphism.tuple())
        uinv = 1/u
        uinv2 = uinv**2
        uinv3 = uinv*uinv2

        x = self.__poly_ring.gen()
        y = self.__xyfield.gen(1)

        self.__posti_ratl_maps = (x - r) * uinv2, (y - s*(x-r) - t) * uinv3

    def __setup_post_isomorphism(self, codomain, model):
        r"""
        Private function to set up the post-isomorphism given the
        codomain.

        EXAMPLES:

        The following examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(7)(1728))
            sage: E2 = EllipticCurve(GF(7), [0,0,0,5,0])
            sage: phi = EllipticCurveIsogeny(E, E((0,0)), E2); phi
            Isogeny of degree 2
             from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7
               to Elliptic Curve defined by y^2 = x^3 + 5*x over Finite Field of size 7
            sage: E3 = EllipticCurve(GF(7), [0,0,0,6,0])
            sage: phi._EllipticCurveIsogeny__setup_post_isomorphism(E3, None)
            sage: phi
            Isogeny of degree 2
             from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7
               to Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7

            sage: EllipticCurveIsogeny(E, E(0,0), model='montgomery')
            Isogeny of degree 2
             from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7
               to Elliptic Curve defined by y^2 = x^3 + x^2 + x over Finite Field of size 7

            sage: R.<x> = QQ[]
            sage: E = EllipticCurve(j=1728)
            sage: f = x^3 - x
            sage: phi = EllipticCurveIsogeny(E, f, model='minimal'); phi
            Isogeny of degree 4
             from Elliptic Curve defined by y^2 = x^3 - x over Rational Field
               to Elliptic Curve defined by y^2 = x^3 - x over Rational Field

            sage: phi = EllipticCurveIsogeny(E, f, model=None)
            sage: phi._EllipticCurveIsogeny__setup_post_isomorphism(None, 'minimal')
            sage: phi
            Isogeny of degree 4
             from Elliptic Curve defined by y^2 = x^3 - x over Rational Field
               to Elliptic Curve defined by y^2 = x^3 - x over Rational Field
        """
        if model is codomain is None:
            return

        oldE2 = self._codomain

        if model is not None:
            if codomain is not None:
                raise ValueError("cannot specify a codomain curve and model name simultaneously")

            from sage.schemes.elliptic_curves.ell_field import compute_model
            codomain = compute_model(oldE2, model)

        else:  # codomain is not None
            if not isinstance(codomain, EllipticCurve_generic):
                raise ValueError("given codomain is not an elliptic curve")

            if not oldE2.is_isomorphic(codomain):
                raise ValueError("given codomain is not isomorphic to the computed codomain")

        post_isom = oldE2.isomorphism_to(codomain)
        self.__set_post_isomorphism(codomain, post_isom)

    ###########################
    # Velu's Formula Functions
    ###########################

    #
    # Setup function for Velu's formula
    #

    def __init_from_kernel_gens(self, kernel_gens):
        r"""
        Private function that initializes the isogeny from a list of
        points which generate the kernel (For Vélu's formulas.)

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: phi = EllipticCurveIsogeny(E, E((0,0))); phi
            Isogeny of degree 2
             from Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7
               to Elliptic Curve defined by y^2 = x^3 + 4*x over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__init_from_kernel_gens([E(0), E((0,0))])

        The following example demonstrates the necessity of avoiding any calls
        to P.order(), since such calls involve factoring the group order which
        could take a long time. ::

            sage: # needs sage.rings.finite_rings
            sage: p = 12 * next_prime(2^180) * next_prime(2^194) - 1
            sage: F = FiniteField(p, proof=False)
            sage: E = EllipticCurve([F(1), F(0)])
            sage: P = E(0).division_points(3)[1]
            sage: EllipticCurveIsogeny(E, P)
            Isogeny of degree 3
             from Elliptic Curve defined by y^2 = x^3 + x
                  over Finite Field of size 461742260113997803268895001173557974278278194575766957660028841364655249961609425998827452443620996655395008156411
               to Elliptic Curve defined by y^2 = x^3 + 80816485163488178037199320944019099858815874115367810482828676054000067654558381377552245721755005198633191074893*x + 301497584865165444049833326660609767433467459033532853758006118022998267706948164646650354324860226263546558337993
                  over Finite Field of size 461742260113997803268895001173557974278278194575766957660028841364655249961609425998827452443620996655395008156411
        """
        if self.__check:
            for P in kernel_gens:
                if not P.has_finite_order():
                    raise ValueError("given kernel contains point of infinite order")

        self.__kernel_mod_sign = {}
        self.__v = self.__w = 0

        # Fast path: The kernel is given by a single generating point.
        if len(kernel_gens) == 1 and kernel_gens[0]:
            self.__init_from_kernel_point(kernel_gens[0])
            return

        # Compute a list of points in the subgroup generated by the
        # points in kernel_gens.  This is very naive: when finite
        # subgroups are implemented better, this could be simplified,
        # but it won't speed things up too much.

        def all_multiples(itr, terminal):
            mult_list = [terminal]
            R = terminal + itr
            while R != terminal:
                mult_list.append(R)
                R += itr
            return mult_list

        kernel_set = {self._domain(0)}
        for P in kernel_gens:
            kernel_set.update(R for Q in tuple(kernel_set)
                                for R in all_multiples(P,Q))

        self._degree = Integer(len(kernel_set))
        self.__kernel_list = list(kernel_set)
        self.__init_from_kernel_list()

    #
    # Precompute the values in Velu's Formula.
    #
    def __update_kernel_data(self, xQ, yQ):
        r"""
        Internal helper function to update some data coming from the
        kernel points of this isogeny when using Vélu's formulas.

        TESTS:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, [P,P]); phi  # implicit doctest
            Isogeny of degree 4
             from Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7
               to Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 7
        """
        a1, a2, a3, a4, _ = self._domain.a_invariants()

        gxQ = (3*xQ + 2*a2)*xQ + a4 - a1*yQ
        gyQ = -2*yQ - a1*xQ - a3

        uQ = gyQ**2

        if 2*yQ == -a1*xQ - a3: # Q is 2-torsion
            vQ = gxQ
        else:                   # Q is not 2-torsion
            vQ = 2*gxQ - a1*gyQ

        self.__kernel_mod_sign[xQ] = yQ, gxQ, gyQ, vQ, uQ

        self.__v += vQ
        self.__w += uQ + xQ*vQ

    def __init_from_kernel_point(self, ker):
        r"""
        Private function with functionality equivalent to
        :meth:`__init_from_kernel_list`, but optimized for when
        the kernel is given by a single point.

        TESTS:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P); phi  # implicit doctest
            Isogeny of degree 4
             from Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7
               to Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 7

        We check that the result is the same as for :meth:`__init_from_kernel_list`::

            sage: psi = EllipticCurveIsogeny(E, [P, P]); psi
            Isogeny of degree 4
             from Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7
               to Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 7
            sage: phi == psi
            True
        """
        self._degree = Integer(1)

        Q, prevQ = ker, self._domain(0)

        while Q and Q != -prevQ:
            self.__update_kernel_data(*Q.xy())

            if Q == -Q:
                self._degree += 1
                break

            prevQ = Q
            Q += ker
            self._degree += 2

    def __init_from_kernel_list(self):
        r"""
        Private function that sorts the list of points in the kernel
        (For Vélu's formulas). Sorts out the 2-torsion points, and
        puts them in a dictionary.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, [P,P]); phi  # implicit doctest
            Isogeny of degree 4
             from Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7
               to Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 7
        """
        for Q in self.__kernel_list:

            if Q.is_zero():
                continue

            xQ, yQ = Q.xy()

            if xQ in self.__kernel_mod_sign:
                continue

            self.__update_kernel_data(xQ, yQ)

    #
    # Velu's formula computing the codomain curve
    #
    def __compute_codomain_via_velu(self):
        r"""
        Private function that computes the codomain via Vélu's
        formulas.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P)
            sage: phi.codomain()
            Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__compute_codomain_via_velu()
            Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 7
        """
        return compute_codomain_formula(self._domain, self.__v, self.__w)

    @staticmethod
    def __velu_sum_helper(xQ, Qvalues, a1, a3, x, y):
        r"""
        Private function for Vélu's formulas, helper function to help
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
            sage: phi = EllipticCurveIsogeny(E, E((0,0)) )
            sage: Qvals = phi._EllipticCurveIsogeny__kernel_mod_sign[0]
            sage: phi._EllipticCurveIsogeny__velu_sum_helper(0, Qvals, 0, 0, F(5), F(5))
            (3, 3)
            sage: R.<x,y> = GF(7)[]
            sage: phi._EllipticCurveIsogeny__velu_sum_helper(0, Qvals, 0, 0, x, y)
            (1/x, y/x^2)
        """
        yQ, gxQ, gyQ, vQ, uQ = Qvalues

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
        F = FractionField(y.parent())
        tY = tY0*F(inv_t1_3) + (tY1 + tY2)*F(inv_t1_2)

        return tX, tY

    def __compute_via_velu_numeric(self, xP, yP):
        r"""
        Private function that sorts the list of points in the kernel
        (for Vélu's formulas). Sorts out the 2-torsion points, and
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
        if xP in self.__kernel_mod_sign:
            return ()

        return self.__compute_via_velu(xP,yP)

    def __compute_via_velu(self, xP, yP):
        r"""
        Private function for Vélu's formulas, to perform the summation.

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

        TESTS:

        Check for :issue:`33214`::

            sage: # needs sage.rings.finite_rings
            sage: z2 = GF(71^2).gen()
            sage: E = EllipticCurve(GF(71^2), [5,5])
            sage: phi = E.isogeny(E.lift_x(0))
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: pre = WeierstrassIsomorphism(None, (z2,7,8,9), E)
            sage: phi = phi * pre
            sage: P = phi.domain()(1, 46*z2+49)
            sage: phi(P)  # indirect doctest
            (33 : 61*z2 + 10 : 1)

        The rational maps are also computed via this code path; check
        that they are plausible (this failed prior to :issue:`33214`)::

            sage: # needs sage.rings.finite_rings
            sage: fx,fy = phi.rational_maps()  # indirect doctest
            sage: R.<x,y> = GF(71^2)[]
            sage: E0, E2 = phi.domain(), phi.codomain()
            sage: eqs = [EE.defining_polynomial()(x,y,1) for EE in (E0,E2)]
            sage: eqs[1](fx,fy).numerator() % eqs[0]
            0
        """
        if self.__pre_isomorphism is None:
            E = self._domain
        else:
            E = self.__pre_isomorphism.codomain()

        a1 = E.a1()
        a3 = E.a3()

        X = xP
        Y = yP

        for xQ, Qvalues in self.__kernel_mod_sign.items():
            tX, tY = self.__velu_sum_helper(xQ, Qvalues, a1, a3, xP, yP)
            X += tX
            Y -= tY

        return X, Y

    def __initialize_rational_maps_via_velu(self):
        r"""
        Private function for Vélu's formulas, helper function to
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
        y = self.__xyfield.gen(1)
        return self.__compute_via_velu(x,y)

    def __init_kernel_polynomial_velu(self):
        r"""
        Private function for Vélu's formulas, helper function to
        initialize the rational maps.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P)
            sage: phi.kernel_polynomial()  # implicit doctest
            x^2 + 2*x + 4
        """
        poly_ring, x = self.__poly_ring.objgen()

        if self.__pre_isomorphism is not None:
            pre_isom = self.__pre_isomorphism
            invX = pre_isom.u**2 * x + pre_isom.r
        else:
            invX = x

        psi = poly_ring.one()
        for xQ in self.__kernel_mod_sign.keys():
            psi *= x - invX(xQ)

        self.__kernel_polynomial = psi

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
            Isogeny of degree 2
             from Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7
               to Elliptic Curve defined by y^2 = x^3 + 4*x over Finite Field of size 7

            sage: phi._EllipticCurveIsogeny__init_from_kernel_polynomial(x)

            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3); phi
            Isogeny of degree 3
             from Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 1 over Finite Field of size 7
               to Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 4*x + 2 over Finite Field of size 7

            sage: phi._EllipticCurveIsogeny__init_from_kernel_polynomial(x+6)
        """
        poly_ring = self.__poly_ring
        E = self._domain

        # Convert to a univariate polynomial, even if it had a
        # bivariate parent, or was given as a list:
        psi = poly_ring(kernel_polynomial)

        self.__kernel_polynomial = psi

        if psi.leading_coefficient() != 1:
            raise ValueError("given kernel polynomial is not monic")

        #
        # Determine if kernel polynomial is entirely 2-torsion
        #
        psi_G = two_torsion_part(E, psi).monic()

        if psi_G.degree() != 0: # even degree case

            psi_quo = psi//psi_G

            if psi_quo.degree() != 0:
                raise NotImplementedError("Kohel's algorithm currently only supports cyclic isogenies (except for [2])")

            phi, omega, v, w, _, d = self.__init_even_kernel_polynomial(E, psi_G)

        else: # odd degree case

            phi, omega, v, w, _, d = self.__init_odd_kernel_polynomial(E, psi)

        #
        # Set up the necessary instance variables
        #

        self.__kernel_polynomial = psi
        self.__inner_kernel_polynomial = psi

        self._degree = Integer(d)  # degree of the isogeny

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

        A tuple (``phi``, ``omega``, ``v``, ``w``, ``n``, ``d``) where:

        - ``phi`` is a univariate polynomial, the numerator of the
          `X`-coordinate of the isogeny;

        - ``omega`` is a bivariate polynomial, the numerator of the
          `Y`-coordinate of the isogeny;

        - ``v``, ``w`` are the Vélu parameters of the isogeny;

        - ``n`` is the degree of ``psi``;

        - ``d`` is the degree of the isogeny.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [-1,0])
            sage: phi = EllipticCurveIsogeny(E, x); phi
            Isogeny of degree 2
             from Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7
               to Elliptic Curve defined by y^2 = x^3 + 4*x over Finite Field of size 7

            sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import two_torsion_part
            sage: psig = two_torsion_part(E,x)
            sage: phi._EllipticCurveIsogeny__init_even_kernel_polynomial(E,psig)
            (x^3 + 6*x, (x^3 + x)*y, 6, 0, 1, 2)

            sage: F = GF(2^4, 'alpha'); R.<x> = F[]                                     # needs sage.rings.finite_rings
            sage: E = EllipticCurve(F, [1,1,0,1,0])                                     # needs sage.rings.finite_rings
            sage: phi = EllipticCurveIsogeny(E, x); phi                                 # needs sage.rings.finite_rings
            Isogeny of degree 2
             from Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + x over Finite Field in alpha of size 2^4
               to Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + 1 over Finite Field in alpha of size 2^4

            sage: psig = two_torsion_part(E,x)                                          # needs sage.rings.finite_rings
            sage: phi._EllipticCurveIsogeny__init_even_kernel_polynomial(E,psig)        # needs sage.rings.finite_rings
            (x^3 + x, (x^3 + x)*y + x^2, 1, 0, 1, 2)

            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: R.<x> = GF(7)[]
            sage: f = x^3 + 6*x^2 + 1
            sage: phi = EllipticCurveIsogeny(E, f); phi
            Isogeny of degree 4
             from Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 1 over Finite Field of size 7
               to Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 2*x + 5 over Finite Field of size 7
            sage: psig = two_torsion_part(E,f)
            sage: psig = two_torsion_part(E,f)
            sage: phi._EllipticCurveIsogeny__init_even_kernel_polynomial(E,psig)
            (x^7 + 5*x^6 + 2*x^5 + 6*x^4 + 3*x^3 + 5*x^2 + 6*x + 3, (x^9 + 4*x^8 + 2*x^7 + 4*x^3 + 2*x^2 + x + 6)*y, 1, 6, 3, 4)
        """
        # check if the polynomial really divides the two_torsion_polynomial
        if self.__check and E.division_polynomial(2, x=self.__poly_ring.gen()) % psi_G != 0:
            raise ValueError(f"the polynomial {psi_G} does not define a finite subgroup of {E}")

        n = psi_G.degree() # 1 or 3
        d = n+1            # 2 or 4

        a1, a2, a3, a4, a6 = E.a_invariants()
        b2, b4, _, _ = E.b_invariants()
        x = self.__poly_ring.gen()
        y = self.__mpoly_ring.gen()

        if n == 1:
            x0 = -psi_G.constant_coefficient()

            # determine y0
            if self.__base_field.characteristic() == 2:
                y0 = (x0**3 + a2*x0**2 + a4*x0 + a6).sqrt()
            else:
                y0 = -(a1*x0 + a3)/2

            v,w = compute_vw_kohel_even_deg1(x0, y0, a1, a2, a4)

            phi = (x*psi_G + v)*psi_G
            omega = (y*psi_G**2 - v*(a1*psi_G + (y - y0)))*psi_G

        elif n == 3:
            s1 = -psi_G[n - 1]
            s2 = psi_G[n - 2]
            s3 = -psi_G[n - 3]

            psi_G_pr = psi_G.derivative()
            psi_G_prpr = psi_G_pr.derivative()

            phi = (psi_G_pr**2) + (-2*psi_G_prpr + (4*x - s1))*psi_G
            phi_pr = phi.derivative(x)

            psi_2 = 2*y + a1*x + a3

            omega = (psi_2*(phi_pr*psi_G - phi*psi_G_pr) - (a1*phi + a3*psi_G)*psi_G)/2

            phi *= psi_G
            omega *= psi_G

            v,w = compute_vw_kohel_even_deg3(b2, b4, s1, s2, s3)

        else:
            raise ValueError(f"input polynomial must have degree 1 or 3, not {n}")

        return phi, omega, v, w, n, d

    def __init_odd_kernel_polynomial(self, E, psi):
        r"""
        Return the isogeny parameters for a cyclic isogeny of odd degree.

        INPUT:

        - ``E`` -- an elliptic curve

        - ``psi`` -- a univariate polynomial over the base field of
          ``E``, assumed to be a kernel polynomial

        OUTPUT:

        A tuple (``phi``, ``omega``, ``v``, ``w``, ``n``, ``d``) where:

        - ``phi`` is a univariate polynomial, the numerator of the
          `X`-coordinate of the isogeny;

        - ``omega`` is a bivariate polynomial, the numerator of the
          `Y`-coordinate of the isogeny;

        - ``v``, ``w`` are the Vélu parameters of the isogeny;

        - ``n`` is the degree of ``psi``;

        - ``d`` is the degree of the isogeny.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3); phi
            Isogeny of degree 3
             from Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 1 over Finite Field of size 7
               to Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 4*x + 2 over Finite Field of size 7

            sage: R.<x> = GF(7)[]
            sage: phi._EllipticCurveIsogeny__init_odd_kernel_polynomial(E, x+6)
            (x^3 + 5*x^2 + 3*x + 2, (x^3 + 4*x^2 + x)*y, 2, 6, 1, 3)

            sage: # needs sage.rings.finite_rings
            sage: F = GF(2^4, 'alpha'); R.<x> = F[]
            sage: alpha = F.gen()
            sage: E = EllipticCurve(F, [1,1,F.gen(),F.gen()^2+1,1])
            sage: f = x + alpha^2 + 1
            sage: phi = EllipticCurveIsogeny(E, f); phi
            Isogeny of degree 3
             from Elliptic Curve defined by y^2 + x*y + alpha*y = x^3 + x^2 + (alpha^2+1)*x + 1 over Finite Field in alpha of size 2^4
               to Elliptic Curve defined by y^2 + x*y + alpha*y = x^3 + x^2 + alpha*x + alpha^3 over Finite Field in alpha of size 2^4

            sage: R.<x> = F[]                                                           # needs sage.rings.finite_rings
            sage: f = x + alpha^2 + 1                                                   # needs sage.rings.finite_rings
            sage: phi._EllipticCurveIsogeny__init_odd_kernel_polynomial(E, f)           # needs sage.rings.finite_rings
            (x^3 + (alpha^2 + 1)*x + alpha^3 + alpha^2 + alpha, (x^3 + (alpha^2 + 1)*x^2 + (alpha^2 + 1)*x + alpha)*y + (alpha^2 + alpha + 1)*x^2 + (alpha^2 + alpha)*x + alpha, alpha^2 + alpha + 1, alpha^3 + alpha^2 + alpha, 1, 3)

            sage: E = EllipticCurve(j=-262537412640768000)
            sage: f = E.isogenies_prime_degree()[0].kernel_polynomial()
            sage: f.degree()
            81
            sage: E.isogeny(kernel=f, check=False)
            Isogeny of degree 163
             from Elliptic Curve defined by y^2 + y = x^3 - 2174420*x + 1234136692 over Rational Field
               to Elliptic Curve defined by y^2 + y = x^3 - 57772164980*x - 5344733777551611 over Rational Field
        """
        n = psi.degree()
        d = 2*n + 1

        # check if the polynomial really divides the torsion polynomial :
        if self.__check:
            from .isogeny_small_degree import is_kernel_polynomial
            if not is_kernel_polynomial(E, d, psi):
                raise ValueError(f"the polynomial {psi} does not define a finite subgroup of {E}")

        b2, b4, b6, _ = E.b_invariants()

        s1 = s2 = s3 = 0
        if 1 <= n:
            s1 = -psi[n-1]
        if 2 <= n:
            s2 = psi[n-2]
        if 3 <= n:
            s3 = -psi[n-3]

        # initializing these allows us to calculate E2.
        v, w = compute_vw_kohel_odd(b2, b4, b6, s1, s2, s3, n)

        # initialize the polynomial temporary variables

        psi_pr = psi.derivative()
        psi_prpr = psi_pr.derivative()

        x = self.__poly_ring.gen()

        phi = (4*x**3 + b2*x**2 + 2*b4*x + b6)*(psi_pr**2 - psi_prpr*psi) \
              - (6*x**2 + b2*x + b4)*psi_pr*psi + (d*x - 2*s1)*psi**2

        phi_pr = phi.derivative(x)

        if self.__base_field.characteristic() != 2:
            omega = self.__compute_omega_fast(E, psi, psi_pr, phi, phi_pr)
        else:
            omega = self.__compute_omega_general(E, psi, psi_pr, phi, phi_pr)

        return phi, omega, v, w, n, d

    #
    # This is the fast omega computation that works when characteristic is not 2
    #
    def __compute_omega_fast(self, E, psi, psi_pr, phi, phi_pr):
        r"""
        Return omega from phi, psi and their derivatives, used when
        the characteristic field is not 2.

        INPUT:

        - ``E`` -- an elliptic curve

        - ``psi``, ``psi_pr``, ``phi``, ``phi_pr`` -- univariate
          polynomials over the base field of ``E``, where ``psi``
          is the kernel polynomial and ``phi`` the numerator of
          the `X`-coordinate of the isogeny, together with their
          derivatives

        OUTPUT:

        A bivariate polynomial ``omega`` giving the numerator of
        the `Y`-coordinate of the isogeny.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3); phi
            Isogeny of degree 3
             from Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 1 over Finite Field of size 7
               to Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 4*x + 2 over Finite Field of size 7

            sage: R.<x,y> = GF(7)[]
            sage: psi = phi._EllipticCurveIsogeny__psi
            sage: psi_pr = psi.derivative()
            sage: fi = phi._EllipticCurveIsogeny__phi
            sage: fi_pr = fi.derivative()
            sage: phi._EllipticCurveIsogeny__compute_omega_fast(E, psi, psi_pr, fi, fi_pr)
            (x^3 + 4*x^2 + x)*y
        """
        a1 = E.a1()
        a3 = E.a3()
        x = self.__poly_ring.gen()
        y = self.__mpoly_ring.gen()

        psi_2 = 2*y + a1*x + a3

        # The formula in Kohel's thesis has some typos:
        # Notably, the first plus sign should be a minus
        # as it is below.
        return phi_pr*psi*psi_2/2 - phi*psi_pr*psi_2 - (a1*phi + a3*psi**2)*psi/2

    def __compute_omega_general(self, E, psi, psi_pr, phi, phi_pr):
        r"""
        Return omega from phi, psi and their derivatives, in any
        characteristic.

        INPUT:

        - ``E`` -- an elliptic curve

        - ``psi``, ``psi_pr``, ``phi``, ``phi_pr`` -- univariate polynomials
          over the base field of ``E``, where ``psi`` is the kernel
          polynomial and ``phi`` the numerator of the `X`-coordinate
          of the isogeny, together with their derivatives

        OUTPUT:

        A bivariate polynomial ``omega`` giving the numerator of
        the `Y`-coordinate of the isogeny.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(2^4, 'alpha'); R.<x> = F[]
            sage: alpha = F.gen()
            sage: E = EllipticCurve(F, [1, 1, F.gen(), F.gen()^2+1, 1])
            sage: f = x + alpha^2 + 1
            sage: phi = EllipticCurveIsogeny(E, f); phi
            Isogeny of degree 3
             from Elliptic Curve defined by y^2 + x*y + alpha*y = x^3 + x^2 + (alpha^2+1)*x + 1 over Finite Field in alpha of size 2^4
               to Elliptic Curve defined by y^2 + x*y + alpha*y = x^3 + x^2 + alpha*x + alpha^3 over Finite Field in alpha of size 2^4

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y> = F[]
            sage: psi = phi._EllipticCurveIsogeny__psi
            sage: psi_pr = psi.derivative()
            sage: fi = phi._EllipticCurveIsogeny__phi
            sage: fi_pr = fi.derivative()
            sage: phi._EllipticCurveIsogeny__compute_omega_general(E, psi, psi_pr, fi, fi_pr)
            (x^3 + (alpha^2 + 1)*x^2 + (alpha^2 + 1)*x + alpha)*y + (alpha^2 + alpha + 1)*x^2 + (alpha^2 + alpha)*x + alpha

        A bug fixed in :issue:`7907`::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(128,'a')
            sage: a = F.gen()
            sage: E = EllipticCurve([1,0,0,0,(a**6+a**4+a**2+a)])
            sage: x = polygen(F)
            sage: ker = (x^6 + (a^6 + a^5 + a^4 + a^3 + a^2 + a)*x^5 + (a^6 + a^5 + a^2 + 1)*x^4
            ....:        + (a^6 + a^5 + a^4 + a^3 + a^2 + 1)*x^3 + (a^6 + a^3 + a)*x^2 + (a^4 + a^3 + 1)*x + a^5 + a^4 + a)
            sage: E.isogeny(ker)
            Isogeny of degree 13
             from Elliptic Curve defined by y^2 + x*y = x^3 + (a^6+a^4+a^2+a) over Finite Field in a of size 2^7
               to Elliptic Curve defined by y^2 + x*y = x^3 + (a^6+a^5+a^4+a^3+a^2+a)*x + (a^5+a^3) over Finite Field in a of size 2^7
        """
        a1, a2, a3, a4, a6 = E.a_invariants()
        b2, b4, _, _ = E.b_invariants()
        x = self.__poly_ring.gen()
        y = self.__mpoly_ring.gen()

        n = psi.degree()
        d = 2 * n + 1

        s1 = -psi[n-1] if n > 0 else 0

        psi_prpr = 0
        cur_x_pow = 1

        # Note: we now get the "derivatives" of psi
        # these are not actually the derivatives
        # furthermore, the formulas in Kohel's
        # thesis are wrong, the correct formulas
        # are coded below

        from sage.arith.misc import binomial

        for j in range(n - 1):
            psi_prpr += binomial(j+2, 2) * psi[j+2] * cur_x_pow
            cur_x_pow = x * cur_x_pow

        psi_prprpr = 0
        cur_x_pow = 1

        for j in range(n - 2):
            psi_prprpr += (3 * binomial(j+3, 3)) * psi[j+3] * cur_x_pow
            cur_x_pow = x * cur_x_pow

        psi_2 = 2 * y + a1 * x + a3

        omega = phi_pr*psi*y - phi*psi_pr*psi_2 \
                + ((a1*x + a3)*(psi_2**2)*(psi_prpr*psi_pr-psi_prprpr*psi)
                  + (a1*psi_2**2 - 3*(a1*x + a3)*(6*x**2 + b2*x + b4))*psi_prpr*psi
                  + (a1*x**3 + 3*a3*x**2 + (2*a2*a3 - a1*a4)*x + (a3*a4 - 2*a1*a6))*psi_pr**2
                  + (-(3*a1*x**2 + 6*a3*x + (-a1*a4 + 2*a2*a3))
                  + (a1*x + a3)*(d*x - 2*s1) )*psi_pr*psi + (a1*s1 + a3*n)*psi**2) * psi

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
            ()
        """
        # first check if the point is in the kernel
        if self.__inner_kernel_polynomial(xP) == 0:
            return ()

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
        omega0 = self.__omega[0]
        omega1 = self.__omega[1]
        b = omega0(xP) + omega1(xP)*yP
        c = self.__psi(xP)
        return a/c**2, b/c**3

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
        return self.__compute_via_kohel(x, y)

    #
    # Kohel's formula computing the codomain curve
    #
    def __compute_codomain_via_kohel(self):
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
            sage: phi._EllipticCurveIsogeny__compute_codomain_via_kohel()
            Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 4*x + 2 over Finite Field of size 7
        """
        return compute_codomain_formula(self._domain, self.__v, self.__w)

    #
    # public isogeny methods
    #

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
            sage: phi = EllipticCurveIsogeny(E, E((0,0)))
            sage: phi.rational_maps()
            ((x^2 + 3)/x, (x^2*y - 3*y)/x^2)
        """
        self.__initialize_rational_maps()
        return tuple(self.__xyfield(f) for f in self.__ratl_maps)

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
            sage: phi = EllipticCurveIsogeny(E, E((0,0)))
            sage: phi.x_rational_map()
            (x^2 + 3)/x
        """
        self.__initialize_rational_maps()
        return self.__ratl_maps[0]

    def scaling_factor(self):
        r"""
        Return the Weierstrass scaling factor associated to this
        elliptic-curve isogeny.

        The scaling factor is the constant `u` (in the base field)
        such that `\varphi^* \omega_2 = u \omega_1`, where
        `\varphi: E_1\to E_2` is this isogeny and `\omega_i` are
        the standard Weierstrass differentials on `E_i` defined by
        `\mathrm dx/(2y+a_1x+a_3)`.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve(GF(257^2), [0,1])
            sage: phi = E.isogeny(E.lift_x(240))
            sage: phi.degree()
            43
            sage: phi.scaling_factor()
            1
            sage: phi.dual().scaling_factor()
            43

        TESTS:

        Check for :issue:`36638`::

            sage: phi.scaling_factor().parent()                                         # needs sage.rings.finite_rings
            Finite Field in z2 of size 257^2

        ALGORITHM: The "inner" isogeny is normalized by construction,
        so we only need to account for the scaling factors of a pre-
        and post-isomorphism.
        """
        sc = self.__base_field.one()
        if self.__pre_isomorphism is not None:
            sc *= self.__pre_isomorphism.scaling_factor()
        if self.__post_isomorphism is not None:
            sc *= self.__post_isomorphism.scaling_factor()
        return sc

    def kernel_polynomial(self):
        r"""
        Return the kernel polynomial of this isogeny.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,0,0,2,0])
            sage: phi = EllipticCurveIsogeny(E, E((0,0)))
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

    def inseparable_degree(self):
        r"""
        Return the inseparable degree of this isogeny.

        Since this class only implements separable isogenies,
        this method always returns one.

        TESTS::

            sage: EllipticCurveIsogeny.inseparable_degree(None)
            1
        """
        return Integer(1)

    @cached_method
    def dual(self):
        r"""
        Return the isogeny dual to this isogeny.

        .. NOTE::

            If `\varphi\colon E \to E'` is the given isogeny and `n`
            is its degree, then the dual is by definition the unique
            isogeny `\hat\varphi\colon E'\to E` such that the
            compositions `\hat\varphi\circ\varphi` and
            `\varphi\circ\hat\varphi` are the multiplication-by-`n`
            maps on `E` and `E'`, respectively.

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

        Inseparable duals should be computed correctly::

            sage: # needs sage.rings.finite_rings
            sage: z2 = GF(71^2).gen()
            sage: E = EllipticCurve(j=57*z2+51)
            sage: E.isogeny(3*E.lift_x(0)).dual()
            Composite morphism of degree 71 = 71*1:
              From: Elliptic Curve defined by y^2 = x^3 + (32*z2+67)*x + (24*z2+37)
                    over Finite Field in z2 of size 71^2
              To:   Elliptic Curve defined by y^2 = x^3 + (41*z2+56)*x + (18*z2+42)
                    over Finite Field in z2 of size 71^2
            sage: E.isogeny(E.lift_x(0)).dual()
            Composite morphism of degree 213 = 71*3:
              From: Elliptic Curve defined by y^2 = x^3 + (58*z2+31)*x + (34*z2+58)
                    over Finite Field in z2 of size 71^2
              To:   Elliptic Curve defined by y^2 = x^3 + (41*z2+56)*x + (18*z2+42)
                    over Finite Field in z2 of size 71^2

        ...even if pre- or post-isomorphisms are present::

            sage: # needs sage.rings.finite_rings
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: phi = E.isogeny(E.lift_x(0))
            sage: pre = ~WeierstrassIsomorphism(phi.domain(), (z2,2,3,4))
            sage: post = WeierstrassIsomorphism(phi.codomain(), (5,6,7,8))
            sage: phi = post * phi * pre
            sage: phi.dual()
            Composite morphism of degree 213 = 71*3:
              From: Elliptic Curve defined
                    by y^2 + 17*x*y + 45*y = x^3 + 30*x^2 + (6*z2+64)*x + (48*z2+65)
                    over Finite Field in z2 of size 71^2
              To:   Elliptic Curve defined
                    by y^2 + (60*z2+22)*x*y + (69*z2+37)*y = x^3 + (32*z2+48)*x^2
                       + (19*z2+58)*x + (56*z2+22)
                    over Finite Field in z2 of size 71^2

        TESTS:

        Test for :issue:`23928`::

            sage: E = EllipticCurve(j=GF(431**2)(4))                                    # needs sage.rings.finite_rings
            sage: phi = E.isogeny(E.lift_x(0))                                          # needs sage.rings.finite_rings
            sage: phi.dual()                                                            # needs sage.rings.finite_rings
            Isogeny of degree 2
             from Elliptic Curve defined by y^2 = x^3 + 427*x over Finite Field in z2 of size 431^2
               to Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 431^2

        Test (for :issue:`7096`)::

            sage: E = EllipticCurve('11a1')
            sage: phi = E.isogeny(E(5,5))
            sage: phi.dual().dual() == phi
            True

            sage: k = GF(103)
            sage: E = EllipticCurve(k,[11,11])
            sage: phi = E.isogeny(E(4,4))
            sage: phi
            Isogeny of degree 5
             from Elliptic Curve defined by y^2 = x^3 + 11*x + 11 over Finite Field of size 103
               to Elliptic Curve defined by y^2 = x^3 + 25*x + 80 over Finite Field of size 103
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: phi = WeierstrassIsomorphism(phi.codomain(),(5,0,1,2)) * phi
            sage: phi.dual().dual() == phi
            True

            sage: E = EllipticCurve(GF(103),[1,0,0,1,-1])
            sage: phi = E.isogeny(E(60,85))
            sage: phi.dual()
            Isogeny of degree 7
             from Elliptic Curve defined by y^2 + x*y = x^3 + 84*x + 34 over Finite Field of size 103
               to Elliptic Curve defined by y^2 + x*y = x^3 + x + 102 over Finite Field of size 103

        Check that :issue:`17293` is fixed::

            sage: # needs sage.rings.number_field
            sage: k.<s> = QuadraticField(2)
            sage: E = EllipticCurve(k, [-3*s*(4 + 5*s), 2*s*(2 + 14*s + 11*s^2)])
            sage: phi = E.isogenies_prime_degree(3)[0]
            sage: (-phi).dual() == -phi.dual()
            True

        Check that :issue:`37168` is fixed::

            sage: R.<x> = GF(23)[]
            sage: F.<a> = FiniteField(23^2, modulus=x^2-x+1)
            sage: E0 = EllipticCurve(F, (0, 1))
            sage: E1 = EllipticCurve(F, (8, 1))
            sage: phi = E0.isogeny(kernel=E0((a, 0)), codomain=E1)
            sage: phi.dual()
            Isogeny of degree 2
             from Elliptic Curve defined by y^2 = x^3 + 8*x + 1 over Finite Field in a of size 23^2
             to Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field in a of size 23^2
        """
        if self.__base_field.characteristic() in (2, 3):
            raise NotImplementedError("computation of dual isogenies not yet implemented in characteristics 2 and 3")

        # trac 7096
        E1, E2pr, _, _ = compute_intermediate_curves(self.codomain(), self.domain())

        F = self.__base_field
        d = self._degree

        from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite

        if F(d) == 0:   # inseparable dual!
            p = F.characteristic()
            k = d.valuation(p)

            from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius
            frob = EllipticCurveHom_frobenius(self._codomain, k)

            dsep = d // p**k
            if dsep > 1:
                #TODO: We could also use resultants here; this is much
                # faster in some cases (but seems worse in general).
                # Presumably there should be a wrapper function that
                # decides on the fly which method to use.
                # Eventually this should become a .separable_part() method.

                f = self.kernel_polynomial()

                psi = self._domain.division_polynomial(p)
                mu_num = self._domain._multiple_x_numerator(p)
                mu_den = self._domain._multiple_x_denominator(p)

                for _ in range(k):
                    f //= f.gcd(psi)
                    S = f.parent().quotient_ring(f)
                    mu = S(mu_num) / S(mu_den)
                    f = mu.minpoly()

                sep = self._domain.isogeny(f, codomain=frob.codomain()).dual()

            else:
                sep = frob.codomain().isomorphism_to(self._domain)

            phi_hat = EllipticCurveHom_composite.from_factors([frob, sep])

            from sage.schemes.elliptic_curves.hom import find_post_isomorphism
            mult = self._domain.scalar_multiplication(d)
            rhs = phi_hat * self
            corr = find_post_isomorphism(mult, rhs)
            return corr * phi_hat

        else:
            # trac 7096
            # this should take care of the case when the isogeny is not normalized.
            u = self.scaling_factor()
            E2 = E2pr.change_weierstrass_model(u/F(d), 0, 0, 0)

            phi_hat = EllipticCurveIsogeny(E1, None, E2, d)
#            assert phi_hat.scaling_factor() == 1

            for pre_iso in self._codomain.isomorphisms(E1):
                for post_iso in E2.isomorphisms(self._domain):
                    sc = u * pre_iso.scaling_factor() * post_iso.scaling_factor()
                    if sc == d:
                        break
                else:
                    continue
                break
            else:
                assert "bug in dual()"

            return phi_hat._compose_with_isomorphism(pre_iso, post_iso)

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
            Isogeny of degree 11
             from Elliptic Curve defined by y^2 + 2*x*y + 2*y = x^3 + 2*x^2 + 6*x + 7 over Finite Field of size 127
               to Elliptic Curve defined by y^2 = x^3 + 37*x + 85 over Finite Field of size 127
            sage: iso2 * phi            # indirect doctest
            Isogeny of degree 11
             from Elliptic Curve defined by y^2 = x^3 + 5*x + 2 over Finite Field of size 127
               to Elliptic Curve defined by y^2 = x^3 + 117*x + 58 over Finite Field of size 127
            sage: iso2 * phi * iso1     # indirect doctest
            Isogeny of degree 11
             from Elliptic Curve defined by y^2 + 2*x*y + 2*y = x^3 + 2*x^2 + 6*x + 7 over Finite Field of size 127
               to Elliptic Curve defined by y^2 = x^3 + 117*x + 58 over Finite Field of size 127

        TESTS:

        We should return ``NotImplemented`` when passed a combination of
        elliptic-curve morphism types that we don't handle here::

            sage: phi._composition_impl(iso1, iso1**-1)
            NotImplemented
        """
        if isinstance(left, WeierstrassIsomorphism) and isinstance(right, EllipticCurveIsogeny):
            return right._compose_with_isomorphism(post_isomorphism=left)

        if isinstance(left, EllipticCurveIsogeny) and isinstance(right, WeierstrassIsomorphism):
            return left._compose_with_isomorphism(pre_isomorphism=right)

        return NotImplemented


def compute_isogeny_bmss(E1, E2, l):
    r"""
    Compute the kernel polynomial of the unique normalized isogeny
    of degree ``l`` between ``E1`` and ``E2``.

    Both curves must be given in short Weierstrass form, and the
    characteristic must be either `0` or no smaller than `4l+4`.

    ALGORITHM: [BMSS2006]_, algorithm *fastElkies'*.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_bmss
        sage: E1 = EllipticCurve(GF(167), [153, 112])
        sage: E2 = EllipticCurve(GF(167), [56, 40])
        sage: compute_isogeny_bmss(E1, E2, 13)
        x^6 + 139*x^5 + 73*x^4 + 139*x^3 + 120*x^2 + 88*x
    """
    # Original author of this function: Rémy Oudompheng.
    # https://github.com/remyoudompheng/isogeny_weber/blob/64289127a337ac1bf258b711e02fea02b7df5275/isogeny_weber/isogenies.py#L272-L332
    # Released under the MIT license: https://github.com/remyoudompheng/isogeny_weber/blob/64289127a337ac1bf258b711e02fea02b7df5275/LICENSE
    # Slightly adjusted for inclusion in the Sage library.
    if E1.a1() or E1.a2() or E1.a3():
        raise ValueError('E1 must be a short Weierstrass curve')
    if E2.a1() or E2.a2() or E2.a3():
        raise ValueError('E2 must be a short Weierstrass curve')
    char = E1.base_ring().characteristic()
    if char != 0 and char < 4*l + 4:
        raise ValueError('characteristic must be at least 4*degree+4')
    Rx, x = E1.base_ring()["x"].objgen()
    # Compute C = 1/(1 + Ax^4 + Bx^6) mod x^4l
    A, B = E1.a4(), E1.a6()
    C = (1 + A * x**4 + B * x**6).inverse_series_trunc(4 * l)
    # Solve differential equation
    # The number of terms doubles at each iteration.
    # S'^2 = G(x,S) = (1 + A2 S^4 + B2 S^6) / (1 + Ax^4 + Bx^6)
    # S = x + O(x^2)
    A2, B2 = E2.a4(), E2.a6()
    S = x + (A2 - A) / 10 * x**5 + (B2 - B) / 14 * x**7
    sprec = 8
    while sprec < 4 * l:
        assert sprec % 2 == 0
        sprec = min(sprec, 2 * l)
        # s1 => s1 + x^k s2
        # 2 s1' s2' - dG/dS(x, s1) s2 = G(x, s1) - s1'2
        s1 = S
        ds1 = s1.derivative()
        s1pows = [1, s1]
        while len(s1pows) < 7:
            s1pows.append(s1._mul_trunc_(s1pows[-1], 2 * sprec))
        GS = C * (1 + A2 * s1pows[4] + B2 * s1pows[6])
        dGS = C * (4 * A2 * s1pows[3] + 6 * B2 * s1pows[5])
        # s2' = (dGS / 2s1') s2 + (G(x, s1) - s1'2) / (2s1')
        denom = (2 * ds1).inverse_series_trunc(2 * sprec)
        a = dGS._mul_trunc_(denom, 2 * sprec)
        b = (GS - ds1**2)._mul_trunc_(denom, 2 * sprec)
        s2 = a.add_bigoh(2 * sprec).solve_linear_de(prec=2 * sprec, b=b, f0=0)
        S = s1 + Rx(s2)
        sprec = 2 * sprec
    # Reconstruct:
    # S = x * T(x^2)
    # Compute U = 1/T^2
    # Reconstruct N(1/x) / D(1/x) = U
    T = Rx([S[2 * i + 1] for i in range(2 * l)])
    U = T._mul_trunc_(T, 2 * l).inverse_series_trunc(2 * l)
    _, Q = Rx(U).rational_reconstruction(x ** (2 * l), l, l)
    Q = Q.add_bigoh((l + 1) // 2)
    if not Q.is_square():
        raise ValueError(f"the two curves are not linked by a cyclic normalized isogeny of degree {l}")
    Q = Q.sqrt()
    ker = Rx(Q).reverse(degree=l//2)
    return ker.monic()


def compute_isogeny_stark(E1, E2, ell):
    r"""
    Return the kernel polynomial of an isogeny of degree ``ell``
    from ``E1`` to ``E2``.

    INPUT:

    - ``E1`` -- domain elliptic curve in short Weierstrass form

    - ``E2`` -- codomain elliptic curve in short Weierstrass form

    - ``ell`` -- the degree of an isogeny from ``E1`` to ``E2``

    OUTPUT: the kernel polynomial of an isogeny from ``E1`` to ``E2``

    .. NOTE::

        If there is no degree-``ell``, cyclic, separable, normalized
        isogeny from ``E1`` to ``E2``, a :exc:`ValueError` will be
        raised.

    ALGORITHM:

    This function uses Stark's algorithm as presented in Section 6.2
    of [BMSS2006]_.

    .. NOTE::

        As published in [BMSS2006]_, the algorithm is incorrect, and a
        correct version (with slightly different notation) can be found
        in [Mo2009]_.  The algorithm originates in [Sta1973]_.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_stark, compute_sequence_of_maps

        sage: E = EllipticCurve(GF(97), [1,0,1,1,0])
        sage: R.<x> = GF(97)[]; f = x^5 + 27*x^4 + 61*x^3 + 58*x^2 + 28*x + 21
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: isom1, isom2, E1pr, E2pr, ker_poly = compute_sequence_of_maps(E, E2, 11)
        sage: compute_isogeny_stark(E1pr, E2pr, 11)
        x^10 + 37*x^9 + 53*x^8 + 66*x^7 + 66*x^6 + 17*x^5 + 57*x^4 + 6*x^3 + 89*x^2 + 53*x + 8

        sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
        sage: R.<x> = GF(37)[]
        sage: f = (x + 14) * (x + 30)
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_isogeny_stark(E, E2, 5)
        x^4 + 14*x^3 + x^2 + 34*x + 21
        sage: f**2
        x^4 + 14*x^3 + x^2 + 34*x + 21

        sage: E = EllipticCurve(QQ, [0,0,0,1,0])
        sage: R.<x> = QQ[]
        sage: f = x
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_isogeny_stark(E, E2, 2)
        x

    TESTS:

    Check for :issue:`21883`::

        sage: E1 = EllipticCurve([0,1])
        sage: E2 = EllipticCurve([0,-27])
        sage: E1.isogeny(None, E2, degree=3)
        Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field to Elliptic Curve defined by y^2 = x^3 - 27 over Rational Field
    """
    K = E1.base_field()
    R, x = PolynomialRing(K, 'x').objgen()

    wp1 = E1.weierstrass_p(prec=4*ell+4)  #BMSS2006 claim 2*ell is enough, but it is not M09
    wp2 = E2.weierstrass_p(prec=4*ell+4)

    # viewed them as power series in Z = z^2
    Z = LaurentSeriesRing(K, 'Z').gen()
    pe1 = pe2 = 1/Z
    for i in range(2*ell + 1):
        pe1 += wp1[2*i] * Z**i
        pe2 += wp2[2*i] * Z**i
    pe1 = pe1.add_bigoh(2*ell+3)
    pe2 = pe2.add_bigoh(2*ell+3)

    n = 1
    q = [R.one(), R.zero()]
    T = pe2

    while q[n].degree() < ell-1:
        n += 1
        a_n = 0
        r = -T.valuation()
        while 0 <= r:
            t_r = T[-r]
            a_n = a_n + t_r * x**r
            T = T - t_r*pe1**r
            r = -T.valuation()

        q_n = a_n*q[n-1] + q[n-2]
        q.append(q_n)

        if n == ell+1 or T == 0:
            if T == 0 or T.valuation() < 2:
                raise ValueError(f"the two curves are not linked by a cyclic normalized isogeny of degree {ell}")
            break

        T = 1/T

    qn = q[n]
    qn /= qn.leading_coefficient()
    return qn


def compute_isogeny_kernel_polynomial(E1, E2, ell, algorithm=None):
    r"""
    Return the kernel polynomial of a cyclic, separable, normalized
    isogeny of degree ``ell`` from ``E1`` to ``E2``.

    INPUT:

    - ``E1`` -- domain elliptic curve in short Weierstrass form

    - ``E2`` -- codomain elliptic curve in short Weierstrass form

    - ``ell`` -- the degree of an isogeny from ``E1`` to ``E2``

    - ``algorithm`` -- ``None`` (default, choose automatically) or
      ``'bmss'`` (:func:`compute_isogeny_bmss`) or
      ``'stark'`` (:func:`compute_isogeny_stark`)

    OUTPUT: the kernel polynomial of an isogeny from ``E1`` to ``E2``

    .. NOTE::

        If there is no degree-``ell``, cyclic, separable, normalized
        isogeny from ``E1`` to ``E2``, a :exc:`ValueError` will be
        raised.

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

        sage: # needs sage.rings.number_field
        sage: R.<x> = QQ[]
        sage: K.<i> = NumberField(x^2 + 1)
        sage: E = EllipticCurve(K, [0,0,0,1,0])
        sage: E2 = EllipticCurve(K, [0,0,0,16,0])
        sage: compute_isogeny_kernel_polynomial(E, E2, 4)
        x^3 + x

    TESTS:

    Check that :meth:`Polynomial.radical` is doing the right thing for us::

        sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
        sage: R.<x> = GF(37)[]
        sage: f = (x + 10) * (x + 12) * (x + 16)
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_stark
        sage: ker_poly = compute_isogeny_stark(E, E2, 7); ker_poly
        x^6 + 2*x^5 + 20*x^4 + 11*x^3 + 36*x^2 + 35*x + 16
        sage: ker_poly.factor()
        (x + 10)^2 * (x + 12)^2 * (x + 16)^2
        sage: poly = ker_poly.radical(); poly
        x^3 + x^2 + 28*x + 33
        sage: poly.factor()
        (x + 10) * (x + 12) * (x + 16)
    """
    if algorithm is None:
        char = E1.base_ring().characteristic()
        if char != 0 and char < 4*ell + 4:
            raise NotImplementedError(f'no algorithm for computing kernel polynomial from domain and codomain is implemented for degree {ell} and characteristic {char}')
        algorithm = 'stark' if ell < 10 else 'bmss'

    if algorithm == 'bmss':
        return compute_isogeny_bmss(E1, E2, ell)
    if algorithm == 'stark':
        return compute_isogeny_stark(E1, E2, ell).radical()

    raise NotImplementedError(f'unknown algorithm {algorithm}')


def compute_intermediate_curves(E1, E2):
    r"""
    Return intermediate curves and isomorphisms.

    .. NOTE::

        This is used to compute `\wp` functions from the short
        Weierstrass model more easily.

    .. WARNING::

        The base field must be of characteristic not equal to `2` or `3`.

    INPUT:

    - ``E1``, ``E2`` -- elliptic curves

    OUTPUT:

    A tuple (``pre_isomorphism``, ``post_isomorphism``,
    ``intermediate_domain``, ``intermediate_codomain``) where:

    - ``intermediate_domain`` is a short Weierstrass curve isomorphic
      to ``E1``;

    - ``intermediate_codomain`` is a short Weierstrass curve isomorphic
      to ``E2``;

    - ``pre_isomorphism`` is a normalized isomorphism from ``E1`` to
      ``intermediate_domain``;

    - ``post_isomorphism`` is a normalized isomorphism from
      ``intermediate_codomain`` to ``E2``.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_intermediate_curves
        sage: E = EllipticCurve(GF(83), [1,0,1,1,0])
        sage: R.<x> = GF(83)[]; f = x + 24
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_intermediate_curves(E, E2)
        (Elliptic Curve defined by y^2 = x^3 + 62*x + 74 over Finite Field of size 83,
         Elliptic Curve defined by y^2 = x^3 + 65*x + 69 over Finite Field of size 83,
         Elliptic-curve morphism:
          From: Elliptic Curve defined by y^2 + x*y + y = x^3 + x
                over Finite Field of size 83
          To:   Elliptic Curve defined by y^2 = x^3 + 62*x + 74
                over Finite Field of size 83
          Via:  (u,r,s,t) = (1, 76, 41, 3),
         Elliptic-curve morphism:
          From: Elliptic Curve defined by y^2 = x^3 + 65*x + 69
                over Finite Field of size 83
          To:   Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + 16
                over Finite Field of size 83
          Via:  (u,r,s,t) = (1, 7, 42, 42))

        sage: # needs sage.rings.number_field
        sage: R.<x> = QQ[]
        sage: K.<i> = NumberField(x^2 + 1)
        sage: E = EllipticCurve(K, [0,0,0,1,0])
        sage: E2 = EllipticCurve(K, [0,0,0,16,0])
        sage: compute_intermediate_curves(E, E2)
        (Elliptic Curve defined by y^2 = x^3 + x
          over Number Field in i with defining polynomial x^2 + 1,
         Elliptic Curve defined by y^2 = x^3 + 16*x
          over Number Field in i with defining polynomial x^2 + 1,
         Elliptic-curve endomorphism of Elliptic Curve defined by y^2 = x^3 + x
          over Number Field in i with defining polynomial x^2 + 1
          Via:  (u,r,s,t) = (1, 0, 0, 0),
         Elliptic-curve endomorphism of Elliptic Curve defined by y^2 = x^3 + 16*x
          over Number Field in i with defining polynomial x^2 + 1
          Via:  (u,r,s,t) = (1, 0, 0, 0))
    """
    if E1.base_ring().characteristic() in (2, 3):
        raise NotImplementedError("compute_intermediate_curves is only defined for characteristics not 2 or 3")

    # We cannot just use
    # E1w = E1.short_weierstrass_model()
    # E2w = E2.short_weierstrass_model()
    # as the resulting isomorphisms would not be normalised (u=1)

    c4, c6 = E1.c_invariants()
    E1w = EllipticCurve([0, 0, 0, -c4/48, -c6/864])
    c4, c6 = E2.c_invariants()
    E2w = EllipticCurve([0, 0, 0, -c4/48, -c6/864])

    # We cannot even just use pre_iso = E1.isomorphism_to(E1w) since
    # it may have u=-1; similarly for E2

    urst = [w for w in _isomorphisms(E1, E1w) if w[0] == 1][0]
    pre_iso = WeierstrassIsomorphism(E1, urst, E1w)
    urst = [w for w in _isomorphisms(E2w, E2) if w[0] == 1][0]
    post_iso = WeierstrassIsomorphism(E2w, urst, E2)
    return E1w, E2w, pre_iso, post_iso


def compute_sequence_of_maps(E1, E2, ell):
    r"""
    Return intermediate curves, isomorphisms and kernel polynomial.

    INPUT:

    - ``E1``, ``E2`` -- elliptic curves

    - ``ell`` -- a prime such that there is a degree-``ell`` separable
      normalized isogeny from ``E1`` to ``E2``

    OUTPUT:

    A tuple (``pre_isom``, ``post_isom``, ``E1pr``, ``E2pr``,
    ``ker_poly``) where:

    - ``E1pr`` is an elliptic curve in short Weierstrass form
      isomorphic to ``E1``;

    - ``E2pr`` is an elliptic curve in short Weierstrass form
      isomorphic to ``E2``;

    - ``pre_isom`` is a normalized isomorphism from ``E1`` to ``E1pr``;

    - ``post_isom`` is a normalized isomorphism from ``E2pr`` to ``E2``;

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
          From: Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20
                over Rational Field
          To:   Elliptic Curve defined by y^2 = x^3 - 31/3*x - 2501/108
                over Rational Field
          Via:  (u,r,s,t) = (1, 1/3, 0, -1/2),
         Elliptic-curve morphism:
          From: Elliptic Curve defined by y^2 = x^3 - 23461/3*x - 28748141/108
                over Rational Field
          To:   Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580
                over Rational Field
          Via:  (u,r,s,t) = (1, -1/3, 0, 1/2),
         Elliptic Curve defined by y^2 = x^3 - 31/3*x - 2501/108 over Rational Field,
         Elliptic Curve defined by y^2 = x^3 - 23461/3*x - 28748141/108 over Rational Field,
         x^2 - 61/3*x + 658/9)

        sage: # needs sage.rings.number_field
        sage: K.<i> = NumberField(x^2 + 1)
        sage: E = EllipticCurve(K, [0,0,0,1,0])
        sage: E2 = EllipticCurve(K, [0,0,0,16,0])
        sage: compute_sequence_of_maps(E, E2, 4)
        (Elliptic-curve endomorphism of Elliptic Curve defined by y^2 = x^3 + x
          over Number Field in i with defining polynomial x^2 + 1
          Via:  (u,r,s,t) = (1, 0, 0, 0),
         Elliptic-curve endomorphism of Elliptic Curve defined by y^2 = x^3 + 16*x
          over Number Field in i with defining polynomial x^2 + 1
          Via:  (u,r,s,t) = (1, 0, 0, 0),
         Elliptic Curve defined by y^2 = x^3 + x
          over Number Field in i with defining polynomial x^2 + 1,
         Elliptic Curve defined by y^2 = x^3 + 16*x
          over Number Field in i with defining polynomial x^2 + 1,
         x^3 + x)

        sage: E = EllipticCurve(GF(97), [1,0,1,1,0])
        sage: R.<x> = GF(97)[]; f = x^5 + 27*x^4 + 61*x^3 + 58*x^2 + 28*x + 21
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_sequence_of_maps(E, E2, 11)
        (Elliptic-curve morphism:
          From: Elliptic Curve defined by y^2 + x*y + y = x^3 + x
                over Finite Field of size 97
          To:   Elliptic Curve defined by y^2 = x^3 + 52*x + 31
                over Finite Field of size 97
          Via:  (u,r,s,t) = (1, 8, 48, 44),
         Elliptic-curve morphism:
          From: Elliptic Curve defined by y^2 = x^3 + 41*x + 66
                over Finite Field of size 97
          To:   Elliptic Curve defined by y^2 + x*y + y = x^3 + 87*x + 26
                over Finite Field of size 97
          Via:  (u,r,s,t) = (1, 89, 49, 49),
         Elliptic Curve defined by y^2 = x^3 + 52*x + 31 over Finite Field of size 97,
         Elliptic Curve defined by y^2 = x^3 + 41*x + 66 over Finite Field of size 97,
         x^5 + 67*x^4 + 13*x^3 + 35*x^2 + 77*x + 69)
    """
    E1pr, E2pr, pre_isom, post_isom = compute_intermediate_curves(E1, E2)

    ker_poly = compute_isogeny_kernel_polynomial(E1pr, E2pr, ell)

    return pre_isom, post_isom, E1pr, E2pr, ker_poly

# Utility functions for manipulating isogeny degree matrices


def fill_isogeny_matrix(M):
    """
    Return a filled isogeny matrix giving all degrees from one giving only prime degrees.

    INPUT:

    - ``M`` -- a square symmetric matrix whose off-diagonal `i`, `j`
      entry is either a prime `l` if the `i`-th and `j`-th curves
      have an `l`-isogeny between them, otherwise `0`

    OUTPUT:

    A square matrix with entries `1` on the diagonal, and in
    general the `i`, `j` entry is `d>0` if `d` is the minimal degree
    of an isogeny from the `i`-th to the `j`-th curve.

    EXAMPLES::

        sage: M = Matrix([[0, 2, 3, 3, 0, 0], [2, 0, 0, 0, 3, 3], [3, 0, 0, 0, 2, 0],
        ....:             [3, 0, 0, 0, 0, 2], [0, 3, 2, 0, 0, 0], [0, 3, 0, 2, 0, 0]]); M
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
    from sage.matrix.constructor import Matrix
    from sage.rings.infinity import Infinity

    n = M.nrows()
    M0 = copy(M)
    for i in range(n):
        M0[i,i] = 1

    def fix(d):
        return d if d != 0 else Infinity

    def fix2(d):
        return d if d != Infinity else 0

    def pr(M1, M2):
        return Matrix([[fix2(min([fix(M1[i,k]*M2[k,j]) for k in range(n)])) for i in range(n)] for j in range(n)])

    M1 = M0
    M2 = pr(M0, M1)
    while M1 != M2:
        M1 = M2
        M2 = pr(M0, M1)
    return M1


def unfill_isogeny_matrix(M):
    """
    Reverses the action of ``fill_isogeny_matrix``.

    INPUT:

    - ``M`` -- a square symmetric matrix of integers

    OUTPUT:

    A square symmetric matrix obtained from ``M`` by
    replacing non-prime entries with `0`.

    EXAMPLES::

        sage: M = Matrix([[0, 2, 3, 3, 0, 0], [2, 0, 0, 0, 3, 3], [3, 0, 0, 0, 2, 0],
        ....:             [3, 0, 0, 0, 0, 2], [0, 3, 2, 0, 0, 0], [0, 3, 0, 2, 0, 0]]); M
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
        M1[i,i] = zero
        for j in range(i):
            if not M1[i,j].is_prime():
                M1[i,j] = zero
                M1[j,i] = zero
    return M1
