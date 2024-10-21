# sage_setup: distribution = sagemath-schemes
r"""
Mestre's algorithm

This file contains functions that:

- create hyperelliptic curves from the Igusa-Clebsch invariants (over
  `\QQ` and finite fields)
- create Mestre's conic from the Igusa-Clebsch invariants

AUTHORS:

- Florian Bouyer
- Marco Streng
"""
#*****************************************************************************
#       Copyright (C) 2011, 2012, 2013
#                  Florian Bouyer <f.j.s.c.bouyer@gmail.com>
#                  Marco Streng <marco.streng@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.constructor import Matrix
from sage.schemes.plane_conics.constructor import Conic
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve


def HyperellipticCurve_from_invariants(i, reduced=True, precision=None,
                                       algorithm='default'):
    r"""
    Return a hyperelliptic curve with the given Igusa-Clebsch invariants up to
    scaling.

    The output is a curve over the field in which the Igusa-Clebsch invariants
    are given. The output curve is unique up to isomorphism over the algebraic
    closure. If no such curve exists over the given field, then raise a
    :exc:`ValueError`.

    INPUT:

    - ``i`` -- list or tuple of length 4 containing the four Igusa-Clebsch
      invariants: I2,I4,I6,I10
    - ``reduced`` -- boolean (default: ``True``); if ``True``, tries to reduce
      the polynomial defining the hyperelliptic curve using the function
      :func:`reduce_polynomial` (see the :func:`reduce_polynomial`
      documentation for more details).
    - ``precision`` -- integer (default: ``None``); which precision for real and
      complex numbers should the reduction use. This only affects the
      reduction, not the correctness. If ``None``, the algorithm uses the default
      53 bit precision.
    - ``algorithm`` -- ``'default'`` or ``'magma'``. If set to ``'magma'``, uses
      Magma to parameterize Mestre's conic (needs Magma to be installed)

    OUTPUT: a hyperelliptic curve object

    EXAMPLES:

    Examples over the rationals::

        sage: HyperellipticCurve_from_invariants([3840,414720,491028480,2437709561856])
        Traceback (most recent call last):
        ...
        NotImplementedError: Reduction of hyperelliptic curves not yet implemented.
        See issues #14755 and #14756.

        sage: HyperellipticCurve_from_invariants([3840,414720,491028480,2437709561856], reduced=False)
        Hyperelliptic Curve over Rational Field defined by
        y^2 = -46656*x^6 + 46656*x^5 - 19440*x^4 + 4320*x^3 - 540*x^2 + 4410*x - 1

        sage: HyperellipticCurve_from_invariants([21, 225/64, 22941/512, 1])
        Traceback (most recent call last):
        ...
        NotImplementedError: Reduction of hyperelliptic curves not yet implemented.
        See issues #14755 and #14756.

    An example over a finite field::

        sage: H = HyperellipticCurve_from_invariants([GF(13)(1), 3, 7, 5]); H
        Hyperelliptic Curve over Finite Field of size 13 defined by ...
        sage: H.igusa_clebsch_invariants()
        (4, 9, 6, 11)

    An example over a number field::

        sage: K = QuadraticField(353, 'a')                                              # needs sage.rings.number_field
        sage: H = HyperellipticCurve_from_invariants([21, 225/64, 22941/512, 1],        # needs sage.rings.number_field
        ....:                                        reduced=false)
        sage: f = K['x'](H.hyperelliptic_polynomials()[0])                              # needs sage.rings.number_field

    If the Mestre Conic defined by the Igusa-Clebsch invariants has no rational
    points, then there exists no hyperelliptic curve over the base field with
    the given invariants.::

        sage: HyperellipticCurve_from_invariants([1,2,3,4])
        Traceback (most recent call last):
        ...
        ValueError: No such curve exists over Rational Field as there are
        no rational points on Projective Conic Curve over Rational Field defined by
        -2572155000*u^2 - 317736000*u*v + 1250755459200*v^2 + 2501510918400*u*w
        + 39276887040*v*w + 2736219686912*w^2

    Mestre's algorithm only works for generic curves of genus two, so another
    algorithm is needed for those curves with extra automorphism. See also
    :issue:`12199`::

        sage: P.<x> = QQ[]
        sage: C = HyperellipticCurve(x^6 + 1)
        sage: i = C.igusa_clebsch_invariants()
        sage: HyperellipticCurve_from_invariants(i)
        Traceback (most recent call last):
        ...
        TypeError: F (=0) must have degree 2


    Igusa-Clebsch invariants also only work over fields of characteristic
    different from 2, 3, and 5, so another algorithm will be needed for fields
    of those characteristics. See also :issue:`12200`::

        sage: P.<x> = GF(3)[]
        sage: HyperellipticCurve(x^6 + x + 1).igusa_clebsch_invariants()
        Traceback (most recent call last):
        ...
        NotImplementedError: Invariants of binary sextics/genus 2 hyperelliptic curves
        not implemented in characteristics 2, 3, and 5
        sage: HyperellipticCurve_from_invariants([GF(5)(1), 1, 0, 1])
        Traceback (most recent call last):
        ...
        ZeroDivisionError: inverse of Mod(0, 5) does not exist

    ALGORITHM:

    This is Mestre's algorithm [Mes1991]_. Our implementation is based on the
    formulae on page 957 of [LY2001]_, cross-referenced with [Wam1999b]_ to
    correct typos.

    First construct Mestre's conic using the :func:`Mestre_conic` function.
    Parametrize the conic if possible.
    Let `f_1, f_2, f_3` be the three coordinates of the parametrization of the
    conic by the projective line, and change them into one variable by letting
    `F_i = f_i(t, 1)`. Note that each `F_i` has degree at most 2.

    Then construct a sextic polynomial
    `f = \sum_{0<=i,j,k<=3}{c_{ijk}*F_i*F_j*F_k}`,
    where `c_{ijk}` are defined as rational functions in the invariants
    (see the source code for detailed formulae for `c_{ijk}`).
    The output is the hyperelliptic curve `y^2 = f`.
    """
    from sage.structure.sequence import Sequence
    i = Sequence(i)
    k = i.universe()
    try:
        k = k.fraction_field()
    except (TypeError, AttributeError, NotImplementedError):
        pass

    MConic, x, y, z = Mestre_conic(i, xyz=True)
    if k.is_finite():
        reduced = False

    t = k['t'].gen()

    if algorithm == 'magma':
        from sage.interfaces.magma import magma
        if MConic.has_rational_point(algorithm='magma'):
            parametrization = magma(MConic).Parametrization().DefiningPolynomials().sage()
            F1, F2, F3 = (p(t, 1) for p in parametrization)
        else:
            raise ValueError(f"No such curve exists over {k} as there are no "
                             f"rational points on {MConic}")
    else:
        if MConic.has_rational_point():
            parametrization = MConic.parametrization(morphism=False)[0]
            F1, F2, F3 = (p(t, 1) for p in parametrization)
        else:
            raise ValueError(f"No such curve exists over {k} as there are no "
                             f"rational points on {MConic}")

    # setting the cijk from Mestre's algorithm
    c111 = 12*x*y - 2*y/3 - 4*z
    c112 = -18*x**3 - 12*x*y - 36*y**2 - 2*z
    c113 = -9*x**3 - 36*x**2*y - 4*x*y - 6*x*z - 18*y**2
    c122 = c113
    c123 = -54*x**4 - 36*x**2*y - 36*x*y**2 - 6*x*z - 4*y**2 - 24*y*z
    c133 = -27*x**4/2 - 72*x**3*y - 6*x**2*y - 9*x**2*z - 39*x*y**2 - \
           36*y**3 - 2*y*z
    c222 = -27*x**4 - 18*x**2*y - 6*x*y**2 - 8*y**2/3 + 2*y*z
    c223 = 9*x**3*y - 27*x**2*z + 6*x*y**2 + 18*y**3 - 8*y*z
    c233 = -81*x**5/2 - 27*x**3*y - 9*x**2*y**2 - 4*x*y**2 + 3*x*y*z - 6*z**2
    c333 = 27*x**4*y/2 - 27*x**3*z/2 + 9*x**2*y**2 + 3*x*y**3 - 6*x*y*z + \
           4*y**3/3 - 10*y**2*z

    # writing out the hyperelliptic curve polynomial
    f = c111*F1**3 + c112*F1**2*F2 + c113*F1**2*F3 + c122*F1*F2**2 + \
        c123*F1*F2*F3 + c133*F1*F3**2 + c222*F2**3 + c223*F2**2*F3 + \
        c233*F2*F3**2 + c333*F3**3

    try:
        f = f * f.denominator()  # clear the denominator
    except (AttributeError, TypeError):
        pass

    if reduced:
        raise NotImplementedError("Reduction of hyperelliptic curves not "
                                  "yet implemented. "
                                  "See issues #14755 and #14756.")

    return HyperellipticCurve(f)


def Mestre_conic(i, xyz=False, names='u,v,w'):
    r"""
    Return the conic equation from Mestre's algorithm given the Igusa-Clebsch
    invariants.

    It has a rational point if and only if a hyperelliptic curve
    corresponding to the invariants exists.

    INPUT:

    - ``i`` -- list or tuple of length 4 containing the four Igusa-Clebsch
      invariants: I2, I4, I6, I10
    - ``xyz`` -- boolean (default: ``False``); if ``True``, the algorithm also
      returns three invariants `x`,`y`,`z` used in Mestre's algorithm
    - ``names`` -- (default: ``'u,v,w'``) the variable names for the conic

    OUTPUT: a Conic object

    EXAMPLES:

    A standard example::

        sage: Mestre_conic([1,2,3,4])
        Projective Conic Curve over Rational Field defined by
        -2572155000*u^2 - 317736000*u*v + 1250755459200*v^2 + 2501510918400*u*w
        + 39276887040*v*w + 2736219686912*w^2

    Note that the algorithm works over number fields as well::

        sage: x = polygen(ZZ, 'x')
        sage: k = NumberField(x^2 - 41, 'a')                                            # needs sage.rings.number_field
        sage: a = k.an_element()                                                        # needs sage.rings.number_field
        sage: Mestre_conic([1, 2 + a, a, 4 + a])                                        # needs sage.rings.number_field
        Projective Conic Curve over Number Field in a with defining polynomial x^2 - 41
         defined by (-801900000*a + 343845000)*u^2 + (855360000*a + 15795864000)*u*v
          + (312292800000*a + 1284808579200)*v^2 + (624585600000*a + 2569617158400)*u*w
          + (15799910400*a + 234573143040)*v*w + (2034199306240*a + 16429854656512)*w^2

    And over finite fields::

        sage: Mestre_conic([GF(7)(10), GF(7)(1), GF(7)(2), GF(7)(3)])
        Projective Conic Curve over Finite Field of size 7
        defined by -2*u*v - v^2 - 2*u*w + 2*v*w - 3*w^2

    An example with ``xyz``::

        sage: Mestre_conic([5,6,7,8], xyz=True)
        (Projective Conic Curve over Rational Field
          defined by -415125000*u^2 + 608040000*u*v + 33065136000*v^2
                     + 66130272000*u*w + 240829440*v*w + 10208835584*w^2,
         232/1125, -1072/16875, 14695616/2109375)

    ALGORITHM:

    The formulas are taken from pages 956 - 957 of [LY2001]_ and based on pages
    321 and 332 of [Mes1991]_.

    See the code or [LY2001]_ for the detailed formulae defining x, y, z and L.
    """
    from sage.structure.sequence import Sequence
    k = Sequence(i).universe()
    try:
        k = k.fraction_field()
    except (TypeError, AttributeError, NotImplementedError):
        pass

    I2, I4, I6, I10 = i

    #Setting x,y,z as in Mestre's algorithm (Using Lauter and Yang's formulas)
    x = 8*(1 + 20*I4/(I2**2))/225
    y = 16*(1 + 80*I4/(I2**2) - 600*I6/(I2**3))/3375
    z = -64*(-10800000*I10/(I2**5) - 9 - 700*I4/(I2**2) + 3600*I6/(I2**3) +
              12400*I4**2/(I2**4) - 48000*I4*I6/(I2**5))/253125

    L = Matrix([[x+6*y     , 6*x**2+2*y         , 2*z                      ],
                [6*x**2+2*y, 2*z                , 9*x**3 + 4*x*y + 6*y**2  ],
                [2*z       , 9*x**3+4*x*y+6*y**2, 6*x**2*y + 2*y**2 + 3*x*z]])

    try:
        L = L*L.denominator()  # clears the denominator
    except (AttributeError, TypeError):
        pass

    u, v, w = PolynomialRing(k, names).gens()
    MConic = Conic(k, L, names)
    if xyz:
        return MConic, x, y, z
    return MConic
