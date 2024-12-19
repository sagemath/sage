r"""
    Jacobians of hyperelliptic curves

    In Sage, a hyperelliptic curve `H` of genus `g` is always
    specified by an (affine) equation in Weierstrass form

    .. MATH::

        H : y^2 + h(x) y = f(x),

    for some polynomials `h` and `f`.

    Elements of the Jacobian of such curves can be identified
    with equivalence classes of divisors. In the following, we
    denote by `\infty_+`, \infty_-` the points at infinity of `H`,
    and we set

    ..MATH::
        D_\infty =
        \lceil g/2 \rceil \infty_+ + \lfloor g/2 \rfloor \infty_-.

    Here, `\infty_- = \infty_+`, if `H` is ramified.
    Unless the genus `g` is odd and `H` is inert, the divisor
    `D_\infty` is rational. In these cases, any element on
    the Jacobian admits a unique representative of the form

    ..MATH::

        [P_1 + ... + P_r + n \cdot \infty_+ + m\cdot \infty_- - D_\infty],

    with `n` and `m` non-negative integers and `P_1 + ... + P_r`
    an affine and reduced  divisor on `H`.

    This module implements the arithemtic for Jacobians of
    hyperelliptic curves. Elements of the Jacobian are represented
    by tuples of the form `(u, v : n)`, where
    - (u,v) is the Mumford representative of the divisor `P_1 + ... + P_r`,
    - n is the coefficient of `\infty_+`

    We note that `m = g - \deg(u) - n` and is therefore omitted in
    the description. Similarly, if `H` ramified or inert,
    then `n` can be deduced from `\deg(u)` and `g`. In these cases,
    `n` is omitted in the description as well.


    EXAMPLES:

    We construct the Jacobian of a hyperelliptic curve with affine equation
    `y^2 + (x^3 + x + 1) y  = 2*x^5 + 4*x^4 + x^3 - x` over the rationals.
    This curve has two points at infinity::

        sage: R.<x> = QQ[]
        sage: H = HyperellipticCurveSmoothModel(2*x^5 + 4*x^4 + x^3 - x, x^3 + x + 1)
        sage: J = Jacobian(H); J
        Jacobian of Hyperelliptic Curve over Rational Field defined by y^2 + (x^3 + x + 1)*y = 2*x^5 + 4*x^4 + x^3 - x

    The points `P = (0, 0)` and `Q = (-1, -1)` are on `H`. We construct the
    element `D_1 = [P - Q] = [P + (-Q) - D_\infty`] on the Jacobian::

        sage: P = H.point([0, 0])
        sage: Q = H.point([-1, -1])
        sage: D1 = J(P,Q); D1
        (x^2 + x, -2*x : 0)

    Elements of the Jacobian can also be constructed by directly providing
    the Mumford representation::

        sage: D1 == J(x^2 + x, -2*x, 0)
        True

    We can also embed single points into the Jacobian. Below we construct
    `D_2 = [P - P_0]`, where `P_0` is the distinguished point of `H`
    (by default one of the points at infinity)::

        sage: D2 = J(P); D2
        (x, 0 : 0)
        sage: P0 = H.distinguished_point(); P0
        (1 : 0 : 0)
        sage: D2 == J(P, P0)
        True

    We may add elements, or multiply by integers::

        sage: 2*D1
        (x, -1 : 1)
        sage: D1 + D2
        (x^2 + x, -1 : 0)
        sage: -D2
        (x, -1 : 1)

    Note that the neutral element is given by `[D_\infty - D_\infty]`,
    in particular `n = 1`::

        sage: J.zero()
        (1, 0 : 1)

    There are two more elements of the Jacobian that are only supported
    at infinity: `[\infty_+ - \infty_-]` and `[\infty_- - \infty_+]`::

        sage: [P_plus, P_minus] = H.points_at_infinity()
        sage: P_plus == P0
        True
        sage: J(P_plus,P_minus)
        (1, 0 : 2)
        sage: J(P_minus, P_plus)
        (1, 0 : 0)

    Now, we consider the Jacobian of a hyperelliptic curve with only one
    point at infinity, defined over a finite field::

        sage: K = FiniteField(7)
        sage: R.<x> = K[]
        sage: H = HyperellipticCurveSmoothModel(x^7 + 3*x + 2)
        sage: J = Jacobian(H); J
        Jacobian of Hyperelliptic Curve over Finite Field of size 7 defined by y^2 = x^7 + 3*x + 2

    Elements on the Jacobian can be constructed as before. But the value
    `n` is not used here, since there is only one point at infinity::

        sage: P = H.point([3, 0])
        sage: Q = H.point([5, 1])
        sage: D1 = J(P,Q); D1
        (x^2 + 6*x + 1, 3*x + 5)
        sage: D2 = J(x^3 + 3*x^2 + 4*x + 3, 2*x^2 + 4*x)
        sage: D1 + D2
        (x^3 + 2, 4)

    Over finite fields, we may also construct random elements and
    compute the order of the Jacobian::

        sage: J.random_element() #random
        (x^3 + x^2 + 4*x + 5, 3*x^2 + 3*x)
        sage: J.order()
        344

    Note that arithmetic on the Jacobian is not implemented if the
    underlying hyperelliptic curve is inert (i.e. has no points at
    infinity) and the genus is odd::

        sage: R.<x> = GF(13)[]
        sage: H = HyperellipticCurveSmoothModel(x^8+1,x^4+1)
        sage: J = Jacobian(H); J
        Jacobian of Hyperelliptic Curve over Finite Field of size 13 defined by y^2 + (x^4 + 1)*y = x^8 + 1


    TODO:
    - finish example for the inert case.

    AUTHORS:

    - David Kohel (2006): initial version
    - Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2024): adapt to smooth model
"""

# ****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#                2024 Sabrina Kunzweiler, Gareth Ma, Giacomo Pope
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.schemes.hyperelliptic_curves_smooth_model import (
    jacobian_homset_inert,
    jacobian_homset_ramified,
    jacobian_homset_split,
    jacobian_morphism,
)
from sage.schemes.jacobians.abstract_jacobian import Jacobian_generic


class HyperellipticJacobian_generic(Jacobian_generic):
    """
    This is the base class for Jacobians of hyperelliptic curves.
    """

    def dimension(self):
        """
        Return the dimension of this Jacobian.
        """
        return Integer(self.curve().genus())

    def point(self, *mumford, check=True, **kwargs):
        try:
            return self.point_homset()(*mumford, check=check)
        except AttributeError:
            raise ValueError("Arguments must determine a valid Mumford divisor.")

    def _point_homset(self, *args, **kwds):
        # TODO: make a constructor for this??
        H = self.curve()
        if H.is_ramified():
            return jacobian_homset_ramified.HyperellipticJacobianHomsetRamified(
                *args, **kwds
            )
        elif H.is_split():
            return jacobian_homset_split.HyperellipticJacobianHomsetSplit(*args, **kwds)
        return jacobian_homset_inert.HyperellipticJacobianHomsetInert(*args, **kwds)

    def _point(self, *args, **kwds):
        H = self.curve()
        if H.is_ramified():
            return jacobian_morphism.MumfordDivisorClassFieldRamified(*args, **kwds)
        elif H.is_split():
            return jacobian_morphism.MumfordDivisorClassFieldSplit(*args, **kwds)
        return jacobian_morphism.MumfordDivisorClassFieldInert(*args, **kwds)

    @cached_method
    def order(self):
        return self.point_homset().order()

    def count_points(self, *args, **kwds):
        return self.point_homset().count_points(*args, **kwds)

    def lift_u(self, *args, **kwds):
        return self.point_homset().lift_u(*args, **kwds)

    def random_element(self, *args, **kwds):
        return self.point_homset().random_element(*args, **kwds)

    def points(self, *args, **kwds):
        return self.point_homset().points(*args, **kwds)

    def list(self):
        return self.point_homset().points()

    def __iter__(self):
        yield from self.list()

    rational_points = points
