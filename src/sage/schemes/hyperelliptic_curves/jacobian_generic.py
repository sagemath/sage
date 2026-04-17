r"""
Jacobian of a general hyperelliptic curve

AUTHORS:

- David Kohel (2006): initial version
- Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2024): adapt to smooth model
"""

# ****************************************************************************
#       Copyright (C) 2025 Sabrina Kunzweiler <sabrina.kunzweiler@math.u-bordeaux.fr>
#                     2025 Gareth Ma <grhkm21@gmail.com>
#                     2025 Giacomo Pope <giacomopope@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.schemes.hyperelliptic_curves import (
    jacobian_homset_inert,
    jacobian_homset_ramified,
    jacobian_homset_split,
    jacobian_morphism,
)
from sage.schemes.jacobians.abstract_jacobian import Jacobian_generic


class HyperellipticJacobian_generic(Jacobian_generic):
    r"""
    This is the base class for Jacobians of hyperelliptic curves.

    We represent elements of the Jacobian by tuples of the form
    `(u, v : n)`, where

    - `(u,v)` is the Mumford representative of a divisor `P_1 + ... + P_r`,

    - `n` is a non-negative integer

    This tuple represents the equivalence class

    .. MATH::

        [P_1 + ... + P_r + n \cdot \infty_+ + m\cdot \infty_- - D_\infty],

    where  `m = g - \deg(u) - n`, and `\infty_+`, `\infty_-` are the
    points at infinity of the hyperelliptic curve,

    .. MATH::
        D_\infty =
        \lceil g/2 \rceil \infty_+ + \lfloor g/2 \rfloor \infty_-.

    Here, `\infty_- = \infty_+`, if the hyperelliptic curve is ramified.

    Such a representation exists and is unique, unless the genus `g` is odd
    and the curve is inert.

    If the hyperelliptic curve is ramified or inert, then `n` can be deduced
    from `\deg(u)` and `g`. In these cases, `n` is omitted in the description.
    """

    def dimension(self) -> Integer:
        r"""
        Return the dimension of this Jacobian.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x^2, x^4+1); H
            Hyperelliptic Curve over Rational Field defined by y^2 + (x^4 + 1)*y = x^2
            sage: J = Jacobian(H)
            sage: J.dimension()
            3
        """
        return self.curve().genus()

    def is_parent_of(self, element):
        r"""
        Return whether ``self`` is the parent of ``element``.

        TESTS::

            sage: R.<x> = GF(19)[]
            sage: H = HyperellipticCurve(x^5 + x, x^2 + 1)
            sage: J = H.jacobian()
            sage: J.is_parent_of(J.zero())
            True
        """
        from sage.structure.element import parent as get_parent
        p = get_parent(element)
        return p == self or p == self.point_homset()

    def point(self, *mumford, check=True, **kwargs):
        r"""
        Return a point on the Jacobian, given:

        1. No arguments or the integer `0`; return `0 \in J`;

        2. A point `P` on `J = Jac(C)`, return `P`;

        3. A point `P` on the curve `H` such that `J = Jac(H)`;
           return `[P - P_0]`, where `P_0` is the distinguished point of `H`.
           By default, `P_0 = \infty`;

        4. Two points `P, Q` on the curve `H` such that `J = Jac(H)`;
           return `[P - Q]`;

        5. Polynomials `(u, v)` such that `v^2 + hv - f \equiv 0 \pmod u`;
           return `[(u(x), y - v(x))]`.

        .. SEEALSO::

            :mod:`sage.schemes.hyperelliptic_curves.jacobian_homset_generic`.
        """
        try:
            return self.point_homset()(*mumford, check=check)
        except AttributeError:
            raise ValueError("Arguments must determine a valid Mumford divisor.")

    def _point_homset(self, *args, **kwds):
        r"""
        Create the Hom-Set of the Jacobian according to the type of `self`.
        """
        # TODO: make a constructor for this??
        H = self.curve()
        if H.is_ramified():
            return jacobian_homset_ramified.HyperellipticJacobianHomsetRamified(
                *args, **kwds
            )
        if H.is_split():
            return jacobian_homset_split.HyperellipticJacobianHomsetSplit(*args, **kwds)
        return jacobian_homset_inert.HyperellipticJacobianHomsetInert(*args, **kwds)

    def _point(self, *args, **kwds):
        H = self.curve()
        if H.is_ramified():
            return jacobian_morphism.MumfordDivisorClassFieldRamified(*args, **kwds)
        if H.is_split():
            return jacobian_morphism.MumfordDivisorClassFieldSplit(*args, **kwds)
        return jacobian_morphism.MumfordDivisorClassFieldInert(*args, **kwds)

    @cached_method
    def order(self):
        r"""
        Compute the order of the Jacobian.

        EXAMPLES::

            sage: R.<x> = GF(3663031)[]
            sage: HyperellipticCurve(x^5 + 1758294*x^4 + 1908793*x^3 + 3033920*x^2 + 3445698*x + 3020661).jacobian().cardinality()
            13403849798842

        .. SEEALSO::

            :meth:`sage.schemes.hyperelliptic_curves.jacobian_homset_generic.order`.
        """
        return self.point_homset().order()

    cardinality = order

    def count_points(self, *args, **kwds):
        r"""
        .. SEEALSO::

            :meth:`sage.schemes.hyperelliptic_curves.jacobian_homset_generic.count_points`.
        """
        return self.point_homset().count_points(*args, **kwds)

    def lift_u(self, *args, **kwds):
        r"""
        Return one or all points with given `u`-coordinate.

        .. SEEALSO::

            :meth:`sage.schemes.hyperelliptic_curves.jacobian_homset_generic.lift_u`.
        """
        return self.point_homset().lift_u(*args, **kwds)

    def random_element(self, *args, **kwds):
        r"""
        Return a random element of the Jacobian.

        .. SEEALSO::

            :meth:`sage.schemes.hyperelliptic_curves.jacobian_homset_generic.random_element`.
        """
        return self.point_homset().random_element(*args, **kwds)

    def some_elements(self) -> list[jacobian_morphism.MumfordDivisorClassField]:
        return [self.zero()] + [self.random_element(), self.random_element(), self.random_element()]

    def points(self, *args, **kwds):
        r"""
        Return all points on the Jacobian.

        .. SEEALSO::

            :meth:`sage.schemes.hyperelliptic_curves.jacobian_homset_generic.points`.
        """

        return self.point_homset().points(*args, **kwds)

    def list(self):
        r"""
        Return all rational elements of the Jacobian.

        .. SEEALSO::

            :meth:`sage.schemes.hyperelliptic_curves.jacobian_homset_generic.points`.
        """

        return self.point_homset().points()

    def __iter__(self):
        r"""
        Return an iterator over the elements of the Jacobian.
        """
        yield from self.point_homset().points()

    rational_points = points

    ####################################################################
    # Some properties of geometric Endomorphism ring and algebra
    ####################################################################

    @lazy_attribute
    def _have_established_geometrically_trivial(self):
        r"""
        Initialize the flag which determines whether or not we have
        already established if the geometric endomorphism ring is
        trivial.

        This is related to the warning at the top of the
        ``jacobian_endomorphism_utils.py`` module.

        INPUT:

        - ``self`` -- the Jacobian

        OUTPUT: the boolean ``False``; this will be updated by other methods

        EXAMPLES:

        This is LMFDB curve 262144.d.524288.2::

            sage: R.<x> = QQ[]
            sage: f = x^5 + x^4 + 4*x^3 + 8*x^2 + 5*x + 1
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: J._have_established_geometrically_trivial
            False
        """
        return False

    @lazy_attribute
    def _have_established_geometrically_field(self):
        r"""
        Initialize the flag which determines whether or not we have
        already established if the geometric endomorphism ring is
        trivial.

        This is related to the warning at the top of the
        ``jacobian_endomorphism_utils.py`` module.

        INPUT:

        - ``self`` -- the Jacobian

        OUTPUT: the boolean ``False``; this will be updated by other methods

        EXAMPLES:

        This is LMFDB curve 262144.d.524288.2::

            sage: R.<x> = QQ[]
            sage: f = x^5 + x^4 + 4*x^3 + 8*x^2 + 5*x + 1
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: J._have_established_geometrically_field
            False
        """
        return False

    def geometric_endomorphism_algebra_is_field(self, B=200, proof=False):
        r"""
        Return whether the geometric endomorphism algebra is a field.

        This implies that the Jacobian of the curve is geometrically
        simple. It is based on Algorithm 4.10 from [Lom2019]_

        INPUT:

        - ``B`` -- (default: 200) the bound which appears in the statement of
          the algorithm from [Lom2019]_

        - ``proof`` -- boolean (default: ``False``); whether or not to insist
          on a provably correct answer. This is related to the warning in the
          docstring of this module: if this function returns ``False``, then
          strictly speaking this has not been proven to be ``False`` until one
          has exhibited a non-trivial endomorphism, which these methods are not
          designed to carry out. If one is convinced that this method should
          return ``True``, but it is returning ``False``, then this can be
          exhibited by increasing `B`.

        OUTPUT:

        Boolean indicating whether or not the geometric endomorphism
        algebra is a field.

        EXAMPLES:

        This is LMFDB curve 262144.d.524288.2 which has QM. Although its
        Jacobian is geometrically simple, the geometric endomorphism algebra
        is not a field::

            sage: R.<x> = QQ[]
            sage: f = x^5 + x^4 + 4*x^3 + 8*x^2 + 5*x + 1
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: J.geometric_endomorphism_algebra_is_field()
            False

        This is LMFDB curve 50000.a.200000.1::

            sage: f = 8*x^5 + 1
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: J.geometric_endomorphism_algebra_is_field()
            True
        """
        from sage.interfaces.genus2reduction import genus2reduction

        from .jacobian_endomorphism_utils import get_is_geom_field

        if self._have_established_geometrically_field:
            return True
        C = self.curve()
        if C.genus() != 2:
            raise NotImplementedError(
                "Current implementation requires the curve to be of genus 2"
            )
        if C.base_ring() != QQ:
            raise NotImplementedError(
                "Current implementation requires the curve to be defined over the rationals"
            )
        f, h = C.hyperelliptic_polynomials()
        if h != 0:
            raise NotImplementedError(
                "Current implementation requires the curve to be in the form y^2 = f(x)"
            )
        red_data = genus2reduction(0, f)
        cond_C = red_data.conductor  # WARNING: this is only the prime_to_2 conductor.
        bad_primes = cond_C.prime_divisors()
        self._bad_primes = bad_primes

        is_abs_simp, is_def_geom_trivial = get_is_geom_field(f, C, bad_primes, B)

        if is_def_geom_trivial:
            self._have_established_geometrically_trivial = True
        if is_abs_simp:
            self._have_established_geometrically_field = True
            return True
        if proof:
            raise NotImplementedError(
                "Rigorous computation of lower bounds of endomorphism algebras has not yet been implemented."
            )
        return False

    def geometric_endomorphism_ring_is_ZZ(self, B=200, proof=False):
        r"""
        Return whether the geometric endomorphism ring of ``self`` is the
        integer ring `\ZZ`.

        INPUT:

        - ``B`` -- (default: 200) the bound which appears in the statement of
          the algorithm from [Lom2019]_

        - ``proof`` -- boolean (default: ``False``); whether or not to insist
          on a provably correct answer. This is related to the warning in the
          module docstring of `jacobian_endomorphisms.py`: if this function
          returns ``False``, then strictly speaking this has not been proven to
          be ``False`` until one has exhibited a non-trivial endomorphism,
          which the methods in that module are not designed to carry out. If
          one is convinced that this method should return ``True``, but it is
          returning ``False``, then this can be exhibited by increasing `B`.

        OUTPUT:

        Boolean indicating whether or not the geometric endomorphism
        ring is isomorphic to the integer ring.

        EXAMPLES:

        This is LMFDB curve 603.a.603.2::

            sage: R.<x> = QQ[]
            sage: f = 4*x^5 + x^4 - 4*x^3 + 2*x^2 + 4*x + 1
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: J.geometric_endomorphism_ring_is_ZZ()
            True

        This is LMFDB curve 1152.a.147456.1 whose geometric endomorphism ring
        is isomorphic to the group of 2x2 matrices over `\QQ`::

            sage: f = x^6 - 2*x^4 + 2*x^2 - 1
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: J.geometric_endomorphism_ring_is_ZZ()
            False

        This is LMFDB curve 20736.k.373248.1 whose geometric endomorphism ring
        is isomorphic to the group of 2x2 matrices over a CM field::

            sage: f = x^6 + 8
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: J.geometric_endomorphism_ring_is_ZZ()
            False

        This is LMFDB curve 708.a.181248.1::

            sage: R.<x> = QQ[]
            sage: f = -3*x^6 - 16*x^5 + 36*x^4 + 194*x^3 - 164*x^2 - 392*x - 143
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: J.geometric_endomorphism_ring_is_ZZ()
            True

        This is LMFDB curve 10609.a.10609.1 whose geometric endomorphism ring
        is an order in a real quadratic field::

            sage: f = x^6 + 2*x^4 + 2*x^3 + 5*x^2 + 6*x + 1
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: J.geometric_endomorphism_ring_is_ZZ()
            False

        This is LMFDB curve 160000.c.800000.1 whose geometric endomorphism ring
        is an order in a CM field::

            sage: f = x^5 - 1
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: J.geometric_endomorphism_ring_is_ZZ()
            False

        This is LMFDB curve 262144.d.524288.2 whose geometric endomorphism ring
        is an order in a quaternion algebra::

            sage: f = x^5 + x^4 + 4*x^3 + 8*x^2 + 5*x + 1
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: J.geometric_endomorphism_ring_is_ZZ()
            False

        This is LMFDB curve 578.a.2312.1 whose geometric endomorphism ring
        is `\QQ \times \QQ`::

            sage: f = 4*x^5 - 7*x^4 + 10*x^3 - 7*x^2 + 4*x
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: J.geometric_endomorphism_ring_is_ZZ()
            False
        """
        from .jacobian_endomorphism_utils import is_geom_trivial_when_field

        if self._have_established_geometrically_trivial:
            return True
        is_abs_simple = self.geometric_endomorphism_algebra_is_field(B=B, proof=proof)
        if self._have_established_geometrically_trivial:
            return True
        if is_abs_simple and is_geom_trivial_when_field(self.curve(), self._bad_primes):
            return True
        if proof:
            raise NotImplementedError(
                "Rigorous computation of lower bounds of endomorphism rings has not yet been implemented."
            )
        return False
