r"""
Points on projective schemes

This module implements scheme morphism for points on weighted projective schemes.

AUTHORS:

- Gareth Ma (2024): initial version, based on unweighted version.
"""

# ****************************************************************************
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2024 Gareth Ma <grhkm21@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy

from sage.categories.integral_domains import IntegralDomains
from sage.rings.fraction_field import FractionField
from sage.schemes.generic.morphism import (
    SchemeMorphism,
    SchemeMorphism_point,
    is_SchemeMorphism,
)
from sage.structure.richcmp import op_EQ, op_NE, richcmp
from sage.structure.sequence import Sequence

# --------------------
# Projective varieties
# --------------------


class SchemeMorphism_point_weighted_projective_ring(SchemeMorphism_point):
    """
    A rational point of projective space over a ring.

    INPUT:

    -  ``X`` -- a homset of a subscheme of an ambient projective space over a ring `K`.

    - ``v`` -- a list or tuple of coordinates in `K`.

    - ``check`` -- boolean (default:``True``). Whether to check the input for consistency.

    EXAMPLES::

        sage: P = ProjectiveSpace(2, ZZ)
        sage: P(2,3,4)
        (2 : 3 : 4)
    """

    def __init__(self, X, v, check=True):
        """
        The Python constructor.

        EXAMPLES::

            TODO
        """
        SchemeMorphism.__init__(self, X)

        if check:
            # check parent
            from sage.schemes.hyperelliptic_curves_smooth_model.weighted_projective_homset import (
                SchemeHomset_points_weighted_projective_ring,
            )
            if not isinstance(X, SchemeHomset_points_weighted_projective_ring):
                raise TypeError(f"ambient space {X} must be a weighted projective space")

            from sage.rings.ring import CommutativeRing
            from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
            d = X.codomain().ambient_space().ngens()
            if isinstance(v, SchemeMorphism) or isinstance(v, EllipticCurvePoint_field):
                v = list(v)
            else:
                try:
                    if isinstance(v.parent(), CommutativeRing):
                        v = [v]
                except AttributeError:
                    pass
            if not isinstance(v, (list, tuple)):
                raise TypeError("argument v (= %s) must be a scheme point, list, or tuple" % str(v))
            if len(v) != d and len(v) != d-1:
                raise TypeError("v (=%s) must have %s components" % (v, d))

            R = X.value_ring()
            v = Sequence(v, R)
            if len(v) == d-1:     # very common special case
                v.append(R.one())

            if R in IntegralDomains():
                # Over integral domains, any tuple with at least one
                # non-zero coordinate is a valid projective point.
                if not any(v):
                    raise ValueError(f"{v} does not define a valid projective "
                                     "point since all entries are zero")
            # Over rings with zero divisors, a more careful check
            # is required: We test whether the coordinates of the
            # point generate the unit ideal. See #31576.
            elif 1 not in R.ideal(v):
                raise ValueError(f"{v} does not define a valid projective point "
                                 "since it is a multiple of a zero divisor")

            X.extended_codomain()._check_satisfies_equations(v)

        self._coords = tuple(v)
        self._normalized = False

    def _repr_(self):
        return "({})".format(" : ".join(map(repr, self._coords)))

    def _richcmp_(self, other, op):
        """
        Test the projective equality of two points.

        INPUT:

        - ``right`` -- a point on projective space

        OUTPUT:

        Boolean

        EXAMPLES::

            TODO
        """
        assert isinstance(other, SchemeMorphism_point)

        if self.codomain() != other.codomain():
            return op == op_NE

        n = len(self._coords)
        if op in [op_EQ, op_NE]:
            b = all(self[i] * other[j] == self[j] * other[i]
                    for i in range(n) for j in range(i + 1, n))
            return b == (op == op_EQ)
        return richcmp(self._coords, other._coords, op)

    def __hash__(self):
        """
        Compute the hash value of this point.

        If the base ring has a fraction field, normalize the point in
        the fraction field and then hash so that equal points have
        equal hash values. If the base ring is not an integral domain,
        return the hash of the parent.

        OUTPUT: Integer.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: hash(P([1, 1])) == hash(P.point([2, 2], False))
            True

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<w> = NumberField(x^2 + 3)
            sage: O = K.maximal_order()
            sage: P.<x,y> = ProjectiveSpace(O, 1)
            sage: hash(P([1 + w, 2])) == hash(P([2, 1 - w]))
            True

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(Zmod(10), 1)
            sage: Q.<x,y> = ProjectiveSpace(Zmod(10), 1)
            sage: hash(P([2, 5])) == hash(Q([2, 5]))
            True
            sage: hash(P([2, 5])) == hash(P([2, 5]))
            True
            sage: hash(P([3, 7])) == hash(P([2, 5]))
            True
        """
        R = self.codomain().base_ring()
        #if there is a fraction field normalize the point so that
        #equal points have equal hash values
        if R in IntegralDomains():
            P = self.change_ring(FractionField(R))
            P.normalize_coordinates()
            return hash(tuple(P))
        #if there is no good way to normalize return
        #a constant value
        return hash(self.codomain())


class SchemeMorphism_point_weighted_projective_field(SchemeMorphism_point_weighted_projective_ring):
    """
    A rational point of projective space over a field.

    INPUT:

    - ``X`` -- a homset of a subscheme of an ambient projective space
       over a field `K`.

    - ``v`` -- a list or tuple of coordinates in `K`.

    - ``check`` -- boolean (default:``True``). Whether to
      check the input for consistency.

    EXAMPLES::

        sage: # needs sage.rings.real_mpfr
        sage: P = ProjectiveSpace(3, RR)
        sage: P(2, 3, 4, 5)
        (0.400000000000000 : 0.600000000000000 : 0.800000000000000 : 1.00000000000000)
    """

    def __init__(self, X, v, check=True):
        """
        The Python constructor.

        See :class:`SchemeMorphism_point_projective_ring` for details.

        This function still normalizes points so that the rightmost non-zero coordinate is 1.
        This is to maintain functionality with current
        implementations of curves in projectives space (plane, conic, elliptic, etc).
        The :class:`SchemeMorphism_point_projective_ring` is for general use.

        EXAMPLES::

            sage: P = ProjectiveSpace(2, QQ)
            sage: P(2, 3/5, 4)
            (1/2 : 3/20 : 1)

        ::

            sage: P = ProjectiveSpace(3, QQ)
            sage: P(0, 0, 0, 0)
            Traceback (most recent call last):
            ...
            ValueError: [0, 0, 0, 0] does not define a valid projective point since all entries are zero

        ::

            sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: X = P.subscheme([x^2 - y*z])
            sage: X([2, 2, 2])
            (1 : 1 : 1)

        ::

            sage: P = ProjectiveSpace(1, GF(7))
            sage: Q = P([2, 1])
            sage: Q[0].parent()
            Finite Field of size 7

        ::

            sage: P = ProjectiveSpace(QQ, 1)
            sage: P.point(Infinity)
            (1 : 0)
            sage: P(infinity)
            (1 : 0)

        ::

            sage: P = ProjectiveSpace(QQ, 2)
            sage: P(infinity)
            Traceback (most recent call last):
            ...
            ValueError: +Infinity not well defined in dimension > 1
            sage: P.point(infinity)
            Traceback (most recent call last):
            ...
            ValueError: +Infinity not well defined in dimension > 1
        """
        SchemeMorphism.__init__(self, X)

        self._normalized = False

        if check:
            from sage.rings.ring import CommutativeRing
            from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
            d = X.codomain().ambient_space().ngens()
            if is_SchemeMorphism(v) or isinstance(v, EllipticCurvePoint_field):
                v = list(v)
            else:
                try:
                    if isinstance(v.parent(), CommutativeRing):
                        v = [v]
                except AttributeError:
                    pass
            if not isinstance(v, (list,tuple)):
                raise TypeError("argument v (= %s) must be a scheme point, list, or tuple" % str(v))
            if len(v) != d and len(v) != d-1:
                raise TypeError("v (=%s) must have %s components" % (v, d))

            R = X.value_ring()
            v = Sequence(v, R)
            if len(v) == d-1:     # very common special case
                v.append(R.one())

            for last in reversed(range(len(v))):
                c = v[last]
                if c.is_one():
                    break
                if c:
                    for j in range(last):
                        v[j] /= c
                    v[last] = R.one()
                    break
            else:
                raise ValueError(f"{v} does not define a valid projective "
                                 "point since all entries are zero")
            self._normalized = True

            X.extended_codomain()._check_satisfies_equations(v)

        self._coords = tuple(v)

    def __hash__(self):
        """
        Computes the hash value of this point.

        OUTPUT: Integer.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: hash(P([1/2, 1])) == hash(P.point([1, 2], False))
            True
        """
        P = copy(self)
        P.normalize_coordinates()
        return hash(tuple(P))

    def normalize_coordinates(self):
        r"""
        Normalizes the point so that the last non-zero coordinate is `1`.

        OUTPUT: None.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5), 2)
            sage: Q = P.point([GF(5)(1), GF(5)(3), GF(5)(0)], False); Q
            (1 : 3 : 0)
            sage: Q.normalize_coordinates(); Q
            (2 : 1 : 0)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme(x^2 - y^2);
            sage: Q = X.point([23, 23, 46], False); Q
            (23 : 23 : 46)
            sage: Q.normalize_coordinates(); Q
            (1/2 : 1/2 : 1)
        """
        if self._normalized:
            return
        for index in reversed(range(len(self._coords))):
            c = self._coords[index]
            if c.is_one():
                break
            if c:
                inv = c.inverse()
                new_coords = [d * inv for d in self._coords[:index]]
                new_coords.append(self.base_ring().one())
                new_coords.extend(self._coords[index+1:])
                self._coords = tuple(new_coords)
                break
        else:
            assert False, 'bug: invalid projective point'
        self._normalized = True


class SchemeMorphism_point_weighted_projective_finite_field(SchemeMorphism_point_weighted_projective_field):

    def __hash__(self):
        r"""
        Returns the integer hash of this point.

        OUTPUT: Integer.

        EXAMPLES: TODO
        """
        p = self.codomain().base_ring().order()
        N = self.codomain().ambient_space().dimension_relative()
        return hash(sum(hash(self[i]) * p**i for i in range(N + 1)))
