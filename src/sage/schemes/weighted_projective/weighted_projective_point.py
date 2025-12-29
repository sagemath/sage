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

from sage.categories.fields import Fields
from sage.rings.fraction_field import FractionField
from sage.misc.misc_c import prod

from sage.schemes.generic.morphism import (SchemeMorphism,
                                           SchemeMorphism_point)
from sage.structure.sequence import Sequence
from sage.structure.richcmp import richcmp, op_EQ, op_NE


# --------------------
# Projective varieties
# --------------------

class SchemeMorphism_point_weighted_projective_ring(SchemeMorphism_point):
    """
    A rational point of weighted projective space over a ring.

    INPUT:

    - ``X`` -- a homset of a subscheme of an ambient weighted projective space over a ring `R`.

    - ``v`` -- a list or tuple of coordinates in `R`.

    - ``check`` -- boolean (default: ``True``); Whether to check the input for consistency.

    EXAMPLES::

        sage: WP = WeightedProjectiveSpace(2, ZZ)
        sage: WP(2, 3, 4)
        (2 : 3 : 4)
    """

    def __init__(self, X, v, check: bool = True):
        """
        The Python constructor.

        EXAMPLES::

            sage: WP = WeightedProjectiveSpace([3, 4, 5], QQ)
            sage: P1 = WP(2, 3, 4); P1
            (2 : 3 : 4)
            sage: P2 = WP(2000, 30000, 400000); P2
            (2000 : 30000 : 400000)
            sage: P1 == P2
            True

        The point constructor normalises coordinates at the last position with
        weight `1`::

            sage: WP = WeightedProjectiveSpace([3, 1, 4, 1], QQ)
            sage: P = WP(8, 3, 16, 2); P
            (1 : 3/2 : 1 : 1)
            sage: P == WP(1, 3 / 2, 1, 1)
            True
        """
        super().__init__(X)

        self._normalized = False

        if check:
            # check parent
            from sage.schemes.weighted_projective.weighted_projective_homset import SchemeHomset_points_weighted_projective_ring
            if not isinstance(X, SchemeHomset_points_weighted_projective_ring):
                raise TypeError(f"ambient space {X} must be a weighted projective space")

            from sage.rings.ring import CommutativeRing
            d = X.codomain().ambient_space().ngens()
            if isinstance(v, SchemeMorphism):
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
            if not R.is_integral_domain():
                raise ValueError("cannot validate point over a ring that is not an integral domain, "
                                 "pass check=False to construct the point")
            v = Sequence(v, R)
            if len(v) == d-1:     # very common special case
                v.append(R.one())

            # (0 : 0 : ... : 0) is not a valid (weighted) projective point
            if not any(v):
                raise ValueError(f"{v} does not define a valid projective "
                                 "point since all entries are zero")

            X.extended_codomain()._check_satisfies_equations(v)

            self._coords = tuple(v)

            # we (try to) normalise coordinates!
            self.normalize_coordinates()

        else:
            self._coords = tuple(v)

    def _repr_(self) -> str:
        return "({})".format(" : ".join(map(repr, self._coords)))

    def _latex_(self) -> str:
        from sage.misc.latex import latex
        return r"\left({}\right)".format(" : ".join(map(latex, self._coords)))

    def _richcmp_(self, other: SchemeMorphism_point, op) -> bool:
        """
        Test the weighted projective equality of two points.

        INPUT:

        - ``other`` -- a point on weighted projective space

        OUTPUT:

        Boolean

        EXAMPLES::

            sage: WP = WeightedProjectiveSpace([3, 4, 5], QQ)
            sage: u = abs(QQ.random_element()) + 1 / 4 # nonzero
            sage: P, Q = WP(2, 3, 4), WP(2 * u^3, 3 * u^4, 4 * u^5)
            sage: P == Q
            True
            sage: P != Q
            False

        ::

            sage: P, Q = WP(2, 3, 4), WP(2, 3, 5)
            sage: P < Q
            True
            sage: P <= Q
            True
            sage: P > Q
            False
            sage: P >= Q
            False

        In the current implementation, points over different base rings
        are never equal (this may change later)::

            sage: WP(2, 3, 4) == WeightedProjectiveSpace([3, 4, 5], ZZ)(2, 3, 4)
            False
        """
        space = self.codomain()
        if space is not other.codomain():
            return op == op_NE

        if op in (op_EQ, op_NE):
            weights = space._weights
            # (other[i] / self[i])^(1 / weight[i]) all equal
            # check weights
            if (weights == other.codomain()._weights) != (op == op_EQ):
                return False

            # check zeros
            b1 = all(c1 == c2
                     for c1, c2 in zip(self._coords, other._coords)
                     if c1 == 0 or c2 == 0)
            if b1 != (op == op_EQ):
                return False

            # check nonzeros
            prod_weights = prod(weights)
            ratio = [(c1 / c2) ** (prod_weights // w)
                     for c1, c2, w in zip(self._coords, other._coords, weights)
                     if c1 != 0 and c2 != 0]
            r0 = ratio[0]
            b2 = all(r == r0 for r in ratio)
            return b2 == (op == op_EQ)

        return richcmp(self._coords, other._coords, op)

    def __hash__(self):
        """
        Compute the hash value of this point.

        We attempt to normalise the coordinates of this point over the field of
        fractions of the base ring. If this is not possible, return the hash of
        the parent. See :meth:`normalize_coordinates` for more details.

        OUTPUT: Integer.

        .. SEEALSO ::

            :meth:`normalize_coordinates`

        EXAMPLES::

            sage: WP = WeightedProjectiveSpace([1, 3, 1], QQ)
            sage: hash(WP(6, 24, 2)) == hash(WP(3, 3, 1))
            True
            sage: WP = WeightedProjectiveSpace([1, 3, 1], ZZ)
            sage: hash(WP(6, 24, 2)) == hash(WP(3, 3, 1))
            True
        """
        R = self.codomain().base_ring()
        P = self.change_ring(FractionField(R))
        P.normalize_coordinates()
        if P._normalized:
            return hash(tuple(P))
        # if there is no good way to normalize return a constant value
        return hash(self.codomain())

    def scale_by(self, t) -> None:
        """
        Scale the coordinates of the point by ``t``.

        A :exc:`TypeError` occurs if the point is not in the
        base_ring of the codomain after scaling.

        INPUT:

        - ``t`` -- a ring element

        OUTPUT: none

        EXAMPLES::

            sage: WP = WeightedProjectiveSpace([3, 4, 5], ZZ)
            sage: P = WP([8, 16, 32]); P
            (8 : 16 : 32)
            sage: P.scale_by(1 / 2); P
            (1 : 1 : 1)
            sage: P.scale_by(1 / 3); P
            Traceback (most recent call last):
            ...
            TypeError: ...
        """
        if t.is_zero():
            raise ValueError("Cannot scale by 0")
        R = self.codomain().base_ring()
        self._coords = tuple([R(u * t**w) for u, w in zip(self._coords, self.codomain()._weights)])
        self._normalized = False

    def normalize_coordinates(self) -> None:
        """
        Normalise coordinates of this weighted projective point if possible.

        Currently, this method checks if (1) the ambient weighted projective
        space is defined over a field and (2) the weight of any index is `1`.
        If so, the last of which is rescaled to `1`.

        EXAMPLES::

            sage: WP = WeightedProjectiveSpace([3, 1, 5], QQ)
            sage: P = WP([8, 16, 32]); P
            (1/512 : 1 : 1/32768)
            sage: P.scale_by(13); P
            (2197/512 : 13 : 371293/32768)
            sage: P.normalize_coordinates(); P
            (1/512 : 1 : 1/32768)

        ::

            sage: WP = WeightedProjectiveSpace([3, 4, 5], ZZ)
            sage: P = WP([8, 16, 32]); P
            (8 : 16 : 32)
            sage: P.normalize_coordinates(); P
            (8 : 16 : 32)
        """
        if self._normalized:
            return

        if self.base_ring() in Fields():
            weights = self.codomain()._weights
            coords = self._coords
            for i in reversed(range(len(coords))):
                w, c = weights[i], coords[i]
                if w.is_one() and not c.is_zero():
                    # we normalise w.r.t this coordinate
                    self.scale_by(~c)
                    self._normalized = True
                    return
