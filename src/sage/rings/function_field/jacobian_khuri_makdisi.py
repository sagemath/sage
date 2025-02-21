r"""
Jacobians in Khuri-Makdisi model

This module implements Jacobian arithmetic by Khuri-Makdisi's algorithms
[Khu2004]_ based on divisor representation by linear spaces.

Jacobian
--------

There are three models for Jacobian arithmetic by Khuri-Makdisi's algorithms.
For each of the models, one should provide a base divisor satisfying certain
degree condition. The following lists the names of the three models and the
corresponding conditions on base divisors. Let `g` be the genus of the function
field.

- ``km_large``: large model; requires an effective divisor of degree at least `2g + 1`

- ``km_medium``: medium model; requires an effective divisor of degree at least `2g + 1`

- ``km_small``: small model; requires an effective divisor of degree at least `g + 1`

To create a Jacobian in this model, specify ``'km_[large|medium|small]'`` as ``model`` and
provide a base divisor satisfying the degree condition.

EXAMPLES:

We construct a function field (of a projective curve) over a finite field::

    sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
    sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
    sage: C.geometric_genus()
    1
    sage: H = C.function(y/x).divisor_of_poles()
    sage: H.degree()
    3

Now we use `H` as base divisor for the large and medium models::

    sage: J_large = C.jacobian(model='km_large', base_div=H)
    sage: J_large
    Jacobian of Projective Plane Curve over Finite Field of size 29
     defined by x^3 - y^2*z + 5*z^3 (Khuri-Makdisi large model)
    sage: J_medium = C.jacobian(model='km_medium', base_div=H)
    sage: J_medium
    Jacobian of Projective Plane Curve over Finite Field of size 29
     defined by x^3 - y^2*z + 5*z^3 (Khuri-Makdisi medium model)

and for the small model, we construct an effective divisor of degree 2::

    sage: B = sum(H.support()[:2])
    sage: B.degree()
    2
    sage: J_small = C.jacobian(model='km_small', base_div=B)
    sage: J_small
    Jacobian of Projective Plane Curve over Finite Field of size 29
     defined by x^3 - y^2*z + 5*z^3 (Khuri-Makdisi small model)

Group of rational points
------------------------

The group of rational points of a Jacobian is created from the Jacobian. A
point of the Jacobian group is represented by a divisor `D - B` where `D` is
an effective divisor of the same degree with the base divisor `B`. The
divisor `D` in turn is determined by a linear subspace of the Riemann-Roch
space associated with certain multiple of `B` (depending on the model). This
allows representing points of Jacobian as matrices once we fix a basis of the
Riemann-Roch space.

EXAMPLES::

    sage: # long time
    sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
    sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
    sage: F = C.function_field()
    sage: H = C.function(y/x).divisor_of_poles()
    sage: J = C.jacobian(model='km_large', base_div=H)
    sage: G = J.group()
    sage: D = C([0,1,0]).place()
    sage: P1 = C([-1,2,1]).place()
    sage: P2 = C([3,7,1]).place()
    sage: p1 = G.point(P1 - D)
    sage: p2 = G.point(P2 - D)
    sage: p1
    Point of Jacobian determined by
    [ 1  0  0  0  0  0  0 12 15]
    [ 0  1  0  0  0  0  0  0 13]
    [ 0  0  1  0  0  0  0  0  2]
    [ 0  0  0  1  0  0  0  0 16]
    [ 0  0  0  0  0  1  0  0 15]
    [ 0  0  0  0  0  0  1  0  1]
    sage: p2
    Point of Jacobian determined by
    [ 1  0  0  0  0  0  0 12  5]
    [ 0  1  0  0  0  0  0  0  2]
    [ 0  0  1  0  0  0  0  0 13]
    [ 0  0  0  1  0  0  0  0  8]
    [ 0  0  0  0  0  1  0  0 10]
    [ 0  0  0  0  0  0  1  0 14]
    sage: p1 + p2
    Point of Jacobian determined by
    [ 1  0  0  0  0 16  0  5  3]
    [ 0  1  0  0  0  6  0  8 16]
    [ 0  0  1  0  0 15  0  3 10]
    [ 0  0  0  1  0  3  0  0  0]
    [ 0  0  0  0  1 12  0 16  8]
    [ 0  0  0  0  0  0  1  3  0]

AUTHORS:

- Kwankyu Lee (2022-01-24): initial version
"""

# ****************************************************************************
#       Copyright (C) 2022 Kwankyu Lee <ekwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.richcmp import op_EQ, richcmp

from sage.categories.map import Map
from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from sage.categories.homset import Hom

from sage.matrix.constructor import matrix

from sage.combinat.integer_vector_weighted import WeightedIntegerVectors

from .place import FunctionFieldPlace
from .divisor import FunctionFieldDivisor

from .jacobian_base import (Jacobian_base,
                            JacobianGroup_base,
                            JacobianGroup_finite_field_base,
                            JacobianPoint_base,
                            JacobianPoint_finite_field_base)


class JacobianPoint(JacobianPoint_base):
    """
    Points of a Jacobian group.

    INPUT:

    - ``parent`` -- Jacobian group

    - ``w`` -- matrix

    EXAMPLES::

        sage: # long time
        sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
        sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
        sage: b = C([0,1,0]).place()
        sage: h = C.function(y/x).divisor_of_poles()
        sage: J = C.jacobian(model='km_large', base_div=h)
        sage: G = J.group()
        sage: pl = C([3,2,1]).place()
        sage: G.point(pl - b)
        Point of Jacobian determined by
        [1 0 0 0 0 0 0 2 3]
        [0 1 0 0 0 0 0 0 3]
        [0 0 1 0 0 0 0 0 1]
        [0 0 0 1 0 0 0 0 5]
        [0 0 0 0 0 1 0 0 5]
        [0 0 0 0 0 0 1 0 4]
    """
    def __init__(self, parent, w):
        """
        Initialize.

        TESTS::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: pl = C([3,2,1]).place()
            sage: p = G.point(pl - b)
            sage: TestSuite(p).run(skip=['_test_category','_test_pickling'])
        """
        super().__init__(parent)
        w.set_immutable()
        self._w = w

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: G.zero()
            Point of Jacobian determined by
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 1 0 0 0]
            [0 0 0 0 0 0 0 1 0]
        """
        return f'Point of Jacobian determined by \n{self._w}'

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: F = C.function_field()
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: zero = G.zero()
            sage: {zero: 1}
            {Point of Jacobian determined by
             [1 0 0 0 0 0 0 0 0]
             [0 1 0 0 0 0 0 0 0]
             [0 0 1 0 0 0 0 0 0]
             [0 0 0 0 1 0 0 0 0]
             [0 0 0 0 0 1 0 0 0]
             [0 0 0 0 0 0 0 1 0]: 1}
        """
        return hash(self._w)

    def _richcmp_(self, other, op):
        """
        Compare ``self`` with ``other`` with respect to operator ``op``.

        EXAMPLES::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: b = C([0,1,0]).place()
            sage: pl1 = C([-1,2,1]).place()
            sage: pl2 = C([3,2,1]).place()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: p1 = G.point(pl1 - b)
            sage: p2 = G.point(pl2 - b)
            sage: p1 == p1
            True
            sage: p1 != p2
            True
            sage: p1 > p1
            False
            sage: p1 > p2
            True
            sage: p1 < p2
            False
        """
        if op is op_EQ:
            km = self.parent()._km
            return km.equal(self._w, other._w)
        else:
            return richcmp(self._w, other._w, op)

    def _add_(self, other):
        """
        Return the sum of ``self`` and ``other``.

        EXAMPLES::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: b = C([0,1,0]).place()
            sage: pl1 = C([-1,2,1]).place()
            sage: pl2 = C([3,2,1]).place()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: p1 = G.point(pl1 - b)
            sage: p2 = G.point(pl2 - b)
            sage: p1 + p2
            Point of Jacobian determined by
            [1 0 0 0 0 0 3 1 1]
            [0 1 0 0 0 0 3 0 2]
            [0 0 1 0 0 0 0 3 0]
            [0 0 0 1 0 0 0 0 3]
            [0 0 0 0 1 0 4 0 3]
            [0 0 0 0 0 1 6 5 2]
            sage: p1 + p2 == p2 + p1
            True
        """
        G = self.parent()
        km = G._km
        return G.element_class(self.parent(), km.add(self._w, other._w))

    def _neg_(self):
        """
        Return the negative of this point.

        EXAMPLES::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: F = C.function_field()
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: b = F.get_place(1)
            sage: p = C([-1,2,1]).place()
            sage: pt = G.point(p - b)
            sage: -pt
            Point of Jacobian determined by
            [1 0 0 0 0 0 1 6 0]
            [0 1 0 0 0 0 2 6 3]
            [0 0 1 0 0 0 4 1 1]
            [0 0 0 1 0 0 1 4 1]
            [0 0 0 0 1 0 5 2 6]
            [0 0 0 0 0 1 3 1 3]
            sage: -(-pt) == pt
            True
        """
        G = self.parent()
        km = G._km
        return G.element_class(self.parent(), km.negate(self._w))

    def _rmul_(self, n):
        """
        Return the ``n``-th multiple of this point.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: b = C([0,1,0]).place()
            sage: pl = C([-1,2,1]).place()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: p = G.point(pl - b)
            sage: 10*(10*p) == 100*p
            True
        """
        return self.multiple(n)

    def multiple(self, n):
        """
        Return the ``n``-th multiple of this point.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: b = C([0,1,0]).place()
            sage: pl = C([-1,2,1]).place()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: p = G.point(pl - b)
            sage: p.multiple(100)
            Point of Jacobian determined by
            [1 0 0 0 0 2 0 1 1]
            [0 1 0 0 0 5 0 1 6]
            [0 0 1 0 0 2 0 6 3]
            [0 0 0 1 0 1 0 0 0]
            [0 0 0 0 1 5 0 1 4]
            [0 0 0 0 0 0 1 1 0]
        """
        G = self.parent()
        km = G._km
        return G.element_class(self.parent(), km.multiple(self._w, n))

    def addflip(self, other):
        """
        Return the addflip of this and ``other`` point.

        The addflip of two points is by definition the negative of the sum of
        the points. This operation is faster than addition in Jacobian
        arithmetic by Khuri-Makdisi algorithms.

        EXAMPLES::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: b = C([0,1,0]).place()
            sage: pl1 = C([-1,2,1]).place()
            sage: pl2 = C([3,2,1]).place()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: p1 = G.point(pl1 - b)
            sage: p2 = G.point(pl2 - b)
            sage: p1.addflip(p2)
            Point of Jacobian determined by
            [1 0 0 0 0 0 0 2 6]
            [0 1 0 0 0 0 0 0 3]
            [0 0 1 0 0 0 0 0 4]
            [0 0 0 1 0 0 0 0 3]
            [0 0 0 0 0 1 0 0 5]
            [0 0 0 0 0 0 1 0 2]
            sage: _ == -(p1 + p2)
            True
        """
        G = self.parent()
        km = G._km
        return G.element_class(self.parent(), km.addflip(self._w, other._w))

    def defining_matrix(self):
        """
        Return the matrix whose row span determines the effective divisor
        representing this point.

        EXAMPLES::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: b = C([0,1,0]).place()
            sage: pl = C([-1,2,1]).place()
            sage: p = G.point(pl - b)
            sage: p.defining_matrix()
            [1 0 0 0 0 0 0 2 5]
            [0 1 0 0 0 0 0 0 3]
            [0 0 1 0 0 0 0 0 2]
            [0 0 0 1 0 0 0 0 6]
            [0 0 0 0 0 1 0 0 5]
            [0 0 0 0 0 0 1 0 1]
        """
        return self._w

    def divisor(self):
        """
        Return the divisor representing this point.

        EXAMPLES::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: F = C.function_field()
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: b = F.get_place(1)
            sage: p = C([-1,2,1]).place()
            sage: pt = G.point(p - b)
            sage: G.point(pt.divisor()) == pt
            True

        ALGORITHM: Lemma 2.1 of [Khu2004]_.
        """
        G = self.parent()
        F = G._function_field
        data = [G._from_L(f).divisor() for f in self._w.rows()]
        supp = set()
        for d in data:
            supp.update(d.support())
        supp = list(supp)
        d = F.divisor_group().zero()
        for p in supp:
            d += min(d.valuation(p) for d in data) * p.divisor()
        return d + G._div_L - G._base_div


class JacobianPoint_finite_field(JacobianPoint, JacobianPoint_finite_field_base):
    pass


class JacobianGroupEmbedding(Map):
    """
    Embeddings between Jacobian groups.

    INPUT:

    - ``base_group`` -- Jacobian group over a base field

    - ``extension_group`` -- Jacobian group over an extension field

    EXAMPLES::

        sage: # long time
        sage: k = GF(5)
        sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
        sage: C = Curve(x^3 + z^3 - y^2*z, P2)
        sage: h = C.function(y/x).divisor_of_poles()
        sage: J = C.jacobian(model='km_large', base_div=h)
        sage: G1 = J.group()
        sage: K = k.extension(2)
        sage: G2 = J.group(K)
        sage: G2.coerce_map_from(G1)
        Jacobian group embedding map:
          From: Group of rational points of Jacobian
                over Finite Field of size 5 (Khuri-Makdisi large model)
          To:   Group of rational points of Jacobian
                over Finite Field in z2 of size 5^2 (Khuri-Makdisi large model)
    """
    def __init__(self, base_group, extension_group):
        """
        Initialize.

        TESTS::

            sage: # long time
            sage: k = GF(5)
            sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
            sage: C = Curve(x^3 + z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G1 = J.group()
            sage: K = k.extension(2)
            sage: G2 = J.group(K)
            sage: map = G2.coerce_map_from(G1)
            sage: TestSuite(map).run(skip=['_test_category', '_test_pickling'])
        """
        F_ext = extension_group._function_field
        K_ext = F_ext.constant_base_field()

        self._K_ext = K_ext

        Map.__init__(self, Hom(base_group, extension_group, CommutativeAdditiveGroups()))

    def _repr_type(self):
        """
        Return string representation of ``self``.

        TESTS::

            sage: # long time
            sage: k = GF(5)
            sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
            sage: C = Curve(x^3 + z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G1 = J.group()
            sage: K = k.extension(2)
            sage: G2 = J.group(K)
            sage: G2.coerce_map_from(G1)  # indirect doctest
            Jacobian group embedding map:
              From: Group of rational points of Jacobian
                    over Finite Field of size 5 (Khuri-Makdisi large model)
              To:   Group of rational points of Jacobian
                    over Finite Field in z2 of size 5^2 (Khuri-Makdisi large model)
        """
        return 'Jacobian group embedding'

    def _call_(self, x):
        """
        Conorm map from F to F_ext.

        TESTS::

            sage: # long time
            sage: k = GF(5)
            sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
            sage: C = Curve(x^3 + z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G1 = J.group()
            sage: K = k.extension(2)
            sage: G2 = J.group(K)
            sage: m = G2.coerce_map_from(G1)
            sage: m(G1.zero()) == G2.zero()
            True
        """
        w = (x._w).change_ring(self._K_ext)

        return self.codomain().element_class(self.codomain(), w)


class JacobianGroup(UniqueRepresentation, JacobianGroup_base):
    """
    Groups of rational points of a Jacobian.

    INPUT:

    - ``parent`` -- a Jacobian

    - ``function_field`` -- a function field

    - ``base_div`` -- an effective divisor of the function field

    EXAMPLES::

        sage: # long time
        sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
        sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
        sage: h = C.function(y/x).divisor_of_poles()
        sage: J = C.jacobian(model='km_large', base_div=h)
        sage: J.group()
        Group of rational points of Jacobian
         over Finite Field of size 7 (Khuri-Makdisi large model)
    """
    Element = JacobianPoint
    _embedding_map_class = JacobianGroupEmbedding

    def __init__(self, parent, function_field, base_div):
        """
        Initialize.

        TESTS::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
        """
        super().__init__(parent, function_field, base_div)

        D0 = base_div

        self._base_div_degree = base_div.degree()
        self._V_cache = 10*[None]

        V_cache = self._V_cache

        def V(n):
            if n in V_cache:
                return V_cache[n]

            Vn, from_Vn, to_Vn = (n * D0).function_space()
            V_cache[n] = (Vn, from_Vn, to_Vn)

            return Vn, from_Vn, to_Vn

        def mu(n, m, i, j):
            Vn, from_Vn, to_Vn = V(n)
            Vm, from_Vm, to_Vm = V(m)
            Vnm, from_Vnm, to_Vnm = V(n + m)
            return to_Vnm(from_Vn(Vn.gen(i)) * from_Vm(Vm.gen(j)))

        model = parent._model

        if model == 'large':
            div_L = 3 * D0
            L, from_L, to_L = V(3)
            from sage.rings.function_field.khuri_makdisi import KhuriMakdisi_large as KM
        elif model == 'medium':
            div_L = 2 * D0
            L, from_L, to_L = V(2)
            from sage.rings.function_field.khuri_makdisi import KhuriMakdisi_medium as KM
        elif model == 'small':
            div_L = 3 * D0
            L, from_L, to_L = V(3)
            from sage.rings.function_field.khuri_makdisi import KhuriMakdisi_small as KM

        self._div_L = div_L
        self._L = L
        self._from_L = from_L
        self._to_L = to_L

        w0 = self.point(function_field.divisor_group().zero()).defining_matrix()
        km = KM(lambda n: V(n)[0], mu, w0, D0.degree(), self._genus)

        self._km = km
        self._w0 = w0

        self._base_place = None

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: J.group()
            Group of rational points of Jacobian
             over Finite Field of size 7 (Khuri-Makdisi large model)
        """
        r = super()._repr_()
        return r + f' (Khuri-Makdisi {self._parent._model} model)'

    def _wd_from_divisor(self, x):
        """
        Return the matrix representing the divisor ``x``.

        INPUT:

        - ``x`` -- an effective divisor

        TESTS:

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_medium', base_div=h)
            sage: G = J.group()
            sage: P = C([-1,2,1]).place()
            sage: G._wd_from_divisor(2*P)
            [1 0 0 0 4 0]
            [0 1 0 0 4 6]
            [0 0 1 0 2 1]
            [0 0 0 1 1 6]
        """
        WD = (self._div_L - x).basis_function_space()
        wd = matrix([self._to_L(f) for f in WD])
        wd.echelonize()
        return wd

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        If ``x`` is an effective divisor, then it is assumed to have the same
        degree with the base divisor.

        TESTS::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: G(0)
            Point of Jacobian determined by
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 1 0 0 0]
            [0 0 0 0 0 0 0 1 0]
            sage: b = C([0,1,0]).place()
            sage: J.set_base_place(b)
            sage: p = C([-1,2,1]).place()
            sage: G(p)
            Point of Jacobian determined by
            [1 0 0 0 0 0 0 2 5]
            [0 1 0 0 0 0 0 0 3]
            [0 0 1 0 0 0 0 0 2]
            [0 0 0 1 0 0 0 0 6]
            [0 0 0 0 0 1 0 0 5]
            [0 0 0 0 0 0 1 0 1]
        """
        if x == 0:
            return self.zero()

        if isinstance(x, FunctionFieldPlace):
            if (self._base_place is not None
                and x in self._function_field.place_set()
                and x.degree() == 1):
                x = x - self._base_place
            else:
                x = x.divisor()

        if (isinstance(x, FunctionFieldDivisor)
            and x in self._function_field.divisor_group()):
            if x.degree() == 0:
                return self.point(x)
            if x.is_effective():
                if x.degree() != self._base_div_degree:
                    raise ValueError(f"effective divisor is not of degree {self._base_div_degree}")
                wd = self._wd_from_divisor(x)
                return self.element_class(self, wd)

        raise ValueError(f"cannot construct a point from {x}")

    def point(self, divisor):
        """
        Return the point represented by the divisor.

        INPUT:

        - ``divisor`` -- a divisor of degree zero

        EXAMPLES::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: b = C([0,1,0]).place()
            sage: p = C([-1,2,1]).place()
            sage: G.point(p - b)
            Point of Jacobian determined by
            [1 0 0 0 0 0 0 2 5]
            [0 1 0 0 0 0 0 0 3]
            [0 0 1 0 0 0 0 0 2]
            [0 0 0 1 0 0 0 0 6]
            [0 0 0 0 0 1 0 0 5]
            [0 0 0 0 0 0 1 0 1]
        """
        if divisor.degree() != 0:
            raise ValueError('divisor not of degree zero')

        c = divisor + self._base_div
        f = c.basis_function_space()[0]
        d = f.divisor() + c

        wd = self._wd_from_divisor(d)
        return self.element_class(self, wd)

    @cached_method
    def zero(self):
        """
        Return the zero element of this group.

        EXAMPLES::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: G.zero()
            Point of Jacobian determined by
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 1 0 0 0]
            [0 0 0 0 0 0 0 1 0]
        """
        return self.element_class(self, self._w0)


class JacobianGroup_finite_field(JacobianGroup, JacobianGroup_finite_field_base):
    """
    Jacobian groups of function fields over finite fields.

    INPUT:

    - ``parent`` -- a Jacobian

    - ``function_field`` -- a function field

    - ``base_div`` -- an effective divisor of the function field

    EXAMPLES::

        sage: # long time
        sage: k = GF(7)
        sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
        sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
        sage: h = C.function(y/x).divisor_of_poles()
        sage: J = C.jacobian(model='km_large', base_div=h)
        sage: G1 = J.group()
        sage: K = k.extension(2)
        sage: G2 = J.group(K)
        sage: G2.coerce_map_from(G1)
        Jacobian group embedding map:
          From: Group of rational points of Jacobian
                over Finite Field of size 7 (Khuri-Makdisi large model)
          To:   Group of rational points of Jacobian
                over Finite Field in z2 of size 7^2 (Khuri-Makdisi large model)
    """
    Element = JacobianPoint_finite_field

    def __init__(self, parent, function_field, base_div):
        """
        Initialize.

        TESTS::

            sage: # long time
            sage: k = GF(7)
            sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: TestSuite(G).run(skip=['_test_elements', '_test_pickling'])
        """
        super().__init__(parent, function_field, base_div)

        F = self._function_field
        K = F.constant_base_field()

        r = self._parent._function_field.constant_base_field().degree()
        frob_K = K.frobenius_endomorphism(r)

        self._frobenius_of_constant_field = frob_K

    def __iter__(self):
        """
        Return generator of points of this group.

        TESTS::

            sage: k = GF(7)
            sage: A.<x,y> = AffineSpace(k,2)
            sage: C = Curve(y^2 + x^3 + 2*x + 1).projective_closure()
            sage: b = C([0,0,1]).place()
            sage: J = C.jacobian(model='km_large', base_div=3*b)
            sage: G = J.group()
            sage: len([pt for pt in G])  # long time
            11
        """
        d0 = self._base_div.degree()
        F = self._function_field

        zero_divisor = self._km.zero_divisor()
        deg = 1
        support = []
        degrees = []
        multiples = []
        lst = []

        places_infinite = F.places_infinite()
        generators = [iter(places_infinite)]
        while True:
            while True:
                try:
                    new_pl = next(generators[-1])
                    break
                except StopIteration:
                    if deg > d0:
                        return
                    generators.append(F._places_finite(deg))
                    deg += 1
            multiples.append((d0 + 1)*[None])
            wn = self._wd_from_divisor(new_pl.divisor())
            dn = new_pl.degree()
            wr = zero_divisor
            dr = 0
            for r in range(1, d0 // dn + 1):
                wr = self._km.add_divisor(wr, wn, dr, dn)
                multiples[-1][r] = wr
                dr += dn
                for weights in WeightedIntegerVectors(d0 - dr, degrees):
                    d = dr
                    wD = wr
                    for i in range(len(support)):
                        w = weights[i]
                        if w > 0:
                            m = w * degrees[i]
                            wD = self._km.add_divisor(wD, multiples[i][w], d, m)
                            d += m
                    pt = self.element_class(self, wD)
                    if pt not in lst:
                        lst.append(pt)
                        yield pt
            support.append(new_pl)
            degrees.append(new_pl.degree())

    def _frobenius_on(self, pt):
        """
        Return the image of ``pt`` acted by the Frobenius automorphism.

        INPUT:

        - ``pt`` -- a point of ``self``

        TESTS::

            sage: # long time
            sage: k = GF(7)
            sage: A.<x,y> = AffineSpace(k,2)
            sage: C = Curve(y^2 + x^3 + 2*x + 1).projective_closure()
            sage: b = C([0,0,1]).place()
            sage: J = C.jacobian(model='km_large', base_div=3*b)
            sage: G1 = J.group()
            sage: K = k.extension(3)
            sage: G3 = J.group(K)
            sage: pt = G3.get_points(12)[-2]  # expected to be a point rational over K
            sage: pt.frobenius().frobenius().frobenius() == pt  # indirect doctest
            True
        """
        w = pt._w.apply_morphism(self._frobenius_of_constant_field)
        return self.element_class(self, w)


class Jacobian(UniqueRepresentation, Jacobian_base):
    """
    Jacobians implemented by Khuri-Makdisi's algorithms.

    EXAMPLES::

        sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
        sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
        sage: C.jacobian(model='km')
        Jacobian of Projective Plane Curve over Finite Field of size 7
         defined by x^3 - y^2*z - 2*z^3 (Khuri-Makdisi large model)
    """
    def __init__(self, function_field, base_div, model, **kwds):
        """
        Initialize.

        TESTS::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: J = C.jacobian(model='km_large')
            sage: TestSuite(J).run(skip=['_test_elements', '_test_pickling'])

        ::

            sage: # long time
            sage: J = C.jacobian(model='km_unknown')
            Traceback (most recent call last):
            ...
            ValueError: unknown model
        """
        super().__init__(function_field, base_div, **kwds)

        if model not in ['large', 'medium', 'small']:
            raise ValueError('unknown model')

        self._model = model

        if function_field.constant_base_field().is_finite():
            self._group_class = JacobianGroup_finite_field
        else:
            self._group_class = JacobianGroup

    def _repr_(self):
        """
        Return the string representation of ``self``.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: C.jacobian(model='km_large')
            Jacobian of Projective Plane Curve over Finite Field of size 7
             defined by x^3 - y^2*z - 2*z^3 (Khuri-Makdisi large model)
        """
        r = super()._repr_()
        return r + f' (Khuri-Makdisi {self._model} model)'
