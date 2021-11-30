r"""
Set of homomorphisms between two toric varieties

For schemes `X` and `Y`, this module implements the set of morphisms
`Hom(X,Y)`. This is done by
:class:`~sage.schemes.generic.homset.SchemeHomset_generic`.

As a special case, the Hom-sets can also represent the points of a
scheme. Recall that the `K`-rational points of a scheme `X` over `k`
can be identified with the set of morphisms `Spec(K) \to X`. In Sage,
the rational points are implemented by such scheme morphisms. This is
done by :class:`~sage.schemes.generic.homset.SchemeHomset_points` and
its subclasses.

.. note::

    You should not create the Hom-sets manually. Instead, use the
    :meth:`~sage.structure.parent.Hom` method that is inherited by all
    schemes.

AUTHORS:

- Volker Braun (2012-02-18): Initial version

EXAMPLES:

Here is a simple example, the projection of
`\mathbb{P}^1\times\mathbb{P}^1\to \mathbb{P}^1` ::

    sage: P1xP1 = toric_varieties.P1xP1()
    sage: P1 = toric_varieties.P1()
    sage: hom_set = P1xP1.Hom(P1);  hom_set
    Set of morphisms
      From: 2-d CPR-Fano toric variety covered by 4 affine patches
      To:   1-d CPR-Fano toric variety covered by 2 affine patches

In terms of the fan, we can define this morphism by the projection
onto the first coordinate. The Hom-set can construct the morphism from
the projection matrix alone::

    sage: hom_set(matrix([[1],[0]]))
    Scheme morphism:
      From: 2-d CPR-Fano toric variety covered by 4 affine patches
      To:   1-d CPR-Fano toric variety covered by 2 affine patches
      Defn: Defined by sending Rational polyhedral fan in 2-d lattice N
            to Rational polyhedral fan in 1-d lattice N.
    sage: _.as_polynomial_map()
    Scheme morphism:
      From: 2-d CPR-Fano toric variety covered by 4 affine patches
      To:   1-d CPR-Fano toric variety covered by 2 affine patches
      Defn: Defined on coordinates by sending [s : t : x : y] to
            [s : t]

In the case of toric algebraic schemes (defined by polynomials in
toric varieties), this module defines the underlying morphism of the
ambient toric varieties::

    sage: P1xP1.inject_variables()
    Defining s, t, x, y
    sage: S = P1xP1.subscheme([s*x-t*y])
    sage: type(S.Hom(S))
    <class 'sage.schemes.toric.homset.SchemeHomset_toric_variety_with_category'>

Finally, you can have morphisms defined through homogeneous
coordinates where the codomain is not implemented as a toric variety::

    sage: P2_toric.<x,y,z> = toric_varieties.P2()
    sage: P2_native.<u,v,w> = ProjectiveSpace(QQ, 2)
    sage: toric_to_native = P2_toric.Hom(P2_native);  toric_to_native
    Set of morphisms
      From: 2-d CPR-Fano toric variety covered by 3 affine patches
      To:   Projective Space of dimension 2 over Rational Field
    sage: type(toric_to_native)
    <class 'sage.schemes.toric.homset.SchemeHomset_toric_variety_with_category'>
    sage: toric_to_native([x^2, y^2, z^2])
    Scheme morphism:
      From: 2-d CPR-Fano toric variety covered by 3 affine patches
      To:   Projective Space of dimension 2 over Rational Field
      Defn: Defined on coordinates by sending [x : y : z] to
            (x^2 : y^2 : z^2)

    sage: native_to_toric = P2_native.Hom(P2_toric);  native_to_toric
    Set of morphisms
      From: Projective Space of dimension 2 over Rational Field
      To:   2-d CPR-Fano toric variety covered by 3 affine patches
    sage: type(native_to_toric)
    <class 'sage.schemes.generic.homset.SchemeHomset_generic_with_category'>
    sage: native_to_toric([u^2, v^2, w^2])
    Scheme morphism:
      From: Projective Space of dimension 2 over Rational Field
      To:   2-d CPR-Fano toric variety covered by 3 affine patches
      Defn: Defined on coordinates by sending (u : v : w) to
            [u^2 : v^2 : w^2]
"""



#*****************************************************************************
#       Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.categories.finite_fields import FiniteFields
from sage.rings.integer_ring import ZZ

from sage.structure.element import is_Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.geometry.fan_morphism import FanMorphism

from sage.schemes.generic.homset import (SchemeHomset_generic,
                                         SchemeHomset_points)


class SchemeHomset_toric_variety(SchemeHomset_generic):
    """
    Set of homomorphisms between two toric varieties.

    EXAMPLES::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: P1 = toric_varieties.P1()
        sage: hom_set = P1xP1.Hom(P1);  hom_set
        Set of morphisms
          From: 2-d CPR-Fano toric variety covered by 4 affine patches
          To:   1-d CPR-Fano toric variety covered by 2 affine patches
        sage: type(hom_set)
        <class 'sage.schemes.toric.homset.SchemeHomset_toric_variety_with_category'>

        sage: hom_set(matrix([[1],[0]]))
        Scheme morphism:
          From: 2-d CPR-Fano toric variety covered by 4 affine patches
          To:   1-d CPR-Fano toric variety covered by 2 affine patches
          Defn: Defined by sending Rational polyhedral fan in 2-d lattice N
                to Rational polyhedral fan in 1-d lattice N.
    """

    def __init__(self, X, Y, category=None, check=True, base=ZZ):
        """
        The Python constructor.

        INPUT:

        The same as for any homset, see
        :mod:`~sage.categories.homset`.

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: hom_set = P1xP1.Hom(P1);  hom_set
            Set of morphisms
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   1-d CPR-Fano toric variety covered by 2 affine patches

        An integral matrix defines a fan morphism, since we think of
        the matrix as a linear map on the toric lattice. This is why
        we need to ``register_conversion`` in the constructor
        below. The result is::

            sage: hom_set(matrix([[1],[0]]))
            Scheme morphism:
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   1-d CPR-Fano toric variety covered by 2 affine patches
              Defn: Defined by sending Rational polyhedral fan in 2-d lattice N
                    to Rational polyhedral fan in 1-d lattice N.
        """
        SchemeHomset_generic.__init__(self, X, Y, category=category, check=check, base=base)
        from sage.schemes.toric.variety import is_ToricVariety
        if is_ToricVariety(X) and is_ToricVariety(Y):
            self.register_conversion(MatrixSpace(ZZ, X.fan().dim(), Y.fan().dim()))

    def _element_constructor_(self, x, check=True):
        """
        Construct a scheme morphism.

        INPUT:

        - `x` -- anything that defines a morphism of toric
          varieties. A matrix, fan morphism, or a list or tuple of
          homogeneous polynomials that define a morphism.

        - ``check`` -- boolean (default: ``True``) passed onto
          functions called by this to be more careful about input
          argument type checking

        OUTPUT:

        The morphism of toric varieties determined by ``x``.

        EXAMPLES:

        First, construct from fan morphism::

            sage: dP8.<t,x0,x1,x2> = toric_varieties.dP8()
            sage: P2.<y0,y1,y2> = toric_varieties.P2()
            sage: hom_set = dP8.Hom(P2)

            sage: fm = FanMorphism(identity_matrix(2), dP8.fan(), P2.fan())
            sage: hom_set(fm)     # calls hom_set._element_constructor_()
            Scheme morphism:
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   2-d CPR-Fano toric variety covered by 3 affine patches
              Defn: Defined by sending Rational polyhedral fan in 2-d lattice N
                    to Rational polyhedral fan in 2-d lattice N.

        A matrix will automatically be converted to a fan morphism::

            sage: hom_set(identity_matrix(2))
            Scheme morphism:
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   2-d CPR-Fano toric variety covered by 3 affine patches
              Defn: Defined by sending Rational polyhedral fan in 2-d lattice N
                    to Rational polyhedral fan in 2-d lattice N.

        Alternatively, one can use homogeneous polynomials to define morphisms::

            sage: P2.inject_variables()
            Defining y0, y1, y2
            sage: dP8.inject_variables()
            Defining t, x0, x1, x2
            sage: hom_set([x0,x1,x2])
            Scheme morphism:
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   2-d CPR-Fano toric variety covered by 3 affine patches
              Defn: Defined on coordinates by sending [t : x0 : x1 : x2] to
                    [x0 : x1 : x2]

        A morphism of the coordinate ring will also work::

            sage: ring_hom = P2.coordinate_ring().hom([x0,x1,x2], dP8.coordinate_ring())
            sage: ring_hom
            Ring morphism:
              From: Multivariate Polynomial Ring in y0, y1, y2 over Rational Field
              To:   Multivariate Polynomial Ring in t, x0, x1, x2 over Rational Field
              Defn: y0 |--> x0
                    y1 |--> x1
                    y2 |--> x2
            sage: hom_set(ring_hom)
            Scheme morphism:
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   2-d CPR-Fano toric variety covered by 3 affine patches
              Defn: Defined on coordinates by sending [t : x0 : x1 : x2] to
                    [x0 : x1 : x2]
        """
        from sage.schemes.toric.morphism import SchemeMorphism_polynomial_toric_variety
        if isinstance(x, (list, tuple)):
            return SchemeMorphism_polynomial_toric_variety(self, x, check=check)

        from sage.categories.map import Map
        from sage.categories.all import Rings
        if isinstance(x, Map) and x.category_for().is_subcategory(Rings()):
            # x is a morphism of Rings
            assert x.domain() is self.codomain().coordinate_ring()
            assert x.codomain() is self.domain().coordinate_ring()
            return SchemeMorphism_polynomial_toric_variety(self, x.im_gens(), check=check)

        if is_Matrix(x):
            x = FanMorphism(x, self.domain().fan(), self.codomain().fan())
        if isinstance(x, FanMorphism):
            if x.is_dominant():
                from sage.schemes.toric.morphism import SchemeMorphism_fan_toric_variety_dominant
                return SchemeMorphism_fan_toric_variety_dominant(self, x, check=check)
            else:
                from sage.schemes.toric.morphism import SchemeMorphism_fan_toric_variety
                return SchemeMorphism_fan_toric_variety(self, x, check=check)

        raise TypeError("x must be a fan morphism or a list/tuple of polynomials")


    def _an_element_(self):
        """
        Construct a sample morphism.

        OUTPUT:

        An element of the homset.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: homset = P2.Hom(P2)
            sage: homset.an_element()   # indirect doctest
            Scheme endomorphism of 2-d CPR-Fano toric variety covered by 3 affine patches
              Defn: Defined by sending Rational polyhedral fan in 2-d lattice N to
                    Rational polyhedral fan in 2-d lattice N.
        """
        from sage.matrix.constructor import zero_matrix
        zero = zero_matrix(self.domain().dimension_relative(),
                           self.codomain().dimension_relative())
        return self(zero)

class SchemeHomset_points_toric_base(SchemeHomset_points):
    """
    Base class for homsets with toric ambient spaces

    INPUT:

    - same as for :class:`SchemeHomset_points`.

    OUTPUT:

    A scheme morphism of type
    :class:`SchemeHomset_points_toric_base`.

    EXAMPLES::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: P1xP1(QQ)
        Set of rational points of 2-d CPR-Fano toric variety
        covered by 4 affine patches

    TESTS::

        sage: import sage.schemes.toric.homset as HOM
        sage: HOM.SchemeHomset_points_toric_base(Spec(QQ), P1xP1)
        Set of rational points of 2-d CPR-Fano toric variety covered by 4 affine patches
    """

    def is_finite(self):
        """
        Return whether there are finitely many points.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: P2.point_set().is_finite()
            False
            sage: P2.change_ring(GF(7)).point_set().is_finite()
            True
        """
        variety = self.codomain()
        return variety.dimension() == 0 or variety.base_ring().is_finite()

    def _naive_enumerator(self, ring=None):
        """
        The naive enumerator over points of the toric variety.

        INPUT:

        - ``ring`` -- a ring (optional; defaults to the base ring of
          the toric variety). The ring over which the points are
          considered.

        OUTPUT:

        A :class:`sage.schemes.toric.points.NaiveFinitePointEnumerator`
        instance that can be used to iterate over the points.

        EXAMPLES::

            sage: P123 = toric_varieties.P2_123(base_ring=GF(3))
            sage: point_set = P123.point_set()
            sage: next(iter(point_set._naive_enumerator()))
            (0, 0, 1)
            sage: next(iter(point_set))
            [0 : 0 : 1]
        """
        from sage.schemes.toric.points import NaiveFinitePointEnumerator
        variety = self.codomain()
        if ring is None:
            ring = variety.base_ring()
        return NaiveFinitePointEnumerator(variety.fan(), ring)

    def _finite_field_enumerator(self, finite_field=None):
        """
        Fast enumerator for points of the toric variety.

        INPUT:

        - ``finite_field`` -- a finite field (optional; defaults to
          the base ring of the toric variety). The ring over which the
          points are considered.

        OUTPUT:

        A
        :class:`sage.schemes.toric.points.FiniteFieldPointEnumerator`
        instance that can be used to iterate over the points.

        EXAMPLES::

            sage: P123 = toric_varieties.P2_123(base_ring=GF(3))
            sage: point_set = P123.point_set()
            sage: next(iter(point_set._finite_field_enumerator()))
            (0, 0, 1)
            sage: next(iter(point_set))
            [0 : 0 : 1]
        """
        from sage.schemes.toric.points import FiniteFieldPointEnumerator
        variety = self.codomain()
        if finite_field is None:
            finite_field = variety.base_ring()
        if finite_field not in FiniteFields():
            raise ValueError('not a finite field')
        return FiniteFieldPointEnumerator(variety.fan(), finite_field)

    def _enumerator(self):
        """
        Return the most suitable enumerator for points.

        OUTPUT:

        An iterable that yields the coordinates of all points as
        tuples.

        EXAMPLES::

            sage: P123 = toric_varieties.P2_123(base_ring=GF(3))
            sage: point_set = P123.point_set()
            sage: point_set._enumerator()
            <sage.schemes.toric.points.FiniteFieldPointEnumerator object at 0x...>
        """
        ring = self.domain().base_ring()
        if ring in FiniteFields():
            return self._finite_field_enumerator()
        elif ring.is_finite():
            return self._naive_enumerator()
        else:
            from sage.schemes.toric.points import InfinitePointEnumerator
            return InfinitePointEnumerator(self.codomain().fan(), ring)


class SchemeHomset_points_toric_field(SchemeHomset_points_toric_base):
    r"""
    Set of rational points of a toric variety.

    You should not use this class directly. Instead, use the
    :meth:`~sage.schemes.generic.scheme.Scheme.point_set` method to
    construct the point set of a toric variety.

    INPUT:

    - same as for :class:`~sage.schemes.generic.homset.SchemeHomset_points`.

    OUTPUT:

    A scheme morphism of type
    :class:`SchemeHomset_points_toric_field`.

    EXAMPLES::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: P1xP1.point_set()
        Set of rational points of 2-d CPR-Fano toric variety
        covered by 4 affine patches
        sage: P1xP1(QQ)
        Set of rational points of 2-d CPR-Fano toric variety
        covered by 4 affine patches

    The quotient `\mathbb{P}^2 / \ZZ_3` over `GF(7)` by the diagonal
    action. This is tricky because the base field has a 3-rd root of
    unity::

        sage: fan = NormalFan(ReflexivePolytope(2, 0))
        sage: X = ToricVariety(fan, base_field=GF(7))
        sage: point_set = X.point_set()
        sage: point_set.cardinality()
        21
        sage: sorted(X.point_set().list())
        [[0 : 0 : 1], [0 : 1 : 0], [0 : 1 : 1], [0 : 1 : 3],
         [1 : 0 : 0], [1 : 0 : 1], [1 : 0 : 3], [1 : 1 : 0],
         [1 : 1 : 1], [1 : 1 : 2], [1 : 1 : 3], [1 : 1 : 4],
         [1 : 1 : 5], [1 : 1 : 6], [1 : 3 : 0], [1 : 3 : 1],
         [1 : 3 : 2], [1 : 3 : 3], [1 : 3 : 4], [1 : 3 : 5],
         [1 : 3 : 6]]

    As for a non-compact example, the blow-up of the plane is the line
    bundle $O_{\mathbf{P}^1}(-1)$. Its point set is the Cartesian
    product of the points on the base $\mathbf{P}^1$ with the points
    on the fiber::

        sage: fan = Fan([Cone([(1,0), (1,1)]), Cone([(1,1), (0,1)])])
        sage: blowup_plane = ToricVariety(fan, base_ring=GF(3))
        sage: point_set = blowup_plane.point_set()
        sage: sorted(point_set.list())
        [[0 : 1 : 0], [0 : 1 : 1], [0 : 1 : 2],
         [1 : 0 : 0], [1 : 0 : 1], [1 : 0 : 2],
         [1 : 1 : 0], [1 : 1 : 1], [1 : 1 : 2],
         [1 : 2 : 0], [1 : 2 : 1], [1 : 2 : 2]]

    Toric varieties with torus factors (that is, where the fan is not
    full-dimensional) also work::

        sage: F_times_Fstar = ToricVariety(Fan([Cone([(1,0)])]), base_field=GF(3))
        sage: sorted(F_times_Fstar.point_set().list())
        [[0 : 1], [0 : 2], [1 : 1], [1 : 2], [2 : 1], [2 : 2]]

    TESTS::

        sage: import sage.schemes.toric.homset as HOM
        sage: HOM.SchemeHomset_points_toric_field(Spec(QQ), P1xP1)
        Set of rational points of 2-d CPR-Fano toric variety covered by 4 affine patches
    """

    def cardinality(self):
        r"""
        Return the number of points of the toric variety.

        OUTPUT:

        An integer or infinity. The cardinality of the set of points.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: V = ToricVariety(FaceFan(o))
            sage: V.change_ring(GF(2)).point_set().cardinality()
            27
            sage: V.change_ring(GF(8, "a")).point_set().cardinality()
            729
            sage: V.change_ring(GF(101)).point_set().cardinality()
            1061208

        For non-smooth varieties over finite fields, the homogeneous
        rescalings are solved. This is somewhat slower::

            sage: fan = NormalFan(ReflexivePolytope(2, 0))
            sage: X = ToricVariety(fan, base_field=GF(7))
            sage: X.point_set().cardinality()
            21

        Fulton's formula does not apply since the variety is not
        smooth. And, indeed, naive application gives a different
        result::

            sage: q = X.base_ring().order()
            sage: n = X.dimension()
            sage: d = map(len, fan().cones())
            sage: sum(dk * (q-1)**(n-k) for k, dk in enumerate(d))
            57

        Over infinite fields the number of points is not very tricky::

            sage: V.count_points()
            +Infinity

        ALGORITHM:

        Uses the formula in Fulton [Ful1993]_, section 4.5.

        AUTHORS:

        - Beth Malmskog (2013-07-14)

        - Adriana Salerno (2013-07-14)

        - Yiwei She (2013-07-14)

        - Christelle Vincent (2013-07-14)

        - Ursula Whitcher (2013-07-14)
        """
        variety = self.codomain()
        if not variety.base_ring().is_finite():
            if variety.dimension_relative() == 0:
                return ZZ.one()
            else:
                from sage.rings.infinity import Infinity
                return Infinity
        if not variety.is_smooth():
            try:
                return self._enumerator().cardinality()
            except AttributeError:
                return super(SchemeHomset_points_toric_field, self).cardinality()
        q = variety.base_ring().order()
        n = variety.dimension()
        d = map(len, variety.fan().cones())
        return sum(dk * (q-1)**(n-k) for k, dk in enumerate(d))

    def __iter__(self):
        """
        Iterate over the points of the variety.

        OUTPUT:

        Iterator over points.

        EXAMPLES::

            sage: P123 = toric_varieties.P2_123(base_ring=GF(3))
            sage: point_set = P123.point_set()
            sage: next(iter(point_set.__iter__()))
            [0 : 0 : 1]
            sage: next(iter(point_set))  # syntactic sugar
            [0 : 0 : 1]
        """
        for pt in self._enumerator():
            yield self(pt)


class SchemeHomset_points_subscheme_toric_field(SchemeHomset_points_toric_base):

    def _enumerator(self):
        """
        Return the most suitable enumerator for points.

        OUTPUT:

        An iterable that yields the coordinates of all points as
        tuples.

        EXAMPLES::

            sage: P123 = toric_varieties.P2_123(base_ring=GF(3))
            sage: point_set = P123.point_set()
            sage: point_set._enumerator()
            <sage.schemes.toric.points.FiniteFieldPointEnumerator object at 0x...>
        """
        ambient = super(
            SchemeHomset_points_subscheme_toric_field, self
        )._enumerator()
        ring = self.domain().base_ring()
        if ring in FiniteFields():
            from sage.schemes.toric.points import FiniteFieldSubschemePointEnumerator
            Enumerator = FiniteFieldSubschemePointEnumerator
        else:
            from sage.schemes.toric.points import NaiveSubschemePointEnumerator
            Enumerator = NaiveSubschemePointEnumerator
        return Enumerator(self.codomain().defining_polynomials(), ambient)

    def __iter__(self):
        """
        Iterate over the points of the variety.

        OUTPUT:

        Iterator over points.

        EXAMPLES::

            sage: P2.<x,y,z> = toric_varieties.P2(base_ring=GF(5))
            sage: cubic = P2.subscheme([x^3 + y^3 + z^3])
            sage: list(cubic.point_set())
            [[0 : 1 : 4], [1 : 0 : 4], [1 : 4 : 0], [1 : 1 : 2], [1 : 2 : 1], [1 : 3 : 3]]
            sage: cubic.point_set().cardinality()
            6
        """
        for p in self._enumerator():
            yield self(p)

    def cardinality(self):
        """
        Return the number of points of the toric variety.

        OUTPUT:

        An integer or infinity. The cardinality of the set of points.

        EXAMPLES::

            sage: P2.<x,y,z> = toric_varieties.P2(base_ring=GF(5))
            sage: cubic = P2.subscheme([x^3 + y^3 + z^3])
            sage: list(cubic.point_set())
            [[0 : 1 : 4], [1 : 0 : 4], [1 : 4 : 0], [1 : 1 : 2], [1 : 2 : 1], [1 : 3 : 3]]
            sage: cubic.point_set().cardinality()
            6
        """
        try:
            return self._enumerator().cardinality()
        except AttributeError:
            return super(SchemeHomset_points_subscheme_toric_field, self).cardinality()


