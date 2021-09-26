# -*- coding: utf-8 -*-
r"""
Methods of constructing simplicial sets

This implements various constructions on simplicial sets:
subsimplicial sets, pullbacks, products, pushouts, quotients, wedges,
disjoint unions, smash products, cones, and suspensions. The best way
to access these is with methods attached to simplicial sets
themselves, as in the following.

EXAMPLES::

    sage: K = simplicial_sets.Simplex(1)
    sage: square = K.product(K)

    sage: K = simplicial_sets.Simplex(1)
    sage: endpoints = K.n_skeleton(0)
    sage: circle = K.quotient(endpoints)

The mapping cone of a morphism of simplicial sets is constructed as a
pushout::

    sage: eta = simplicial_sets.HopfMap()
    sage: CP2 = eta.mapping_cone()
    sage: type(CP2)
    <class 'sage.topology.simplicial_set_constructions.PushoutOfSimplicialSets_finite_with_category'>

See the main documentation for simplicial sets, as well as for the
classes for pushouts, pullbacks, etc., for more details.

Many of the classes defined here inherit from
:class:`sage.structure.unique_representation.UniqueRepresentation`. This
means that they produce identical output if given the same input, so
for example, if ``K`` is a simplicial set, calling ``K.suspension()``
twice returns the same result both times::

    sage: CP2.suspension() is CP2.suspension()
    True

So on one hand, a command like ``simplicial_sets.Sphere(2)``
constructs a distinct copy of a 2-sphere each time it is called; on
the other, once you have constructed a 2-sphere, then constructing its
cone, its suspension, its product with another simplicial set, etc.,
will give you the same result each time::

    sage: simplicial_sets.Sphere(2) == simplicial_sets.Sphere(2)
    False
    sage: S2 = simplicial_sets.Sphere(2)
    sage: S2.product(S2) == S2.product(S2)
    True
    sage: S2.disjoint_union(CP2, S2) == S2.disjoint_union(CP2, S2)
    True

AUTHORS:

- John H. Palmieri (2016-07)
"""
#*****************************************************************************
#  Copyright (C) 2016 John H. Palmieri <palmieri at math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#
#*****************************************************************************

import itertools

from sage.graphs.graph import Graph
from sage.misc.latex import latex
from sage.sets.set import Set
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from .simplicial_set import AbstractSimplex, \
    SimplicialSet_arbitrary, SimplicialSet_finite, \
    standardize_degeneracies, face_degeneracies
from .simplicial_set_examples import Empty, Point

from sage.misc.lazy_import import lazy_import
lazy_import('sage.categories.simplicial_sets', 'SimplicialSets')

########################################################################
# classes which inherit from SimplicialSet_arbitrary

# Note: many of the classes below for infinite simplicial sets have an
# attribute '_n_skeleton'. This is used to cache the highest
# dimensional skeleton calculated so far for this simplicial set,
# along with its dimension, so for example, the starting value is
# often (-1, Empty()): the (-1)-skeleton is the empty simplicial
# set. It gets used and updated in the n_skeleton method.

class SubSimplicialSet(SimplicialSet_finite, UniqueRepresentation):
    @staticmethod
    def __classcall__(self, data, ambient=None):
        """
        Convert ``data`` from a dict to a tuple.

        TESTS::

            sage: from sage.topology.simplicial_set_constructions import SubSimplicialSet
            sage: K = simplicial_sets.Simplex(2)
            sage: e = K.n_cells(1)[0]
            sage: A = SubSimplicialSet({e: K.faces(e)}, ambient=K)
            sage: B = SubSimplicialSet({e: list(K.faces(e))}, ambient=K)
            sage: A == B
            True
        """
        L = []
        for x in data:
            if data[x] is None:
                L.append((x, None))
            else:
                L.append((x, tuple(data[x])))
        return super(SubSimplicialSet, self).__classcall__(self, tuple(L), ambient)

    def __init__(self, data, ambient=None):
        r"""
        Return a finite simplicial set as a subsimplicial set of another
        simplicial set.

        This keeps track of the ambient simplicial set and the
        inclusion map from the subcomplex into it.

        INPUT:

        - ``data`` -- the data defining the subset: a dictionary where
          the keys are simplices from the ambient simplicial set and
          the values are their faces.

        - ``ambient`` -- the ambient simplicial set. If omitted, use
          the same simplicial set as the subset and the ambient
          complex.

        EXAMPLES::

            sage: S3 = simplicial_sets.Sphere(3)
            sage: K = simplicial_sets.KleinBottle()
            sage: X = S3.disjoint_union(K)
            sage: Y = X.structure_map(0).image() # the S3 summand
            sage: Y.inclusion_map()
            Simplicial set morphism:
              From: Simplicial set with 2 non-degenerate simplices
              To:   Disjoint union: (S^3 u Klein bottle)
              Defn: [v_0, sigma_3] --> [v_0, sigma_3]
            sage: Y.ambient_space()
            Disjoint union: (S^3 u Klein bottle)

        TESTS::

            sage: T = simplicial_sets.Torus()
            sage: latex(T.n_skeleton(2))
            S^{1} \times S^{1}

            sage: T.n_skeleton(1).n_skeleton(1) == T.n_skeleton(1)
            True

            sage: T.n_skeleton(1) is T.n_skeleton(1)
            True
        """
        data = dict(data)
        if ambient is None:
            ambient = self
        if (ambient.is_pointed()
            and hasattr(ambient, '_basepoint')
            and ambient.base_point() in data):
            SimplicialSet_finite.__init__(self, data, base_point=ambient.base_point())
        else:
            SimplicialSet_finite.__init__(self, data)
        if self == ambient:
            if hasattr(ambient, '__custom_name'):
                self.rename(str(ambient))
            self._latex_name = latex(ambient)
        # When constructing the inclusion map, we do not need to check
        # the validity of the morphism, and more importantly, we
        # cannot check it in the infinite case: the appropriate data
        # may not have yet been constructed. So use "check=False".
        self._inclusion = self.Hom(ambient)({x:x for x in data}, check=False)

    def inclusion_map(self):
        r"""
        Return the inclusion map from this subsimplicial set into its
        ambient space.

        EXAMPLES::

            sage: RP6 = simplicial_sets.RealProjectiveSpace(6)
            sage: K = RP6.n_skeleton(2)
            sage: K.inclusion_map()
            Simplicial set morphism:
              From: Simplicial set with 3 non-degenerate simplices
              To:   RP^6
              Defn: [1, f, f * f] --> [1, f, f * f]

        `RP^6` itself is constructed as a subsimplicial set of
        `RP^\infty`::

            sage: latex(RP6.inclusion_map())
            RP^{6} \to RP^{\infty}
        """
        return self._inclusion

    def ambient_space(self):
        """
        Return the simplicial set of which this is a subsimplicial set.

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: eight = T.wedge_as_subset()
            sage: eight
            Simplicial set with 3 non-degenerate simplices
            sage: eight.fundamental_group()
            Finitely presented group < e0, e1 |  >
            sage: eight.ambient_space()
            Torus
        """
        return self._inclusion.codomain()


class PullbackOfSimplicialSets(SimplicialSet_arbitrary, UniqueRepresentation):
    @staticmethod
    def __classcall_private__(self, maps=None):
        """
        TESTS::

            sage: from sage.topology.simplicial_set_constructions import PullbackOfSimplicialSets
            sage: S2 = simplicial_sets.Sphere(2)
            sage: one = S2.Hom(S2).identity()
            sage: PullbackOfSimplicialSets([one, one]) == PullbackOfSimplicialSets((one, one))
            True
        """
        if maps:
            return super(PullbackOfSimplicialSets, self).__classcall__(self, tuple(maps))
        return super(PullbackOfSimplicialSets, self).__classcall__(self)

    def __init__(self, maps=None):
        r"""
        Return the pullback obtained from the morphisms ``maps``.

        INPUT:

        - ``maps`` -- a list or tuple of morphisms of simplicial sets

        If only a single map `f: X \to Y` is given, then return
        `X`. If no maps are given, return the one-point simplicial
        set. Otherwise, given a simplicial set `Y` and maps `f_i: X_i
        \to Y` for `0 \leq i \leq m`, construct the pullback `P`: see
        :wikipedia:`Pullback_(category_theory)`. This is constructed
        as pullbacks of sets for each set of `n`-simplices, so `P_n`
        is the subset of the product `\prod (X_i)_n` consisting of
        those elements `(x_i)` for which `f_i(x_i) = f_j(x_j)` for all
        `i`, `j`.

        This is pointed if the maps `f_i` are.

        EXAMPLES:

        The pullback of a quotient map by a subsimplicial set and the
        base point map gives a simplicial set isomorphic to the
        original subcomplex::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: K = RP5.quotient(RP5.n_skeleton(2))
            sage: X = K.pullback(K.quotient_map(), K.base_point_map())
            sage: X.homology() == RP5.n_skeleton(2).homology()
            True

        Pullbacks of identity maps::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: one = S2.Hom(S2).identity()
            sage: P = S2.pullback(one, one)
            sage: P.homology()
            {0: 0, 1: 0, 2: Z}

        The pullback is constructed in terms of the product -- of
        course, the product is a special case of the pullback -- and
        the simplices are named appropriately::

            sage: P.nondegenerate_simplices()
            [(v_0, v_0), (sigma_2, sigma_2)]
        """
        # Import this here to prevent circular imports.
        from sage.topology.simplicial_set_morphism import SimplicialSetMorphism
        if maps and any(not isinstance(f, SimplicialSetMorphism) for f in maps):
            raise ValueError('the maps must be morphisms of simplicial sets')

        Cat = SimplicialSets()
        if maps:
            if all(f.domain().is_finite() for f in maps):
                Cat = Cat.Finite()
        if all(f.is_pointed() for f in maps):
            Cat = Cat.Pointed()
        Parent.__init__(self, category=Cat)
        self._maps = maps
        self._n_skeleton = (-1, Empty())

    def n_skeleton(self, n):
        """
        Return the `n`-skeleton of this simplicial set.

        That is, the simplicial set generated by all nondegenerate
        simplices of dimension at most `n`.

        INPUT:

        - ``n`` -- the dimension

        The `n`-skeleton of the pullback is computed as the pullback
        of the `n`-skeleta of the component simplicial sets.

        EXAMPLES::

            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: one = Hom(B,B).identity()
            sage: c = Hom(B,B).constant_map()
            sage: P = B.pullback(one, c)
            sage: P.n_skeleton(2)
            Pullback of maps:
              Simplicial set endomorphism of Simplicial set with 3 non-degenerate simplices
                Defn: Identity map
              Simplicial set endomorphism of Simplicial set with 3 non-degenerate simplices
                Defn: Constant map at 1
            sage: P.n_skeleton(3).homology()
            {0: 0, 1: C2, 2: 0, 3: Z}
        """
        if self.is_finite():
            maps = self._maps
            if maps:
                codomain = SimplicialSet_finite.n_skeleton(maps[0].codomain(), n)
                domains = [SimplicialSet_finite.n_skeleton(f.domain(), n) for f in maps]
                new_maps = [f.n_skeleton(n, d, codomain) for (f, d) in zip(maps, domains)]
                return PullbackOfSimplicialSets_finite(new_maps)
            return PullbackOfSimplicialSets_finite(maps)
        start, skel = self._n_skeleton
        if start == n:
            return skel
        elif start > n:
            return skel.n_skeleton(n)
        ans = PullbackOfSimplicialSets_finite([f.n_skeleton(n) for f in self._maps])
        self._n_skeleton = (n, ans)
        return ans

    def defining_map(self, i):
        r"""
        Return the `i`-th map defining the pullback.

        INPUT:

        - ``i`` -- integer

        If this pullback was constructed as ``Y.pullback(f_0, f_1, ...)``,
        this returns `f_i`.

        EXAMPLES::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: K = RP5.quotient(RP5.n_skeleton(2))
            sage: Y = K.pullback(K.quotient_map(), K.base_point_map())
            sage: Y.defining_map(1)
            Simplicial set morphism:
              From: Point
              To:   Quotient: (RP^5/Simplicial set with 3 non-degenerate simplices)
              Defn: Constant map at *
            sage: Y.defining_map(0).domain()
            RP^5
        """
        return self._maps[i]

    def _repr_(self):
        """
        Print representation

        EXAMPLES::

            sage: S3 = simplicial_sets.Sphere(3)
            sage: c = Hom(S3,S3).constant_map()
            sage: one = Hom(S3, S3).identity()
            sage: S3.pullback(c, one)
            Pullback of maps:
              Simplicial set endomorphism of S^3
                Defn: Constant map at v_0
              Simplicial set endomorphism of S^3
                Defn: Identity map
        """
        if not self._maps:
            return 'Point'
        s = 'Pullback of maps:'
        for f in self._maps:
            t = '\n' + str(f)
            s += t.replace('\n', '\n  ')
        return s


class PullbackOfSimplicialSets_finite(PullbackOfSimplicialSets, SimplicialSet_finite):
    """
    The pullback of finite simplicial sets obtained from ``maps``.

    When the simplicial sets involved are all finite, there are more
    methods available to the resulting pullback, as compared to case
    when some of the components are infinite: the structure maps from
    the pullback and the pullback's universal property: see
    :meth:`structure_map` and :meth:`universal_property`.
    """
    @staticmethod
    def __classcall_private__(self, maps=None):
        """
        TESTS::

            sage: from sage.topology.simplicial_set_constructions import PullbackOfSimplicialSets_finite
            sage: S2 = simplicial_sets.Sphere(2)
            sage: one = S2.Hom(S2).identity()
            sage: PullbackOfSimplicialSets_finite([one, one]) == PullbackOfSimplicialSets_finite((one, one))
            True
        """
        if maps:
            return super(PullbackOfSimplicialSets_finite, self).__classcall__(self, tuple(maps))
        return super(PullbackOfSimplicialSets_finite, self).__classcall__(self)

    def __init__(self, maps=None):
        r"""
        Return the pullback obtained from the morphisms ``maps``.

        See :class:`PullbackOfSimplicialSets` for more information.

        INPUT:

        - ``maps`` -- a list or tuple of morphisms of simplicial sets

        EXAMPLES::

            sage: eta = simplicial_sets.HopfMap()
            sage: S3 = eta.domain()
            sage: S2 = eta.codomain()
            sage: c = Hom(S2,S2).constant_map()
            sage: S2.pullback(eta, c).is_finite()
            True

            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: one = Hom(B,B).identity()
            sage: c = Hom(B,B).constant_map()
            sage: B.pullback(one, c).is_finite()
            False

        TESTS::

            sage: P = simplicial_sets.Point()
            sage: P.pullback(P.constant_map(), P.constant_map())
            Pullback of maps:
              Simplicial set endomorphism of Point
                Defn: Identity map
              Simplicial set endomorphism of Point
                Defn: Identity map
        """
        # Import this here to prevent circular imports.
        from sage.topology.simplicial_set_morphism import SimplicialSetMorphism
        if maps and any(not isinstance(f, SimplicialSetMorphism) for f in maps):
            raise ValueError('the maps must be morphisms of simplicial sets')
        if not maps:
            star = AbstractSimplex(0, name='*')
            SimplicialSet_finite.__init__(self, {star: None}, base_point=star, name='Point')
            self._maps = ()
            self._translation = {}
            return
        if len(maps) == 1:
            f = maps[0]
            if f.is_pointed():
                SimplicialSet_finite.__init__(self, f.domain().face_data(),
                                       base_point=f.domain().base_point())
            else:
                SimplicialSet_finite.__init__(self, f.domain().face_data())
            self._maps = (f,)
            return
        codomain = maps[0].codomain()
        if any(codomain != f.codomain() for f in maps[1:]):
            raise ValueError('the codomains of the maps must be equal')
        # Now construct the pullback by constructing the product and only
        # keeping the appropriate simplices.
        domains = [f.domain() for f in maps]
        nondegen = [X.nondegenerate_simplices() for X in domains]
        data_factors = [X.face_data() for X in domains]
        # data: dictionary to construct the new simplicial set.
        data = {}
        # translate: keep track of the nondegenerate simplices in the
        # new simplicial set for computing faces: keys are tuples of
        # pairs (sigma, degens), with sigma a nondegenerate simplex in
        # one of the factors, degens the tuple of applied
        # degeneracies. The associated value is the actual simplex in
        # the product.
        translate = {}
        for simplices in itertools.product(*nondegen):
            dims =  [_.dimension() for _ in simplices]
            dim_max = max(dims)
            sum_dims = sum(dims)
            for d in range(dim_max, sum_dims + 1):
                S = Set(range(d))
                # Is there a way to speed up the following? Given the
                # tuple dims=(n_1, n_2, ..., n_k) and given d between
                # max(dims) and sum(dims), we are trying to construct
                # k-tuples of subsets (D_1, D_2, ..., D_k) of range(d)
                # such that the intersection of all of the D_i's is
                # empty.
                for I in itertools.product(*[S.subsets(d - _) for _ in dims]):
                    if set.intersection(*[set(_) for _ in I]):
                        # To get a nondegenerate face, can't have a
                        # degeneracy in common for all the factors.
                        continue
                    degens = [tuple(sorted(_, reverse=True)) for _ in I]

                    sigma = simplices[0].apply_degeneracies(*degens[0])
                    target = maps[0](sigma)
                    if any(target != f(tau.apply_degeneracies(*degen))
                               for (f, tau, degen) in zip(maps[1:], simplices[1:], degens[1:])):
                        continue

                    simplex_factors = tuple(zip(simplices, tuple(degens)))
                    s = '(' + ', '.join(['{}'.format(_[0].apply_degeneracies(*_[1]))
                                         for _ in simplex_factors]) + ')'
                    ls = '(' + ', '.join(['{}'.format(latex(_[0].apply_degeneracies(*_[1])))
                                          for _ in simplex_factors]) + ')'
                    simplex = AbstractSimplex(d, name=s, latex_name=ls)
                    translate[simplex_factors] = simplex
                    # Now compute the faces of simplex.
                    if d == 0:
                        # It's a vertex, so it has no faces.
                        faces = None
                    else:
                        faces = []
                        for i in range(d+1):
                            # Compute d_i on simplex.
                            #
                            # face_degens: tuple of pairs (J, t): J is the
                            # list of degeneracies to apply to the
                            # corresponding entry in simplex_factors, t is
                            # the face map to apply.
                            face_degens = [face_degeneracies(i, _) for _ in degens]
                            face_factors = []
                            new_degens = []
                            for x, Face, face_dict in zip(simplices, face_degens, data_factors):
                                J = Face[0]
                                t = Face[1]
                                if t is None:
                                    face_factors.append(x.nondegenerate())
                                else:
                                    underlying = face_dict[x][t]
                                    temp_degens = underlying.degeneracies()
                                    underlying = underlying.nondegenerate()
                                    J = standardize_degeneracies(*(J + list(temp_degens)))
                                    face_factors.append(underlying)
                                new_degens.append(J)

                            # By the simplicial identities, s_{i_1}
                            # s_{i_2} ... s_{i_n} z (if decreasing) is in
                            # the image of s_{i_k} for each k.
                            #
                            # So find the intersection K of each J, the
                            # degeneracies applied to left_face and
                            # right_face. Then the face will be s_{K}
                            # (s_{J'_L} left_face, s_{J'_R} right_face),
                            # where you get J'_L from J_L by pulling out K
                            # from J_L.
                            #
                            # J'_L is obtained as follows: for each j in
                            # J_L, decrease j by q if q = #{x in K: x < j}
                            K = set.intersection(*[set(J) for J in new_degens])

                            face_degens = []
                            for J in new_degens:
                                new_J = []
                                for j in J:
                                    if j not in K:
                                        q = len([x for x in K if x < j])
                                        new_J.append(j - q)
                                face_degens.append(tuple(new_J))
                            K = sorted(K, reverse=True)
                            underlying_face = translate[tuple(zip(tuple(face_factors), tuple(face_degens)))]
                            faces.append(underlying_face.apply_degeneracies(*K))
                        data[simplex] = faces

        if all(f.is_pointed() for f in maps):
            basept = translate[tuple([(sset.base_point(), ()) for sset in domains])]
            if not data:
                data = {basept: None}
            SimplicialSet_finite.__init__(self, data, base_point=basept)
        else:
            SimplicialSet_finite.__init__(self, data)
        self._maps = maps
        # self._translation: tuple converted from dict. keys: tuples
        # of pairs (sigma, degens), with sigma a nondegenerate simplex
        # in one of the factors, degens the tuple of applied
        # degeneracies. The associated value is the actual simplex in
        # the product.
        self._translation = tuple(translate.items())

    def structure_map(self, i):
        r"""
        Return the `i`-th projection map of the pullback.

        INPUT:

        - ``i`` -- integer

        If this pullback `P` was constructed as ``Y.pullback(f_0, f_1,
        ...)``, where `f_i: X_i \to Y`, then there are structure maps
        `\bar{f}_i: P \to X_i`. This method constructs `\bar{f}_i`.

        EXAMPLES::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: K = RP5.quotient(RP5.n_skeleton(2))
            sage: Y = K.pullback(K.quotient_map(), K.base_point_map())
            sage: Y.structure_map(0)
            Simplicial set morphism:
              From: Pullback of maps:
              Simplicial set morphism:
                From: RP^5
                To:   Quotient: (RP^5/Simplicial set with 3 non-degenerate simplices)
                Defn: [1, f, f * f, f * f * f, f * f * f * f, f * f * f * f * f] --> [*, s_0 *, s_1 s_0 *, f * f * f, f * f * f * f, f * f * f * f * f]
              Simplicial set morphism:
                From: Point
                To:   Quotient: (RP^5/Simplicial set with 3 non-degenerate simplices)
                Defn: Constant map at *
              To:   RP^5
              Defn: [(1, *), (f, s_0 *), (f * f, s_1 s_0 *)] --> [1, f, f * f]
            sage: Y.structure_map(1).codomain()
            Point

        These maps are also accessible via ``projection_map``::

            sage: Y.projection_map(1).codomain()
            Point
        """
        if len(self._maps) == 1:
            return self.Hom(self).identity()
        f = {}
        for x in self._translation:
            f[x[1]] = x[0][i][0].apply_degeneracies(*x[0][i][1])
        codomain = self.defining_map(i).domain()
        return self.Hom(codomain)(f)

    projection_map = structure_map

    def universal_property(self, *maps):
        r"""
        Return the map induced by ``maps``.

        INPUT:

        - ``maps`` -- maps from a simplicial set `Z` to the "factors"
          `X_i` forming the pullback.

        If the pullback `P` is formed by maps `f_i: X_i \to Y`, then
        given maps `g_i: Z \to X_i` such that `f_i g_i = f_j g_j` for
        all `i`, `j`, then there is a unique map `g: Z \to P` making
        the appropriate diagram commute. This constructs that map.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = S1.product(S1)
            sage: K = T.factor(0, as_subset=True)
            sage: f = S1.Hom(T)({S1.n_cells(0)[0]:K.n_cells(0)[0], S1.n_cells(1)[0]:K.n_cells(1)[0]})
            sage: P = S1.product(T)
            sage: P.universal_property(S1.Hom(S1).identity(), f)
            Simplicial set morphism:
              From: S^1
              To:   S^1 x S^1 x S^1
              Defn: [v_0, sigma_1] --> [(v_0, (v_0, v_0)), (sigma_1, (sigma_1, s_0 v_0))]
        """
        if len(self._maps) != len(maps):
            raise ValueError('wrong number of maps specified')
        if len(self._maps) == 1:
            return maps[0]
        domain = maps[0].domain()
        if any(g.domain() != domain for g in maps[1:]):
            raise ValueError('the maps do not all have the same codomain')
        composite = self._maps[0] * maps[0]
        if any(f*g != composite for f,g in zip(self._maps[1:], maps[1:])):
            raise ValueError('the maps are not compatible')
        data = {}
        translate = dict(self._translation)
        for sigma in domain.nondegenerate_simplices():
            target = tuple([(f(sigma).nondegenerate(), tuple(f(sigma).degeneracies()))
                               for f in maps])
            # If there any degeneracies in common, remove them: the
            # dictionary "translate" has nondegenerate simplices as
            # its keys.
            in_common = set.intersection(*[set(_[1]) for _ in target])
            if in_common:
                target = tuple((tau, tuple(sorted(set(degens).difference(in_common),
                                                  reverse=True)))
                               for tau, degens in target)
            in_common = sorted(in_common, reverse=True)
            data[sigma] = translate[target].apply_degeneracies(*in_common)
        return domain.Hom(self)(data)

class Factors(object):
    """
    Classes which inherit from this should define a ``_factors``
    attribute for their instances, and this class accesses that
    attribute. This is used by :class:`ProductOfSimplicialSets`,
    :class:`WedgeOfSimplicialSets`, and
    :class:`DisjointUnionOfSimplicialSets`.
    """
    def factors(self):
        """
        Return the factors involved in this construction of simplicial sets.

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: S3 = simplicial_sets.Sphere(3)
            sage: S2.wedge(S3).factors() == (S2, S3)
            True
            sage: S2.product(S3).factors()[0]
            S^2
        """
        return self._factors

    def factor(self, i):
        r"""
        Return the $i$-th factor of this construction of simplicial sets.

        INPUT:

        - ``i`` -- integer, the index of the factor

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: S3 = simplicial_sets.Sphere(3)
            sage: K = S2.disjoint_union(S3)
            sage: K.factor(0)
            S^2
            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: X = B.wedge(S3, B)
            sage: X.factor(1)
            S^3
            sage: X.factor(2)
            Classifying space of Multiplicative Abelian group isomorphic to C2
        """
        return self.factors()[i]


class ProductOfSimplicialSets(PullbackOfSimplicialSets, Factors):
    @staticmethod
    def __classcall__(cls, factors=None):
        """
        TESTS::

            sage: from sage.topology.simplicial_set_constructions import ProductOfSimplicialSets
            sage: S2 = simplicial_sets.Sphere(2)
            sage: ProductOfSimplicialSets([S2, S2]) == ProductOfSimplicialSets((S2, S2))
            True
        """
        if factors:
            return super(ProductOfSimplicialSets, cls).__classcall__(cls, factors=tuple(factors))
        return super(ProductOfSimplicialSets, cls).__classcall__(cls)

    def __init__(self, factors=None):
        r"""
        Return the product of simplicial sets.

        INPUT:

        - ``factors`` -- a list or tuple of simplicial sets

        Return the product of the simplicial sets in ``factors``.

        If `X` and `Y` are simplicial sets, then their product `X
        \times Y` is defined to be the simplicial set with
        `n`-simplices `X_n \times Y_n`.  Therefore the simplices in
        the product have the form `(s_I \sigma, s_J \tau)`, where `s_I
        = s_{i_1} ... s_{i_p}` and `s_J = s_{j_1} ... s_{j_q}` are
        composites of degeneracy maps, written in decreasing order.
        Such a simplex is nondegenerate if the indices `I` and `J` are
        disjoint. Therefore if `\sigma` and `\tau` are nondegenerate
        simplices of dimensions `m` and `n`, in the product they will
        lead to nondegenerate simplices up to dimension `m+n`, and no
        further.

        This extends in the more or less obvious way to products with
        more than two factors: with three factors, a simplex `(s_I
        \sigma, s_J \tau, s_K \rho)` is nondegenerate if `I \cap J
        \cap K` is empty, etc.

        If a simplicial set is constructed as a product, the factors
        are recorded and are accessible via the method
        :meth:`Factors.factors`. If it is constructed as a product and then
        copied, this information is lost.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: e = AbstractSimplex(1, name='e')
            sage: X = SimplicialSet({e: (v, w)})
            sage: square = X.product(X)

        ``square`` is now the standard triangulation of the square: 4
        vertices, 5 edges (the four on the border plus the diagonal),
        2 triangles::

            sage: square.f_vector()
            [4, 5, 2]

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = S1.product(S1)
            sage: T.homology(reduced=False)
            {0: Z, 1: Z x Z, 2: Z}

        Since ``S1`` is pointed, so is ``T``::

            sage: S1.is_pointed()
            True
            sage: S1.base_point()
            v_0
            sage: T.is_pointed()
            True
            sage: T.base_point()
            (v_0, v_0)

            sage: S2 = simplicial_sets.Sphere(2)
            sage: S3 = simplicial_sets.Sphere(3)
            sage: Z = S2.product(S3)
            sage: Z.homology()
            {0: 0, 1: 0, 2: Z, 3: Z, 4: 0, 5: Z}

        Products involving infinite simplicial sets::

            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: B.rename('RP^oo')
            sage: X = B.product(B)
            sage: X
            RP^oo x RP^oo
            sage: X.n_cells(1)
            [(f, f), (f, s_0 1), (s_0 1, f)]
            sage: X.homology(range(3), base_ring=GF(2))
            {0: Vector space of dimension 0 over Finite Field of size 2,
             1: Vector space of dimension 2 over Finite Field of size 2,
             2: Vector space of dimension 3 over Finite Field of size 2}
            sage: Y = B.product(S2)
            sage: Y.homology(range(5), base_ring=GF(2))
            {0: Vector space of dimension 0 over Finite Field of size 2,
             1: Vector space of dimension 1 over Finite Field of size 2,
             2: Vector space of dimension 2 over Finite Field of size 2,
             3: Vector space of dimension 2 over Finite Field of size 2,
             4: Vector space of dimension 2 over Finite Field of size 2}
        """
        PullbackOfSimplicialSets.__init__(self, [space.constant_map()
                                                 for space in factors])
        self._factors = factors

    def n_skeleton(self, n):
        """
        Return the `n`-skeleton of this simplicial set.

        That is, the simplicial set generated by all nondegenerate
        simplices of dimension at most `n`.

        INPUT:

        - ``n`` -- the dimension

        In the finite case, this returns the ordinary `n`-skeleton. In
        the infinite case, it computes the `n`-skeleton of the product
        of the `n`-skeleta of the factors.

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: S3 = simplicial_sets.Sphere(3)
            sage: S2.product(S3).n_skeleton(2)
            Simplicial set with 2 non-degenerate simplices
            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: X = B.product(B)
            sage: X.n_skeleton(2)
            Simplicial set with 13 non-degenerate simplices
        """
        n_skel = SimplicialSet_finite.n_skeleton
        if self.is_finite():
            n_skel = SimplicialSet_finite.n_skeleton
            return n_skel(self, n)
        start, skel = self._n_skeleton
        if start == n:
            return skel
        elif start > n:
            return skel.n_skeleton(n)
        ans = n_skel(ProductOfSimplicialSets_finite([X.n_skeleton(n) for X in self._factors]), n)
        self._n_skeleton = (n, ans)
        return ans

    def factor(self, i, as_subset=False):
        r"""
        Return the $i$-th factor of the product.

        INPUT:

        - ``i`` -- integer, the index of the factor

        - ``as_subset`` -- boolean, optional (default ``False``)

        If ``as_subset`` is ``True``, return the $i$-th factor as a
        subsimplicial set of the product, identifying it with its
        product with the base point in each other factor. As a
        subsimplicial set, it comes equipped with an inclusion
        map. This option will raise an error if any factor does not
        have a base point.

        If ``as_subset`` is ``False``, return the $i$-th factor in
        its original form as a simplicial set.

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: S3 = simplicial_sets.Sphere(3)
            sage: K = S2.product(S3)
            sage: K.factor(0)
            S^2

            sage: K.factor(0, as_subset=True)
            Simplicial set with 2 non-degenerate simplices
            sage: K.factor(0, as_subset=True).homology()
            {0: 0, 1: 0, 2: Z}

            sage: K.factor(0) is S2
            True
            sage: K.factor(0, as_subset=True) is S2
            False
        """
        if as_subset:
            if any(not _.is_pointed() for _ in self.factors()):
                raise ValueError('"as_subset=True" is only valid '
                                 'if each factor is pointed')

            basept_factors = [sset.base_point() for sset in self.factors()]
            basept_factors = basept_factors[:i] + basept_factors[i+1:]
            to_factors = dict((v,k) for k,v in self._translation)
            simps = []
            for x in self.nondegenerate_simplices():
                simplices = [sigma[0] for sigma in to_factors[x]]
                if simplices[:i] + simplices[i+1:] == basept_factors:
                    simps.append(x)
            return self.subsimplicial_set(simps)
        return self.factors()[i]

    def _repr_(self):
        """
        Print representation

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: K = simplicial_sets.KleinBottle()
            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: S2.product(S2)
            S^2 x S^2
            sage: S2.product(K, B)
            S^2 x Klein bottle x Classifying space of Multiplicative Abelian group isomorphic to C2
        """
        return ' x '.join([str(X) for X in self._factors])

    def _latex_(self):
        r"""
        LaTeX representation

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: latex(S2.product(S2))
            S^{2} \times S^{2}
            sage: RPoo = simplicial_sets.RealProjectiveSpace(Infinity)
            sage: latex(S2.product(RPoo, S2))
            S^{2} \times RP^{\infty} \times S^{2}
        """
        return ' \\times '.join([latex(X) for X in self._factors])


class ProductOfSimplicialSets_finite(ProductOfSimplicialSets, PullbackOfSimplicialSets_finite):
    r"""
    The product of finite simplicial sets.

    When the factors are all finite, there are more methods available
    for the resulting product, as compared to products with infinite
    factors: projection maps, the wedge as a subcomplex, and the fat
    wedge as a subcomplex. See :meth:`projection_map`,
    :meth:`wedge_as_subset`, and :meth:`fat_wedge_as_subset`
    """
    def __init__(self, factors=None):
        r"""
        Return the product of finite simplicial sets.

        See :class:`ProductOfSimplicialSets` for more information.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(1)
            sage: X = SimplicialSet({e: (v, v)})
            sage: W = X.product(X, X)
            sage: W.homology()
            {0: 0, 1: Z x Z x Z, 2: Z x Z x Z, 3: Z}
            sage: W.is_pointed()
            False

            sage: X = X.set_base_point(v)
            sage: w = AbstractSimplex(0, name='w')
            sage: f = AbstractSimplex(1)
            sage: Y = SimplicialSet({f: (v,w)}, base_point=w)
            sage: Z = Y.product(X)
            sage: Z.is_pointed()
            True
            sage: Z.base_point()
            (w, v)
        """
        PullbackOfSimplicialSets_finite.__init__(self, [space.constant_map()
                                                 for space in factors])
        self._factors = tuple([f.domain() for f in self._maps])

    def projection_map(self, i):
        """
        Return the map projecting onto the $i$-th factor.

        INPUT:

        - ``i`` -- integer, the index of the projection map

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: f_0 = T.projection_map(0)
            sage: f_1 = T.projection_map(1)
            sage: m_0 = f_0.induced_homology_morphism().to_matrix(1) # matrix in dim 1
            sage: m_1 = f_1.induced_homology_morphism().to_matrix(1)
            sage: m_0.rank()
            1
            sage: m_0 == m_1
            False
        """
        return self.structure_map(i)

    def wedge_as_subset(self):
        """
        Return the wedge as a subsimplicial set of this product of pointed
        simplicial sets.

        This will raise an error if any factor is not pointed.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(1, name='e')
            sage: w = AbstractSimplex(0, name='w')
            sage: f = AbstractSimplex(1, name='f')
            sage: X = SimplicialSet({e: (v, v)}, base_point=v)
            sage: Y = SimplicialSet({f: (w, w)}, base_point=w)
            sage: P = X.product(Y)
            sage: W = P.wedge_as_subset()
            sage: W.nondegenerate_simplices()
            [(v, w), (e, s_0 w), (s_0 v, f)]
            sage: W.homology()
            {0: 0, 1: Z x Z}
        """
        basept_factors = [sset.base_point() for sset in self.factors()]
        to_factors = dict((v,k) for k,v in self._translation)
        simps = []
        for x in self.nondegenerate_simplices():
            simplices = to_factors[x]
            not_base_pt = 0
            for sigma, star in zip(simplices, basept_factors):
                if not_base_pt > 1:
                    continue
                if sigma[0].nondegenerate() != star:
                    not_base_pt += 1
            if not_base_pt <= 1:
                simps.append(x)
        return self.subsimplicial_set(simps)

    def fat_wedge_as_subset(self):
        """
        Return the fat wedge as a subsimplicial set of this product of
        pointed simplicial sets.

        The fat wedge consists of those terms where at least one
        factor is the base point. Thus with two factors this is the
        ordinary wedge, but with more factors, it is larger.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: X = S1.product(S1, S1)
            sage: W = X.fat_wedge_as_subset()
            sage: W.homology()
            {0: 0, 1: Z x Z x Z, 2: Z x Z x Z}
        """
        basept_factors = [sset.base_point() for sset in self.factors()]
        to_factors = {v: k for k, v in self._translation}
        simps = []
        for x in self.nondegenerate_simplices():
            simplices = to_factors[x]
            combined = zip(simplices, basept_factors)
            if any(sigma[0] == pt for (sigma, pt) in combined):
                simps.append(x)
        return self.subsimplicial_set(simps)


class PushoutOfSimplicialSets(SimplicialSet_arbitrary, UniqueRepresentation):
    @staticmethod
    def __classcall_private__(cls, maps=None, vertex_name=None):
        """
        TESTS::

            sage: from sage.topology.simplicial_set_constructions import PushoutOfSimplicialSets
            sage: S2 = simplicial_sets.Sphere(2)
            sage: one = S2.Hom(S2).identity()
            sage: PushoutOfSimplicialSets([one, one]) == PushoutOfSimplicialSets((one, one))
            True
        """
        if maps:
            return super(PushoutOfSimplicialSets, cls).__classcall__(cls, maps=tuple(maps),
                                                                      vertex_name=vertex_name)
        return super(PushoutOfSimplicialSets, cls).__classcall__(cls, vertex_name=vertex_name)

    def __init__(self, maps=None, vertex_name=None):
        r"""
        Return the pushout obtained from the morphisms ``maps``.

        INPUT:

        - ``maps`` -- a list or tuple of morphisms of simplicial sets
        - ``vertex_name`` -- optional, default ``None``

        If only a single map `f: X \to Y` is given, then return
        `Y`. If no maps are given, return the empty simplicial
        set. Otherwise, given a simplicial set `X` and maps `f_i: X
        \to Y_i` for `0 \leq i \leq m`, construct the pushout `P`: see
        :wikipedia:`Pushout_(category_theory)`. This is constructed as
        pushouts of sets for each set of `n`-simplices, so `P_n` is
        the disjoint union of the sets `(Y_i)_n`, with elements
        `f_i(x)` identified for `n`-simplex `x` in `X`.

        Simplices in the pushout are given names as follows: if a
        simplex comes from a single `Y_i`, it inherits its
        name. Otherwise it must come from a simplex (or several) in
        `X`, and then it inherits one of those names, and it should be
        the first alphabetically. For example, if vertices `v`, `w`,
        and `z` in `X` are glued together, then the resulting vertex
        in the pushout will be called `v`.

        Base points are taken care of automatically: if each of the
        maps `f_i` is pointed, so is the pushout. If `X` is a point or
        if `X` is nonempty and any of the spaces `Y_i` is a point, use
        those for the base point. In all of these cases, if
        ``vertex_name`` is ``None``, generate the name of the base
        point automatically; otherwise, use ``vertex_name`` for its
        name.

        In all other cases, the pushout is not pointed.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: a = AbstractSimplex(0, name='a')
            sage: b = AbstractSimplex(0, name='b')
            sage: c = AbstractSimplex(0, name='c')
            sage: e0 = AbstractSimplex(1, name='e_0')
            sage: e1 = AbstractSimplex(1, name='e_1')
            sage: e2 = AbstractSimplex(1, name='e_2')
            sage: X = SimplicialSet({e2: (b, a)})
            sage: Y0 = SimplicialSet({e2: (b,a), e0: (c,b), e1: (c,a)})
            sage: Y1 = simplicial_sets.Simplex(0)
            sage: f0_data = {a:a, b:b, e2: e2}
            sage: v = Y1.n_cells(0)[0]
            sage: f1_data = {a:v, b:v, e2:v.apply_degeneracies(0)}
            sage: f0 = X.Hom(Y0)(f0_data)
            sage: f1 = X.Hom(Y1)(f1_data)
            sage: P = X.pushout(f0, f1)
            sage: P.nondegenerate_simplices()
            [a, c, e_0, e_1]

        There are defining maps `f_i: X \to Y_i` and structure maps
        `\bar{f}_i: Y_i \to P`; the latter are only implemented in
        Sage when each `Y_i` is finite. ::

            sage: P.defining_map(0) == f0
            True
            sage: P.structure_map(1)
            Simplicial set morphism:
              From: 0-simplex
              To:   Pushout of maps:
              Simplicial set morphism:
                From: Simplicial set with 3 non-degenerate simplices
                To:   Simplicial set with 6 non-degenerate simplices
                Defn: [a, b, e_2] --> [a, b, e_2]
              Simplicial set morphism:
                From: Simplicial set with 3 non-degenerate simplices
                To:   0-simplex
                Defn: Constant map at (0,)
              Defn: Constant map at a
            sage: P.structure_map(0).domain() == Y0
            True
            sage: P.structure_map(0).codomain() == P
            True

        An inefficient way of constructing a suspension for an
        unpointed set: take the pushout of two copies of the inclusion
        map `X \to CX`::

            sage: T = simplicial_sets.Torus()
            sage: T = T.unset_base_point()
            sage: CT = T.cone()
            sage: inc = CT.base_as_subset().inclusion_map()
            sage: P = T.pushout(inc, inc)
            sage: P.homology()
            {0: 0, 1: 0, 2: Z x Z, 3: Z}
            sage: len(P.nondegenerate_simplices())
            20

        It is more efficient to construct the suspension as the
        quotient `CX/X`::

            sage: len(CT.quotient(CT.base_as_subset()).nondegenerate_simplices())
            8

        It is more efficient still if the original simplicial set has
        a base point::

            sage: T = simplicial_sets.Torus()
            sage: len(T.suspension().nondegenerate_simplices())
            6

            sage: S1 = simplicial_sets.Sphere(1)
            sage: pt = simplicial_sets.Point()
            sage: bouquet = pt.pushout(S1.base_point_map(), S1.base_point_map(), S1.base_point_map())
            sage: bouquet.homology(1)
            Z x Z x Z
        """
        # Import this here to prevent circular imports.
        from sage.topology.simplicial_set_morphism import SimplicialSetMorphism
        if maps and any(not isinstance(f, SimplicialSetMorphism) for f in maps):
            raise ValueError('the maps must be morphisms of simplicial sets')
        Cat = SimplicialSets()
        if maps:
            if all(f.codomain().is_finite() for f in maps):
                Cat = Cat.Finite()
            if all(f.is_pointed() for f in maps):
                Cat = Cat.Pointed()
            Parent.__init__(self, category=Cat)
        self._maps = maps
        self._n_skeleton = (-1, Empty())
        self._vertex_name = vertex_name

    def n_skeleton(self, n):
        """
        Return the `n`-skeleton of this simplicial set.

        That is, the simplicial set generated by all nondegenerate
        simplices of dimension at most `n`.

        INPUT:

        - ``n`` -- the dimension

        The `n`-skeleton of the pushout is computed as the pushout
        of the `n`-skeleta of the component simplicial sets.

        EXAMPLES::

            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: K = B.n_skeleton(3)
            sage: Q = K.pushout(K.inclusion_map(), K.constant_map())
            sage: Q.n_skeleton(5).homology()
            {0: 0, 1: 0, 2: 0, 3: 0, 4: Z, 5: Z}

        Of course, computing the `n`-skeleton and then taking homology
        need not yield the same answer as asking for homology through
        dimension `n`, since the latter computation will use the
        `(n+1)`-skeleton::

            sage: Q.homology(range(6))
            {0: 0, 1: 0, 2: 0, 3: 0, 4: Z, 5: C2}
        """
        if self.is_finite():
            maps = self._maps
            if maps:
                domain = SimplicialSet_finite.n_skeleton(maps[0].domain(), n)
                codomains = [SimplicialSet_finite.n_skeleton(f.codomain(), n) for f in maps]
                new_maps = [f.n_skeleton(n, domain, c) for (f, c) in zip(maps, codomains)]
                return PushoutOfSimplicialSets_finite(new_maps,
                                                      vertex_name=self._vertex_name)
            return PushoutOfSimplicialSets_finite(maps,
                                                  vertex_name=self._vertex_name)
        start, skel = self._n_skeleton
        if start == n:
            return skel
        elif start > n:
            return skel.n_skeleton(n)
        ans = PushoutOfSimplicialSets_finite([f.n_skeleton(n) for f in self._maps],
                                             vertex_name=self._vertex_name)
        self._n_skeleton = (n, ans)
        return ans

    def defining_map(self, i):
        r"""
        Return the `i`-th map defining the pushout.

        INPUT:

        - ``i`` -- integer

        If this pushout was constructed as ``X.pushout(f_0, f_1, ...)``,
        this returns `f_i`.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = simplicial_sets.Torus()
            sage: X = S1.wedge(T) # a pushout
            sage: X.defining_map(0)
            Simplicial set morphism:
              From: Point
              To:   S^1
              Defn: Constant map at v_0
            sage: X.defining_map(1).domain()
            Point
            sage: X.defining_map(1).codomain()
            Torus
        """
        return self._maps[i]

    def _repr_(self):
        """
        Print representation

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: S3 = simplicial_sets.Sphere(3)
            sage: pt = simplicial_sets.Point()
            sage: pt.pushout(S2.base_point_map(), S3.base_point_map())
            Pushout of maps:
              Simplicial set morphism:
                From: Point
                To:   S^2
                Defn: Constant map at v_0
              Simplicial set morphism:
                From: Point
                To:   S^3
                Defn: Constant map at v_0
        """
        if not self._maps:
            return 'Empty simplicial set'
        s = 'Pushout of maps:'
        for f in self._maps:
            t = '\n' + str(f)
            s += t.replace('\n', '\n  ')
        return s


class PushoutOfSimplicialSets_finite(PushoutOfSimplicialSets, SimplicialSet_finite):
    """
    The pushout of finite simplicial sets obtained from ``maps``.

    When the simplicial sets involved are all finite, there are more
    methods available to the resulting pushout, as compared to case
    when some of the components are infinite: the structure maps to the
    pushout and the pushout's universal property: see
    :meth:`structure_map` and :meth:`universal_property`.
    """
    @staticmethod
    def __classcall_private__(cls, maps=None, vertex_name=None):
        """
        TESTS::

            sage: from sage.topology.simplicial_set_constructions import PushoutOfSimplicialSets_finite
            sage: S2 = simplicial_sets.Sphere(2)
            sage: one = S2.Hom(S2).identity()
            sage: PushoutOfSimplicialSets_finite([one, one]) == PushoutOfSimplicialSets_finite((one, one))
            True
        """
        if maps:
            return super(PushoutOfSimplicialSets_finite, cls).__classcall__(cls, maps=tuple(maps),
                                                                      vertex_name=vertex_name)
        return super(PushoutOfSimplicialSets_finite, cls).__classcall__(cls, vertex_name=vertex_name)

    def __init__(self, maps=None, vertex_name=None):
        r"""
        Return the pushout obtained from the morphisms ``maps``.

        See :class:`PushoutOfSimplicialSets` for more information.

        INPUT:

        - ``maps`` -- a list or tuple of morphisms of simplicial sets
        - ``vertex_name`` -- optional, default ``None``

        EXAMPLES::

            sage: from sage.topology.simplicial_set_constructions import PushoutOfSimplicialSets_finite
            sage: T = simplicial_sets.Torus()
            sage: S2 = simplicial_sets.Sphere(2)
            sage: PushoutOfSimplicialSets_finite([T.base_point_map(), S2.base_point_map()]).n_cells(0)[0]
            *
            sage: PushoutOfSimplicialSets_finite([T.base_point_map(), S2.base_point_map()], vertex_name='v').n_cells(0)[0]
            v
        """
        # Import this here to prevent circular imports.
        from sage.topology.simplicial_set_morphism import SimplicialSetMorphism
        if maps and any(not isinstance(f, SimplicialSetMorphism) for f in maps):
            raise ValueError('the maps must be morphisms of simplicial sets')
        if not maps:
            SimplicialSet_finite.__init__(self, {})
            self._maps = ()
            self._structure = ()
            return
        domain = maps[0].domain()
        if len(maps) == 1:
            # f: X --> Y
            f = maps[0]
            codomain = f.codomain()
            if f.is_pointed():
                base_point=codomain.base_point()
                if vertex_name is not None:
                    base_point.rename(vertex_name)
                SimplicialSet_finite.__init__(self, codomain.face_data(),
                                       base_point=base_point)
            elif len(domain.nondegenerate_simplices()) == 1:
                # X is a point.
                base_point = f(domain().n_cells(0)[0])
                if vertex_name is not None:
                    base_point.rename(vertex_name)
                SimplicialSet_finite.__init__(self, codomain.face_data(),
                                       base_point=base_point)
            elif len(codomain.nondegenerate_simplices()) == 1:
                # Y is a point.
                base_point = codomain.n_cells(0)[0]
                if vertex_name is not None:
                    base_point.rename(vertex_name)
                SimplicialSet_finite.__init__(self, codomain.face_data(),
                                       base_point=base_point)
            else:
                SimplicialSet_finite.__init__(self, codomain.face_data())
            self._maps = (f,)
            self._structure = (f,)
            return
        if any(domain != f.domain() for f in maps[1:]):
            raise ValueError('the domains of the maps must be equal')
        # Data to define the pushout:
        data = {}
        codomains = [f.codomain() for f in maps]
        # spaces: indexed list of spaces. Entries are of the form
        # (space, int) where int=-1 for the domain, and for the
        # codomains, int is the corresponding index.
        spaces = [(Y,i-1) for (i,Y) in enumerate([domain] + codomains)]
        # Dictionaries to translate from simplices in domain,
        # codomains to simplices in the pushout. The keys are of the
        # form (space, int). int=-1 for the domain, and for the
        # codomains, int is the corresponding index.
        _to_P = {Y:{} for Y in spaces}
        max_dim = max(Y.dimension() for Y in codomains)
        for n in range(1 + max_dim):
            # Now we impose an equivalence relation on the simplices,
            # setting x equivalent to f_i(x) for each simplex x in X
            # and each defining map f_i. We do this by constructing a
            # graph and finding its connected components: the vertices
            # of the graph are the n-cells of X and the Y_i, and
            # there are edges from x to f_i(x).
            vertices = []
            for (Y,i) in spaces:
                vertices.extend([(cell,i) for cell in Y.n_cells(n)])
            edges = []
            for x in domain.n_cells(n):
                edges.extend([[(x,-1), (f(x),i)] for (i,f) in enumerate(maps)])
            G = Graph([vertices, edges], format='vertices_and_edges')
            data[n] = [set(_) for _ in G.connected_components()]
        # data is now a dictionary indexed by dimension, and data[n]
        # consists of sets of n-simplices of the domain and the
        # codomains, each set an equivalence class of n-simplices
        # under the gluing. So if any element of one of those sets is
        # degenerate, we can throw the whole thing away. Otherwise, we
        # can choose a representative to compute the faces.
        simplices = {}
        for dim in sorted(data):
            for s in data[dim]:
                degenerate = any(sigma[0].is_degenerate() for sigma in s)
                if degenerate:
                    # Identify the degeneracies involved.
                    degens = []
                    for (sigma, j) in s:
                        if len(sigma.degeneracies()) > len(degens):
                            degens = sigma.degeneracies()
                            space = spaces[j+1]
                            old = _to_P[space][sigma.nondegenerate()]
                    for (sigma,j) in s:
                        # Now update the _to_P[space] dictionaries.
                        space = spaces[j+1]
                        _to_P[space][sigma] = old.apply_degeneracies(*degens)
                else: # nondegenerate
                    if len(s) == 1:
                        name = str(list(s)[0][0])
                        latex_name = latex(list(s)[0][0])
                    else:
                        # Choose a name from a simplex in domain.
                        for (sigma,j) in sorted(s):
                            if j == -1:
                                name = str(sigma)
                                latex_name = latex(sigma)
                                break
                    new = AbstractSimplex(dim, name=name,
                                          latex_name=latex_name)
                    if dim == 0:
                        faces = None
                    for (sigma,j) in s:
                        space = spaces[j+1]
                        _to_P[space][sigma] = new
                        if dim > 0:
                            faces = [_to_P[space][tau.nondegenerate()].apply_degeneracies(*tau.degeneracies())
                                     for tau in space[0].faces(sigma)]
                    simplices[new] = faces

        some_Y_is_pt = False
        if len(domain.nondegenerate_simplices()) > 1:
            # Only investigate this if X is not empty and not a point.
            for (Y,i) in spaces:
                if len(Y.nondegenerate_simplices()) == 1:
                    some_Y_is_pt = True
                    break
        if len(domain.nondegenerate_simplices()) == 1:
            # X is a point.
            base_point = _to_P[(domain,-1)][domain.n_cells(0)[0]]
            if vertex_name is not None:
                base_point.rename(vertex_name)
            SimplicialSet_finite.__init__(self, simplices, base_point=base_point)
        elif some_Y_is_pt:
            # We found (Y,i) above.
            base_point = _to_P[(Y,i)][Y.n_cells(0)[0]]
            if vertex_name is not None:
                base_point.rename(vertex_name)
            SimplicialSet_finite.__init__(self, simplices, base_point=base_point)
        elif all(f.is_pointed() for f in maps):
            pt = _to_P[(codomains[0],0)][codomains[0].base_point()]
            if any(_to_P[(Y,i)][Y.base_point()] != pt for (Y,i) in spaces[2:]):
                raise ValueError('something unexpected went wrong '
                                 'with base points')
            base_point = _to_P[(domain,-1)][domain.base_point()]
            if vertex_name is not None:
                base_point.rename(vertex_name)
            SimplicialSet_finite.__init__(self, simplices, base_point=base_point)
        else:
            SimplicialSet_finite.__init__(self, simplices)
        # The relevant maps:
        self._maps = maps
        self._structure = tuple([Y.Hom(self)(_to_P[(Y,i)])
                               for (Y,i) in spaces[1:]])
        self._vertex_name = vertex_name

    def structure_map(self, i):
        r"""
        Return the $i$-th structure map of the pushout.

        INPUT:

        - ``i`` -- integer

        If this pushout `Z` was constructed as ``X.pushout(f_0, f_1, ...)``,
        where `f_i: X \to Y_i`, then there are structure maps
        `\bar{f}_i: Y_i \to Z`. This method constructs `\bar{f}_i`.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = simplicial_sets.Torus()
            sage: X = S1.disjoint_union(T) # a pushout
            sage: X.structure_map(0)
            Simplicial set morphism:
              From: S^1
              To:   Disjoint union: (S^1 u Torus)
              Defn: [v_0, sigma_1] --> [v_0, sigma_1]
            sage: X.structure_map(1).domain()
            Torus
            sage: X.structure_map(1).codomain()
            Disjoint union: (S^1 u Torus)
        """
        return self._structure[i]

    def universal_property(self, *maps):
        r"""
        Return the map induced by ``maps``

        INPUT:

        - ``maps`` -- maps "factors" `Y_i` forming the pushout to a
          fixed simplicial set `Z`.

        If the pushout `P` is formed by maps `f_i: X \to Y_i`, then
        given maps `g_i: Y_i \to Z` such that `g_i f_i = g_j f_j` for
        all `i`, `j`, then there is a unique map `g: P \to Z` making
        the appropriate diagram commute. This constructs that map.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: x = AbstractSimplex(0, name='x')
            sage: evw = AbstractSimplex(1, name='vw')
            sage: evx = AbstractSimplex(1, name='vx')
            sage: ewx = AbstractSimplex(1, name='wx')
            sage: X = SimplicialSet({evw: (w, v), evx: (x, v)})
            sage: Y_0 = SimplicialSet({evw: (w, v), evx: (x, v), ewx: (x, w)})
            sage: Y_1 = SimplicialSet({evx: (x, v)})

            sage: f_0 = Hom(X, Y_0)({v:v, w:w, x:x, evw:evw, evx:evx})
            sage: f_1 = Hom(X, Y_1)({v:v, w:v, x:x, evw:v.apply_degeneracies(0), evx:evx})
            sage: P = X.pushout(f_0, f_1)

            sage: one = Hom(Y_1, Y_1).identity()
            sage: g = Hom(Y_0, Y_1)({v:v, w:v, x:x, evw:v.apply_degeneracies(0), evx:evx, ewx:evx})
            sage: P.universal_property(g, one)
            Simplicial set morphism:
              From: Pushout of maps:
              Simplicial set morphism:
                From: Simplicial set with 5 non-degenerate simplices
                To:   Simplicial set with 6 non-degenerate simplices
                Defn: [v, w, x, vw, vx] --> [v, w, x, vw, vx]
              Simplicial set morphism:
                From: Simplicial set with 5 non-degenerate simplices
                To:   Simplicial set with 3 non-degenerate simplices
                Defn: [v, w, x, vw, vx] --> [v, v, x, s_0 v, vx]
              To:   Simplicial set with 3 non-degenerate simplices
              Defn: [v, x, vx, wx] --> [v, x, vx, vx]
        """
        codomain = maps[0].codomain()
        if any(g.codomain() != codomain for g in maps[1:]):
            raise ValueError('the maps do not all have the same codomain')
        composite = maps[0] * self._maps[0]
        if any(g*f != composite for g,f in zip(maps[1:], self._maps[1:])):
            raise ValueError('the maps are not compatible')
        data = {}
        for i,g in enumerate(maps):
            f_i_dict = self.structure_map(i)._dictionary
            for sigma in f_i_dict:
                tau = f_i_dict[sigma]
                # For sigma_i in Y_i, define the map G by
                # G(\bar{f}_i)(sigma_i) = g_i(sigma_i).
                if tau not in data:
                    data[tau] = g(sigma)
        return self.Hom(codomain)(data)


class QuotientOfSimplicialSet(PushoutOfSimplicialSets):
    def __init__(self, inclusion, vertex_name='*'):
        r"""
        Return the quotient of a simplicial set by a subsimplicial set.

        INPUT:

        - ``inclusion`` -- inclusion map of a subcomplex (=
          subsimplicial set) of a simplicial set
        - ``vertex_name`` -- optional, default ``'*'``

        A subcomplex `A` comes equipped with the inclusion map `A \to
        X` to its ambient complex `X`, and this constructs the
        quotient `X/A`, collapsing `A` to a point. The resulting point
        is called ``vertex_name``, which is ``'*'`` by default.

        When the simplicial sets involved are finite, there is a
        :meth:`QuotientOfSimplicialSet_finite.quotient_map` method available.

        EXAMPLES::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP2 = RP5.n_skeleton(2)
            sage: RP5_2 = RP5.quotient(RP2)
            sage: RP5_2
            Quotient: (RP^5/Simplicial set with 3 non-degenerate simplices)
            sage: RP5_2.quotient_map()
            Simplicial set morphism:
              From: RP^5
              To:   Quotient: (RP^5/Simplicial set with 3 non-degenerate simplices)
              Defn: [1, f, f * f, f * f * f, f * f * f * f, f * f * f * f * f] --> [*, s_0 *, s_1 s_0 *, f * f * f, f * f * f * f, f * f * f * f * f]
        """
        subcomplex = inclusion.domain()
        PushoutOfSimplicialSets.__init__(self, [inclusion,
                                                subcomplex.constant_map()],
                                         vertex_name=vertex_name)

        ambient = inclusion.codomain()
        if ambient.is_pointed() and ambient.is_finite():
            if ambient.base_point() not in subcomplex:
                self._basepoint = self.structure_map(0)(ambient.base_point())

    def ambient(self):
        """
        Return the ambient space.

        That is, if this quotient is `K/L`, return `K`.

        EXAMPLES::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP2 = RP5.n_skeleton(2)
            sage: RP5_2 = RP5.quotient(RP2)
            sage: RP5_2.ambient()
            RP^5

            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: K = B.n_skeleton(3)
            sage: Q = B.quotient(K)
            sage: Q.ambient()
            Classifying space of Multiplicative Abelian group isomorphic to C2
        """
        return self._maps[0].codomain()

    def subcomplex(self):
        """
        Return the subcomplex space associated to this quotient.

        That is, if this quotient is `K/L`, return `L`.

        EXAMPLES::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP2 = RP5.n_skeleton(2)
            sage: RP5_2 = RP5.quotient(RP2)
            sage: RP5_2.subcomplex()
            Simplicial set with 3 non-degenerate simplices

            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: K = B.n_skeleton(3)
            sage: Q = B.quotient(K)
            sage: Q.subcomplex()
            Simplicial set with 4 non-degenerate simplices
        """
        return self._maps[0].domain()

    def n_skeleton(self, n):
        """
        Return the `n`-skeleton of this simplicial set.

        That is, the simplicial set generated by all nondegenerate
        simplices of dimension at most `n`.

        INPUT:

        - ``n`` -- the dimension

        The `n`-skeleton of the quotient is computed as the quotient
        of the `n`-skeleta.

        EXAMPLES::

            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: K = B.n_skeleton(3)
            sage: Q = B.quotient(K)
            sage: Q.n_skeleton(6)
            Quotient: (Simplicial set with 7 non-degenerate simplices/Simplicial set with 4 non-degenerate simplices)
            sage: Q.n_skeleton(6).homology()
            {0: 0, 1: 0, 2: 0, 3: 0, 4: Z, 5: C2, 6: 0}
        """
        if self.is_finite():
            ambient = SimplicialSet_finite.n_skeleton(self.ambient(), n)
            subcomplex = SimplicialSet_finite.n_skeleton(self.subcomplex(), n)
            subcomplex = ambient.subsimplicial_set(subcomplex.nondegenerate_simplices())
            return QuotientOfSimplicialSet_finite(subcomplex.inclusion_map(),
                                                  vertex_name=self._vertex_name)
        start, skel = self._n_skeleton
        if start == n:
            return skel
        elif start > n:
            return skel.n_skeleton(n)
        ambient = self.ambient().n_skeleton(n)
        subcomplex = ambient.subsimplicial_set(self.subcomplex().nondegenerate_simplices(n))
        ans = QuotientOfSimplicialSet_finite(subcomplex.inclusion_map(),
                                             vertex_name=self._vertex_name)
        self._n_skeleton = (n, ans)
        return ans

    def _repr_(self):
        """
        Print representation

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: T.quotient(T.n_skeleton(1))
            Quotient: (Torus/Simplicial set with 4 non-degenerate simplices)
        """
        return 'Quotient: ({}/{})'.format(self.ambient(), self.subcomplex())

    def _latex_(self):
        r"""
        LaTeX representation

        EXAMPLES::

            sage: RPoo = simplicial_sets.RealProjectiveSpace(Infinity)
            sage: RP3 = RPoo.n_skeleton(3)
            sage: RP3.rename_latex('RP^{3}')
            sage: latex(RPoo.quotient(RP3))
            RP^{\infty} / RP^{3}
        """
        return '{} / {}'.format(latex(self.ambient()), latex(self.subcomplex()))


class QuotientOfSimplicialSet_finite(QuotientOfSimplicialSet,
                                     PushoutOfSimplicialSets_finite):
    """
    The quotient of finite simplicial sets.

    When the simplicial sets involved are finite, there is a
    :meth:`quotient_map` method available.
    """
    def __init__(self, inclusion, vertex_name='*'):
        r"""
        Return the quotient of a simplicial set by a subsimplicial set.

        See :class:`QuotientOfSimplicialSet` for more information.

        EXAMPLES::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP2 = RP5.n_skeleton(2)
            sage: RP5_2 = RP5.quotient(RP2)
            sage: RP5_2
            Quotient: (RP^5/Simplicial set with 3 non-degenerate simplices)
            sage: RP5_2.quotient_map()
            Simplicial set morphism:
              From: RP^5
              To:   Quotient: (RP^5/Simplicial set with 3 non-degenerate simplices)
              Defn: [1, f, f * f, f * f * f, f * f * f * f, f * f * f * f * f] --> [*, s_0 *, s_1 s_0 *, f * f * f, f * f * f * f, f * f * f * f * f]
        """
        subcomplex = inclusion.domain()
        PushoutOfSimplicialSets_finite.__init__(self, [inclusion,
                                                       subcomplex.constant_map()],
                                                vertex_name=vertex_name)
        ambient = inclusion.codomain()
        if ambient.is_pointed():
            if ambient.base_point() not in subcomplex:
                self._basepoint = self.structure_map(0)(ambient.base_point())

    def quotient_map(self):
        """
        Return the quotient map from the original simplicial set to the
        quotient.

        EXAMPLES::

            sage: K = simplicial_sets.Simplex(1)
            sage: S1 = K.quotient(K.n_skeleton(0))
            sage: q = S1.quotient_map()
            sage: q
            Simplicial set morphism:
              From: 1-simplex
              To:   Quotient: (1-simplex/Simplicial set with 2 non-degenerate simplices)
              Defn: [(0,), (1,), (0, 1)] --> [*, *, (0, 1)]
            sage: q.domain() == K
            True
            sage: q.codomain() == S1
            True
        """
        return self.structure_map(0)


class SmashProductOfSimplicialSets_finite(QuotientOfSimplicialSet_finite,
                                          Factors):
    @staticmethod
    def __classcall__(cls, factors=None):
        """
        TESTS::

            sage: from sage.topology.simplicial_set_constructions import SmashProductOfSimplicialSets_finite as Smash
            sage: S2 = simplicial_sets.Sphere(2)
            sage: Smash([S2, S2]) == Smash((S2, S2))
            True
        """
        if factors:
            return super(SmashProductOfSimplicialSets_finite, cls).__classcall__(cls, factors=tuple(factors))
        return super(SmashProductOfSimplicialSets_finite, cls).__classcall__(cls)

    def __init__(self, factors=None):
        r"""
        Return the smash product of finite pointed simplicial sets.

        INPUT:

        - ``factors`` -- a list or tuple of simplicial sets

        Return the smash product of the simplicial sets in
        ``factors``: the smash product `X \wedge Y` is defined to be
        the quotient `(X \times Y) / (X \vee Y)`, where `X \vee Y` is
        the wedge sum.

        Each element of ``factors`` must be finite and pointed. (As of
        July 2016, constructing the wedge as a subcomplex of the
        product is only possible in Sage for finite simplicial sets.)

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: S2 = simplicial_sets.Sphere(2)
            sage: T.smash_product(S2).homology() == T.suspension(2).homology()
            True
        """
        if any(not space.is_pointed() for space in factors):
            raise ValueError('the simplicial sets must be pointed')
        prod = ProductOfSimplicialSets_finite(factors)
        wedge = prod.wedge_as_subset()
        QuotientOfSimplicialSet_finite.__init__(self, wedge.inclusion_map())
        self._factors = factors

    def _repr_(self):
        """
        Print representation

        EXAMPLES::

            sage: RP4 = simplicial_sets.RealProjectiveSpace(4)
            sage: S1 = simplicial_sets.Sphere(1)
            sage: S1.smash_product(RP4, S1)
            Smash product: (S^1 ^ RP^4 ^ S^1)
        """
        s = 'Smash product: ('
        s += ' ^ '.join([str(X) for X in self._factors])
        s += ')'
        return s

    def _latex_(self):
        r"""
        LaTeX representation

        EXAMPLES::

            sage: RP4 = simplicial_sets.RealProjectiveSpace(4)
            sage: S1 = simplicial_sets.Sphere(1)
            sage: latex(S1.smash_product(RP4, S1))
            S^{1} \wedge RP^{4} \wedge S^{1}
        """
        return ' \\wedge '.join([latex(X) for X in self._factors])


class WedgeOfSimplicialSets(PushoutOfSimplicialSets, Factors):
    @staticmethod
    def __classcall__(cls, factors=None):
        """
        TESTS::

            sage: from sage.topology.simplicial_set_constructions import WedgeOfSimplicialSets
            sage: S2 = simplicial_sets.Sphere(2)
            sage: WedgeOfSimplicialSets([S2, S2]) == WedgeOfSimplicialSets((S2, S2))
            True
        """
        if factors:
            return super(WedgeOfSimplicialSets, cls).__classcall__(cls, factors=tuple(factors))
        return super(WedgeOfSimplicialSets, cls).__classcall__(cls)

    def __init__(self, factors=None):
        r"""
        Return the wedge sum of pointed simplicial sets.

        INPUT:

        - ``factors`` -- a list or tuple of simplicial sets

        Return the wedge of the simplicial sets in ``factors``: the
        wedge sum `X \vee Y` is formed by taking the disjoint
        union of `X` and `Y` and identifying their base points. In
        this construction, the new base point is renamed '*'.

        The wedge comes equipped with maps to and from each factor, or
        actually, maps from each factor, and maps to simplicial sets
        isomorphic to each factor. The codomains of the latter maps
        are quotients of the wedge, not identical to the original
        factors.

        EXAMPLES::

            sage: CP2 = simplicial_sets.ComplexProjectiveSpace(2)
            sage: K = simplicial_sets.KleinBottle()
            sage: W = CP2.wedge(K)
            sage: W.homology()
            {0: 0, 1: Z x C2, 2: Z, 3: 0, 4: Z}

            sage: W.inclusion_map(1)
            Simplicial set morphism:
              From: Klein bottle
              To:   Wedge: (CP^2 v Klein bottle)
              Defn: [Delta_{0,0}, Delta_{1,0}, Delta_{1,1}, Delta_{1,2}, Delta_{2,0}, Delta_{2,1}] --> [*, Delta_{1,0}, Delta_{1,1}, Delta_{1,2}, Delta_{2,0}, Delta_{2,1}]

            sage: W.projection_map(0).domain()
            Wedge: (CP^2 v Klein bottle)
            sage: W.projection_map(0).codomain() # copy of CP^2
            Quotient: (Wedge: (CP^2 v Klein bottle)/Simplicial set with 6 non-degenerate simplices)
            sage: W.projection_map(0).codomain().homology()
            {0: 0, 1: 0, 2: Z, 3: 0, 4: Z}

        An error occurs if any of the factors is not pointed::

            sage: CP2.wedge(simplicial_sets.Simplex(1))
            Traceback (most recent call last):
            ...
            ValueError: the simplicial sets must be pointed
        """
        if any(not space.is_pointed() for space in factors):
            raise ValueError('the simplicial sets must be pointed')
        PushoutOfSimplicialSets.__init__(self, [space.base_point_map()
                                                for space in factors])
        if factors:
            vertices = PushoutOfSimplicialSets_finite([space.n_skeleton(0).base_point_map()
                                                       for space in factors])
            self._basepoint = vertices.base_point()
        self.base_point().rename('*')
        self._factors = factors

    summands = Factors.factors
    summand = Factors.factor

    def _repr_(self):
        """
        Print representation

        EXAMPLES::

            sage: K = simplicial_sets.KleinBottle()
            sage: K.wedge(K, K)
            Wedge: (Klein bottle v Klein bottle v Klein bottle)
        """
        s = 'Wedge: ('
        s += ' v '.join([str(X) for X in self._factors])
        s += ')'
        return s

    def _latex_(self):
        r"""
        LaTeX representation

        EXAMPLES::

            sage: RP4 = simplicial_sets.RealProjectiveSpace(4)
            sage: S1 = simplicial_sets.Sphere(1)
            sage: latex(S1.wedge(RP4, S1))
            S^{1} \vee RP^{4} \vee S^{1}
        """
        return ' \\vee '.join([latex(X) for X in self._factors])


class WedgeOfSimplicialSets_finite(WedgeOfSimplicialSets, PushoutOfSimplicialSets_finite):
    """
    The wedge sum of finite pointed simplicial sets.
    """
    def __init__(self, factors=None):
        r"""
        Return the wedge sum of finite pointed simplicial sets.

        INPUT:

        - ``factors`` -- a tuple of simplicial sets

        If there are no factors, a point is returned.

        See :class:`WedgeOfSimplicialSets` for more information.

        EXAMPLES::

            sage: from sage.topology.simplicial_set_constructions import WedgeOfSimplicialSets_finite
            sage: K = simplicial_sets.Simplex(3)
            sage: WedgeOfSimplicialSets_finite((K,K))
            Traceback (most recent call last):
            ...
            ValueError: the simplicial sets must be pointed
        """
        if not factors:
            # An empty wedge is a point, constructed as a pushout.
            PushoutOfSimplicialSets_finite.__init__(self, [Point().identity()])
        else:
            if any(not space.is_pointed() for space in factors):
                raise ValueError('the simplicial sets must be pointed')
            PushoutOfSimplicialSets_finite.__init__(self, [space.base_point_map()
                                                           for space in factors])
        self.base_point().rename('*')
        self._factors = factors

    def inclusion_map(self, i):
        """
        Return the inclusion map of the $i$-th factor.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: S2 = simplicial_sets.Sphere(2)
            sage: W = S1.wedge(S2, S1)
            sage: W.inclusion_map(1)
            Simplicial set morphism:
              From: S^2
              To:   Wedge: (S^1 v S^2 v S^1)
              Defn: [v_0, sigma_2] --> [*, sigma_2]
            sage: W.inclusion_map(0).domain()
            S^1
            sage: W.inclusion_map(2).domain()
            S^1
        """
        return self.structure_map(i)

    def projection_map(self, i):
        """
        Return the projection map onto the $i$-th factor.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: S2 = simplicial_sets.Sphere(2)
            sage: W = S1.wedge(S2, S1)
            sage: W.projection_map(1)
            Simplicial set morphism:
              From: Wedge: (S^1 v S^2 v S^1)
              To:   Quotient: (Wedge: (S^1 v S^2 v S^1)/Simplicial set with 3 non-degenerate simplices)
              Defn: [*, sigma_1, sigma_1, sigma_2] --> [*, s_0 *, s_0 *, sigma_2]
            sage: W.projection_map(1).image().homology(1)
            0
            sage: W.projection_map(1).image().homology(2)
            Z
        """
        m = len(self._factors)
        simplices = ([self.inclusion_map(j).image().nondegenerate_simplices()
                      for j in range(i)]
                     + [self.inclusion_map(j).image().nondegenerate_simplices()
                        for j in range(i+1,m)])
        return self.quotient(list(itertools.chain(*simplices))).quotient_map()


class DisjointUnionOfSimplicialSets(PushoutOfSimplicialSets, Factors):
    @staticmethod
    def __classcall__(cls, factors=None):
        """
        TESTS::

            sage: from sage.topology.simplicial_set_constructions import DisjointUnionOfSimplicialSets
            sage: from sage.topology.simplicial_set_examples import Empty
            sage: S2 = simplicial_sets.Sphere(2)
            sage: DisjointUnionOfSimplicialSets([S2, S2]) == DisjointUnionOfSimplicialSets((S2, S2))
            True
            sage: DisjointUnionOfSimplicialSets([S2, Empty(), S2, Empty()]) == DisjointUnionOfSimplicialSets((S2, S2))
            True
        """
        if factors:
            # Discard any empty factors.
            factors = [F for F in factors if F != Empty()]
        if factors:
            return super(DisjointUnionOfSimplicialSets, cls).__classcall__(cls, factors=tuple(factors))
        return super(DisjointUnionOfSimplicialSets, cls).__classcall__(cls)

    def __init__(self, factors=None):
        r"""
        Return the disjoint union of simplicial sets.

        INPUT:

        - ``factors`` -- a list or tuple of simplicial sets

        Discard any factors which are empty and return the disjoint
        union of the remaining simplicial sets in ``factors``.  The
        disjoint union comes equipped with a map from each factor, as
        long as all of the factors are finite.

        EXAMPLES::

            sage: CP2 = simplicial_sets.ComplexProjectiveSpace(2)
            sage: K = simplicial_sets.KleinBottle()
            sage: W = CP2.disjoint_union(K)
            sage: W.homology()
            {0: Z, 1: Z x C2, 2: Z, 3: 0, 4: Z}

            sage: W.inclusion_map(1)
            Simplicial set morphism:
              From: Klein bottle
              To:   Disjoint union: (CP^2 u Klein bottle)
              Defn: [Delta_{0,0}, Delta_{1,0}, Delta_{1,1}, Delta_{1,2}, Delta_{2,0}, Delta_{2,1}] --> [Delta_{0,0}, Delta_{1,0}, Delta_{1,1}, Delta_{1,2}, Delta_{2,0}, Delta_{2,1}]
        """
        PushoutOfSimplicialSets.__init__(self, [space._map_from_empty_set()
                                                for space in factors])
        self._factors = factors
        self._n_skeleton = (-1, Empty())

    def n_skeleton(self, n):
        """
        Return the `n`-skeleton of this simplicial set.

        That is, the simplicial set generated by all nondegenerate
        simplices of dimension at most `n`.

        INPUT:

        - ``n`` -- the dimension

        The `n`-skeleton of the disjoint union is computed as the
        disjoint union of the `n`-skeleta of the component simplicial
        sets.

        EXAMPLES::

            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: T = simplicial_sets.Torus()
            sage: X = B.disjoint_union(T)
            sage: X.n_skeleton(3).homology()
            {0: Z, 1: Z x Z x C2, 2: Z, 3: Z}
        """
        if self.is_finite():
            return DisjointUnionOfSimplicialSets_finite(tuple([X.n_skeleton(n)
                                                               for X in self._factors]))
        start, skel = self._n_skeleton
        if start == n:
            return skel
        elif start > n:
            return skel.n_skeleton(n)
        ans = DisjointUnionOfSimplicialSets_finite(tuple([X.n_skeleton(n)
                                                          for X in self._factors]))
        self._n_skeleton = (n, ans)
        return ans

    summands = Factors.factors
    summand = Factors.factor

    def _repr_(self):
        """
        Print representation

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: RP3 = simplicial_sets.RealProjectiveSpace(3)
            sage: T.disjoint_union(T, RP3)
            Disjoint union: (Torus u Torus u RP^3)
        """
        s = 'Disjoint union: ('
        s += ' u '.join([str(X) for X in self._factors])
        s += ')'
        return s

    def _latex_(self):
        r"""
        LaTeX representation

        EXAMPLES::

            sage: RP4 = simplicial_sets.RealProjectiveSpace(4)
            sage: S1 = simplicial_sets.Sphere(1)
            sage: latex(S1.disjoint_union(RP4, S1))
            S^{1} \amalg RP^{4} \amalg S^{1}
        """
        return ' \\amalg '.join([latex(X) for X in self._factors])


class DisjointUnionOfSimplicialSets_finite(DisjointUnionOfSimplicialSets,
                                           PushoutOfSimplicialSets_finite):
    """
    The disjoint union of finite simplicial sets.
    """
    def __init__(self, factors=None):
        r"""
        Return the disjoint union of finite simplicial sets.

        INPUT:

        - ``factors`` -- a tuple of simplicial sets

        Return the disjoint union of the simplicial sets in
        ``factors``.  The disjoint union comes equipped with a map
        from each factor. If there are no factors, this returns the
        empty simplicial set.

        EXAMPLES::

            sage: from sage.topology.simplicial_set_constructions import DisjointUnionOfSimplicialSets_finite
            sage: from sage.topology.simplicial_set_examples import Empty
            sage: S = simplicial_sets.Sphere(4)
            sage: DisjointUnionOfSimplicialSets_finite((S,S,S))
            Disjoint union: (S^4 u S^4 u S^4)
            sage: DisjointUnionOfSimplicialSets_finite([Empty(), Empty()]) == Empty()
            True
        """
        if not factors:
            PushoutOfSimplicialSets_finite.__init__(self)
        else:
            PushoutOfSimplicialSets_finite.__init__(self, [space._map_from_empty_set()
                                                           for space in factors])
        self._factors = factors

    def inclusion_map(self, i):
        """
        Return the inclusion map of the $i$-th factor.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: S2 = simplicial_sets.Sphere(2)
            sage: W = S1.disjoint_union(S2, S1)
            sage: W.inclusion_map(1)
            Simplicial set morphism:
              From: S^2
              To:   Disjoint union: (S^1 u S^2 u S^1)
              Defn: [v_0, sigma_2] --> [v_0, sigma_2]
            sage: W.inclusion_map(0).domain()
            S^1
            sage: W.inclusion_map(2).domain()
            S^1
        """
        return self.structure_map(i)


class ConeOfSimplicialSet(SimplicialSet_arbitrary, UniqueRepresentation):
    def __init__(self, base):
        r"""
        Return the unreduced cone on a finite simplicial set.

        INPUT:

        - ``base`` -- return the cone on this simplicial set.

        Add a point `*` (which will become the base point) and for
        each simplex `\sigma` in ``base``, add both `\sigma` and a
        simplex made up of `*` and `\sigma` (topologically, form the
        join of `*` and `\sigma`).

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(1, name='e')
            sage: X = SimplicialSet({e: (v, v)})
            sage: CX = X.cone() # indirect doctest
            sage: CX.nondegenerate_simplices()
            [*, v, (v,*), e, (e,*)]
            sage: CX.base_point()
            *
        """
        Cat = SimplicialSets().Pointed()
        if base.is_finite():
            Cat = Cat.Finite()
        Parent.__init__(self, category=Cat)
        star = AbstractSimplex(0, name='*')
        self._base = base
        self._basepoint = star
        self._n_skeleton = (-1, Empty())

    def n_skeleton(self, n):
        """
        Return the `n`-skeleton of this simplicial set.

        That is, the simplicial set generated by all nondegenerate
        simplices of dimension at most `n`.

        INPUT:

        - ``n`` -- the dimension

        In the case when the cone is infinite, the `n`-skeleton of the
        cone is computed as the `n`-skeleton of the cone of the
        `n`-skeleton.

        EXAMPLES::

            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: X = B.disjoint_union(B)
            sage: CX = B.cone()
            sage: CX.n_skeleton(3).homology()
            {0: 0, 1: 0, 2: 0, 3: Z}
        """
        if self.is_finite():
            return SimplicialSet_finite.n_skeleton(self, n)
        start, skel = self._n_skeleton
        if start == n:
            return skel
        elif start > n:
            return skel.n_skeleton(n)
        ans = ConeOfSimplicialSet_finite(self._base.n_skeleton(n)).n_skeleton(n)
        self._n_skeleton = (n, ans)
        self._basepoint = ans.base_point()
        return ans

    def _repr_(self):
        """
        Print representation

        EXAMPLES::

            sage: simplicial_sets.Simplex(3).cone()
            Cone of 3-simplex
        """
        return 'Cone of {}'.format(self._base)

    def _latex_(self):
        r"""
        LaTeX representation

        EXAMPLES::

            sage: latex(simplicial_sets.Simplex(3).cone())
            C \Delta^{3}
        """
        return 'C {}'.format(latex(self._base))


class ConeOfSimplicialSet_finite(ConeOfSimplicialSet, SimplicialSet_finite):
    def __init__(self, base):
        r"""
        Return the unreduced cone on a finite simplicial set.

        INPUT:

        - ``base`` -- return the cone on this simplicial set.

        Add a point `*` (which will become the base point) and for
        each simplex `\sigma` in ``base``, add both `\sigma` and a
        simplex made up of `*` and `\sigma` (topologically, form the
        join of `*` and `\sigma`).

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(1, name='e')
            sage: X = SimplicialSet({e: (v, v)})
            sage: CX = X.cone() # indirect doctest
            sage: CX.nondegenerate_simplices()
            [*, v, (v,*), e, (e,*)]
            sage: CX.base_point()
            *
        """
        star = AbstractSimplex(0, name='*')
        data = {}
        data[star] = None
        # Dictionary for translating old simplices to new: keys are
        # old simplices, corresponding value is the new simplex
        # (sigma, *).
        new_simplices = {'cone': star}
        for sigma in base.nondegenerate_simplices():
            new = AbstractSimplex(sigma.dimension()+1,
                                  name='({},*)'.format(sigma),
                                  latex_name='({},*)'.format(latex(sigma)))
            if sigma.dimension() == 0:
                data[sigma] = None
                data[new] = (star, sigma)
            else:
                sigma_faces = base.face_data()[sigma]
                data[sigma] = sigma_faces
                new_faces = [new_simplices[face.nondegenerate()].apply_degeneracies(*face.degeneracies())
                             for face in sigma_faces]
                data[new] = (new_faces + [sigma])
            new_simplices[sigma] = new
        SimplicialSet_finite.__init__(self, data, base_point=star)
        # self._base: original simplicial set.
        self._base = base
        # self._joins: dictionary, each key is a simplex sigma in
        # base, the corresponding value is the new simplex (sigma, *)
        # in the cone. Also, one other key is 'cone', and the value is
        # the cone vertex. This is used in the suspension class to
        # construct the suspension of a morphism. It could be used to
        # construct the cone of a morphism, also, although cones of
        # morphisms are not yet implemented.
        self._joins = new_simplices

    def base_as_subset(self):
        """
        If this is the cone `CX` on `X`, return `X` as a subsimplicial set.

        EXAMPLES::

            sage: X = simplicial_sets.RealProjectiveSpace(4).unset_base_point()
            sage: Y = X.cone()
            sage: Y.base_as_subset()
            Simplicial set with 5 non-degenerate simplices
            sage: Y.base_as_subset() == X
            True
        """
        X = self._base
        return self.subsimplicial_set(X.nondegenerate_simplices())

    def map_from_base(self):
        r"""
        If this is the cone `CX` on `X`, return the inclusion map from `X`.

        EXAMPLES::

            sage: X = simplicial_sets.Simplex(2).n_skeleton(1)
            sage: Y = X.cone()
            sage: Y.map_from_base()
            Simplicial set morphism:
              From: Simplicial set with 6 non-degenerate simplices
              To:   Cone of Simplicial set with 6 non-degenerate simplices
              Defn: [(0,), (1,), (2,), (0, 1), (0, 2), (1, 2)] --> [(0,), (1,), (2,), (0, 1), (0, 2), (1, 2)]
        """
        return self.base_as_subset().inclusion_map()


class ReducedConeOfSimplicialSet(QuotientOfSimplicialSet):
    def __init__(self, base):
        r"""
        Return the reduced cone on a simplicial set.

        INPUT:

        - ``base`` -- return the cone on this simplicial set.

        Start with the unreduced cone: take ``base`` and add a point
        `*` (which will become the base point) and for each simplex
        `\sigma` in ``base``, add both `\sigma` and a simplex made up
        of `*` and `\sigma` (topologically, form the join of `*` and
        `\sigma`).

        Now reduce: take the quotient by the 1-simplex connecting the
        old base point to the new one.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(1, name='e')
            sage: X = SimplicialSet({e: (v, v)})
            sage: X = X.set_base_point(v)
            sage: CX = X.cone()  # indirect doctest
            sage: CX.nondegenerate_simplices()
            [*, e, (e,*)]
        """
        C = ConeOfSimplicialSet(base)
        for t in C.n_cells(1):
            edge_faces = sorted([C.base_point(), base.base_point()])
            if sorted(C.faces(t)) == edge_faces:
                edge = t
                break
        inc = C.subsimplicial_set([edge]).inclusion_map()
        QuotientOfSimplicialSet.__init__(self, inc)
        self._base = base
        self._n_skeleton = (-1, Empty())

    def n_skeleton(self, n):
        """
        Return the `n`-skeleton of this simplicial set.

        That is, the simplicial set generated by all nondegenerate
        simplices of dimension at most `n`.

        INPUT:

        - ``n`` -- the dimension

        In the case when the cone is infinite, the `n`-skeleton of the
        cone is computed as the `n`-skeleton of the cone of the
        `n`-skeleton.

        EXAMPLES::

            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: B.cone().n_skeleton(3).homology()
            {0: 0, 1: 0, 2: 0, 3: Z}
        """
        if self.is_finite():
            return SimplicialSet_finite.n_skeleton(self, n)
        start, skel = self._n_skeleton
        if start == n:
            return skel
        elif start > n:
            return skel.n_skeleton(n)
        ans = ReducedConeOfSimplicialSet_finite(self._base.n_skeleton(n)).n_skeleton(n)
        self._n_skeleton = (n, ans)
        return ans

    def _repr_(self):
        """
        Print representation

        EXAMPLES::

            sage: X = simplicial_sets.Sphere(4)
            sage: X.cone()
            Reduced cone of S^4
        """
        return 'Reduced cone of {}'.format(self._base)

    def _latex_(self):
        r"""
        LaTeX representation

        EXAMPLES::

            sage: latex(simplicial_sets.Sphere(4).cone())
            C S^{4}
        """
        return 'C {}'.format(latex(self._base))


class ReducedConeOfSimplicialSet_finite(ReducedConeOfSimplicialSet,
                                        QuotientOfSimplicialSet_finite):
    def __init__(self, base):
        r"""
        Return the reduced cone on a simplicial set.

        INPUT:

        - ``base`` -- return the cone on this simplicial set.

        Start with the unreduced cone: take ``base`` and add a point
        `*` (which will become the base point) and for each simplex
        `\sigma` in ``base``, add both `\sigma` and a simplex made up
        of `*` and `\sigma` (topologically, form the join of `*` and
        `\sigma`).

        Now reduce: take the quotient by the 1-simplex connecting the
        old base point to the new one.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(1, name='e')
            sage: X = SimplicialSet({e: (v, v)})
            sage: X = X.set_base_point(v)
            sage: CX = X.cone()  # indirect doctest
            sage: CX.nondegenerate_simplices()
            [*, e, (e,*)]
        """
        C = ConeOfSimplicialSet_finite(base)
        edge_faces = sorted([C.base_point(), base.base_point()])
        for t in C.n_cells(1):
            if sorted(C.faces(t)) == edge_faces:
                edge = t
                break
        inc = C.subsimplicial_set([edge]).inclusion_map()
        QuotientOfSimplicialSet_finite.__init__(self, inc)
        self._base = base
        q = self.quotient_map()
        self._joins = {sigma:q(C._joins[sigma]) for sigma in C._joins}

    def map_from_base(self):
        r"""
        If this is the cone `\tilde{C}X` on `X`, return the map from `X`.

        The map is defined to be the composite `X \to CX \to
        \tilde{C}X`.  This is used by the
        :class:`SuspensionOfSimplicialSet_finite` class to construct
        the reduced suspension: take the quotient of the reduced cone
        by the image of `X` therein.

        EXAMPLES::

            sage: S3 = simplicial_sets.Sphere(3)
            sage: CS3 = S3.cone()
            sage: CS3.map_from_base()
            Simplicial set morphism:
              From: S^3
              To:   Reduced cone of S^3
              Defn: [v_0, sigma_3] --> [*, sigma_3]
        """
        quotient_map = self.quotient_map()
        unreduced = quotient_map.domain()
        temp_map = unreduced.map_from_base()
        X = self._base
        incl = X.Hom(unreduced)(temp_map._dictionary)
        return quotient_map * incl


class SuspensionOfSimplicialSet(SimplicialSet_arbitrary, UniqueRepresentation):
    def __init__(self, base):
        r"""
        Return the (reduced) suspension of a simplicial set.

        INPUT:

        - ``base`` -- return the suspension of this simplicial set.

        If this simplicial set ``X=base`` is not pointed, or if it is
        itself an unreduced suspension, return the unreduced
        suspension: the quotient `CX/X`, where `CX` is the (ordinary,
        unreduced) cone on `X`. If `X` is pointed, then use the
        reduced cone instead, and so return the reduced suspension.

        We use `S` to denote unreduced suspension, `\Sigma` for
        reduced suspension.

        EXAMPLES::

            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: B.suspension()
            Sigma(Classifying space of Multiplicative Abelian group isomorphic to C2)
            sage: B.suspension().n_skeleton(3).homology()
            {0: 0, 1: 0, 2: C2, 3: 0}

        If ``X`` is finite, the suspension comes with a quotient map
        from the cone::

            sage: S3 = simplicial_sets.Sphere(3)
            sage: S4 = S3.suspension()
            sage: S4.quotient_map()
            Simplicial set morphism:
              From: Reduced cone of S^3
              To:   Sigma(S^3)
              Defn: [*, sigma_3, (sigma_3,*)] --> [*, s_2 s_1 s_0 *, (sigma_3,*)]

        TESTS::

            sage: S3.suspension() == S3.suspension()
            True
            sage: S3.suspension() == simplicial_sets.Sphere(3).suspension()
            False
            sage: B.suspension() == B.suspension()
            True
        """
        Cat = SimplicialSets()
        if base.is_finite():
            Cat = Cat.Finite()
        reduced = (base.is_pointed()
                   and (not hasattr(base, '_reduced')
                        or (hasattr(base, '_reduced') and base._reduced)))
        if reduced:
            Cat = Cat.Pointed()
        Parent.__init__(self, category=Cat)
        self._reduced = reduced
        self._base = base
        self._n_skeleton = (-1, Empty())

    def n_skeleton(self, n):
        """
        Return the `n`-skeleton of this simplicial set.

        That is, the simplicial set generated by all nondegenerate
        simplices of dimension at most `n`.

        INPUT:

        - ``n`` -- the dimension

        In the case when the suspension is infinite, the `n`-skeleton
        of the suspension is computed as the `n`-skeleton of the
        suspension of the `n`-skeleton.

        EXAMPLES::

            sage: B = simplicial_sets.ClassifyingSpace(groups.misc.MultiplicativeAbelian([2]))
            sage: SigmaB = B.suspension()
            sage: SigmaB.n_skeleton(4).homology(base_ring=GF(2))
            {0: Vector space of dimension 0 over Finite Field of size 2,
             1: Vector space of dimension 0 over Finite Field of size 2,
             2: Vector space of dimension 1 over Finite Field of size 2,
             3: Vector space of dimension 1 over Finite Field of size 2,
             4: Vector space of dimension 1 over Finite Field of size 2}
        """
        if self.is_finite():
            return SimplicialSet_finite.n_skeleton(self, n)
        start, skel = self._n_skeleton
        if start == n:
            return skel
        elif start > n:
            return skel.n_skeleton(n)
        ans = SuspensionOfSimplicialSet_finite(self._base.n_skeleton(n)).n_skeleton(n)
        self._n_skeleton = (n, ans)
        return ans

    def __repr_or_latex__(self, output_type=None):
        r"""
        Print representation, for either :meth:`_repr_` or :meth:`_latex_`.

        INPUT:

        - ``output_type`` -- either ``"latex"`` for LaTeX output or
          anything else for ``str`` output.

        We use `S` to denote unreduced suspension, `\Sigma` for
        reduced suspension.

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: K = T.suspension(10)
            sage: K.__repr_or_latex__()
            'Sigma^10(Torus)'
            sage: K.__repr_or_latex__('latex')
            '\\Sigma^{10}(S^{1} \\times S^{1})'
        """
        latex_output = (output_type == 'latex')
        base = self._base
        if self._reduced:
            # Reduced suspension.
            if latex_output:
                symbol = '\\Sigma'
            else:
                symbol = 'Sigma'
        else:
            # Unreduced suspension.
            symbol = 'S'
        idx = 1
        while isinstance(base, SuspensionOfSimplicialSet):
            idx += 1
            base = base._base
        if latex_output:
            base = latex(base)
            exp = '^{{{}}}'
        else:
            exp = '^{}'
        if idx > 1:
            return ('{}' + exp + '({})').format(symbol, idx, base)
        else:
            return ('{}({})').format(symbol, base)

    def _repr_(self):
        r"""
        Print representation

        We use `S` to denote unreduced suspension, `\Sigma` for
        reduced suspension.

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: S2.suspension(3)
            Sigma^3(S^2)
            sage: K = simplicial_sets.Simplex(2)
            sage: K.suspension(3)
            S^3(2-simplex)
            sage: K.suspension()
            S(2-simplex)
        """
        return self.__repr_or_latex__()

    def _latex_(self):
        r"""
        LaTeX representation

        We use `S` to denote unreduced suspension, `\Sigma` for
        reduced suspension.

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: latex(S2.suspension(3))
            \Sigma^{3}(S^{2})
            sage: K = simplicial_sets.Simplex(2)
            sage: latex(K.suspension(3))
            S^{3}(\Delta^{2})
            sage: latex(K.suspension())
            S(\Delta^{2})
        """
        return self.__repr_or_latex__('latex')


class SuspensionOfSimplicialSet_finite(SuspensionOfSimplicialSet,
                                       QuotientOfSimplicialSet_finite):
    """
    The (reduced) suspension of a finite simplicial set.

    See :class:`SuspensionOfSimplicialSet` for more information.
    """
    def __init__(self, base):
        r"""
        INPUT:

        - ``base`` -- return the suspension of this finite simplicial set.

        See :class:`SuspensionOfSimplicialSet` for more information.

        EXAMPLES::

            sage: X = simplicial_sets.Sphere(3)
            sage: X.suspension(2)
            Sigma^2(S^3)
            sage: Y = X.unset_base_point()
            sage: Y.suspension(2)
            S^2(Simplicial set with 2 non-degenerate simplices)
        """
        self._base = base
        reduced = (base.is_pointed()
                   and (not hasattr(base, '_reduced')
                        or (hasattr(base, '_reduced') and base._reduced)))
        if reduced:
            C = ReducedConeOfSimplicialSet_finite(base)
            subcomplex = C.map_from_base().image()
        else:
            C = ConeOfSimplicialSet_finite(base)
            subcomplex = C.base_as_subset()
        QuotientOfSimplicialSet_finite.__init__(self, subcomplex.inclusion_map())
        self._reduced = reduced
        # self._suspensions: dictionary, each key is a simplex sigma
        # in base, the corresponding value is the new simplex (sigma, *)
        # in S(base). Another key is 'cone', and its value is the cone
        # vertex in C(base). This is used to construct the suspension of a
        # morphism.
        q = self.quotient_map()
        self._suspensions = {sigma: q(C._joins[sigma]) for sigma in C._joins}
