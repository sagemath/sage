# sage_setup: distribution = sagemath-categories
"""
Simplicial Sets
"""
# ****************************************************************************
#  Copyright (C) 2015 John H. Palmieri <palmieri at math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_singleton import Category_singleton
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.homsets import HomsetsCategory
from sage.categories.sets_cat import Sets
from sage.functions.generalized import sign
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ

lazy_import('sage.matrix.constructor', 'matrix')
lazy_import('sage.matrix.special', 'identity_matrix')


class SimplicialSets(Category_singleton):
    r"""
    The category of simplicial sets.

    A simplicial set `X` is a collection of sets `X_i`, indexed by
    the nonnegative integers, together with maps

    .. math::

        d_i: X_n \to X_{n-1}, \quad 0 \leq i \leq n \quad \text{(face maps)} \\
        s_j: X_n \to X_{n+1}, \quad 0 \leq j \leq n \quad \text{(degeneracy maps)}

    satisfying the *simplicial identities*:

    .. math::

        d_i d_j &= d_{j-1} d_i \quad \text{if } i<j \\
        d_i s_j &= s_{j-1} d_i \quad \text{if } i<j \\
        d_j s_j &= 1 = d_{j+1} s_j \\
        d_i s_j &= s_{j} d_{i-1} \quad \text{if } i>j+1 \\
        s_i s_j &= s_{j+1} s_{i} \quad \text{if } i \leq j

    Morphisms are sequences of maps `f_i : X_i \to Y_i` which commute
    with the face and degeneracy maps.

    EXAMPLES::

        sage: from sage.categories.simplicial_sets import SimplicialSets
        sage: C = SimplicialSets(); C
        Category of simplicial sets

    TESTS::

        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.simplicial_sets import SimplicialSets
            sage: SimplicialSets().super_categories()
            [Category of sets]
        """
        return [Sets()]

    class ParentMethods:
        def is_finite(self):
            """
            Return ``True`` if this simplicial set is finite, i.e., has a
            finite number of nondegenerate simplices.

            EXAMPLES::

                sage: simplicial_sets.Torus().is_finite()                               # needs sage.graphs
                True
                sage: C5 = groups.misc.MultiplicativeAbelian([5])                       # needs sage.graphs sage.groups
                sage: simplicial_sets.ClassifyingSpace(C5).is_finite()                  # needs sage.graphs sage.groups
                False
            """
            return SimplicialSets.Finite() in self.categories()

        def is_pointed(self):
            """
            Return ``True`` if this simplicial set is pointed, i.e., has a
            base point.

            EXAMPLES::

                sage: # needs sage.graphs
                sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
                sage: v = AbstractSimplex(0)
                sage: w = AbstractSimplex(0)
                sage: e = AbstractSimplex(1)
                sage: X = SimplicialSet({e: (v, w)})
                sage: Y = SimplicialSet({e: (v, w)}, base_point=w)
                sage: X.is_pointed()
                False
                sage: Y.is_pointed()
                True
            """
            return SimplicialSets.Pointed() in self.categories()

        def set_base_point(self, point):
            """
            Return a copy of this simplicial set in which the base point is
            set to ``point``.

            INPUT:

            - ``point`` -- a 0-simplex in this simplicial set

            EXAMPLES::

                sage: # needs sage.graphs
                sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
                sage: v = AbstractSimplex(0, name='v_0')
                sage: w = AbstractSimplex(0, name='w_0')
                sage: e = AbstractSimplex(1)
                sage: X = SimplicialSet({e: (v, w)})
                sage: Y = SimplicialSet({e: (v, w)}, base_point=w)
                sage: Y.base_point()
                w_0
                sage: X_star = X.set_base_point(w)
                sage: X_star.base_point()
                w_0
                sage: Y_star = Y.set_base_point(v)
                sage: Y_star.base_point()
                v_0

            TESTS::

                sage: X.set_base_point(e)                                               # needs sage.graphs
                Traceback (most recent call last):
                ...
                ValueError: the "point" is not a zero-simplex
                sage: pt = AbstractSimplex(0)                                           # needs sage.graphs
                sage: X.set_base_point(pt)                                              # needs sage.graphs
                Traceback (most recent call last):
                ...
                ValueError: the point is not a simplex in this simplicial set
            """
            from sage.topology.simplicial_set import SimplicialSet
            if point.dimension() != 0:
                raise ValueError('the "point" is not a zero-simplex')
            if point not in self._simplices:
                raise ValueError('the point is not a simplex in this '
                                 'simplicial set')
            return SimplicialSet(self.face_data(), base_point=point)

    class Homsets(HomsetsCategory):
        class Endset(CategoryWithAxiom):
            class ParentMethods:
                def one(self):
                    r"""
                    Return the identity morphism in `\operatorname{Hom}(S, S)`.

                    EXAMPLES::

                        sage: T = simplicial_sets.Torus()                               # needs sage.graphs
                        sage: Hom(T, T).identity()                                      # needs sage.graphs
                        Simplicial set endomorphism of Torus
                          Defn: Identity map
                    """
                    from sage.topology.simplicial_set_morphism import SimplicialSetMorphism
                    return SimplicialSetMorphism(domain=self.domain(),
                                                 codomain=self.codomain(),
                                                 identity=True)

    class Finite(CategoryWithAxiom):
        """
        Category of finite simplicial sets.

        The objects are simplicial sets with finitely many
        non-degenerate simplices.
        """
        pass

    class SubcategoryMethods:
        def Pointed(self):
            """
            A simplicial set is *pointed* if it has a distinguished base
            point.

            EXAMPLES::

                sage: from sage.categories.simplicial_sets import SimplicialSets
                sage: SimplicialSets().Pointed().Finite()
                Category of finite pointed simplicial sets
                sage: SimplicialSets().Finite().Pointed()
                Category of finite pointed simplicial sets
            """
            return self._with_axiom("Pointed")

    class Pointed(CategoryWithAxiom):
        class ParentMethods:
            def base_point(self):
                """
                Return this simplicial set's base point.

                EXAMPLES::

                    sage: # needs sage.graphs
                    sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
                    sage: v = AbstractSimplex(0, name='*')
                    sage: e = AbstractSimplex(1)
                    sage: S1 = SimplicialSet({e: (v, v)}, base_point=v)
                    sage: S1.is_pointed()
                    True
                    sage: S1.base_point()
                    *
                """
                return self._basepoint

            def base_point_map(self, domain=None):
                """
                Return a map from a one-point space to this one, with image the
                base point.

                This raises an error if this simplicial set does not have a
                base point.

                INPUT:

                - ``domain`` -- (default: ``None``) use this to specify a
                  particular one-point space as the domain. The default
                  behavior is to use the
                  :func:`sage.topology.simplicial_set.Point` function to use a
                  standard one-point space.

                EXAMPLES::

                    sage: # needs sage.graphs
                    sage: T = simplicial_sets.Torus()
                    sage: f = T.base_point_map(); f
                    Simplicial set morphism:
                      From: Point
                      To:   Torus
                      Defn: Constant map at (v_0, v_0)
                    sage: S3 = simplicial_sets.Sphere(3)
                    sage: g = S3.base_point_map()
                    sage: f.domain() == g.domain()
                    True
                    sage: RP3 = simplicial_sets.RealProjectiveSpace(3)                  # needs sage.groups
                    sage: temp = simplicial_sets.Simplex(0)
                    sage: pt = temp.set_base_point(temp.n_cells(0)[0])
                    sage: h = RP3.base_point_map(domain=pt)                             # needs sage.groups
                    sage: f.domain() == h.domain()                                      # needs sage.groups
                    False

                    sage: C5 = groups.misc.MultiplicativeAbelian([5])                   # needs sage.graphs sage.groups
                    sage: BC5 = simplicial_sets.ClassifyingSpace(C5)                    # needs sage.graphs sage.groups
                    sage: BC5.base_point_map()                                          # needs sage.graphs sage.groups
                    Simplicial set morphism:
                      From: Point
                      To:   Classifying space of Multiplicative Abelian group isomorphic to C5
                      Defn: Constant map at 1
                """
                from sage.topology.simplicial_set_examples import Point
                if domain is None:
                    domain = Point()
                else:
                    if len(domain._simplices) > 1:
                        raise ValueError('domain has more than one nondegenerate simplex')
                target = self.base_point()
                return domain.Hom(self).constant_map(point=target)

            def fundamental_group(self, simplify=True):
                r"""
                Return the fundamental group of this pointed simplicial set.

                INPUT:

                - ``simplify`` -- boolean (default: ``True``); if
                  ``False``, then return a presentation of the group
                  in terms of generators and relations. If ``True``,
                  the default, simplify as much as GAP is able to.

                Algorithm: we compute the edge-path group -- see
                Section 19 of [Kan1958]_ and
                :wikipedia:`Fundamental_group`. Choose a spanning tree
                for the connected component of the 1-skeleton
                containing the base point, and then the group's
                generators are given by the non-degenerate
                edges. There are two types of relations: `e=1` if `e`
                is in the spanning tree, and for every 2-simplex, if
                its faces are `e_0`, `e_1`, and `e_2`, then we impose
                the relation `e_0 e_1^{-1} e_2 = 1`, where we first
                set `e_i=1` if `e_i` is degenerate.

                EXAMPLES::

                    sage: S1 = simplicial_sets.Sphere(1)                                # needs sage.graphs
                    sage: eight = S1.wedge(S1)                                          # needs sage.graphs
                    sage: eight.fundamental_group()  # free group on 2 generators       # needs sage.graphs sage.groups
                    Finitely presented group < e0, e1 |  >

                The fundamental group of a disjoint union of course depends on
                the choice of base point::

                    sage: T = simplicial_sets.Torus()                                   # needs sage.graphs
                    sage: K = simplicial_sets.KleinBottle()                             # needs sage.graphs
                    sage: X = T.disjoint_union(K)                                       # needs sage.graphs

                    sage: # needs sage.graphs
                    sage: X_0 = X.set_base_point(X.n_cells(0)[0])
                    sage: X_0.fundamental_group().is_abelian()                          # needs sage.groups
                    True
                    sage: X_1 = X.set_base_point(X.n_cells(0)[1])
                    sage: X_1.fundamental_group().is_abelian()                          # needs sage.groups
                    False

                    sage: RP3 = simplicial_sets.RealProjectiveSpace(3)                  # needs sage.graphs sage.groups
                    sage: RP3.fundamental_group()                                       # needs sage.graphs sage.groups
                    Finitely presented group < e | e^2 >

                Compute the fundamental group of some classifying spaces::

                    sage: C5 = groups.misc.MultiplicativeAbelian([5])                   # needs sage.graphs sage.groups
                    sage: BC5 = C5.nerve()                                              # needs sage.graphs sage.groups
                    sage: BC5.fundamental_group()                                       # needs sage.graphs sage.groups
                    Finitely presented group < e0 | e0^5 >

                    sage: # needs sage.graphs sage.groups
                    sage: Sigma3 = groups.permutation.Symmetric(3)
                    sage: BSigma3 = Sigma3.nerve()
                    sage: pi = BSigma3.fundamental_group(); pi
                    Finitely presented group < e1, e2 | e2^2, e1^3, (e2*e1)^2 >
                    sage: pi.order()
                    6
                    sage: pi.is_abelian()
                    False

                The sphere has a trivial fundamental group::

                    sage: S2 = simplicial_sets.Sphere(2)                                # needs sage.graphs
                    sage: S2.fundamental_group()                                        # needs sage.graphs sage.groups
                    Finitely presented group <  |  >
                """
                # Import this here to prevent importing libgap upon startup.
                from sage.groups.free_group import FreeGroup
                if not self.n_cells(1):
                    return FreeGroup([]).quotient([])
                FG = self._universal_cover_dict()[0]
                if simplify:
                    return FG.simplified()
                else:
                    return FG

            def _universal_cover_dict(self):
                r"""
                Return the fundamental group and dictionary sending each edge to
                the corresponding group element

                TESTS::

                    sage: RP2 = simplicial_sets.RealProjectiveSpace(2)                  # needs sage.graphs sage.groups
                    sage: RP2._universal_cover_dict()                                   # needs sage.graphs sage.groups
                    (Finitely presented group < e | e^2 >, {f: e})
                    sage: RP2.nondegenerate_simplices()                                 # needs sage.graphs sage.groups
                    [1, f, f * f]
                """
                from sage.groups.free_group import FreeGroup
                graph = self.graph()
                if not self.is_connected():
                    graph = graph.subgraph(self.base_point())
                edges = [e[2] for e in graph.edges(sort=False)]
                spanning_tree = [e[2] for e in graph.min_spanning_tree()]
                gens = [e for e in edges if e not in spanning_tree]
                gens_dict = dict(zip(gens, range(len(gens))))
                FG = FreeGroup(len(gens), 'e')
                rels = []

                for f in self.n_cells(2):
                    z = {}
                    for i, sigma in enumerate(self.faces(f)):
                        if sigma in spanning_tree:
                            z[i] = FG.one()
                        elif sigma.is_degenerate():
                            z[i] = FG.one()
                        elif sigma in edges:
                            z[i] = FG.gen(gens_dict[sigma])
                        else:
                            # sigma is not in the correct connected component.
                            z[i] = FG.one()
                    rels.append(z[0]*z[1].inverse()*z[2])
                G = FG.quotient(rels)
                char = {g: G.gen(i) for i, g in enumerate(gens)}
                for e in edges:
                    if e not in gens:
                        char[e] = G.one()
                return (G, char)

            def universal_cover_map(self):
                r"""
                Return the universal covering map of the simplicial set.

                It requires the fundamental group to be finite.

                EXAMPLES::

                    sage: RP2 = simplicial_sets.RealProjectiveSpace(2)                  # needs sage.graphs sage.groups
                    sage: phi = RP2.universal_cover_map(); phi                          # needs sage.graphs sage.groups gap_package_polenta
                    Simplicial set morphism:
                      From: Simplicial set with 6 non-degenerate simplices
                      To:   RP^2
                      Defn: [(1, 1), (1, e), (f, 1), (f, e), (f * f, 1), (f * f, e)]
                            --> [1, 1, f, f, f * f, f * f]
                    sage: phi.domain().face_data()                                      # needs sage.graphs sage.groups gap_package_polenta
                        {(1, 1): None,
                         (1, e): None,
                         (f, 1): ((1, e), (1, 1)),
                         (f, e): ((1, 1), (1, e)),
                         (f * f, 1): ((f, e), s_0 (1, 1), (f, 1)),
                         (f * f, e): ((f, 1), s_0 (1, e), (f, e))}
                """
                edges = self.n_cells(1)
                if not edges:
                    return self.identity()
                G, char = self._universal_cover_dict()
                return self.covering_map(char)

            def covering_map(self, character):
                r"""
                Return the covering map associated to a character.

                The character is represented by a dictionary that assigns an
                element of a finite group to each nondegenerate 1-dimensional
                cell. It should correspond to an epimorphism from the fundamental
                group.

                INPUT:

                - ``character`` -- dictionary

                EXAMPLES::

                    sage: # needs sage.graphs sage.groups
                    sage: S1 = simplicial_sets.Sphere(1)
                    sage: S1_ = simplicial_sets.Sphere(1)
                    sage: S1_.n_cells(1)[0].rename("sigma_1'")
                    sage: W = S1.wedge(S1_)
                    sage: G = CyclicPermutationGroup(3)
                    sage: a, b = W.n_cells(1)
                    sage: C = W.covering_map({a : G.gen(0), b : G.one()}); C
                    Simplicial set morphism:
                      From: Simplicial set with 9 non-degenerate simplices
                      To:   Wedge: (S^1 v S^1)
                      Defn: [(*, ()), (*, (1,2,3)), (*, (1,3,2)), (sigma_1', ()),
                             (sigma_1', (1,2,3)), (sigma_1', (1,3,2)), (sigma_1, ()),
                             (sigma_1, (1,2,3)), (sigma_1, (1,3,2))]
                            --> [*, *, *, sigma_1', sigma_1', sigma_1', sigma_1, sigma_1, sigma_1]
                    sage: C.domain()
                    Simplicial set with 9 non-degenerate simplices
                    sage: C.domain().face_data()
                    {(*, ()): None,
                     (*, (1,2,3)): None,
                     (*, (1,3,2)): None,
                     (sigma_1', ()): ((*, ()), (*, ())),
                     (sigma_1', (1,2,3)): ((*, (1,2,3)), (*, (1,2,3))),
                     (sigma_1', (1,3,2)): ((*, (1,3,2)), (*, (1,3,2))),
                     (sigma_1, ()): ((*, (1,2,3)), (*, ())),
                     (sigma_1, (1,2,3)): ((*, (1,3,2)), (*, (1,2,3))),
                     (sigma_1, (1,3,2)): ((*, ()), (*, (1,3,2)))}
                """
                from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
                from sage.topology.simplicial_set_morphism import SimplicialSetMorphism
                char = dict(character.items())
                G = list(char.values())[0].parent()
                if not G.is_finite():
                    raise NotImplementedError("can only compute universal covers of spaces with finite fundamental group")
                cells_dict = {}
                faces_dict = {}

                for s in self.n_cells(0):
                    for g in G:
                        cell = AbstractSimplex(0, name="({}, {})".format(s, g))
                        cells_dict[(s, g)] = cell
                        char[s] = G.one()

                for d in range(1, self.dimension() + 1):
                    for s in self.all_n_simplices(d):
                        if s not in char.keys():
                            if d == 1 and s.is_degenerate():
                                char[s] = G.one()
                            elif s.is_degenerate():
                                if 0 in s.degeneracies():
                                    char[s] = G.one()
                                else:
                                    char[s] = char[s.nondegenerate()]
                            else:
                                char[s] = char[self.face(s, d)]
                        if s.is_nondegenerate():
                            for g in G:
                                cell = AbstractSimplex(d, name="({}, {})".format(s, g))
                                cells_dict[(s, g)] = cell
                                faces = self.faces(s)
                                f0 = faces[0]
                                for h in G:
                                    if h == g*char[s]:
                                        lifted = h
                                        break
                                grelems = [cells_dict[(f0.nondegenerate(), lifted)].apply_degeneracies(*f0.degeneracies())]
                                grelems.extend(cells_dict[(f.nondegenerate(), g)].apply_degeneracies(*f.degeneracies()) for f in faces[1:])
                                faces_dict[cell] = grelems
                cover = SimplicialSet(faces_dict, base_point=cells_dict[(self.base_point(), G.one())])
                cover_map_data = {c: s[0] for (s, c) in cells_dict.items()}
                return SimplicialSetMorphism(data=cover_map_data, domain=cover, codomain=self)

            def cover(self, character):
                r"""
                Return the cover of the simplicial set associated to a character
                of the fundamental group.

                The character is represented by a dictionary, that assigns an
                element of a finite group to each nondegenerate 1-dimensional
                cell. It should correspond to an epimorphism from the fundamental
                group.

                INPUT:

                - ``character`` -- dictionary

                EXAMPLES::

                    sage: # needs sage.graphs sage.groups
                    sage: S1 = simplicial_sets.Sphere(1)
                    sage: S1_ = simplicial_sets.Sphere(1)
                    sage: S1_.n_cells(1)[0].rename("sigma_1'")
                    sage: W = S1.wedge(S1_)
                    sage: G = CyclicPermutationGroup(3)
                    sage: (a, b) = W.n_cells(1)
                    sage: C = W.cover({a : G.gen(0), b : G.gen(0)^2})
                    sage: C.face_data()
                    {(*, ()): None,
                     (*, (1,2,3)): None,
                     (*, (1,3,2)): None,
                     (sigma_1', ()): ((*, (1,3,2)), (*, ())),
                     (sigma_1', (1,2,3)): ((*, ()), (*, (1,2,3))),
                     (sigma_1', (1,3,2)): ((*, (1,2,3)), (*, (1,3,2))),
                     (sigma_1, ()): ((*, (1,2,3)), (*, ())),
                     (sigma_1, (1,2,3)): ((*, (1,3,2)), (*, (1,2,3))),
                     (sigma_1, (1,3,2)): ((*, ()), (*, (1,3,2)))}
                    sage: C.homology(1)                                                 # needs sage.modules
                    Z x Z x Z x Z
                    sage: C.fundamental_group()
                    Finitely presented group < e0, e1, e2, e3 |  >
                """
                return self.covering_map(character).domain()

            def universal_cover(self):
                r"""
                Return the universal cover of the simplicial set.
                The fundamental group must be finite in order to ensure that the
                universal cover is a simplicial set of finite type.

                EXAMPLES::

                    sage: # needs sage.graphs sage.groups
                    sage: RP3 = simplicial_sets.RealProjectiveSpace(3)
                    sage: C = RP3.universal_cover(); C
                    Simplicial set with 8 non-degenerate simplices
                    sage: C.face_data()  # needs gap_package_polenta
                    {(1, 1): None,
                     (1, e): None,
                     (f, 1): ((1, e), (1, 1)),
                     (f, e): ((1, 1), (1, e)),
                     (f * f, 1): ((f, e), s_0 (1, 1), (f, 1)),
                     (f * f, e): ((f, 1), s_0 (1, e), (f, e)),
                     (f * f * f, 1): ((f * f, e), s_0 (f, 1), s_1 (f, 1), (f * f, 1)),
                     (f * f * f, e): ((f * f, 1), s_0 (f, e), s_1 (f, e), (f * f, e))}
                    sage: C.fundamental_group()
                    Finitely presented group <  |  >

                TESTS::

                    sage: RP2 = simplicial_sets.RealProjectiveSpace(2)
                    sage: S3 = simplicial_sets.Sphere(3)
                    sage: X = S3.wedge(RP2)
                    sage: XU = X.universal_cover()
                    sage: [XU.homology(i) for i in range(5)]
                    [0, 0, Z, Z x Z, 0]

                """
                return self.universal_cover_map().domain()

            def _canonical_twisting_operator(self):
                r"""
                The canonical twisting operator corresponds to the abelianization of the fundamental
                group. It assigns each edge to the corresponding element in the algebra of the
                abelianization of the fundamental group, which is a quotient of a Laurent polynomial
                ring.

                EXAMPLES::

                    sage: # needs sage.graphs
                    sage: X = simplicial_sets.Torus()
                    sage: d = X._canonical_twisting_operator()
                    sage: d
                    {(s_0 v_0, sigma_1): f2, (sigma_1, s_0 v_0): f1*f2^-1, (sigma_1, sigma_1): f1}
                    sage: list(d.values())[0].parent()
                    Multivariate Laurent Polynomial Ring in f1, f2 over Integer Ring
                    sage: Y = simplicial_sets.RealProjectiveSpace(2)
                    sage: d2 = Y._canonical_twisting_operator()
                    sage: d2
                    {f: F1bar}
                    sage: list(d2.values())[0].parent()
                    Quotient of Univariate Laurent Polynomial Ring in F1 over Integer Ring by the ideal (-1 + F1^2)
                """
                G, d = self._universal_cover_dict()
                abelG, R, I, images = G.abelianization_to_algebra(ZZ)
                QRP = R.quotient_ring(I)
                res = {}
                for s, el in d.items():
                    res[s] = QRP(prod(images[abs(a)-1]**sign(a) for a in el.Tietze()))
                return res

            def twisted_chain_complex(self, twisting_operator=None, dimensions=None, augmented=False,
                                      cochain=False, verbose=False, subcomplex=None,
                                      check=False):
                r"""
                Return the normalized chain complex twisted by some operator.

                A twisting operator is a map from the set of simplices to some algebra.
                The differentials are then twisted by this operator.

                INPUT:

                - ``twisting_operator`` -- dictionary, associating the twist of each
                  simplex. If it is not given, the canonical one (associated to the
                  laurent polynomial ring abelianization of the fundamental group, ignoring
                  torsion) is used.

                - ``dimensions`` -- if ``None``, compute the chain complex in all
                  dimensions.  If a list or tuple of integers, compute the
                  chain complex in those dimensions, setting the chain groups
                  in all other dimensions to zero.

                - ``augmented`` -- boolean (default: ``False``); if ``True``,
                  return the augmented chain complex (that is, include a class
                  in dimension `-1` corresponding to the empty cell).

                - ``cochain`` -- boolean (default: ``False``); if ``True``,
                  return the cochain complex (that is, the dual of the chain
                  complex).

                - ``verbose`` -- boolean (default: ``False``); ignored

                - ``subcomplex`` -- (default: ``None``) if present,
                  compute the chain complex relative to this subcomplex

                - ``check`` -- boolean (default: ``False``); if ``True``, make
                  sure that the chain complex is actually a chain complex:
                  the differentials are composable and their product is zero.

                The normalized chain complex of a simplicial set is isomorphic
                to the chain complex obtained by modding out by degenerate
                simplices, and the latter is what is actually constructed
                here.

                EXAMPLES::

                    sage: # needs sage.graphs
                    sage: W = simplicial_sets.Sphere(1).wedge(simplicial_sets.Sphere(2))
                    sage: W.nondegenerate_simplices()
                    [*, sigma_1, sigma_2]
                    sage: s1 = W.nondegenerate_simplices()[1]
                    sage: L.<t> = LaurentPolynomialRing(QQ)
                    sage: tw = {s1:t}
                    sage: ChC = W.twisted_chain_complex(tw)
                    sage: ChC.differential(1)
                    [-1 + t]
                    sage: ChC.differential(2)
                    [0]

                ::

                    sage: # needs sage.graphs
                    sage: X = simplicial_sets.Torus()
                    sage: C = X.twisted_chain_complex()
                    sage: C.differential(1)
                    [      f2 - 1 f1*f2^-1 - 1       f1 - 1]
                    sage: C.differential(2)
                    [       1 f1*f2^-1]
                    [      f2        1]
                    [      -1       -1]
                    sage: C.differential(3)
                    []

                ::

                    sage: # needs sage.graphs
                    sage: Y = simplicial_sets.RealProjectiveSpace(2)
                    sage: C = Y.twisted_chain_complex()
                    sage: C.differential(1)
                    [-1 + F1]
                    sage: C.differential(2)
                    [1 + F1]
                    sage: C.differential(3)
                    []
                """
                from sage.homology.chain_complex import ChainComplex
                from sage.structure.element import get_coercion_model
                cm = get_coercion_model()

                if twisting_operator:
                    twop = twisting_operator
                else:
                    di = self._canonical_twisting_operator()
                    if hasattr(list(di.values())[0], "lift"):
                        twop = {a: b.lift() for a, b in di.items()}
                    else:
                        twop = di

                def twist(s):
                    if s in twop:
                        return twop[s]
                    if s.dimension() > 1:
                        return twist(self.face(s,s.dimension()))
                    return 1
                base_ring = cm.common_parent(*twop.values())

                if dimensions is None:
                    if not self.cells():  # Empty
                        if cochain:
                            return ChainComplex({-1: matrix(base_ring, 0, 0)},
                                                degree_of_differential=1)
                        return ChainComplex({0: matrix(base_ring, 0, 0)},
                                            degree_of_differential=-1)
                    dimensions = list(range(self.dimension() + 1))
                else:
                    if not isinstance(dimensions, (list, tuple, range)):
                        dimensions = list(range(dimensions - 1, dimensions + 2))
                    else:
                        dimensions = [n for n in dimensions if n >= 0]
                    if not dimensions:
                        # Return the empty chain complex.
                        if cochain:
                            return ChainComplex(base_ring=base_ring, degree=1)
                        else:
                            return ChainComplex(base_ring=base_ring, degree=-1)

                differentials = {}
                # Convert the tuple self._data to a dictionary indexed by the
                # non-degenerate simplices.
                if subcomplex:
                    X = self.quotient(subcomplex)
                    face_data = X.face_data()
                    nondegens = X.nondegenerate_simplices()
                else:
                    face_data = self.face_data()
                    nondegens = self.nondegenerate_simplices()
                # simplices: dictionary indexed by dimension, values the list
                # of non-degenerate simplices in that dimension.
                simplices = {}
                for sigma in nondegens:
                    if sigma.dimension() in simplices:
                        simplices[sigma.dimension()].append(sigma)
                    else:
                        simplices[sigma.dimension()] = [sigma]
                first = dimensions.pop(0)
                if first in simplices:
                    rank = len(simplices[first])
                    current = sorted(simplices[first])
                else:
                    rank = 0
                    current = []
                if augmented and first == 0:
                    differentials[first-1] = matrix(base_ring, 0, 1)
                    differentials[first] = matrix(base_ring, 1, rank,
                                                  [1] * rank)
                else:
                    differentials[first] = matrix(base_ring, 0, rank)

                for d in dimensions:
                    old_rank = rank
                    faces = {_[1]: _[0] for _ in enumerate(current)}
                    if d in simplices:
                        current = sorted(simplices[d])
                        rank = len(current)
                        # old_rank: number of simplices in dimension d-1.
                        # faces: list of simplices in dimension d-1.
                        # rank: number of simplices in dimension d.
                        # current: list of simplices in dimension d.
                        if not faces:
                            differentials[d] = matrix(base_ring, old_rank, rank)
                        else:
                            matrix_data = {}
                            for col, sigma in enumerate(current):
                                sign = 1
                                twists = len(face_data[sigma]) * [1]
                                twists[0] = twist(sigma)
                                for (ch, tau) in zip(twists, face_data[sigma]):
                                    if tau.is_nondegenerate():
                                        row = faces[tau]
                                        if (row, col) in matrix_data:
                                            matrix_data[(row, col)] += sign*ch
                                        else:
                                            matrix_data[(row, col)] = sign*ch
                                    sign *= -1

                            differentials[d] = matrix(base_ring, old_rank,
                                                      rank, matrix_data, sparse=False)

                    else:
                        rank = 0
                        current = []
                        differentials[d] = matrix(base_ring, old_rank, rank, sparse=False)

                if cochain:
                    new_diffs = {}
                    for d in differentials:
                        new_diffs[d-1] = differentials[d].transpose()
                    return ChainComplex(new_diffs, degree_of_differential=1,
                                        check=check)
                return ChainComplex(differentials, degree_of_differential=-1,
                                    check=check)

            def twisted_homology(self, n, reduced=False):
                r"""
                The `n`-th twisted homology module of the simplicial set with respect to
                the abelianization of the fundamental_group.

                It is a module over a polynomial ring, including relations to make some
                variables the multiplicative inverses of others.

                INPUT:

                - ``n`` -- positive integer

                - ``reduced`` -- boolean (default: ``False``); if set to True,
                  the presentation matrix will be reduced

                EXAMPLES::

                    sage: # needs sage.graphs
                    sage: X = simplicial_sets.Sphere(1).wedge(simplicial_sets.Sphere(2))
                    sage: X.twisted_homology(1)
                    Quotient module by Submodule of Ambient free module of rank 0 over the integral domain Multivariate Polynomial Ring in f1, f1inv over Integer Ring
                    Generated by the rows of the matrix:
                    []
                    sage: X.twisted_homology(2)
                    Quotient module by Submodule of Ambient free module of rank 1 over the integral domain Multivariate Polynomial Ring in f1, f1inv over Integer Ring
                    Generated by the rows of the matrix:
                    [f1*f1inv - 1]

                ::

                    sage: # needs sage.graphs
                    sage: Y = simplicial_sets.Torus()
                    sage: Y.twisted_homology(1)
                    Quotient module by Submodule of Ambient free module of rank 5 over the integral domain Multivariate Polynomial Ring in f1, f1inv, f2, f2inv over Integer Ring
                    Generated by the rows of the matrix:
                    [           1            0            0            0            0]
                    [           0            1            0            0            0]
                    [           0            0            1            0            0]
                    [           0            0            0            1            0]
                    [           0            0            0            0            1]
                    [f1*f1inv - 1            0            0            0            0]
                    [           0 f1*f1inv - 1            0            0            0]
                    [           0            0 f1*f1inv - 1            0            0]
                    [           0            0            0 f1*f1inv - 1            0]
                    [           0            0            0            0 f1*f1inv - 1]
                    [f2*f2inv - 1            0            0            0            0]
                    [           0 f2*f2inv - 1            0            0            0]
                    [           0            0 f2*f2inv - 1            0            0]
                    [           0            0            0 f2*f2inv - 1            0]
                    [           0            0            0            0 f2*f2inv - 1]
                    sage: Y.twisted_homology(2)
                    Quotient module by Submodule of Ambient free module of rank 0 over the integral domain Multivariate Polynomial Ring in f1, f1inv, f2, f2inv over Integer Ring
                    Generated by the rows of the matrix:
                    []
                    sage: Y.twisted_homology(1, reduced=True)
                    Quotient module by Submodule of Ambient free module of rank 5 over the integral domain Multivariate Polynomial Ring in f1, f1inv, f2, f2inv over Integer Ring
                    Generated by the rows of the matrix:
                    [1 0 0 0 0]
                    [0 1 0 0 0]
                    [0 0 1 0 0]
                    [0 0 0 1 0]
                    [0 0 0 0 1]

                TESTS::

                    sage: # needs sage.graphs
                    sage: X = simplicial_sets.PresentationComplex(groups.presentation.FGAbelian((3,2)))
                    sage: TW2 = X.twisted_homology(2, reduced=True)
                    sage: M = TW2.relations_matrix()
                    sage: from sage.libs.singular.function import singular_function
                    sage: vdim = singular_function("vdim")
                    sage: vdim(M.T, ring=M.base_ring())
                    // ** considering the image in Q[...]
                    // ** _ is no standard basis
                    5
                    sage: X.universal_cover().homology(2)
                    Z^5
                    sage: from sage.libs.singular.function import singular_function
                """
                from sage.libs.singular.function import singular_function
                from sage.libs.singular.option import opt_verb
                opt_verb['not_warn_sb'] = True
                singstd = singular_function("std")
                singsyz = singular_function("syz")
                singred = singular_function("reduce")
                singlift = singular_function("lift")
                G, d = self._universal_cover_dict()
                abelG, R, I, images = G.abelianization_to_algebra(ZZ)
                CC = self.twisted_chain_complex()
                M1 = CC.differential(n).T
                M2 = CC.differential(n + 1).T

                def convert_to_polynomial(p):
                    if hasattr(p, "lift"):
                        return p.lift()._as_extended_polynomial()
                    else:
                        return p._as_extended_polynomial()
                M1 = M1.apply_map(convert_to_polynomial)
                M2 = M2.apply_map(convert_to_polynomial)
                RP = R._extended_ring
                IP = RP.ideal([convert_to_polynomial(g) for g in I])
                JP = R._extended_ring_ideal
                GB = (IP+JP).groebner_basis()
                GBI = RP.ideal(GB)

                def reduce_laurent(a):
                    return singred(a, GBI, ring=RP)

                def group_to_polynomial(el, RP):
                    res = RP.one()
                    for a in el.Tietze():
                        if a > 0:
                            res *= RP.gen(2*a-2)
                        else:
                            res *= RP.gen(-2*a-1)
                    return res

                def mkernel(M):
                    if M.nrows() == 0:
                        return matrix(M.base_ring(), 0, 0)
                    if M.ncols() == 0:
                        return M.T
                    res = M
                    n = res.ncols()
                    for g in (IP+JP).gens():
                        res = res.stack(g*identity_matrix(n))
                    syz = matrix(singsyz(res.T, ring=res.base_ring())).T
                    trimmed = syz.T.submatrix(0, 0, syz.ncols(), M.nrows())
                    trimmed = trimmed.apply_map(reduce_laurent)
                    to_delete = [i for (i, r) in enumerate(trimmed.rows()) if not r]
                    return trimmed.delete_rows(to_delete)

                def lift_to_submodule(S, M):
                    if S.nrows() == 0 or S.ncols() == 0:
                        return S
                    res = M
                    for g in GB:
                        res = res.stack(g*identity_matrix(M.ncols()))
                    singres = matrix(singlift(res.T, S.T,ring=res.base_ring()))
                    return singres.submatrix(0, 0, M.nrows(), S.nrows())

                def mgb(M):
                    if M.nrows() == 0 or M.ncols() == 0:
                        return M
                    res = M
                    for g in GB:
                        res = res.stack(g*identity_matrix(M.ncols()))
                    sres = matrix(singstd(res.T, ring=RP))
                    to_delete = [i for i, r in enumerate(sres.apply_map(reduce_laurent)) if not r]
                    return sres.delete_rows(to_delete)
                    M2 = border_matrix(n+1)
                if M1.nrows() == 0:
                    opt_verb.reset_default()
                    return (RP**0).quotient_module([])
                K = mkernel(M1)
                DK = mkernel(K)
                if M2.nrows() > 0:
                    S = lift_to_submodule(M2, K)
                    if S.nrows() > 0 and S.ncols() > 0:
                        res = mgb(DK.stack(S.T))
                    else:
                        res = DK
                else:
                    res = mgb(DK)
                resmat = mgb(res.apply_map(reduce_laurent))
                AM = RP ** K.nrows()
                if resmat.ncols() == 0:
                    SM = AM.submodule([])
                    opt_verb.reset_default()
                    return AM.quotient_module(SM)
                for g in (IP+JP).gens():
                    resmat = resmat.stack(g * identity_matrix(resmat.ncols()))
                if reduced:
                    resmat = matrix(singstd(resmat.T, ring=RP))
                SM = AM.submodule(resmat)
                opt_verb.reset_default()
                return AM.quotient_module(SM)

            def is_simply_connected(self):
                """
                Return ``True`` if this pointed simplicial set is simply connected.

                .. WARNING::

                    Determining simple connectivity is not always
                    possible, because it requires determining when a
                    group, as given by generators and relations, is
                    trivial. So this conceivably may give a false
                    negative in some cases.

                EXAMPLES::

                    sage: # needs sage.graphs sage.groups
                    sage: T = simplicial_sets.Torus()
                    sage: T.is_simply_connected()
                    False
                    sage: T.suspension().is_simply_connected()
                    True
                    sage: simplicial_sets.KleinBottle().is_simply_connected()
                    False

                    sage: # needs sage.graphs
                    sage: S2 = simplicial_sets.Sphere(2)
                    sage: S3 = simplicial_sets.Sphere(3)
                    sage: (S2.wedge(S3)).is_simply_connected()                          # needs sage.groups
                    True
                    sage: X = S2.disjoint_union(S3)
                    sage: X = X.set_base_point(X.n_cells(0)[0])
                    sage: X.is_simply_connected()
                    False

                    sage: C3 = groups.misc.MultiplicativeAbelian([3])                   # needs sage.graphs sage.groups
                    sage: BC3 = simplicial_sets.ClassifyingSpace(C3)                    # needs sage.graphs sage.groups
                    sage: BC3.is_simply_connected()                                     # needs sage.graphs sage.groups
                    False
                """
                if not self.is_connected():
                    return False
                try:
                    if not self.is_pointed():
                        space = self.set_base_point(self.n_cells(0)[0])
                    else:
                        space = self
                    return bool(space.fundamental_group().IsTrivial())
                except AttributeError:
                    try:
                        return space.fundamental_group().order() == 1
                    except (NotImplementedError, RuntimeError):
                        # I don't know of any simplicial sets for which the
                        # code reaches this point, but there are certainly
                        # groups for which these errors are raised. 'IsTrivial'
                        # works for all of the examples I've seen, though.
                        raise ValueError('unable to determine if the fundamental '
                                         'group is trivial')

            def connectivity(self, max_dim=None):
                """
                Return the connectivity of this pointed simplicial set.

                INPUT:

                - ``max_dim`` -- specify a maximum dimension through
                  which to check. This is required if this simplicial
                  set is simply connected and not finite.

                The dimension of the first nonzero homotopy group. If
                simply connected, this is the same as the dimension of
                the first nonzero homology group.

                .. WARNING::

                   See the warning for the :meth:`is_simply_connected` method.

                The connectivity of a contractible space is ``+Infinity``.

                EXAMPLES::

                    sage: # needs sage.graphs sage.groups
                    sage: simplicial_sets.Sphere(3).connectivity()
                    2
                    sage: simplicial_sets.Sphere(0).connectivity()
                    -1
                    sage: K = simplicial_sets.Simplex(4)
                    sage: K = K.set_base_point(K.n_cells(0)[0])
                    sage: K.connectivity()
                    +Infinity
                    sage: X = simplicial_sets.Torus().suspension(2)
                    sage: X.connectivity()
                    2

                    sage: C2 = groups.misc.MultiplicativeAbelian([2])                   # needs sage.graphs sage.groups
                    sage: BC2 = simplicial_sets.ClassifyingSpace(C2)                    # needs sage.graphs sage.groups
                    sage: BC2.connectivity()                                            # needs sage.graphs sage.groups
                    0
                """
                if not self.is_connected():
                    return Integer(-1)
                if not self.is_simply_connected():
                    return Integer(0)
                if max_dim is None:
                    if self.is_finite():
                        max_dim = self.dimension()
                    else:
                        # Note: at the moment, this will never be reached,
                        # because our only examples (so far) of infinite
                        # simplicial sets are not simply connected.
                        raise ValueError('this simplicial set may be infinite, '
                                         'so specify a maximum dimension through '
                                         'which to check')

                H = self.homology(range(2, max_dim + 1))
                for i in range(2, max_dim + 1):
                    if i in H and H[i].order() != 1:
                        return i-1
                return Infinity

        class Finite(CategoryWithAxiom):
            class ParentMethods:

                def unset_base_point(self):
                    """
                    Return a copy of this simplicial set in which the base point has
                    been forgotten.

                    EXAMPLES::

                        sage: # needs sage.graphs
                        sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
                        sage: v = AbstractSimplex(0, name='v_0')
                        sage: w = AbstractSimplex(0, name='w_0')
                        sage: e = AbstractSimplex(1)
                        sage: Y = SimplicialSet({e: (v, w)}, base_point=w)
                        sage: Y.is_pointed()
                        True
                        sage: Y.base_point()
                        w_0
                        sage: Z = Y.unset_base_point()
                        sage: Z.is_pointed()
                        False
                    """
                    from sage.topology.simplicial_set import SimplicialSet
                    return SimplicialSet(self.face_data())

                def fat_wedge(self, n):
                    """
                    Return the `n`-th fat wedge of this pointed simplicial set.

                    This is the subcomplex of the `n`-fold product `X^n`
                    consisting of those points in which at least one
                    factor is the base point. Thus when `n=2`, this is the
                    wedge of the simplicial set with itself, but when `n`
                    is larger, the fat wedge is larger than the `n`-fold
                    wedge.

                    EXAMPLES::

                        sage: # needs sage.graphs
                        sage: S1 = simplicial_sets.Sphere(1)
                        sage: S1.fat_wedge(0)
                        Point
                        sage: S1.fat_wedge(1)
                        S^1
                        sage: S1.fat_wedge(2).fundamental_group()                       # needs sage.groups
                        Finitely presented group < e0, e1 |  >
                        sage: S1.fat_wedge(4).homology()                                # needs sage.modules
                        {0: 0, 1: Z x Z x Z x Z, 2: Z^6, 3: Z x Z x Z x Z}
                    """
                    from sage.topology.simplicial_set_examples import Point
                    if n == 0:
                        return Point()
                    if n == 1:
                        return self
                    return self.product(*[self]*(n-1)).fat_wedge_as_subset()

                def smash_product(self, *others):
                    """
                    Return the smash product of this simplicial set with ``others``.

                    INPUT:

                    - ``others`` -- one or several simplicial sets

                    EXAMPLES::

                        sage: # needs sage.graphs sage.groups
                        sage: S1 = simplicial_sets.Sphere(1)
                        sage: RP2 = simplicial_sets.RealProjectiveSpace(2)
                        sage: X = S1.smash_product(RP2)
                        sage: X.homology(base_ring=GF(2))                               # needs sage.modules
                        {0: Vector space of dimension 0 over Finite Field of size 2,
                         1: Vector space of dimension 0 over Finite Field of size 2,
                         2: Vector space of dimension 1 over Finite Field of size 2,
                         3: Vector space of dimension 1 over Finite Field of size 2}

                        sage: T = S1.product(S1)                                        # needs sage.graphs sage.groups
                        sage: X = T.smash_product(S1)                                   # needs sage.graphs sage.groups
                        sage: X.homology(reduced=False)                                 # needs sage.graphs sage.groups sage.modules
                        {0: Z, 1: 0, 2: Z x Z, 3: Z}
                    """
                    from sage.topology.simplicial_set_constructions import SmashProductOfSimplicialSets_finite
                    return SmashProductOfSimplicialSets_finite((self,) + others)
