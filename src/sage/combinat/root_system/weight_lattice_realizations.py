# sage_setup: distribution = sagemath-modules
"""
Weight lattice realizations
"""
# ****************************************************************************
#       Copyright (C) 2007-2012 Nicolas M. Thiery <nthiery at users.sf.net>
#
#       (with contributions of many others)
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc_c import prod
from sage.categories.category_types import Category_over_base_ring
from sage.sets.family import Family
from .root_lattice_realizations import RootLatticeRealizations
from sage.rings.rational_field import QQ


class WeightLatticeRealizations(Category_over_base_ring):
    r"""
    The category of weight lattice realizations over a given base ring.

    A *weight lattice realization* `L` over a base ring `R` is a free
    module (or vector space if `R` is a field) endowed with an embedding
    of the root lattice of some root system. By restriction, this
    embedding defines an embedding of the root lattice of this root
    system, which makes `L` a root lattice realization.

    Typical weight lattice realizations over `\ZZ` include the weight
    lattice, and ambient lattice. Typical weight lattice realizations
    over `\QQ` include the weight space, and ambient space.

    To describe the embedding, a weight lattice realization must
    implement a method
    :meth:`~RootLatticeRealizations.ParentMethods.fundamental_weight`(i)
    returning for each `i` in the index set the image of the fundamental
    weight `\Lambda_i` under the embedding.

    In order to be a proper root lattice realization, a weight lattice
    realization should also implement the scalar product with the coroot
    lattice; on the other hand, the embedding of the simple roots is
    given for free.

    .. SEEALSO::

        - :class:`~sage.combinat.root_system.root_system.RootSystem`
        - :class:`~sage.combinat.root_system.root_lattice_realizations.RootLatticeRealizations`
        - :class:`~sage.combinat.root_system.weight_space.WeightSpace`
        - :class:`~sage.combinat.root_system.ambient_space.AmbientSpace`

    EXAMPLES:

    Here, we consider the root system of type `A_7`, and embed the weight
    lattice element `x = \Lambda_1 + 2 \Lambda_3` in several root lattice
    realizations::

        sage: R = RootSystem(["A",7])
        sage: Lambda = R.weight_lattice().fundamental_weights()
        sage: x = Lambda[2] + 2 * Lambda[5]

        sage: L = R.weight_space()
        sage: L(x)
        Lambda[2] + 2*Lambda[5]

        sage: L = R.ambient_lattice()
        sage: L(x)
        (3, 3, 2, 2, 2, 0, 0, 0)

    We embed the weight space element `x = \Lambda_1 + 1/2 \Lambda_3` in
    the ambient space::

        sage: Lambda = R.weight_space().fundamental_weights()
        sage: x = Lambda[2] + 1/2 * Lambda[5]

        sage: L = R.ambient_space()
        sage: L(x)
        (3/2, 3/2, 1/2, 1/2, 1/2, 0, 0, 0)

    Of course, one can't embed the weight space in the ambient lattice::

        sage: L = R.ambient_lattice()
        sage: L(x)
        Traceback (most recent call last):
        ...
        TypeError: do not know how to make x (= Lambda[2] + 1/2*Lambda[5])
        an element of self (=Ambient lattice of the Root system of type ['A', 7])

    If `K_1` is a subring of `K_2`, then one could in theory have an
    embedding from the weight space over `K_1` to any weight lattice
    realization over `K_2`; this is not implemented::

        sage: K1 = QQ
        sage: K2 = QQ['q']
        sage: L = R.ambient_space(K2)

        sage: Lambda = R.weight_space(K2).fundamental_weights()
        sage: L(Lambda[1])
        (1, 0, 0, 0, 0, 0, 0, 0)

        sage: Lambda = R.weight_space(K1).fundamental_weights()
        sage: L(Lambda[1])
        Traceback (most recent call last):
        ...
        TypeError: do not know how to make x (= Lambda[1]) an element
        of self (=Ambient space of the Root system of type ['A', 7])
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.combinat.root_system.weight_lattice_realizations import WeightLatticeRealizations
            sage: WeightLatticeRealizations(QQ).super_categories()
            [Category of root lattice realizations over Rational Field]
        """
        return [RootLatticeRealizations(self.base_ring())]

    class ParentMethods:

        @abstract_method
        def fundamental_weight(self, i):
            r"""
            Return the `i`-th fundamental weight.

            INPUT:

            - ``i`` -- an element of the index set

            By a slight notational abuse, for an affine type this method
            should also accept ``'delta'`` as input, and return the image
            of `\delta` of the extended weight lattice in this
            realization.

            This should be overridden by any subclass, and typically
            be implemented as a cached method for efficiency.

            EXAMPLES::

                sage: L = RootSystem(["A",3]).ambient_lattice()
                sage: L.fundamental_weight(1)
                (1, 0, 0, 0)

                sage: L = RootSystem(["A",3,1]).weight_lattice(extended=True)
                sage: L.fundamental_weight(1)
                Lambda[1]
                sage: L.fundamental_weight("delta")
                delta

            TESTS::

                sage: super(sage.combinat.root_system.weight_space.WeightSpace, L).fundamental_weight(1)
                Traceback (most recent call last):
                ...
                NotImplementedError: <abstract method fundamental_weight at ...>
            """

        def is_extended(self):
            """
            Return whether this is a realization of the extended weight lattice.

            .. SEEALSO:: :class:`sage.combinat.root_system.weight_space.WeightSpace`

            EXAMPLES::

                sage: RootSystem(["A",3,1]).weight_lattice().is_extended()
                False
                sage: RootSystem(["A",3,1]).weight_lattice(extended=True).is_extended()
                True

            This method is irrelevant for finite root systems, since the
            weight lattice need not be extended to ensure that the root
            lattice embeds faithfully::

                sage: RootSystem(["A",3]).weight_lattice().is_extended()
                False
            """
            return False

        def __init_extra__(self):
            r"""
            Registers the embedding of the weight lattice into ``self``.

            Also registers the embedding of the weight space over the same
            base field `K` into ``self`` if `K` is not `\ZZ`.

            If ``self`` is a realization of the extended weight lattice,
            then the embeddings from the extended weight space/lattices
            are registered instead.

            EXAMPLES:

            We embed the fundamental weight `\Lambda_1` of the weight
            lattice in the ambient lattice::

                sage: R = RootSystem(["A",3])
                sage: Lambda = R.root_lattice().simple_roots()
                sage: L = R.ambient_space()
                sage: L(Lambda[2])
                (0, 1, -1, 0)

            .. NOTE::

                More examples are given in :class:`WeightLatticeRealizations`;
                The embeddings are systematically tested in
                :meth:`_test_weight_lattice_realization`.
            """
            from sage.rings.integer_ring import ZZ
            from .weight_space import WeightSpace
            K = self.base_ring()
            # If self is the root lattice or the root space, we don't want
            # to register its trivial embedding into itself. This builds
            # the domains from which we want to register an embedding.
            domains = []
            if not isinstance(self, WeightSpace) or K is not ZZ:
                domains.append(self.root_system.weight_lattice(extended=self.is_extended()))
            if not isinstance(self, WeightSpace):
                domains.append(self.root_system.weight_space(K,extended=self.is_extended()))
            # Build and register the embeddings
            for domain in domains:
                domain.module_morphism(self.fundamental_weight,
                                       codomain=self
                                       ).register_as_coercion()

        def _test_weight_lattice_realization(self, **options):
            """
            Run sanity checks on this weight lattice realization.

            - scalar products between the fundamental weights and simple coroots
            - embeddings from the weight lattice and weight space
            - rho, highest_root, ...

            .. SEEALSO:: :class:`TestSuite`

            EXAMPLES::

                sage: RootSystem(['A',3]).weight_lattice()._test_weight_lattice_realization()       # needs sage.graphs
            """
            from sage.rings.integer_ring import ZZ
            tester = self._tester(**options)

            try:
                Lambda = self.fundamental_weights()
                alphacheck = self.simple_coroots()
            except ImportError:
                return

            tester.assertEqual(tuple(Lambda.keys()), self.index_set())

            # Check the consistency between simple_root and simple_roots
            for i in self.index_set():
                tester.assertEqual(self.fundamental_weight(i), Lambda[i])

            # Check the embeddings from:
            # - the weight lattice
            # - the weight space over the same base ring
            #
            # For an affine root system, this will check the embedding of
            # the extended ones, and also of the non extended ones if this
            # realization is not extended
            domains = [self.root_system.weight_space(base_ring, extended=extended)
                       for base_ring in set([ZZ, self.base_ring()])
                       for extended in set([self.cartan_type().is_affine(), self.is_extended()])]
            for domain in domains:
                tester.assertIsNot(self._internal_coerce_map_from(domain), None)
                for i in self.index_set():
                    # This embedding maps fundamental weights to fundamental weights
                    tester.assertEqual(self(domain.fundamental_weight(i)), Lambda[i])
                if self.cartan_type().is_affine():
                    tester.assertEqual(self(domain.null_root()), self.null_root())
                    if self.is_extended():
                        a = self.cartan_type().col_annihilator()
                        # This could be an over specification; we
                        # could imagine realizations of the extended
                        # weight lattice where the null root would not
                        # be a (multiple of) basis element.
                        tester.assertEqual(self.null_root(), self.term("delta", a[0]))
                    for i in self.index_set():
                        # The level of the fundamental weights is consistent
                        tester.assertEqual(domain.fundamental_weight(i).level(), Lambda[i].level())

            # Check that the fundamental weights form the dual basis of the simple coroots
            for i in self.index_set():
                assert (Lambda[i].is_dominant())
                for j in self.index_set():
                    tester.assertEqual(Lambda[j].scalar(alphacheck[i]), (1 if i == j else 0))

            tester.assertTrue(self.rho().is_dominant())
            if self.root_system.is_finite() and self.root_system.is_irreducible():
                try:
                    tester.assertTrue(self.highest_root().is_dominant())
                except ImportError:  # when sage.graphs is missing (Dynkin diagram)
                    pass

        @cached_method
        def fundamental_weights(self):
            r"""
            Return the family `(\Lambda_i)_{i\in I}` of the fundamental weights.

            EXAMPLES::

                sage: e = RootSystem(['A',3]).ambient_lattice()
                sage: f = e.fundamental_weights()
                sage: [f[i] for i in [1,2,3]]
                [(1, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 0)]
            """
            return Family(self.index_set(), self.fundamental_weight)
            # It would be nice to give this family a nice name with
            # ``rename``, but this currently break some doctests.

        @cached_method
        def simple_root(self, i):
            r"""
            Return the `i`-th simple root.

            This default implementation takes the `i`-th simple root in
            the weight lattice and embeds it in ``self``.

            EXAMPLES:

            Since all the weight lattice realizations in Sage currently
            implement a ``simple_root`` method, we have to call this one by
            hand::

                sage: from sage.combinat.root_system.weight_lattice_realizations import WeightLatticeRealizations
                sage: simple_root = WeightLatticeRealizations(QQ).parent_class.simple_root.f
                sage: L = RootSystem("A3").ambient_space()
                sage: simple_root(L, 1)                                                 # needs sage.graphs
                (1, -1, 0, 0)
                sage: simple_root(L, 2)                                                 # needs sage.graphs
                (0, 1, -1, 0)
                sage: simple_root(L, 3)                                                 # needs sage.graphs
                (1, 1, 2, 0)

            Note that this last root differs from the one implemented in
            ``L`` by a multiple of the vector ``(1,1,1,1)``::

                sage: L.simple_roots()                                                  # needs sage.graphs
                Finite family {1: (1, -1, 0, 0), 2: (0, 1, -1, 0), 3: (0, 0, 1, -1)}

            This is a harmless artefact of the `SL` versus `GL`
            interpretation of type `A`; see the thematic tutorial on Lie
            Methods and Related Combinatorics in Sage for details.
            """
            if i not in self.index_set():
                raise ValueError("{} is not in the index set".format(i))
            alphai = self.root_system.weight_lattice().simple_root(i)
            # Note: it would be nicer to just return ``self(alpha[i])``,
            # However the embedding from the weight lattice is defined
            # after the embedding from the root lattice, and the later
            # uses the simple roots. So we compute that embedding by hand.
            Lambda = self.fundamental_weights()
            return self.linear_combination( (Lambda[j], c) for j,c in alphai )

        @cached_method
        def rho(self):
            """
            EXAMPLES::

                sage: RootSystem(['A',3]).ambient_lattice().rho()
                (3, 2, 1, 0)
            """
            return sum(self.fundamental_weights())

        def reduced_word_of_alcove_morphism(self, f):
            r"""
            Return the reduced word of an alcove morphism.

            INPUT:

            - ``f`` -- a linear map from ``self`` to ``self`` which
              preserves alcoves

            Let `A` be the fundamental alcove. This returns a reduced word
            `i_1, \ldots, i_k` such that the affine Weyl group element `w =
            s_{i_1} \circ \cdots \circ s_{i_k}` maps the alcove `f(A)` back
            to `A`. In other words, the alcove walk `i_1, \ldots, i_k` brings
            the fundamental alcove to the corresponding translated alcove.

            Let us throw in a bit of context to explain the main use
            case.  It is customary to realize the alcove picture in
            the coroot or coweight lattice `R^\vee`. The extended
            affine Weyl group is then the group of linear maps on
            `R^\vee` which preserve the alcoves. By
            [Kac "Infinite-dimensional Lie algebra", Proposition 6.5]
            the affine Weyl group is the semidirect product of the
            associated finite Weyl group and the group of translations
            in the coroot lattice (the extended affine Weyl group uses
            the coweight lattice instead). In other words, an element
            of the extended affine Weyl group admits a unique
            decomposition of the form:

            .. MATH:: f = d w ,

            where `w` is in the Weyl group, and `d` is a function which
            maps the fundamental alcove to itself. As `d` permutes the
            walls of the fundamental alcove, it permutes accordingly the
            corresponding simple roots, which induces an automorphism of
            the Dynkin diagram.

            This method returns a reduced word for `w`, whereas the method
            :meth:`dynkin_diagram_automorphism_of_alcove_morphism` returns
            `d` as a permutation of the nodes of the Dynkin diagram.

            Nota bene: recall that the coroot (resp. coweight) lattice is
            implemented as the root (resp weight) lattice of the dual root
            system. Hence, this method is implemented for weight lattice
            realizations, but in practice is most of the time used on the
            dual side.

            EXAMPLES:

            We start with type `A` which is simply laced; hence we do not
            have to worry about the distinction between the weight and
            coweight lattice::

                sage: R = RootSystem(["A",2,1]).weight_lattice()
                sage: alpha = R.simple_roots()                                          # needs sage.graphs
                sage: Lambda = R.fundamental_weights()

            We consider first translations by elements of the root lattice::

                sage: R.reduced_word_of_alcove_morphism(alpha[0].translation)           # needs sage.graphs
                [1, 2, 1, 0]
                sage: R.reduced_word_of_alcove_morphism(alpha[1].translation)           # needs sage.graphs
                [0, 2, 0, 1]
                sage: R.reduced_word_of_alcove_morphism(alpha[2].translation)           # needs sage.graphs
                [0, 1, 0, 2]

            We continue with translations by elements of the classical
            weight lattice, embedded at level `0`:

                sage: omega1 = Lambda[1] - Lambda[0]
                sage: omega2 = Lambda[2] - Lambda[0]

                sage: R.reduced_word_of_alcove_morphism(omega1.translation)             # needs sage.graphs
                [0, 2]
                sage: R.reduced_word_of_alcove_morphism(omega2.translation)             # needs sage.graphs
                [0, 1]

            The following tests ensure that the code agrees with the tables
            in Kashiwara's private notes on affine quantum algebras (2008).

            TESTS::

                sage: # needs sage.graphs
                sage: R = RootSystem(['A',5,1]).weight_lattice()
                sage: alpha = R.simple_roots()
                sage: Lambda = R.fundamental_weights()
                sage: omega1 = Lambda[1] - Lambda[0]
                sage: R.reduced_word_of_alcove_morphism(omega1.translation)
                [0, 5, 4, 3, 2]
                sage: R.reduced_word_of_alcove_morphism(alpha[0].translation)
                [1, 2, 3, 4, 5, 4, 3, 2, 1, 0]

                sage: # needs sage.graphs
                sage: R = RootSystem(['C',3,1]).weight_lattice()
                sage: alpha = R.simple_roots()
                sage: Lambda = R.fundamental_weights()
                sage: omega1 = 2*(Lambda[1] - Lambda[0])
                sage: omega2 = 2*(Lambda[2] - Lambda[0])
                sage: omega3 = Lambda[3] - Lambda[0]
                sage: R.reduced_word_of_alcove_morphism(omega1.translation)
                [0, 1, 2, 3, 2, 1]
                sage: R.reduced_word_of_alcove_morphism(omega2.translation)
                [0, 1, 0, 2, 1, 3, 2, 1, 3, 2]
                sage: R.reduced_word_of_alcove_morphism(omega3.translation)
                [0, 1, 0, 2, 1, 0]

                sage: # needs sage.libs.gap
                sage: W = WeylGroup(['C',3,1])
                sage: s = W.simple_reflections()
                sage: w = s[0]*s[1]*s[2]*s[3]*s[2]
                sage: W.from_reduced_word(R.reduced_word_of_alcove_morphism(omega2.translation)) == w*w                 # needs sage.graphs
                True
                sage: w = s[0]*s[1]*s[2]*s[0]*s[1]*s[0]
                sage: W.from_reduced_word(R.reduced_word_of_alcove_morphism(omega3.translation)) == w                   # needs sage.graphs
                True

                sage: # needs sage.graphs
                sage: R = RootSystem(['D',4,1]).weight_lattice()
                sage: Lambda = R.fundamental_weights()
                sage: omega1 = Lambda[1] - Lambda[0]
                sage: omega2 = Lambda[2] - 2*Lambda[0]
                sage: omega3 = Lambda[3] - Lambda[0]
                sage: omega4 = Lambda[4] - Lambda[0]
                sage: R.reduced_word_of_alcove_morphism(omega1.translation)
                [0, 2, 3, 4, 2, 0]
                sage: R.reduced_word_of_alcove_morphism(omega2.translation)
                [0, 2, 1, 3, 2, 4, 2, 1, 3, 2]
                sage: R.reduced_word_of_alcove_morphism(omega3.translation)
                [0, 2, 1, 4, 2, 0]
                sage: R.reduced_word_of_alcove_morphism(omega4.translation)
                [0, 2, 1, 3, 2, 0]

                sage: # needs sage.libs.gap
                sage: W = WeylGroup(['D',4,1])
                sage: s = W.simple_reflections()
                sage: w = s[0]*s[2]*s[3]*s[4]*s[2]
                sage: w1= s[1]*s[2]*s[3]*s[4]*s[2]
                sage: W.from_reduced_word(R.reduced_word_of_alcove_morphism(omega2.translation)) == w*w1                # needs sage.graphs
                True

                sage: R = RootSystem(['D',5,1]).weight_lattice()
                sage: Lambda = R.fundamental_weights()
                sage: omega1 = Lambda[1] - Lambda[0]
                sage: omega2 = Lambda[2] - 2*Lambda[0]
                sage: R.reduced_word_of_alcove_morphism(omega1.translation)             # needs sage.graphs
                [0, 2, 3, 4, 5, 3, 2, 0]

                sage: # needs sage.libs.gap
                sage: W = WeylGroup(['D',5,1])
                sage: s = W.simple_reflections()
                sage: w = s[0]*s[2]*s[3]*s[4]*s[5]*s[3]*s[2]
                sage: w1= s[1]*s[2]*s[3]*s[4]*s[5]*s[3]*s[2]
                sage: W.from_reduced_word(R.reduced_word_of_alcove_morphism(omega2.translation)) == w*w1                # needs sage.graphs
                True
            """
            return f(self.rho()).reduced_word()

        def dynkin_diagram_automorphism_of_alcove_morphism(self, f):
            r"""
            Return the Dynkin diagram automorphism induced by an alcove morphism.

            INPUT:

            - ``f`` -- a linear map from ``self`` to ``self`` which preserves alcoves

            This method returns the Dynkin diagram automorphism for
            the decomposition `f = d w` (see
            :meth:`reduced_word_of_alcove_morphism`), as a dictionary
            mapping elements of the index set to itself.

            EXAMPLES::

                sage: R = RootSystem(["A",2,1]).weight_lattice()
                sage: alpha = R.simple_roots()                                          # needs sage.graphs
                sage: Lambda = R.fundamental_weights()

            Translations by elements of the root lattice induce a
            trivial Dynkin diagram automorphism::

                sage: # needs sage.graphs sage.libs.gap
                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(alpha[0].translation)
                {0: 0, 1: 1, 2: 2}
                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(alpha[1].translation)
                {0: 0, 1: 1, 2: 2}
                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(alpha[2].translation)
                {0: 0, 1: 1, 2: 2}

            This is no more the case for translations by general
            elements of the (classical) weight lattice at level 0::

                sage: omega1 = Lambda[1] - Lambda[0]
                sage: omega2 = Lambda[2] - Lambda[0]

                sage: # needs sage.graphs sage.libs.gap
                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(omega1.translation)
                {0: 1, 1: 2, 2: 0}
                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(omega2.translation)
                {0: 2, 1: 0, 2: 1}

                sage: # needs sage.graphs sage.libs.gap
                sage: R = RootSystem(['C',2,1]).weight_lattice()
                sage: alpha = R.simple_roots()
                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(alpha[1].translation)
                {0: 2, 1: 1, 2: 0}

                sage: # needs sage.graphs sage.libs.gap
                sage: R = RootSystem(['D',5,1]).weight_lattice()
                sage: Lambda = R.fundamental_weights()
                sage: omega1 = Lambda[1] - Lambda[0]
                sage: omega2 = Lambda[2] - 2*Lambda[0]
                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(omega1.translation)
                {0: 1, 1: 0, 2: 2, 3: 3, 4: 5, 5: 4}
                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(omega2.translation)
                {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}

            Algorithm: computes `w` of the decomposition, and see how
            `f\circ w^{-1}` permutes the simple roots.
            """
            alpha = self.simple_roots()
            w = self.weyl_group().from_reduced_word(self.reduced_word_of_alcove_morphism(f))
            # Now, we have d = f w^-1
            winv = ~w
            assert all(alpha[i].level().is_zero() for i in self.index_set())
            rank_simple_roots = dict( (alpha[i],i) for i in self.index_set())
            permutation = dict()
            for i in self.index_set():
                root = f(winv.action(alpha[i])) # This is d(alpha_i)
                assert root in rank_simple_roots
                permutation[i] = rank_simple_roots[root]
                assert set(permutation.values()), set(self.index_set())
            return permutation

        def reduced_word_of_translation(self, t):
            r"""
            Given an element of the root lattice, this returns a reduced
            word `i_1, \ldots, i_k` such that the Weyl group element `s_{i_1}
            \circ \cdots \circ s_{i_k}` implements the "translation"
            where `x` maps to `x + level(x)*t`. In other words, the alcove walk
            `i_1, \ldots, i_k` brings the fundamental alcove to the
            corresponding translated alcove.

            .. NOTE::

                There are some technical conditions for `t` to actually
                be a translation; those are not tested (TODO: detail).

            EXAMPLES::

                sage: # needs sage.graphs
                sage: R = RootSystem(["A",2,1]).weight_lattice()
                sage: alpha = R.simple_roots()
                sage: R.reduced_word_of_translation(alpha[1])
                [0, 2, 0, 1]
                sage: R.reduced_word_of_translation(alpha[2])
                [0, 1, 0, 2]
                sage: R.reduced_word_of_translation(alpha[0])
                [1, 2, 1, 0]

                sage: R = RootSystem(['D',5,1]).weight_lattice()
                sage: Lambda = R.fundamental_weights()
                sage: omega1 = Lambda[1] - Lambda[0]
                sage: omega2 = Lambda[2] - 2*Lambda[0]
                sage: R.reduced_word_of_translation(omega1)                             # needs sage.graphs
                [0, 2, 3, 4, 5, 3, 2, 0]
                sage: R.reduced_word_of_translation(omega2)                             # needs sage.graphs
                [0, 2, 1, 3, 2, 4, 3, 5, 3, 2, 1, 4, 3, 2]

            A non simply laced case::

                sage: R = RootSystem(["C",2,1]).weight_lattice()
                sage: Lambda = R.fundamental_weights()
                sage: c = R.cartan_type().translation_factors(); c                      # needs sage.graphs
                Finite family {0: 1, 1: 2, 2: 1}
                sage: R.reduced_word_of_translation((Lambda[1]-Lambda[0]) * c[1])       # needs sage.graphs
                [0, 1, 2, 1]
                sage: R.reduced_word_of_translation((Lambda[2]-Lambda[0]) * c[2])       # needs sage.graphs
                [0, 1, 0]

            See also :meth:`_test_reduced_word_of_translation`.

            .. TODO::

                 - Add a picture in the doc
                 - Add a method which, given an element of the classical
                   weight lattice, constructs the appropriate value for t
            """
            return self.reduced_word_of_alcove_morphism(t.translation)

        def _test_reduced_word_of_translation(self, elements=None, **options):
            r"""
            Test the method :meth:`reduced_word_of_translation`.

            INPUT:

            - ``options`` -- any keyword arguments accepted by :meth:`_tester`

            EXAMPLES::

                sage: R = RootSystem(['D',4,1]).weight_lattice()
                sage: R._test_reduced_word_of_translation()                             # needs sage.graphs

            See the documentation for :class:`TestSuite` for more information.
            """
            tester = self._tester(**options)
            if not self.cartan_type().is_affine(): # won't be necessary anymore once root systems are categorified
                return
            try:
                alpha = self.simple_roots()
            except ImportError:  # when sage.graphs is not available
                return
            Lambda = self.fundamental_weights()
            rho = self.rho()
            G = self.dynkin_diagram()
            permutations = []

            # Note: this uses a special set of default elements instead of
            # the usual tester.some_elements(), namely the smallest
            # elements in the weight lattice giving rise to translations
            # preserving the alcoves.
            if elements is None:
                c = self.cartan_type().c()
                elements = [ c[i] * Lambda[i] for i in self.cartan_type().classical().index_set() ]

            # When the null root is zero in this root lattice realization,
            # the roots correspond to the classical roots. We use that to
            # check that w permute the simple roots according to a Dynkin
            # diagram automorphism. This test currently requires the index
            # set to be of the form 0..n
            test_automorphism = self.null_root().is_zero() and set(self.index_set()) == set(i for i in range(len(self.index_set())))
            # dictionary assigning a simple root to its index
            rank_simple_roots = dict( (alpha[i],i) for i in self.index_set() )

            try:
                W = self.weyl_group()
            except ImportError:
                return

            for t in elements:
                t = t - self.base_ring()(t.level()/Lambda[0].level()) * Lambda[0]
                w = W.from_reduced_word(self.reduced_word_of_translation(t))
                if self.null_root().is_zero():
                    # The following formula is only valid when the null root is zero
                    tester.assertEqual(w.action(rho), rho + rho.level()*t)
                    # TODO: fix this formula to take delta into account,
                    # and remove the above condition
                if test_automorphism:
                    permutation = [None for i in self.index_set()]
                    for i in self.index_set():
                        root = w.action(alpha[i])
                        tester.assertIn(root, rank_simple_roots)
                        permutation[i] = rank_simple_roots[root]
                    tester.assertEqual(set(permutation), set(self.index_set()))
                    # It could be nicer to test equality of G and its relabelling
                    for i in self.index_set():
                        for j in self.index_set():
                            tester.assertEqual(G[permutation[i],permutation[j]], G[i,j])
                    permutations.append(permutation)

            if test_automorphism and elements is None: # note: the test on elements is broken
                # Check that, if we start from all fundamental weights, we
                # get the full automorphism group
                # Disabled: this should actually check that one gets all special
                # automorphisms, which are in bijection with the special nodes
                #from sage.groups.perm_gps.permgroup import PermutationGroup
                #P = PermutationGroup([[i+1 for i in permutation] for permutation in permutations])
                #tester.assertEqual(P, G.automorphism_group())
                pass

        def signs_of_alcovewalk(self, walk):
            r"""
            Let walk = `[i_1,\ldots,i_n]` denote an alcove walk starting
            from the fundamental alcove `y_0`, crossing at step 1 the
            wall `i_1`, and so on.

            For each `k`, set `w_k = s_{i_1} \circ s_{i_k}`, and denote
            by `y_k = w_k(y_0)` the alcove reached after `k` steps. Then,
            `y_k` is obtained recursively from `y_{k-1}` by applying the
            following reflection:

            .. MATH::

                  y_k = s_{w_{k-1} \alpha_{i_k}} y_{k-1}.

            The step is said positive if `w_{k-1} \alpha_{i_k}` is a
            negative root (considering `w_{k-1}` as element of the
            classical Weyl group and `\alpha_{i_k}` as a classical
            root) and negative otherwise. The algorithm implemented
            here use the equivalent property::

                .. MATH:: \langle w_{k-1}^{-1} \rho_0, \alpha^\vee_{i_k}\rangle > 0

            Where `\rho_0` is the sum of the classical fundamental
            weights embedded at level 0 in this space (see
            :meth:`rho_classical`), and `\alpha^\vee_{i_k}` is the
            simple coroot associated to `\alpha_{i_k}`.

            This function returns a list of the form `[+1,+1,-1,...]`,
            where the `k`-th entry denotes whether the `k`-th step was
            positive or negative.

            See equation 3.4, of Ram: Alcove walks ..., :arxiv:`math/0601343v1`

            EXAMPLES::

                sage: L = RootSystem(['C',2,1]).weight_lattice()
                sage: L.signs_of_alcovewalk([1,2,0,1,2,1,2,0,1,2])                      # needs sage.libs.gap
                [-1, -1, 1, -1, 1, 1, 1, 1, 1, 1]

                sage: L = RootSystem(['A',2,1]).weight_lattice()
                sage: L.signs_of_alcovewalk([0,1,2,1,2,0,1,2,0,1,2,0])                  # needs sage.libs.gap
                [1, 1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 1]

                sage: L = RootSystem(['B',2,1]).coweight_lattice()
                sage: L.signs_of_alcovewalk([0,1,2,0,1,2])                              # needs sage.libs.gap
                [1, -1, 1, -1, 1, 1]

            .. WARNING::

                This method currently does not work in the weight
                lattice for type BC dual because `\rho_0` does not
                live in this lattice (but an integral multiple of it
                would do the job as well).
            """
            # Below, w is w_{k-1} and we use the fact that, for a root
            # `a` the following are equivalent:
            # - w a is a negative root
            # - <w a, rho^\vee> < 0
            # - <w a^\vee, rho> < 0
            # - <a^\vee, w^-1 rho> < 0
            W = self.weyl_group()
            s = W.simple_reflections()
            alphacheck = self.alphacheck()
            rho0 = self.rho_classical()
            w = W.one()
            signs = []
            for i in walk:
                if (w.action(rho0).scalar(alphacheck[i]) > 0):
                    signs.append(-1)
                else:
                    signs.append(1)
                w = s[i]*w
            return signs

        def rho_classical(self):
            r"""
            Return the embedding at level 0 of `\rho` of the classical lattice.

            EXAMPLES::

                sage: RootSystem(['C',4,1]).weight_lattice().rho_classical()            # needs sage.graphs
                -4*Lambda[0] + Lambda[1] + Lambda[2] + Lambda[3] + Lambda[4]
                sage: L = RootSystem(['D',4,1]).weight_lattice()
                sage: L.rho_classical().scalar(L.null_coroot())                         # needs sage.graphs
                0

            .. WARNING::

                In affine type BC dual, this does not live in the weight lattice::

                    sage: L = CartanType(["BC",2,2]).dual().root_system().weight_space()
                    sage: L.rho_classical()                                             # needs sage.graphs
                    -3/2*Lambda[0] + Lambda[1] + Lambda[2]
                    sage: L = CartanType(["BC",2,2]).dual().root_system().weight_lattice()
                    sage: L.rho_classical()                                             # needs sage.graphs
                    Traceback (most recent call last):
                    ...
                    ValueError: 5 is not divisible by 2
            """
            rho = self.rho()
            Lambda = self.fundamental_weights()
            return rho - Lambda[0] * rho.level() / Lambda[0].level()

        def embed_at_level(self, x, level=1):
            r"""
            Embed the classical weight `x` in the level ``level`` hyperplane.

            This is achieved by translating the straightforward
            embedding of `x` by `c\Lambda_0` for `c` some appropriate
            scalar.

            INPUT:

            - ``x`` -- an element of the corresponding classical weight/ambient lattice
            - ``level`` -- integer or element of the base ring (default: 1)

            EXAMPLES::

                sage: # needs sage.graphs
                sage: L = RootSystem(["B",3,1]).weight_space()
                sage: L0 = L.classical()
                sage: alpha = L0.simple_roots()
                sage: omega = L0.fundamental_weights()
                sage: L.embed_at_level(omega[1], 1)
                Lambda[1]
                sage: L.embed_at_level(omega[2], 1)
                -Lambda[0] + Lambda[2]
                sage: L.embed_at_level(omega[3], 1)
                Lambda[3]
                sage: L.embed_at_level(alpha[1], 1)
                Lambda[0] + 2*Lambda[1] - Lambda[2]
            """
            if not self.classical().is_parent_of(x):
                raise ValueError("x must be an element of the classical type")
            Lambda = self.fundamental_weights()
            result = self.sum_of_terms(x)
            result += Lambda[0] * (level-result.level()) / (Lambda[0].level())
            assert result.level() == level
            return result

        # Should it be a method of highest_weight?
        def weyl_dimension(self, highest_weight):
            r"""
            Return the dimension of the highest weight representation of highest weight ``highest_weight``.

            EXAMPLES::

                sage: RootSystem(['A',3]).ambient_lattice().weyl_dimension([2,1,0,0])
                20
                sage: P = RootSystem(['C',2]).weight_lattice()
                sage: La = P.basis()
                sage: P.weyl_dimension(La[1]+La[2])                                     # needs sage.graphs
                16

                sage: type(RootSystem(['A',3]).ambient_lattice().weyl_dimension([2,1,0,0]))
                <class 'sage.rings.integer.Integer'>
            """
            highest_weight = self(highest_weight)
            if not highest_weight.is_dominant():
                raise ValueError("the highest weight must be dominant")
            rho = self.rho()
            pr = self.coroot_lattice().positive_roots()
            from sage.rings.integer import Integer
            n = prod(((rho+highest_weight).scalar(x) for x in pr), Integer(1))
            d = prod((rho.scalar(x) for x in pr), Integer(1))
            return Integer(n/d)

        @lazy_attribute
        def _inverse_cartan_matrix(self):
            r"""
            Return the inverse Cartan matrix defining ``self``.

            EXAMPLES::

                sage: RootSystem(['A', 3]).ambient_lattice()._inverse_cartan_matrix
                [3/4 1/2 1/4]
                [1/2   1 1/2]
                [1/4 1/2 3/4]
                sage: RootSystem(['G', 2]).weight_lattice()._inverse_cartan_matrix
                [2 3]
                [1 2]
            """
            ret = self.cartan_type().cartan_matrix().inverse()
            ret.set_immutable()
            return ret

        @lazy_attribute
        def _symmetric_form_matrix(self):
            r"""
            Return the matrix for the symmetric form `( | )` in
            the weight lattice basis.

            Let `A` be a symmetrizable Cartan matrix with symmetrizer `D`,.
            This returns the matrix `M^t DA M`, where `M` is dependent upon
            the type given below.

            In finite types, `M` is the inverse of the Cartan matrix.

            In affine types, `M` takes the basis
            `(\Lambda_0, \Lambda_1, \ldots, \Lambda_r, \delta)` to
            `(\alpha_0, \ldots, \alpha_r, \Lambda_0)` where `r` is the
            rank of ``self``.

            This is used in computing the symmetric form for affine
            root systems.

            EXAMPLES::

                sage: P = RootSystem(['B',2]).weight_lattice()
                sage: P._symmetric_form_matrix                                          # needs sage.graphs
                [2 1]
                [1 1]

                sage: P = RootSystem(['C',2]).weight_lattice()
                sage: P._symmetric_form_matrix                                          # needs sage.graphs
                [1 1]
                [1 2]

                sage: P = RootSystem(['C',2,1]).weight_lattice()
                sage: P._symmetric_form_matrix                                          # needs sage.graphs
                [0 0 0 1]
                [0 1 1 1]
                [0 1 2 1]
                [1 1 1 0]

                sage: P = RootSystem(['A',4,2]).weight_lattice()
                sage: P._symmetric_form_matrix                                          # needs sage.graphs
                [  0   0   0 1/2]
                [  0   2   2   1]
                [  0   2   4   1]
                [1/2   1   1   0]
            """
            from sage.matrix.constructor import matrix
            ct = self.cartan_type()
            cm = ct.cartan_matrix()
            if cm.det() != 0:
                diag = matrix.diagonal(cm.symmetrizer())
                return self._inverse_cartan_matrix.transpose() * diag

            if not ct.is_affine():
                raise ValueError("only implemented for affine types when the"
                                 " Cartan matrix is singular")

            r = ct.rank()
            a = ct.a()
            # Determine the change of basis matrix
            # La[0], ..., La[r], delta -> al[0], ..., al[r], La[0]
            M = cm.stack( matrix([1] + [0]*(r-1)) )
            M = matrix.block([[ M, matrix([[1]] + [[0]]*r) ]])
            M = M.inverse()

            if a[0] != 1:
                from sage.rings.rational_field import QQ
                S = matrix([~a[0]]+[0]*(r-1))
                A = cm.symmetrized_matrix().change_ring(QQ).stack(S)
            else:
                A = cm.symmetrized_matrix().stack(matrix([1]+[0]*(r-1)))
            A = matrix.block([[A, matrix([[~a[0]]] + [[0]]*r)]])
            return M.transpose() * A * M

    class ElementMethods:
        def symmetric_form(self, la):
            r"""
            Return the symmetric form of ``self`` with ``la``.

            Return the pairing `( | )` on the weight lattice. See Chapter 6
            in Kac, Infinite Dimensional Lie Algebras for more details.

            .. WARNING::

                For affine root systems, if you are not working in the
                extended weight lattice/space, this may return incorrect
                results.

            EXAMPLES::

                sage: # needs sage.graphs
                sage: P = RootSystem(['C',2]).weight_lattice()
                sage: al = P.simple_roots()
                sage: al[1].symmetric_form(al[1])
                2
                sage: al[1].symmetric_form(al[2])
                -2
                sage: al[2].symmetric_form(al[1])
                -2
                sage: Q = RootSystem(['C',2]).root_lattice()
                sage: alQ = Q.simple_roots()
                sage: all(al[i].symmetric_form(al[j]) == alQ[i].symmetric_form(alQ[j])
                ....:     for i in P.index_set() for j in P.index_set())
                True

                sage: # needs sage.graphs
                sage: P = RootSystem(['C',2,1]).weight_lattice(extended=True)
                sage: al = P.simple_roots()
                sage: al[1].symmetric_form(al[1])
                2
                sage: al[1].symmetric_form(al[2])
                -2
                sage: al[1].symmetric_form(al[0])
                -2
                sage: al[0].symmetric_form(al[1])
                -2
                sage: Q = RootSystem(['C',2,1]).root_lattice()
                sage: alQ = Q.simple_roots()
                sage: all(al[i].symmetric_form(al[j]) == alQ[i].symmetric_form(alQ[j])
                ....:     for i in P.index_set() for j in P.index_set())
                True
                sage: La = P.basis()
                sage: [La['delta'].symmetric_form(al) for al in P.simple_roots()]
                [0, 0, 0]
                sage: [La[0].symmetric_form(al) for al in P.simple_roots()]
                [1, 0, 0]

                sage: P = RootSystem(['C',2,1]).weight_lattice()
                sage: Q = RootSystem(['C',2,1]).root_lattice()
                sage: al = P.simple_roots()                                             # needs sage.graphs
                sage: alQ = Q.simple_roots()                                            # needs sage.graphs
                sage: all(al[i].symmetric_form(al[j]) == alQ[i].symmetric_form(alQ[j])  # needs sage.graphs
                ....:     for i in P.index_set() for j in P.index_set())
                True

            The result of `(\Lambda_0 | \alpha_0)` should be `1`, however we
            get `0` because we are not working in the extended weight
            lattice::

                sage: La = P.basis()
                sage: [La[0].symmetric_form(al) for al in P.simple_roots()]             # needs sage.graphs
                [0, 0, 0]

            TESTS:

            We check that `A_{2n}^{(2)}` has 3 different root lengths::

                sage: P = RootSystem(['A',4,2]).weight_lattice()
                sage: al = P.simple_roots()                                             # needs sage.graphs
                sage: [al[i].symmetric_form(al[i]) for i in P.index_set()]              # needs sage.graphs
                [2, 4, 8]

            Check that :issue:`31410` is fixed, and the symmetric form
            computed on the weight space is the same as the symmetric
            form computed on the root space::

                sage: def s1(ct):
                ....:    L = RootSystem(ct).weight_space()
                ....:    P = L.positive_roots()
                ....:    rho = L.rho()
                ....:    return [beta.symmetric_form(rho) for beta in P]

                sage: def s2(ct):
                ....:    R = RootSystem(ct).root_space()
                ....:    P = R.positive_roots()
                ....:    rho = 1/2*sum(P)
                ....:    return [beta.symmetric_form(rho) for beta in P]

                sage: all(s1(ct) == s2(ct)                                              # needs sage.graphs
                ....:     for ct in CartanType.samples(finite=True, crystallographic=True))
                True
            """
            P = self.parent()
            ct = P.cartan_type()
            sym = P._symmetric_form_matrix

            if ct.is_finite():
                iset = P.index_set()
            else:
                iset = P.index_set() + ('delta',)

            return sum(cl*sym[iset.index(ml),iset.index(mr)]*cr
                       for ml, cl in self for mr, cr in la)

        #    # This should be in a method to_weight_lattice()
        #    alphac = self.simple_coroots()
        #    Lambda = self.fundamental_weights()
        #    assert t == self.plus(t.scalar(alphac[i]) * Lambda[i] for i in self.index_set())
        #    t = self.plus( t.scalar(alphac[i]) * c[i] * Lambda[i] for i in self.index_set())

        def to_weight_space(self, base_ring=None):
            r"""
            Map ``self`` to the weight space.

            .. WARNING::

                Implemented for finite Cartan type.

            EXAMPLES::

                sage: b = CartanType(['B',2]).root_system().ambient_space().from_vector(vector([1,-2])); b
                (1, -2)
                sage: b.to_weight_space()
                3*Lambda[1] - 4*Lambda[2]
                sage: b = CartanType(['B',2]).root_system().ambient_space().from_vector(vector([1/2,0])); b
                (1/2, 0)
                sage: b.to_weight_space()
                1/2*Lambda[1]
                sage: b.to_weight_space(ZZ)
                Traceback (most recent call last):
                ...
                TypeError: no conversion of this rational to integer
                sage: b = CartanType(['G',2]).root_system().ambient_space().from_vector(vector([4,-5,1])); b
                (4, -5, 1)
                sage: b.to_weight_space()
                -6*Lambda[1] + 5*Lambda[2]
            """
            L = self.parent()
            if base_ring is None:
                base_ring = L.base_ring()

            wt_space = L.root_system.weight_space(base_ring)
            simple_coroots = L.simple_coroots()
            return wt_space.sum_of_terms(((i, base_ring(self.scalar(ac)))
                                          for i, ac in simple_coroots.items()),
                                         distinct=True)

        @cached_method
        def _to_root_vector(self):
            r"""
            Helper method to express ``self`` as a linear combination
            of simple roots.

            OUTPUT:

            A list with entries in `\QQ` representing ``self`` as a linear
            combination of simple roots.

            EXAMPLES::

                sage: L = RootSystem(['A', 3]).ambient_space()
                sage: e = L.basis()
                sage: (e[0] + 3*e[3])._to_root_vector()  # not in the root space
                sage: (e[0] - e[1])._to_root_vector()
                (1, 0, 0)
                sage: (e[0] + 2*e[1] - 3*e[2])._to_root_vector()
                (1, 3, 0)
            """
            v = self.to_vector().change_ring(QQ)
            al = [a.to_vector() for a in self.parent().simple_roots()]
            b = v.parent().linear_dependence([v] + al)
            if len(b) != 1 or b[0] == 0:
                return None
            b = b[0]  # Get the actual vector that gives the linear dependency
            return b[1:].change_ring(QQ) / -b[0]
