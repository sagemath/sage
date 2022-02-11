r"""
Finite posets

Here is some terminology used in this file:

- An *order filter* (or *upper set*) of a poset `P` is a subset `S` of `P`
  such that if `x \leq y` and `x\in S` then `y\in S`.

- An *order ideal* (or *lower set*) of a poset `P` is a subset `S` of `P`
  such that if `x \leq y` and `y\in S` then `x\in S`.
"""
#*****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.categories.category_with_axiom import CategoryWithAxiom

class FinitePosets(CategoryWithAxiom):
    r"""
    The category of finite posets i.e. finite sets with a partial
    order structure.

    EXAMPLES::

        sage: FinitePosets()
        Category of finite posets
        sage: FinitePosets().super_categories()
        [Category of posets, Category of finite sets]
        sage: FinitePosets().example()
        NotImplemented

    .. SEEALSO:: :class:`~sage.categories.posets.Posets`, :func:`Poset`

    TESTS::

        sage: C = FinitePosets()
        sage: C is Posets().Finite()
        True
        sage: TestSuite(C).run()

    """

    class ParentMethods:

        ##########################################################################
        # Properties of this poset

        def is_lattice(self):
            r"""
            Return whether the poset is a lattice.

            A poset is a lattice if all pairs of elements have
            both a least upper bound ("join") and a greatest lower bound
            ("meet") in the poset.

            EXAMPLES::

                sage: P = Poset([[1, 3, 2], [4], [4, 5, 6], [6], [7], [7], [7], []])
                sage: P.is_lattice()
                True

                sage: P = Poset([[1, 2], [3], [3], []])
                sage: P.is_lattice()
                True

                sage: P = Poset({0: [2, 3], 1: [2, 3]})
                sage: P.is_lattice()
                False

                sage: P = Poset({1: [2, 3, 4], 2: [5, 6], 3: [5, 7], 4: [6, 7], 5: [8, 9],
                ....:            6: [8, 10], 7: [9, 10], 8: [11], 9: [11], 10: [11]})
                sage: P.is_lattice()
                False

            TESTS::

                sage: P = Poset()
                sage: P.is_lattice()
                True

            .. SEEALSO::

                - Weaker properties: :meth:`~sage.combinat.posets.posets.FinitePoset.is_join_semilattice`,
                  :meth:`~sage.combinat.posets.posets.FinitePoset.is_meet_semilattice`
            """
            return (self.cardinality() == 0 or
                     (self.has_bottom() and self.is_join_semilattice()))

        def is_self_dual(self):
            r"""
            Return whether the poset is *self-dual*.

            A poset is self-dual if it is isomorphic to its dual poset.

            EXAMPLES::

                sage: P = Poset({1: [3, 4], 2: [3, 4]})
                sage: P.is_self_dual()
                True

                sage: P = Poset({1: [2, 3]})
                sage: P.is_self_dual()
                False

            TESTS::

                sage: P = Poset()
                sage: P.is_self_dual()
                True

            .. SEEALSO::

                - Stronger properties: :meth:`~sage.combinat.posets.lattices.FiniteLatticePoset.is_orthocomplemented` (for lattices)
                - Other: :meth:`~sage.combinat.posets.posets.FinitePoset.dual`
            """
            # Two quick checks before full isomorphic test.
            if sorted(self._hasse_diagram.in_degree()) != sorted(self._hasse_diagram.out_degree()):
                return False
            levels_orig = [len(x) for x in self._hasse_diagram.level_sets()]
            dual_poset_hasse = self._hasse_diagram.reverse()
            levels_dual = [len(x) for x in dual_poset_hasse.level_sets()]
            if levels_orig != levels_dual:
                return False
            return self._hasse_diagram.is_isomorphic(dual_poset_hasse)

        ##########################################################################
        # Properties of morphisms

        def is_poset_isomorphism(self, f, codomain):
            r"""
            Return whether `f` is an isomorphism of posets from
            ``self`` to ``codomain``.

            INPUT:

            - ``f`` -- a function from ``self`` to ``codomain``
            - ``codomain`` -- a poset

            EXAMPLES:

            We build the poset `D` of divisors of 30, and check that
            it is isomorphic to the boolean lattice `B` of the subsets
            of `\{2,3,5\}` ordered by inclusion, via the reverse
            function `f: B \to D, b \mapsto \prod_{x\in b} x`::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: B = Poset(([frozenset(s) for s in Subsets([2,3,5])], attrcall("issubset")))
                sage: def f(b): return D(prod(b))
                sage: B.is_poset_isomorphism(f, D)
                True

            On the other hand, `f` is not an isomorphism to the chain
            of divisors of 30, ordered by usual comparison::

                sage: P = Poset((divisors(30), operator.le))
                sage: def f(b): return P(prod(b))
                sage: B.is_poset_isomorphism(f, P)
                False

            A non surjective case::

                sage: B = Poset(([frozenset(s) for s in Subsets([2,3])], attrcall("issubset")))
                sage: def f(b): return D(prod(b))
                sage: B.is_poset_isomorphism(f, D)
                False

            A non injective case::

                sage: B = Poset(([frozenset(s) for s in Subsets([2,3,5,6])], attrcall("issubset")))
                sage: def f(b): return D(gcd(prod(b), 30))
                sage: B.is_poset_isomorphism(f, D)
                False

            .. note:: since ``D`` and ``B`` are not facade posets, ``f`` is
               responsible for the conversions between integers and subsets to
               elements of ``D`` and ``B`` and back.

            .. SEEALSO:: :meth:`FiniteLatticePosets.ParentMethods.is_lattice_morphism`
            """
            image = set(f(x) for x in self)
            if len(image) != self.cardinality():
                # Not injective
                return False
            if len(image) != codomain.cardinality():
                # Not surjective
                return False
            for x in self:
                if set(f(y) for y in self.upper_covers(x)) != set(codomain.upper_covers(f(x))):
                    return False
            return True

        def is_poset_morphism(self, f, codomain):
            r"""
            Return whether `f` is a morphism of posets from ``self``
            to ``codomain``, that is

            .. MATH::

                x\leq y \Longrightarrow f(x) \leq f(y)

            for all `x` and `y` in ``self``.

            INPUT:

            - ``f`` -- a function from ``self`` to ``codomain``
            - ``codomain`` -- a poset

            EXAMPLES:

            We build the boolean lattice of the subsets of
            `\{2,3,5,6\}` and the lattice of divisors of `30`, and
            check that the map `b \mapsto \gcd(\prod_{x\in b} x, 30)`
            is a morphism of posets::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: B = Poset(([frozenset(s) for s in Subsets([2,3,5,6])], attrcall("issubset")))
                sage: def f(b): return D(gcd(prod(b), 30))
                sage: B.is_poset_morphism(f, D)
                True

            .. note:: since ``D`` and ``B`` are not facade posets, ``f`` is responsible
               for the conversions between integers and subsets to elements of
               ``D`` and ``B`` and back.

            `f` is also a morphism of posets to the chain of divisors
            of 30, ordered by usual comparison::

                sage: P = Poset((divisors(30), operator.le))
                sage: def f(b): return P(gcd(prod(b), 30))
                sage: B.is_poset_morphism(f, P)
                True

            FIXME: should this be ``is_order_preserving_morphism``?

            .. SEEALSO:: :meth:`is_poset_isomorphism`

            TESTS:

            Base cases::

                sage: P = posets.ChainPoset(2)
                sage: Q = posets.AntichainPoset(2)
                sage: f = lambda x: 1-x
                sage: P.is_poset_morphism(f, P)
                False
                sage: P.is_poset_morphism(f, Q)
                False
                sage: Q.is_poset_morphism(f, Q)
                True
                sage: Q.is_poset_morphism(f, P)
                True

                sage: P = Poset(); P
                Finite poset containing 0 elements
                sage: P.is_poset_morphism(f, P)
                True

            """
            for x in self:
                for y in self.upper_covers(x):
                    if not codomain.is_lequal(f(x),f(y)):
                        return False
            return True

        ##########################################################################
        # About order ideals, order filters and the like

        def order_ideal_generators(self, ideal, direction='down'):
            r"""
            Return the antichain of (minimal) generators of the order
            ideal (resp. order filter) ``ideal``.

            INPUT:

            - ``ideal`` -- an order ideal `I` (resp. order filter)
              of ``self``, as a list (or iterable); this should be
              an order ideal if ``direction`` is set to ``'down'``,
              and an order filter if ``direction`` is set to
              ``'up'``.
            - ``direction`` -- ``'up'`` or ``'down'`` (default:
              ``'down'``).

            The antichain of (minimal) generators of an order ideal
            `I` in a poset `P` is the set of all minimal elements of
            `P`. In the case of an order filter, the definition is
            similar, but with "maximal" used instead of "minimal".

            EXAMPLES:

            We build the boolean lattice of all subsets of `\{1,2,3\}`
            ordered by inclusion, and compute an order ideal there::

                sage: P = Poset((Subsets([1,2,3]), attrcall("issubset")))
                sage: I = P.order_ideal([Set([1,2]), Set([2,3]), Set([1])])
                sage: sorted(sorted(p) for p in I)
                [[], [1], [1, 2], [2], [2, 3], [3]]

            Then, we retrieve the generators of this ideal::

                sage: gen = P.order_ideal_generators(I)
                sage: sorted(sorted(p) for p in gen)
                [[1, 2], [2, 3]]

            If ``direction`` is 'up', then this instead computes
            the minimal generators for an order filter::

                sage: I = P.order_filter([Set([1,2]), Set([2,3]), Set([1])])
                sage: sorted(sorted(p) for p in I)
                [[1], [1, 2], [1, 2, 3], [1, 3], [2, 3]]
                sage: gen = P.order_ideal_generators(I, direction='up')
                sage: sorted(sorted(p) for p in gen)
                [[1], [2, 3]]

            Complexity: `O(n+m)` where `n` is the cardinality of `I`,
            and `m` the number of upper covers of elements of `I`.
            """
            if direction == 'down':
                covers = self.upper_covers
            else:
                covers = self.lower_covers
            ideal_as_set = set(ideal)
            from sage.sets.set import Set
            return Set(x for x in ideal if all(y not in ideal_as_set
                                               for y in covers(x)))

        def order_filter_generators(self, filter):
            r"""
            Generators for an order filter

            INPUT:

            - ``filter`` -- an order filter of ``self``, as a list (or iterable)

            EXAMPLES::

                sage: P = Poset((Subsets([1,2,3]), attrcall("issubset")))
                sage: I = P.order_filter([Set([1,2]), Set([2,3]), Set([1])])
                sage: sorted(sorted(p) for p in I)
                [[1], [1, 2], [1, 2, 3], [1, 3], [2, 3]]
                sage: gen = P.order_filter_generators(I)
                sage: sorted(sorted(p) for p in gen)
                [[1], [2, 3]]

            .. SEEALSO:: :meth:`order_ideal_generators`
            """
            return self.order_ideal_generators(filter, direction='up')

        def order_ideal_complement_generators(self, antichain, direction='up'):
            r"""
            Return the Panyushev complement of the antichain
            ``antichain``.

            Given an antichain `A` of a poset `P`, the Panyushev
            complement of `A` is defined to be the antichain consisting
            of the minimal elements of the order filter `B`, where `B`
            is the (set-theoretic) complement of the order ideal of
            `P` generated by `A`.

            Setting the optional keyword variable ``direction`` to
            ``'down'`` leads to the inverse Panyushev complement being
            computed instead of the Panyushev complement. The inverse
            Panyushev complement of an antichain `A` is the antichain
            whose Panyushev complement is `A`. It can be found as the
            antichain consisting of the maximal elements of the order
            ideal `C`, where `C` is the (set-theoretic) complement of
            the order filter of `P` generated by `A`.

            :meth:`panyushev_complement` is an alias for this method.

            Panyushev complementation is related (actually, isomorphic)
            to rowmotion (:meth:`rowmotion`).

            INPUT:

            - ``antichain`` -- an antichain of ``self``, as a list (or
              iterable), or, more generally, generators of an order ideal
              (resp. order filter)
            - ``direction`` -- 'up' or 'down' (default: 'up')

            OUTPUT:

            - the generating antichain of the complement order filter
              (resp. order ideal) of the order ideal (resp. order filter)
              generated by the antichain ``antichain``

            EXAMPLES::

                sage: P = Poset( ( [1,2,3], [ [1,3], [2,3] ] ) )
                sage: P.order_ideal_complement_generators([1])
                {2}
                sage: P.order_ideal_complement_generators([3])
                set()
                sage: P.order_ideal_complement_generators([1,2])
                {3}
                sage: P.order_ideal_complement_generators([1,2,3])
                set()

                sage: P.order_ideal_complement_generators([1], direction="down")
                {2}
                sage: P.order_ideal_complement_generators([3], direction="down")
                {1, 2}
                sage: P.order_ideal_complement_generators([1,2], direction="down")
                set()
                sage: P.order_ideal_complement_generators([1,2,3], direction="down")
                set()

            .. WARNING::

                This is a brute force implementation, building the
                order ideal generated by the antichain, and searching
                for order filter generators of its complement
            """
            if direction == 'up':
                I = self.order_ideal(antichain)
            else:
                I = self.order_filter(antichain)
            I_comp = set(self).difference(I)
            return set(self.order_ideal_generators(I_comp, direction = direction))

        panyushev_complement = order_ideal_complement_generators

        def rowmotion(self, order_ideal):
            r"""
            The image of the order ideal ``order_ideal`` under rowmotion
            in ``self``.

            Rowmotion on a finite poset `P` is an automorphism of the set
            `J(P)` of all order ideals of `P`. One way to define it is as
            follows: Given an order ideal `I \in J(P)`, we let `F` be the
            set-theoretic complement of `I` in `P`. Furthermore we let
            `A` be the antichain consisting of all minimal elements of
            `F`. Then, the rowmotion of `I` is defined to be the order
            ideal of `P` generated by the antichain `A` (that is, the
            order ideal consisting of each element of `P` which has some
            element of `A` above it).

            Rowmotion is related (actually, isomorphic) to Panyushev
            complementation (:meth:`panyushev_complement`).

            INPUT:

            - ``order_ideal`` -- an order ideal of ``self``, as a set

            OUTPUT:

            - the image of ``order_ideal`` under rowmotion, as a set again

            EXAMPLES::

                sage: P = Poset( {1: [2, 3], 2: [], 3: [], 4: [8], 5: [], 6: [5], 7: [1, 4], 8: []} )
                sage: I = Set({2, 6, 1, 7})
                sage: P.rowmotion(I)
                {1, 3, 4, 5, 6, 7}

                sage: P = Poset( {} )
                sage: I = Set({})
                sage: P.rowmotion(I)
                {}
            """
            result = order_ideal
            for i in reversed(self.linear_extension()):
                result = self.order_ideal_toggle(result, i)
            return result

        def birational_free_labelling(self, linear_extension=None,
                                      prefix='x', base_field=None,
                                      reduced=False, addvars=None,
                                      labels=None,
                                      min_label=None,
                                      max_label=None):
            r"""
            Return the birational free labelling of ``self``.

            Let us hold back defining this, and introduce birational
            toggles and birational rowmotion first. These notions have
            been introduced in [EP2013]_ as generalizations of the notions
            of toggles (:meth:`~sage.categories.posets.Posets.ParentMethods.order_ideal_toggle`)
            and :meth:`rowmotion <rowmotion>` on order ideals of a finite poset. They
            have been studied further in [GR2013]_.

            Let `\mathbf{K}` be a field, and `P` be a finite poset. Let
            `\widehat{P}` denote the poset obtained from `P` by adding a
            new element `1` which is greater than all existing elements
            of `P`, and a new element `0` which is smaller than all
            existing elements of `P` and `1`. Now, a `\mathbf{K}`-*labelling
            of* `P` will mean any function from `\widehat{P}` to `\mathbf{K}`.
            The image of an element `v` of `\widehat{P}` under this labelling
            will be called the *label* of this labelling at `v`. The set
            of all `\mathbf{K}`-labellings of `P` is clearly
            `\mathbf{K}^{\widehat{P}}`.

            For any `v \in P`, we now define a rational map
            `T_v : \mathbf{K}^{\widehat{P}} \dashrightarrow
            \mathbf{K}^{\widehat{P}}` as follows: For every `f \in
            \mathbf{K}^{\widehat{P}}`, the image `T_v f` should send every
            element `u \in \widehat{P}` distinct from `v` to `f(u)` (so the
            labels at all `u \neq v` don't change), while `v` is sent to

            .. MATH::

                \frac{1}{f(v)} \cdot
                \frac{\sum_{u \lessdot v} f(u)}
                {\sum_{u \gtrdot v} \frac{1}{f(u)}}

            (both sums are over all `u \in \widehat{P}` satisfying the
            respectively given conditions). Here, `\lessdot` and `\gtrdot`
            mean (respectively) "covered by" and "covers", interpreted with
            respect to the poset `\widehat{P}`. This rational map `T_v`
            is an involution and is called the *(birational)* `v`-*toggle*; see
            :meth:`birational_toggle` for its implementation.

            Now, *birational rowmotion* is defined as the composition
            `T_{v_1} \circ T_{v_2} \circ \cdots \circ T_{v_n}`, where
            `(v_1, v_2, \ldots, v_n)` is a linear extension of `P`
            (written as a linear ordering of the elements of `P`). This
            is a rational map
            `\mathbf{K}^{\widehat{P}} \dashrightarrow \mathbf{K}^{\widehat{P}}`
            which does not depend on the choice of the linear extension;
            it is denoted by `R`. See :meth:`birational_rowmotion` for
            its implementation.

            The definitions of birational toggles and birational
            rowmotion extend to the case of `\mathbf{K}` being any semifield
            rather than necessarily a field (although it becomes less
            clear what constitutes a rational map in this generality).
            The most useful case is that of the :class:`tropical semiring
            <sage.rings.semirings.tropical_semiring.TropicalSemiring>`,
            in which case birational rowmotion relates to classical
            constructions such as promotion of rectangular semistandard
            Young tableaux (page 5 of [EP2013b]_ and future work, via the
            related notion of birational *promotion*) and rowmotion on
            order ideals of the poset ([EP2013]_).

            The *birational free labelling* is a special labelling
            defined for every finite poset `P` and every linear extension
            `(v_1, v_2, \ldots, v_n)` of `P`. It is given by sending
            every element `v_i` in `P` to `x_i`, sending the element `0`
            of `\widehat{P}` to `a`, and sending the element `1` of
            `\widehat{P}` to `b`, where the ground field `\mathbf{K}` is the
            field of rational functions in `n+2` indeterminates
            `a, x_1, x_2, \ldots, x_n, b` over `\mathbb Q`.

            In Sage, a labelling `f` of a poset `P` is encoded as a
            `4`-tuple `(\mathbf{K}, d, u, v)`, where `\mathbf{K}` is the
            ground field of the labelling (i. e., its target), `d` is the
            dictionary containing the values of `f` at the elements of
            `P` (the keys being the respective elements of `P`), `u`
            is the label of `f` at `0`, and `v` is the label of `f` at
            `1`.

            .. WARNING::

                The dictionary `d` is labelled by the elements of `P`.
                If `P` is a poset with ``facade`` option set to
                ``False``, these might not be what they seem to be!
                (For instance, if
                ``P == Poset({1: [2, 3]}, facade=False)``, then the
                value of `d` at `1` has to be accessed by ``d[P(1)]``, not
                by ``d[1]``.)

            .. WARNING::

                Dictionaries are mutable. They do compare correctly,
                but are not hashable and need to be cloned to avoid
                spooky action at a distance. Be careful!

            INPUT:

            - ``linear_extension`` -- (default: the default linear
              extension of ``self``) a linear extension of ``self``
              (as a linear extension or as a list), or more generally
              a list of all elements of all elements of ``self`` each
              occurring exactly once

            - ``prefix`` -- (default: ``'x'``) the prefix to name
              the indeterminates corresponding to the elements of
              ``self`` in the labelling (so, setting it to
              ``'frog'`` will result in these indeterminates being
              called ``frog1, frog2, ..., frogn`` rather than
              ``x1, x2, ..., xn``).

            - ``base_field`` -- (default: ``QQ``) the base field to
              be used instead of `\QQ` to define the rational
              function field over; this is not going to be the base
              field of the labelling, because the latter will have
              indeterminates adjoined!

            - ``reduced`` -- (default: ``False``) if set to
              ``True``, the result will be the *reduced* birational
              free labelling, which differs from the regular one by
              having `0` and `1` both sent to `1` instead of `a` and
              `b` (the indeterminates `a` and `b` then also won't
              appear in the ground field)

            - ``addvars`` -- (default: ``''``) a string containing
              names of extra variables to be adjoined to the ground
              field (these don't have an effect on the labels)

            - ``labels`` -- (default: ``'x'``) Either a function
              that takes an element of the poset and returns a name
              for the indeterminate corresponding to that element,
              or a string containing a comma-separated list of
              indeterminates that will be assigned to elements in
              the order of ``linear_extension``. If the
              list contains more indeterminates than needed, the
              excess will be ignored. If it contains too few, then
              the needed indeterminates will be constructed from
              ``prefix``.

            - ``min_label`` -- (default: ``'a'``) a string to be
              used as the label for the element `0` of `\widehat{P}`

            - ``max_label`` -- (default: ``'b'``) a string to be
              used as the label for the element `1` of `\widehat{P}`

            OUTPUT:

            The birational free labelling of the poset ``self`` and the
            linear extension ``linear_extension``. Or, if ``reduced``
            is set to ``True``, the reduced birational free labelling.

            EXAMPLES:

            We construct the birational free labelling on a simple
            poset::

                sage: P = Poset({1: [2, 3]})
                sage: l = P.birational_free_labelling(); l
                (Fraction Field of Multivariate Polynomial Ring in a, x1, x2, x3, b over Rational Field,
                 {...},
                 a,
                 b)
                sage: sorted(l[1].items())
                [(1, x1), (2, x2), (3, x3)]

                sage: l = P.birational_free_labelling(linear_extension=[1, 3, 2]); l
                (Fraction Field of Multivariate Polynomial Ring in a, x1, x2, x3, b over Rational Field,
                 {...},
                 a,
                 b)
                sage: sorted(l[1].items())
                [(1, x1), (2, x3), (3, x2)]

                sage: l = P.birational_free_labelling(linear_extension=[1, 3, 2], reduced=True, addvars="spam, eggs"); l
                (Fraction Field of Multivariate Polynomial Ring in x1, x2, x3, spam, eggs over Rational Field,
                 {...},
                 1,
                 1)
                sage: sorted(l[1].items())
                [(1, x1), (2, x3), (3, x2)]

                sage: l = P.birational_free_labelling(linear_extension=[1, 3, 2], prefix="wut", reduced=True, addvars="spam, eggs"); l
                (Fraction Field of Multivariate Polynomial Ring in wut1, wut2, wut3, spam, eggs over Rational Field,
                 {...},
                 1,
                 1)
                sage: sorted(l[1].items())
                [(1, wut1), (2, wut3), (3, wut2)]

                sage: l = P.birational_free_labelling(linear_extension=[1, 3, 2], reduced=False, addvars="spam, eggs"); l
                (Fraction Field of Multivariate Polynomial Ring in a, x1, x2, x3, b, spam, eggs over Rational Field,
                 {...},
                 a,
                 b)
                sage: sorted(l[1].items())
                [(1, x1), (2, x3), (3, x2)]
                sage: l[1][2]
                x3

            Illustrating labelling with a function::

                sage: P = posets.ChainPoset(2).product(posets.ChainPoset(2))
                sage: l = P.birational_free_labelling(labels=lambda e : 'x_' + str(e[0]) + str(e[1]))
                sage: sorted(l[1].items())
                [((0, 0), x_00), ((0, 1), x_01), ((1, 0), x_10), ((1, 1), x_11)]
                sage: l[2]
                a

            The same, but with ``min_label`` and ``max_label`` provided::

                sage: P = posets.ChainPoset(2).product(posets.ChainPoset(2))
                sage: l = P.birational_free_labelling(labels=lambda e : 'x_' + str(e[0]) + str(e[1]), min_label="lambda", max_label="mu")
                sage: sorted(l[1].items())
                [((0, 0), x_00), ((0, 1), x_01), ((1, 0), x_10), ((1, 1), x_11)]
                sage: l[2]
                lambda
                sage: l[3]
                mu

            Illustrating labelling with a comma separated list of labels::

                sage: l = P.birational_free_labelling(labels='w,x,y,z')
                sage: sorted(l[1].items())
                [((0, 0), w), ((0, 1), x), ((1, 0), y), ((1, 1), z)]
                sage: l = P.birational_free_labelling(labels='w,x,y,z,m')
                sage: sorted(l[1].items())
                [((0, 0), w), ((0, 1), x), ((1, 0), y), ((1, 1), z)]
                sage: l = P.birational_free_labelling(labels='w')
                sage: sorted(l[1].items())
                [((0, 0), w), ((0, 1), x1), ((1, 0), x2), ((1, 1), x3)]

            Illustrating the warning about facade::

                sage: P = Poset({1: [2, 3]}, facade=False)
                sage: l = P.birational_free_labelling(linear_extension=[1, 3, 2], reduced=False, addvars="spam, eggs"); l
                (Fraction Field of Multivariate Polynomial Ring in a, x1, x2, x3, b, spam, eggs over Rational Field,
                 {...},
                 a,
                 b)
                sage: l[1][2]
                Traceback (most recent call last):
                ...
                KeyError: 2
                sage: l[1][P(2)]
                x3

            Another poset::

                sage: P = posets.SSTPoset([2,1])
                sage: lext = sorted(P)
                sage: l = P.birational_free_labelling(linear_extension=lext, addvars="ohai")
                sage: l
                (Fraction Field of Multivariate Polynomial Ring in a, x1, x2, x3, x4, x5, x6, x7, x8, b, ohai over Rational Field,
                 {...},
                 a,
                 b)
                sage: sorted(l[1].items())
                [([[1, 1], [2]], x1), ([[1, 1], [3]], x2), ([[1, 2], [2]], x3), ([[1, 2], [3]], x4),
                 ([[1, 3], [2]], x5), ([[1, 3], [3]], x6), ([[2, 2], [3]], x7), ([[2, 3], [3]], x8)]

            See :meth:`birational_rowmotion`, :meth:`birational_toggle` and
            :meth:`birational_toggles` for more substantial examples of what
            one can do with the birational free labelling.

            TESTS:

            The ``linear_extension`` keyword does not have to be given an
            actual linear extension::

                sage: P = posets.ChainPoset(2).product(posets.ChainPoset(3))
                sage: P
                Finite lattice containing 6 elements
                sage: lex = [(1,0),(0,0),(1,1),(0,1),(1,2),(0,2)]
                sage: l = P.birational_free_labelling(linear_extension=lex,
                ....:                                 prefix="u", reduced=True)
                sage: l
                (Fraction Field of Multivariate Polynomial Ring in u1, u2, u3, u4, u5, u6 over Rational Field,
                 {...},
                 1,
                 1)
                sage: sorted(l[1].items())
                [((0, 0), u2),
                 ((0, 1), u4),
                 ((0, 2), u6),
                 ((1, 0), u1),
                 ((1, 1), u3),
                 ((1, 2), u5)]

            For comparison, the standard linear extension::

                sage: l = P.birational_free_labelling(prefix="u", reduced=True); l
                (Fraction Field of Multivariate Polynomial Ring in u1, u2, u3, u4, u5, u6 over Rational Field,
                 {...},
                 1,
                 1)
                sage: sorted(l[1].items())
                [((0, 0), u1),
                 ((0, 1), u2),
                 ((0, 2), u3),
                 ((1, 0), u4),
                 ((1, 1), u5),
                 ((1, 2), u6)]

            If you want your linear extension to be tested for being a
            linear extension, just call the ``linear_extension`` method
            on the poset::

                sage: lex = [(0,0),(0,1),(1,0),(1,1),(0,2),(1,2)]
                sage: l = P.birational_free_labelling(linear_extension=P.linear_extension(lex),
                ....:                                 prefix="u", reduced=True)
                sage: l
                (Fraction Field of Multivariate Polynomial Ring in u1, u2, u3, u4, u5, u6 over Rational Field,
                 {...},
                 1,
                 1)
                sage: sorted(l[1].items())
                [((0, 0), u1),
                 ((0, 1), u2),
                 ((0, 2), u5),
                 ((1, 0), u3),
                 ((1, 1), u4),
                 ((1, 2), u6)]

            Nonstandard base field::

                sage: P = Poset({1: [3], 2: [3,4]})
                sage: lex = [1, 2, 4, 3]
                sage: l = P.birational_free_labelling(linear_extension=lex,
                ....:                                 prefix="aaa",
                ....:                                 base_field=Zmod(13))
                sage: l
                (Fraction Field of Multivariate Polynomial Ring in a, aaa1, aaa2, aaa3, aaa4, b over Ring of integers modulo 13,
                 {...},
                 a,
                 b)
                sage: l[1][4]
                aaa3

            The empty poset::

                sage: P = Poset({})
                sage: P.birational_free_labelling(reduced=False, addvars="spam, eggs")
                (Fraction Field of Multivariate Polynomial Ring in a, b, spam, eggs over Rational Field,
                 {},
                 a,
                 b)
                sage: P.birational_free_labelling(reduced=True, addvars="spam, eggs")
                (Fraction Field of Multivariate Polynomial Ring in spam, eggs over Rational Field,
                 {},
                 1,
                 1)
                sage: P.birational_free_labelling(reduced=True)
                (Multivariate Polynomial Ring in no variables over Rational Field,
                 {},
                 1,
                 1)
                sage: P.birational_free_labelling(prefix="zzz")
                (Fraction Field of Multivariate Polynomial Ring in a, b over Rational Field,
                 {},
                 a,
                 b)
                sage: P.birational_free_labelling(labels="x,y,z", min_label="spam", max_label="eggs")
                (Fraction Field of Multivariate Polynomial Ring in spam, eggs over Rational Field,
                 {},
                 spam,
                 eggs)
            """
            if base_field is None:
                from sage.rings.rational_field import QQ
                base_field = QQ
            if linear_extension is None:
                linear_extension = self.linear_extension()
            n = self.cardinality()
            label_list = []
            if labels:
                if callable(labels):
                    label_list = [labels(e) for e in linear_extension]
                else:
                    label_list = labels.split(',')
                    if len(label_list) > n:
                        label_list = label_list[:n]
                    elif len(label_list) < n:
                        label_list += [prefix + str(i) for i in range(1, n + 1 - len(label_list))]
            else:
                label_list = [prefix + str(i) for i in range(1, n + 1)]
            if not reduced:
                if min_label is None:
                    min_label = 'a'
                if max_label is None:
                    max_label = 'b'
                label_list = [min_label] + label_list + [max_label]
            if addvars:
                label_list += addvars.split(',')
            varstring = ','.join(label_list)
            varnum = len(label_list)

            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            PR = PolynomialRing(base_field, varstring, varnum)
            # Now, ``PR`` is the polynomial ring in `n + 2` indeterminates
            # (or more, if ``addvars`` was set; or less, if ``reduced`` is
            # ``True``) over ``base_field``.
            # The first `n + 2` of these indeterminates are named
            # ``a, x1, x2, ..., xn, b`` (if ``reduced`` is ``False``).
            # These will label the vertices of `\widehat{P}`.
            if reduced:
                xs = tuple(PR.gens()[: n])
            else:
                xs = tuple(PR.gens()[1 : n + 1])
            # So ``xs`` is the list ``[x1, x2, ..., xn]``.
            if not reduced:
                a = PR.gens()[0]
                b = PR.gens()[n + 1]
            else:
                a = PR.one()
                b = PR.one()
            # So ``a`` and ``b`` are the labels at `0` and `1`.
            FF = PR.fraction_field()
            # ``FF`` is the field of rational functions.
            dct = {self(p): xs[i] for (i, p) in enumerate(linear_extension)}
            return (FF, dct, a, b)

        def birational_toggle(self, v, labelling):
            r"""
            Return the result of applying the birational `v`-toggle `T_v`
            to the `\mathbf{K}`-labelling ``labelling`` of the poset ``self``.

            See the documentation of :meth:`birational_free_labelling`
            for a definition of this toggle and of `\mathbf{K}`-labellings as
            well as an explanation of how `\mathbf{K}`-labellings are to be
            encoded to be understood by Sage. This implementation allows
            `\mathbf{K}` to be a semifield, not just a field. The birational
            `v`-toggle is only a rational map, so an exception (most
            likely, ``ZeroDivisionError``) will be thrown if the
            denominator is zero.

            INPUT:

            - ``v`` -- an element of ``self`` (must have ``self`` as
              parent if ``self`` is a ``facade=False`` poset)

            - ``labelling`` -- a `\mathbf{K}`-labelling of ``self`` in the
              sense as defined in the documentation of
              :meth:`birational_free_labelling`

            OUTPUT:

            The `\mathbf{K}`-labelling `T_v f` of ``self``, where `f` is
            ``labelling``.

            EXAMPLES:

            Let us start with the birational free labelling of the
            "V"-poset (the three-element poset with Hasse diagram looking
            like a "V")::

                sage: V = Poset({1: [2, 3]})
                sage: s = V.birational_free_labelling(); s
                (Fraction Field of Multivariate Polynomial Ring in a, x1, x2, x3, b over Rational Field,
                 {...},
                 a,
                 b)
                sage: sorted(s[1].items())
                [(1, x1), (2, x2), (3, x3)]

            The image of `s` under the `1`-toggle `T_1` is::

                sage: s1 = V.birational_toggle(1, s); s1
                (Fraction Field of Multivariate Polynomial Ring in a, x1, x2, x3, b over Rational Field,
                 {...},
                 a,
                 b)
                sage: sorted(s1[1].items())
                [(1, a*x2*x3/(x1*x2 + x1*x3)), (2, x2), (3, x3)]

            Now let us apply the `2`-toggle `T_2` (to the old ``s``)::

                sage: s2 = V.birational_toggle(2, s); s2
                (Fraction Field of Multivariate Polynomial Ring in a, x1, x2, x3, b over Rational Field,
                 {...},
                 a,
                 b)
                sage: sorted(s2[1].items())
                [(1, x1), (2, x1*b/x2), (3, x3)]

            On the other hand, we can also apply `T_2` to the image of `s`
            under `T_1`::

                sage: s12 = V.birational_toggle(2, s1); s12
                (Fraction Field of Multivariate Polynomial Ring in a, x1, x2, x3, b over Rational Field,
                 {...},
                 a,
                 b)
                sage: sorted(s12[1].items())
                [(1, a*x2*x3/(x1*x2 + x1*x3)), (2, a*x3*b/(x1*x2 + x1*x3)), (3, x3)]

            Each toggle is an involution::

                sage: all( V.birational_toggle(i, V.birational_toggle(i, s)) == s
                ....:      for i in V )
                True

            We can also start with a less generic labelling::

                sage: t = (QQ, {1: 3, 2: 6, 3: 7}, 2, 10)
                sage: t1 = V.birational_toggle(1, t); t1
                (Rational Field, {...}, 2, 10)
                sage: sorted(t1[1].items())
                [(1, 28/13), (2, 6), (3, 7)]
                sage: t13 = V.birational_toggle(3, t1); t13
                (Rational Field, {...}, 2, 10)
                sage: sorted(t13[1].items())
                [(1, 28/13), (2, 6), (3, 40/13)]

            However, labellings have to be sufficiently generic, lest
            denominators vanish::

                sage: t = (QQ, {1: 3, 2: 5, 3: -5}, 1, 15)
                sage: t1 = V.birational_toggle(1, t)
                Traceback (most recent call last):
                ...
                ZeroDivisionError: rational division by zero

            We don't get into zero-division issues in the tropical
            semiring (unless the zero of the tropical semiring appears
            in the labelling)::

                sage: TT = TropicalSemiring(QQ)
                sage: t = (TT, {1: TT(2), 2: TT(4), 3: TT(1)}, TT(6), TT(0))
                sage: t1 = V.birational_toggle(1, t); t1
                (Tropical semiring over Rational Field, {...}, 6, 0)
                sage: sorted(t1[1].items())
                [(1, 8), (2, 4), (3, 1)]
                sage: t12 = V.birational_toggle(2, t1); t12
                (Tropical semiring over Rational Field, {...}, 6, 0)
                sage: sorted(t12[1].items())
                [(1, 8), (2, 4), (3, 1)]
                sage: t123 = V.birational_toggle(3, t12); t123
                (Tropical semiring over Rational Field, {...}, 6, 0)
                sage: sorted(t123[1].items())
                [(1, 8), (2, 4), (3, 7)]

            We turn to more interesting posets. Here is the `6`-element
            poset arising from the weak order on `S_3`::

                sage: P = posets.SymmetricGroupWeakOrderPoset(3)
                sage: sorted(list(P))
                ['123', '132', '213', '231', '312', '321']
                sage: t = (TT, {'123': TT(4), '132': TT(2), '213': TT(3), '231': TT(1), '321': TT(1), '312': TT(2)}, TT(7), TT(1))
                sage: t1 = P.birational_toggle('123', t); t1
                (Tropical semiring over Rational Field, {...}, 7, 1)
                sage: sorted(t1[1].items())
                [('123', 6), ('132', 2), ('213', 3), ('231', 1), ('312', 2), ('321', 1)]
                sage: t13 = P.birational_toggle('213', t1); t13
                (Tropical semiring over Rational Field, {...}, 7, 1)
                sage: sorted(t13[1].items())
                [('123', 6), ('132', 2), ('213', 4), ('231', 1), ('312', 2), ('321', 1)]

            Let us verify on this example some basic properties of
            toggles. First of all, again let us check that `T_v` is an
            involution for every `v`::

                sage: all( P.birational_toggle(v, P.birational_toggle(v, t)) == t
                ....:      for v in P )
                True

            Furthermore, two toggles `T_v` and `T_w` commute unless
            one of `v` or `w` covers the other::

                sage: all( P.covers(v, w) or P.covers(w, v)
                ....:      or P.birational_toggle(v, P.birational_toggle(w, t))
                ....:         == P.birational_toggle(w, P.birational_toggle(v, t))
                ....:      for v in P for w in P )
                True

            TESTS:

            Setting ``facade`` to ``False`` does not break
            ``birational_toggle``::

                sage: P = Poset({'x': ['y', 'w'], 'y': ['z'], 'w': ['z']}, facade=False)
                sage: lex = ['x', 'y', 'w', 'z']
                sage: t = P.birational_free_labelling(linear_extension=lex)
                sage: all( P.birational_toggle(v, P.birational_toggle(v, t)) == t
                ....:      for v in P )
                True
                sage: t4 = P.birational_toggle(P('z'), t); t4
                (Fraction Field of Multivariate Polynomial Ring in a, x1, x2, x3, x4, b over Rational Field,
                 {...},
                 a,
                 b)
                sage: t4[1][P('x')]
                x1
                sage: t4[1][P('y')]
                x2
                sage: t4[1][P('w')]
                x3
                sage: t4[1][P('z')]
                (x2*b + x3*b)/x4

            The one-element poset::

                sage: P = Poset({8: []})
                sage: t = P.birational_free_labelling()
                sage: t8 = P.birational_toggle(8, t); t8
                (Fraction Field of Multivariate Polynomial Ring in a, x1, b over Rational Field,
                 {...},
                 a,
                 b)
                sage: t8[1][8]
                a*b/x1
            """
            FF = labelling[0]       # base field
            a = labelling[2]        # label at `0 \in \widehat{P}`
            b = labelling[3]
            newdict = labelling[1].copy()
            # Construct the harmonic sum ``x`` of the labels at the
            # elements covering ``v``:
            uppers = self.upper_covers(v)
            if len(uppers) == 0:
                x = FF.one() / b
            else:
                x = FF.sum(FF.one() / newdict[j] for j in uppers)
                # ``FF.sum``, not ``sum``, see trac #15591.
            x = FF.one() / x
            # Construct the sum ``y`` of the labels at the elements
            # covered by ``v``:
            lowers = self.lower_covers(v)
            if len(lowers) == 0:
                y = a
            else:
                y = FF.sum(newdict[j] for j in lowers)
            # Now, transform the label at v:
            newdict[v] = x * y / newdict[v]
            return (FF, newdict, a, b)

        def birational_toggles(self, vs, labelling):
            r"""
            Return the result of applying a sequence of birational
            toggles (specified by ``vs``) to the `\mathbf{K}`-labelling
            ``labelling`` of the poset ``self``.

            See the documentation of :meth:`birational_free_labelling`
            for a definition of birational toggles and `\mathbf{K}`-labellings
            and for an explanation of how `\mathbf{K}`-labellings are to be
            encoded to be understood by Sage. This implementation allows
            `\mathbf{K}` to be a semifield, not just a field. The birational
            `v`-toggle is only a rational map, so an exception (most
            likely, ``ZeroDivisionError``) will be thrown if the
            denominator is zero.

            INPUT:

            - ``vs`` -- an iterable comprising elements of ``self``
              (which must have ``self`` as parent if ``self`` is a
              ``facade=False`` poset)

            - ``labelling`` -- a `\mathbf{K}`-labelling of ``self`` in the
              sense as defined in the documentation of
              :meth:`birational_free_labelling`

            OUTPUT:

            The `\mathbf{K}`-labelling `T_{v_n} T_{v_{n-1}} \cdots T_{v_1} f`
            of ``self``, where `f` is ``labelling`` and
            `(v_1, v_2, \ldots, v_n)` is ``vs`` (written as list).

            EXAMPLES::

                sage: P = posets.SymmetricGroupBruhatOrderPoset(3)
                sage: sorted(list(P))
                ['123', '132', '213', '231', '312', '321']
                sage: TT = TropicalSemiring(ZZ)
                sage: t = (TT, {'123': TT(4), '132': TT(2), '213': TT(3), '231': TT(1), '321': TT(1), '312': TT(2)}, TT(7), TT(1))
                sage: tA = P.birational_toggles(['123', '231', '312'], t); tA
                (Tropical semiring over Integer Ring, {...}, 7, 1)
                sage: sorted(tA[1].items())
                [('123', 6), ('132', 2), ('213', 3), ('231', 2), ('312', 1), ('321', 1)]
                sage: tAB = P.birational_toggles(['132', '213', '321'], tA); tAB
                (Tropical semiring over Integer Ring, {...}, 7, 1)
                sage: sorted(tAB[1].items())
                [('123', 6), ('132', 6), ('213', 5), ('231', 2), ('312', 1), ('321', 1)]

                sage: P = Poset({1: [2, 3], 2: [4], 3: [4]})
                sage: Qx = PolynomialRing(QQ, 'x').fraction_field()
                sage: x = Qx.gen()
                sage: t = (Qx, {1: 1, 2: x, 3: (x+1)/x, 4: x^2}, 1, 1)
                sage: t1 = P.birational_toggles((i for i in range(1, 5)), t); t1
                (Fraction Field of Univariate Polynomial Ring in x over Rational Field,
                 {...},
                 1,
                 1)
                sage: sorted(t1[1].items())
                [(1, (x^2 + x)/(x^2 + x + 1)), (2, (x^3 + x^2)/(x^2 + x + 1)), (3, x^4/(x^2 + x + 1)), (4, 1)]
                sage: t2 = P.birational_toggles(reversed(range(1, 5)), t)
                sage: sorted(t2[1].items())
                [(1, 1/x^2), (2, (x^2 + x + 1)/x^4), (3, (x^2 + x + 1)/(x^3 + x^2)), (4, (x^2 + x + 1)/x^3)]

            Facade set to ``False`` works::

                sage: P = Poset({'x': ['y', 'w'], 'y': ['z'], 'w': ['z']}, facade=False)
                sage: lex = ['x', 'y', 'w', 'z']
                sage: t = P.birational_free_labelling(linear_extension=lex)
                sage: sorted(P.birational_toggles([P('x'), P('y')], t)[1].items())
                [(x, a*x2*x3/(x1*x2 + x1*x3)), (y, a*x3*x4/(x1*x2 + x1*x3)), (w, x3), (z, x4)]
            """
            l = labelling
            for v in vs:
                l = self.birational_toggle(v, l)
            return l

        def birational_rowmotion(self, labelling):
            r"""
            Return the result of applying birational rowmotion to the
            `\mathbf{K}`-labelling ``labelling`` of the poset ``self``.

            See the documentation of :meth:`birational_free_labelling`
            for a definition of birational rowmotion and
            `\mathbf{K}`-labellings and for an explanation of how
            `\mathbf{K}`-labellings are to be encoded to be understood
            by Sage. This implementation allows `\mathbf{K}` to be a
            semifield, not just a field. Birational rowmotion is only a
            rational map, so an exception (most likely, ``ZeroDivisionError``)
            will be thrown if the denominator is zero.

            INPUT:

            - ``labelling`` -- a `\mathbf{K}`-labelling of ``self`` in the
              sense as defined in the documentation of
              :meth:`birational_free_labelling`

            OUTPUT:

            The image of the `\mathbf{K}`-labelling `f` under birational
            rowmotion.

            EXAMPLES::

                sage: P = Poset({1: [2, 3], 2: [4], 3: [4]})
                sage: lex = [1, 2, 3, 4]
                sage: t = P.birational_free_labelling(linear_extension=lex); t
                (Fraction Field of Multivariate Polynomial Ring in a, x1, x2, x3, x4, b over Rational Field,
                 {...},
                 a,
                 b)
                sage: sorted(t[1].items())
                [(1, x1), (2, x2), (3, x3), (4, x4)]
                sage: t = P.birational_rowmotion(t); t
                (Fraction Field of Multivariate Polynomial Ring in a, x1, x2, x3, x4, b over Rational Field,
                 {...},
                 a,
                 b)
                sage: sorted(t[1].items())
                [(1, a*b/x4), (2, (x1*x2*b + x1*x3*b)/(x2*x4)),
                 (3, (x1*x2*b + x1*x3*b)/(x3*x4)), (4, (x2*b + x3*b)/x4)]

            A result of [GR2013]_ states that applying birational rowmotion
            `n+m` times to a `\mathbf{K}`-labelling `f` of the poset
            `[n] \times [m]` gives back `f`. Let us check this::

                sage: def test_rectangle_periodicity(n, m, k):
                ....:     P = posets.ChainPoset(n).product(posets.ChainPoset(m))
                ....:     t0 = P.birational_free_labelling(P)
                ....:     t = t0
                ....:     for i in range(k):
                ....:         t = P.birational_rowmotion(t)
                ....:     return t == t0
                sage: test_rectangle_periodicity(2, 2, 4)
                True
                sage: test_rectangle_periodicity(2, 2, 2)
                False
                sage: test_rectangle_periodicity(2, 3, 5)  # long time
                True

            While computations with the birational free labelling quickly
            run out of memory due to the complexity of the rational
            functions involved, it is computationally cheap to check
            properties of birational rowmotion on examples in the tropical
            semiring::

                sage: def test_rectangle_periodicity_tropical(n, m, k):
                ....:     P = posets.ChainPoset(n).product(posets.ChainPoset(m))
                ....:     TT = TropicalSemiring(ZZ)
                ....:     t0 = (TT, {v: TT(floor(random()*100)) for v in P}, TT(0), TT(124))
                ....:     t = t0
                ....:     for i in range(k):
                ....:         t = P.birational_rowmotion(t)
                ....:     return t == t0
                sage: test_rectangle_periodicity_tropical(7, 6, 13)
                True

            Tropicalization is also what relates birational rowmotion to
            classical rowmotion on order ideals. In fact, if `T` denotes
            the :class:`tropical semiring
            <sage.rings.semirings.tropical_semiring.TropicalSemiring>` of
            `\ZZ` and `P` is a finite poset, then we can define an embedding
            `\phi` from the set `J(P)` of all order ideals of `P` into the
            set `T^{\widehat{P}}` of all `T`-labellings of `P` by sending
            every `I \in J(P)` to the indicator function of `I` extended by
            the value `1` at the element `0` and the value `0` at the
            element `1`. This map `\phi` has the property that
            `R \circ \phi = \phi \circ r`, where `R` denotes birational
            rowmotion, and `r` denotes :meth:`classical rowmotion <rowmotion>`
            on `J(P)`. An example::

                sage: P = posets.IntegerPartitions(5)
                sage: TT = TropicalSemiring(ZZ)
                sage: def indicator_labelling(I):
                ....:     # send order ideal `I` to a `T`-labelling of `P`.
                ....:     dct = {v: TT(v in I) for v in P}
                ....:     return (TT, dct, TT(1), TT(0))
                sage: all(indicator_labelling(P.rowmotion(I))
                ....:     == P.birational_rowmotion(indicator_labelling(I))
                ....:     for I in P.order_ideals_lattice(facade=True))
                True

            TESTS:

            Facade set to false works::

                sage: P = Poset({1: [2, 3], 2: [4], 3: [4]}, facade=False)
                sage: lex = [1, 2, 3, 4]
                sage: t = P.birational_free_labelling(linear_extension=lex); t
                (Fraction Field of Multivariate Polynomial Ring in a, x1, x2, x3, x4, b over Rational Field,
                 {...},
                 a,
                 b)
                sage: t = P.birational_rowmotion(t); t
                (Fraction Field of Multivariate Polynomial Ring in a, x1, x2, x3, x4, b over Rational Field,
                 {...},
                 a,
                 b)
                sage: t[1][P(2)]
                (x1*x2*b + x1*x3*b)/(x2*x4)
                sage: t = P.birational_rowmotion(t)
                sage: t[1][P(2)]
                a*b/x3
            """
            l = labelling
            for v in reversed(self.linear_extension()):
                l = self.birational_toggle(v, l)
            return l

        def panyushev_orbits(self, element_constructor = set):
            r"""
            Return the Panyushev orbits of antichains in ``self``.

            The Panyushev orbit of an antichain is its orbit under
            Panyushev complementation (see
            :meth:`panyushev_complement`).

            INPUT:

            - ``element_constructor`` (defaults to ``set``) -- a type
              constructor (``set``, ``tuple``, ``list``, ``frozenset``,
              ``iter``, etc.) which is to be applied to the antichains
              before they are returned.

            OUTPUT:

            - the partition of the set of all antichains of ``self`` into
              orbits under Panyushev complementation. This is returned as
              a list of lists ``L`` such that for each ``L`` and ``i``,
              cyclically:
              ``self.order_ideal_complement_generators(L[i]) == L[i+1]``.
              The entries ``L[i]`` are sets by default, but depending on
              the optional keyword variable ``element_constructors``
              they can also be tuples, lists etc.

            EXAMPLES::

                sage: P = Poset( ( [1,2,3], [ [1,3], [2,3] ] ) )
                sage: orb = P.panyushev_orbits()
                sage: sorted(sorted(o) for o in orb)
                [[set(), {1, 2}, {3}], [{2}, {1}]]
                sage: orb = P.panyushev_orbits(element_constructor=list)
                sage: sorted(sorted(o) for o in orb)
                [[[], [1, 2], [3]], [[1], [2]]]
                sage: orb = P.panyushev_orbits(element_constructor=frozenset)
                sage: sorted(sorted(o) for o in orb)
                [[frozenset(), frozenset({1, 2}), frozenset({3})],
                 [frozenset({2}), frozenset({1})]]
                sage: orb = P.panyushev_orbits(element_constructor=tuple)
                sage: sorted(sorted(o) for o in orb)
                [[(), (1, 2), (3,)], [(1,), (2,)]]
                sage: P = Poset( {} )
                sage: P.panyushev_orbits()
                [[set()]]
            """
            # TODO: implement a generic function taking a set and
            # bijections on this set, and returning the orbits.
            AC = set(self.antichains(element_constructor = frozenset))
            orbits = []
            while AC:
                A = AC.pop()
                orbit = [ A ]
                while True:
                    A = frozenset(self.order_ideal_complement_generators(A))
                    if A not in AC:
                        break
                    orbit.append( A )
                    AC.remove( A )
                orbits.append([element_constructor(_) for _ in orbit])
            return orbits

        def rowmotion_orbits(self, element_constructor = set):
            r"""
            Return the rowmotion orbits of order ideals in ``self``.

            The rowmotion orbit of an order ideal is its orbit under
            rowmotion (see :meth:`rowmotion`).

            INPUT:

            - ``element_constructor`` (defaults to ``set``) -- a type
              constructor (``set``, ``tuple``, ``list``, ``frozenset``,
              ``iter``, etc.) which is to be applied to the antichains
              before they are returned.

            OUTPUT:

            - the partition of the set of all order ideals of ``self``
              into orbits under rowmotion. This is returned as
              a list of lists ``L`` such that for each ``L`` and ``i``,
              cyclically: ``self.rowmotion(L[i]) == L[i+1]``.
              The entries ``L[i]`` are sets by default, but depending on
              the optional keyword variable ``element_constructors``
              they can also be tuples, lists etc.

            EXAMPLES::

                sage: P = Poset( {1: [2, 3], 2: [], 3: [], 4: [2]} )
                sage: sorted(len(o) for o in P.rowmotion_orbits())
                [3, 5]
                sage: orb = P.rowmotion_orbits(element_constructor=list)
                sage: sorted(sorted(e) for e in orb)
                [[[], [4, 1], [4, 1, 2, 3]], [[1], [1, 3], [4], [4, 1, 2], [4, 1, 3]]]
                sage: orb = P.rowmotion_orbits(element_constructor=tuple)
                sage: sorted(sorted(e) for e in orb)
                [[(), (4, 1), (4, 1, 2, 3)], [(1,), (1, 3), (4,), (4, 1, 2), (4, 1, 3)]]
                sage: P = Poset({})
                sage: P.rowmotion_orbits(element_constructor=tuple)
                [[()]]
            """
            pan_orbits = self.panyushev_orbits(element_constructor = list)
            return [[element_constructor(self.order_ideal(oideal)) for oideal in orbit] for orbit in pan_orbits]

        def rowmotion_orbits_plots(self):
            r"""
            Return plots of the rowmotion orbits of order ideals in ``self``.

            The rowmotion orbit of an order ideal is its orbit under
            rowmotion (see :meth:`rowmotion`).

            EXAMPLES::

                sage: P = Poset( {1: [2, 3], 2: [], 3: [], 4: [2]} )
                sage: P.rowmotion_orbits_plots()
                Graphics Array of size 2 x 5
                sage: P = Poset({})
                sage: P.rowmotion_orbits_plots()
                Graphics Array of size 1 x 1

            """
            from sage.plot.plot import graphics_array
            plot_of_orb_plots=[]
            max_orbit_size = 0
            for orb in self.rowmotion_orbits():
                orb_plots=[]
                if len(orb) > max_orbit_size:
                    max_orbit_size = len(orb)
                for oi in orb:
                    oiplot = self.order_ideal_plot(oi)
                    orb_plots.append(oiplot)
                plot_of_orb_plots.append(orb_plots)
            return graphics_array(plot_of_orb_plots, ncols = max_orbit_size)


        def toggling_orbits(self, vs, element_constructor = set):
            r"""
            Return the orbits of order ideals in ``self`` under the
            operation of toggling the vertices ``vs[0], vs[1], ...``
            in this order.

            See :meth:`~sage.categories.posets.Posets.ParentMethods.order_ideal_toggle` for a definition of toggling.

            .. WARNING::

                The orbits are those under the composition of toggles,
                *not* under the single toggles themselves. Thus, for
                example, if ``vs == [1,2]``, then the orbits have the
                form `(I, T_2 T_1 I, T_2 T_1 T_2 T_1 I, \ldots)`
                (where `I` denotes an order ideal and `T_i` means
                toggling at `i`) rather than
                `(I, T_1 I, T_2 T_1 I, T_1 T_2 T_1 I, \ldots)`.

            INPUT:

            - ``vs``: a list (or other iterable) of elements of ``self``
              (but since the output depends on the order, sets should
              not be used as ``vs``).

            OUTPUT:

            - a partition of the order ideals of ``self``, as a list of
              sets ``L`` such that for each ``L`` and ``i``, cyclically:
              ``self.order_ideal_toggles(L[i], vs) == L[i+1]``.

            EXAMPLES::

                sage: P = Poset( {1: [2, 4], 2: [], 3: [4], 4: []} )
                sage: sorted(len(o) for o in P.toggling_orbits([1, 2]))
                [2, 3, 3]
                sage: P = Poset( {1: [3], 2: [1, 4], 3: [], 4: [3]} )
                sage: sorted(len(o) for o in P.toggling_orbits((1, 2, 4, 3)))
                [3, 3]
            """
            # TODO: implement a generic function taking a set and
            # bijections on this set, and returning the orbits.
            OI = set(self.order_ideals_lattice(facade=True))
            orbits = []
            while OI:
                A = OI.pop()
                orbit = [ A ]
                while True:
                    A = self.order_ideal_toggles(A, vs)
                    if A not in OI:
                        break
                    orbit.append( A )
                    OI.remove( A )
                orbits.append([element_constructor(_) for _ in orbit])
            return orbits

        def toggling_orbits_plots(self, vs):
            r"""
            Return plots of the orbits of order ideals in ``self`` under the
            operation of toggling the vertices ``vs[0], vs[1], ...``
            in this order.

            See :meth:`toggling_orbits` for more information.

            EXAMPLES::

                sage: P = Poset( {1: [2, 3], 2: [], 3: [], 4: [2]} )
                sage: P.toggling_orbits_plots([1,2,3,4])
                Graphics Array of size 2 x 5
                sage: P = Poset({})
                sage: P.toggling_orbits_plots([])
                Graphics Array of size 1 x 1

            """
            from sage.plot.plot import graphics_array
            plot_of_orb_plots=[]
            max_orbit_size = 0
            for orb in self.toggling_orbits(vs):
                orb_plots=[]
                if len(orb) > max_orbit_size:
                    max_orbit_size = len(orb)
                for oi in orb:
                    oiplot = self.order_ideal_plot(oi)
                    orb_plots.append(oiplot)
                plot_of_orb_plots.append(orb_plots)
            return graphics_array(plot_of_orb_plots, ncols = max_orbit_size)

        def panyushev_orbit_iter(self, antichain, element_constructor=set, stop=True, check=True):
            r"""
            Iterate over the Panyushev orbit of an antichain
            ``antichain`` of ``self``.

            The Panyushev orbit of an antichain is its orbit under
            Panyushev complementation (see
            :meth:`panyushev_complement`).

            INPUT:

            - ``antichain`` -- an antichain of ``self``, given as an
              iterable.

            - ``element_constructor`` (defaults to ``set``) -- a type
              constructor (``set``, ``tuple``, ``list``, ``frozenset``,
              ``iter``, etc.) which is to be applied to the antichains
              before they are yielded.

            - ``stop`` -- a Boolean (default: ``True``) determining
              whether the iterator should stop once it completes its
              cycle (this happens when it is set to ``True``) or go on
              forever (this happens when it is set to ``False``).

            - ``check`` -- a Boolean (default: ``True``) determining
              whether ``antichain`` should be checked for being an
              antichain.

            OUTPUT:

            - an iterator over the orbit of the antichain ``antichain``
              under Panyushev complementation. This iterator `I` has the
              property that ``I[0] == antichain`` and each `i` satisfies
              ``self.order_ideal_complement_generators(I[i]) == I[i+1]``,
              where ``I[i+1]`` has to be understood as ``I[0]`` if it is
              undefined.
              The entries ``I[i]`` are sets by default, but depending on
              the optional keyword variable ``element_constructors``
              they can also be tuples, lists etc.

            EXAMPLES::

                sage: P = Poset( ( [1,2,3], [ [1,3], [2,3] ] ) )
                sage: list(P.panyushev_orbit_iter(set([1, 2])))
                [{1, 2}, {3}, set()]
                sage: list(P.panyushev_orbit_iter([1, 2]))
                [{1, 2}, {3}, set()]
                sage: list(P.panyushev_orbit_iter([2, 1]))
                [{1, 2}, {3}, set()]
                sage: list(P.panyushev_orbit_iter(set([1, 2]), element_constructor=list))
                [[1, 2], [3], []]
                sage: list(P.panyushev_orbit_iter(set([1, 2]), element_constructor=frozenset))
                [frozenset({1, 2}), frozenset({3}), frozenset()]
                sage: list(P.panyushev_orbit_iter(set([1, 2]), element_constructor=tuple))
                [(1, 2), (3,), ()]

                sage: P = Poset( {} )
                sage: list(P.panyushev_orbit_iter([]))
                [set()]

                sage: P = Poset({ 1: [2, 3], 2: [4], 3: [4], 4: [] })
                sage: Piter = P.panyushev_orbit_iter([2], stop=False)
                sage: next(Piter)
                {2}
                sage: next(Piter)
                {3}
                sage: next(Piter)
                {2}
                sage: next(Piter)
                {3}
            """
            # TODO: implement a generic function taking a set and
            # bijections on this set, and returning an orbit of a given
            # element.
            if check:
                if not self.is_antichain_of_poset(antichain):
                    raise ValueError("the given antichain is not an antichain")
            starter = set(antichain)     # sanitize input
            yield element_constructor(starter)
            next = starter
            if stop:
                while True:
                    next = self.order_ideal_complement_generators(next)
                    if next == starter:
                        break
                    yield element_constructor(next)
            else:
                while True:
                    next = self.order_ideal_complement_generators(next)
                    yield element_constructor(next)

        def rowmotion_orbit_iter(self, oideal, element_constructor=set, stop=True, check=True):
            r"""
            Iterate over the rowmotion orbit of an order ideal
            ``oideal`` of ``self``.

            The rowmotion orbit of an order ideal is its orbit under
            rowmotion (see :meth:`rowmotion`).

            INPUT:

            - ``oideal`` -- an order ideal of ``self``, given as an
              iterable.

            - ``element_constructor`` (defaults to ``set``) -- a type
              constructor (``set``, ``tuple``, ``list``, ``frozenset``,
              ``iter``, etc.) which is to be applied to the order
              ideals before they are yielded.

            - ``stop`` -- a Boolean (default: ``True``) determining
              whether the iterator should stop once it completes its
              cycle (this happens when it is set to ``True``) or go on
              forever (this happens when it is set to ``False``).

            - ``check`` -- a Boolean (default: ``True``) determining
              whether ``oideal`` should be checked for being an
              order ideal.

            OUTPUT:

            - an iterator over the orbit of the order ideal ``oideal``
              under rowmotion. This iterator `I` has the property that
              ``I[0] == oideal`` and that every `i` satisfies
              ``self.rowmotion(I[i]) == I[i+1]``, where ``I[i+1]`` has
              to be understood as ``I[0]`` if it is undefined.
              The entries ``I[i]`` are sets by default, but depending on
              the optional keyword variable ``element_constructors``
              they can also be tuples, lists etc.

            EXAMPLES::

                sage: P = Poset( ( [1,2,3], [ [1,3], [2,3] ] ) )
                sage: list(P.rowmotion_orbit_iter(set([1, 2])))
                [{1, 2}, {1, 2, 3}, set()]
                sage: list(P.rowmotion_orbit_iter([1, 2]))
                [{1, 2}, {1, 2, 3}, set()]
                sage: list(P.rowmotion_orbit_iter([2, 1]))
                [{1, 2}, {1, 2, 3}, set()]
                sage: list(P.rowmotion_orbit_iter(set([1, 2]), element_constructor=list))
                [[1, 2], [1, 2, 3], []]
                sage: list(P.rowmotion_orbit_iter(set([1, 2]), element_constructor=frozenset))
                [frozenset({1, 2}), frozenset({1, 2, 3}), frozenset()]
                sage: list(P.rowmotion_orbit_iter(set([1, 2]), element_constructor=tuple))
                [(1, 2), (1, 2, 3), ()]

                sage: P = Poset( {} )
                sage: list(P.rowmotion_orbit_iter([]))
                [set()]

                sage: P = Poset({ 1: [2, 3], 2: [4], 3: [4], 4: [] })
                sage: Piter = P.rowmotion_orbit_iter([1, 2, 3], stop=False)
                sage: next(Piter)
                {1, 2, 3}
                sage: next(Piter)
                {1, 2, 3, 4}
                sage: next(Piter)
                set()
                sage: next(Piter)
                {1}
                sage: next(Piter)
                {1, 2, 3}

                sage: P = Poset({ 1: [4], 2: [4, 5], 3: [5] })
                sage: list(P.rowmotion_orbit_iter([1, 2], element_constructor=list))
                [[1, 2], [1, 2, 3, 4], [2, 3, 5], [1], [2, 3], [1, 2, 3, 5], [1, 2, 4], [3]]
            """
            # TODO: implement a generic function taking a set and
            # bijections on this set, and returning an orbit of a given
            # element.
            if check:
                if not self.is_order_ideal(oideal):
                    raise ValueError("the given order ideal is not an order ideal")
            starter = set(oideal)     # sanitize input
            yield element_constructor(starter)
            next = starter
            if stop:
                while True:
                    next = self.rowmotion(next)
                    if next == starter:
                        break
                    yield element_constructor(next)
            else:
                while True:
                    next = self.rowmotion(next)
                    yield element_constructor(next)

        def toggling_orbit_iter(self, vs, oideal, element_constructor=set, stop=True, check=True):
            r"""
            Iterate over the orbit of an order ideal ``oideal`` of
            ``self`` under the operation of toggling the vertices
            ``vs[0], vs[1], ...`` in this order.

            See :meth:`~sage.categories.posets.Posets.ParentMethods.order_ideal_toggle` for a definition of toggling.

            .. WARNING::

                The orbit is that under the composition of toggles,
                *not* under the single toggles themselves. Thus, for
                example, if ``vs == [1,2]``, then the orbit has the
                form `(I, T_2 T_1 I, T_2 T_1 T_2 T_1 I, \ldots)`
                (where `I` denotes ``oideal`` and `T_i` means
                toggling at `i`) rather than
                `(I, T_1 I, T_2 T_1 I, T_1 T_2 T_1 I, \ldots)`.

            INPUT:

            - ``vs``: a list (or other iterable) of elements of ``self``
              (but since the output depends on the order, sets should
              not be used as ``vs``).

            - ``oideal`` -- an order ideal of ``self``, given as an
              iterable.

            - ``element_constructor`` (defaults to ``set``) -- a type
              constructor (``set``, ``tuple``, ``list``, ``frozenset``,
              ``iter``, etc.) which is to be applied to the order
              ideals before they are yielded.

            - ``stop`` -- a Boolean (default: ``True``) determining
              whether the iterator should stop once it completes its
              cycle (this happens when it is set to ``True``) or go on
              forever (this happens when it is set to ``False``).

            - ``check`` -- a Boolean (default: ``True``) determining
              whether ``oideal`` should be checked for being an
              order ideal.

            OUTPUT:

            - an iterator over the orbit of the order ideal ``oideal``
              under toggling the vertices in the list ``vs`` in this
              order. This iterator `I` has the property that
              ``I[0] == oideal`` and that every `i` satisfies
              ``self.order_ideal_toggles(I[i], vs) == I[i+1]``, where
              ``I[i+1]`` has to be understood as ``I[0]`` if it is
              undefined.
              The entries ``I[i]`` are sets by default, but depending on
              the optional keyword variable ``element_constructors``
              they can also be tuples, lists etc.

            EXAMPLES::

                sage: P = Poset( ( [1,2,3], [ [1,3], [2,3] ] ) )
                sage: list(P.toggling_orbit_iter([1, 3, 1], set([1, 2])))
                [{1, 2}]
                sage: list(P.toggling_orbit_iter([1, 2, 3], set([1, 2])))
                [{1, 2}, set(), {1, 2, 3}]
                sage: list(P.toggling_orbit_iter([3, 2, 1], set([1, 2])))
                [{1, 2}, {1, 2, 3}, set()]
                sage: list(P.toggling_orbit_iter([3, 2, 1], set([1, 2]), element_constructor=list))
                [[1, 2], [1, 2, 3], []]
                sage: list(P.toggling_orbit_iter([3, 2, 1], set([1, 2]), element_constructor=frozenset))
                [frozenset({1, 2}), frozenset({1, 2, 3}), frozenset()]
                sage: list(P.toggling_orbit_iter([3, 2, 1], set([1, 2]), element_constructor=tuple))
                [(1, 2), (1, 2, 3), ()]
                sage: list(P.toggling_orbit_iter([3, 2, 1], [2, 1], element_constructor=tuple))
                [(1, 2), (1, 2, 3), ()]

                sage: P = Poset( {} )
                sage: list(P.toggling_orbit_iter([], []))
                [set()]

                sage: P = Poset({ 1: [2, 3], 2: [4], 3: [4], 4: [] })
                sage: Piter = P.toggling_orbit_iter([1, 2, 4, 3], [1, 2, 3], stop=False)
                sage: next(Piter)
                {1, 2, 3}
                sage: next(Piter)
                {1}
                sage: next(Piter)
                set()
                sage: next(Piter)
                {1, 2, 3}
                sage: next(Piter)
                {1}
            """
            # TODO: implement a generic function taking a set and
            # bijections on this set, and returning an orbit of a given
            # element.
            if check:
                if not self.is_order_ideal(oideal):
                    raise ValueError("the given order ideal is not an order ideal")
            starter = set(oideal)     # sanitize input
            yield element_constructor(starter)
            next = starter
            if stop:
                while True:
                    next = self.order_ideal_toggles(next, vs)
                    if next == starter:
                        break
                    yield element_constructor(next)
            else:
                while True:
                    next = self.order_ideal_toggles(next, vs)
                    yield element_constructor(next)

        def order_ideals_lattice(self, as_ideals=True, facade=None):
            r"""
            Return the lattice of order ideals of a poset ``self``,
            ordered by inclusion.

            The lattice of order ideals of a poset `P` is usually
            denoted by `J(P)`. Its underlying set is the set of order
            ideals of `P`, and its partial order is given by
            inclusion.

            The order ideals of `P` are in a canonical bijection
            with the antichains of `P`. The bijection maps every
            order ideal to the antichain formed by its maximal
            elements. By setting the ``as_ideals`` keyword variable to
            ``False``, one can make this method apply this bijection
            before returning the lattice.

            INPUT:

            - ``as_ideals`` -- Boolean, if ``True`` (default) returns
              a poset on the set of order ideals, otherwise on the set
              of antichains
            - ``facade`` -- Boolean or ``None`` (default). Whether to
              return a facade lattice or not. By default return facade
              lattice if the poset is a facade poset.

            EXAMPLES::

                sage: P = posets.PentagonPoset()
                sage: P.cover_relations()
                [[0, 1], [0, 2], [1, 4], [2, 3], [3, 4]]
                sage: J = P.order_ideals_lattice(); J
                Finite lattice containing 8 elements
                sage: sorted(sorted(e) for e in J)
                 [[], [0], [0, 1], [0, 1, 2], [0, 1, 2, 3], [0, 1, 2, 3, 4], [0, 2], [0, 2, 3]]

            As a lattice on antichains::

                sage: J2 = P.order_ideals_lattice(False); J2
                Finite lattice containing 8 elements
                sage: sorted(J2)
                [(), (0,), (1,), (1, 2), (1, 3), (2,), (3,), (4,)]

            TESTS::

                sage: J = posets.DiamondPoset(4, facade = True).order_ideals_lattice(); J
                Finite lattice containing 6 elements
                sage: sorted(sorted(e) for e in J)
                [[], [0], [0, 1], [0, 1, 2], [0, 1, 2, 3], [0, 2]]
                sage: sorted(sorted(sorted(e) for e in c) for c in J.cover_relations())
                [[[], [0]], [[0], [0, 1]], [[0], [0, 2]], [[0, 1], [0, 1, 2]], [[0, 1, 2], [0, 1, 2, 3]], [[0, 1, 2], [0, 2]]]

                sage: P = Poset({1:[2]})
                sage: J_facade = P.order_ideals_lattice()
                sage: J_nonfacade = P.order_ideals_lattice(facade=False)
                sage: type(J_facade[0]) == type(J_nonfacade[0])
                False
            """
            from sage.combinat.posets.lattices import LatticePoset
            if facade is None:
                facade = self._is_facade
            if as_ideals:
                from sage.misc.call import attrcall
                from sage.sets.set import Set
                ideals = [Set(self.order_ideal(antichain))
                          for antichain in self.antichains()]
                return LatticePoset((ideals, attrcall("issubset")),
                                    facade=facade)
            else:
                from sage.misc.cachefunc import cached_function
                antichains = [tuple(a) for a in self.antichains()]
                @cached_function
                def is_above(a, xb):
                    return any(self.is_lequal(xa, xb) for xa in a)
                def compare(a, b):
                    return all(is_above(a, xb) for xb in b)
                return LatticePoset((antichains, compare), facade=facade)

        @abstract_method(optional = True)
        def antichains(self):
            r"""
            Return all antichains of ``self``.

            EXAMPLES::

                sage: A = posets.PentagonPoset().antichains(); A
                Set of antichains of Finite lattice containing 5 elements
                sage: list(A)
                [[], [0], [1], [1, 2], [1, 3], [2], [3], [4]]
            """

        def directed_subsets(self, direction):
            r"""
            Return the order filters (resp. order ideals) of ``self``, as lists.

            If ``direction`` is 'up', returns the order filters (upper sets).

            If ``direction`` is 'down', returns the order ideals (lower sets).

            INPUT:

            - ``direction`` -- 'up' or 'down'

            EXAMPLES::

                sage: P = Poset((divisors(12), attrcall("divides")), facade=True)
                sage: A = P.directed_subsets('up')
                sage: sorted(list(A))
                [[], [1, 2, 4, 3, 6, 12], [2, 4, 3, 6, 12], [2, 4, 6, 12], [3, 6, 12], [4, 3, 6, 12], [4, 6, 12], [4, 12], [6, 12], [12]]

            TESTS::

                sage: list(Poset().directed_subsets('up'))
                [[]]
            """
            if direction != 'up' and direction != 'down':
                raise ValueError("Direction must be either 'up' or 'down'.")
            return self.antichains().map(lambda elements: self.directed_subset(elements, direction))
