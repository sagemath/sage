# sage_setup: distribution = sagemath-categories
r"""
Filtered Modules With Basis

A *filtered module with basis* over a ring `R` means
(for the purpose of this code) a filtered `R`-module `M`
with filtration `(F_i)_{i \in I}` (typically `I = \NN`)
endowed with a basis `(b_j)_{j \in J}` of `M` and a partition
`J = \bigsqcup_{i \in I} J_i` of the set `J` (it is allowed
that some `J_i` are empty) such that for every `n \in I`,
the subfamily `(b_j)_{j \in U_n}`, where
`U_n = \bigcup_{i \leq n} J_i`, is a basis of the
`R`-submodule `F_n`.

For every `i \in I`, the `R`-submodule of `M` spanned by
`(b_j)_{j \in J_i}` is called the `i`-*th graded component*
(aka the `i`-*th homogeneous component*) of the filtered
module with basis `M`; the elements of this submodule are
referred to as *homogeneous elements of degree* `i`.

See the class documentation
:class:`~sage.categories.filtered_modules_with_basis.FilteredModulesWithBasis`
for further details.
"""
#*****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.filtered_modules import FilteredModulesCategory
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.subobjects import SubobjectsCategory
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring


class FilteredModulesWithBasis(FilteredModulesCategory):
    r"""
    The category of filtered modules with a distinguished basis.

    A *filtered module with basis* over a ring `R` means
    (for the purpose of this code) a filtered `R`-module `M`
    with filtration `(F_i)_{i \in I}` (typically `I = \NN`)
    endowed with a basis `(b_j)_{j \in J}` of `M` and a partition
    `J = \bigsqcup_{i \in I} J_i` of the set `J` (it is allowed
    that some `J_i` are empty) such that for every `n \in I`,
    the subfamily `(b_j)_{j \in U_n}`, where
    `U_n = \bigcup_{i \leq n} J_i`, is a basis of the
    `R`-submodule `F_n`.

    For every `i \in I`, the `R`-submodule of `M` spanned by
    `(b_j)_{j \in J_i}` is called the `i`-*th graded component*
    (aka the `i`-*th homogeneous component*) of the filtered
    module with basis `M`; the elements of this submodule are
    referred to as *homogeneous elements of degree* `i`.
    The `R`-module `M` is the direct sum of its `i`-th graded
    components over all `i \in I`, and thus becomes a graded
    `R`-module with basis.
    Conversely, any graded `R`-module with basis canonically
    becomes a filtered `R`-module with basis (by defining
    `F_n = \bigoplus_{i \leq n} G_i` where `G_i` is the `i`-th
    graded component, and defining `J_i` as the indexing set
    of the basis of the `i`-th graded component). Hence, the
    notion of a filtered `R`-module with basis is equivalent
    to the notion of a graded `R`-module with basis.

    However, the *category* of filtered `R`-modules with basis is not
    the category of graded `R`-modules with basis. Indeed, the *morphisms*
    of filtered `R`-modules with basis are defined to be morphisms of
    `R`-modules which send each `F_n` of the domain to the corresponding
    `F_n` of the target; in contrast, the morphisms of graded `R`-modules
    with basis must preserve each homogeneous component. Also,
    the notion of a filtered algebra with basis differs from
    that of a graded algebra with basis.

    .. NOTE::

        Currently, to make use of the functionality of this class,
        an instance of ``FilteredModulesWithBasis`` should fulfill
        the contract of a :class:`CombinatorialFreeModule` (most
        likely by inheriting from it). It should also have the
        indexing set `J` encoded as its ``_indices`` attribute,
        and ``_indices.subset(size=i)`` should yield the subset
        `J_i` (as an iterable). If the latter conditions are not
        satisfied, then :meth:`basis` must be overridden.

    .. NOTE::

        One should implement a ``degree_on_basis`` method in the parent
        class in order to fully utilize the methods of this category.
        This might become a required abstract method in the future.

    EXAMPLES::

        sage: C = ModulesWithBasis(ZZ).Filtered(); C
        Category of filtered modules with basis over Integer Ring
        sage: sorted(C.super_categories(), key=str)
        [Category of filtered modules over Integer Ring,
         Category of modules with basis over Integer Ring]
        sage: C is ModulesWithBasis(ZZ).Filtered()
        True

    TESTS::

        sage: C = ModulesWithBasis(ZZ).Filtered()
        sage: TestSuite(C).run()
        sage: C = ModulesWithBasis(QQ).Filtered()
        sage: TestSuite(C).run()
    """
    class ParentMethods:
        # TODO: which syntax do we prefer?
        # A.basis(degree = 3)
        # A.basis().subset(degree=3)

        # This is related to the following design question:
        # If F = (f_i)_{i\in I} is a family, should ``F.subset(degree = 3)``
        # be the elements of F of degree 3 or those whose index is of degree 3?

        @cached_method
        def basis(self, d=None):
            r"""
            Return the basis for (the ``d``-th homogeneous component
            of) ``self``.

            INPUT:

            - ``d`` -- (default: ``None``) nonnegative integer
              or ``None``

            OUTPUT:

            If ``d`` is ``None``, returns the basis of the module.
            Otherwise, returns the basis of the homogeneous component
            of degree ``d`` (i.e., the subfamily of the basis of the
            whole module which consists only of the basis vectors
            lying in `F_d \setminus \bigcup_{i<d} F_i`).

            The basis is always returned as a family.

            EXAMPLES::

                sage: A = ModulesWithBasis(ZZ).Filtered().example()
                sage: A.basis(4)
                Lazy family (Term map from Partitions to An example of a
                 filtered module with basis: the free module on partitions
                 over Integer Ring(i))_{i in Partitions of the integer 4}

            Without arguments, the full basis is returned::

                sage: A.basis()
                Lazy family (Term map from Partitions to An example of a
                 filtered module with basis: the free module on partitions
                 over Integer Ring(i))_{i in Partitions}
                sage: A.basis()
                Lazy family (Term map from Partitions to An example of a
                 filtered module with basis: the free module on partitions
                 over Integer Ring(i))_{i in Partitions}

            Checking this method on a filtered algebra. Note that this
            will typically raise a :exc:`NotImplementedError` when this
            feature is not implemented. ::

                sage: A = AlgebrasWithBasis(ZZ).Filtered().example()
                sage: A.basis(4)
                Traceback (most recent call last):
                ...
                NotImplementedError: infinite set

            Without arguments, the full basis is returned::

                sage: A.basis()
                Lazy family (Term map from Free abelian monoid indexed by
                 {'x', 'y', 'z'} to An example of a filtered algebra with
                 basis: the universal enveloping algebra of Lie algebra
                 of RR^3 with cross product over Integer Ring(i))_{i in
                 Free abelian monoid indexed by {'x', 'y', 'z'}}

            An example with a graded algebra::

                sage: E.<x,y> = ExteriorAlgebra(QQ)
                sage: E.basis()
                Lazy family (Term map from Subsets of {0,1} to
                 The exterior algebra of rank 2 over Rational Field(i))_{i in
                 Subsets of {0,1}}
            """
            if d is None:
                from sage.sets.family import Family
                return Family(self._indices, self.monomial)
            else:
                return self.homogeneous_component_basis(d)

        # TODO: Change `list(self._indices)` to `self._indices` and move
        #   this fallback to the category of finite-dimensional filtered
        #   modules with basis when it is implemented and the MRO issues
        #   are fixed (see trac #19397)
        def homogeneous_component_basis(self, d):
            """
            Return a basis for the ``d``-th homogeneous component of ``self``.

            EXAMPLES::

                sage: A = GradedModulesWithBasis(ZZ).example()
                sage: A.homogeneous_component_basis(4)
                Lazy family (Term map
                 from Partitions
                   to An example of a graded module with basis: the free module
                      on partitions over Integer Ring(i))_{i in Partitions of the integer 4}

                sage: # needs sage.modules
                sage: cat = GradedModulesWithBasis(ZZ)
                sage: C = CombinatorialFreeModule(ZZ, ['a', 'b'], category=cat)
                sage: C.degree_on_basis = lambda x: 1 if x == 'a' else 2
                sage: C.homogeneous_component_basis(1)
                Finite family {'a': B['a']}
                sage: C.homogeneous_component_basis(2)
                Finite family {'b': B['b']}
            """
            from sage.sets.family import Family
            try:
                S = self._indices.subset(size=d)
            except (AttributeError, ValueError, TypeError):
                S = [i for i in list(self._indices) if self.degree_on_basis(i) == d]
            return Family(S, self.monomial)

        def homogeneous_component(self, d):
            """
            Return the ``d``-th homogeneous component of ``self``.

            EXAMPLES::

                sage: A = GradedModulesWithBasis(ZZ).example()
                sage: A.homogeneous_component(4)
                Degree 4 homogeneous component of An example of a graded module
                 with basis: the free module on partitions over Integer Ring
            """
            from sage.categories.modules_with_basis import ModulesWithBasis
            from sage.categories.filtered_algebras import FilteredAlgebras
            if self.base_ring() in FilteredAlgebras:
                raise NotImplementedError("this is only a natural module over"
                                          " the degree 0 component of the filtered"
                                          " algebra and coordinate rings are not"
                                          " yet implemented for submodules")
            category = ModulesWithBasis(self.category().base_ring())
            M = self.submodule(self.homogeneous_component_basis(d),
                               category=category,
                               already_echelonized=True)
            M.rename("Degree {} homogeneous component of {}".format(d, self))
            return M

        def hilbert_series(self, prec=None):
            r"""
            Return the Hilbert series of ``self``.

            Let `R` be a commutative ring (with unit). Let
            `M = \bigcup_{n=0}^{\infty} M_n` be a filtered `R`-module.
            The *Hilbert series* of `M` is the formal power series

            .. MATH::

                H(t) = \sum_{n=0}^{\infty} \ell(M_n / M_{n-1}) t^n,

            where `\ell(N)` is the *length* of `N`, which is the
            longest chain of submodules (over `R`), and by convention
            `M_{-1} = \{0\}`. By the assumptions of the category,
            `M_n / M_{n-1}` is a free `R`-module, and so `\ell(M_n / M_{n-1})`
            is equal to the rank of `M_n / M_{n-1}`.

            INPUT:

            - ``prec`` -- (default: `\infty`) the precision

            OUTPUT:

            If the precision is finite, then this returns an element in the
            :class:`PowerSeriesRing` over ``ZZ``. Otherwise it returns an
            element in the :class:`LazyPowerSeriesRing` over ``ZZ``.

            EXAMPLES::

                sage: A = GradedModulesWithBasis(ZZ).example()
                sage: A.hilbert_series()
                1 + t + 2*t^2 + 3*t^3 + 5*t^4 + 7*t^5 + 11*t^6 + O(t^7)
                sage: A.hilbert_series(10)
                1 + t + 2*t^2 + 3*t^3 + 5*t^4 + 7*t^5 + 11*t^6 + 15*t^7 + 22*t^8 + 30*t^9 + O(t^10)
            """
            from sage.rings.integer_ring import ZZ
            if prec is None:
                from sage.rings.lazy_series_ring import LazyPowerSeriesRing
                R = LazyPowerSeriesRing(ZZ, 't')
                return R(lambda n: self.homogeneous_component_basis(n).cardinality())
            from sage.rings.power_series_ring import PowerSeriesRing
            R = PowerSeriesRing(ZZ, 't')
            elt = R([self.homogeneous_component_basis(n).cardinality() for n in range(prec)])
            return elt.O(prec)

        def graded_algebra(self):
            r"""
            Return the associated graded module to ``self``.

            See :class:`~sage.algebras.associated_graded.AssociatedGradedAlgebra`
            for the definition and the properties of this.

            If the filtered module ``self`` with basis is called `A`,
            then this method returns `\operatorname{gr} A`. The method
            :meth:`to_graded_conversion` returns the canonical
            `R`-module isomorphism `A \to \operatorname{gr} A` induced
            by the basis of `A`, and the method
            :meth:`from_graded_conversion` returns the inverse of this
            isomorphism. The method :meth:`projection` projects
            elements of `A` onto `\operatorname{gr} A` according to
            their place in the filtration on `A`.

            .. WARNING::

                When not overridden, this method returns the default
                implementation of an associated graded module --
                namely, ``AssociatedGradedAlgebra(self)``, where
                ``AssociatedGradedAlgebra`` is
                :class:`~sage.algebras.associated_graded.AssociatedGradedAlgebra`.
                But some instances of :class:`FilteredModulesWithBasis`
                override this method, as the associated graded module
                often is (isomorphic) to a simpler object (for instance,
                the associated graded module of a graded module can be
                identified with the graded module itself). Generic code
                that uses associated graded modules (such as the code
                of the :meth:`induced_graded_map` method below) should
                make sure to only communicate with them via the
                :meth:`to_graded_conversion`,
                :meth:`from_graded_conversion` and
                :meth:`projection` methods (in particular,
                do not expect there to be a conversion from ``self``
                to ``self.graded_algebra()``; this currently does not
                work for Clifford algebras). Similarly, when
                overriding :meth:`graded_algebra`, make sure to
                accordingly redefine these three methods, unless their
                definitions below still apply to your case (this will
                happen whenever the basis of your :meth:`graded_algebra`
                has the same indexing set as ``self``, and the partition
                of this indexing set according to degree is the same as
                for ``self``).

            EXAMPLES::

                sage: A = ModulesWithBasis(ZZ).Filtered().example()
                sage: A.graded_algebra()
                Graded Module of An example of a filtered module with basis:
                 the free module on partitions over Integer Ring
            """
            from sage.algebras.associated_graded import AssociatedGradedAlgebra
            return AssociatedGradedAlgebra(self)

        # Maps

        def to_graded_conversion(self):
            r"""
            Return the canonical `R`-module isomorphism
            `A \to \operatorname{gr} A` induced by the basis of `A`
            (where `A = ` ``self``).

            This is an isomorphism of `R`-modules. See
            the class documentation :class:`AssociatedGradedAlgebra`.

            .. SEEALSO::

                :meth:`from_graded_conversion`

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: A = Modules(QQ).WithBasis().Filtered().example()
                sage: p = -2 * A.an_element(); p
                -4*P[] - 4*P[1] - 6*P[2]
                sage: q = A.to_graded_conversion()(p); q
                -4*Bbar[[]] - 4*Bbar[[1]] - 6*Bbar[[2]]
                sage: q.parent() is A.graded_algebra()
                True
            """
            base_one = self.base_ring().one()
            return self.module_morphism(diagonal=lambda x: base_one,
                                        codomain=self.graded_algebra())

        def from_graded_conversion(self):
            r"""
            Return the inverse of the canonical `R`-module isomorphism
            `A \to \operatorname{gr} A` induced by the basis of `A`
            (where `A = ` ``self``). This inverse is an isomorphism
            `\operatorname{gr} A \to A`.

            This is an isomorphism of `R`-modules. See
            the class documentation :class:`AssociatedGradedAlgebra`.

            .. SEEALSO::

                :meth:`to_graded_conversion`

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: A = Modules(QQ).WithBasis().Filtered().example()
                sage: p = -2 * A.an_element(); p
                -4*P[] - 4*P[1] - 6*P[2]
                sage: q = A.to_graded_conversion()(p); q
                -4*Bbar[[]] - 4*Bbar[[1]] - 6*Bbar[[2]]
                sage: A.from_graded_conversion()(q) == p
                True
                sage: q.parent() is A.graded_algebra()
                True
            """
            base_one = self.base_ring().one()
            return self.graded_algebra().module_morphism(diagonal=lambda x: base_one,
                                                         codomain=self)

        def projection(self, i):
            r"""
            Return the `i`-th projection `p_i : F_i \to G_i` (in the
            notations of the class documentation
            :class:`AssociatedGradedAlgebra`, where `A = ` ``self``).

            This method actually does not return the map `p_i` itself,
            but an extension of `p_i` to the whole `R`-module `A`.
            This extension is the composition of the `R`-module
            isomorphism `A \to \operatorname{gr} A` with the canonical
            projection of the graded `R`-module `\operatorname{gr} A`
            onto its `i`-th graded component `G_i`. The codomain of
            this map is `\operatorname{gr} A`, although its actual
            image is `G_i`. The map `p_i` is obtained from this map
            by restricting its domain to `F_i` and its image to `G_i`.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: A = Modules(ZZ).WithBasis().Filtered().example()
                sage: p = -2 * A.an_element(); p
                -4*P[] - 4*P[1] - 6*P[2]
                sage: q = A.projection(2)(p); q
                -6*Bbar[[2]]
                sage: q.parent() is A.graded_algebra()
                True
                sage: A.projection(3)(p)
                0
            """
            base_zero = self.base_ring().zero()
            base_one = self.base_ring().one()
            grA = self.graded_algebra()
            proj = lambda x: (base_one if self.degree_on_basis(x) == i
                              else base_zero)
            return self.module_morphism(diagonal=proj, codomain=grA)

        def induced_graded_map(self, other, f):
            r"""
            Return the graded linear map between the associated graded
            modules of ``self`` and ``other`` canonically induced by
            the filtration-preserving map ``f : self -> other``.

            Let `A` and `B` be two filtered modules with basis, and let
            `(F_i)_{i \in I}` and `(G_i)_{i \in I}` be their
            filtrations. Let `f : A \to B` be a linear map which
            preserves the filtration (i.e., satisfies `f(F_i) \subseteq
            G_i` for all `i \in I`). Then, there is a canonically
            defined graded linear map
            `\operatorname{gr} f : \operatorname{gr} A \to
            \operatorname{gr} B` which satisfies

            .. MATH::

                (\operatorname{gr} f) (p_i(a)) = p_i(f(a))
                \qquad \text{for all } i \in I \text{ and } a \in F_i ,

            where the `p_i` on the left hand side is the canonical
            projection from `F_i` onto the `i`-th graded component
            of `\operatorname{gr} A`, while the `p_i` on the right
            hand side is the canonical projection from `G_i` onto
            the `i`-th graded component of `\operatorname{gr} B`.

            INPUT:

            - ``other`` -- a filtered algebra with basis

            - ``f`` -- a filtration-preserving linear map from ``self``
              to ``other`` (can be given as a morphism or as a function)

            OUTPUT: the graded linear map `\operatorname{gr} f`

            EXAMPLES:

            **Example 1.**

            We start with the free `\QQ`-module with basis the set of all
            partitions::

                sage: A = Modules(QQ).WithBasis().Filtered().example(); A               # needs sage.combinat sage.modules
                An example of a filtered module with basis: the free module
                 on partitions over Rational Field
                sage: M = A.indices(); M                                                # needs sage.combinat sage.modules
                Partitions
                sage: p1, p2, p21, p321 = [A.basis()[Partition(i)]                      # needs sage.combinat sage.modules
                ....:                      for i in [[1], [2], [2,1], [3,2,1]]]

            Let us define a map from ``A`` to itself which acts on the
            basis by sending every partition `\lambda` to the sum of
            the conjugates of all partitions `\mu` for which
            `\lambda / \mu` is a horizontal strip::

                sage: # needs sage.combinat sage.modules
                sage: def map_on_basis(lam):
                ....:     def mus(k):
                ....:         return lam.remove_horizontal_border_strip(k)
                ....:     return A.sum_of_monomials([Partition(mu).conjugate()
                ....:                                for k in range(sum(lam) + 1)
                ....:                                for mu in mus(k)])
                sage: f = A.module_morphism(on_basis=map_on_basis,
                ....:                       codomain=A)
                sage: f(p1)
                P[] + P[1]
                sage: f(p2)
                P[] + P[1] + P[1, 1]
                sage: f(p21)
                P[1] + P[1, 1] + P[2] + P[2, 1]
                sage: f(p21 - p1)
                -P[] + P[1, 1] + P[2] + P[2, 1]
                sage: f(p321)
                P[2, 1] + P[2, 1, 1] + P[2, 2] + P[2, 2, 1]
                 + P[3, 1] + P[3, 1, 1] + P[3, 2] + P[3, 2, 1]

            We now compute `\operatorname{gr} f` ::

                sage: # needs sage.combinat sage.modules
                sage: grA = A.graded_algebra(); grA
                Graded Module of An example of a filtered module with basis:
                 the free module on partitions over Rational Field
                sage: pp1, pp2, pp21, pp321 = [A.to_graded_conversion()(i)
                ....:                          for i in [p1, p2, p21, p321]]
                sage: pp2 + 4 * pp21
                Bbar[[2]] + 4*Bbar[[2, 1]]
                sage: grf = A.induced_graded_map(A, f); grf
                Generic endomorphism of Graded Module of
                 An example of a filtered module with basis:
                  the free module on partitions over Rational Field
                sage: grf(pp1)
                Bbar[[1]]
                sage: grf(pp2 + 4 * pp21)
                Bbar[[1, 1]] + 4*Bbar[[2, 1]]

            **Example 2.**

            We shall now construct `\operatorname{gr} f` for a
            different map `f` out of the same ``A``; the new map
            `f` will lead into a graded algebra already, namely into
            the algebra of symmetric functions::

                sage: # needs sage.combinat sage.modules
                sage: h = SymmetricFunctions(QQ).h()
                sage: def map_on_basis(lam):  # redefining map_on_basis
                ....:     def mus(k):
                ....:         return lam.remove_horizontal_border_strip(k)
                ....:     return h.sum_of_monomials([Partition(mu).conjugate()
                ....:                                for k in range(sum(lam) + 1)
                ....:                                for mu in mus(k)])
                sage: f = A.module_morphism(on_basis=map_on_basis,
                ....:                       codomain=h)  # redefining f
                sage: f(p1)
                h[] + h[1]
                sage: f(p2)
                h[] + h[1] + h[1, 1]
                sage: f(A.zero())
                0
                sage: f(p2 - 3*p1)
                -2*h[] - 2*h[1] + h[1, 1]

            The algebra ``h`` of symmetric functions in the `h`-basis
            is already graded, so its associated graded algebra is
            implemented as itself::

                sage: # needs sage.combinat sage.modules
                sage: grh = h.graded_algebra(); grh is h
                True
                sage: grf = A.induced_graded_map(h, f); grf
                Generic morphism:
                  From: Graded Module of An example of a filtered
                   module with basis: the free module on partitions
                   over Rational Field
                  To:   Symmetric Functions over Rational Field
                   in the homogeneous basis
                sage: grf(pp1)
                h[1]
                sage: grf(pp2)
                h[1, 1]
                sage: grf(pp321)
                h[3, 2, 1]
                sage: grf(pp2 - 3*pp1)
                -3*h[1] + h[1, 1]
                sage: grf(pp21)
                h[2, 1]
                sage: grf(grA.zero())
                0

            **Example 3.**

            After having had a graded module as the codomain, let us try to
            have one as the domain instead. Our new ``f`` will go from ``h``
            to ``A``::

                sage: # needs sage.combinat sage.modules
                sage: def map_on_basis(lam):  # redefining map_on_basis
                ....:     def mus(k):
                ....:         return lam.remove_horizontal_border_strip(k)
                ....:     return A.sum_of_monomials([Partition(mu).conjugate()
                ....:                                for k in range(sum(lam) + 1)
                ....:                                for mu in mus(k)])
                sage: f = h.module_morphism(on_basis=map_on_basis,
                ....:                       codomain=A)  # redefining f
                sage: f(h[1])
                P[] + P[1]
                sage: f(h[2])
                P[] + P[1] + P[1, 1]
                sage: f(h[1, 1])
                P[1] + P[2]
                sage: f(h[2, 2])
                P[1, 1] + P[2, 1] + P[2, 2]
                sage: f(h[3, 2, 1])
                P[2, 1] + P[2, 1, 1] + P[2, 2] + P[2, 2, 1]
                 + P[3, 1] + P[3, 1, 1] + P[3, 2] + P[3, 2, 1]
                sage: f(h.one())
                P[]
                sage: grf = h.induced_graded_map(A, f); grf
                Generic morphism:
                  From: Symmetric Functions over Rational Field
                   in the homogeneous basis
                  To:   Graded Module of An example of a filtered
                   module with basis: the free module on partitions
                   over Rational Field
                sage: grf(h[1])
                Bbar[[1]]
                sage: grf(h[2])
                Bbar[[1, 1]]
                sage: grf(h[1, 1])
                Bbar[[2]]
                sage: grf(h[2, 2])
                Bbar[[2, 2]]
                sage: grf(h[3, 2, 1])
                Bbar[[3, 2, 1]]
                sage: grf(h.one())
                Bbar[[]]

            **Example 4.**

            The construct `\operatorname{gr} f` also makes sense when `f`
            is a filtration-preserving map between graded modules. ::

                sage: # needs sage.combinat sage.modules
                sage: def map_on_basis(lam):  # redefining map_on_basis
                ....:     def mus(k):
                ....:         return lam.remove_horizontal_border_strip(k)
                ....:     return h.sum_of_monomials([Partition(mu).conjugate()
                ....:                                for k in range(sum(lam) + 1)
                ....:                                for mu in mus(k)])
                sage: f = h.module_morphism(on_basis=map_on_basis,
                ....:                       codomain=h)  # redefining f
                sage: f(h[1])
                h[] + h[1]
                sage: f(h[2])
                h[] + h[1] + h[1, 1]
                sage: f(h[1, 1])
                h[1] + h[2]
                sage: f(h[2, 1])
                h[1] + h[1, 1] + h[2] + h[2, 1]
                sage: f(h.one())
                h[]
                sage: grf = h.induced_graded_map(h, f); grf
                Generic endomorphism of
                 Symmetric Functions over Rational Field in the homogeneous basis
                sage: grf(h[1])
                h[1]
                sage: grf(h[2])
                h[1, 1]
                sage: grf(h[1, 1])
                h[2]
                sage: grf(h[2, 1])
                h[2, 1]
                sage: grf(h.one())
                h[]
            """
            grA = self.graded_algebra()
            grB = other.graded_algebra()
            from sage.categories.graded_modules_with_basis import GradedModulesWithBasis
            cat = GradedModulesWithBasis(self.base_ring())
            from_gr = self.from_graded_conversion()

            def on_basis(m):
                i = grA.degree_on_basis(m)
                lifted_img_of_m = f(from_gr(grA.monomial(m)))
                return other.projection(i)(lifted_img_of_m)
            return grA.module_morphism(on_basis=on_basis,
                                       codomain=grB, category=cat)
            # If we could assume that the projection of the basis
            # element of ``self`` indexed by an index ``m`` is the
            # basis element of ``grA`` indexed by ``m``, then this
            # could go faster:
            #
            # def on_basis(m):
            #     i = grA.degree_on_basis(m)
            #     return grB.projection(i)(f(self.monomial(m)))
            # return grA.module_morphism(on_basis=on_basis,
            #                            codomain=grB, category=cat)
            #
            # But this assumption might come back to bite us in the
            # ass one day. What do you think?

    class ElementMethods:

        def is_homogeneous(self):
            r"""
            Return whether the element ``self`` is homogeneous.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: A = ModulesWithBasis(ZZ).Filtered().example()
                sage: x = A(Partition((3,2,1)))
                sage: y = A(Partition((4,4,1)))
                sage: z = A(Partition((2,2,2)))
                sage: (3*x).is_homogeneous()
                True
                sage: (x - y).is_homogeneous()
                False
                sage: (x+2*z).is_homogeneous()
                True

            Here is an example with a graded algebra::

                sage: # needs sage.combinat sage.modules
                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: (x, y) = (S[2], S[3])
                sage: (3*x).is_homogeneous()
                True
                sage: (x^3 - y^2).is_homogeneous()
                True
                sage: ((x + y)^2).is_homogeneous()
                False

            Let us now test a filtered algebra (but remember that the
            notion of homogeneity now depends on the choice of a
            basis, or at least on a definition of homogeneous
            components)::

                sage: A = AlgebrasWithBasis(QQ).Filtered().example()
                sage: x,y,z = A.algebra_generators()
                sage: (x*y).is_homogeneous()
                True
                sage: (y*x).is_homogeneous()
                False
                sage: A.one().is_homogeneous()
                True
                sage: A.zero().is_homogeneous()
                True
                sage: (A.one()+x).is_homogeneous()
                False
            """
            degree_on_basis = self.parent().degree_on_basis
            degree = None
            for m in self.support():
                if degree is None:
                    degree = degree_on_basis(m)
                else:
                    if degree != degree_on_basis(m):
                        return False
            return True

        @abstract_method(optional=True)
        def degree_on_basis(self, m):
            r"""
            Return the degree of the basis element indexed by ``m``
            in ``self``.

            EXAMPLES::

                sage: A = GradedModulesWithBasis(QQ).example()                          # needs sage.combinat sage.modules
                sage: A.degree_on_basis(Partition((2,1)))                               # needs sage.combinat sage.modules
                3
                sage: A.degree_on_basis(Partition((4,2,1,1,1,1)))                       # needs sage.combinat sage.modules
                10
            """

        def homogeneous_degree(self):
            r"""
            The degree of a nonzero homogeneous element ``self`` in the
            filtered module.

            .. NOTE::

                This raises an error if the element is not homogeneous.
                To compute the maximum of the degrees of the homogeneous
                summands of a (not necessarily homogeneous) element, use
                :meth:`maximal_degree` instead.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: A = ModulesWithBasis(ZZ).Filtered().example()
                sage: x = A(Partition((3,2,1)))
                sage: y = A(Partition((4,4,1)))
                sage: z = A(Partition((2,2,2)))
                sage: x.degree()
                6
                sage: (x + 2*z).degree()
                6
                sage: (y - x).degree()
                Traceback (most recent call last):
                ...
                ValueError: element is not homogeneous

            An example in a graded algebra::

                sage: # needs sage.combinat sage.modules
                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: (x, y) = (S[2], S[3])
                sage: x.homogeneous_degree()
                2
                sage: (x^3 + 4*y^2).homogeneous_degree()
                6
                sage: ((1 + x)^3).homogeneous_degree()
                Traceback (most recent call last):
                ...
                ValueError: element is not homogeneous

            Let us now test a filtered algebra (but remember that the
            notion of homogeneity now depends on the choice of a
            basis)::

                sage: A = AlgebrasWithBasis(QQ).Filtered().example()
                sage: x,y,z = A.algebra_generators()
                sage: (x*y).homogeneous_degree()
                2
                sage: (y*x).homogeneous_degree()
                Traceback (most recent call last):
                ...
                ValueError: element is not homogeneous
                sage: A.one().homogeneous_degree()
                0

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()                      # needs sage.combinat sage.modules
                sage: S.zero().degree()                                                 # needs sage.combinat sage.modules
                Traceback (most recent call last):
                ...
                ValueError: the zero element does not have a well-defined degree
            """
            if not self.support():
                raise ValueError("the zero element does not have a well-defined degree")
            if not self.is_homogeneous():
                raise ValueError("element is not homogeneous")
            return self.parent().degree_on_basis(self.leading_support())

        # default choice for degree; will be overridden as necessary
        degree = homogeneous_degree

        def maximal_degree(self):
            """
            The maximum of the degrees of the homogeneous components
            of ``self``.

            This is also the smallest `i` such that ``self`` belongs
            to `F_i`. Hence, it does not depend on the basis of the
            parent of ``self``.

            .. SEEALSO:: :meth:`homogeneous_degree`

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: A = ModulesWithBasis(ZZ).Filtered().example()
                sage: x = A(Partition((3,2,1)))
                sage: y = A(Partition((4,4,1)))
                sage: z = A(Partition((2,2,2)))
                sage: x.maximal_degree()
                6
                sage: (x + 2*z).maximal_degree()
                6
                sage: (y - x).maximal_degree()
                9
                sage: (3*z).maximal_degree()
                6

            Now, we test this on a graded algebra::

                sage: # needs sage.combinat sage.modules
                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: (x, y) = (S[2], S[3])
                sage: x.maximal_degree()
                2
                sage: (x^3 + 4*y^2).maximal_degree()
                6
                sage: ((1 + x)^3).maximal_degree()
                6

            Let us now test a filtered algebra::

                sage: # needs sage.combinat sage.modules
                sage: A = AlgebrasWithBasis(QQ).Filtered().example()
                sage: x,y,z = A.algebra_generators()
                sage: (x*y).maximal_degree()
                2
                sage: (y*x).maximal_degree()
                2
                sage: A.one().maximal_degree()
                0
                sage: A.zero().maximal_degree()
                Traceback (most recent call last):
                ...
                ValueError: the zero element does not have a well-defined degree
                sage: (A.one()+x).maximal_degree()
                1

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()                      # needs sage.combinat sage.modules
                sage: S.zero().degree()                                                 # needs sage.combinat sage.modules
                Traceback (most recent call last):
                ...
                ValueError: the zero element does not have a well-defined degree
            """
            if self.is_zero():
                raise ValueError("the zero element does not have a well-defined degree")
            degree_on_basis = self.parent().degree_on_basis
            return max(degree_on_basis(m) for m in self.support())

        def homogeneous_component(self, n):
            r"""
            Return the homogeneous component of degree ``n`` of the
            element ``self``.

            Let `m` be an element of a filtered `R`-module `M` with
            basis. Then, `m` can be uniquely written in the form
            `m = \sum_{i \in I} m_i`, where each `m_i` is a
            homogeneous element of degree `i`. For `n \in I`, we
            define the homogeneous component of degree `n` of the
            element `m` to be `m_n`.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: A = ModulesWithBasis(ZZ).Filtered().example()
                sage: x = A.an_element(); x
                2*P[] + 2*P[1] + 3*P[2]
                sage: x.homogeneous_component(-1)
                0
                sage: x.homogeneous_component(0)
                2*P[]
                sage: x.homogeneous_component(1)
                2*P[1]
                sage: x.homogeneous_component(2)
                3*P[2]
                sage: x.homogeneous_component(3)
                0

                sage: # needs sage.combinat sage.modules
                sage: A = ModulesWithBasis(ZZ).Graded().example()
                sage: x = A.an_element(); x
                2*P[] + 2*P[1] + 3*P[2]
                sage: x.homogeneous_component(-1)
                0
                sage: x.homogeneous_component(0)
                2*P[]
                sage: x.homogeneous_component(1)
                2*P[1]
                sage: x.homogeneous_component(2)
                3*P[2]
                sage: x.homogeneous_component(3)
                0

                sage: A = AlgebrasWithBasis(ZZ).Filtered().example()
                sage: G = A.algebra_generators()
                sage: g = A.an_element() - 2 * G['x'] * G['y']; g
                U['x']^2*U['y']^2*U['z']^3 - 2*U['x']*U['y']
                 + 2*U['x'] + 3*U['y'] + 1
                sage: g.homogeneous_component(-1)
                0
                sage: g.homogeneous_component(0)
                1
                sage: g.homogeneous_component(2)
                -2*U['x']*U['y']
                sage: g.homogeneous_component(5)
                0
                sage: g.homogeneous_component(7)
                U['x']^2*U['y']^2*U['z']^3
                sage: g.homogeneous_component(8)
                0

            TESTS:

            Check that this really returns ``A.zero()`` and not a plain ``0``::

                sage: # needs sage.combinat sage.modules
                sage: A = ModulesWithBasis(ZZ).Filtered().example()
                sage: x = A.an_element()
                sage: x.homogeneous_component(3).parent() is A
                True
            """
            degree_on_basis = self.parent().degree_on_basis
            return self.parent().sum_of_terms((i, c)
                                              for (i, c) in self
                                              if degree_on_basis(i) == n)

        def truncate(self, n):
            """
            Return the sum of the homogeneous components of degree
            strictly less than ``n`` of ``self``.

            See :meth:`homogeneous_component` for the notion of a
            homogeneous component.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: A = ModulesWithBasis(ZZ).Filtered().example()
                sage: x = A.an_element(); x
                2*P[] + 2*P[1] + 3*P[2]
                sage: x.truncate(0)
                0
                sage: x.truncate(1)
                2*P[]
                sage: x.truncate(2)
                2*P[] + 2*P[1]
                sage: x.truncate(3)
                2*P[] + 2*P[1] + 3*P[2]

                sage: # needs sage.combinat sage.modules
                sage: A = ModulesWithBasis(ZZ).Graded().example()
                sage: x = A.an_element(); x
                2*P[] + 2*P[1] + 3*P[2]
                sage: x.truncate(0)
                0
                sage: x.truncate(1)
                2*P[]
                sage: x.truncate(2)
                2*P[] + 2*P[1]
                sage: x.truncate(3)
                2*P[] + 2*P[1] + 3*P[2]

                sage: A = AlgebrasWithBasis(ZZ).Filtered().example()
                sage: G = A.algebra_generators()
                sage: g = A.an_element() - 2 * G['x'] * G['y']; g
                U['x']^2*U['y']^2*U['z']^3 - 2*U['x']*U['y']
                 + 2*U['x'] + 3*U['y'] + 1
                sage: g.truncate(-1)
                0
                sage: g.truncate(0)
                0
                sage: g.truncate(2)
                2*U['x'] + 3*U['y'] + 1
                sage: g.truncate(3)
                -2*U['x']*U['y'] + 2*U['x'] + 3*U['y'] + 1
                sage: g.truncate(5)
                -2*U['x']*U['y'] + 2*U['x'] + 3*U['y'] + 1
                sage: g.truncate(7)
                -2*U['x']*U['y'] + 2*U['x'] + 3*U['y'] + 1
                sage: g.truncate(8)
                U['x']^2*U['y']^2*U['z']^3 - 2*U['x']*U['y']
                 + 2*U['x'] + 3*U['y'] + 1

            TESTS:

            Check that this really return ``A.zero()`` and not a plain ``0``::

                sage: # needs sage.combinat sage.modules
                sage: A = ModulesWithBasis(ZZ).Filtered().example()
                sage: x = A.an_element()
                sage: x.truncate(0).parent() is A
                True
            """
            degree_on_basis = self.parent().degree_on_basis
            return self.parent().sum_of_terms((i, c) for (i, c) in self
                                              if degree_on_basis(i) < n)

    class Subobjects(SubobjectsCategory):
        class ParentMethods:
            def degree_on_basis(self, m):
                r"""
                Return the degree of the basis element indexed by ``m``
                in ``self``.

                EXAMPLES::

                    sage: # needs sage.modules
                    sage: E.<x,y,z> = ExteriorAlgebra(QQ)
                    sage: S = E.submodule([x + y, x*y - y*z, y])
                    sage: B = S.basis()
                    sage: [B[0].lift(), B[1].lift(), B[2].lift()]
                    [x, y, x*y - y*z]
                    sage: S.degree_on_basis(0)
                    1
                    sage: S.degree_on_basis(1)
                    1
                    sage: S.degree_on_basis(2)
                    2
                """
                return self.basis()[m].lift().degree()

        class ElementMethods:
            def degree(self):
                r"""
                Return the degree of ``self``.

                EXAMPLES::

                    sage: # needs sage.modules
                    sage: E.<x,y,z> = ExteriorAlgebra(QQ)
                    sage: S = E.submodule([x + y, x*y - y*z, y])
                    sage: B = S.basis()
                    sage: [B[0].lift(), B[1].lift(), B[2].lift()]
                    [x, y, x*y - y*z]
                    sage: B[0].degree()
                    1
                    sage: B[1].degree()
                    1
                    sage: (B[0] + 3*B[1]).degree()
                    1

                The degree of inhomogeneous elements is not defined
                (following the behavior of the exterior algebra)::

                    sage: (B[0] + B[2]).degree()                                        # needs sage.modules
                    Traceback (most recent call last):
                    ...
                    ValueError: element is not homogeneous

                We can still get the maximal degree::

                    sage: (B[0] + B[2]).maximal_degree()                                # needs sage.modules
                    2
                """
                return self.lift().degree()

            def maximal_degree(self):
                r"""
                The maximum of the degrees of the homogeneous components
                of ``self``.

                This is also the smallest `i` such that ``self`` belongs
                to `F_i`. Hence, it does not depend on the basis of the
                parent of ``self``.

                .. SEEALSO:: :meth:`homogeneous_degree`

                EXAMPLES::

                    sage: # needs sage.modules
                    sage: E.<x,y,z> = ExteriorAlgebra(QQ)
                    sage: F = E.submodule([x + 1, x*y - 1])
                    sage: B = F.basis()
                    sage: [B[0].lift(), B[1].lift()]
                    [-x*y + 1, x*y + x]
                    sage: B[0].maximal_degree()
                    2
                    sage: B[1].maximal_degree()
                    2
                """
                return self.lift().maximal_degree()

    class FiniteDimensional(CategoryWithAxiom_over_base_ring):
        class ParentMethods:
            def hilbert_series(self, prec=None):
                r"""
                Return the Hilbert series of ``self`` as a polynomial.

                Let `R` be a commutative ring (with unit). Let
                `M = \bigcup_{n=0}^{\infty} M_n` be a filtered `R`-module.
                The *Hilbert series* of `M` is the formal power series

                .. MATH::

                    H(t) = \sum_{n=0}^{\infty} \ell(M_n / M_{n-1}) t^n,

                where `\ell(N)` is the *length* of `N`, which is the
                longest chain of submodules (over `R`), and by convention
                `M_{-1} = \{0\}`. By the assumptions of the category,
                `M_n / M_{n-1}` is a free `R`-module, and so
                `\ell(M_n / M_{n-1})` is equal to the rank of `M_n / M_{n-1}`.

                EXAMPLES::

                    sage: OS = hyperplane_arrangements.braid(3).orlik_solomon_algebra(QQ)
                    sage: OS.hilbert_series()
                    2*t^2 + 3*t + 1

                    sage: OS = matroids.Uniform(5, 3).orlik_solomon_algebra(ZZ)
                    sage: OS.hilbert_series()
                    t^3 + 3*t^2 + 3*t + 1

                    sage: OS = matroids.PG(2, 3).orlik_solomon_algebra(ZZ['x','y'])
                    sage: OS.hilbert_series()
                    27*t^3 + 39*t^2 + 13*t + 1
                """
                from collections import defaultdict
                from sage.rings.integer_ring import ZZ
                from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                PR = PolynomialRing(ZZ, 't')
                dims = defaultdict(ZZ)
                for b in self.basis():
                    dims[b.homogeneous_degree()] += 1
                return PR(dims)
