# sage_setup: distribution = sagemath-categories
r"""
Filtered Algebras With Basis

A filtered algebra with basis over a commutative ring `R`
is a filtered algebra over `R` endowed with the structure
of a filtered module with basis (with the same underlying
filtered-module structure). See
:class:`~sage.categories.filtered_algebras.FilteredAlgebras` and
:class:`~sage.categories.filtered_modules_with_basis.FilteredModulesWithBasis`
for these two notions.
"""
#*****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.filtered_modules import FilteredModulesCategory


class FilteredAlgebrasWithBasis(FilteredModulesCategory):
    """
    The category of filtered algebras with a distinguished
    homogeneous basis.

    A filtered algebra with basis over a commutative ring `R`
    is a filtered algebra over `R` endowed with the structure
    of a filtered module with basis (with the same underlying
    filtered-module structure). See
    :class:`~sage.categories.filtered_algebras.FilteredAlgebras` and
    :class:`~sage.categories.filtered_modules_with_basis.FilteredModulesWithBasis`
    for these two notions.

    EXAMPLES::

        sage: C = AlgebrasWithBasis(ZZ).Filtered(); C
        Category of filtered algebras with basis over Integer Ring
        sage: sorted(C.super_categories(), key=str)
        [Category of algebras with basis over Integer Ring,
         Category of filtered algebras over Integer Ring,
         Category of filtered modules with basis over Integer Ring]

    TESTS::

        sage: TestSuite(C).run()
    """
    class ParentMethods:
        def graded_algebra(self):
            r"""
            Return the associated graded algebra to ``self``.

            See :class:`~sage.algebras.associated_graded.AssociatedGradedAlgebra`
            for the definition and the properties of this.

            If the filtered algebra ``self`` with basis is called `A`,
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
                implementation of an associated graded algebra --
                namely, ``AssociatedGradedAlgebra(self)``, where
                ``AssociatedGradedAlgebra`` is
                :class:`~sage.algebras.associated_graded.AssociatedGradedAlgebra`.
                But many instances of :class:`FilteredAlgebrasWithBasis`
                override this method, as the associated graded algebra
                often is (isomorphic) to a simpler object (for instance,
                the associated graded algebra of a graded algebra can be
                identified with the graded algebra itself). Generic code
                that uses associated graded algebras (such as the code
                of the :meth:`induced_graded_map` method below) should
                make sure to only communicate with them via the
                :meth:`to_graded_conversion`,
                :meth:`from_graded_conversion`, and
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

            .. TODO::

                Maybe the thing about the conversion from ``self``
                to ``self.graded_algebra()`` on the Clifford at least
                could be made to work? (I would still warn the user
                against ASSUMING that it must work -- as there is
                probably no way to guarantee it in all cases, and
                we shouldn't require users to mess with
                element constructors.)

            EXAMPLES::

                sage: A = AlgebrasWithBasis(ZZ).Filtered().example()
                sage: A.graded_algebra()
                Graded Algebra of An example of a filtered algebra with basis:
                 the universal enveloping algebra of
                 Lie algebra of RR^3 with cross product over Integer Ring
            """
            from sage.algebras.associated_graded import AssociatedGradedAlgebra
            return AssociatedGradedAlgebra(self)

        # Maps

        def to_graded_conversion(self):
            r"""
            Return the canonical `R`-module isomorphism
            `A \to \operatorname{gr} A` induced by the basis of `A`
            (where `A = ` ``self``).

            This is an isomorphism of `R`-modules, not of algebras. See
            the class documentation :class:`AssociatedGradedAlgebra`.

            .. SEEALSO::

                :meth:`from_graded_conversion`

            EXAMPLES::

                sage: A = Algebras(QQ).WithBasis().Filtered().example()
                sage: p = A.an_element() + A.algebra_generators()['x'] + 2; p
                U['x']^2*U['y']^2*U['z']^3 + 3*U['x'] + 3*U['y'] + 3
                sage: q = A.to_graded_conversion()(p); q
                bar(U['x']^2*U['y']^2*U['z']^3) + 3*bar(U['x'])
                 + 3*bar(U['y']) + 3*bar(1)
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

            This is an isomorphism of `R`-modules, not of algebras. See
            the class documentation :class:`AssociatedGradedAlgebra`.

            .. SEEALSO::

                :meth:`to_graded_conversion`

            EXAMPLES::

                sage: A = Algebras(QQ).WithBasis().Filtered().example()
                sage: p = A.an_element() + A.algebra_generators()['x'] + 2; p
                U['x']^2*U['y']^2*U['z']^3 + 3*U['x'] + 3*U['y'] + 3
                sage: q = A.to_graded_conversion()(p)
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

                sage: A = Algebras(QQ).WithBasis().Filtered().example()
                sage: p = A.an_element() + A.algebra_generators()['x'] + 2; p
                U['x']^2*U['y']^2*U['z']^3 + 3*U['x'] + 3*U['y'] + 3
                sage: q = A.projection(7)(p); q
                bar(U['x']^2*U['y']^2*U['z']^3)
                sage: q.parent() is A.graded_algebra()
                True
                sage: A.projection(8)(p)
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
            algebras of ``self`` and ``other`` canonically induced by
            the filtration-preserving map ``f : self -> other``.

            Let `A` and `B` be two filtered algebras with basis, and let
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

            We start with the universal enveloping algebra of the
            Lie algebra `\RR^3` (with the cross product serving as
            Lie bracket)::

                sage: A = AlgebrasWithBasis(QQ).Filtered().example(); A
                An example of a filtered algebra with basis: the
                 universal enveloping algebra of Lie algebra of RR^3
                 with cross product over Rational Field
                sage: M = A.indices(); M
                Free abelian monoid indexed by {'x', 'y', 'z'}
                sage: x,y,z = [A.basis()[M.gens()[i]] for i in "xyz"]

            Let us define a stupid filtered map from ``A`` to
            itself::

                sage: def map_on_basis(m):
                ....:     d = m.dict()
                ....:     i = d.get('x', 0); j = d.get('y', 0); k = d.get('z', 0)
                ....:     g = (y ** (i+j)) * (z ** k)
                ....:     if i > 0:
                ....:         g += i * (x ** (i-1)) * (y ** j) * (z ** k)
                ....:     return g
                sage: f = A.module_morphism(on_basis=map_on_basis,
                ....:                       codomain=A)
                sage: f(x)
                U['y'] + 1
                sage: f(x*y*z)
                U['y']^2*U['z'] + U['y']*U['z']
                sage: f(x*x*y*z)
                U['y']^3*U['z'] + 2*U['x']*U['y']*U['z']
                sage: f(A.one())
                1
                sage: f(y*z)
                U['y']*U['z']

            (There is nothing here that is peculiar to this
            universal enveloping algebra; we are only using its
            module structure, and we could just as well be using
            a polynomial algebra in its stead.)

            We now compute `\operatorname{gr} f` ::

                sage: grA = A.graded_algebra(); grA
                Graded Algebra of An example of a filtered algebra with
                 basis: the universal enveloping algebra of Lie algebra
                 of RR^3 with cross product over Rational Field
                sage: xx, yy, zz = [A.to_graded_conversion()(i) for i in [x, y, z]]
                sage: xx+yy*zz
                bar(U['y']*U['z']) + bar(U['x'])
                sage: grf = A.induced_graded_map(A, f); grf
                Generic endomorphism of Graded Algebra of An example
                 of a filtered algebra with basis: the universal
                 enveloping algebra of Lie algebra of RR^3 with cross
                 product over Rational Field
                sage: grf(xx)
                bar(U['y'])
                sage: grf(xx*yy*zz)
                bar(U['y']^2*U['z'])
                sage: grf(xx*xx*yy*zz)
                bar(U['y']^3*U['z'])
                sage: grf(grA.one())
                bar(1)
                sage: grf(yy*zz)
                bar(U['y']*U['z'])
                sage: grf(yy*zz-2*yy)
                bar(U['y']*U['z']) - 2*bar(U['y'])

            **Example 2.**

            We shall now construct `\operatorname{gr} f` for a
            different map `f` out of the same ``A``; the new map
            `f` will lead into a graded algebra already, namely into
            the algebra of symmetric functions::

                sage: # needs sage.combinat sage.modules
                sage: h = SymmetricFunctions(QQ).h()
                sage: def map_on_basis(m):  # redefining map_on_basis
                ....:     d = m.dict()
                ....:     i = d.get('x', 0); j = d.get('y', 0); k = d.get('z', 0)
                ....:     g = (h[1] ** i) * (h[2] ** (j // 2) * (h[3] ** (k // 3)))
                ....:     g += i * (h[1] ** (i+j+k))
                ....:     return g
                sage: f = A.module_morphism(on_basis=map_on_basis,
                ....:                       codomain=h)  # redefining f
                sage: f(x)
                2*h[1]
                sage: f(y)
                h[]
                sage: f(z)
                h[]
                sage: f(y**2)
                h[2]
                sage: f(x**2)
                3*h[1, 1]
                sage: f(x*y*z)
                h[1] + h[1, 1, 1]
                sage: f(x*x*y*y*z)
                2*h[1, 1, 1, 1, 1] + h[2, 1, 1]
                sage: f(A.one())
                h[]

            The algebra ``h`` of symmetric functions in the `h`-basis
            is already graded, so its associated graded algebra is
            implemented as itself::

                sage: # needs sage.combinat sage.modules
                sage: grh = h.graded_algebra(); grh is h
                True
                sage: grf = A.induced_graded_map(h, f); grf
                Generic morphism:
                  From: Graded Algebra of An example of a filtered
                   algebra with basis: the universal enveloping
                   algebra of Lie algebra of RR^3 with cross
                   product over Rational Field
                  To:   Symmetric Functions over Rational Field
                   in the homogeneous basis
                sage: grf(xx)
                2*h[1]
                sage: grf(yy)
                0
                sage: grf(zz)
                0
                sage: grf(yy**2)
                h[2]
                sage: grf(xx**2)
                3*h[1, 1]
                sage: grf(xx*yy*zz)
                h[1, 1, 1]
                sage: grf(xx*xx*yy*yy*zz)
                2*h[1, 1, 1, 1, 1]
                sage: grf(grA.one())
                h[]

            **Example 3.**

            After having had a graded algebra as the codomain, let us try to
            have one as the domain instead. Our new ``f`` will go from ``h``
            to ``A``::

                sage: # needs sage.combinat sage.modules
                sage: def map_on_basis(lam):  # redefining map_on_basis
                ....:     return x ** (sum(lam)) + y ** (len(lam))
                sage: f = h.module_morphism(on_basis=map_on_basis,
                ....:                       codomain=A)  # redefining f
                sage: f(h[1])
                U['x'] + U['y']
                sage: f(h[2])
                U['x']^2 + U['y']
                sage: f(h[1, 1])
                U['x']^2 + U['y']^2
                sage: f(h[2, 2])
                U['x']^4 + U['y']^2
                sage: f(h[3, 2, 1])
                U['x']^6 + U['y']^3
                sage: f(h.one())
                2
                sage: grf = h.induced_graded_map(A, f); grf
                Generic morphism:
                  From: Symmetric Functions over Rational Field
                   in the homogeneous basis
                  To:   Graded Algebra of An example of a filtered
                   algebra with basis: the universal enveloping
                   algebra of Lie algebra of RR^3 with cross
                   product over Rational Field
                sage: grf(h[1])
                bar(U['x']) + bar(U['y'])
                sage: grf(h[2])
                bar(U['x']^2)
                sage: grf(h[1, 1])
                bar(U['x']^2) + bar(U['y']^2)
                sage: grf(h[2, 2])
                bar(U['x']^4)
                sage: grf(h[3, 2, 1])
                bar(U['x']^6)
                sage: grf(h.one())
                2*bar(1)

            **Example 4.**

            The construct `\operatorname{gr} f` also makes sense when `f`
            is a filtration-preserving map between graded algebras. ::

                sage: # needs sage.combinat sage.modules
                sage: def map_on_basis(lam):  # redefining map_on_basis
                ....:     return h[lam] + h[len(lam)]
                sage: f = h.module_morphism(on_basis=map_on_basis,
                ....:                       codomain=h)  # redefining f
                sage: f(h[1])
                2*h[1]
                sage: f(h[2])
                h[1] + h[2]
                sage: f(h[1, 1])
                h[1, 1] + h[2]
                sage: f(h[2, 1])
                h[2] + h[2, 1]
                sage: f(h.one())
                2*h[]
                sage: grf = h.induced_graded_map(h, f); grf
                Generic endomorphism of Symmetric Functions over Rational
                 Field in the homogeneous basis
                sage: grf(h[1])
                2*h[1]
                sage: grf(h[2])
                h[2]
                sage: grf(h[1, 1])
                h[1, 1] + h[2]
                sage: grf(h[2, 1])
                h[2, 1]
                sage: grf(h.one())
                2*h[]

            **Example 5.**

            For another example, let us compute `\operatorname{gr} f` for a
            map `f` between two Clifford algebras::

                sage: # needs sage.modules
                sage: Q = QuadraticForm(ZZ, 2, [1,2,3])
                sage: B = CliffordAlgebra(Q, names=['u','v']); B
                The Clifford algebra of the Quadratic form in 2
                 variables over Integer Ring with coefficients:
                [ 1 2 ]
                [ * 3 ]
                sage: m = Matrix(ZZ, [[1, 2], [1, -1]])
                sage: f = B.lift_module_morphism(m, names=['x','y'])
                sage: A = f.domain(); A
                The Clifford algebra of the Quadratic form in 2
                 variables over Integer Ring with coefficients:
                [ 6 0 ]
                [ * 3 ]
                sage: x, y = A.gens()
                sage: f(x)
                u + v
                sage: f(y)
                2*u - v
                sage: f(x**2)
                6
                sage: f(x*y)
                -3*u*v + 3
                sage: grA = A.graded_algebra(); grA
                The exterior algebra of rank 2 over Integer Ring
                sage: A.to_graded_conversion()(x)
                x
                sage: A.to_graded_conversion()(y)
                y
                sage: A.to_graded_conversion()(x*y)
                x*y
                sage: u = A.to_graded_conversion()(x*y+1); u
                x*y + 1
                sage: A.from_graded_conversion()(u)
                x*y + 1
                sage: A.projection(2)(x*y+1)
                x*y
                sage: A.projection(1)(x+2*y-2)
                x + 2*y
                sage: grf = A.induced_graded_map(B, f); grf
                Generic morphism:
                  From: The exterior algebra of rank 2 over Integer Ring
                  To:   The exterior algebra of rank 2 over Integer Ring
                sage: grf(A.to_graded_conversion()(x))
                u + v
                sage: grf(A.to_graded_conversion()(y))
                2*u - v
                sage: grf(A.to_graded_conversion()(x**2))
                6
                sage: grf(A.to_graded_conversion()(x*y))
                -3*u*v
                sage: grf(grA.one())
                1
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
        pass
