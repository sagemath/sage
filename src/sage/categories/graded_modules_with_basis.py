# sage_setup: distribution = sagemath-categories
r"""
Graded modules with basis
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.quotients import QuotientsCategory


class GradedModulesWithBasis(GradedModulesCategory):
    """
    The category of graded modules with a distinguished basis.

    EXAMPLES::

        sage: C = GradedModulesWithBasis(ZZ); C
        Category of graded modules with basis over Integer Ring
        sage: sorted(C.super_categories(), key=str)
        [Category of filtered modules with basis over Integer Ring,
         Category of graded modules over Integer Ring]
        sage: C is ModulesWithBasis(ZZ).Graded()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    class ParentMethods:
        def degree_negation(self, element):
            r"""
            Return the image of ``element`` under the degree negation
            automorphism of the graded module ``self``.

            The degree negation is the module automorphism which scales
            every homogeneous element of degree `k` by `(-1)^k` (for all
            `k`). This assumes that the module ``self`` is `\ZZ`-graded.

            INPUT:

            - ``element`` -- element of the module ``self``

            EXAMPLES::

                sage: E.<a,b> = ExteriorAlgebra(QQ)                                     # needs sage.modules
                sage: E.degree_negation((1 + a) * (1 + b))                              # needs sage.modules
                a*b - a - b + 1
                sage: E.degree_negation(E.zero())                                       # needs sage.modules
                0

                sage: P = GradedModulesWithBasis(ZZ).example(); P                       # needs sage.combinat sage.modules
                An example of a graded module with basis:
                 the free module on partitions over Integer Ring
                sage: pbp = lambda x: P.basis()[Partition(list(x))]
                sage: p = pbp([3,1]) - 2 * pbp([2]) + 4 * pbp([1])                      # needs sage.combinat sage.modules
                sage: P.degree_negation(p)                                              # needs sage.combinat sage.modules
                -4*P[1] - 2*P[2] + P[3, 1]
            """
            base_one = self.base_ring().one()
            base_minusone = - base_one
            diag = lambda x: (base_one if self.degree_on_basis(x) % 2 == 0
                              else base_minusone)
            return self.sum_of_terms([(key, diag(key) * value)
                                      for key, value in
                                      element.monomial_coefficients(copy=False).items()])

        def submodule(self, gens, check=True, already_echelonized=False,
                      unitriangular=False, support_order=None, category=None,
                      *args, **opts):
            r"""
            Return the submodule spanned by a finite set of elements.

            INPUT:

            - ``gens`` -- list or family of elements of ``self``
            - ``check`` -- boolean (default: ``True``); whether to verify that
              the elements of ``gens`` are in ``self``
            - ``already_echelonized`` -- boolean (default: ``False``); whether
               the elements of ``gens`` are already in (not necessarily
               reduced) echelon form
            - ``unitriangular`` -- boolean (default: ``False``); whether
              the lift morphism is unitriangular
            - ``support_order`` -- (optional) either something that can
              be converted into a tuple or a key function
            - ``category`` -- (optional) the category of the submodule

            If ``already_echelonized`` is ``False``, then the
            generators are put in reduced echelon form using
            :meth:`echelonize`, and reindexed by `0,1,...`.

            .. WARNING::

                At this point, this method only works for finite
                dimensional submodules and if matrices can be
                echelonized over the base ring.

            If in addition ``unitriangular`` is ``True``, then
            the generators are made such that the coefficients of
            the pivots are 1, so that lifting map is unitriangular.

            The basis of the submodule uses the same index set as the
            generators, and the lifting map sends `y_i` to `gens[i]`.

            .. SEEALSO::

                 - :meth:`ModulesWithBasis.FiniteDimensional.ParentMethods.quotient_module`
                 - :class:`sage.modules.with_basis.subquotient.SubmoduleWithBasis`

            EXAMPLES:

            A graded submodule of a graded module generated by homogeneous
            elements is naturally graded::

                sage: # needs sage.combinat sage.modules
                sage: E.<x,y,z> = ExteriorAlgebra(QQ)
                sage: S = E.submodule([x + y, x*y - y*z])
                sage: S.category()
                Join of
                 Category of graded vector spaces with basis over Rational Field and
                 Category of subobjects of filtered modules with basis over Rational Field and
                 Category of finite dimensional filtered modules with basis over Rational Field
                sage: S.basis()[0].degree()
                1
                sage: S.basis()[1].degree()
                2

            We check on the echelonized basis::

                sage: Sp = E.submodule([1, x + y + 5, x*y - y*z + x + y - 2])           # needs sage.combinat sage.modules
                sage: Sp.category()                                                     # needs sage.combinat sage.modules
                Join of
                 Category of graded vector spaces with basis over Rational Field and
                 Category of subobjects of filtered modules with basis over Rational Field and
                 Category of finite dimensional filtered modules with basis over Rational Field

            If it is generated by inhomogeneous elements, then it is
            filtered by default::

                sage: F = E.submodule([x + y*z, x*z + y*x])                             # needs sage.combinat sage.modules
                sage: F.category()                                                      # needs sage.combinat sage.modules
                Join of
                 Category of subobjects of filtered modules with basis over Rational Field and
                 Category of finite dimensional filtered modules with basis over Rational Field and
                 Category of filtered vector spaces with basis over Rational Field

            If ``category`` is specified, then it does not give any extra
            structure to the submodule (we can think of this as applying
            the forgetful functor)::

                sage: # needs sage.combinat sage.modules
                sage: SM = E.submodule([x + y, x*y - y*z],
                ....:                  category=ModulesWithBasis(QQ))
                sage: SM.category()
                Join of
                 Category of finite dimensional vector spaces with basis over Rational Field and
                 Category of subobjects of sets
                sage: FM = E.submodule([x + 1, x*y - x*y*z],
                ....:                  category=ModulesWithBasis(QQ))
                sage: FM.category()
                Join of
                 Category of finite dimensional vector spaces with basis over Rational Field and
                 Category of subobjects of sets

            If we have specified that this is a graded submodule of a graded
            module, then the echelonized elements must be homogeneous::

                sage: Cat = ModulesWithBasis(QQ).Graded().Subobjects()
                sage: E.submodule([x + y, x*y - 1], category=Cat)                       # needs sage.combinat sage.modules
                Traceback (most recent call last):
                ...
                ValueError: all of the generators must be homogeneous
                sage: E.submodule([x + y, x*y - x - y], category=Cat)                   # needs sage.combinat sage.modules
                Free module generated by {0, 1} over Rational Field
            """
            # Make sure gens consists of elements of ``self``
            from sage.sets.family import Family, AbstractFamily
            if isinstance(gens, AbstractFamily):
                gens = gens.map(self)
            elif isinstance(gens, dict):
                gens = Family(gens.keys(), gens.__getitem__)
            else:
                gens = [self(y) for y in gens]
            support_order = self._compute_support_order(gens, support_order)
            if not already_echelonized:
                gens = self.echelon_form(gens, unitriangular, order=support_order)

            GMod = GradedModulesWithBasis(self.category().base_ring())
            if category is None:
                if all(g.is_homogeneous() for g in gens):
                    category = GMod.Subobjects()
            elif (category.is_subcategory(GMod.Subobjects())
                  and not all(g.is_homogeneous() for g in gens)):
                raise ValueError("all of the generators must be homogeneous")

            from sage.modules.with_basis.subquotient import SubmoduleWithBasis
            return SubmoduleWithBasis(gens, ambient=self,
                                      support_order=support_order,
                                      unitriangular=unitriangular,
                                      category=category, *args, **opts)

        def quotient_module(self, submodule, check=True, already_echelonized=False, category=None):
            r"""
            Construct the quotient module ``self`` / ``submodule``.

            INPUT:

            - ``submodule`` -- a submodule with basis of ``self``, or
              something that can be turned into one via
              ``self.submodule(submodule)``
            - ``check``, ``already_echelonized`` -- passed down to
              :meth:`ModulesWithBasis.ParentMethods.submodule`
            - ``category`` -- (optional) the category of the quotient module

            .. WARNING::

                At this point, this only supports quotients by free
                submodules admitting a basis in unitriangular echelon
                form. In this case, the quotient is also a free
                module, with a basis consisting of the retract of a
                subset of the basis of ``self``.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: E.<x,y,z> = ExteriorAlgebra(QQ)
                sage: S = E.submodule([x + y, x*y - y*z, y])
                sage: Q = E.quotient_module(S)
                sage: Q.category()
                Join of
                 Category of quotients of graded modules with basis over Rational Field and
                 Category of graded vector spaces with basis over Rational Field and
                 Category of finite dimensional filtered modules with basis over Rational Field

            .. SEEALSO::

                 - :meth:`Modules.WithBasis.ParentMethods.submodule`
                 - :meth:`Rings.ParentMethods.quotient`
                 - :class:`sage.modules.with_basis.subquotient.QuotientModuleWithBasis`
            """
            from sage.modules.with_basis.subquotient import SubmoduleWithBasis, QuotientModuleWithBasis
            if not isinstance(submodule, SubmoduleWithBasis):
                submodule = self.submodule(submodule, check=check,
                                           unitriangular=True,
                                           already_echelonized=already_echelonized)

            GMod = GradedModulesWithBasis(self.category().base_ring())
            if category is None:
                if all(g.is_homogeneous() for g in submodule.basis()):
                    category = GMod.Quotients()
            elif (category.is_subcategory(GMod.Quotients())
                  and not all(g.is_homogeneous() for g in submodule.basis())):
                raise ValueError("all of the basis elements must be homogeneous")

            return QuotientModuleWithBasis(submodule, category=category)

    class ElementMethods:
        def degree_negation(self):
            r"""
            Return the image of ``self`` under the degree negation
            automorphism of the graded module to which ``self`` belongs.

            The degree negation is the module automorphism which scales
            every homogeneous element of degree `k` by `(-1)^k` (for all
            `k`). This assumes that the module to which ``self`` belongs
            (that is, the module ``self.parent()``) is `\ZZ`-graded.

            EXAMPLES::

                sage: E.<a,b> = ExteriorAlgebra(QQ)                                     # needs sage.modules
                sage: ((1 + a) * (1 + b)).degree_negation()                             # needs sage.modules
                a*b - a - b + 1
                sage: E.zero().degree_negation()                                        # needs sage.modules
                0

                sage: P = GradedModulesWithBasis(ZZ).example(); P                       # needs sage.combinat sage.modules
                An example of a graded module with basis:
                 the free module on partitions over Integer Ring
                sage: pbp = lambda x: P.basis()[Partition(list(x))]
                sage: p = pbp([3,1]) - 2 * pbp([2]) + 4 * pbp([1])                      # needs sage.combinat sage.modules
                sage: p.degree_negation()                                               # needs sage.combinat sage.modules
                -4*P[1] - 2*P[2] + P[3, 1]
            """
            return self.parent().degree_negation(self)

    class Quotients(QuotientsCategory):
        class ParentMethods:
            def degree_on_basis(self, m):
                r"""
                Return the degree of the basis element indexed by ``m``
                in ``self``.

                EXAMPLES::

                    sage: # needs sage.combinat sage.modules
                    sage: E.<x,y,z> = ExteriorAlgebra(QQ)
                    sage: S = E.submodule([x + y, x*y - y*z, y])
                    sage: Q = E.quotient_module(S)
                    sage: B = Q.basis()
                    sage: [B[i].lift() for i in Q.indices()]
                    [1, z, x*z, y*z, x*y*z]
                    sage: [Q.degree_on_basis(i) for i in Q.indices()]
                    [0, 1, 2, 2, 3]
                """
                return self.basis()[m].lift().degree()

        class ElementMethods:
            def degree(self):
                r"""
                Return the degree of ``self``.

                EXAMPLES::

                    sage: # needs sage.combinat sage.modules
                    sage: E.<x,y,z> = ExteriorAlgebra(QQ)
                    sage: S = E.submodule([x + y, x*y - y*z, y])
                    sage: Q = E.quotient_module(S)
                    sage: B = Q.basis()
                    sage: [B[i].lift() for i in Q.indices()]
                    [1, z, x*z, y*z, x*y*z]
                    sage: [B[i].degree() for i in Q.indices()]
                    [0, 1, 2, 2, 3]
                """
                return self.lift().degree()
