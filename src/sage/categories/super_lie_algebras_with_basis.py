r"""
Super Lie algebras with basis
"""
# ****************************************************************************
#  Copyright (C) 2023 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.modules import Modules
from sage.categories.super_modules import SuperModulesCategory
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring

class SuperLieAlgebrasWithBasis(SuperModulesCategory):
    """
    The category of super Lie algebras with a distinguished basis.

    EXAMPLES::

        sage: C = LieAlgebras(ZZ).WithBasis().Super(); C
        Category of super Lie algebras with basis over Integer Ring
        sage: TestSuite(C).run()
    """
    def extra_super_categories(self):
        """
        EXAMPLES::

            sage: C = LieAlgebras(ZZ).WithBasis().Super()
            sage: sorted(C.super_categories(), key=str) # indirect doctest
            [Category of graded modules with basis over Integer Ring,
             Category of super Lie algebras over Integer Ring,
             Category of super modules with basis over Integer Ring]
        """
        return [Modules(self.base()).WithBasis().Graded()]

    class ParentMethods:
        # TODO: Move this to the more general SuperLieAlgebras.ParentMethods
        #   once that category has been created.
        def bracket(self, lhs, rhs):
            r"""
            Return the super Lie bracket ``[lhs, rhs]`` after coercing ``lhs``
            and ``rhs`` into elements of ``self``.

            If ``lhs`` and ``rhs`` are super Lie algebras, then this constructs
            the product space, and if only one of them is a super Lie algebra,
            then it constructs the corresponding ideal.

            EXAMPLES::

                sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z':1}}, (1,1,2))
                sage: L.bracket(x, x + y)
                z
                sage: L.bracket(x, 0)
                0
                sage: L.bracket(0, x)
                0

            Constructing the product space::

                sage: Z = L.bracket(L, L); Z
                Traceback (most recent call last):
                ...
                NotImplementedError: sub Lie superalgebras not yet implemented

            Constructing ideals::

                sage: L.bracket(x, L)
                Traceback (most recent call last):
                ...
                NotImplementedError: sub Lie superalgebras not yet implemented
                sage: L.bracket(L, y)
                Traceback (most recent call last):
                ...
                NotImplementedError: sub Lie superalgebras not yet implemented
            """
            if lhs in SuperLieAlgebrasWithBasis:
                raise NotImplementedError("sub Lie superalgebras not yet implemented")
                #if rhs in LieAlgebras:
                #    return lhs.product_space(rhs)
                #return lhs.ideal(rhs)
            elif rhs in SuperLieAlgebrasWithBasis:
                raise NotImplementedError("sub Lie superalgebras not yet implemented")
                #return rhs.ideal(lhs)
            return self(lhs)._bracket_(self(rhs))

    class FiniteDimensional(CategoryWithAxiom_over_base_ring):
        """
        The category of finite dimensional super Lie algebras with basis.

        EXAMPLES::

            sage: C = LieAlgebras(ZZ).WithBasis().Super().FiniteDimensional()
            sage: C is LieAlgebras(ZZ).WithBasis().FiniteDimensional().Super()
            True
        """
        class ParentMethods:
            def module(self, sparse=False):
                r"""
                Return an (ungraded) module corresponding to ``self``.

                EXAMPLES::

                    sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z':1}}, (1,1,2))
                    sage: L.module(True)
                    Sparse vector space of dimension 3 over Rational Field
                    sage: L.module(False)
                    Vector space of dimension 3 over Rational Field
                """
                from sage.modules.free_module import FreeModule
                return FreeModule(self.base_ring(), self.dimension(), sparse=sparse)

            def graded_morphism(self, on_generators, codomain=None, base_map=None, check=True):
                r"""
                Return a graded morphim of super Lie algebras.

                EXAMPLES::

                    sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z': 1}}, (1,1,2))
                    sage: phi = L.graded_morphism({x:x, y:y, z:z})
                    sage: phi
                    super Lie algebra endomorphism of Super Lie algebra generated
                     by (x, y, z) in degrees (1, 1, 2) over Rational Field
                      Defn: x |--> x
                            y |--> y
                            z |--> z
                    sage: psi = L.graded_morphism({x:x, y:y}, L)
                    sage: psi
                    super Lie algebra endomorphism of Super Lie algebra generated
                     by (x, y, z) in degrees (1, 1, 2) over Rational Field
                      Defn: x |--> x
                            y |--> y
                            z |--> z
                """
                from sage.algebras.lie_algebras.morphism import LieAlgebraMorphism_from_generators
                return LieAlgebraMorphism_from_generators(on_generators, domain=self,
                                                          codomain=codomain, base_map=base_map, check=check)
