"""
Lie Conformal Algebras With Basis

AUTHORS:

- Reimundo Heluani (2019-10-05): Initial implementation.
"""

#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.graded_lie_conformal_algebras import GradedLieConformalAlgebrasCategory
from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.super_modules import SuperModulesCategory

class LieConformalAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    The category of Lie conformal algebras with basis.

    EXAMPLES::

        sage: LieConformalAlgebras(QQbar).WithBasis()
        Category of Lie conformal algebras with basis over Algebraic Field
    """
    class Super(SuperModulesCategory):
        """
        The category of super Lie conformal algebras with basis.

        EXAMPLES::

            sage: LieConformalAlgebras(AA).WithBasis().Super()
            Category of super Lie conformal algebras with basis over Algebraic Real Field
        """
        class ParentMethods:

            def _even_odd_on_basis(self, m):
                """
                Return the parity of the basis element indexed by ``m``.

                OUTPUT:

                ``0`` if ``m`` is for an even element or ``1`` if ``m``
                is for an odd element.

                EXAMPLES::

                    sage: V = lie_conformal_algebras.NeveuSchwarz(QQ)
                    sage: B = V._indices
                    sage: V._even_odd_on_basis(B(('G',1)))
                    1
                """
                return self._parity[self.monomial((m[0],0))]

        class Graded(GradedLieConformalAlgebrasCategory):
            """
            The category of H-graded super Lie conformal algebras with basis.

            EXAMPLES::

                sage: LieConformalAlgebras(QQbar).WithBasis().Super().Graded()
                Category of H-graded super Lie conformal algebras with basis over Algebraic Field
            """

    class Graded(GradedLieConformalAlgebrasCategory):
        """
        The category of H-graded Lie conformal algebras with basis.

        EXAMPLES::

            sage: LieConformalAlgebras(QQbar).WithBasis().Graded()
            Category of H-graded Lie conformal algebras with basis over Algebraic Field
        """

    class FinitelyGeneratedAsLambdaBracketAlgebra(CategoryWithAxiom_over_base_ring):
        """
        The category of finitely generated Lie conformal
        algebras with basis.

        EXAMPLES::

            sage: C = LieConformalAlgebras(QQbar)
            sage: C.WithBasis().FinitelyGenerated()
            Category of finitely generated Lie conformal algebras with basis over Algebraic Field
            sage: C.WithBasis().FinitelyGenerated() is C.FinitelyGenerated().WithBasis()
            True
        """
        class Super(SuperModulesCategory):
            """
            The category of super finitely generated Lie conformal
            algebras with basis.

            EXAMPLES::

                sage: LieConformalAlgebras(AA).WithBasis().FinitelyGenerated().Super()
                Category of super finitely generated Lie conformal algebras with basis over Algebraic Real Field
            """
            class Graded(GradedModulesCategory):
                """
                The category of H-graded super finitely generated Lie
                conformal algebras with basis.

                EXAMPLES::

                    sage: C = LieConformalAlgebras(QQbar).WithBasis().FinitelyGenerated()
                    sage: C.Graded().Super()
                    Category of H-graded super finitely generated Lie conformal algebras with basis over Algebraic Field
                    sage: C.Graded().Super() is C.Super().Graded()
                    True
                """
                def _repr_object_names(self):
                    """
                    The names of the objects of ``self``.

                    EXAMPLES::

                        sage: C = LieConformalAlgebras(QQbar).WithBasis().FinitelyGenerated()
                        sage: C.Super().Graded()
                        Category of H-graded super finitely generated Lie conformal algebras with basis over Algebraic Field
                    """
                    return "H-graded {}".format(self.base_category()._repr_object_names())

        class Graded(GradedLieConformalAlgebrasCategory):
            """
            The category of H-graded finitely generated Lie conformal
            algebras with basis.

            EXAMPLES::

                sage: LieConformalAlgebras(QQbar).WithBasis().FinitelyGenerated().Graded()
                Category of H-graded finitely generated Lie conformal algebras with basis over Algebraic Field
            """
