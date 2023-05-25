"""
Finitely Generated Lie Conformal Algebras

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
from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.super_modules import SuperModulesCategory
from sage.categories.lie_conformal_algebras import LieConformalAlgebras

class FinitelyGeneratedLieConformalAlgebras(CategoryWithAxiom_over_base_ring):
    """
    The category of finitely generated Lie conformal algebras.

    EXAMPLES::

        sage: LieConformalAlgebras(QQbar).FinitelyGenerated()
        Category of finitely generated lie conformal algebras over Algebraic Field
    """
    _base_category_class_and_axiom = (LieConformalAlgebras, "FinitelyGeneratedAsLambdaBracketAlgebra")

    class ParentMethods:

        def some_elements(self):
            """
            Some elements of this Lie conformal algebra.

            Returns a list with elements containing at least the
            generators.

            EXAMPLES::

                sage: V = lie_conformal_algebras.Affine(QQ, 'A1', names=('e', 'h', 'f'))
                sage: V.some_elements()
                [e, h, f, K, ...]
                sage: all(v.parent() is V for v in V.some_elements())
                True
            """
            S = list(self.gens())
            from sage.misc.misc import some_tuples
            for x,y in some_tuples(S, 2, 0, max_samples=self.ngens()):
                S.append(x.T() + 2*y.T(2))
            return S

    class Super(SuperModulesCategory):
        """
        The category of super finitely generated Lie conformal algebras.

        EXAMPLES::

            sage: LieConformalAlgebras(AA).FinitelyGenerated().Super()
            Category of super finitely generated lie conformal algebras over Algebraic Real Field
        """
        class Graded(GradedModulesCategory):
            """
            The category of H-graded super finitely generated Lie conformal algebras.

            EXAMPLES::

                sage: LieConformalAlgebras(QQbar).FinitelyGenerated().Super().Graded()
                Category of H-graded super finitely generated lie conformal algebras over Algebraic Field
            """
            def _repr_object_names(self):
                """
                The names of the objects of ``self``.

                EXAMPLES::

                    sage: C = LieConformalAlgebras(QQbar).FinitelyGenerated()
                    sage: C.Super().Graded()
                    Category of H-graded super finitely generated lie conformal algebras over Algebraic Field
                """
                return "H-graded {}".format(self.base_category()._repr_object_names())

    class Graded(GradedModulesCategory):
        """
        The category of H-graded finitely generated Lie conformal algebras.

        EXAMPLES::

            sage: LieConformalAlgebras(QQbar).FinitelyGenerated().Graded()
            Category of H-graded finitely generated lie conformal algebras over Algebraic Field
        """
        def _repr_object_names(self):
            """
            The names of the objects of ``self``.

            EXAMPLES::

                sage: LieConformalAlgebras(QQbar).WithBasis().FinitelyGenerated().Graded()
                Category of H-graded finitely generated Lie conformal algebras with basis over Algebraic Field
            """
            return "H-graded {}".format(self.base_category()._repr_object_names())
