# sage_setup: distribution = sagemath-categories
"""
Finitely Generated Lambda bracket Algebras

AUTHORS:

- Reimundo Heluani (2020-08-21): Initial implementation.
"""

#******************************************************************************
#       Copyright (C) 2020 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.lambda_bracket_algebras import LambdaBracketAlgebras


class FinitelyGeneratedLambdaBracketAlgebras(CategoryWithAxiom_over_base_ring):
    """
    The category of finitely generated lambda bracket algebras.

    EXAMPLES::

        sage: from sage.categories.lambda_bracket_algebras import LambdaBracketAlgebras
        sage: LambdaBracketAlgebras(QQbar).FinitelyGenerated()                          # needs sage.rings.number_field
        Category of finitely generated lambda bracket algebras over Algebraic Field
    """
    _base_category_class_and_axiom = (LambdaBracketAlgebras, "FinitelyGeneratedAsLambdaBracketAlgebra")

    class ParentMethods:
        def ngens(self):
            r"""
            The number of generators of this Lie conformal algebra.

            EXAMPLES::

                sage: Vir = lie_conformal_algebras.Virasoro(QQ)                         # needs sage.combinat sage.modules
                sage: Vir.ngens()                                                       # needs sage.combinat sage.modules
                2

                sage: V = lie_conformal_algebras.Affine(QQ, 'A2')                       # needs sage.combinat sage.modules
                sage: V.ngens()                                                         # needs sage.combinat sage.modules
                9
            """
            return len(self.gens())

        def gen(self, i):
            r"""
            The ``i``-th generator of this Lie conformal algebra.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: V = lie_conformal_algebras.Affine(QQ, 'A1')
                sage: V.gens()
                (B[alpha[1]], B[alphacheck[1]], B[-alpha[1]], B['K'])
                sage: V.gen(0)
                B[alpha[1]]
                sage: V.1
                B[alphacheck[1]]
            """
            return self.gens()[i]

        def some_elements(self):
            """
            Some elements of this Lie conformal algebra.

            This method returns a list with elements containing at
            least the generators.

            EXAMPLES::

                sage: V = lie_conformal_algebras.Affine(QQ, 'A1',                       # needs sage.combinat sage.modules
                ....:                                   names=('e', 'h', 'f'))
                sage: V.some_elements()                                                 # needs sage.combinat sage.modules
                [e, h, f, K, ...]
                sage: all(v.parent() is V for v in V.some_elements())                   # needs sage.combinat sage.modules
                True
            """
            S = list(self.gens())
            from sage.misc.misc import some_tuples
            for x,y in some_tuples(S, 2, 0, max_samples=self.ngens()):
                S.append(x.T() + 2*y.T(2))
            return S

    class Graded(GradedModulesCategory):
        """
        The category of H-graded finitely generated Lie conformal algebras.

        EXAMPLES::

            sage: LieConformalAlgebras(QQbar).FinitelyGenerated().Graded()              # needs sage.rings.number_field
            Category of H-graded finitely generated Lie conformal algebras
             over Algebraic Field
        """
        def _repr_object_names(self):
            """
            The names of the objects of ``self``.

            EXAMPLES::

                sage: LieConformalAlgebras(QQbar).WithBasis().FinitelyGenerated().Graded()          # needs sage.rings.number_field
                Category of H-graded finitely generated Lie conformal algebras with basis over Algebraic Field
            """
            return "H-graded {}".format(self.base_category()._repr_object_names())
