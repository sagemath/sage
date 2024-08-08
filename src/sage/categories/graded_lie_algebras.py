# sage_setup: distribution = sagemath-categories
r"""
Graded Lie Algebras

AUTHORS:

- Eero Hakavuori (2018-08-16): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Eero Hakavuori <eero.hakavuori@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.graded_modules import GradedModulesCategory


class GradedLieAlgebras(GradedModulesCategory):
    r"""
    Category of graded Lie algebras.

    TESTS::

        sage: C = LieAlgebras(QQ).Graded()
        sage: TestSuite(C).run()
    """
    class SubcategoryMethods:
        def Stratified(self):
            r"""
            Return the full subcategory of stratified objects of ``self``.

            A Lie algebra is stratified if it is graded and generated as a
            Lie algebra by its component of degree one.

            EXAMPLES::

                sage: LieAlgebras(QQ).Graded().Stratified()
                Category of stratified Lie algebras over Rational Field
            """
            return self._with_axiom("Stratified")

    class Stratified(CategoryWithAxiom_over_base_ring):
        r"""
        Category of stratified Lie algebras.

        A graded Lie algebra `L = \bigoplus_{k=1}^M L_k` (where
        possibly `M = \infty`) is called *stratified* if it is generated
        by `L_1`; in other words, we have `L_{k+1} = [L_1, L_k]`.

        TESTS::

            sage: C = LieAlgebras(QQ).Graded().Stratified()
            sage: TestSuite(C).run()
        """
        class FiniteDimensional(CategoryWithAxiom_over_base_ring):
            r"""
            Category of finite dimensional stratified Lie algebras.

            EXAMPLES::

                sage: LieAlgebras(QQ).Graded().Stratified().FiniteDimensional()
                Category of finite dimensional stratified Lie algebras over Rational Field

            TESTS::

                sage: C = LieAlgebras(QQ).Graded().Stratified().FiniteDimensional()
                sage: TestSuite(C).run()
            """
            def extra_super_categories(self):
                """
                Implement the fact that a finite dimensional stratified Lie
                algebra is nilpotent.

                EXAMPLES::

                    sage: C = LieAlgebras(QQ).Graded().Stratified().FiniteDimensional()
                    sage: C.extra_super_categories()
                    [Category of nilpotent Lie algebras over Rational Field]
                    sage: C is C.Nilpotent()
                    True
                    sage: C.is_subcategory(LieAlgebras(QQ).Nilpotent())
                    True
                """
                from sage.categories.lie_algebras import LieAlgebras
                return [LieAlgebras(self.base_ring()).Nilpotent()]
