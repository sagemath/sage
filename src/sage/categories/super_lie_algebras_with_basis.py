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
        Category of super lie algebras with basis over Integer Ring
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

    class FiniteDimensional(CategoryWithAxiom_over_base_ring):
        """
        The category of finite dimensional super Lie algebras with basis.

        EXAMPLES::

            sage: C = LieAlgebras(ZZ).WithBasis().Super().FiniteDimensional()
            sage: C is LieAlgebras(ZZ).WithBasis().FiniteDimensional().Super()
            True
        """
