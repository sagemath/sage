# sage_setup: distribution = sagemath-categories
r"""
Kac-Moody Algebras

AUTHORS:

- Travis Scrimshaw (07-15-2017): Initial implementation
"""

# ****************************************************************************
#       Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.cachefunc import cached_method
from sage.categories.category_types import Category_over_base_ring
from sage.categories.lie_algebras import LieAlgebras


class KacMoodyAlgebras(Category_over_base_ring):
    """
    Category of Kac-Moody algebras.
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.kac_moody_algebras import KacMoodyAlgebras
            sage: KacMoodyAlgebras(QQ).super_categories()
            [Category of Lie algebras over Rational Field]
        """
        return [LieAlgebras(self.base_ring())]

    def example(self, n=2):
        """
        Return an example of a Kac-Moody algebra as per
        :meth:`Category.example <sage.categories.category.Category.example>`.

        EXAMPLES::

            sage: from sage.categories.kac_moody_algebras import KacMoodyAlgebras
            sage: KacMoodyAlgebras(QQ).example()                                        # needs sage.combinat sage.modules
            Lie algebra of ['A', 2] in the Chevalley basis

        We can specify the rank of the example::

            sage: KacMoodyAlgebras(QQ).example(4)                                       # needs sage.combinat sage.modules
            Lie algebra of ['A', 4] in the Chevalley basis
        """
        from sage.algebras.lie_algebras.classical_lie_algebra import LieAlgebraChevalleyBasis
        return LieAlgebraChevalleyBasis(self.base_ring(), ['A', n])

    class ParentMethods:
        def cartan_type(self):
            """
            Return the Cartan type of ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['A', 2])                          # needs sage.combinat sage.modules
                sage: L.cartan_type()                                                   # needs sage.combinat sage.modules
                ['A', 2]
            """
            return self._cartan_type

        def weyl_group(self):
            """
            Return the Weyl group of ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['A', 2])                          # needs sage.combinat sage.modules
                sage: L.weyl_group()                                                    # needs sage.combinat sage.modules
                Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)
            """
            from sage.combinat.root_system.weyl_group import WeylGroup
            return WeylGroup(self.cartan_type())
