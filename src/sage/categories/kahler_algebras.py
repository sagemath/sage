r"""
Kähler Algebras

AUTHORS:

- Shriya M
"""

from sage.categories.category_types import Category_over_base_ring
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.misc.abstract_method import abstract_method
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.misc.cachefunc import cached_method


# ****************************************************************************
#       Copyright (C) 2024 Shriya M <25shriya at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
class KahlerAlgebras(Category_over_base_ring):
    r"""
    The category of graded algebras satisfying the Kähler package.

    EXAMPLES::

        sage: from sage.categories.kahler_algebras import KahlerAlgebras

        sage: C = KahlerAlgebras(QQ); C
        Category of kahler algebras over Rational Field
        sage: sorted(C.super_categories(), key=str)
        [Category of graded algebras with basis over Rational Field]

    TESTS::

        sage: C = KahlerAlgebras(QQ)
        sage: TestSuite(C).run()
    """
    def super_categories(self):
        return [GradedAlgebrasWithBasis(self.base_ring())]

    class ParentMethods:
        @abstract_method
        def poincare_pairing():
            pass

        @cached_method
        def _top_degree(self):
            r"""
            Return the top degree of the Kähler algebra.

            EXAMPLES::

                sage: ch = matroids.Uniform(4,6).chow_ring(QQ, False)
                sage: ch._top_degree()
                3
                sage: ch = matroids.Wheel(3).chow_ring(QQ, 'atom-free')
                sage: ch._top_degree()
                2
            """
            return max([b.degree() for b in self.basis()])

  # check all methods with matroids with loops/parallel elements
  # add properties and Kahler algebra def. Be as detailed as possible
  # issue with category

        @abstract_method
        def lefschetz_element():
            pass

        def hodge_riemann_relations(self, k, r):
            r"""
            Return the quadratic form for the corresponding k (< r/2) for the
            Kähler algebra.

            EXAMPLES::

                sage: ch = matroids.Uniform(4,6).chow_ring(QQ, False)
                sage: ch.hodge_riemann_relations(1, ch.matroid().rank() - 1)
                Quadratic form in 36 variables over Rational Field with coefficients:
                [ 3 -1 -1 3 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 3 ]
                [ * 3 -1 3 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 3 ]
                [ * * 3 3 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 3 ]
                [ * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * 3 -1 3 -1 3 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 3 ]
                [ * * * * * 3 3 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 3 ]
                [ * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * 3 3 3 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 ]
                [ * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * 3 -1 3 -1 -1 3 -1 -1 3 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 ]
                [ * * * * * * * * * * * 3 3 -1 3 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 ]
                [ * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * 3 3 3 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 3 ]
                [ * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * 3 3 3 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 3 ]
                [ * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * 3 -1 3 -1 3 -1 -1 -1 -1 3 -1 -1 -1 3 -1 3 ]
                [ * * * * * * * * * * * * * * * * * * * * * 3 3 -1 -1 3 -1 -1 3 -1 -1 -1 3 -1 -1 3 ]
                [ * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * 3 3 3 -1 3 -1 -1 -1 3 -1 -1 -1 3 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * 3 3 3 3 -1 -1 -1 -1 3 3 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 3 3 3 3 3 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 ]
                sage: ch.hodge_riemann_relations(3, ch.matroid().rank() - 1)
                Traceback (most recent call last):
                ...
                ValueError: k must be less than r < 2
            """
            if k > (r/2):
                raise ValueError("k must be less than r < 2")
            basis_k = []
            lefschetz_el = self.lefschetz_element()
            for b in self.basis():
                if b.homogeneous_degree() == k:
                    basis_k.append(b)
            coeff = []
            for i,el in enumerate(basis_k):
                for j in range(i, len(basis_k)):
                    coeff.append((el * (lefschetz_el ** (r-(2*k)) * basis_k[j])).degree())
            return QuadraticForm(self.base_ring(), len(basis_k), coeff)