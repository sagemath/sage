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

                sage: E.<a,b> = ExteriorAlgebra(QQ)
                sage: E.degree_negation((1 + a) * (1 + b))
                a*b - a - b + 1
                sage: E.degree_negation(E.zero())
                0

                sage: P = GradedModulesWithBasis(ZZ).example(); P
                An example of a graded module with basis: the free module on partitions over Integer Ring
                sage: pbp = lambda x: P.basis()[Partition(list(x))]
                sage: p = pbp([3,1]) - 2 * pbp([2]) + 4 * pbp([1])
                sage: P.degree_negation(p)
                -4*P[1] - 2*P[2] + P[3, 1]
            """
            base_one = self.base_ring().one()
            base_minusone = - base_one
            diag = lambda x: (base_one if self.degree_on_basis(x) % 2 == 0
                              else base_minusone)
            return self.sum_of_terms([(key, diag(key) * value)
                                      for key, value in
                                      element.monomial_coefficients(copy=False).items()])

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

                sage: E.<a,b> = ExteriorAlgebra(QQ)
                sage: ((1 + a) * (1 + b)).degree_negation()
                a*b - a - b + 1
                sage: E.zero().degree_negation()
                0

                sage: P = GradedModulesWithBasis(ZZ).example(); P
                An example of a graded module with basis: the free module on partitions over Integer Ring
                sage: pbp = lambda x: P.basis()[Partition(list(x))]
                sage: p = pbp([3,1]) - 2 * pbp([2]) + 4 * pbp([1])
                sage: p.degree_negation()
                -4*P[1] - 2*P[2] + P[3, 1]
            """
            return self.parent().degree_negation(self)

