r"""
Category of Kähler Algebras.

AUTHORS:

- Shriya M
"""

from sage.categories.category_types import Category_over_base_ring
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.misc.abstract_method import abstract_method
from sage.quadratic_forms.quadratic_form import QuadraticForm


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
        def poincare_pairing(self, el1, el2, r):
            r"""
            Return the Poincaré pairing of any two elements of the
            Kähler algebra.

            EXAMPLES::

                sage: ch = matroids.Wheel(3).chow_ring(QQ, True, 'atom-free')
                sage: A0, A1, A2, A3, A4, A5, A013, A025, A04, A124, A15, A23, A345, A012345 = ch.gens()
                sage: u = ch(-1/6*A2*A012345 + 41/48*A012345^2); u
                -1/6*A2*A012345 + 41/48*A012345^2
                sage: v = ch(-A345^2 - 1/4*A345); v
                -A345^2 - 1/4*A345
                sage: ch.poincare_pairing(v, u, ch.matroid().rank())
                3
            """
            hom_components1 = el1.lift().homogeneous_components()
            hom_components2 = el2.lift().homogeneous_components()
            new_el = self.base_ring().zero()
            for i in hom_components1:
                for j in hom_components2:
                    if i == r - j:
                        new_el += hom_components1[i] * hom_components2[j]
                    #  the 'else' case is new_el += self.base_ring().zero()
            return new_el.degree()

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
            if k < (r/2):
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
            else:
                raise ValueError("k must be less than r < 2")