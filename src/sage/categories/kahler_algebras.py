r"""
Category of Kahler Algebras.

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

        def hodge_riemann_relations(self, k, lefschetz_el, r):
            basis_k = self.basis(d=k)
            coeff = []
            for el, i in enumerate(basis_k):
                for j in range(i+1, len(basis_k)):
                    coeff.append((el*(lefschetz_el**(r-(2*k))*basis_k[j])).degree())
            return QuadraticForm(self.base_ring(), len(basis_k), coeff)