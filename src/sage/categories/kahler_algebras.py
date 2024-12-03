r"""
Category of Kahler Algebras.

AUTHORS:

- Shriya M
"""

from sage.categories.category_types import Category_over_base_ring
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis

class KahlerAlgebras(Category_over_base_ring):
    class ParentMethods:
        def super_categories(self):
            return GradedAlgebrasWithBasis
        def poincare_pairing(self, a, b, r):
            if (a.homogeneous_degree() <= (r/2)) & (b.homogeneous_degree() == (r - a.homogeneous_degree())):
                el = a*b
                return el.degree()
                



            
    
    