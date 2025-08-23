"""
Vectors over cartesian product of rings
"""
from operator import truediv

from sage.modules.free_module_element import FreeModuleElement_generic_dense
from sage.structure.element import Vector


class Vector_cartesian_product_dense(FreeModuleElement_generic_dense):
    def __truediv__(self, other):
        """
        Division of the vector ``self`` by the scalar, vector or matrix ``other``.

        TESTS::

        Test if :issue:`40626` is fixed::

            sage: R = cartesian_product([QQ, QQ])
            sage: b = vector([R(2)])
            sage: b / R(2) # vector-by-scalar
            ((1, 1))
            sage: c = vector([R(1)])
            sage: b / c # vector-by-vector
            (2, 2)
            sage: d = vector([R(0)])
            sage: b / d # vector-by-vector
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero vector
            sage: A = matrix([R(1) / R(2)])
            sage: b / A # vector-by-matrix
            ((4, 4))

            sage: R2 = cartesian_product([ZZ, R])
            sage: v = vector([R2.one()])
            sage: two = R2(2)
            sage: v / (v / two)
            (2, (2, 2))
        """
        if isinstance(other, Vector_cartesian_product_dense):
            base_ring = self.base_ring()

            # Cartesian product ring may not admit a basis, therefore division is
            # instead being performed component-wise and then stitched back
            right_base_ring = other.base_ring()

            n = len(base_ring._sets)
            left_lsts, right_lsts = [[] for _ in range(n)], [[] for _ in range(n)]
            m = len(self)
            assert m == len(other), "sizes of vectors are different"

            # Split vectors into component vectors
            for i in range(m):
                for j in range(n):
                    left_lsts[j].append(self[i][j])
                    right_lsts[j].append(other[i][j])

            result = []
            for j in range(n):
                result.append((base_ring._sets[j]**m)(left_lsts[j]) / (right_base_ring._sets[j]**m)(right_lsts[j]))

            # Convert result to cartesian product
            return base_ring._cartesian_product_of_elements(result)

        # fallback
        return Vector.__truediv__(self, other)
