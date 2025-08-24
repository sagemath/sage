"""
Vectors over cartesian product of rings
"""
from sage.modules.free_module_element cimport FreeModuleElement_generic_dense
from sage.structure.element cimport Vector


cdef class Vector_cartesian_product(FreeModuleElement_generic_dense):

    def _to_cartesian_factors(self):
        return [self.apply_map(lambda x: x[i]) for i in range(len(self.base_ring().cartesian_factors()))]

    def __truediv__(self, other):
        """
        Division of the vector ``self`` by the scalar, vector or matrix ``other``.

        TESTS:

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
        if isinstance(other, Vector_cartesian_product):
            m = len(self)
            assert m == len(other), "sizes of vectors are different"

            base_ring = self.base_ring()

            # Cartesian product ring may not admit a basis, therefore division is
            # instead being performed component-wise and then stitched back
            self_facs = self._to_cartesian_factors()
            other_facs = other._to_cartesian_factors()

            result = [(self_facs[j]) / (other_facs[j]) for j in range(len(self_facs))]

            # Convert result to cartesian product
            return base_ring._cartesian_product_of_elements(result)

        # fallback
        return Vector.__truediv__(self, other)
