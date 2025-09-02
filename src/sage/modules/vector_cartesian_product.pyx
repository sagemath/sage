"""
Vectors over cartesian product of rings
"""
from sage.modules.free_module_element cimport FreeModuleElement_generic_dense
from sage.structure.coerce cimport py_scalar_to_element
from sage.structure.element cimport Vector

from sage.categories.cartesian_product import cartesian_product
from sage.sets.cartesian_product import CartesianProduct

cdef class Vector_cartesian_product(FreeModuleElement_generic_dense):
    def _to_cartesian_factors(self):
        return (self.apply_map(lambda x: x[i]) for i in range(len(self.base_ring().cartesian_factors())))

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

            sage: A = matrix([R(1)])
            sage: b / A # vector-by-matrix
            ((2, 2))

        Test for conversion between ZZ and QQ as component rings::

            sage: R = cartesian_product([ZZ, ZZ])
            sage: v = vector([R.one()])
            sage: two = R(2)
            sage: w = (v / two); w
            ((1/2, 1/2))
            sage: v / w
            (2, 2)
        """
        other = py_scalar_to_element(other)

        # vector-by-vector division
        if isinstance(other, Vector_cartesian_product):
            m = len(self)
            assert m == len(other), "sizes of vectors are different"

            # Cartesian product ring may not admit a basis, therefore division is
            # instead being performed component-wise and then stitched back
            self_facs = self._to_cartesian_factors()
            other_facs = other._to_cartesian_factors()

            result = [v1 / v2 for v1, v2 in zip(self_facs, other_facs)]

            # Convert result to cartesian product
            return cartesian_product(result)

        # vector-by-scalar division
        if isinstance(other.parent(), CartesianProduct):
            # __invert__ method does not allow inverting eg 2 in ZZ
            inverted = cartesian_product([~x for x in other.cartesian_factors()])

            return self * inverted

        # fallback
        return Vector.__truediv__(self, other)
