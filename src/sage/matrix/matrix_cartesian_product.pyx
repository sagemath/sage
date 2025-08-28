from .matrix_generic_dense cimport Matrix_generic_dense
from sage.structure.element cimport Matrix
from sage.categories.cartesian_product import cartesian_product

cdef class Matrix_cartesian_product(Matrix_generic_dense):
    """
    A class for matrices over the cartesian product of multiple rings.
    """
    def _to_cartesian_product_of_matrices(self):
        """
        Convert the matrix to a cartesian product of matrices.
        """
        return cartesian_product([
            self.apply_map(lambda x: x[i])
            for i in range(len(self.base_ring().cartesian_factors()))])

    @staticmethod
    def _from_cartesian_product_of_matrices(a):
        """
        Create a ``Matrix_cartesian_product`` from a cartesian product of matrices.
        """
        from sage.matrix.constructor import matrix
        return matrix(a[0].nrows(), a[0].ncols(), list(map(cartesian_product, zip(*(c.list() for c in a)))))

    cdef Matrix _matrix_times_matrix_(self, Matrix right):
        """
        EXAMPLES::

            sage: R = cartesian_product([ZZ, ZZ])
            sage: matrix.identity(R, 2)                                                 # needs sage.modules
            [(1, 1) (0, 0)]
            [(0, 0) (1, 1)]
            sage: type(matrix.identity(R, 2))                                           # needs sage.modules
            <class 'sage.matrix.matrix_cartesian_product.Matrix_cartesian_product'>
            sage: matrix.random(R, 100) * matrix.random(R, 100)  # should finish quickly            # needs sage.modules
            100 x 100 dense matrix over The Cartesian product of 2 copies of Integer Ring (use...)
        """
        return self._from_cartesian_product_of_matrices(self._to_cartesian_product_of_matrices() *
                                                        right._to_cartesian_product_of_matrices())

    def __invert__(self):
        """
        EXAMPLES::

            sage: R = cartesian_product([GF(next_prime(2^80)), GF(next_prime(2^81))])  # P[singular matrix] is negligible
            sage: ~matrix.random(R, 100)  # should finish in <2s                        # needs sage.modules
            100 x 100 dense matrix over The Cartesian product of (Finite Field of..., Finite Field of...) (use...)
        """
        return self._from_cartesian_product_of_matrices(~self._to_cartesian_product_of_matrices())
