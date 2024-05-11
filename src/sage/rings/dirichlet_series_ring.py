from sage.structure.parent import Parent
from sage.rings.fast_arith import prime_range
from sage.rings.ring import CommutativeRing
from sage.rings.dirichlet_series_ring_element import DirichletSeries_dense, DirichletSeries_sparse

class DirichletSeriesRing(CommutativeRing, Parent):
    """
    A ring of Dirichlet series over a base ring, truncated to a fixed precision.
    
    EXAMPLES::

        sage: R = DirichletSeriesRing(ZZ, 10)
        sage: u = R([1,3,1]); u
        1 + 3*2^-s + 3^-s + O(10^-s)
        sage: v = 1/u; v
        1 - 3*2^-s - 3^-s + 9*4^-s + 6*6^-s - 27*8^-s + 9^-s + O(10^-s)
        sage: u*v
        1 + O(10^-s)
    """
    def __init__(self, base_ring, precision, sparse=True):
        """
        Create a Dirichlet series ring.
        
        EXAMPLES::

            sage: R = DirichletSeriesRing(ZZ, 10, sparse=False)
            sage: S = DirichletSeriesRing(ZZ, 10, sparse=True)
            sage: R.has_coerce_map_from(ZZ)
            True
            sage: R.has_coerce_map_from(S)
            False
            sage: S.has_coerce_map_from(R)
            True
        """
        self.Element = DirichletSeries_sparse if sparse else DirichletSeries_dense
        CommutativeRing.__init__(self, base_ring, names=None, category=base_ring.category())
        Parent.__init__(self, base_ring, names=None, category=base_ring.category())
        self.__precision = precision
        self.__is_sparse = sparse
        
    def is_sparse(self):
        """
        Return `True` if this ring uses sparse internal representation.
 
        EXAMPLES::
        
            sage: R = DirichletSeriesRing(ZZ, 10, sparse=False)
            sage: S = DirichletSeriesRing(ZZ, 10, sparse=True)
            sage: R.is_sparse()
            False
            sage: S.is_sparse()
            True
        """
        return self.__is_sparse
        
    def precision(self):
        """
        Return the specified precision for this ring.

        EXAMPLES::

            sage: R = DirichletSeriesRing(ZZ, 10)
            sage: R.precision()
            10
        """
        return self.__precision
        
    def _coerce_map_from_(self, S):
        """
        Implement coercion.

        EXAMPLES::
        
            sage: R = DirichletSeriesRing(ZZ, 10, sparse=False)
            sage: S = DirichletSeriesRing(ZZ, 10, sparse=True)
            sage: R.has_coerce_map_from(ZZ)
            True
            sage: R.has_coerce_map_from(S)
            False
            sage: S.has_coerce_map_from(R)
            True
        """
        base_ring = self.base_ring()
        if base_ring.has_coerce_map_from(S):
            return True
        if isinstance(S, DirichletSeriesRing) and base_ring.has_coerce_map_from(S.base_ring()) and self.precision() <= S.precision() and (self.is_sparse() is True or S.is_sparse() is False):
            return True

    def euler_product(self, factor):
        """
        EXAMPLES:

        Construct the Riemann zeta function as a formal Dirichlet series::

             sage: R = DirichletSeriesRing(ZZ, 10)
             sage: R.euler_product(lambda p: 1 - R({p:1}))
             1 + 2^-s + 3^-s + 4^-s + 5^-s + 6^-s + 7^-s + 8^-s + 9^-s + O(10^-s)
             sage: R.euler_product({p: 1 - R({p:1}) for p in prime_range(25)})
             1 + 2^-s + 3^-s + 4^-s + 5^-s + 6^-s + 7^-s + 8^-s + 9^-s + O(10^-s)
        """
        ans = self.one()
        if isinstance(factor, dict):
            for p in prime_range(self.precision()):
                ans /= factor[p]
        else:
            for p in prime_range(self.precision()):
                ans /= factor(p)
        return ans

