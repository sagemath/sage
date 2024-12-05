r"""
Rings of formal Dirichlet series

For any ring `R` and any positive integer `n`, we obtain a ring of formal Dirichlet
series over `R` truncated to precision `n` by considering the formal expressions of
the form `\sum_{i=1}^{n-1} a_i i^{-s}` where the `a_i` are elements of `R`, with the
addition and multiplication rules implied by the syntax.

Sage provides dense and sparse implementations of fixed-precision Dirichlet series
over any Sage base ring.

EXAMPLES::

    sage: R = DirichletSeriesRing(ZZ, 10)
    sage: u = R([1,3,1]); u
    1 + 3*2^-s + 3^-s + O(10^-s)
    sage: v = 1/u; v
    1 - 3*2^-s - 3^-s + 9*4^-s + 6*6^-s - 27*8^-s + 9^-s + O(10^-s)
    sage: u*v
    1 + O(10^-s)

AUTHORS:

- Kiran Kedlaya (2024-08-01): initial version

"""

# ****************************************************************************
#       Copyright (C) 2024 Kiran S. Kedlaya <kedlaya@ucsd.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.misc_c import prod
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

    def _repr_(self) -> str:
        """
        Return a string representation.

        EXAMPLES::

            sage: R = DirichletSeriesRing(ZZ, 10)
            sage: R
            Dirichlet Series Ring over Integer Ring with fixed precision 10
        """
        return "Dirichlet Series Ring over {} with fixed precision {}".format(self.base_ring(), self.__precision)

    def is_sparse(self) -> bool:
        """
        Return ``True`` if this ring uses sparse internal representation.

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

    def euler_product(self, Lpoly):
        """
        Construct an Euler product out of `L`-polynomials.

        EXAMPLES:

        Construct the Riemann zeta function as a formal Dirichlet series::

             sage: R = DirichletSeriesRing(ZZ, 10)
             sage: P.<T> = ZZ[]
             sage: R.euler_product({p: 1-T for p in prime_range(25)})
             1 + 2^-s + 3^-s + 4^-s + 5^-s + 6^-s + 7^-s + 8^-s + 9^-s + O(10^-s)
        """
        return prod(1 / self({p**i: j for i, j in Lpoly[p].dict().items()})
            for p in prime_range(self.precision()))
