# sage.doctest: needs sage.combinat sage.modules
"""
Freely Generated Lie Conformal Algebras

AUTHORS:

- Reimundo Heluani (2019-08-09): Initial implementation
"""

#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .lie_conformal_algebra_with_basis import LieConformalAlgebraWithBasis
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.categories.cartesian_product import cartesian_product
from sage.rings.integer import Integer
from sage.sets.family import Family
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets


class FreelyGeneratedLieConformalAlgebra(LieConformalAlgebraWithBasis):
    """
    Base class for a central extension of a freely generated Lie
    conformal algebra.

    This class provides minimal functionality, it sets up the
    family of Lie conformal algebra generators.

    .. NOTE::

        We now only accept direct sums of free modules plus
        some central generators `C_i` such that `TC_i = 0`.
    """
    def __init__(self, R, index_set=None, central_elements=None, category=None,
                 element_class=None, prefix=None, **kwds):
        """
        Initialize ``self``.

        TESTS::

            sage: V = lie_conformal_algebras.Virasoro(QQ)
            sage: TestSuite(V).run()
        """
        self._generators = Family(index_set)
        E = cartesian_product([index_set, NonNegativeIntegers()])
        if central_elements is not None:
            self._generators = DisjointUnionEnumeratedSets([index_set,
                                                    Family(central_elements)])
            E = DisjointUnionEnumeratedSets((cartesian_product([
                Family(central_elements), {Integer(0)}]), E))

        super().__init__(R, basis_keys=E, element_class=element_class,
                         category=category, prefix=prefix, **kwds)

        if central_elements is not None:
            self._central_elements = Family(central_elements)
        else:
            self._central_elements = ()

    def lie_conformal_algebra_generators(self):
        """
        The generators of this Lie conformal algebra.

        OUTPUT: a (possibly infinite) family of generators (as an
        `R[T]`-module) of this Lie conformal algebra.

        EXAMPLES::

            sage: Vir = lie_conformal_algebras.Virasoro(QQ)
            sage: Vir.lie_conformal_algebra_generators()
            (L, C)
            sage: V = lie_conformal_algebras.Affine(QQ,'A1')
            sage: V.lie_conformal_algebra_generators()
            (B[alpha[1]], B[alphacheck[1]], B[-alpha[1]], B['K'])
        """
        F = Family(self._generators,
                   lambda i: self.monomial((i, Integer(0))),
                   name="generator map")
        from sage.categories.sets_cat import Sets
        if F in Sets().Finite():
            return tuple(F)
        return F

    def central_elements(self):
        """
        The central generators of this Lie conformal algebra.

        EXAMPLES::

            sage: Vir = lie_conformal_algebras.Virasoro(QQ)
            sage: Vir.central_elements()
            (C,)
            sage: V = lie_conformal_algebras.Affine(QQ, 'A1')
            sage: V.central_elements()
            (B['K'],)
        """
        return Family(self._central_elements,
                      lambda i: self.monomial((i, Integer(0))),
                      name="central_element map")
