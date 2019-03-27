"""
Benkart-Kang-Kashiwara crystals for the general-linear Lie superalgebra
"""

#*****************************************************************************
#       Copyright (C) 2017 Franco Saliola <saliola@gmail.com>
#                     2017 Travis Scrimshaw <tcscrims at gmail.com>
#                     2017 Anne Schilling <anne@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.categories.regular_supercrystals import RegularSuperCrystals
from sage.combinat.partition import _Partitions
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.crystals.letters import CrystalOfBKKLetters
from sage.combinat.crystals.tensor_product import CrystalOfWords
from sage.combinat.crystals.tensor_product_element import CrystalOfBKKTableauxElement


class CrystalOfBKKTableaux(CrystalOfWords):
    r"""
    Crystal of tableaux for type `A(m|n)`.

    This is an implementation of the tableaux model of the
    Benkart-Kang-Kashiwara crystal [BKK2000]_ for the Lie
    superalgebra `\mathfrak{gl}(m+1,n+1)`.

    INPUT:

    - ``ct`` -- a super Lie Cartan type of type `A(m|n)`
    - ``shape`` -- shape specifying the highest weight; this should be
      a partition contained in a hook of height `n+1` and width `m+1`

    EXAMPLES::

        sage: T = crystals.Tableaux(['A', [1,1]], shape = [2,1])
        sage: T.cardinality()
        20
    """
    @staticmethod
    def __classcall_private__(cls, ct, shape):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: crystals.Tableaux(['A', [1, 2]], shape=[2,1])
            Crystal of BKK tableaux of shape [2, 1] of gl(2|3)
            sage: crystals.Tableaux(['A', [1, 1]], shape=[3,3,3])
            Traceback (most recent call last):
            ...
            ValueError: invalid hook shape
        """
        ct = CartanType(ct)
        shape = _Partitions(shape)
        if len(shape) > ct.m + 1 and shape[ct.m+1] > ct.n + 1:
            raise ValueError("invalid hook shape")
        return super(CrystalOfBKKTableaux, cls).__classcall__(cls, ct, shape)

    def __init__(self, ct, shape):
        r"""
        Initialize ``self``.

        TESTS::

            sage: T = crystals.Tableaux(['A', [1,1]], shape = [2,1]); T
            Crystal of BKK tableaux of shape [2, 1] of gl(2|2)
            sage: TestSuite(T).run()
        """
        self._shape = shape
        self._cartan_type = ct
        m = ct.m + 1
        C = CrystalOfBKKLetters(ct)
        tr = shape.conjugate()
        mg = []
        for i,col_len in enumerate(tr):
            for j in range(col_len - m):
                mg.append(C(i+1))
            for j in range(max(0, m - col_len), m):
                mg.append(C(-j-1))
        mg = list(reversed(mg))
        Parent.__init__(self, category=RegularSuperCrystals())
        self.module_generators = (self.element_class(self, mg),)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: crystals.Tableaux(['A', [1, 2]], shape=[2,1])
            Crystal of BKK tableaux of shape [2, 1] of gl(2|3)
        """
        m = self._cartan_type.m + 1
        n = self._cartan_type.n + 1
        return "Crystal of BKK tableaux of shape {} of gl({}|{})".format(self.shape(), m, n)

    def shape(self):
        r"""
        Return the shape of ``self``.

        EXAMPLES::

            sage: T = crystals.Tableaux(['A', [1, 2]], shape=[2,1])
            sage: T.shape()
            [2, 1]
        """
        return self._shape

    def genuine_highest_weight_vectors(self, index_set=None):
        """
        Return a tuple of genuine highest weight elements.

        A *fake highest weight vector* is one which is annihilated by
        `e_i` for all `i` in the index set, but whose weight is not
        bigger in dominance order than all other elements in the
        crystal. A *genuine highest weight vector* is a highest
        weight element that is not fake.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A', [1,1]], shape=[3,2,1])
            sage: B.genuine_highest_weight_vectors()
            ([[-2, -2, -2], [-1, -1], [1]],)
            sage: B.highest_weight_vectors()
            ([[-2, -2, -2], [-1, -1], [1]],
             [[-2, -2, -2], [-1, 2], [1]],
             [[-2, -2, 2], [-1, -1], [1]])
        """
        if index_set is None or index_set == self.index_set():
            return self.module_generators
        return super(CrystalOfBKKTableaux, self).genuine_highest_weight_vectors(index_set)

    class Element(CrystalOfBKKTableauxElement):
        pass

