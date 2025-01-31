r"""
Ribbon Shaped Tableaux
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.skew_tableau import SkewTableau, SkewTableaux, StandardSkewTableaux
from sage.combinat.tableau import Tableaux
from sage.combinat.permutation import descents_composition_first, descents_composition_list, descents_composition_last
from sage.rings.integer import Integer
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.sets_cat import Sets


class RibbonShapedTableau(SkewTableau):
    r"""
    A ribbon shaped tableau.

    For the purposes of this class, a ribbon shaped tableau is a skew
    tableau whose shape is a skew partition which:

    - has at least one cell in row `1`;

    - has at least one cell in column `1`;

    - has exactly one cell in each of `q` consecutive diagonals, for
      some nonnegative integer `q`.

    A ribbon is given by a list of the rows from top to bottom.

    EXAMPLES::

        sage: x = RibbonShapedTableau([[None, None, None, 2, 3], [None, 1, 4, 5], [3, 2]]); x
        [[None, None, None, 2, 3], [None, 1, 4, 5], [3, 2]]
        sage: x.pp()
          .  .  .  2  3
          .  1  4  5
          3  2
        sage: x.shape()
        [5, 4, 2] / [3, 1]

    The entries labeled by ``None`` correspond to the inner partition.
    Using ``None`` is optional; the entries will be shifted accordingly.  ::

        sage: x = RibbonShapedTableau([[2,3],[1,4,5],[3,2]]); x.pp()
          .  .  .  2  3
          .  1  4  5
          3  2

    TESTS::

        sage: r = RibbonShapedTableau([[1], [2,3], [4, 5, 6]])
        sage: r.to_permutation()
        [4, 5, 6, 2, 3, 1]

        sage: RibbonShapedTableau([[1,2],[3,4]]).evaluation()
        [1, 1, 1, 1]
    """
    @staticmethod
    def __classcall_private__(cls, rows):
        r"""
        Return a ribbon shaped tableau object.

        EXAMPLES::

            sage: RibbonShapedTableau([[2,3],[1,4,5]])
            [[None, None, 2, 3], [1, 4, 5]]

        TESTS::

            sage: RibbonShapedTableau([4,5])
            Traceback (most recent call last):
            ...
            TypeError: rows must be lists of positive integers

            sage: RibbonShapedTableau([[2,3],[-4,5]])
            Traceback (most recent call last):
            ...
            TypeError: r must be a list of positive integers
        """
        try:
            r = [tuple(r) for r in rows]
        except TypeError:
            raise TypeError("rows must be lists of positive integers")
        if not r:
            return StandardRibbonShapedTableaux()(r)
        if all(j is None or (isinstance(j, (int, Integer)) and j > 0)
               for i in r for j in i):
            return StandardRibbonShapedTableaux()(r)
        raise TypeError("r must be a list of positive integers")

    def __init__(self, parent, t):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: R = RibbonShapedTableau([[2,3],[1,4,5]])
            sage: TestSuite(R).run()
        """
        if not isinstance(t, SkewTableau):
            # scrubbing None
            t = [[i for i in row if i is not None] for row in t]

            st = []
            space_count = 0
            for row in reversed(t):
                st.append([None] * space_count + row)
                space_count += len(row) - 1
            st.reverse()
            t = st
        else:
            t = list(t)

        SkewTableau.__init__(self, parent, t)

    def height(self):
        """
        Return the height of ``self``.

        The height is given by the number of rows in the outer partition.

        EXAMPLES::

            sage: RibbonShapedTableau([[2,3],[1,4,5]]).height()
            2
        """
        return len(self.outer_shape())

    def spin(self):
        """
        Return the spin of ``self``.

        EXAMPLES::

            sage: RibbonShapedTableau([[2,3],[1,4,5]]).spin()
            1/2
        """
        return Integer(self.height() - 1) / 2

    def width(self):
        """
        Return the width of the ribbon.

        This is given by the length of the longest row in the outer partition.

        EXAMPLES::

            sage: RibbonShapedTableau([[2,3],[1,4,5]]).width()
            4
            sage: RibbonShapedTableau([]).width()
            0
        """
        return len(self[0]) if self else 0


class RibbonShapedTableaux(SkewTableaux):
    """
    The set of all ribbon shaped tableaux.
    """
    @staticmethod
    def __classcall_private__(cls, shape=None, **kwds):
        """
        Normalize input to ensure a unique representation and pick the correct
        class based on input.

        The ``shape`` parameter is currently ignored.

        EXAMPLES::

            sage: S1 = RibbonShapedTableaux([4, 2, 2, 1])
            sage: S2 = RibbonShapedTableaux((4, 2, 2, 1))
            sage: S1 is S2
            True
        """
        # if shape is not None:
        #     from sage.combinat.partition import Partition
        #     return RibbonShapedTableaux_shape(Partition(shape))

        # Otherwise arg0 takes the place of the category in pickling
        return super().__classcall__(cls, **kwds)

    def __init__(self, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = RibbonShapedTableaux()                                            # needs sage.graphs
            sage: TestSuite(S).run()                                                    # needs sage.graphs
        """
        if category is None:
            category = Sets()

        SkewTableaux.__init__(self, category=category)

    def _repr_(self):
        """
        TESTS::

            sage: repr(RibbonShapedTableaux())    # indirect doctest
            'Ribbon shaped tableaux'
        """
        return "Ribbon shaped tableaux"

    Element = RibbonShapedTableau
    options = Tableaux.options

    def from_shape_and_word(self, shape, word):
        """
        Return the ribbon corresponding to the given ribbon shape and word.

        EXAMPLES::

            sage: RibbonShapedTableaux().from_shape_and_word([1,3],[1,3,3,7])
            [[None, None, 1], [3, 3, 7]]
        """
        pos = 0
        r = []
        for l in shape:
            r.append(word[pos:pos + l])
            pos += l
        return self.element_class(self, r)


class StandardRibbonShapedTableaux(StandardSkewTableaux):
    """
    The set of all standard ribbon shaped tableaux.

    INPUT:

    - ``shape`` -- (optional) the composition shape of the rows
    """
    @staticmethod
    def __classcall_private__(cls, shape=None, **kwds):
        """
        Normalize input to ensure a unique representation and pick the correct
        class based on input.

        EXAMPLES::

            sage: S1 = StandardRibbonShapedTableaux([4, 2, 2, 1])
            sage: S2 = StandardRibbonShapedTableaux((4, 2, 2, 1))
            sage: S1 is S2
            True
        """
        if shape is not None:
            from sage.combinat.partition import Partition
            return StandardRibbonShapedTableaux_shape(Partition(shape))

        # Otherwise arg0 takes the place of the category in pickling
        return super().__classcall__(cls, **kwds)

    def __init__(self, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: # needs sage.graphs sage.modules
            sage: S = StandardRibbonShapedTableaux()
            sage: TestSuite(S).run()
        """
        if category is None:
            category = InfiniteEnumeratedSets()

        StandardSkewTableaux.__init__(self, category=category)

    def _repr_(self):
        """
        TESTS::

            sage: repr(StandardRibbonShapedTableaux())    # indirect doctest
            'Standard ribbon shaped tableaux'
        """
        return "Standard ribbon shaped tableaux"

    def __iter__(self):
        """
        Iterate through ``self``.

        EXAMPLES::

            sage: # needs sage.graphs sage.modules
            sage: it = StandardRibbonShapedTableaux().__iter__()
            sage: [next(it) for x in range(10)]
            [[],
             [[1]],
             [[1, 2]],
             [[1], [2]],
             [[1, 2, 3]],
             [[None, 1], [2, 3]],
             [[None, 2], [1, 3]],
             [[1], [2], [3]],
             [[1, 2, 3, 4]],
             [[None, None, 1], [2, 3, 4]]]
        """
        from sage.combinat.partition import _Partitions
        for p in _Partitions:
            for r in StandardRibbonShapedTableaux_shape(p):
                yield self.element_class(self, r)

    Element = RibbonShapedTableau
    options = Tableaux.options

    def from_shape_and_word(self, shape, word):
        """
        Return the ribbon corresponding to the given ribbon shape and word.

        EXAMPLES::

            sage: StandardRibbonShapedTableaux().from_shape_and_word([2,3],[1,2,3,4,5])
            [[None, None, 1, 2], [3, 4, 5]]
        """
        pos = 0
        r = []
        for l in shape:
            r.append(word[pos:pos + l])
            pos += l
        return self.element_class(self, r)

    def from_permutation(self, p):
        """
        Return a standard ribbon of size ``len(p)`` from a permutation ``p``. The
        lengths of each row are given by the distance between the descents
        of the permutation ``p``.

        EXAMPLES::

            sage: import sage.combinat.ribbon_shaped_tableau as rst
            sage: [StandardRibbonShapedTableaux().from_permutation(p)
            ....:  for p in Permutations(3)]
            [[[1, 2, 3]],
             [[None, 2], [1, 3]],
             [[1, 3], [2]],
             [[None, 1], [2, 3]],
             [[1, 2], [3]],
             [[1], [2], [3]]]
        """
        if not p:
            return self.element_class(self, [])

        comp = p.descents()

        if not comp:
            return self.element_class(self, [p[:]])

        r = [[p[j] for j in range(comp[0])]]
        r.extend([p[j] for j in range(comp[i], comp[i + 1])]
                 for i in range(len(comp) - 1))
        r.append([p[j] for j in range(comp[-1], len(p))])
        r.reverse()
        return self.element_class(self, r)


class StandardRibbonShapedTableaux_shape(StandardRibbonShapedTableaux):
    """
    Class of standard ribbon shaped tableaux of ribbon shape ``shape``.

    EXAMPLES::

        sage: StandardRibbonShapedTableaux([2,2])
        Standard ribbon shaped tableaux of shape [2, 2]
        sage: StandardRibbonShapedTableaux([2,2]).first()
        [[None, 2, 4], [1, 3]]
        sage: StandardRibbonShapedTableaux([2,2]).last()
        [[None, 1, 2], [3, 4]]

        sage: # needs sage.graphs sage.modules
        sage: StandardRibbonShapedTableaux([2,2]).cardinality()
        5
        sage: StandardRibbonShapedTableaux([2,2]).list()
        [[[None, 1, 3], [2, 4]],
         [[None, 1, 2], [3, 4]],
         [[None, 2, 3], [1, 4]],
         [[None, 2, 4], [1, 3]],
         [[None, 1, 4], [2, 3]]]
        sage: StandardRibbonShapedTableaux([3,2,2]).cardinality()
        155
    """
    @staticmethod
    def __classcall_private__(cls, shape):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: S = StandardRibbonShapedTableaux([2,2])
            sage: S2 = StandardRibbonShapedTableaux((2,2))
            sage: S is S2
            True
        """
        return super(StandardRibbonShapedTableaux, cls).__classcall__(cls, tuple(shape))

    def __init__(self, shape):
        """
        TESTS::

            sage: S = StandardRibbonShapedTableaux([2,2])
            sage: TestSuite(S).run()                                                    # needs sage.graphs
        """
        self.shape = shape
        StandardRibbonShapedTableaux.__init__(self, FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: StandardRibbonShapedTableaux([2,2])
            Standard ribbon shaped tableaux of shape [2, 2]
        """
        return "Standard ribbon shaped tableaux of shape %s" % list(self.shape)

    def first(self):
        """
        Return the first standard ribbon of ``self``.

        EXAMPLES::

            sage: StandardRibbonShapedTableaux([2,2]).first()
            [[None, 2, 4], [1, 3]]
        """
        return self.from_permutation(descents_composition_first(self.shape))

    def last(self):
        """
        Return the last standard ribbon of ``self``.

        EXAMPLES::

            sage: StandardRibbonShapedTableaux([2,2]).last()
            [[None, 1, 2], [3, 4]]
        """
        return self.from_permutation(descents_composition_last(self.shape))

    def __iter__(self):
        """
        An iterator for the standard ribbon of ``self``.

        EXAMPLES::

            sage: [t for t in StandardRibbonShapedTableaux([2,2])]                      # needs sage.graphs
            [[[None, 1, 3], [2, 4]],
             [[None, 1, 2], [3, 4]],
             [[None, 2, 3], [1, 4]],
             [[None, 2, 4], [1, 3]],
             [[None, 1, 4], [2, 3]]]
        """
        for p in descents_composition_list(self.shape):
            yield self.from_permutation(p)


class Ribbon_class(RibbonShapedTableau):
    """
    This exists solely for unpickling ``Ribbon_class`` objects.
    """

    def __setstate__(self, state):
        r"""
        Unpickle old ``Ribbon_class`` objects.

        EXAMPLES::

            sage: loads(b'x\x9ck`J.NLO\xd5K\xce\xcfM\xca\xccK,\xd1+\xcaLJ\xca\xcf\xe3\n\x02S\xf1\xc99\x89\xc5\xc5\\\x85\x8c\x9a\x8d\x85L\xb5\x85\xcc\x1a\xa1\xac\xf1\x19\x89\xc5\x19\x85,~@VNfqI!kl!\x9bFl!\xbb\x06\xc4\x9c\xa2\xcc\xbc\xf4b\xbd\xcc\xbc\x92\xd4\xf4\xd4"\xae\xdc\xc4\xec\xd4x\x18\xa7\x90#\x94\xd1\xb05\xa8\x903\x03\xc80\x022\xb8Rc\x0b\xb95@<c \x8f\x07\xc40\x012xSSK\x93\xf4\x00l\x811\x17')
            [[None, 1, 2], [3, 4]]
            sage: loads(dumps( RibbonShapedTableau([[3,2,1], [1,1]]) ))  # indirect doctest
            [[None, 3, 2, 1], [1, 1]]
        """
        self.__class__ = RibbonShapedTableau
        self.__init__(RibbonShapedTableaux(), state['_list'])


from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.ribbon', 'Ribbon_class', Ribbon_class)
register_unpickle_override('sage.combinat.ribbon', 'StandardRibbons_shape', StandardRibbonShapedTableaux)
