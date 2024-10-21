# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs lrcalc_python
r"""
Littlewood-Richardson tableaux

A semistandard tableau is Littlewood-Richardson with respect to
the sequence of partitions `(\mu^{(1)},\ldots,\mu^{(k)})` if,
when restricted to each alphabet `\{|\mu^{(1)}|+\cdots+|\mu^{(i-1)}|+1,
\ldots, |\mu^{(1)}|+\cdots+|\mu^{(i)}|-1\}`, is Yamanouchi.

AUTHORS:

- Maria Gillespie, Jake Levinson, Anne Schilling (2016): initial version
"""

#*****************************************************************************
#       Copyright (C) 2016 Maria Gillespie
#                          Anne Schilling <anne at math.ucdavis.edu>
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
#                  http://www.gnu.org/licenses/
#****************************************************************************

from itertools import zip_longest

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.tableau import SemistandardTableau, SemistandardTableaux
from sage.combinat.partition import Partition, Partitions


class LittlewoodRichardsonTableau(SemistandardTableau):
    r"""
    A semistandard tableau is Littlewood-Richardson with respect to
    the sequence of partitions `(\mu^{(1)}, \ldots, \mu^{(k)})` if,
    when restricted to each alphabet `\{|\mu^{(1)}|+\cdots+|\mu^{(i-1)}|+1,
    \ldots, |\mu^{(1)}|+\cdots+|\mu^{(i)}|-1\}`, is Yamanouchi.

    INPUT:

    - ``t`` -- Littlewood-Richardson tableau; the input is supposed to be
      a list of lists specifying the rows of the tableau

    EXAMPLES::

        sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableau
        sage: LittlewoodRichardsonTableau([[1,1,3],[2,3],[4]], [[2,1],[2,1]])
        [[1, 1, 3], [2, 3], [4]]
    """
    @staticmethod
    def __classcall_private__(cls, t, weight):
        r"""
        Implement the shortcut ``LittlewoodRichardsonTableau(t, weight)`` to
        ``LittlewoodRichardsonTableaux(shape , weight)(t)``
        where ``shape`` is the shape of the tableau.

        TESTS::

            sage: LR = LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            sage: t = LR([[1, 1, 3], [2, 3], [4]])
            sage: t.check()
            sage: type(t)
            <class 'sage.combinat.lr_tableau.LittlewoodRichardsonTableaux_with_category.element_class'>
            sage: TestSuite(t).run()
            sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableau
            sage: LittlewoodRichardsonTableau([[1,1,3],[2,3],[4]], [[2,1],[2,1]])
            [[1, 1, 3], [2, 3], [4]]
        """
        if isinstance(t, cls):
            return t
        tab = SemistandardTableau(list(t))
        shape = tab.shape()
        return LittlewoodRichardsonTableaux(shape, weight)(t)

    def __init__(self, parent, t):
        r"""
        Initialize ``self``.

        TESTS::

            sage: LR = LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            sage: t = LR([[1, 1, 3], [2, 3], [4]])
            sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableau
            sage: s = LittlewoodRichardsonTableau([[1,1,3],[2,3],[4]], [[2,1],[2,1]])
            sage: s == t
            True
            sage: type(t)
            <class 'sage.combinat.lr_tableau.LittlewoodRichardsonTableaux_with_category.element_class'>
            sage: t.parent()
            Littlewood-Richardson Tableaux of shape [3, 2, 1] and weight ([2, 1], [2, 1])
            sage: TestSuite(t).run()
        """
        self._shape = parent._shape
        self._weight = parent._weight
        super().__init__(parent, list(t))

    def check(self):
        r"""
        Check that ``self`` is a valid Littlewood-Richardson tableau.

        EXAMPLES::

            sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableau
            sage: t = LittlewoodRichardsonTableau([[1,1,3],[2,3],[4]], [[2,1],[2,1]])
            sage: t.check()

        TESTS::

            sage: LR = LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            sage: LR([[1, 1, 2], [3, 3], [4]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 1, 2], [3, 3], [4]] is not an element of
             Littlewood-Richardson Tableaux of shape [3, 2, 1] and weight ([2, 1], [2, 1])
            sage: LR([[1, 1, 2, 3], [3], [4]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 1, 2, 3], [3], [4]] is not an element of
             Littlewood-Richardson Tableaux of shape [3, 2, 1] and weight ([2, 1], [2, 1])
            sage: LR([[1, 1, 3], [3, 3], [4]])
            Traceback (most recent call last):
            ...
            ValueError: weight of the parent does not agree with the weight of the tableau
        """
        super().check()
        if not [i for a in self.parent()._weight for i in a] == self.weight():
            raise ValueError("weight of the parent does not agree "
                             "with the weight of the tableau")
        if not self.shape() == self.parent()._shape:
            raise ValueError("shape of the parent does not agree "
                             "with the shape of the tableau")


class LittlewoodRichardsonTableaux(SemistandardTableaux):
    r"""
    Littlewood-Richardson tableaux.

    A semistandard tableau `t` is *Littlewood-Richardson* with respect to
    the sequence of partitions `(\mu^{(1)}, \ldots, \mu^{(k)})` (called
    the weight) if `t` is Yamanouchi when restricted to each alphabet
    `\{|\mu^{(1)}| + \cdots + |\mu^{(i-1)}| + 1, \ldots,
    |\mu^{(1)}| + \cdots + |\mu^{(i)}| - 1\}`.

    INPUT:

    - ``shape`` -- the shape of the Littlewood-Richardson tableaux
    - ``weight`` -- the weight is a sequence of partitions

    EXAMPLES::

        sage: LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
        Littlewood-Richardson Tableaux of shape [3, 2, 1] and weight ([2, 1], [2, 1])
    """
    @staticmethod
    def __classcall_private__(cls, shape, weight):
        r"""
        Straighten arguments before unique representation.

        TESTS::

            sage: LR = LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            sage: TestSuite(LR).run()
            sage: LittlewoodRichardsonTableaux([3,2,1],[[2,1]])
            Traceback (most recent call last):
            ...
            ValueError: the sizes of shapes and sequence of weights do not match
        """
        shape = Partition(shape)
        weight = tuple(Partition(a) for a in weight)
        if shape.size() != sum(a.size() for a in weight):
            raise ValueError("the sizes of shapes and sequence of weights do not match")
        return super().__classcall__(cls, shape, weight)

    def __init__(self, shape, weight):
        r"""
        Initialize the parent class of Littlewood-Richardson tableaux.

        INPUT:

        - ``shape`` -- the shape of the Littlewood-Richardson tableaux
        - ``weight`` -- the weight is a sequence of partitions

        TESTS::

            sage: LR = LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            sage: TestSuite(LR).run()
        """
        self._shape = shape
        self._weight = weight
        self._heights = [a.length() for a in self._weight]
        super().__init__(category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            Littlewood-Richardson Tableaux of shape [3, 2, 1] and weight ([2, 1], [2, 1])
        """
        return "Littlewood-Richardson Tableaux of shape %s and weight %s" % (self._shape, self._weight)

    def __iter__(self):
        r"""
        TESTS::

            sage: LR = LittlewoodRichardsonTableaux([3,2,1], [[2,1],[2,1]])
            sage: LR.list()
            [[[1, 1, 3], [2, 3], [4]], [[1, 1, 3], [2, 4], [3]]]
        """
        from sage.libs.lrcalc.lrcalc import lrskew
        if not self._weight:
            yield self.element_class(self, [])
            return

        for nu in Partitions(self._shape.size() - self._weight[-1].size(),
                             outer=self._shape):
            for s in lrskew(self._shape, nu, weight=self._weight[-1]):
                for t in LittlewoodRichardsonTableaux(nu, self._weight[:-1]):
                    shift = sum(a.length() for a in self._weight[:-1])
                    yield self.element_class(self, _tableau_join(t, s, shift=shift))

    def __contains__(self, t):
        """
        Check if ``t`` is contained in ``self``.

        TESTS::

            sage: LR = LittlewoodRichardsonTableaux([3,2,1], [[2,1],[2,1]])
            sage: SST = SemistandardTableaux([3,2,1], [2,1,2,1])
            sage: [t for t in SST if t in LR]
            [[[1, 1, 3], [2, 3], [4]], [[1, 1, 3], [2, 4], [3]]]
            sage: [t for t in SST if t in LR] == LR.list()
            True

            sage: LR = LittlewoodRichardsonTableaux([3,2,1], [[2,1],[2,1]])
            sage: T = [[1,1,3], [2,3], [4]]
            sage: T in LR
            True
        """
        return (SemistandardTableaux.__contains__(self, t)
                and is_littlewood_richardson(t, self._heights))

    Element = LittlewoodRichardsonTableau

#### common or global functions related to LR tableaux


def is_littlewood_richardson(t, heights):
    r"""
    Return whether semistandard tableau ``t`` is Littleword-Richardson
    with respect to ``heights``.

    A tableau is Littlewood-Richardson with respect to ``heights`` given
    by `(h_1, h_2, \ldots)` if each subtableau with respect to the
    alphabets `\{1, 2, \ldots, h_1\}`, `\{h_1+1, \ldots, h_1+h_2\}`,
    etc. is Yamanouchi.

    EXAMPLES::

        sage: from sage.combinat.lr_tableau import is_littlewood_richardson
        sage: t = Tableau([[1,1,2,3,4],[2,3,3],[3]])
        sage: is_littlewood_richardson(t,[2,2])
        False
        sage: t = Tableau([[1,1,3],[2,3],[4,4]])
        sage: is_littlewood_richardson(t,[2,2])
        True
        sage: t = Tableau([[7],[8]])
        sage: is_littlewood_richardson(t,[2,3,3])
        False
        sage: is_littlewood_richardson([[2],[3]],[3,3])
        False
    """
    from sage.combinat.words.word import Word
    partial = [sum(heights[i] for i in range(j)) for j in range(len(heights)+1)]
    try:
        w = t.to_word()
    except AttributeError:  # Not an instance of Tableau
        w = sum(reversed(t), [])
    for i in range(len(heights)):
        subword = Word([j for j in w if partial[i]+1 <= j <= partial[i+1]],
                       alphabet=list(range(partial[i]+1,partial[i+1]+1)))
        if not subword.is_yamanouchi():
            return False
    return True


def _tableau_join(t1, t2, shift=0):
    """
    Join semistandard tableau ``t1`` with semistandard tableau ``t2``
    shifted by ``shift``.

    Concatenate the rows of ``t1`` and ``t2``, dropping any ``None``'s
    from ``t2``. This method is intended for the case when the outer
    shape of ``t1`` is equal to the inner shape of ``t2``.

    EXAMPLES::

        sage: from sage.combinat.lr_tableau import _tableau_join
        sage: _tableau_join([[1,2]],[[None,None,2],[3]],shift=5)
        [[1, 2, 7], [8]]
    """
    return [list(row1) + [e2+shift for e2 in row2 if e2 is not None]
            for (row1, row2) in zip_longest(t1, t2, fillvalue=[])]
