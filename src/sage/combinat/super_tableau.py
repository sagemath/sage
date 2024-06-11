r"""
Super Tableaux

AUTHORS:

- Matthew Lancellotti (2007): initial version

- Chaman Agrawal (2019-07-23): Modify standard and semistandard tableaux for
  super tableaux.
"""

# ****************************************************************************
#       Copyright (C) 2019 Matthew Lancellotti <mareoraft at gmail.com>
#                     2019 Chaman Agrawal <chaman.ag at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations

from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.rings.integer_ring import ZZ
from sage.sets.family import Family
from sage.structure.parent import Parent
from sage.rings.integer import Integer
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.shifted_primed_tableau import PrimedEntry
from sage.combinat.tableau import (Tableau, Tableaux, SemistandardTableaux,
                                   StandardTableaux)


class SemistandardSuperTableau(Tableau):
    """
    A semistandard super tableau.

    A semistandard super tableau is a tableau with primed positive integer entries.
    As defined in [Muth2019]_, a semistandard super tableau weakly increases along
    the rows and down the columns. Also, the letters of even parity (unprimed)
    strictly increases down the columns, and letters of oddd parity (primed)
    strictly increases along the rows. Note that Sage uses the English convention
    for partitions and tableaux; the longer rows are displayed on top.

    INPUT:

    - ``t`` -- a tableau, a list of iterables, or an empty list

    EXAMPLES::

        sage: t = SemistandardSuperTableau([['1p',2,"3'"],[2,3]]); t
        [[1', 2, 3'], [2, 3]]
        sage: t.shape()
        [3, 2]
        sage: t.pp() # pretty printing
        1' 2 3'
        2 3
        sage: t = Tableau([["1p",2],[2]])
        sage: s = SemistandardSuperTableau(t); s
        [[1', 2], [2]]
        sage: SemistandardSuperTableau([]) # The empty tableau
        []

    TESTS::

        sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
        sage: t = Tableaux()([[1,1],[2]])
        sage: s = SemistandardSuperTableaux()([[PrimedEntry(1),PrimedEntry(1)],
        ....:                                   [PrimedEntry(2)]])
        sage: s == t
        True
        sage: s.parent()
        Semistandard super tableaux
        sage: r = SemistandardSuperTableaux()(t); r.parent()
        Semistandard super tableaux
        sage: isinstance(r, Tableau)
        True
        sage: s2 = SemistandardSuperTableaux()([[PrimedEntry(1),
        ....:       PrimedEntry(1)], [PrimedEntry(2)]])
        sage: s2 == s
        True
        sage: s2.parent()
        Semistandard super tableaux
    """
    @staticmethod
    def __classcall_private__(cls, t):
        r"""
        This ensures that a SemistandardSuperTableau is only ever constructed
        as an element_class call of an appropriate parent.

        TESTS::

            sage: t = SemistandardSuperTableau([[1,1],[2]])
            sage: TestSuite(t).run()
            sage: t.parent()
            Semistandard super tableaux
            sage: t.category()
            Category of elements of Semistandard super tableaux
            sage: type(t)
            <class 'sage.combinat.super_tableau.SemistandardSuperTableaux_all_with_category.element_class'>
        """
        if isinstance(t, cls):
            return t

        # We must verify ``t`` is a list of iterables, and also
        # normalize it to be a list of tuples.
        try:
            t = [tuple(it) for it in t]
        except TypeError:
            raise ValueError("a tableau must be a list of iterables")
        SST = SemistandardSuperTableaux_all()
        return SST.element_class(SST, t)

    def __init__(self, parent, t, check=True, preprocessed=False):
        r"""
        Initialize a semistandard super tableau for given tableau ``t``.

        TESTS::

            sage: s = SemistandardSuperTableau([[1,"2'","3'",3], [2,"3'"]])
            sage: t = SemistandardSuperTableaux()([[1,"2p","3p",3], [2,"3p"]])
            sage: s == t
            True
            sage: s.parent()
            Semistandard super tableaux
            sage: TestSuite(t).run()
            sage: r = SemistandardSuperTableaux()(s); r.parent()
            Semistandard super tableaux
            sage: s is t  # identical super tableaux are distinct objects
            False
        """
        if not preprocessed:
            t = self._preprocess(t)
        super().__init__(parent, t, check=check)

    @staticmethod
    def _preprocess(t):
        """
        Preprocessing list ``T`` to initialize the tableau.

        The output is a list of rows as tuples, with entries being
        ``PrimedEntry`` instances.

        TESTS::

            sage: SemistandardSuperTableau._preprocess([["2'", "3p", 3.5]])
            [[2', 3', 4']]
            sage: SemistandardSuperTableau._preprocess([[None]])
            [[None]]
            sage: SemistandardSuperTableau._preprocess([])
            []
        """
        if isinstance(t, SemistandardSuperTableau):
            return t
        # Preprocessing list t for primes and other symbols
        t = [[PrimedEntry(entry) if entry is not None else entry for entry in row]
             for row in t]
        while t and not t[-1]:
            t = t[:-1]
        return t

    def check(self):
        """
        Check that ``self`` is a valid semistandard super tableau.

        TESTS::

            sage: SemistandardSuperTableau([[1,2,3],[1]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the unprimed entries of each column must
            be strictly increasing

            sage: SemistandardSuperTableau([[1,2,1]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the entries in each row of a semistandard super
            tableau must be weakly increasing

            sage: SemistandardSuperTableau([[0,1]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the entries of a semistandard super tableau must be
            nonnegative primed integers
        """
        super().check()
        for row in self:
            if not all(isinstance(c, PrimedEntry) and c > 0 for c in row):
                raise ValueError("the entries of a semistandard super tableau"
                                 " must be nonnegative primed integers")
            if any(row[c] > row[c + 1] for c in range(len(row) - 1)):
                raise ValueError("the entries in each row of a semistandard"
                                 " super tableau must be weakly increasing")

        if self:
            for row, next in zip(self, self[1:]):
                # Check that letters are weakly increasing down columns
                if any(row[c] > next[c] for c in range(len(next))):
                    raise ValueError("the entries of each column of a "
                                     "semistandard super tableau must be "
                                     "weakly increasing")
                # Check that unprimed letters are column strict
                if not all(row[c] < next[c]
                           for c in range(len(next))
                           if (row[c].is_unprimed() or next[c].is_unprimed())):
                    raise ValueError("the unprimed entries of each column"
                                     " must be strictly increasing")

            # Check that primed letters are row strict
            for row in self:
                if not all(row[c] < row[c + 1]
                           for c in range(len(row) - 1)
                           if (row[c].is_primed() or row[c + 1].is_primed())):
                    raise ValueError("the primed entries in each row must be"
                                     " strictly increasing")


class StandardSuperTableau(SemistandardSuperTableau):
    r"""
    A standard super tableau.

    A standard super tableau is a semistandard super tableau whose entries
    are in bijection with positive primed integers `1', 1, 2' \ldots n`.

    For more information refer [Muth2019]_.

    INPUT:

    - ``t`` -- a Tableau, a list of iterables, or an empty list

    EXAMPLES::

        sage: t = StandardSuperTableau([["1'",1,"2'",2,"3'"],[3,"4'"]]); t
        [[1', 1, 2', 2, 3'], [3, 4']]
        sage: t.shape()
        [5, 2]
        sage: t.pp() # pretty printing
        1' 1 2' 2 3'
        3 4'
        sage: t.is_standard()
        True
        sage: StandardSuperTableau([]) # The empty tableau
        []

    TESTS::

        sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
        sage: t = Tableaux()([[PrimedEntry('1p'), PrimedEntry('2p')],
        ....:                [PrimedEntry(1)]])
        sage: s = StandardSuperTableaux()([['1p','2p'],[1]])
        sage: s == t
        True
        sage: s.parent()
        Standard super tableaux
        sage: r = StandardSuperTableaux()([]); r.parent()
        Standard super tableaux
        sage: isinstance(r, Tableau)
        True
    """
    @staticmethod
    def __classcall_private__(self, t):
        r"""
        This ensures that a :class:`StandardSuperTableau` is only ever
        constructed as an ``element_class`` call of an appropriate parent.

        TESTS::

            sage: t = StandardSuperTableau([['1p','2p'],[1]])
            sage: TestSuite(t).run()
            sage: t.parent()
            Standard super tableaux
            sage: type(t)
            <class 'sage.combinat.super_tableau.StandardSuperTableaux_all_with_category.element_class'>
        """
        if isinstance(t, StandardSuperTableau):
            return t

        SST = StandardSuperTableaux_all()
        return SST.element_class(SST, t)

    def check(self):
        r"""
        Check that ``self`` is a standard tableau.

        TESTS::

            sage: StandardSuperTableau([[1,2,3],[4,5]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the entries in a standard tableau must be in
            bijection with 1',1,2',2,...,n

            sage: StandardSuperTableau([[1,3,2]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the entries in each row of a semistandard super
            tableau must be weakly increasing
        """
        super().check()
        # t is semistandard so we only need to check
        # that its entries are in bijection with {1',1,2', 2, ..., n}
        flattened_list = [i for row in self for i in row]
        a = PrimedEntry('1p')
        primed_list = []
        for i in range(len(flattened_list)):
            primed_list.append(a)
            a = a.increase_half()

        if sorted(flattened_list) != primed_list:
            raise ValueError("the entries in a standard tableau must be in"
                             " bijection with 1',1,2',2,...,n")

    def is_standard(self) -> bool:
        """
        Return ``True`` since ``self`` is a standard super tableau.

        EXAMPLES::

            sage: StandardSuperTableau([['1p', 1], ['2p', 2]]).is_standard()
            True
        """
        return True


################################
# Semi-standard Super tableaux #
################################
class SemistandardSuperTableaux(SemistandardTableaux):
    r"""
    The set of semistandard super tableaux.

    A semistandard super tableau is a tableau with primed positive integer entries.
    As defined in [Muth2019]_, a semistandard super tableau weakly increases along
    the rows and down the columns. Also, the letters of even parity (unprimed)
    strictly increases down the columns, and letters of oddd parity (primed)
    strictly increases along the rows. Note that Sage uses the English convention
    for partitions and tableaux; the longer rows are displayed on top.

    EXAMPLES::

        sage: SST = SemistandardSuperTableaux(); SST
        Semistandard super tableaux
    """
    @staticmethod
    def __classcall_private__(cls):
        r"""
        Normalize and process input to return the correct parent and
        ensure a unique representation.

        TESTS::

            sage: SemistandardSuperTableaux()
            Semistandard super tableaux
        """
        return SemistandardSuperTableaux_all()

    Element = SemistandardSuperTableau

    def __contains__(self, x) -> bool:
        """
        Return ``True`` if ``t`` can be interpreted as a
        :class:`SemistandardSuperTableau`.

        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: T = sage.combinat.super_tableau.SemistandardSuperTableaux_all()
            sage: [[1,2],[2]] in T
            True
            sage: [[PrimedEntry('1p'),PrimedEntry(2)],[PrimedEntry(2)]] in T
            True
            sage: [] in T
            True
            sage: Tableau([[1]]) in T
            True
            sage: StandardSuperTableau([[1]]) in T
            Traceback (most recent call last):
            ...
            ValueError: the entries in a standard tableau must be in bijection
            with 1',1,2',2,...,n

            sage: [[1,2],[1]] in T
            False
            sage: [[1,1],[5]] in T
            True
            sage: [[1,3,2]] in T
            False
        """
        if isinstance(x, SemistandardSuperTableau):
            return True
        elif Tableaux.__contains__(self, x):
            x = SemistandardSuperTableau._preprocess(x)
            for row in x:
                if any(row[c] > row[c + 1] for c in range(len(row) - 1)):
                    return False
                if not all(row[c] < row[c + 1]
                           for c in range(len(row) - 1)
                           if (row[c].is_primed() or row[c + 1].is_primed())):
                    return False
            for row, next in zip(x, x[1:]):
                if any(row[c] > next[c] for c in range(len(next))):
                    return False
                if not all(row[c] < next[c]
                           for c in range(len(next))
                           if (row[c].is_unprimed() or next[c].is_unprimed())):
                    return False
            return True
        else:
            return False


class SemistandardSuperTableaux_all(SemistandardSuperTableaux):
    """
    All semistandard super tableaux.
    """

    def __init__(self):
        r"""
        Initialize the class of all semistandard super tableaux.

        TESTS::

            sage: from sage.combinat.super_tableau import SemistandardSuperTableaux_all
            sage: SST = SemistandardSuperTableaux_all(); SST
            Semistandard super tableaux
        """
        Parent.__init__(self, category=InfiniteEnumeratedSets())
        SemistandardSuperTableaux.__init__(self)

    def _repr_(self) -> str:
        """
        TESTS::

            sage: repr(SemistandardSuperTableaux())
            'Semistandard super tableaux'
        """
        return "Semistandard super tableaux"


################################
# Standard Super tableaux #
################################
class StandardSuperTableaux(SemistandardSuperTableaux, Parent):
    r"""
    The set of standard super tableaux.

    A standard super tableau is a tableau whose entries are primed positive
    integers, which are strictly increasing in rows and down columns and
    contains each letters from 1',1,2'...n exactly once.

    For more information refer [Muth2019]_.

    INPUT:

    - ``n`` -- a nonnegative integer or a partition

    EXAMPLES::

        sage: SST = StandardSuperTableaux()
        sage: SST
        Standard super tableaux
        sage: SST([["1'",1,"2'",2,"3'"],[3,"4'"]])
        [[1', 1, 2', 2, 3'], [3, 4']]
        sage: SST = StandardSuperTableaux(3)
        sage: SST
        Standard super tableaux of size 3
        sage: SST.first()
        [[1', 1, 2']]
        sage: SST.last()
        [[1'], [1], [2']]
        sage: SST.cardinality()
        4
        sage: SST.list()
        [[[1', 1, 2']], [[1', 2'], [1]], [[1', 1], [2']], [[1'], [1], [2']]]
        sage: SST = StandardSuperTableaux([3,2])
        sage: SST
        Standard super tableaux of shape [3, 2]

    TESTS::

        sage: StandardSuperTableaux()([])
        []
        sage: SST = StandardSuperTableaux([3,2]); SST
        Standard super tableaux of shape [3, 2]
        sage: SST.first()
        [[1', 2', 3'], [1, 2]]
        sage: SST.last()
        [[1', 1, 2'], [2, 3']]
        sage: SST.cardinality()
        5
        sage: SST.cardinality() == StandardTableaux([3,2]).cardinality()
        True
        sage: SST.list()
        [[[1', 2', 3'], [1, 2]],
         [[1', 1, 3'], [2', 2]],
         [[1', 2', 2], [1, 3']],
         [[1', 1, 2], [2', 3']],
         [[1', 1, 2'], [2, 3']]]
        sage: TestSuite(SST).run()
    """
    @staticmethod
    def __classcall_private__(cls, n=None):
        r"""
        This class returns the appropriate parent based on arguments.
        See the documentation for :class:`StandardTableaux` for more
        information.

        TESTS::

            sage: SST = StandardSuperTableaux(); SST
            Standard super tableaux
            sage: StandardSuperTableaux(3)
            Standard super tableaux of size 3
            sage: StandardSuperTableaux([2,2])
            Standard super tableaux of shape [2, 2]
            sage: StandardSuperTableaux(-1)
            Traceback (most recent call last):
            ...
            ValueError: the argument must be a nonnegative integer or a
            partition
            sage: StandardSuperTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: the argument must be a nonnegative integer or a
            partition
        """
        from sage.combinat.partition import _Partitions
        from sage.combinat.skew_partition import SkewPartitions

        if n is None:
            return StandardSuperTableaux_all()

        elif n in _Partitions:
            return StandardSuperTableaux_shape(_Partitions(n))

        elif n in SkewPartitions():
            raise NotImplementedError("standard super tableau for skew "
                                      "partitions is not implemented yet")

        if not isinstance(n, (int, Integer)) or n < 0:
            raise ValueError("the argument must be a nonnegative integer"
                             " or a partition")

        return StandardSuperTableaux_size(n)

    Element = StandardSuperTableau

    def __contains__(self, x) -> bool:
        """
        Return ``True`` if ``t`` can be interpreted as a
        :class:`StandardSuperTableau`.

        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: T = sage.combinat.super_tableau.StandardSuperTableaux_all()
            sage: [[0.5,1],[1.5]] in T
            True
            sage: [[PrimedEntry('1p'),PrimedEntry('2p')],[PrimedEntry(1)]] in T
            True
            sage: [] in T
            True
            sage: Tableau([['1p']]) in T
            True
            sage: StandardSuperTableau([['1p']]) in T
            True
            sage: [[1,2],[1]] in T
            False
            sage: [[1,1],[5]] in T
            False
            sage: [[1,3,2]] in T
            False
        """
        if isinstance(x, StandardSuperTableau):
            return True
        elif Tableaux.__contains__(self, x):
            x = SemistandardSuperTableau._preprocess(x)
            flattened_list = [i for row in x for i in row]
            a = PrimedEntry('1p')
            primed_list = []
            for i in range(len(flattened_list)):
                primed_list.append(a)
                a = a.increase_half()
            # return True
            return sorted(flattened_list) == primed_list and (x or
                    (all(row[i] < row[i + 1] for row in x for i in range(len(row) - 1)) and
                     all(x[r][c] < x[r + 1][c] for r in range(len(x) - 1)
                         for c in range(len(x[r + 1])))))
        else:
            return False


class StandardSuperTableaux_all(StandardSuperTableaux,
                                DisjointUnionEnumeratedSets):
    """
    All standard super tableaux.
    """

    def __init__(self):
        r"""
        Initialize the class of all standard super tableaux.

        TESTS::

            sage: from sage.combinat.super_tableau import StandardSuperTableaux_all
            sage: SST = StandardSuperTableaux_all(); SST
            Standard super tableaux
            sage: TestSuite(SST).run()
        """
        StandardSuperTableaux.__init__(self)
        DisjointUnionEnumeratedSets.__init__(self,
                                             Family(NonNegativeIntegers(),
                                                    StandardSuperTableaux_size),
                                             facade=True, keepkey=False)

    def _repr_(self) -> str:
        """
        TESTS::

            sage: repr(StandardSuperTableaux())
            'Standard super tableaux'
        """
        return "Standard super tableaux"


class StandardSuperTableaux_size(StandardSuperTableaux,
                                 DisjointUnionEnumeratedSets):
    """
    Standard super tableaux of fixed size `n`.

    EXAMPLES::

        sage: [ t for t in StandardSuperTableaux(1) ]
        [[[1']]]
        sage: [ t for t in StandardSuperTableaux(2) ]
        [[[1', 1]], [[1'], [1]]]
        sage: [ t for t in StandardSuperTableaux(3) ]
        [[[1', 1, 2']], [[1', 2'], [1]], [[1', 1], [2']], [[1'], [1], [2']]]
        sage: StandardSuperTableaux(4)[:]
        [[[1', 1, 2', 2]],
         [[1', 2', 2], [1]],
         [[1', 1, 2], [2']],
         [[1', 1, 2'], [2]],
         [[1', 2'], [1, 2]],
         [[1', 1], [2', 2]],
         [[1', 2], [1], [2']],
         [[1', 2'], [1], [2]],
         [[1', 1], [2'], [2]],
         [[1'], [1], [2'], [2]]]
    """

    def __init__(self, n):
        r"""
        Initialize the class of all standard super tableaux of size ``n``.

        TESTS::

            sage: TestSuite( StandardSuperTableaux(4) ).run()
        """
        StandardSuperTableaux.__init__(self)
        from sage.combinat.partition import Partitions_n
        DisjointUnionEnumeratedSets.__init__(self,
                                             Family(Partitions_n(n),
                                                    StandardSuperTableaux_shape),
                                             category=FiniteEnumeratedSets(),
                                             facade=True, keepkey=False)
        self.size = Integer(n)

    def _repr_(self) -> str:
        """
        TESTS::

            sage: StandardSuperTableaux(3)
            Standard super tableaux of size 3
        """
        return "Standard super tableaux of size %s" % self.size

    def __contains__(self, x) -> bool:
        """
        TESTS::

            sage: SST3 = StandardSuperTableaux(3)
            sage: all(st in SST3 for st in SST3)
            True
            sage: SST4 = StandardSuperTableaux(4)
            sage: [x for x in SST4 if x in SST3]
            []
            sage: 1 in StandardSuperTableaux(4)
            False
        """
        return (StandardSuperTableaux.__contains__(self, x) and
                sum(map(len, x)) == self.size)

    def cardinality(self):
        r"""
        Return the number of all standard super tableaux of size ``n``.

        The standard super tableaux of size `n` are in bijection with the
        corresponding standard tableaux (under the alphabet relabeling). Refer
        :class:`sage.combinat.tableau.StandardTableaux_size` for more details.

        EXAMPLES::

            sage: StandardSuperTableaux(3).cardinality()
            4
            sage: ns = [1,2,3,4,5,6]
            sage: sts = [StandardSuperTableaux(n) for n in ns]
            sage: all(st.cardinality() == len(st.list()) for st in sts)
            True
            sage: StandardSuperTableaux(50).cardinality()  # long time
            27886995605342342839104615869259776

        TESTS::

            sage: def cardinality_using_hook_formula(n):
            ....:     c = 0
            ....:     for p in Partitions(n):
            ....:         c += StandardSuperTableaux(p).cardinality()
            ....:     return c
            sage: all(cardinality_using_hook_formula(i) ==
            ....:       StandardSuperTableaux(i).cardinality()
            ....:       for i in range(10))
            True
        """
        return StandardTableaux(self.size).cardinality()


class StandardSuperTableaux_shape(StandardSuperTableaux):
    """
    Standard super tableaux of a fixed shape `p`.
    """

    def __init__(self, p):
        r"""
        Initialize the class of all standard super tableaux of a given shape.

        TESTS::

            sage: TestSuite( StandardSuperTableaux([2,2,1]) ).run()
        """
        super().__init__(category=FiniteEnumeratedSets())
        StandardSuperTableaux.__init__(self)
        self.shape = p

    def __contains__(self, x) -> bool:
        """
        EXAMPLES::

            sage: ST = StandardSuperTableaux([2,1,1])
            sage: all(st in ST for st in ST)
            True
            sage: len([x for x in StandardSuperTableaux(4) if x in ST])
            3
            sage: ST.cardinality()
            3
            sage: 1 in StandardSuperTableaux([2,1,1])
            False
        """
        return (StandardSuperTableaux.__contains__(self, x) and
                [len(w) for w in x] == self.shape)

    def _repr_(self) -> str:
        """
        TESTS::

            sage: repr(StandardSuperTableaux([2,1,1]))
            'Standard super tableaux of shape [2, 1, 1]'
        """
        return "Standard super tableaux of shape %s" % str(self.shape)

    def cardinality(self):
        r"""
        Return the number of standard super tableaux of given shape.

        The standard super tableaux of a fixed shape `p` are in bijection with
        the corresponding standard tableaux (under the alphabet relabeling).
        Refer :class:`sage.combinat.tableau.StandardTableaux_shape` for more
        details.

        EXAMPLES::

            sage: StandardSuperTableaux([3,2,1]).cardinality()
            16
            sage: StandardSuperTableaux([2,2]).cardinality()
            2
            sage: StandardSuperTableaux([5]).cardinality()
            1
            sage: StandardSuperTableaux([6,5,5,3]).cardinality()
            6651216
            sage: StandardSuperTableaux([]).cardinality()
            1
        """
        pi = self.shape
        return StandardTableaux(pi).cardinality()

    def __iter__(self):
        r"""
        An iterator for the standard super tableaux associated to the
        shape `p` of ``self``.

        EXAMPLES::

            sage: [t for t in StandardSuperTableaux([2,2])]
            [[[1', 2'], [1, 2]], [[1', 1], [2', 2]]]
            sage: [t for t in StandardSuperTableaux([3,2])]
            [[[1', 2', 3'], [1, 2]],
             [[1', 1, 3'], [2', 2]],
             [[1', 2', 2], [1, 3']],
             [[1', 1, 2], [2', 3']],
             [[1', 1, 2'], [2, 3']]]
            sage: st = StandardSuperTableaux([2,1])
            sage: st[0].parent() is st
            True
        """
        pi = self.shape
        for tableau in StandardTableaux(pi):
            yield self.element_class(self, [[PrimedEntry(ZZ(val) / 2) for val in row]
                                            for row in tableau])
