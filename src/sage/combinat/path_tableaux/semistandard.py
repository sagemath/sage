# sage_setup: distribution = sagemath-combinat
r"""
Semistandard Tableaux
=====================

This is an implementation of the abstract base class
:class:`sage.combinat.path_tableaux.path_tableau.PathTableau`.

This implementation is for semistandard tableaux, represented as a chain of partitions
(essentially, the Gelfand-Tsetlin pattern).
This generalises the jeu de taquin operations of rectification, promotion, evacuation from
standard tableaux to semistandard tableaux. The local rule is the Bender-Knuth involution.

EXAMPLES::

    sage: pt = path_tableaux.SemistandardPathTableau([[], [3], [3,2], [3,3,1],
    ....:                                             [3,3,2,1], [4,3,3,1,0]])
    sage: pt.promotion()
    [(), (2,), (3, 1), (3, 2, 1), (4, 3, 1, 0), (4, 3, 3, 1, 0)]
    sage: pt.evacuation()
    [(), (2,), (4, 0), (4, 2, 0), (4, 3, 1, 0), (4, 3, 3, 1, 0)]

    sage: pt = path_tableaux.SemistandardPathTableau([[], [3], [3,2], [3,3,1],
    ....:                                             [3,3,2,1], [9/2,3,3,1,0]])
    sage: pt.promotion()
    [(), (2,), (3, 1), (3, 2, 1), (9/2, 3, 1, 0), (9/2, 3, 3, 1, 0)]
    sage: pt.evacuation()
    [(), (5/2,), (9/2, 0), (9/2, 2, 0), (9/2, 3, 1, 0), (9/2, 3, 3, 1, 0)]

    sage: pt = path_tableaux.SemistandardPathTableau([[], [3], [4,2], [5,4,1]])
    sage: path_tableaux.CylindricalDiagram(pt)
    [       (),      (3,),    (4, 2), (5, 4, 1)]
    [         ,        (),      (3,),    (5, 2), (5, 4, 1)]
    [         ,          ,        (),      (4,),    (4, 3), (5, 4, 1)]
    [         ,          ,          ,        (),      (3,),    (5, 1), (5, 4, 1)]

    sage: pt2 = path_tableaux.SemistandardPathTableau([[3,2], [3,3,1],
    ....:                                              [3,3,2,1], [4,3,3,1,0]])
    sage: pt1 = path_tableaux.SemistandardPathTableau([[], [3], [3,2]])
    sage: pt1.commutor(pt2)
    ([(), (2,), (2, 2), (4, 2, 0)], [(4, 2, 0), (4, 3, 2, 0), (4, 3, 3, 1, 0)])
    sage: pt1.commutor(pt2,verbose=True)
    [(3, 2), (3, 3, 1), (3, 3, 2, 1), (4, 3, 3, 1, 0)]
    [(3,), (3, 2), (3, 2, 2), (4, 3, 2, 0)]
    [(), (2,), (2, 2), (4, 2, 0)]
    ([(), (2,), (2, 2), (4, 2, 0)], [(4, 2, 0), (4, 3, 2, 0), (4, 3, 3, 1, 0)])

    sage: st = SkewTableau([[None, None, None, 4, 4, 5, 6, 7], [None, 2, 4, 6, 7, 7, 7],
    ....:                   [None, 4, 5, 8, 8, 9], [None, 6, 7, 10], [None, 8, 8, 11],
    ....:                   [None], [4]])
    sage: pt = path_tableaux.SemistandardPathTableau(st)
    sage: bk = [SkewTableau(st.bender_knuth_involution(i+1)) for i in range(10)]
    sage: lr = [pt.local_rule(i+1) for i in range(10)]
    sage: [r.to_tableau() for r in lr] == bk
    True

TESTS::

    sage: pt = path_tableaux.SemistandardPathTableau([[], [3], [3,2], [3,3,1], [3,3,2,1], [4,3,3,1,0]])
    sage: TestSuite(pt).run()

    sage: pt = path_tableaux.SemistandardPathTableau([[], [3], [3,2], [3,3,1], [7/2,3,2,1], [4,3,3,1,0]])
    sage: TestSuite(pt).run()
    Failure in _test_jdt_promotion:
    Traceback (most recent call last):
    ...
    The following tests failed: _test_jdt_promotion

    sage: pt = path_tableaux.SemistandardPathTableau([[3,2], [3,3,1], [3,3,2,1], [4,3,3,1,0]])
    sage: pt.promotion()
    [(3, 2), (3, 2, 2), (4, 3, 2, 0), (4, 3, 3, 1, 0)]

AUTHORS:

- Bruce Westbury (2020): initial version
"""

#*****************************************************************************
#       Copyright (C) 2020 Bruce Westbury <bruce.westbury@gmail.com>,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.path_tableaux.path_tableau import PathTableau, PathTableaux
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.skew_tableau import SkewTableau, SkewTableaux
from sage.combinat.tableau import Tableau
from sage.combinat.gelfand_tsetlin_patterns import GelfandTsetlinPattern
from sage.combinat.partition import _Partitions
from sage.rings.semirings.non_negative_integer_semiring import NN

###############################################################################


class SemistandardPathTableau(PathTableau):
    r"""
    An instance is a sequence of lists. Usually the entries will be nonnegative integers
    in which case this is the chain of partitions of a (skew) semistandard tableau.
    In general the entries are elements of an ordered abelian group; each list is weakly
    decreasing and successive lists are interleaved.

    INPUT:

    Can be any of the following

    * a sequence of partitions
    * a sequence of lists/tuples
    * a semistandard tableau
    * a semistandard skew tableau
    * a Gelfand-Tsetlin pattern

    EXAMPLES::

        sage: path_tableaux.SemistandardPathTableau([[], [2], [2,1]])
        [(), (2,), (2, 1)]

        sage: gt = GelfandTsetlinPattern([[2,1], [2]])
        sage: path_tableaux.SemistandardPathTableau(gt)
        [(), (2,), (2, 1)]

        sage: st = SemistandardTableau([[1,1], [2]])
        sage: path_tableaux.SemistandardPathTableau(st)
        [(), (2,), (2, 1)]

        sage: st = SkewTableau([[1,1], [2]])
        sage: path_tableaux.SemistandardPathTableau(st)
        [(), (2,), (2, 1)]

        sage: st = SkewTableau([[None,1,1], [2]])
        sage: path_tableaux.SemistandardPathTableau(st)
        [(1,), (3, 0), (3, 1, 0)]

        sage: path_tableaux.SemistandardPathTableau([[], [5/2], [7/2,2]])
        [(), (5/2,), (7/2, 2)]

        sage: path_tableaux.SemistandardPathTableau([[], [2.5], [3.5,2]])
        [(), (2.50000000000000,), (3.50000000000000, 2)]
    """

    @staticmethod
    def __classcall_private__(cls, st, check=True):
        r"""
        Ensure that a tableau is only ever constructed as an
        ``element_class`` call of an appropriate parent.

        EXAMPLES::

            sage: t = path_tableaux.SemistandardPathTableau([[], [2]])
            sage: t.parent()
            <sage.combinat.path_tableaux.semistandard.SemistandardPathTableaux_with_category object at ...>
        """
        return SemistandardPathTableaux()(st, check=check)

    def __init__(self, parent, st, check=True):
        r"""
        Initialize a semistandard tableau.

        TESTS::

            sage: path_tableaux.SemistandardPathTableau([(), 3, (3, 2)])
            Traceback (most recent call last):
            ...
            ValueError: [(), 3, (3, 2)] is not a sequence of lists
        """
        w = None

        if isinstance(st, SemistandardPathTableau):
            w = list(st)

        elif isinstance(st, GelfandTsetlinPattern):
            w = list(st)
            w.reverse()
            w = [(),*w]

        elif isinstance(st, (Tableau,SkewTableau)):
            w = st.to_chain()

        elif isinstance(st, (list,tuple)):
            if any(not isinstance(a,(list,tuple)) for a in st):
                raise ValueError(f"{st} is not a sequence of lists")
            w = st

        else:
            raise ValueError(f"invalid input {st} is of type {type(st)}")

        # Pad with zeroes, if necessary
        m = max(len(a)-i for i,a in enumerate(w))
        w = [list(a)+[0]*(m+i-len(a)) for i,a in enumerate(w)]
        # Convert to immutable
        w = tuple([tuple(a) for a in w])

        PathTableau.__init__(self, parent, w, check=check)

    def check(self):
        """
        Check that ``self`` is a valid path.

        EXAMPLES::

            sage: path_tableaux.SemistandardPathTableau([[], [3], [2,2]])       # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: [(), (3,), (2, 2)] does not satisfy
            the required inequalities in row 1

            sage: path_tableaux.SemistandardPathTableau([[], [3/2], [2,5/2]])   # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: [(), (3/2,), (2, 5/2)] does not satisfy
            the required inequalities in row 1


        TESTS::

            sage: path_tableaux.SemistandardPathTableau([[], [2], [1,2]])
            Traceback (most recent call last):
            ...
            ValueError: [(), (2,), (1, 2)] does not satisfy the required inequalities in row 1

            sage: path_tableaux.SemistandardPathTableau([[], [2], [1,2]], check=False)
            [(), (2,), (1, 2)]
        """
        for i in range(1,len(self)-1):
            if not all(r >= s for r,s in zip(self[i+1],self[i])):
                raise ValueError(f"{self} does not satisfy the required inequalities in row {i}")
            if not all(r >= s for r,s in zip(self[i],self[i+1][1:])):
                raise ValueError(f"{self} does not satisfy the required inequalities in row {i}")

    def size(self):
        r"""
        Return the size or length of ``self``.

        EXAMPLES::

            sage: path_tableaux.SemistandardPathTableau([[], [3], [3,2], [3,3,1], [3,3,2,1]]).size()
            5
        """
        return len(self)

    def is_skew(self):
        """
        Return ``True`` if ``self`` is skew.

        EXAMPLES::

            sage: path_tableaux.SemistandardPathTableau([[], [2]]).is_skew()
            False
            sage: path_tableaux.SemistandardPathTableau([[2,1]]).is_skew()
            True
        """
        return bool(self[0])

    def is_integral(self) -> bool:
        """
        Return ``True`` if all entries are nonnegative integers.

        EXAMPLES::

            sage: path_tableaux.SemistandardPathTableau([[], [3], [3,2]]).is_integral()
            True
            sage: path_tableaux.SemistandardPathTableau([[], [5/2], [7/2,2]]).is_integral()
            False
            sage: path_tableaux.SemistandardPathTableau([[], [3], [3,-2]]).is_integral()
            False
        """
        return all(i in NN for a in self for i in a)

    def local_rule(self, i):
        r"""
        This is the Bender-Knuth involution.

        This is implemented by toggling the entries of the `i`-th list.
        The allowed range for ``i`` is ``0 < i < len(self)-1`` so any row except
        the first and last can be changed.

        EXAMPLES::

            sage: pt = path_tableaux.SemistandardPathTableau([[], [3], [3,2],
            ....:                                             [3,3,1], [3,3,2,1]])
            sage: pt.local_rule(1)
            [(), (2,), (3, 2), (3, 3, 1), (3, 3, 2, 1)]
            sage: pt.local_rule(2)
            [(), (3,), (3, 2), (3, 3, 1), (3, 3, 2, 1)]
            sage: pt.local_rule(3)
            [(), (3,), (3, 2), (3, 2, 2), (3, 3, 2, 1)]

        TESTS::

            sage: pt = path_tableaux.SemistandardPathTableau([[], [3], [3,2],
            ....:                                             [3,3,1], [3,3,2,1]])
            sage: pt.local_rule(0)
            Traceback (most recent call last):
            ...
            ValueError: 0 is not defined on [(), (3,), (3, 2), (3, 3, 1), (3, 3, 2, 1)]
            sage: pt.local_rule(4)
            Traceback (most recent call last):
            ...
            ValueError: 4 is not defined on [(), (3,), (3, 2), (3, 3, 1), (3, 3, 2, 1)]
        """
        def toggle(i, j):
            """
            Return the toggle of entry 'self[i][j]'.
            """

            if j == 0:
                left = self[i+1][0]
            else:
                left = min(self[i+1][j], self[i-1][j-1])
            if j == len(self[i])-1:
                right = self[i+1][j+1]
            else:
                right = max(self[i+1][j+1], self[i-1][j])

            return left + right - self[i][j]

        if not 0 < i < self.size()-1:
            raise ValueError(f"{i} is not defined on {self}")

        with self.clone() as result:
            result[i] = tuple([toggle(i,k) for k in range(len(self[i]))])

        return result

    def rectify(self, inner=None, verbose=False):
        """
        Rectify ``self``.

        This gives the usual rectification of a skew standard tableau and gives a
        generalisation to skew semistandard tableaux. The usual construction uses
        jeu de taquin but here we use the Bender-Knuth involutions.

        EXAMPLES::

            sage: st = SkewTableau([[None, None, None, 4], [None, None, 1, 6],
            ....:                   [None, None, 5], [2, 3]])
            sage: path_tableaux.SemistandardPathTableau(st).rectify()
            [(), (1,), (1, 1), (2, 1, 0), (3, 1, 0, 0), (3, 2, 0, 0, 0), (4, 2, 0, 0, 0, 0)]
            sage: path_tableaux.SemistandardPathTableau(st).rectify(verbose=True)
            [[(3, 2, 2), (3, 3, 2, 0), (3, 3, 2, 1, 0), (3, 3, 2, 2, 0, 0),
              (4, 3, 2, 2, 0, 0, 0), (4, 3, 3, 2, 0, 0, 0, 0), (4, 4, 3, 2, 0, 0, 0, 0, 0)],
             [(3, 2), (3, 3, 0), (3, 3, 1, 0), (3, 3, 2, 0, 0), (4, 3, 2, 0, 0, 0),
              (4, 3, 3, 0, 0, 0, 0), (4, 4, 3, 0, 0, 0, 0, 0)],
             [(3,), (3, 1), (3, 1, 1), (3, 2, 1, 0), (4, 2, 1, 0, 0), (4, 3, 1, 0, 0, 0),
              (4, 4, 1, 0, 0, 0, 0)],
             [(), (1,), (1, 1), (2, 1, 0), (3, 1, 0, 0), (3, 2, 0, 0, 0), (4, 2, 0, 0, 0, 0)]]

        TESTS::

            sage: S = SemistandardSkewTableaux([[5,3,3], [3,1]], [3,2,2])
            sage: LHS = [path_tableaux.SemistandardPathTableau(st.rectify()) for st in S]
            sage: RHS = [path_tableaux.SemistandardPathTableau(st).rectify() for st in S]
            sage: LHS == RHS
            True

            sage: st = SkewTableau([[None, None, None, 4], [None, None, 1, 6],
            ....:                   [None, None, 5], [2, 3]])
            sage: pt = path_tableaux.SemistandardPathTableau(st)
            sage: SP = [path_tableaux.SemistandardPathTableau(it)
            ....:       for it in StandardTableaux([3,2,2])]
            sage: len(set(pt.rectify(inner=ip) for ip in SP))
            1
        """
        if not self.is_skew():
            return self

        n = len(self)
        pp = self[0]
        P = self.parent()

        if inner is None:
            initial = [pp[:r] for r in range(len(pp))]
        elif _Partitions(inner[-1]) == _Partitions(pp):
            initial = list(inner)[:-1]
        else:
            raise ValueError(f"the final shape{inner[-1]} must agree with the initial shape {pp}")

        r = len(initial)
        path = P.element_class(P, initial + list(self))
        if verbose:
            rect = [self]

        for i in range(r):
            for j in range(n-1):
                path = path.local_rule(r+j-i)
            if verbose:
                rect.append(P.element_class(P, list(path)[r-i-1:r+n-i-1]))

        if verbose:
            return rect
        else:
            return P.element_class(P, list(path)[:n])

    @combinatorial_map(name='to semistandard tableau')
    def to_tableau(self):
        r"""
        Convert ``self`` to a :class:`SemistandardTableau`.

        The :class:`SemistandardSkewTableau` is not implemented so this returns a :class:`SkewTableau`

        EXAMPLES::

            sage: pt = path_tableaux.SemistandardPathTableau([[], [3], [3,2], [3,3,1],
            ....:                                             [3,3,2,1], [4,3,3,1,0]])
            sage: pt.to_tableau()
            [[1, 1, 1, 5], [2, 2, 3], [3, 4, 5], [4]]

        TESTS::

            sage: SST = SemistandardTableaux(shape=[5,5,3], eval=[2,2,3,4,2])
            sage: all(st == path_tableaux.SemistandardPathTableau(st).to_tableau()      # needs sage.modules
            ....:     for st in SST)
            True
        """
        from sage.combinat.tableau import from_chain

        if not self.is_integral():
            raise ValueError(f"{self} must have all entries nonnegative integers")

        lt = [[i for i in a if i > 0] for a in self]
        if self.is_skew():
            return SkewTableaux().from_chain(lt)
        else:
            return from_chain(lt)

    @combinatorial_map(name='to Gelfand-Tsetlin pattern')
    def to_pattern(self):
        r"""
        Convert ``self`` to a Gelfand-Tsetlin pattern.

        EXAMPLES::

            sage: pt = path_tableaux.SemistandardPathTableau([[], [3], [3,2], [3,3,1],
            ....:                                             [3,3,2,1], [4,3,3,1]])
            sage: pt.to_pattern()
            [[4, 3, 3, 1, 0], [3, 3, 2, 1], [3, 3, 1], [3, 2], [3]]

        TESTS::

            sage: pt = path_tableaux.SemistandardPathTableau([[3,2], [3,3,1], [3,3,2,1], [4,3,3,1]])
            sage: pt.to_pattern()
            Traceback (most recent call last):
            ...
            ValueError: [(3, 2), (3, 3, 1), (3, 3, 2, 1), (4, 3, 3, 1, 0)] cannot be a skew tableau

            sage: GT = GelfandTsetlinPatterns(top_row=[5,5,3])
            sage: all(gt == path_tableaux.SemistandardPathTableau(gt).to_pattern() for gt in GT)
            True

            sage: GT = GelfandTsetlinPatterns(top_row=[5,5,3])
            sage: all(gt.to_tableau() == path_tableaux.SemistandardPathTableau(gt).to_tableau() for gt in GT)
            True
        """
        if self.is_skew():
            raise ValueError(f"{self} cannot be a skew tableau")

        lt = list(self)
        lt.reverse()
        if not lt[-1]:
            lt.pop()

        return GelfandTsetlinPattern([list(a) for a in lt])

    def _test_jdt_promotion(self, **options):
        """
        Check that promotion agrees with :meth:`Tableau.promotion_inverse`
        constructed using jeu de taquin.

        TESTS::

            sage: pt = path_tableaux.SemistandardPathTableau([(), (1,), (2,1), (4,2),
            ....:                                             (4,3,1), (4,3,3)])
            sage: pt._test_jdt_promotion()

            sage: pt = path_tableaux.SemistandardPathTableau([(), (1,), (2,1), (4,2),
            ....:                                             (4,3,1), (9/2,3,3)])
            sage: pt._test_jdt_promotion()
            Traceback (most recent call last):
            ...
            ValueError: [(), (1,), (2, 1), (4, 2, 0), (4, 3, 1, 0), (9/2, 3, 3, 0, 0)]
            must have all entries nonnegative integers
        """
        if not self.is_integral():
            raise ValueError(f"{self} must have all entries nonnegative integers")

        tester = self._tester(**options)
        LHS = self.promotion().to_tableau()
        RHS = self.to_tableau().promotion_inverse(len(self)-2)
        tester.assertEqual(LHS,RHS)


class SemistandardPathTableaux(PathTableaux):
    """
    The parent class for :class:`SemistandardPathTableau`.
    """

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: path_tableaux.SemistandardPathTableaux()._an_element_()
            [(), (2,), (2, 1)]
        """
        return SemistandardPathTableau([[], [2], [2,1]])

    Element = SemistandardPathTableau
