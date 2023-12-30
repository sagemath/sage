r"""
Tabloids

AUTHORS:

- Sacha Goldman, Travis Scrimshaw (2023): initial version
"""

from sage.combinat.partition import Partition
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.arith.misc import multinomial
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.list_clone import ClonableArray
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.typeset.ascii_art import ascii_art
from sage.rings.integer_ring import ZZ

class Tabloid(ClonableArray):
    r"""
    A tabloid.

    A tabloid is an ordered set partition of `\{1, 2, \ldots, n\}` such that
    the sizes of the sets are weakly decreasing (that is, they form an
    integer partition).
    """
    def __init__(self, parent, tableaux, check=True):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: T = Tabloids([5, 3, 1])
            sage: A = T([[1,2,3,9,7], [4,5,8], [6]])
            sage: TestSuite(A).run()
        """
        super().__init__(parent, map(frozenset, tableaux), check=check)

    def check(self):
        r"""
        Check that ``self`` is a valid tabloid.

        EXAMPLES::

            sage: T = Tabloids([3, 2, 1])
            sage: A = T([[1,2,3], [4,5], [6]])
            sage: A.check()
            sage: T([[1,2,3,4], [5], [6]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: incorrect shape
            sage: T([[1,2,3], [4,'a'], [6]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: invalid entry in the row frozenset({...})
            sage: T([[1,2,3], [4,3], [6]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: invalid entry in the row frozenset({3, 4})
        """
        la = self.parent()._shape
        if [len(row) for row in self] != la:
            raise ValueError("incorrect shape")
        X = set(range(1, sum(la)+1))
        for row in self:
            ell = len(X)
            X -= row
            if len(X) != ell - len(row):
                raise ValueError(f"invalid entry in the row {row}")
        assert not X

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: T = Tabloids([5, 3, 1])
            sage: T([[1,2,3,9,7], [4,5,8], [6]])
            [{1, 2, 3, 7, 9}, {4, 5, 8}, {6}]
        """
        ret = "["
        ret += ", ".join("{" + ", ".join(str(val) for val in sorted(row)) + "}" for row in self)
        return ret + "]"

    def _ascii_art_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: T = Tabloids([5, 3, 1])
            sage: A = T([[1,2,3,9,7], [4,5,8], [6]])
            sage: ascii_art(A)
            {1, 2, 3, 7, 9}
            {4, 5, 8}
            {6}
        """
        ret = "\n".join("{" + ", ".join(str(val) for val in sorted(row)) + "}" for row in self)
        return ascii_art(ret)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self`` as a Young diagram.

        EXAMPLES::

            sage: T = Tabloids([5, 3, 1])
            sage: A = T([[1,2,3,9,7], [4,5,8], [6]])
            sage: latex(A)
            {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{5}c}\cline{1-5}
            \lr{1}&\lr{2}&\lr{3}&\lr{7}&\lr{9}\\\cline{1-5}
            \lr{4}&\lr{5}&\lr{8}\\\cline{1-3}
            \lr{6}\\\cline{1-1}
            \end{array}$}
            }

            sage: T = Tabloids([])
            sage: latex(T([]))
            {\emptyset}
        """
        if not self:
            return "{\\emptyset}"
        from sage.combinat.output import tex_from_array
        A = list(map(sorted, self))
        ret = str(tex_from_array(A))
        return ret.replace("|", "")

    def symmetric_group_action(self, permutation):
        r"""
        Return the left action of ``permutation`` on ``self``.

        EXAMPLES::

            sage: T = Tabloids([3, 2, 1])
            sage: A = T([[1,2,6], [3,5], [4]]); ascii_art(A)
            {1, 2, 6}
            {3, 5}
            {4}
            sage: sigma = Permutation([1,6,2,3,5,4])
            sage: ascii_art(A.symmetric_group_action(sigma))
            {1, 4, 6}
            {2, 5}
            {3}
            sage: sigma = SymmetricGroup(6)([(1,2),(3,4,5)])
            sage: ascii_art(A.symmetric_group_action(sigma))
            {1, 2, 6}
            {3, 4}
            {5}
        """
        P = self.parent()
        S = SymmetricGroup(sum(P._shape))
        permutation = S(permutation)
        return P.element_class(P, [[permutation(val) for val in row] for row in self])

    def _acted_upon_(self, sigma, self_on_left=True):
        r"""
        Return the action of ``sigma`` on ``self``.

        EXAMPLES::

            sage: T = Tabloids([2, 2, 1])
            sage: A = T([[1,2], [3,5], [4]]); ascii_art(A)
            {1, 2}
            {3, 5}
            {4}
            sage: ascii_art(Permutation([1,2,3,5,4]) * A)
            {1, 2}
            {3, 4}
            {5}
            sage: sigma = SymmetricGroup(5)([(1,2),(3,4,5)])
            sage: ascii_art(sigma * A)
            {1, 2}
            {3, 4}
            {5}
        """
        P = self.parent()
        S = SymmetricGroup(sum(P._shape))
        if sigma in S:
            return self.symmetric_group_action(sigma)
        return super()._acted_upon_(sigma, self_on_left)


class Tabloids(UniqueRepresentation, Parent):
    r"""
    Tabloids of a fixed shape.
    """
    @staticmethod
    def __classcall_private__(cls, partition):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: T1 = Tabloids([2, 1])
            sage: T2 = Tabloids((2, 1))
            sage: T3 = Tabloids(Partition([2, 1]))
            sage: T1 is T2 and T2 is T3
            True
        """
        partition = Partition(partition)
        return super().__classcall__(cls, partition)

    def __init__(self, partition):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: T = Tabloids([2, 1])
            sage: TestSuite(T).run()

            sage: T = Tabloids([])
            sage: TestSuite(T).run()

            sage: T = Tabloids([1])
            sage: TestSuite(T).run()

            sage: T = Tabloids([10])
            sage: TestSuite(T).run()

            sage: T = Tabloids([1, 1, 1, 1])
            sage: TestSuite(T).run()
        """
        self._shape = partition
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: T = Tabloids([2,1])
            sage: list(T)
            [[{1, 2}, {3}], [{1, 3}, {2}], [{2, 3}, {1}]]

            sage: T = Tabloids([])
            sage: list(T)
            [[]]
            sage: T = Tabloids([1])
            sage: list(T)
            [[{1}]]
            sage: T = Tabloids([3])
            sage: list(T)
            [[{1, 2, 3}]]
            sage: T = Tabloids([1,1,1])
            sage: list(T)
            [[{1}, {2}, {3}],
             [{3}, {1}, {2}],
             [{2}, {3}, {1}],
             [{1}, {3}, {2}],
             [{3}, {2}, {1}],
             [{2}, {1}, {3}]]
        """
        n = sum(self._shape)
        # trivial cases
        if n == 0:
            yield self.element_class(self, [], check=False)
            return
        if n == 1:
            yield self.element_class(self, [[1]], check=False)
            return

        ell = len(self._shape)
        if ell == n:  # single column trivial case
            for sigma in SymmetricGroup(n):
                yield self.element_class(self, [[val] for val in sigma.tuple()], check=False)
            return
        if ell == 1:  # single row trivial case
            yield self.element_class(self, [range(1,n+1)], check=False)
            return

        pos = [-1] * n
        i = 0
        cur = [[] for _ in range(ell)]
        while True:
            if i == n:
                yield self.element_class(self, cur, check=False)
                i -= 1
                cur[pos[i]].pop()
            pos[i] += 1
            while pos[i] < ell and self._shape[pos[i]] - len(cur[pos[i]]) == 0:
                pos[i] += 1
            if pos[i] == ell: # backtrack
                pos[i] = -1
                i -= 1
                if i < 0:
                    break
                cur[pos[i]].pop()
                continue
            cur[pos[i]].append(i+1)
            i += 1

    def from_tableau(self, T):
        r"""
        Construct a tabloid from the tableau ``T``.

        EXAMPLES::

            sage: T = Tabloids([3, 1])
            sage: T.from_tableau([[1,3,4], [2]])
            [{1, 3, 4}, {2}]
        """
        return self.element_class(self, T)

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: T = Tabloids([3,2,1])
            sage: T.cardinality()
            60
        """
        return ZZ(multinomial(list(self._shape)))

    Element = Tabloid
