r"""
Path tableaux

This is an abstract base class for using local rules to construct
rectification and the action of the cactus group [Wes2017]_.

This is a construction of the Henriques-Kamnitzer construction of
the action of the cactus group on tensor powers of a crystal. This is
also a generalisation of the Fomin growth rules, which are a version of
the operations on standard tableaux which were previously constructed
using jeu de taquin.

The basic operations are rectification, evacuation and promotion.
Rectification of standard skew tableaux agrees with the rectification
by jeu de taquin as does evacuation. Promotion agrees with promotion
by jeu de taquin on rectangular tableaux but in general they are different.

REFERENCES:

- [Wes2017]_

AUTHORS:

- Bruce Westbury (2018): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Bruce Westbury <bruce.westbury@gmail.com>,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.abstract_method import abstract_method
from sage.categories.sets_cat import Sets
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
#from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_graded
from sage.structure.sage_object import SageObject
from sage.structure.list_clone import ClonableArray
from sage.misc.latex import latex
#from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
#from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets


class PathTableau(ClonableArray, metaclass=InheritComparisonClasscallMetaclass):
    r"""
    This is the abstract base class for a path tableau.
    """
    @abstract_method
    def local_rule(self, i):
        r"""
        This is the abstract local rule defined in any coboundary category.

        This has input a list of objects. This method first takes
        the list of objects of length three consisting of the `(i-1)`-st,
        `i`-th and `(i+1)`-term and applies the rule. It then replaces
        the `i`-th object  by the object returned by the rule.

        EXAMPLES::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t.local_rule(3)
            [0, 1, 2, 1, 2, 1, 0]
        """

    ################################# Book Keeping ############################

    def size(self):
        r"""
        Return the size or length of ``self``.

        EXAMPLES::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t.size()
            7
        """
        return len(self)

    def initial_shape(self):
        r"""
        Return the initial shape of ``self``.

        EXAMPLES::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t.initial_shape()
            0
        """
        return self[0]

    def final_shape(self):
        r"""
        Return the final shape of ``self``.

        EXAMPLES::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t.final_shape()
            0
        """
        return self[-1]

    ############################# Jeu de taquin ###############################

    def promotion(self):
        r"""
        Return the promotion operator applied to ``self``.

        EXAMPLES::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t.promotion()
            [0, 1, 2, 1, 0, 1, 0]
        """
        with self.clone() as result:
            for i in range(1,self.size()-1):
                result = result.local_rule(i)

        return result

    def evacuation(self):
        r"""
        Return the evacuation operator applied to ``self``.

        EXAMPLES::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t.evacuation()
            [0, 1, 2, 3, 2, 1, 0]
        """
        if self.size() < 3:
            return self

        L = list(self)
        result = []
        P = self.parent()
        for i in range(self.size()):
            L = list(P(L).promotion())
            result.append( L.pop() )
        result.reverse()
        return P(result)

    def commutor(self, other, verbose=False):
        r"""
        Return the commutor of ``self`` with ``other``.

        If ``verbose=True`` then the function will print
        the rectangle.

        EXAMPLES::

            sage: t1 = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t2 = path_tableaux.DyckPath([0,1,2,1,0])
            sage: t1.commutor(t2)
            ([0, 1, 2, 1, 0], [0, 1, 2, 3, 2, 1, 0])
            sage: t1.commutor(t2,verbose=True)
            [0, 1, 2, 1, 0]
            [1, 2, 3, 2, 1]
            [2, 3, 4, 3, 2]
            [3, 4, 5, 4, 3]
            [2, 3, 4, 3, 2]
            [1, 2, 3, 2, 1]
            [0, 1, 2, 1, 0]
            ([0, 1, 2, 1, 0], [0, 1, 2, 3, 2, 1, 0])

        TESTS::

            sage: t1 = path_tableaux.DyckPath([])
            sage: t2 = path_tableaux.DyckPath([0,1,2,1,0])
            sage: t1.commutor(t2)
            Traceback (most recent call last):
            ...
            ValueError: this requires nonempty lists
            sage: t1 = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t2 = path_tableaux.DyckPath([])
            sage: t1.commutor(t2)
            Traceback (most recent call last):
            ...
            ValueError: this requires nonempty lists
            sage: t1 = path_tableaux.DyckPath([0,1,2,3,2,1])
            sage: t2 = path_tableaux.DyckPath([0,1,2,1,0])
            sage: t1.commutor(t2)
            Traceback (most recent call last):
            ...
            ValueError: [0, 1, 2, 3, 2, 1], [0, 1, 2, 1, 0] is not a composable pair
        """
        n = self.size()
        m = len(other)
        if n == 0 or m == 0:
            raise ValueError("this requires nonempty lists")
        if n == 1 or m == 1:
            return (other, self)

        P = self.parent()

        row = list(other)
        col = list(self)
        if col[-1] != row[0]:
            raise ValueError("%s, %s is not a composable pair" % (self,other))

        path = P(col + row[1:])

        for i in range(1,n):
            if verbose:
                print(path[n-i:n+m-i])
            for j in range(m-1):
                path = path.local_rule(n+j-i)
        if verbose:
            print(path[:m])

        return (P(path[:m]), P(path[m-1:]))

    def cactus(self, i, j):
        r"""
        Return the action of the generator `s_{i,j}` of the cactus
        group on ``self``.

        INPUT:

        - ``i`` -- positive integer
        - ``j`` -- positive integer weakly greater than `i`

        EXAMPLES::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t.cactus(1,5)
            [0, 1, 0, 1, 2, 1, 0]

            sage: t.cactus(1,6)
            [0, 1, 2, 1, 0, 1, 0]

            sage: t.cactus(1,7) == t.evacuation()
            True
            sage: t.cactus(1,7).cactus(1,6) == t.promotion()
            True

        TESTS::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t.cactus(1,8)
            Traceback (most recent call last):
            ...
            ValueError: integers out of bounds
            sage: t.cactus(0,3)
            Traceback (most recent call last):
            ...
            ValueError: integers out of bounds
        """
        if not 0 < i <= j <= self.size():
            raise ValueError("integers out of bounds")

        if i == j:
            return self

        if i == 1:
            h = list(self)[:j]
            t = list(self)[j:]
            T = self.parent()(h)
            L = list(T.evacuation()) + t
            return self.parent()(L)

        return self.cactus(1,j).cactus(1,j-i+1).cactus(1,j)

    ########################### Visualisation and checking ####################

    def _test_involution_rule(self, **options):
        """
        Check that the local rule gives an involution.

        TESTS::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t._test_involution_rule()
        """
        tester = self._tester(**options)
        for i in range(self.size()-2):
            tester.assertEqual(self.local_rule(i+1).local_rule(i + 1), self)

    def _test_involution_cactus(self, **options):
        """
        Check that the cactus group generators are involutions.

        TESTS::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t._test_involution_cactus()
        """
        tester = self._tester(**options)
        for i in range(2, self.size()+1):
            tester.assertEqual(self.cactus(1,i).cactus(1,i), self)

    def _test_promotion(self, **options):
        """
        Check that promotion can be expressed in terms of the cactus generators.

        TESTS::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t._test_promotion()
        """
        tester = self._tester(**options)
        n = self.size()
        tester.assertEqual(self.cactus(1,n-1).cactus(1,n).promotion(), self)

    def _test_commutation(self, **options):
        """
        Check the commutation relations in the presentation of the cactus group.

        TESTS::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t._test_commutation()
        """
        from itertools import combinations
        tester = self._tester(**options)

        n = self.size()
        if n < 5:
            return
        for i,j,r,s in combinations(range(1,n+1), 4):
            lhs = self.cactus(i, j).cactus(r, s)
            rhs = self.cactus(r, s).cactus(i, j)
            tester.assertEqual(lhs, rhs)

    def _test_coboundary(self, **options):
        """
        Check the coboundary relations in the presentation of the cactus group.

        TESTS::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t._test_coboundary()
        """
        from itertools import combinations
        tester = self._tester(**options)

        n = self.size()
        if n < 4:
            return
        for i,j,r,s in combinations(range(1,n+3), 4):
            lhs = self.cactus(i, s-2).cactus(j-1, r-1)
            rhs = self.cactus(i+s-r-1, i+s-j-1).cactus(i, s-2)
            tester.assertEqual(lhs, rhs)

    def orbit(self):
        r"""
        Return the orbit of ``self`` under the action of the cactus group.

        EXAMPLES::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: t.orbit()
            {[0, 1, 0, 1, 0, 1, 0],
             [0, 1, 0, 1, 2, 1, 0],
             [0, 1, 2, 1, 0, 1, 0],
             [0, 1, 2, 1, 2, 1, 0],
             [0, 1, 2, 3, 2, 1, 0]}
        """
        orb = set()
        rec = set([self])
        while rec:
            new = set()
            for a in rec:
                for i in range(2, self.size()):
                    b = a.cactus(1, i)
                    if (b not in orb) and (b not in rec):
                        new.add(b)
            orb = orb.union(rec)
            rec = new.copy()

        return orb

    def dual_equivalence_graph(self):
        r"""
        Return the graph with vertices the orbit of ``self``
        and edges given by the action of the cactus group generators.

        In most implementations the generators `s_{i,i+1}` will act
        as the identity operators. The usual dual equivalence graphs
        are given by replacing the label `i,i+2` by `i` and removing
        edges with other labels.

        EXAMPLES::

            sage: s = path_tableaux.DyckPath([0,1,2,3,2,3,2,1,0])
            sage: s.dual_equivalence_graph().adjacency_matrix()                         # needs sage.graphs sage.modules
            [0 1 1 1 0 1 0 1 1 0 0 0 0 0]
            [1 0 1 1 1 1 1 0 1 0 0 1 1 0]
            [1 1 0 1 1 1 0 1 0 1 1 1 0 0]
            [1 1 1 0 1 0 1 1 1 1 0 1 1 0]
            [0 1 1 1 0 0 1 0 0 1 1 0 1 1]
            [1 1 1 0 0 0 1 1 1 1 1 0 1 0]
            [0 1 0 1 1 1 0 1 0 1 1 1 0 1]
            [1 0 1 1 0 1 1 0 1 1 1 1 1 0]
            [1 1 0 1 0 1 0 1 0 1 0 1 1 0]
            [0 0 1 1 1 1 1 1 1 0 0 1 1 1]
            [0 0 1 0 1 1 1 1 0 0 0 1 1 1]
            [0 1 1 1 0 0 1 1 1 1 1 0 1 1]
            [0 1 0 1 1 1 0 1 1 1 1 1 0 1]
            [0 0 0 0 1 0 1 0 0 1 1 1 1 0]
            sage: s = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: s.dual_equivalence_graph().edges(sort=True)                           # needs sage.graphs
            [([0, 1, 0, 1, 0, 1, 0], [0, 1, 0, 1, 2, 1, 0], '4,7'),
             ([0, 1, 0, 1, 0, 1, 0], [0, 1, 2, 1, 0, 1, 0], '2,5'),
             ([0, 1, 0, 1, 0, 1, 0], [0, 1, 2, 1, 2, 1, 0], '2,7'),
             ([0, 1, 0, 1, 2, 1, 0], [0, 1, 2, 1, 0, 1, 0], '2,6'),
             ([0, 1, 0, 1, 2, 1, 0], [0, 1, 2, 1, 2, 1, 0], '1,4'),
             ([0, 1, 0, 1, 2, 1, 0], [0, 1, 2, 3, 2, 1, 0], '2,7'),
             ([0, 1, 2, 1, 0, 1, 0], [0, 1, 2, 1, 2, 1, 0], '4,7'),
             ([0, 1, 2, 1, 0, 1, 0], [0, 1, 2, 3, 2, 1, 0], '3,7'),
             ([0, 1, 2, 1, 2, 1, 0], [0, 1, 2, 3, 2, 1, 0], '3,6')]
        """
        from sage.graphs.graph import Graph
        from itertools import combinations

        G = Graph()
        orb = self.orbit()

        for a in orb:
            for i,j in combinations(range(1,self.size()+1),2):
                b = a.cactus(i,j)
                if a != b:
                    G.add_edge(a,b,"%d,%d" % (i,j))
        return G


class PathTableaux(UniqueRepresentation,Parent):
    """
    The abstract parent class for PathTableau.
    """

    def __init__(self):
        """
        Initialize ``self``.

        TESTS::

            sage: t = path_tableaux.DyckPaths()
            sage: TestSuite(t).run()

            sage: f = path_tableaux.FriezePatterns(QQ)
            sage: TestSuite(f).run()
        """
        Parent.__init__(self, category=Sets())

    def _element_constructor_(self, *args, **kwds):
        r"""
        Construct an object as an element of ``self``, if possible.

        TESTS::

            sage: path_tableaux.DyckPath([0,1,2,1,0])  # indirect doctest
            [0, 1, 2, 1, 0]
        """
        return self.element_class(self, *args, **kwds)


class CylindricalDiagram(SageObject):
    r"""
    Cylindrical growth diagrams.

    EXAMPLES::

        sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
        sage: path_tableaux.CylindricalDiagram(t)
        [0, 1, 2, 3, 2, 1, 0]
        [ , 0, 1, 2, 1, 0, 1, 0]
        [ ,  , 0, 1, 0, 1, 2, 1, 0]
        [ ,  ,  , 0, 1, 2, 3, 2, 1, 0]
        [ ,  ,  ,  , 0, 1, 2, 1, 0, 1, 0]
        [ ,  ,  ,  ,  , 0, 1, 0, 1, 2, 1, 0]
        [ ,  ,  ,  ,  ,  , 0, 1, 2, 3, 2, 1, 0]
    """

    def __init__(self, T):
        """
        Initialize ``self`` from the
        :class:`~sage.combinat.path_tableaux.path_tableau.PathTableau`
        object ``T``.

        TESTS::

            sage: T = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: D = path_tableaux.CylindricalDiagram(T)
            sage: TestSuite(D).run()

            sage: path_tableaux.CylindricalDiagram(2)
            Traceback (most recent call last):
            ...
            ValueError: 2 must be a path tableau
        """
        if not isinstance(T, PathTableau):
            raise ValueError('{0} must be a path tableau'.format(str(T)))
        n = len(T)
        result = [[None]*(2*n-1)] * n
        for i in range(n):
            result[i] = [""]*i + list(T)
            T = T.promotion()

        self.path_tableau = T
        self.diagram = result

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: cd = path_tableaux.CylindricalDiagram(path_tableaux.DyckPath([0,1,2,1,2,1,0]))
            sage: repr(cd) == cd._repr_() # indirect test
            True

            sage: cd = path_tableaux.CylindricalDiagram(path_tableaux.FriezePattern([1,3,4,5,1]))
            sage: repr(cd) == cd._repr_() # indirect test
            True

            sage: print(path_tableaux.DyckPath([0,1,2,1,2,1,0])) # indirect test
            [0, 1, 2, 1, 2, 1, 0]

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: path_tableaux.CylindricalDiagram(t)
            [0, 1, 2, 3, 2, 1, 0]
            [ , 0, 1, 2, 1, 0, 1, 0]
            [ ,  , 0, 1, 0, 1, 2, 1, 0]
            [ ,  ,  , 0, 1, 2, 3, 2, 1, 0]
            [ ,  ,  ,  , 0, 1, 2, 1, 0, 1, 0]
            [ ,  ,  ,  ,  , 0, 1, 0, 1, 2, 1, 0]
            [ ,  ,  ,  ,  ,  , 0, 1, 2, 3, 2, 1, 0]
        """
        data = [[str(x) for x in row] for row in self.diagram]
        if not data[0]:
            data[0] = ['']  # Put sometime there
        max_width = max(max(len(x) for x in row) for row in data if row)
        return '\n'.join('[' + ', '.join(' '*(max_width-len(x)) + x for x in row)
                         + ']' for row in data)

    def __eq__(self, other):
        """
        Check equality.

        EXAMPLES::

            sage: t1 = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: T1 = path_tableaux.CylindricalDiagram(t1)
            sage: t2 = path_tableaux.DyckPath([0,1,2,1,2,1,0])
            sage: T2 = path_tableaux.CylindricalDiagram(t2)
            sage: T1 == T2
            False
            sage: T1 == path_tableaux.CylindricalDiagram(t1)
            True
        """
        return isinstance(other, CylindricalDiagram) and self.diagram == other.diagram

    def __ne__(self, other):
        """
        Check inequality.

        EXAMPLES::

            sage: t1 = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: T1 = path_tableaux.CylindricalDiagram(t1)
            sage: t2 = path_tableaux.DyckPath([0,1,2,1,2,1,0])
            sage: T2 = path_tableaux.CylindricalDiagram(t2)
            sage: T1 != T2
            True
            sage: T1 != path_tableaux.CylindricalDiagram(t1)
            False
        """
        return not (self == other)

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: latex(path_tableaux.CylindricalDiagram(t))
            \begin{array}{ccccccccccccc}
            0 & 1 & 2 & 3 & 2 & 1 & 0\\
             & 0 & 1 & 2 & 1 & 0 & 1 & 0\\
             &  & 0 & 1 & 0 & 1 & 2 & 1 & 0\\
             &  &  & 0 & 1 & 2 & 3 & 2 & 1 & 0\\
             &  &  &  & 0 & 1 & 2 & 1 & 0 & 1 & 0\\
             &  &  &  &  & 0 & 1 & 0 & 1 & 2 & 1 & 0\\
             &  &  &  &  &  & 0 & 1 & 2 & 3 & 2 & 1 & 0
             \end{array}

            sage: t = path_tableaux.FriezePattern([1,3,4,5,1])
            sage: latex(path_tableaux.CylindricalDiagram(t))
            \begin{array}{ccccccccccccc}
            0 & 1 & 3 & 4 & 5 & 1 & 0\\
             & 0 & 1 & \frac{5}{3} & \frac{7}{3} & \frac{2}{3} & 1 & 0\\
             &  & 0 & 1 & 2 & 1 & 3 & 1 & 0\\
             &  &  & 0 & 1 & 1 & 4 & \frac{5}{3} & 1 & 0\\
             &  &  &  & 0 & 1 & 5 & \frac{7}{3} & 2 & 1 & 0\\
             &  &  &  &  & 0 & 1 & \frac{2}{3} & 1 & 1 & 1 & 0\\
             &  &  &  &  &  & 0 & 1 & 3 & 4 & 5 & 1 & 0
             \end{array}
        """
        D = self.diagram
        m = len(D[-1])
        result = "\\begin{array}{"+"c"*m + "}\n"
        result += "\\\\ \n".join( " & ".join(latex(a) for a in x) for x in D )
        result += "\n \\end{array}\n"
        return result

    def __len__(self):
        r"""
        Return the length of ``self``.

        TESTS::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: len(path_tableaux.CylindricalDiagram(t))
            7
        """
        return len(self.diagram)

    def _ascii_art_(self):
        r"""
        Return an ascii art representation of ``self``.

        TESTS::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: ascii_art(path_tableaux.CylindricalDiagram(t))
             0 1 2 3 2 1 0
               0 1 2 1 0 1 0
                 0 1 0 1 2 1 0
                   0 1 2 3 2 1 0
                     0 1 2 1 0 1 0
                       0 1 0 1 2 1 0
                         0 1 2 3 2 1 0

            sage: t = path_tableaux.FriezePattern([1,3,4,5,1])
            sage: ascii_art(path_tableaux.CylindricalDiagram(t))
               0   1   3   4   5   1   0
                   0   1 5/3 7/3 2/3   1   0
                       0   1   2   1   3   1   0
                           0   1   1   4 5/3   1   0
                               0   1   5 7/3   2   1   0
                                   0   1 2/3   1   1   1   0
                                       0   1   3   4   5   1   0
        """
        from sage.typeset.ascii_art import ascii_art
        from sage.misc.misc_c import prod
        data = [[ascii_art(x) for x in row] for row in self.diagram]
        if not data[0]:
            data[0] = [ascii_art('')]  # Put sometime there
        max_width = max(max(len(x) for x in row) for row in data if row)
        return prod((sum((ascii_art(' '*(max_width-len(x)+1)) + x for x in row), ascii_art(''))
                    for row in data), ascii_art(''))

    def _unicode_art_(self):
        r"""
        Return a unicode art representation of ``self``.

        TESTS::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: unicode_art(path_tableaux.CylindricalDiagram(t))
             0 1 2 3 2 1 0
               0 1 2 1 0 1 0
                 0 1 0 1 2 1 0
                   0 1 2 3 2 1 0
                     0 1 2 1 0 1 0
                       0 1 0 1 2 1 0
                         0 1 2 3 2 1 0

            sage: t = path_tableaux.FriezePattern([1,3,4,5,1])
            sage: unicode_art(path_tableaux.CylindricalDiagram(t))
               0   1   3   4   5   1   0
                   0   1 5/3 7/3 2/3   1   0
                       0   1   2   1   3   1   0
                           0   1   1   4 5/3   1   0
                               0   1   5 7/3   2   1   0
                                   0   1 2/3   1   1   1   0
                                       0   1   3   4   5   1   0
        """
        from sage.typeset.unicode_art import unicode_art
        from sage.misc.misc_c import prod
        data = [[unicode_art(x) for x in row] for row in self.diagram]
        if not data[0]:
            data[0] = [unicode_art('')]  # Put sometime there
        max_width = max(max(len(x) for x in row) for row in data if row)
        return prod((sum((unicode_art(' '*(max_width-len(x)+1)) + x for x in row), unicode_art(''))
                    for row in data), unicode_art(''))

    def pp(self):
        r"""
        A pretty print utility method.

        EXAMPLES::

            sage: t = path_tableaux.DyckPath([0,1,2,3,2,1,0])
            sage: path_tableaux.CylindricalDiagram(t).pp()
            0 1 2 3 2 1 0
              0 1 2 1 0 1 0
                0 1 0 1 2 1 0
                  0 1 2 3 2 1 0
                    0 1 2 1 0 1 0
                      0 1 0 1 2 1 0
                        0 1 2 3 2 1 0

            sage: t = path_tableaux.FriezePattern([1,3,4,5,1])
            sage: path_tableaux.CylindricalDiagram(t).pp()
              0   1   3   4   5   1   0
                  0   1 5/3 7/3 2/3   1   0
                      0   1   2   1   3   1   0
                          0   1   1   4 5/3   1   0
                              0   1   5 7/3   2   1   0
                                  0   1 2/3   1   1   1   0
                                      0   1   3   4   5   1   0
        """
        data = [[str(x) for x in row] for row in self.diagram]
        if not data[0]:
            data[0] = ['']  # Put sometime there
        max_width = max(max(len(x) for x in row) for row in data if row)
        print('\n'.join(' '.join(' '*(max_width-len(x)) + x for x in row)
                        for row in data))
