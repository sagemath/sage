r"""
Perfect matchings

A perfect matching of a set `S` is a partition into 2-element sets. If `S` is
the set `\{1,...,n\}`, it is equivalent to fixpoint-free involutions. These
simple combinatorial objects appear in different domains such as combinatoric
of orthogonal polynomials and of the hyperoctaedral groups (see [MV]_, [McD]_
and also [CM]_):

AUTHOR:

    - Valentin Feray, 2010 : initial version

EXAMPLES:

    Create a perfect matching::

        sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')]);m
        [('a', 'e'), ('b', 'c'), ('d', 'f')]

    Count its crossings, if the ground set is totally ordered::

        sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
        [(1, 3), (2, 8), (4, 7), (5, 6)]
        sage: n.number_of_crossings()
        1

    List the perfect matchings of a given ground set::

        sage: PerfectMatchings(4).list()
        [[(1, 2), (3, 4)], [(1, 3), (2, 4)], [(1, 4), (2, 3)]]

REFERENCES:

    .. [MV] combinatorics of orthogonal polynomials (A. de Medicis et
       X.Viennot, Moments des q-polynomes de Laguerre et la bijection de
       Foata-Zeilberger, Adv. Appl. Math., 15 (1994), 262-304)

    .. [McD] combinatorics of hyperoctahedral group, double coset algebra and
       zonal polynomials (I. G. Macdonald, Symmetric functions and Hall
       polynomials, Oxford University Press, second edition, 1995, chapter
       VII).

    .. [CM] Benoit Collins, Sho Matsumoto, On some properties of
       orthogonal Weingarten functions, :arxiv:`0903.5143`.
"""
#*****************************************************************************
#       Copyright (C) 2010 Valentin Feray <feray@labri.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
# python3
from __future__ import division, print_function


from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.element_wrapper import ElementWrapper
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.misc.flatten import flatten
from sage.combinat.permutation import Permutation
from sage.sets.set import Set
from sage.combinat.partition import Partition
from sage.misc.misc_c import prod
from sage.matrix.constructor import Matrix
from sage.combinat.combinatorial_map import combinatorial_map


class PerfectMatching(ElementWrapper):
    r"""
    Class of perfect matching.

    An instance of the class can be created from a list of pairs or from a
    fixed point-free involution as follows::

        sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')]);m
        [('a', 'e'), ('b', 'c'), ('d', 'f')]
        sage: n = PerfectMatching([3,8,1,7,6,5,4,2]);n
        [(1, 3), (2, 8), (4, 7), (5, 6)]
        sage: isinstance(m,PerfectMatching)
        True

    The parent, which is the set of perfect matchings of the ground set, is
    automatically created::

        sage: n.parent()
        Set of perfect matchings of {1, 2, 3, 4, 5, 6, 7, 8}

    If the ground set is ordered, one can, for example, ask if the matching is
    non crossing::

        sage: PerfectMatching([(1, 4), (2, 3), (5, 6)]).is_non_crossing()
        True

    TESTS::

        sage: m = PerfectMatching([]); m
        []
        sage: m.parent()
        Set of perfect matchings of {}
    """
    #the data structure of the element is a list (accessible via x.value)
    wrapped_class = list
    __lt__ = ElementWrapper._lt_by_value
    #During the creation of the instance of the class, the function
    #__classcall_private__ will be called instead of __init__ directly.
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, p):
        r"""
        This function tries to recognize the input (it can be either a list or
        a tuple of pairs, or a fix-point free involution given as a list or as
        a permutation), constructs the parent (enumerated set of
        PerfectMatchings of the ground set) and calls the __init__ function to
        construct our object.

        EXAMPLES::

            sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')]);m
            [('a', 'e'), ('b', 'c'), ('d', 'f')]
            sage: isinstance(m,PerfectMatching)
            True
            sage: n = PerfectMatching([3, 8, 1, 7, 6, 5, 4, 2]);n
            [(1, 3), (2, 8), (4, 7), (5, 6)]
            sage: n.parent()
            Set of perfect matchings of {1, 2, 3, 4, 5, 6, 7, 8}
            sage: PerfectMatching([(1, 4), (2, 3), (5, 6)]).is_non_crossing()
            True

        The function checks that the given list or permutation is a valid perfect
        matching (i.e. a list of pairs with pairwise disjoint elements  or a
        fixpoint-free involution) and raises a ValueError otherwise::

            sage: PerfectMatching([(1, 2, 3), (4, 5)])
            Traceback (most recent call last):
            ...
            ValueError: [(1, 2, 3), (4, 5)] is not a valid perfect matching: all elements of the list must be pairs

        If you know your datas are in a good format, use directly
        ``PerfectMatchings(objects)(data)``.

        TESTS::

             sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')])
             sage: TestSuite(m).run()
             sage: m = PerfectMatching([])
             sage: TestSuite(m).run()
             sage: PerfectMatching(6)
             Traceback (most recent call last):
             ...
             ValueError: cannot convert p (= 6) to a PerfectMatching
             sage: PerfectMatching([(1,2,3)])
             Traceback (most recent call last):
             ...
             ValueError: [(1, 2, 3)] is not a valid perfect matching:
             all elements of the list must be pairs
             sage: PerfectMatching([(1,1)])
             Traceback (most recent call last):
             ...
             ValueError: [(1, 1)] is not a valid perfect matching:
             there are some repetitions
             sage: PerfectMatching(Permutation([4,2,1,3]))
             Traceback (most recent call last):
             ...
             ValueError: The permutation p (= [4, 2, 1, 3]) is not a fixed point free involution
        """
        # we have to extract from the argument p the set of objects of the
        # matching and the list of pairs.
        # First case: p is a list (resp tuple) of lists (resp tuple).
        if (isinstance(p, list) or isinstance(p, tuple)) and (
                all([isinstance(x, list) or isinstance(x, tuple) for x in p])):
            objects = Set(flatten(p))
            data = [tuple(_) for _ in p]
            #check if the data are correct
            if not all([len(t) == 2 for t in data]):
                raise ValueError("%s is not a valid perfect matching:\n"
                                 "all elements of the list must be pairs" % p)
            if len(objects) < 2*len(data):
                raise ValueError("%s is not a valid perfect matching:\n"
                                 "there are some repetitions" % p)
        # Second case: p is a permutation or a list of integers, we have to
        # check if it is a fix-point-free involution.
        elif ((isinstance(p, list) and
               all(((isinstance(x, Integer) or isinstance(x, int)) for x in p)))
              or isinstance(p, Permutation)):
            p = Permutation(p)
            n = len(p)
            if not(p.cycle_type() == [2 for i in range(n//2)]):
                raise ValueError("The permutation p (= %s) is not a "
                                 "fixed point free involution" % p)
            objects = Set(range(1, n+1))
            data = p.to_cycles()
        # Third case: p is already a perfect matching, we return p directly
        elif isinstance(p, PerfectMatching):
            return p
        else:
            raise ValueError("cannot convert p (= %s) to a PerfectMatching" % p)
        # Finally, we create the parent and the element using the element
        # class of the parent. Note: as this function is private, when we
        # create an object via parent.element_class(...), __init__ is directly
        # executed and we do not have an infinite loop.
        return PerfectMatchings(objects)(data)

    def __iter__(self):
        r"""
        Iterate over the edges of the matching ``self``.

        The edges are yielded as 2-tuples. Neither the elements of these
        tuples nor the tuples are necessarily sorted in any predictable
        way.

        EXAMPLES::

            sage: list(PerfectMatching([('a','e'),('b','c'),('d','f')]))
            [('a', 'e'), ('b', 'c'), ('d', 'f')]
            sage: list(PerfectMatchings(2)[0])
            [(1, 2)]
            sage: list(PerfectMatching([3,8,1,7,6,5,4,2]))
            [(1, 3), (2, 8), (4, 7), (5, 6)]
        """
        return iter(self.value)

    def _repr_(self):
        r"""
        Return a string representation of the matching ``self``.

        EXAMPLES::

            sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')]);m
            [('a', 'e'), ('b', 'c'), ('d', 'f')]
            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]);n
            [(1, 3), (2, 8), (4, 7), (5, 6)]
        """
        return '%s' % self.value

    def _latex_(self):
        r"""
        A latex representation of ``self`` using the tikzpicture package.

        EXAMPLES::

            sage: P = PerfectMatching([(1,3),(2,5),(4,6)])
            sage: latex(P)  # optional - dot2tex; random
            \begin{tikzpicture}
            ...
            \end{tikzpicture}

        TESTS:

        Above we added ``random`` since warnings might be displayed
        once. The second time, there should be no warnings::

            sage: print(P._latex_())  # optional - dot2tex
            \begin{tikzpicture}
            ...
            \end{tikzpicture}
        """
        G = self.to_graph()
        G.set_pos(G.layout_circular())
        G.set_latex_options(
            vertex_size=0.4,
            edge_thickness=0.04,
        )
        return G._latex_()

    def __hash__(self):
        r"""
        Returns the hash of ``self`` using the tupled value.

        EXAMPLES::

            sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')])
            sage: m.__hash__() #random
            1053935254331348997
            sage: hash(m) #indirect doctest #random

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2])
            sage: hash(n) #indirect doctest #random
            8097274995140737937
        """
        return hash(tuple(self.value))

    def __eq__(self, other):
        r"""
        Compares two perfect matchings

        EXAMPLES::

            sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')])
            sage: n = PerfectMatching([('c','b'),('d','f'),('e','a')])
            sage: n == m
            True
            sage: n == PerfectMatching([('a','b'),('d','f'),('e','c')])
            False

        """
        try:
            if other.parent() != self.parent():
                return False
        except AttributeError:
            return False
        return Set([Set(_) for _ in self.value]) == Set([Set(_) for _ in other.value])

    def size(self):
        r"""

        Returns the size of the perfect matching ``self``, i.e. the number of
        elements in the ground set.

        EXAMPLES::

            sage: m = PerfectMatching([(-3, 1), (2, 4), (-2, 7)]); m.size()
            6
        """
        return 2*len(self.value)

    def partner(self, x):
        r"""
        Returns the element in the same pair than ``x`` in the matching ``self``.

        EXAMPLES::

            sage: m = PerfectMatching([(-3, 1), (2, 4), (-2, 7)]); m.partner(4)
            2
            sage: n = PerfectMatching([('c','b'),('d','f'),('e','a')])
            sage: n.partner('c')
            'b'
        """
        for i in range(self.size()):
            if self.value[i][0] == x:
                return self.value[i][1]
            if self.value[i][1] == x:
                return self.value[i][0]
        raise ValueError("%s in not an element of the %s" % (x, self))

    def conjugate_by_permutation(self, p):
        r"""
        Returns the conjugate of the perfect matching ``self`` by the
        permutation ``p`` of the ground set.

        EXAMPLE::

            sage: m = PerfectMatching([(1,4),(2,6),(3,5)])
            sage: m.conjugate_by_permutation(Permutation([4,1,5,6,3,2]))
            [(4, 6), (1, 2), (5, 3)]

        TEST::

            sage: PerfectMatching([]).conjugate_by_permutation(Permutation([]))
            []
        """
        return self.parent()([tuple(map(p, t)) for t in self.value])

    def loops_iterator(self, other=None):
        r"""
        INPUT:

             - ``other`` -- a perfect matching of the same set of ``self``.
               (if the second argument is empty, the method :meth:`an_element` is
               called on the parent of the first)

        OUTPUT:

            If we draw the two perfect matchings simultaneously as edges of a
            graph, the graph obtained is a union of cycles of even lengths.
            The function returns an iterator for these cycles (each cycle is
            given as a list).

        EXAMPLES::

            sage: o = PerfectMatching([(1, 7), (2, 4), (3, 8), (5, 6)])
            sage: p = PerfectMatching([(1, 6), (2, 7), (3, 4), (5, 8)])
            sage: it = o.loops_iterator(p)
            sage: next(it)
            [1, 7, 2, 4, 3, 8, 5, 6]
            sage: next(it)
            Traceback (most recent call last):
            ...
            StopIteration
        """
        if other is None:
            other = self.parent().an_element()
        elif self.parent() != other.parent():
            s = "%s is not a matching of the ground set of %s" % (other, self)
            raise ValueError(s)
        remain = flatten(self.value)
        while len(remain) > 0:
            a = remain.pop(0)
            b = self.partner(a)
            remain.remove(b)
            loop = [a, b]
            c = other.partner(b)
            while c != a:
                b = self.partner(c)
                remain.remove(c)
                loop.append(c)
                remain.remove(b)
                loop.append(b)
                c = other.partner(b)
            yield loop

    def loops(self, other=None):
        r"""
        INPUT:

             - ``other`` -- a perfect matching of the same set of ``self``.
               (if the second argument is empty, the method :meth:`an_element` is
               called on the parent of the first)

        OUTPUT:

            If we draw the two perfect matchings simultaneously as edges of a
            graph, the graph obtained is a union of cycles of even lengths.
            The function returns the list of these cycles (each cycle is given
            as a list).

        EXAMPLES::

            sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')])
            sage: n = PerfectMatching([('a','b'),('d','f'),('e','c')])
            sage: m.loops(n)
            [['a', 'e', 'c', 'b'], ['d', 'f']]

            sage: o = PerfectMatching([(1, 7), (2, 4), (3, 8), (5, 6)])
            sage: p = PerfectMatching([(1, 6), (2, 7), (3, 4), (5, 8)])
            sage: o.loops(p)
            [[1, 7, 2, 4, 3, 8, 5, 6]]
        """
        return list(self.loops_iterator(other))

    def loop_type(self, other=None):
        r"""
        INPUT:

             - ``other`` -- a perfect matching of the same set of ``self``.
               (if the second argument is empty, the method :meth:`an_element` is
               called on the parent of the first)

        OUTPUT:

            If we draw the two perfect matchings simultaneously as edges of a
            graph, the graph obtained is a union of cycles of even
            lengths. The function returns the ordered list of the semi-length
            of these cycles (considered as a partition)

        EXAMPLES::

            sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')])
            sage: n = PerfectMatching([('a','b'),('d','f'),('e','c')])
            sage: m.loop_type(n)
            [2, 1]

        TESTS::

            sage: m = PerfectMatching([]); m.loop_type()
            []
        """
        return Partition(reversed(
                sorted([len(l)//2 for l in self.loops_iterator(other)])))

    def number_of_loops(self, other=None):
        r"""
        INPUT:

            - ``other`` -- a perfect matching of the same set of ``self``.
              (if the second argument is empty, the method :meth:`an_element` is
              called on the parent of the first)

        OUTPUT:

            If we draw the two perfect matchings simultaneously as edges of a
            graph, the graph obtained is a union of cycles of even lengths.
            The function returns their numbers.

        EXAMPLES::

            sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')])
            sage: n = PerfectMatching([('a','b'),('d','f'),('e','c')])
            sage: m.number_of_loops(n)
            2
        """
        c = Integer(0)
        one = Integer(1)
        for _ in self.loops_iterator(other):
            c += one
        return c

    def crossings_iterator(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns an iterator over the pairs of
            crossing lines (as a line correspond to a pair, the iterator
            produces pairs of pairs).

        EXAMPLES::

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (2, 8), (4, 7), (5, 6)]
            sage: it = n.crossings_iterator();
            sage: next(it)
            ((1, 3), (2, 8))
            sage: next(it)
            Traceback (most recent call last):
            ...
            StopIteration
        """
        x = self.value[:]
        if len(x) == 0:
            return
        (i, j) = x.pop(0)
        for (a, b) in x:
            # if (i<a<j<b) or (i<b<j<a) or (j<a<i<b) or (j<b<i<a) or (
            #        a<i<b<j) or (a<j<b<i) or (b<i<a<j) or (b<j<a<i):
            labij = sorted([a, b, i, j])
            posij = sorted([labij.index(i), labij.index(j)])
            if posij == [0, 2] or posij == [1, 3]:
                yield ((i, j), (a, b))
        for cr in PerfectMatchings(flatten(x))(x).crossings_iterator():
            yield cr

    def crossings(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns the list of the pairs of
            crossing lines (as a line correspond to a pair, it returns a list
            of pairs of pairs).

        EXAMPLES::

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (2, 8), (4, 7), (5, 6)]
            sage: n.crossings()
            [((1, 3), (2, 8))]

        TESTS::

            sage: m = PerfectMatching([]); m.crossings()
            []
        """
        return list(self.crossings_iterator())

    def number_of_crossings(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns the number the pairs of crossing
            lines.

        EXAMPLES::

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (2, 8), (4, 7), (5, 6)]
            sage: n.number_of_crossings()
            1
        """
        c = Integer(0)
        one = Integer(1)
        for _ in self.crossings_iterator():
            c += one
        return c

    def is_non_crossing(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns ``True`` if the picture obtained
            this way has no crossings.

        EXAMPLES::

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (2, 8), (4, 7), (5, 6)]
            sage: n.is_non_crossing()
            False
            sage: PerfectMatching([(1, 4), (2, 3), (5, 6)]).is_non_crossing()
            True
        """
        it = self.crossings_iterator()
        try:
            next(it)
        except StopIteration:
            return True
        else:
            return False

    def nestings_iterator(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns an iterator over the pairs of
            nesting lines (as a line correspond to a pair, the iterator
            produces pairs of pairs).

        EXAMPLES::

            sage: n = PerfectMatching([(1, 6), (2, 7), (3, 5), (4, 8)])
            sage: it = n.nestings_iterator();
            sage: next(it)
            ((1, 6), (3, 5))
            sage: next(it)
            ((2, 7), (3, 5))
            sage: next(it)
            Traceback (most recent call last):
            ...
            StopIteration
        """
        x = self.value[:]
        if len(x) == 0:
            return
        (i, j) = x.pop(0)
        for (a, b) in x:
            # if (i<a<j<b) or (i<b<j<a) or (j<a<i<b) or (j<b<i<a) or (
            #        a<i<b<j) or (a<j<b<i) or (b<i<a<j) or (b<j<a<i):
            labij = sorted([a, b, i, j])
            posij = sorted([labij.index(i), labij.index(j)])
            if posij == [0, 3] or posij == [1, 2]:
                yield ((i, j), (a, b))
        for nest in PerfectMatchings(flatten(x))(x).nestings_iterator():
            yield nest

    def nestings(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns the list of the pairs of
            nesting lines (as a line correspond to a pair, it returns a list
            of pairs of pairs).

        EXAMPLES::

            sage: m = PerfectMatching([(1, 6), (2, 7), (3, 5), (4, 8)])
            sage: m.nestings()
            [((1, 6), (3, 5)), ((2, 7), (3, 5))]

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (2, 8), (4, 7), (5, 6)]
            sage: n.nestings()
            [((2, 8), (4, 7)), ((2, 8), (5, 6)), ((4, 7), (5, 6))]

        TESTS::

            sage: m = PerfectMatching([]); m.nestings()
            []
        """
        return list(self.nestings_iterator())

    def number_of_nestings(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns the number the pairs of nesting
            lines.

        EXAMPLES::

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (2, 8), (4, 7), (5, 6)]
            sage: n.number_of_nestings()
            3
        """
        c = Integer(0)
        one = Integer(1)
        for _ in self.nestings_iterator():
            c += one
        return c

    def is_non_nesting(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns ``True`` if the picture obtained
            this way has no nestings.

        EXAMPLES::

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (2, 8), (4, 7), (5, 6)]
            sage: n.is_non_nesting()
            False
            sage: PerfectMatching([(1, 3), (2, 5), (4, 6)]).is_non_nesting()
            True
        """
        it = self.nestings_iterator()
        try:
            next(it)
        except StopIteration:
            return True
        else:
            return False

    def Weingarten_function(self, d, other=None):
        r"""
        Returns the Weingarten function of two pairings.

        This function is the value of some integrals over the orthogonal
        groups `O_N`.  With the convention of [CM]_, the method returns
        `Wg^{O(d)}(other,self)`.

        EXAMPLES::

            sage: var('N')
            N
            sage: m = PerfectMatching([(1,3),(2,4)])
            sage: n = PerfectMatching([(1,2),(3,4)])
            sage: factor(m.Weingarten_function(N,n))
            -1/((N + 2)*(N - 1)*N)
        """
        if other is None:
            other = self.parent().an_element()
        W = self.parent().Weingarten_matrix(d)
        return W[other.rank()][self.rank()]

    def to_graph(self):
        r"""
        Returns the graph corresponding to the perfect matching.

        OUTPUT:

            The realization of ``self`` as a graph.

        EXAMPLES::

            sage: PerfectMatching([[1,3], [4,2]]).to_graph().edges(labels=False)
            [(1, 3), (2, 4)]
            sage: PerfectMatching([[1,4], [3,2]]).to_graph().edges(labels=False)
            [(1, 4), (2, 3)]
            sage: PerfectMatching([]).to_graph().edges(labels=False)
            []
        """
        from sage.graphs.graph import Graph
        G = Graph()
        for a, b in self.value:
            G.add_edge((a, b))
        return G

    @combinatorial_map(name='to permutation')
    def to_permutation(self):
        r"""
        Returns the permutation corresponding to the perfect matching.

        OUTPUT:

            The realization of ``self`` as a permutation.

        EXAMPLES::

            sage: PerfectMatching([[1,3], [4,2]]).to_permutation()
            [3, 4, 1, 2]
            sage: PerfectMatching([[1,4], [3,2]]).to_permutation()
            [4, 3, 2, 1]
            sage: PerfectMatching([]).to_permutation()
            []
        """
        from sage.combinat.permutation import Permutation
        return Permutation(self.value)

    def to_non_crossing_set_partition(self):
        r"""
        Returns the noncrossing set partition (on half as many elements) 
        corresponding to the perfect matching if the perfect matching is 
        noncrossing, and otherwise gives an error.

        OUTPUT:

            The realization of ``self`` as a noncrossing set partition.

        EXAMPLES::

            sage: PerfectMatching([[1,3], [4,2]]).to_non_crossing_set_partition()
            Traceback (most recent call last):
            ...
            ValueError: matching must be non-crossing
            sage: PerfectMatching([[1,4], [3,2]]).to_non_crossing_set_partition()
            {{1, 2}}
            sage: PerfectMatching([]).to_non_crossing_set_partition()
            {}
        """
        from sage.combinat.set_partition import SetPartition        
        if not self.is_non_crossing():
            raise ValueError("matching must be non-crossing")
        else:
            perm = self.to_permutation()
            perm2 = Permutation([perm[2 * i] // 2
                                 for i in range(len(perm) // 2)])
        return SetPartition(perm2.cycle_tuples())


class PerfectMatchings(UniqueRepresentation, Parent):
    r"""
    Class of perfect matchings of a ground set. At the creation, the set
    can be given as any iterable object. If the argument is an integer `n`, it
    will be transformed into `[1 .. n]`::

        sage: M = PerfectMatchings(6);M
        Set of perfect matchings of {1, 2, 3, 4, 5, 6}
        sage: PerfectMatchings([-1, -3, 1, 2])
        Set of perfect matchings of {1, 2, -3, -1}

    One can ask for the list, the cardinality or an element of a set of
    perfect matching::

        sage: PerfectMatchings(4).list()
        [[(1, 2), (3, 4)], [(1, 3), (2, 4)], [(1, 4), (2, 3)]]
        sage: PerfectMatchings(8).cardinality()
        105
        sage: M = PerfectMatchings(('a', 'e', 'b', 'f', 'c', 'd'))
        sage: M.an_element()
        [('a', 'b'), ('f', 'e'), ('c', 'd')]
        sage: all([PerfectMatchings(i).an_element() in PerfectMatchings(i)
        ...        for i in range(2,11,2)])
        True

    TESTS::

        sage: PerfectMatchings(5).list()
        []
        sage: TestSuite(PerfectMatchings(6)).run()
        sage: TestSuite(PerfectMatchings([])).run()
    """

    @staticmethod
    def _parse_input(objects):
        r"""
        This function tries to recognize the argument and to transform into a
        set. It is not meant to be called manually, but only as the first of
        the creation of an enumerated set of ``PerfectMatchings``.

        EXAMPLES::

            sage: PerfectMatchings._parse_input(4)
            {1, 2, 3, 4}
            sage: PerfectMatchings._parse_input(['a','b','c','e'])
            {'a', 'c', 'b', 'e'}
        """
        # if the argument is a python int n, we replace it by the list [1 .. n]
        if isinstance(objects, int):
            objects = range(1, objects+1)
        # same thing if the argument is a sage integer.
        elif isinstance(objects, Integer):
            objects = range(1, objects+1)
        # Finally, if it is iterable, we return the corresponding set.
        # Note that it is important to return a hashable object here (in
        # particular, NOT A LIST), see comment below.
        if not hasattr(objects, '__iter__'):
            raise ValueError("do not know how to construct a set of matchings from %s (it must be iterable)")
        return Set(objects)

    @staticmethod
    def __classcall__(cls, objects):
        r"""
        This function is called automatically when the user want to
        create an enumerated set of PerfectMatchings.

        EXAMPLES::

            sage: M = PerfectMatchings(6);M
            Set of perfect matchings of {1, 2, 3, 4, 5, 6}
            sage: PerfectMatchings([-1, -3, 1, 2])
            Set of perfect matchings of {1, 2, -3, -1}

        If one has already created a set of perfect matchings of the same set,
        it does not create a new object, but returns the already existing
        one::

            sage: N = PerfectMatchings((2, 3, 5, 4, 1, 6))
            sage: N is M
            True
        """
        #we call the constructor of an other class, which will
        #    - check if the object has already been constructed (so the
        # second argument, i.e. the output of _parse_input, must be hashable)
        #    - look for a place in memory and call the __init__ function
        return super(PerfectMatchings, cls).__classcall__(
            cls, cls._parse_input(objects))

    def __init__(self, objects):
        r"""
        See :meth:`__classcall__`

        TEST::

            sage: M = PerfectMatchings(6)
            sage: TestSuite(M).run()
        """
        self._objects = objects
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        r"""
        Returns a description of ``self``.

        EXAMPLES::

            sage: PerfectMatchings([-1, -3, 1, 2])
            Set of perfect matchings of {1, 2, -3, -1}
        """
        return "Set of perfect matchings of %s" % self._objects

    def __iter__(self):
        r"""
        Returns an iterator for the elements of ``self``.

        EXAMPLES::

            sage: PerfectMatchings(4).list()
            [[(1, 2), (3, 4)], [(1, 3), (2, 4)], [(1, 4), (2, 3)]]
        """
        if len(self._objects) == 0:
            yield self([])
        elif len(self._objects) == 1:
            pass
        else:
            l = list(self._objects)
            a = l.pop(0)
            for i in range(len(l)):
                obj_rest = l[:]
                b = obj_rest.pop(i)
                for p in PerfectMatchings(obj_rest):
                    yield self([(a, b)]+p.value)

    def __contains__(self, x):
        r"""
        Tests if ``x`` is an element of ``self``.

        EXAMPLES::

            sage: m = PerfectMatching([(1,2),(4,3)])
            sage: m in PerfectMatchings(4)
            True
            sage: m in PerfectMatchings((0, 1, 2, 3))
            False
            sage: all([m in PerfectMatchings(6) for m in PerfectMatchings(6)])
            True

        Note that the class of ``x`` does not need to be ``PerfectMatching``:
        if the data defines a perfect matching of the good set, the function
        returns ``True``::

            sage: [(1, 4), (2, 3)] in PerfectMatchings(4)
            True
            sage: [(1, 3, 6), (2, 4), (5,)] in PerfectMatchings(6)
            False
            sage: [('a', 'b'), ('a', 'c')] in PerfectMatchings(
            ...        ('a', 'b', 'c', 'd'))
            False
        """
        if not isinstance(x, PerfectMatching):
            try:
                x = PerfectMatching(x)
            except ValueError:
                return False
        return x.parent() is self

    def cardinality(self):
        r"""
        Returns the cardinality of the set of perfect matching ``self``

        This is `1*3*5*...*(2n-1)`, where `2n` is the size of the ground set.

        EXAMPLES::

            sage: PerfectMatchings(8).cardinality()
            105
        """
        n = len(self._objects)
        if n % 2 == 1:
            return 0
        else:
            return Integer(prod(i for i in range(n) if i % 2 == 1))

    def random_element(self):
        r"""
        Returns a random element of self.

        EXAMPLES::

            sage: M = PerfectMatchings(('a', 'e', 'b', 'f', 'c', 'd'))
            sage: M.an_element()
            [('a', 'b'), ('f', 'e'), ('c', 'd')]
            sage: all([PerfectMatchings(2*i).an_element() in PerfectMatchings(2*i)
            ...        for i in range(2,11,2)])
            True

        TESTS::

            sage: p = PerfectMatchings(13).random_element()
            Traceback (most recent call last):
            ...
            ValueError: there is no perfect matching on an odd number of elements

        """
        n = len(self._objects)

        if n % 2 == 1:
            raise ValueError("there is no perfect matching on an odd number of elements")

        k = n//2

        from sage.combinat.permutation import Permutations
        p = Permutations(n).random_element()

        return self([(self._objects[p[2*i]-1], self._objects[p[2*i+1]-1]) for i in range(k)])

    an_element = random_element

    @cached_method
    def Weingarten_matrix(self, N):
        r"""
        Returns the Weingarten matrix corresponding to the set of
        PerfectMatchings ``self``.

        It is a useful theoretical tool to compute polynomial integral
        over the orthogonal group `O_N` (see [CM]_).

        EXAMPLES::

            sage: M = PerfectMatchings(4).Weingarten_matrix(var('N'))
            sage: N*(N-1)*(N+2)*M.apply_map(factor)
            [N + 1    -1    -1]
            [   -1 N + 1    -1]
            [   -1    -1 N + 1]
        """
        G = Matrix([[N**(p1.number_of_loops(p2)) for p1 in self]
                    for p2 in self])
        return G**(-1)

    Element = PerfectMatching
