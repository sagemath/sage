# -*- coding: utf-8 -*-
r"""
Linear Extensions of Posets
"""
#*****************************************************************************
#       Copyright (C) 2012 Anne Schilling <anne at math.ucdavis.edu>
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

from sage.rings.rational_field import QQ
from sage.categories.posets import Posets
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.graphs.digraph import DiGraph
import sage.graphs.linearextensions
from sage.combinat.posets.hasse_diagram import HasseDiagram
from sage.combinat.posets.posets import Poset
from sage.combinat.posets.elements import PosetElement
from sage.combinat.permutation import Permutation
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.graphs.dot2tex_utils import have_dot2tex
from sage.structure.list_clone import ClonableArray

class LinearExtensionOfPoset(ClonableArray):
    r"""
    A linear extension of a finite poset `P` of size `n` is a total
    ordering `\pi := \pi_0 \pi_1 \ldots \pi_{n-1}` of its elements
    such that `i<j` whenever `\pi_i < \pi_j` in the poset `P`.

    When the elements of `P` are indexed by `\{1,2,\ldots,n\}`, `\pi`
    denotes a permutation of the elements of `P` in one-line notation.

    INPUT:

    - ``linear_extension`` -- a list of the elements of `P`
    - ``poset`` -- the underlying poset `P`

    .. seealso:: :class:`Poset`, :class:`LinearExtensionsOfPosets`

    EXAMPLES::

        sage: P = Poset(([1,2,3,4], [[1,3],[1,4],[2,3]]), linear_extension=True, facade = False)
        sage: p = P.linear_extension([1,4,2,3]); p
        [1, 4, 2, 3]
        sage: p.parent()
        The set of all linear extensions of Finite poset containing 4 elements
        sage: p[0], p[1], p[2], p[3]
        (1, 4, 2, 3)

    Following Schützenberger and later Haiman and
    Malvenuto-Reutenauer, Stanley [Stanley2009]_ defined a promotion
    and evacuation operator on any finite poset `P` using operators
    `\tau_i` on the linear extensions of `P`::

        sage: p.promotion()
        [1, 2, 3, 4]
        sage: Q = p.promotion().to_poset()
        sage: Q.cover_relations()
        [[1, 3], [1, 4], [2, 3]]
        sage: Q == P
        True

        sage: p.promotion(3)
        [1, 4, 2, 3]
        sage: Q = p.promotion(3).to_poset()
        sage: Q == P
        False
        sage: Q.cover_relations()
        [[1, 2], [1, 4], [3, 4]]

    REFERENCES:

        .. [Stanley2009] Richard Stanley,
           *Promotion and evacuation*,
           Electron. J. Combin. 16 (2009), no. 2, Special volume in honor of Anders Björner,
           Research Paper 9, 24 pp.
    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, linear_extension, poset):
        r"""
        Implements the shortcut ``LinearExtensionOfPoset(linear_extension, poset)`` to ``LinearExtensionsOfPoset(poset)(linear_extension)``

        INPUT:

        - ``linear_extension`` -- a list of elements of ``poset``
        - ``poset`` -- a finite poset

        .. todo:: check whether this method is still useful

        TESTS::

            sage: from sage.combinat.posets.linear_extensions import LinearExtensionOfPoset
            sage: P = Poset(([1,2,3,4], [[1,3],[1,4],[2,3]]))
            sage: p = LinearExtensionOfPoset([1,4,2,3], P)
            sage: p.parent()
            The set of all linear extensions of Finite poset containing 4 elements
            sage: type(p)
            <class 'sage.combinat.posets.linear_extensions.LinearExtensionsOfPoset_with_category.element_class'>
            sage: p.poset()
            Finite poset containing 4 elements
            sage: TestSuite(p).run()

            sage: LinearExtensionOfPoset([4,3,2,1], P)
            Traceback (most recent call last):
            ...
            ValueError: [4, 3, 2, 1] is not a linear extension of Finite poset containing 4 elements
        """
        if isinstance(linear_extension, cls):
            return linear_extension
        return LinearExtensionsOfPoset(poset)(linear_extension)

    def check(self):
        r"""
        Checks whether ``self`` is indeed a linear extension of the underlying poset.

        TESTS::

            sage: P = Poset(([1,2,3,4], [[1,3],[1,4],[2,3]]))
            sage: P.linear_extension([1,4,2,3])
            [1, 4, 2, 3]
            sage: P.linear_extension([4,3,2,1])
            Traceback (most recent call last):
            ...
            ValueError: [4, 3, 2, 1] is not a linear extension of Finite poset containing 4 elements
        """
        P = self.parent().poset()
        if not P.is_linear_extension(self):
            raise ValueError, "%s is not a linear extension of %s"%(self, P)

    def poset(self):
        r"""
        Returns the underlying original poset.

        EXAMPLES::

            sage: P = Poset(([1,2,3,4], [[1,2],[2,3],[1,4]]))
            sage: p = P.linear_extension([1,2,4,3])
            sage: p.poset()
            Finite poset containing 4 elements
        """
        return self.parent().poset()

    def _latex_(self):
        r"""
        Returns the latex string for ``self``.

        EXAMPLES::

            sage: P = Poset(([1,2,3,4], [[1,3],[1,4],[2,3]]))
            sage: p = P.linear_extension([1,2,3,4])
            sage: p._latex_()
            '\\mathtt{(1, 2, 3, 4)}'
        """
        return "\\mathtt{"+str(tuple(self))+"}"

    def to_poset(self):
        r"""
        Returns the poset associated to the linear extension ``self``.

        This method returns the poset obtained from the original poset
        `P` by relabelling the 'i'-th element of ``self`` to the
        `i`-th element of the original poset, while keeping the linear
        extension of the original poset.

        For a poset with default linear extension `1,\dots,n`,
        ``self`` can be interpreted as a permutation, and the
        relabelling is done according to the inverse of this
        permutation.

        EXAMPLES::

            sage: P = Poset(([1,2,3,4], [[1,2],[1,3],[3,4]]), facade = False)
            sage: p = P.linear_extension([1,3,4,2])
            sage: Q = p.to_poset(); Q
            Finite poset containing 4 elements
            sage: P == Q
            False

        The default linear extension remains the same::

            sage: list(P)
            [1, 2, 3, 4]
            sage: list(Q)
            [1, 2, 3, 4]

        But the relabelling can be seen on cover relations::

            sage: P.cover_relations()
            [[1, 2], [1, 3], [3, 4]]
            sage: Q.cover_relations()
            [[1, 2], [1, 4], [2, 3]]

            sage: p = P.linear_extension([1,2,3,4])
            sage: Q = p.to_poset()
            sage: P == Q
            True
        """
        P = self.parent().poset()
        old = [P.unwrap(x) for x in self]
        new = [P.unwrap(x) for x in P]
        relabelling = dict(zip(old,new))
        return P.relabel(relabelling).with_linear_extension(new)

    def tau(self, i):
        r"""
        Returns the operator `\tau_i` on linear extensions ``self`` of a poset.

        INPUT:

        - `i` -- an integer between `1` and `n-1`, where `n` is the cardinality of the poset.

        The operator `\tau_i` on a linear extension `\pi` of a poset
        `P` interchanges positions `i` and `i+1` if the result is
        again a linear extension of `P`, and otherwise acts
        trivially. For more details, see [Stanley2009]_.

        EXAMPLES::

            sage: P = Poset(([1,2,3,4], [[1,3],[1,4],[2,3]]), linear_extension=True)
            sage: L = P.linear_extensions()
            sage: l = L.an_element(); l
            [1, 2, 3, 4]
            sage: l.tau(1)
            [2, 1, 3, 4]
            sage: for p in L:
            ...       for i in range(1,4):
            ...           print i, p, p.tau(i)
            ...
            1 [1, 2, 3, 4] [2, 1, 3, 4]
            2 [1, 2, 3, 4] [1, 2, 3, 4]
            3 [1, 2, 3, 4] [1, 2, 4, 3]
            1 [1, 2, 4, 3] [2, 1, 4, 3]
            2 [1, 2, 4, 3] [1, 4, 2, 3]
            3 [1, 2, 4, 3] [1, 2, 3, 4]
            1 [1, 4, 2, 3] [1, 4, 2, 3]
            2 [1, 4, 2, 3] [1, 2, 4, 3]
            3 [1, 4, 2, 3] [1, 4, 2, 3]
            1 [2, 1, 3, 4] [1, 2, 3, 4]
            2 [2, 1, 3, 4] [2, 1, 3, 4]
            3 [2, 1, 3, 4] [2, 1, 4, 3]
            1 [2, 1, 4, 3] [1, 2, 4, 3]
            2 [2, 1, 4, 3] [2, 1, 4, 3]
            3 [2, 1, 4, 3] [2, 1, 3, 4]

        TESTS::

            sage: type(l.tau(1))
            <class 'sage.combinat.posets.linear_extensions.LinearExtensionsOfPoset_with_category.element_class'>
            sage: l.tau(2) == l
            True
        """
        P = self.poset()
        a = self[i-1]
        b = self[i  ]
        if P.lt(a,b) or P.lt(b,a):
            return self
        with self.clone() as q:
                q[i-1] = b
                q[i  ] = a
        return q

    def promotion(self, i=1):
        r"""
        Computes the (generalized) promotion on the linear extension of a poset.

        INPUT:

        - `i` -- an integer between `1` and `n-1`, where `n` is the cardinality of the poset (default: `1`)

        The `i`-th generalized promotion operator `\partial_i` on a linear extension
        `\pi` is defined as `\pi \tau_i \tau_{i+1} \cdots \tau_{n-1}`, where `n` is the
        size of the linear extension (or size of the underlying poset).

        For more details see [Stanley2009]_.

        .. seealso:: :meth:`tau`, :meth:`evacuation`

        EXAMPLES::

            sage: P = Poset(([1,2,3,4,5,6,7], [[1,2],[1,4],[2,3],[2,5],[3,6],[4,7],[5,6]]))
            sage: p = P.linear_extension([1,2,3,4,5,6,7])
            sage: q = p.promotion(4); q
            [1, 2, 3, 5, 6, 4, 7]
            sage: p.to_poset() == q.to_poset()
            False
            sage: p.to_poset().is_isomorphic(q.to_poset())
            True
        """
        for j in range(i,len(self)):
            self = self.tau(j)
        return self

    def evacuation(self):
        r"""
        Computes evacuation on the linear extension of a poset.

        Evacuation on a linear extension `\pi` of length `n` is defined as
        `\pi (\tau_1 \cdots \tau_{n-1}) (\tau_1 \cdots \tau_{n-2}) \cdots (\tau_1)`.
        For more details see [Stanley2009]_.

        .. seealso:: :meth:`tau`, :meth:`promotion`

        EXAMPLES::

            sage: P = Poset(([1,2,3,4,5,6,7], [[1,2],[1,4],[2,3],[2,5],[3,6],[4,7],[5,6]]))
            sage: p = P.linear_extension([1,2,3,4,5,6,7])
            sage: p.evacuation()
            [1, 4, 2, 3, 7, 5, 6]
            sage: p.evacuation().evacuation() == p
            True
        """
        for i in reversed(range(1,len(self)+1)):
            for j in range(1,i):
                self = self.tau(j)
        return self

class LinearExtensionsOfPoset(UniqueRepresentation, Parent):
    """
    The set of all linear extensions of a finite poset

    INPUT:

    - ``poset`` -- a poset `P` of size `n`
    - ``facade`` -- a boolean (default: False)

    .. seealso::

        - :meth:`sage.combinat.posets.posets.FinitePoset.linear_extensions`
        - :class:`sage.graphs.linearextensions.LinearExtensions`

    EXAMPLES::

        sage: elms = [1,2,3,4]
        sage: rels = [[1,3],[1,4],[2,3]]
        sage: P = Poset((elms, rels), linear_extension=True)
        sage: L = P.linear_extensions(); L
        The set of all linear extensions of Finite poset containing 4 elements
        sage: L.cardinality()
        5
        sage: L.list()
        [[1, 2, 3, 4], [1, 2, 4, 3], [1, 4, 2, 3], [2, 1, 3, 4], [2, 1, 4, 3]]
        sage: L.an_element()
        [1, 2, 3, 4]
        sage: L.poset()
        Finite poset containing 4 elements
    """

    @staticmethod
    def __classcall_private__(cls, poset, facade=False):
        r"""
        Straighten arguments before unique representation.

        TESTS::

            sage: from sage.combinat.posets.linear_extensions import LinearExtensionsOfPoset
            sage: P = Poset(([1,2],[[1,2]]))
            sage: L = LinearExtensionsOfPoset(P)
            sage: type(L)
            <class 'sage.combinat.posets.linear_extensions.LinearExtensionsOfPoset_with_category'>
            sage: L is LinearExtensionsOfPoset(P,facade=False)
            True
        """
        return super(LinearExtensionsOfPoset, cls).__classcall__(cls, poset, facade=facade)

    def __init__(self, poset, facade):
        """
        TESTS::

            sage: from sage.combinat.posets.linear_extensions import LinearExtensionsOfPoset
            sage: P = Poset(([1,2,3],[[1,2],[1,3]]))
            sage: L = P.linear_extensions()
            sage: L is LinearExtensionsOfPoset(P)
            True
            sage: L._poset is P
            True
            sage: L._linear_extensions_of_hasse_diagram
            Linear extensions of Hasse diagram of a poset containing 3 elements
            sage: TestSuite(L).run()

            sage: P = Poset((divisors(15), attrcall("divides")))
            sage: L = P.linear_extensions()
            sage: TestSuite(L).run()

            sage: P = Poset((divisors(15), attrcall("divides")), facade=True)
            sage: L = P.linear_extensions()
            sage: TestSuite(L).run()

            sage: L = P.linear_extensions(facade = True)
            sage: TestSuite(L).run(skip="_test_an_element")
        """
        self._poset = poset
        self._linear_extensions_of_hasse_diagram = sage.graphs.linearextensions.LinearExtensions(poset._hasse_diagram)
        self._is_facade = facade
        if facade:
            facade = (list,)
        Parent.__init__(self, category = FiniteEnumeratedSets(), facade=facade)

    def _repr_(self):
        """
        TESTS::

            sage: P = Poset(([1,2,3],[[1,2],[1,3]]))
            sage: P.linear_extensions()
            The set of all linear extensions of Finite poset containing 3 elements
        """
        return "The set of all linear extensions of %s"%(self._poset)

    def poset(self):
        r"""
        Returns the underlying original poset.

        EXAMPLES::

            sage: P = Poset(([1,2,3,4], [[1,2],[2,3],[1,4]]))
            sage: L = P.linear_extensions()
            sage: L.poset()
            Finite poset containing 4 elements
        """
        return self._poset

    def __iter__(self):
        r"""
        Iterates through the linear extensions of the underlying poset.

        EXAMPLES::

            sage: elms = [1,2,3,4]
            sage: rels = [[1,3],[1,4],[2,3]]
            sage: P = Poset((elms, rels), linear_extension=True)
            sage: L = P.linear_extensions()
            sage: list(L)
            [[1, 2, 3, 4], [1, 2, 4, 3], [1, 4, 2, 3], [2, 1, 3, 4], [2, 1, 4, 3]]
        """
        vertex_to_element = self._poset._vertex_to_element
        for lin_ext in self._linear_extensions_of_hasse_diagram:
            yield self._element_constructor_(map(vertex_to_element,lin_ext))

    def __contains__(self, obj):
        """
        Membership testing

        EXAMPLES::

            sage: P = Poset((divisors(12), attrcall("divides")), facade=True, linear_extension=True)
            sage: P.list()
            [1, 2, 3, 4, 6, 12]
            sage: L = P.linear_extensions()
            sage: L([1, 2, 4, 3, 6, 12]) in L
            True
            sage: [1, 2, 4, 3, 6, 12] in L
            False

            sage: L = P.linear_extensions(facade=True)
            sage: [1, 2, 4, 3, 6, 12] in L
            True
            sage: [1, 3, 2, 6, 4, 12] in L
            True
            sage: [1, 3, 6, 2, 4, 12] in L
            False

            sage: [p for p in Permutations(list(P)) if list(p) in L]
            [[1, 2, 3, 4, 6, 12], [1, 2, 3, 6, 4, 12], [1, 2, 4, 3, 6, 12], [1, 3, 2, 4, 6, 12], [1, 3, 2, 6, 4, 12]]

        """
        if not self._is_facade:
            return super(LinearExtensionsOfPoset, self).__contains__(obj)
        return (isinstance(obj, (list, tuple)) and
                self.poset().is_linear_extension(obj))

    def markov_chain_digraph(self, action = 'promotion', labeling = 'identity'):
        r"""
        Returns the digraph of the action of generalized promotion or tau on ``self``

        INPUT:

        - ``action`` -- 'promotion' or 'tau' (default: 'promotion')
        - ``labeling`` -- 'identity' or 'source' (default: 'identity')

        .. todo::

            - generalize this feature by accepting a family of operators as input
            - move up in some appropriate category

        This method creates a graph with vertices being the linear extensions of a given finite
        poset and an edge from `\pi` to `\pi'` if `\pi' = \pi \partial_i` where `\partial_i` is
        the promotion operator (see :meth:`promotion`) if ``action`` is set to ``promotion``
        and `\tau_i` (see :meth:`tau`) if ``action`` is set to ``tau``. The label of the edge
        is `i` (resp. `\pi_i`) if ``labeling`` is set to ``identity`` (resp. ``source``).

        EXAMPLES::

            sage: P = Poset(([1,2,3,4], [[1,3],[1,4],[2,3]]), linear_extension = True)
            sage: L = P.linear_extensions()
            sage: G = L.markov_chain_digraph(); G
            Looped multi-digraph on 5 vertices
            sage: sorted(G.vertices(), key = repr)
            [[1, 2, 3, 4], [1, 2, 4, 3], [1, 4, 2, 3], [2, 1, 3, 4], [2, 1, 4, 3]]
            sage: sorted(G.edges(), key = repr)
            [([1, 2, 3, 4], [1, 2, 3, 4], 4), ([1, 2, 3, 4], [1, 2, 4, 3], 2), ([1, 2, 3, 4], [1, 2, 4, 3], 3),
            ([1, 2, 3, 4], [2, 1, 4, 3], 1), ([1, 2, 4, 3], [1, 2, 3, 4], 3), ([1, 2, 4, 3], [1, 2, 4, 3], 4),
            ([1, 2, 4, 3], [1, 4, 2, 3], 2), ([1, 2, 4, 3], [2, 1, 3, 4], 1), ([1, 4, 2, 3], [1, 2, 3, 4], 1),
            ([1, 4, 2, 3], [1, 2, 3, 4], 2), ([1, 4, 2, 3], [1, 4, 2, 3], 3), ([1, 4, 2, 3], [1, 4, 2, 3], 4),
            ([2, 1, 3, 4], [1, 2, 4, 3], 1), ([2, 1, 3, 4], [2, 1, 3, 4], 4), ([2, 1, 3, 4], [2, 1, 4, 3], 2),
            ([2, 1, 3, 4], [2, 1, 4, 3], 3), ([2, 1, 4, 3], [1, 4, 2, 3], 1), ([2, 1, 4, 3], [2, 1, 3, 4], 2),
            ([2, 1, 4, 3], [2, 1, 3, 4], 3), ([2, 1, 4, 3], [2, 1, 4, 3], 4)]

            sage: G = L.markov_chain_digraph(labeling = 'source')
            sage: sorted(G.vertices(), key = repr)
            [[1, 2, 3, 4], [1, 2, 4, 3], [1, 4, 2, 3], [2, 1, 3, 4], [2, 1, 4, 3]]
            sage: sorted(G.edges(), key = repr)
            [([1, 2, 3, 4], [1, 2, 3, 4], 4), ([1, 2, 3, 4], [1, 2, 4, 3], 2), ([1, 2, 3, 4], [1, 2, 4, 3], 3),
            ([1, 2, 3, 4], [2, 1, 4, 3], 1), ([1, 2, 4, 3], [1, 2, 3, 4], 4), ([1, 2, 4, 3], [1, 2, 4, 3], 3),
            ([1, 2, 4, 3], [1, 4, 2, 3], 2), ([1, 2, 4, 3], [2, 1, 3, 4], 1), ([1, 4, 2, 3], [1, 2, 3, 4], 1),
            ([1, 4, 2, 3], [1, 2, 3, 4], 4), ([1, 4, 2, 3], [1, 4, 2, 3], 2), ([1, 4, 2, 3], [1, 4, 2, 3], 3),
            ([2, 1, 3, 4], [1, 2, 4, 3], 2), ([2, 1, 3, 4], [2, 1, 3, 4], 4), ([2, 1, 3, 4], [2, 1, 4, 3], 1),
            ([2, 1, 3, 4], [2, 1, 4, 3], 3), ([2, 1, 4, 3], [1, 4, 2, 3], 2), ([2, 1, 4, 3], [2, 1, 3, 4], 1),
            ([2, 1, 4, 3], [2, 1, 3, 4], 4), ([2, 1, 4, 3], [2, 1, 4, 3], 3)]

        The edges of the graph are by default colored using blue for
        edge 1, red for edge 2, green for edge 3, and yellow for edge 4::

            sage: view(G) #optional - dot2tex graphviz

        Alternatively, one may get the graph of the action of the ``tau`` operator::

            sage: G = L.markov_chain_digraph(action='tau'); G
            Looped multi-digraph on 5 vertices
            sage: sorted(G.vertices(), key = repr)
            [[1, 2, 3, 4], [1, 2, 4, 3], [1, 4, 2, 3], [2, 1, 3, 4], [2, 1, 4, 3]]
            sage: sorted(G.edges(), key = repr)
            [([1, 2, 3, 4], [1, 2, 3, 4], 2), ([1, 2, 3, 4], [1, 2, 4, 3], 3), ([1, 2, 3, 4], [2, 1, 3, 4], 1),
            ([1, 2, 4, 3], [1, 2, 3, 4], 3), ([1, 2, 4, 3], [1, 4, 2, 3], 2), ([1, 2, 4, 3], [2, 1, 4, 3], 1),
            ([1, 4, 2, 3], [1, 2, 4, 3], 2), ([1, 4, 2, 3], [1, 4, 2, 3], 1), ([1, 4, 2, 3], [1, 4, 2, 3], 3),
            ([2, 1, 3, 4], [1, 2, 3, 4], 1), ([2, 1, 3, 4], [2, 1, 3, 4], 2), ([2, 1, 3, 4], [2, 1, 4, 3], 3),
            ([2, 1, 4, 3], [1, 2, 4, 3], 1), ([2, 1, 4, 3], [2, 1, 3, 4], 3), ([2, 1, 4, 3], [2, 1, 4, 3], 2)]
            sage: view(G) #optional - dot2tex graphviz

        .. seealso:: :meth:`markov_chain_transition_matrix`, :meth:`promotion`, :meth:`tau`

        TESTS::

            sage: P = Poset(([1,2,3,4], [[1,3],[1,4],[2,3]]), linear_extension = True, facade = True)
            sage: L = P.linear_extensions()
            sage: G = L.markov_chain_digraph(labeling = 'source'); G
            Looped multi-digraph on 5 vertices
        """
        d = dict([x,dict([y,[]] for y in self)] for x in self)
        if action == 'promotion':
            R = range(self.poset().cardinality())
        else:
            R = range(self.poset().cardinality()-1)
        if labeling == 'source':
            for x in self:
                for i in R:
                    child = getattr(x, action)(i+1)
                    d[x][child]+=[self.poset().unwrap(x[i])]
        else:
            for x in self:
                for i in R:
                    child = getattr(x, action)(i+1)
                    d[x][child]+=[i+1]
        G = DiGraph(d)
        if have_dot2tex():
            G.set_latex_options(format="dot2tex", edge_labels = True, color_by_label = {1:"blue", 2:"red", 3:"green", 4:"yellow"})
            #G.set_latex_options(format="dot2tex", edge_labels = True, color_by_label = {1:"green", 2:"blue", 3:"brown", 4:"red"})
        return G

    def markov_chain_transition_matrix(self, action = 'promotion', labeling = 'identity'):
        r"""
        Returns the transition matrix of the Markov chain for the action of generalized promotion or tau on ``self``

        INPUT:

        - ``action`` -- 'promotion' or 'tau' (default: 'promotion')
        - ``labeling`` -- 'identity' or 'source' (default: 'identity')

        This method yields the transition matrix of the Markov chain defined by the action of the generalized
        promotion operator `\partial_i` (resp. `\tau_i`) on the set of linear extensions of a finite poset.
        Here the transition from the linear extension `\pi` to `\pi'`, where `\pi' = \pi \partial_i`
        (resp. `\pi'= \pi \tau_i`) is counted with weight `x_i` (resp. `x_{\pi_i}` if ``labeling`` is set to ``source``).

        EXAMPLES::

            sage: P = Poset(([1,2,3,4], [[1,3],[1,4],[2,3]]), linear_extension = True)
            sage: L = P.linear_extensions()
            sage: L.markov_chain_transition_matrix()
            [-x0 - x1 - x2            x2       x0 + x1             0             0]
            [      x1 + x2 -x0 - x1 - x2             0            x0             0]
            [            0            x1      -x0 - x1             0            x0]
            [            0            x0             0 -x0 - x1 - x2       x1 + x2]
            [           x0             0             0       x1 + x2 -x0 - x1 - x2]

            sage: L.markov_chain_transition_matrix(labeling = 'source')
            [-x0 - x1 - x2            x3       x0 + x3             0             0]
            [      x1 + x2 -x0 - x1 - x3             0            x1             0]
            [            0            x1      -x0 - x3             0            x1]
            [            0            x0             0 -x0 - x1 - x2       x0 + x3]
            [           x0             0             0       x0 + x2 -x0 - x1 - x3]

            sage: L.markov_chain_transition_matrix(action = 'tau')
            [     -x0 - x2            x2             0            x0             0]
            [           x2 -x0 - x1 - x2            x1             0            x0]
            [            0            x1           -x1             0             0]
            [           x0             0             0      -x0 - x2            x2]
            [            0            x0             0            x2      -x0 - x2]

            sage: L.markov_chain_transition_matrix(action = 'tau', labeling = 'source')
            [     -x0 - x2            x3             0            x1             0]
            [           x2 -x0 - x1 - x3            x3             0            x1]
            [            0            x1           -x3             0             0]
            [           x0             0             0      -x1 - x2            x3]
            [            0            x0             0            x2      -x1 - x3]

        .. seealso:: :meth:`markov_chain_digraph`, :meth:`promotion`, :meth:`tau`

        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.matrix.constructor import matrix
        L = self.list()
        n = self.poset().cardinality()
        R = PolynomialRing(QQ, 'x', n)
        x = [R.gen(i) for i in range(n)]
        l = self.cardinality()
        M = dict([(i,j),0] for i in range(l) for j in range(l))
        if labeling == 'source':
            for i in range(l):
                perm = [ self.poset().unwrap(k) for k in L[i] ]
                for j in range(n-1):
                    p = getattr(L[i], action)(j+1)
                    M[(L.index(p),i)] += x[perm[j]-1]
        else:
            for i in range(l):
                for j in range(n-1):
                    p = getattr(L[i], action)(j+1)
                    M[(L.index(p),i)] += x[j]
        for i in range(l):
            M[(i,i)] += -sum(M[(j,i)] for j in range(l))
        return matrix(l,l,lambda x,y : M[(x,y)])

    def _element_constructor_(self, lst, check=True):
        r"""
        Constructor for elements of this class.

        TESTS::

            sage: P = Poset(([1,2,3,4], [[1,2],[1,4],[2,3]]))
            sage: L = P.linear_extensions()
            sage: x = L._element_constructor_([1,2,4,3]); x
            [1, 2, 4, 3]
            sage: x.parent() is L
            True

            sage: L._element_constructor_([4,3,2,1])
            Traceback (most recent call last):
            ...
            ValueError: [4, 3, 2, 1] is not a linear extension of Finite poset containing 4 elements
            sage: L._element_constructor_([4,3,2,1],check=False)
            [4, 3, 2, 1]
        """
        if not isinstance(lst, (list, tuple)):
            raise TypeError, "Input should be a list or tuple."
        lst = map(self._poset, lst)
        if self._is_facade:
            return lst
        else:
            return self.element_class(self, lst, check)

    Element = LinearExtensionOfPoset
