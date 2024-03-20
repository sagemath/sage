# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Blob Algebras

AUTHORS:

- Travis Scrimshaw (2020-05-16): Initial version
"""

# ****************************************************************************
#       Copyright (C) 2020 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element, get_coercion_model
from sage.structure.richcmp import richcmp
#from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.cachefunc import cached_method
from sage.combinat.subset import powerset
from sage.arith.misc import binomial
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.algebras import Algebras
from sage.combinat.diagram_algebras import (TemperleyLiebDiagrams, diagram_latex,
                                            TL_diagram_ascii_art)
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.dyck_word import DyckWords

#@add_metaclass(InheritComparisonClasscallMetaclass)


class BlobDiagram(Element):
    r"""
    A blob diagram.

    A blob diagram consists of a perfect matching of the set
    `\{1, \ldots, n\} \sqcup \{-1, \ldots, -n\}` such that the result
    is a noncrossing matching (a :class:`Temperley-Lieb diagram
    <sage.combinat.diagram_algebras.TemperleyLiebDiagram>`), divided
    into two sets of pairs: one for the pairs with blobs and one for
    those without. The blobed pairs must either be either the leftmost
    propagating strand or to the left of it and not nested.
    """

    def __init__(self, parent, marked, unmarked):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.combinat.blob_algebra import BlobDiagrams
            sage: BD4 = BlobDiagrams(4)
            sage: B = BD4([[1,-3]], [[2,-4], [3,4], [-1,-2]])
            sage: TestSuite(B).run()
        """
        Element.__init__(self, parent)
        self.marked = tuple(sorted([tuple(sorted(pair)) for pair in marked]))
        self.unmarked = tuple(sorted([tuple(sorted(pair)) for pair in unmarked]))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.blob_algebra import BlobDiagrams
            sage: BD4 = BlobDiagrams(4)
            sage: BD4([[1,-3]], [[2,-4], [3,4], [-1,-2]])
            ({{-3, 1}}, {{-4, 2}, {-2, -1}, {3, 4}})
        """
        return '({{{}}}, {{{}}})'.format(', '.join('{' + repr(X)[1:-1] + '}'
                                                   for X in self.marked),
                                         ', '.join('{' + repr(X)[1:-1] + '}'
                                                   for X in self.unmarked))

    def __hash__(self):
        r"""
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.combinat.blob_algebra import BlobDiagrams
            sage: BD4 = BlobDiagrams(4)
            sage: B = BD4([[1,-3]], [[2,-4], [3,4], [-1,-2]])
            sage: hash(B) in [hash(D) for D in BD4]
            True
            sage: len(set([hash(D) for D in BD4])) == len(BD4)
            True
        """
        return hash((self.marked, self.unmarked))

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` to ``other`` with operation ``op``.

        EXAMPLES::

            sage: from sage.combinat.blob_algebra import BlobDiagrams
            sage: BD4 = BlobDiagrams(4)
            sage: B = BD4([[1,-3]], [[2,-4], [3,4], [-1,-2]])
            sage: any(B == D for D in BD4)
            True
            sage: B2 = BD4([], [[1,-3], [2,-4], [3,4], [-1,-2]])
            sage: B == B2
            False
            sage: B != B2
            True
            sage: sorted(BlobDiagrams(3))
            [({}, {{-3, -2}, {-1, 1}, {2, 3}}),
             ({}, {{-3, -2}, {-1, 3}, {1, 2}}),
             ({}, {{-3, 1}, {-2, -1}, {2, 3}}),
             ({}, {{-3, 3}, {-2, -1}, {1, 2}}),
             ({}, {{-3, 3}, {-2, 2}, {-1, 1}}),
             ({{-3, 1}}, {{-2, -1}, {2, 3}}),
             ({{-3, 3}}, {{-2, -1}, {1, 2}}),
             ({{-2, -1}}, {{-3, 1}, {2, 3}}),
             ({{-2, -1}}, {{-3, 3}, {1, 2}}),
             ({{-1, 1}}, {{-3, -2}, {2, 3}}),
             ({{-1, 1}}, {{-3, 3}, {-2, 2}}),
             ({{-1, 3}}, {{-3, -2}, {1, 2}}),
             ({{1, 2}}, {{-3, -2}, {-1, 3}}),
             ({{1, 2}}, {{-3, 3}, {-2, -1}}),
             ({{-3, 1}, {-2, -1}}, {{2, 3}}),
             ({{-3, 3}, {-2, -1}}, {{1, 2}}),
             ({{-3, 3}, {1, 2}}, {{-2, -1}}),
             ({{-2, -1}, {1, 2}}, {{-3, 3}}),
             ({{-1, 3}, {1, 2}}, {{-3, -2}}),
             ({{-3, 3}, {-2, -1}, {1, 2}}, {})]
        """
        return richcmp((len(self.marked), self.marked, self.unmarked),
                       (len(other.marked), other.marked, other.unmarked),
                       op)

    def temperley_lieb_diagram(self):
        r"""
        Return the Temperley-Lieb diagram corresponding to ``self``.

        EXAMPLES::

            sage: from sage.combinat.blob_algebra import BlobDiagrams
            sage: BD4 = BlobDiagrams(4)
            sage: B = BD4([[1,-3]], [[2,-4], [3,4], [-1,-2]])
            sage: B.temperley_lieb_diagram()
            {{-4, 2}, {-3, 1}, {-2, -1}, {3, 4}}
        """
        return self.parent()._TL_diagrams(self.marked + self.unmarked)


class BlobDiagrams(Parent, UniqueRepresentation):
    r"""
    The set of all blob diagrams.
    """

    def __init__(self, n):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.combinat.blob_algebra import BlobDiagrams
            sage: BD4 = BlobDiagrams(4)
            sage: TestSuite(BD4).run()
        """
        self._n = n
        self._TL_diagrams = TemperleyLiebDiagrams(n)
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.blob_algebra import BlobDiagrams
            sage: BlobDiagrams(4)
            Blob diagrams of order 4
        """
        return "Blob diagrams of order {}".format(self._n)

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: from sage.combinat.blob_algebra import BlobDiagrams
            sage: BD4 = BlobDiagrams(4)
            sage: BD4.cardinality()
            70
        """
        return binomial(2*self._n, self._n)

    def order(self):
        r"""
        Return the order of ``self``.

        EXAMPLES::

            sage: from sage.combinat.blob_algebra import BlobDiagrams
            sage: BD4 = BlobDiagrams(4)
            sage: BD4.order()
            4
        """
        return self._n

    @cached_method
    def base_set(self):
        r"""
        Return the base set of ``self``.

        EXAMPLES::

            sage: from sage.combinat.blob_algebra import BlobDiagrams
            sage: BD4 = BlobDiagrams(4)
            sage: sorted(BD4.base_set())
            [-4, -3, -2, -1, 1, 2, 3, 4]
        """
        return frozenset(range(1,self._n+1)).union(range(-self._n,0))

    def _element_constructor_(self, marked, unmarked=None):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: from sage.combinat.blob_algebra import BlobDiagrams
            sage: BD4 = BlobDiagrams(4)
            sage: BD4([[1,-3]], [[-1,-2], [2,3], [-4,4]])
            ({{-3, 1}}, {{-4, 4}, {-2, -1}, {2, 3}})
            sage: BD4([[(1,-3)], ([-1,-2], (2,3), [-4,4])])
            ({{-3, 1}}, {{-4, 4}, {-2, -1}, {2, 3}})
        """
        if unmarked is None:
            marked, unmarked = marked
        ret = self.element_class(self, marked, unmarked)
        if ret not in self:
            raise ValueError("not a blob diagram of order {}".format(self._n))
        return ret

    def __contains__(self, X):
        r"""
        Check if ``X`` is contained in ``self``.

        EXAMPLES::

            sage: from sage.combinat.blob_algebra import BlobDiagrams
            sage: BD4 = BlobDiagrams(4)
            sage: BD4([[1,-3], [-1,-2]], [[2,-4], [3,4]])  # indirect doctest
            ({{-3, 1}, {-2, -1}}, {{-4, 2}, {3, 4}})
            sage: BD4([[1,4], [-1,-2], [-3,-4]], [[2,3]])  # indirect doctest
            ({{-4, -3}, {-2, -1}, {1, 4}}, {{2, 3}})

            sage: BD4([[1,-2], [-1,-3]], [[2,-4], [3,4]])  # crossing strands
            Traceback (most recent call last):
            ...
            ValueError: not a blob diagram of order 4
            sage: BD4([[1,-4], [-1,-2]], [[2,-3], [3,4]])  # crossing strands
            Traceback (most recent call last):
            ...
            ValueError: not a blob diagram of order 4
            sage: BD4([[1,-2], [-1,-3]], [[3,-4], [2,4]])  # crossing strands
            Traceback (most recent call last):
            ...
            ValueError: not a blob diagram of order 4
            sage: BD4([[1,-3], [-1,-2], [3,4]], [[2,-4]])  # trapped blob cup
            Traceback (most recent call last):
            ...
            ValueError: not a blob diagram of order 4
            sage: BD4([[-1,3], [1,2], [-3,-4]], [[-2,4]])  # trapped blob cap
            Traceback (most recent call last):
            ...
            ValueError: not a blob diagram of order 4
            sage: BD4([[1,4], [-1,-2], [-3,-4], [2,3]], [])  # nested blob cup
            Traceback (most recent call last):
            ...
            ValueError: not a blob diagram of order 4
            sage: BD4([[-1,-4], [1,2], [3,4], [-2,-3]], [])  # nested blob cap
            Traceback (most recent call last):
            ...
            ValueError: not a blob diagram of order 4
            sage: BD4([[3,-3]], [[1,-1],[2,-2],[4,-4]])  # trapped propagating line
            Traceback (most recent call last):
            ...
            ValueError: not a blob diagram of order 4
        """
        if not isinstance(X, BlobDiagram):
            return False
        # Check that it is a Temperley-Lieb diagram
        TL = X.marked + X.unmarked  # the TL diagram
        if TL not in self._TL_diagrams:
            return False
        # Check left escaping
        for x, y in X.marked:
            if x > 0:  # Must be a cup
                for P in TL:
                    if P[1] < 0:  # P is a cap
                        continue
                    if P[1] < x:
                        if P[0] < 0:  # A propagating line to the left
                            return False
                    else:  # Note that P[1] != x
                        if 0 < P[0] < x:  # A nesting line
                            return False
            elif y < 0:  # Must be a cap
                for P in TL:
                    if P[0] > 0:  # P is a cup
                        continue
                    if P[0] > y:
                        if P[1] > 0:  # A propagating line to the left
                            return False
                    else:  # Note that P[0] != y
                        if 0 > P[1] > y:  # A nesting line
                            return False
            else:  # Must be a propagating line
                if any(P[0] < 0 and P[1] > 0 and P[1] < y for P in TL):
                    return False
        return True

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: from sage.combinat.blob_algebra import BlobDiagrams
            sage: BD3 = BlobDiagrams(3)
            sage: sorted(BD3)
            [({}, {{-3, -2}, {-1, 1}, {2, 3}}),
             ({}, {{-3, -2}, {-1, 3}, {1, 2}}),
             ({}, {{-3, 1}, {-2, -1}, {2, 3}}),
             ({}, {{-3, 3}, {-2, -1}, {1, 2}}),
             ({}, {{-3, 3}, {-2, 2}, {-1, 1}}),
             ({{-3, 1}}, {{-2, -1}, {2, 3}}),
             ({{-3, 3}}, {{-2, -1}, {1, 2}}),
             ({{-2, -1}}, {{-3, 1}, {2, 3}}),
             ({{-2, -1}}, {{-3, 3}, {1, 2}}),
             ({{-1, 1}}, {{-3, -2}, {2, 3}}),
             ({{-1, 1}}, {{-3, 3}, {-2, 2}}),
             ({{-1, 3}}, {{-3, -2}, {1, 2}}),
             ({{1, 2}}, {{-3, -2}, {-1, 3}}),
             ({{1, 2}}, {{-3, 3}, {-2, -1}}),
             ({{-3, 1}, {-2, -1}}, {{2, 3}}),
             ({{-3, 3}, {-2, -1}}, {{1, 2}}),
             ({{-3, 3}, {1, 2}}, {{-2, -1}}),
             ({{-2, -1}, {1, 2}}, {{-3, 3}}),
             ({{-1, 3}, {1, 2}}, {{-3, -2}}),
             ({{-3, 3}, {-2, -1}, {1, 2}}, {})]
        """
        for D in DyckWords(self._n):
            markable = set()
            unmarked = []
            unpaired = []
            # Determine the pairing and which pairings are markable
            for i,d in enumerate(D):
                if i >= self._n:
                    i = -2*self._n + i
                else:
                    i += 1
                if d == 1:
                    unpaired.append(i)
                else: # d == 0
                    m = unpaired.pop()
                    if not unpaired:
                        markable.add((m, i))
                    else:
                        unmarked.append((m, i))
            for X in powerset(markable):
                yield self.element_class(self, X, unmarked + list(markable.difference(X)))

    Element = BlobDiagram


class BlobAlgebra(CombinatorialFreeModule):
    r"""
    The blob algebra.

    The *blob algebra* (also known as the Temperley-Lieb algebra of type `B`
    in [ILZ2018]_, but is a quotient of the Temperley-Lieb algebra of type `B`
    defined in [Graham1985]_) is a diagram-type algebra introduced in
    [MS1994]_ whose basis consists of :class:`Temperley-Lieb diagrams
    <sage.combinat.diagram_algebras.TemperleyLiebDiagram>`, noncrossing
    perfect matchings, that may contain blobs on strands that can be
    deformed so that the blob touches the left side (which we can think of
    as a frozen pole).

    The form we give here has 3 parameters, the natural one from the
    :class:`Temperley-Lieb algebra <sage.combinat.diagram_algebras.TemperleyLiebAlgebra>`,
    one for the idempotent relation, and one for a loop with a blob.

    INPUT:

    - ``k`` -- the order
    - ``q1`` -- the loop parameter
    - ``q2`` -- the idempotent parameter
    - ``q3`` -- the blob loop parameter

    EXAMPLES::

        sage: R.<q,r,s> = ZZ[]
        sage: B4 = algebras.Blob(4, q, r, s)
        sage: B = sorted(B4.basis())
        sage: B[14]
        B({{-4, -3}}, {{-2, -1}, {1, 2}, {3, 4}})
        sage: B[40]
        B({{3, 4}}, {{-4, -3}, {-2, -1}, {1, 2}})
        sage: B[14] * B[40]
        q*r*s*B({}, {{-4, -3}, {-2, -1}, {1, 2}, {3, 4}})

    REFERENCES:

    - [MS1994]_
    - [ILZ2018]_
    """
    @staticmethod
    def __classcall_private__(cls, k, q1, q2, q3, base_ring=None, prefix='B'):
        r"""
        Normalize input to ensure a unique representation.

        TESTS::

            sage: R.<q,r,s> = ZZ[]
            sage: B3 = algebras.Blob(3, q, r, s)
            sage: Bp = algebras.Blob(3, q, r, s, R, prefix='B')
            sage: B3 is Bp
            True
        """
        if base_ring is None:
            base_ring = get_coercion_model().common_parent(q1, q2, q3)
        q1 = base_ring(q1)
        q2 = base_ring(q2)
        q3 = base_ring(q3)
        return super().__classcall__(cls, k, q1, q2, q3, base_ring, prefix)

    def __init__(self, k, q1, q2, q3, base_ring, prefix):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R.<q,r,s> = ZZ[]
            sage: B4 = algebras.Blob(4, q, r, s)
            sage: TestSuite(B4).run()

            sage: B3 = algebras.Blob(3, q, r, s)
            sage: B = list(B3.basis())
            sage: TestSuite(B3).run(elements=B)  # long time
        """
        self._q1 = q1
        self._q2 = q2
        self._q3 = q3
        diagrams = BlobDiagrams(k)
        cat = Algebras(base_ring.category()).FiniteDimensional().WithBasis()
        CombinatorialFreeModule.__init__(self, base_ring, diagrams, category=cat,
                                         prefix=prefix, bracket=False)

    def _ascii_art_term(self, diagram):
        r"""
        Return an ascii art representation of ``diagram``.

        EXAMPLES::

            sage: R.<q,r,s> = ZZ[]
            sage: B2 = algebras.Blob(2, q, r, s)
            sage: x = B2.an_element()
            sage: ascii_art(x)  # indirect doctest
               o o      o o      o o
            2* `-` + 3* `-` + 2* `0`
               .-.      .0.      .-.
               o o      o o      o o
        """
        return TL_diagram_ascii_art(diagram.marked+diagram.unmarked, use_unicode=False,
                                    blobs=diagram.marked)

    def _unicode_art_term(self, diagram):
        r"""
        Return a unicode art representation of ``diagram``.

        EXAMPLES::

            sage: R.<q,r,s> = ZZ[]
            sage: B2 = algebras.Blob(2, q, r, s)
            sage: x = B2.an_element()
            sage: unicode_art(x)  # indirect doctest
               ⚬ ⚬      ⚬ ⚬      ⚬ ⚬
            2* ╰─╯ + 3* ╰─╯ + 2* ╰●╯
               ╭─╮      ╭●╮      ╭─╮
               ⚬ ⚬      ⚬ ⚬      ⚬ ⚬
        """
        return TL_diagram_ascii_art(diagram.marked+diagram.unmarked, use_unicode=True,
                                    blobs=diagram.marked)

    def _latex_term(self, diagram):
        r"""
        Return a latex representation of ``diagram``.

        EXAMPLES::

            sage: R.<q,r,s> = ZZ[]
            sage: B2 = algebras.Blob(2, q, r, s)
            sage: latex(B2.an_element())  # indirect doctest
            2 \begin{tikzpicture}[scale = 0.5,thick, baseline={(0,-1ex/2)}]
            \tikzstyle{vertex} = [shape = circle, minimum size = 7pt, inner sep = 1pt]
            \node[vertex] (G--2) at (1.5, -1) [shape = circle, draw] {};
            \node[vertex] (G--1) at (0.0, -1) [shape = circle, draw] {};
            \node[vertex] (G-1) at (0.0, 1) [shape = circle, draw] {};
            \node[vertex] (G-2) at (1.5, 1) [shape = circle, draw] {};
            \draw[] (G--2) .. controls +(-0.5, 0.5) and +(0.5, 0.5) .. (G--1);
            \draw[] (G-1) .. controls +(0.5, -0.5) and +(-0.5, -0.5) .. (G-2);
            \end{tikzpicture}
             + 3 \begin{tikzpicture}[scale = 0.5,thick, baseline={(0,-1ex/2)}]
            \tikzstyle{vertex} = [shape = circle, minimum size = 7pt, inner sep = 1pt]
            \node[vertex] (G--2) at (1.5, -1) [shape = circle, draw] {};
            \node[vertex] (G--1) at (0.0, -1) [shape = circle, draw] {};
            \node[vertex] (G-1) at (0.0, 1) [shape = circle, draw] {};
            \node[vertex] (G-2) at (1.5, 1) [shape = circle, draw] {};
            \draw[blue,very thick] (G--2) .. controls +(-0.5, 0.5) and +(0.5, 0.5) .. node[midway,circle,fill,scale=0.6] {} (G--1);
            \draw[] (G-1) .. controls +(0.5, -0.5) and +(-0.5, -0.5) .. (G-2);
            \end{tikzpicture}
             + 2 \begin{tikzpicture}[scale = 0.5,thick, baseline={(0,-1ex/2)}]
            \tikzstyle{vertex} = [shape = circle, minimum size = 7pt, inner sep = 1pt]
            \node[vertex] (G-1) at (0.0, 1) [shape = circle, draw] {};
            \node[vertex] (G-2) at (1.5, 1) [shape = circle, draw] {};
            \node[vertex] (G--2) at (1.5, -1) [shape = circle, draw] {};
            \node[vertex] (G--1) at (0.0, -1) [shape = circle, draw] {};
            \draw[blue,very thick] (G-1) .. controls +(0.5, -0.5) and +(-0.5, -0.5) .. node[midway,circle,fill,scale=0.6] {} (G-2);
            \draw[] (G--2) .. controls +(-0.5, 0.5) and +(0.5, 0.5) .. (G--1);
            \end{tikzpicture}
        """
        def edge_options(P):
            if P[1] < P[0]:
                P = [P[1], P[0]]
            if tuple(P) in diagram.marked:
                return 'blue,very thick'
            return ''

        def edge_additions(P):
            if P[1] < P[0]:
                P = [P[1], P[0]]
            if tuple(P) in diagram.marked:
                return 'node[midway,circle,fill,scale=0.6] {} '
            return ''
        return diagram_latex(diagram.marked+diagram.unmarked,
                             edge_options=edge_options,
                             edge_additions=edge_additions)

    def order(self):
        r"""
        Return the order of ``self``.

        The order of a partition algebra is defined as half of the number
        of nodes in the diagrams.

        EXAMPLES::

            sage: R.<q,r,s> = ZZ[]
            sage: B4 = algebras.Blob(4, q, r, s)
            sage: B4.order()
            4
        """
        return self._indices.order()

    @cached_method
    def one_basis(self):
        r"""
        Return the index of the basis element `1`.

        EXAMPLES::

            sage: R.<q,r,s> = ZZ[]
            sage: B4 = algebras.Blob(4, q, r, s)
            sage: B4.one_basis()
            ({}, {{-4, 4}, {-3, 3}, {-2, 2}, {-1, 1}})
        """
        B = self._indices
        return B.element_class(B, [], [[i, -i] for i in range(1, self.order()+1)])

    def product_on_basis(self, top, bot):
        r"""
        Return the product of the basis elements indexed by ``top``
        and ``bot``.

        EXAMPLES::

            sage: R.<q,r,s> = ZZ[]
            sage: B4 = algebras.Blob(4, q, r, s)
            sage: B = B4.basis()
            sage: BD = sorted(B.keys())
            sage: BD[14]
            ({{-4, -3}}, {{-2, -1}, {1, 2}, {3, 4}})
            sage: BD[40]
            ({{3, 4}}, {{-4, -3}, {-2, -1}, {1, 2}})
            sage: B4.product_on_basis(BD[14], BD[40])
            q*r*s*B({}, {{-4, -3}, {-2, -1}, {1, 2}, {3, 4}})
            sage: all(len((x*y).support()) == 1 for x in B for y in B)
            True
        """
        ret_lists = [[], []]
        coeff = self.base_ring().one()
        top_marked = set(top.marked)
        top_unmarked = set(top.unmarked)
        bot_marked = set(bot.marked)
        bot_unmarked = set(bot.unmarked)

        for top_set, is_unmarked in [(top_marked, 0), (top_unmarked, 1)]:
            while top_set:
                # We are starting a new strand
                cur, stop = top_set.pop()  # note that cur < stop
                unmarked = is_unmarked
                #print(top_set, unmarked, cur, stop)
                if cur > 0:  # Both are anchored to the top
                    ret_lists[unmarked].append((cur, stop))
                    continue
                anchored = bool(stop > 0)  # Possibly only stop is anchored

                # Follow the path from cur until we either reach stop or
                #   we break out of the loop because both ends are anchored
                while anchored or cur != stop:
                    #print(anchored, unmarked, cur, stop)
                    cur = -cur  # Move cur to the bottom diagram
                    for X in bot_marked:
                        if cur in X:
                            if unmarked:
                                unmarked = 0
                            else:
                                coeff *= self._q2
                            prev = cur
                            cur = X[1-X.index(prev)]
                            bot_marked.remove(X)
                            break
                    for X in bot_unmarked:
                        if cur in X:
                            prev = cur
                            cur = X[1-X.index(prev)]
                            bot_unmarked.remove(X)
                            break
                    if cur < 0:  # cur is anchored at the bottom
                        if anchored:
                            ret_lists[unmarked].append((stop, cur))
                            break
                        else:
                            anchored = True
                            stop, cur = cur, stop # stop is now anchored to the bottom
                            continue
                    cur = -cur  # bring cur back to the top diagram
                    for X in top_marked:
                        if cur in X:
                            if unmarked:
                                unmarked = 0
                            else:
                                coeff *= self._q2
                            prev = cur
                            cur = X[1-X.index(prev)]
                            top_marked.remove(X)
                            break
                    for X in top_unmarked:
                        if cur in X:
                            prev = cur
                            cur = X[1-X.index(prev)]
                            top_unmarked.remove(X)
                            break
                    if cur > 0:  # cur is anchored at the top
                        if anchored:
                            ret_lists[unmarked].append((stop, cur))
                            break
                        else:
                            anchored = True
                            stop, cur = cur, stop # stop is now anchored to the top
                if cur == stop:  # We have found a (marked) loop
                    if unmarked:
                        coeff *= self._q1
                    else:
                        coeff *= self._q3
        # Everything remaining in the bottom sets are just anchored
        #   at the bottom, (i.e., are of the form {-i, -j}).
        ret_lists[0].extend(bot_marked)
        ret_lists[1].extend(bot_unmarked)

        if coeff == 0:
            return self.zero()
        diagram = self._indices.element_class(self._indices, ret_lists[0], ret_lists[1])
        return self._from_dict({diagram: coeff}, remove_zeros=False)
