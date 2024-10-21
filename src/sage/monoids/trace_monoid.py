# sage_setup: distribution = sagemath-combinat
r"""
Module of trace monoids (free partially commutative monoids).

EXAMPLES:

We first create a trace monoid::

    sage: from sage.monoids.trace_monoid import TraceMonoid
    sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a'))); M
    Trace monoid on 3 generators ([a], [b], [c]) with independence relation {{a, c}}

Different elements can be equal because of the partially
commutative multiplication::

    sage: c * a * b == a * c * b
    True

We check that it is a monoid::

    sage: M in Monoids()
    True

REFERENCES:

- :wikipedia:`Trace_monoid`

- https://ncatlab.org/nlab/show/trace+monoid

AUTHORS:

- Pavlo Tokariev (2019-05-31): initial version
"""
# ****************************************************************************
#       Copyright (C) 2019 Pavlo Tokariev <pavlo.tokariev@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from itertools import repeat, chain, product

from sage.combinat.words.alphabet import Alphabet
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.monoids.free_monoid import FreeMonoid
from sage.monoids.monoid import Monoid_class
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.infinity import infinity
from sage.structure.element import MonoidElement
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.unique_representation import UniqueRepresentation

lazy_import('sage.graphs.digraph', 'DiGraph')
lazy_import('sage.graphs.graph', 'Graph')


class TraceMonoidElement(ElementWrapper, MonoidElement):
    r"""
    Element of a trace monoid, also known as a trace.

    Elements of trace monoid is actually a equivalence classes
    of related free monoid over some equivalence relation
    that in the case is presented as independence relation.

    .. RUBRIC:: Representative

    We transform each trace to its lexicographic form for the
    representative in the ambient free monoid. This is also used
    for comparisons.

    EXAMPLES::

        sage: from sage.monoids.trace_monoid import TraceMonoid
        sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
        sage: M.<a,b,c,d> = TraceMonoid(I=I)
        sage: x = b * a * d * a * c * b
        sage: x^3
        [b*a^2*d*b^2*c*a^2*d*b^2*c*a^2*d*b*c]
        sage: x^0
        1
        sage: x.lex_normal_form()
        b*a^2*d*b*c
        sage: x.foata_normal_form()
        (b, a*d, a, b*c)
    """
    def _repr_(self) -> str:
        """
        Textual representation of ``self``.

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: a * b
            [a*b]
            sage: b * a
            [b*a]
            sage: d * a
            [a*d]
        """
        if self == self.parent().one():
            return "1"
        return f"[{self.value}]"

    def _richcmp_(self, other, op) -> bool:
        r"""
        Compare two traces by their lexicographic normal forms.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: a^2 > a
            True
            sage: a*b < b*a
            True
            sage: a * c * b == a * b * c
            True
        """
        return self.value._richcmp_(other.value, op)

    def lex_normal_form(self):
        r"""
        Return the lexicographic normal form of ``self``.

        OUTPUT: a free monoid element

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: (a*b).lex_normal_form()
            a*b
            sage: (b*a).lex_normal_form()
            b*a
            sage: (d*a).lex_normal_form()
            a*d
        """
        return self.value

    def foata_normal_form(self) -> tuple:
        r"""
        Return the Foata normal form of ``self``.

        OUTPUT: tuple of free monoid elements

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b * a * d * a * c * b
            sage: x.foata_normal_form()
            (b, a*d, a, b*c)
        """
        return self.parent()._compute_foata_normal_form(self.value)

    def _mul_(self, other):
        r"""
        Concatenate one equivalence class with another.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: a * b * c == a * c * b
            True
        """
        return self.parent(self.value * other.value)

    def _flat_elements(self):
        r"""
        Return flatten list of generator numbers representing the trace.

        OUTPUT: list of generator indexes

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b * a^3 * d * a * c * b^2
            sage: x._flat_elements()
            [b, a, a, a, a, d, b, b, c]
        """
        return [g for g, times in self.value for _ in range(times)]

    @cached_method
    def dependence_graph(self):
        r"""
        Return dependence graph of the trace.

        It is a directed graph where all dependent (non-commutative)
        generators are connected by edges which
        direction depend on the generator position in the trace.

        OUTPUT: directed graph of generator indexes

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b * a * d * a * c * b
            sage: x.dependence_graph()                                                  # needs sage.graphs
            Digraph on 6 vertices
        """
        elements = self._flat_elements()
        independence = self.parent()._independence
        graph = {}

        for i, e in enumerate(elements):
            edges = [(v, i) for v in graph
                     if (e, elements[v]) not in independence]
            graph[i] = []
            for v1, v2 in edges:
                graph[v1].append(v2)

        return DiGraph(graph)

    @cached_method
    def hasse_diagram(self, algorithm='naive'):
        r"""
        Return Hasse diagram of the trace.

        Hasse diagram is a dependence graph without transitive edges.

        INPUT:

        - ``algorithm`` -- string (default: ``'naive'``); defines algorithm
          that will be used to compute Hasse diagram; there are two
          variants: ``'naive'`` and ``'min'``.

        OUTPUT: directed graph of generator indexes

        .. SEEALSO::

            :meth:`~sage.monoids.trace_monoid.TraceMonoidElement.naive_hasse_digram`,
            :meth:`~sage.monoids.trace_monoid.TraceMonoidElement.min_hasse_diagram`.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b * a * d * a * c * b
            sage: x.hasse_diagram()                                                     # needs sage.graphs
            Digraph on 6 vertices

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b * a * d * a * c * b
            sage: x.hasse_diagram(algorithm='naive') == x.hasse_diagram(algorithm='min')            # needs sage.graphs
            True
            sage: y = b * a^3 * d * a * c * b^2
            sage: y.hasse_diagram(algorithm='naive') == y.hasse_diagram(algorithm='min')            # needs sage.graphs
            True
        """
        if algorithm == "naive":
            return self.naive_hasse_diagram()
        elif algorithm == "min":
            return self.min_hasse_diagram()
        raise ValueError("`alg` option must be `naive` "
                         f"or `min`, got `{algorithm}`.")

    def min_hasse_diagram(self):
        r"""
        Return Hasse diagram of the trace.

        OUTPUT: directed graph of generator indexes

        .. SEEALSO::

            :meth:`~sage.monoids.trace_monoid.TraceMonoidElement.hasse_digram`,
            :meth:`~sage.monoids.trace_monoid.TraceMonoidElement.naive_hasse_diagram`.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b * a * d * a * c * b
            sage: x.min_hasse_diagram()                                                 # needs sage.graphs
            Digraph on 6 vertices
        """
        elements = self._flat_elements()
        elements.reverse()
        independence = self.parent()._independence
        reachable = {}
        mini = set()
        graph = DiGraph({})

        for i, x in enumerate(elements):
            reachable[i] = set()
            front = mini.copy()
            while front:
                used = set()
                for j in list(front):
                    y = elements[j]
                    if (x, y) not in independence:
                        graph.add_edge(i, j)
                        reachable[i].add(j)
                        reachable[i].update(reachable[j])
                        if j in mini:
                            mini.remove(j)
                        used.add(j)
                forbidden = set(chain.from_iterable(reachable[v] for v in used))
                front = {dest for _, dest in graph.outgoing_edges(front, labels=False)}
                front = front - forbidden

            mini.add(i)

        length = len(elements)
        graph.relabel(length - 1 - i for i in range(length))
        return graph

    def naive_hasse_diagram(self):
        r"""
        Return Hasse diagram of ``self``.

        ALGORITHM:

        In loop check for every two pair of edges if they
        have common vertex, remove their transitive edge.

        OUTPUT: directed graph of generator indexes

        .. SEEALSO::

            :meth:`~sage.monoids.trace_monoid.TraceMonoidElement.hasse_digram`,
            :meth:`~sage.monoids.trace_monoid.TraceMonoidElement.min_hasse_diagram`.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b * a * d * a * c * b
            sage: x.naive_hasse_diagram()                                               # needs sage.graphs
            Digraph on 6 vertices
        """
        d = self.dependence_graph()
        h = d.copy()

        d_edges = d.edges(sort=False)
        for e1 in d_edges:
            for e2 in d_edges:
                if e1[1] == e2[0]:
                    h.delete_edge((e1[0], e2[1]))

        return h

    def alphabet(self):
        r"""
        Return alphabet of ``self``.

        OUTPUT: a set of free monoid generators

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b*a*d*a*c*b
            sage: x.alphabet()
            {b, a, d, c}
        """
        return Alphabet([g for g, _ in self.value])

    def projection(self, letters):
        r"""
        Return a trace that formed from ``self`` by erasing ``letters``.

        INPUT:

        - ``letters`` -- set of generators; defines set of letters that will be
          used to filter the trace

        OUTPUT: a trace

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c,d> = FreeMonoid()
            sage: I = ((a,d), (d,a), (b,c), (c,b))
            sage: M.<ac,bc,cc,dc> = TraceMonoid(F, I=I)
            sage: x = M(b*a*d*a*c*b)
            sage: x.projection({a,b})
            [b*a^2*b]
            sage: x.projection({b,d,c})
            [b*d*b*c]
        """
        P = self.parent()
        base = P._free_monoid
        return P(base.prod(x for x in self._flat_elements() if x in letters))

    def multiplicative_order(self):
        r"""
        Return the multiplicative order of ``self``, which is `\infty`
        for any element not the identity.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: a.multiplicative_order()
            +Infinity
            sage: M.one().multiplicative_order()
            1
        """
        if self.value.is_one():
            return ZZ.one()
        return infinity


class TraceMonoid(UniqueRepresentation, Monoid_class):
    r"""
    Return a free partially commuting monoid (trace monoid) on `n` generators
    over independence relation `I`.

    We construct a trace monoid by specifying:

    - a free monoid and independence relation
    - or generator names and independence relation,
      FreeMonoid is constructed automatically then.

    INPUT:

    - ``M`` -- a free monoid

    - ``I`` -- commutation relation between generators
      (or their names if the ``names`` are given)

    - ``names`` -- names of generators

    EXAMPLES::

        sage: from sage.monoids.trace_monoid import TraceMonoid
        sage: F = TraceMonoid(names=('a', 'b', 'c'), I={('a','c'), ('c','a')}); F
        Trace monoid on 3 generators ([a], [b], [c]) with independence relation {{a, c}}
        sage: x = F.gens()
        sage: x[0]*x[1]**5 * (x[0]*x[2])
        [a*b^5*a*c]

        sage: from sage.monoids.trace_monoid import TraceMonoid
        sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
        sage: latex(M)
        \langle a, b, c \mid ac=ca \rangle

    TESTS::

        sage: from sage.monoids.trace_monoid import TraceMonoid
        sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
        sage: M.number_of_words(3) == len(M.words(3))                                   # needs sage.graphs
        True
    """
    Element = TraceMonoidElement

    @staticmethod
    def __classcall_private__(cls, M=None, I=frozenset(), names=None):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M1.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: M2.<a,b,c> = TraceMonoid(I=[('a','c')])
            sage: M3 = TraceMonoid(I=[{'a','c'}], names=('a', 'b', 'c'))
            sage: M1 is M2 and M2 is M3
            True
        """
        if not M:
            if names:
                M = FreeMonoid(names=names)
            else:
                raise ValueError("names must be provided")
        elif not names:
            names = [str(g) for g in M.gens()]
        names = tuple(names)

        rels = set()
        gen_from_str = {names[i]: gen for i, gen in enumerate(M.gens())}
        for x, y in I:
            try:
                if isinstance(x, str):
                    x = gen_from_str[x]
                x = M(x)
                if isinstance(y, str):
                    y = gen_from_str[y]
                y = M(y)
                if x == y:
                    raise ValueError
            except (TypeError, ValueError):
                raise ValueError("invalid relation defined")
            rels.add((x, y))
            rels.add((y, x))
        I = frozenset(rels)

        return super().__classcall__(cls, M, I, names)

    def __init__(self, M, I, names) -> None:
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: TestSuite(M).run()
        """
        self._free_monoid = M
        self._independence = I
        Monoid_class.__init__(self, names=names)

    def ngens(self):
        """
        Return the number of generators of ``self``.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: M.ngens()
            3
        """
        return self._free_monoid.ngens()

    def one(self):
        """
        Return the neutral element of ``self``.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: M.one()
            1
        """
        return self.element_class(self, self._free_monoid.one())

    def gen(self, i=0):
        """
        Return the `i`-th generator of the monoid.

        INPUT:

        - ``i`` -- integer (default: 0)

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: M.gen(1)
            [b]
            sage: M.gen(4)
            Traceback (most recent call last):
            ...
            IndexError: argument i (= 4) must be between 0 and 2
        """
        return self.element_class(self, self._free_monoid.gen(i))

    def cardinality(self):
        """
        Return the cardinality of ``self``, which is infinite except for
        the trivial monoid.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: M.cardinality()
            +Infinity
        """
        return self._free_monoid.cardinality()

    def _compute_dependence_stack(self, x):
        r"""
        Return generator stacks formed from trace
        subelements with respect to non-commutativity.

        OUTPUT: used generators and list of stacks as tuple

        ALGORITHM:

        Let `x` be a word of monoid; we scan `x` from right to left;
        when processing a letter `a` it is pushed on its stack and a
        marker is pushed on the stack of all the letters `b` ( `b \neq a` )
        which do not commute with `a`.

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c,d> = FreeMonoid()
            sage: I = (('ac','dc'), ('dc','ac'), ('bc','cc'), ('cc','bc'))
            sage: M.<ac,bc,cc,dc> = TraceMonoid(F, I=I)
            sage: x = b*a*d*a*c*b
            sage: M._compute_dependence_stack(x)
            ({a, b, c, d},
             {a: [False, False, True, True, False],
              b: [True, False, False, False, True],
              c: [True, False, False, False],
              d: [False, False, True, False]})
        """
        independence = self._independence
        generators_set = {e for e, _ in x}
        stacks = dict(sorted((g, []) for g in generators_set))
        for generator, times in reversed(list(x)):
            stacks[generator].extend(repeat(True, times))
            for other_gen in generators_set:
                if other_gen == generator:
                    continue
                if (generator, other_gen) not in independence:
                    stacks[other_gen].extend(repeat(False, times))
        return generators_set, stacks

    @cached_method
    def _compute_lex_normal_form(self, x):
        r"""
        Return lexicographic normal form of the free monoid
        element in free monoid terms.

        OUTPUT: trace monoid element

        ALGORITHM:

        Take among the letters being on the top of some stack that
        letter `a` being minimal with respect to the given lexicographic
        ordering. We pop a marker from each stack corresponding to a
        letter `b` ( `b \neq a` ) which does not commute with `a`. We repeat
        this loop until all stacks are empty.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c,d> = FreeMonoid()
            sage: I = ((a,d), (d,a), (b,c), (c,b))
            sage: M.<ac,bc,cc,dc> = TraceMonoid(F, I=I)
            sage: M._compute_lex_normal_form(c*a*c*b*a^2)
            c*a*b*c*a^2
        """
        if not x._element_list:
            return x
        generators_set, stacks = self._compute_dependence_stack(x)
        independence = self._independence

        elements = []
        while any(stacks.values()):
            for generator, g_stack in stacks.items():
                if g_stack and g_stack[-1]:
                    g_stack.pop()
                    elements.append(generator)
                    for other_gen in generators_set:
                        if (other_gen != generator
                                and (generator, other_gen) not in independence):
                            stacks[other_gen].pop()
                    break

        return prod(elements)

    @cached_method
    def _compute_foata_normal_form(self, x) -> tuple:
        r"""
        Return Foata normal form of the monoid element.

        OUTPUT: tuple of steps

        ALGORITHM:

        Within a loop we form the set using letters being
        on the top of stacks; arranging the letters in the lexicographic
        order yields a step of the Foata normal form;
        This loop is repeated until all stacks are empty.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c,d> = FreeMonoid()
            sage: I = ((a,d), (d,a), (b,c), (c,b))
            sage: M.<ac,bc,cc,dc> = TraceMonoid(F, I=I)
            sage: x = b*a*d*a*c*b
            sage: M._compute_foata_normal_form(x)
            (b, a*d, a, b*c)
            sage: y = b*a*a*d*b*a*b*c^2*a
            sage: M._compute_foata_normal_form(y)
            (b, a*d, a, b, a, b*c, c, a)
        """
        if not x._element_list:
            return ()

        generators_set, stacks = self._compute_dependence_stack(x)
        independence = self._independence

        steps = []
        while any(stacks.values()):
            step = []
            for generator, g_stack in stacks.items():
                if g_stack and g_stack[-1]:
                    g_stack.pop()
                    step.append(generator)

            for g in step:
                for other_gen in generators_set:
                    if other_gen != g and (g, other_gen) not in independence:
                        stacks[other_gen].pop()

            steps.append(step)

        return tuple(prod(step) for step in steps)

    def _element_constructor_(self, x):
        """
        Return ``x`` coerced into this trace monoid.

        One can create a free monoid element from the integer 1,
        free monoid elements of the same generators as internal one,
        and coerce everything that can coerce free monoid.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c,d> = FreeMonoid()
            sage: I = ((a,d), (d,a), (b,c), (c,b))
            sage: M.<ac,bc,cc,dc> = TraceMonoid(F, I=I)
            sage: x = b*a*d*a*c*b
            sage: M(x)
            [b*a^2*d*b*c]
        """
        x = self._compute_lex_normal_form(self._free_monoid(x))
        return self.element_class(self, x)

    @cached_method
    def independence(self):
        r"""
        Return independence relation over the monoid.

        OUTPUT: set of commuting generator pairs

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c> = FreeMonoid()
            sage: I = frozenset(((a,c), (c,a)))
            sage: M.<ac,bc,cc> = TraceMonoid(F, I=I)
            sage: M.independence() == frozenset([frozenset([a,c])])
            True
        """
        return frozenset(map(frozenset, self._independence))

    @cached_method
    def dependence(self):
        r"""
        Return dependence relation over the monoid.

        OUTPUT: set of non-commuting generator pairs

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: sorted(M.dependence())
            [(a, a), (a, b), (b, a), (b, b), (b, c), (c, b), (c, c)]
        """
        return frozenset(pair for pair in product(self._free_monoid.gens(), repeat=2)
                         if pair not in self._independence)

    @cached_method
    def dependence_graph(self):
        r"""
        Return graph of dependence relation.

        OUTPUT: dependence graph with generators as vertices

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c> = FreeMonoid()
            sage: M.<ai,bi,ci> = TraceMonoid(F, I=((a,c), (c,a)))
            sage: M.dependence_graph() == Graph({a:[a,b], b:[b], c:[c,b]})              # needs sage.graphs
            True
        """
        return Graph({frozenset((e1, e2)) if e1 != e2 else (e1, e2)
                      for e1, e2 in self.dependence()}, loops=True,
                     format='list_of_edges',
                     immutable=True)

    @cached_method
    def independence_graph(self):
        r"""
        Return the digraph of independence relations.

        OUTPUT: independence graph with generators as vertices

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c> = FreeMonoid()
            sage: M.<ai,bi,ci> = TraceMonoid(F, I=((a,c), (c,a)))
            sage: M.independence_graph() == Graph({a:[c], b:[], c:[]})                  # needs sage.graphs
            True
        """
        verts = list(self._free_monoid.gens())
        edges = list(map(list, self.independence()))
        return Graph([verts, edges], immutable=True)

    @cached_method
    def dependence_polynomial(self, t=None):
        r"""
        Return dependence polynomial.

        The polynomial is defined as follows: `\sum{i}{(-1)^i c_i t^i}`,
        where `c_i` equals to number of full subgraphs
        of size `i` in the independence graph.

        OUTPUT: a rational function in ``t`` with coefficients in the integer ring

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: M.dependence_polynomial()                                             # needs sage.graphs
            1/(2*t^2 - 4*t + 1)
        """
        if t is None:
            R = PolynomialRing(ZZ, 't')
            t = R.gen()
        clique_seq = self.independence_graph().clique_polynomial().coefficients()
        return ~sum((-1)**i * coeff * (t**i)
                    for i, coeff in enumerate(clique_seq))

    @cached_method
    def number_of_words(self, length):
        r"""
        Return number of unique words of defined length.

        INPUT:

        - ``length`` -- integer; defines size of words what number should be computed

        OUTPUT: words number as integer

        EXAMPLES:

        Get number of words of size 3 ::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: M.number_of_words(3)                                                  # needs sage.graphs
            48
        """
        psr = PowerSeriesRing(ZZ, default_prec=length + 1)
        return psr(self.dependence_polynomial()).coefficients()[length]

    @cached_method
    def words(self, length):
        r"""
        Return all lexicographic forms of defined length.

        INPUT:

        - ``length`` -- integer; defines size of words

        OUTPUT: set of traces of size ``length``

        EXAMPLES:

        All words of size 2::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: sorted(M.words(2))
            [[a^2], [a*b], [a*c], [a*d], [b*a], [b^2], [b*c],
             [b*d], [c*a], [c^2], [c*d], [d*b], [d*c], [d^2]]

        Get number of words of size 3::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: len(M.words(3))
            48

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','b'), ('b','a'), ('b', 'c'), ('c', 'b')))
            sage: for i in range(10):                                                   # needs sage.graphs
            ....:    assert len(M.words(i)) == M.number_of_words(i)
        """
        if length < 0:
            raise ValueError("bad length of words; expected zero or positive number")
        if length == 0:
            return frozenset([self.one()])
        if length == 1:
            return frozenset(self.gens())

        return frozenset([word * suffix for word in self.words(length - 1)
                          for suffix in self.gens()
                          if not ((list(word.value)[-1][0], suffix.value) in self._independence
                                  and list(word.value)[-1][0] > suffix.value)])

    def _sorted_independence(self) -> list:
        r"""
        Return independence relation over the monoid.

        OUTPUT: sorted list of sorted commuting generator pairs

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c> = FreeMonoid()
            sage: I = frozenset(((a,c), (c,a)))
            sage: M.<ac,bc,cc> = TraceMonoid(F, I=I)
            sage: M._sorted_independence()
            [[a, c]]
        """
        return sorted(sorted(x_y)
                      for x_y in self.independence())

    def _repr_(self) -> str:
        r"""
        Textual representation of trace monoids.

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I); M
            Trace monoid on 4 generators ([a], [b], [c], [d])
             with independence relation {{a, d}, {b, c}}
        """
        return ("Trace monoid on {!s} generators {!s} "
                "with independence relation {{{}}}").format(self.ngens(), self.gens(),
                                                            ", ".join(f"{{{x}, {y}}}"
                                                                      for (x, y) in self._sorted_independence()))

    def _latex_(self) -> str:
        r"""
        LaTeX representation of trace monoids.

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I); latex(M)
            \langle a, b, c, d \mid ad=da,bc=cb \rangle
        """
        return "\\langle {} \\mid {} \\rangle".format(
            repr(self._free_monoid.gens())[1:-1],
            ",".join(
                f"{v1!r}{v2!r}={v2!r}{v1!r}"
                for v1, v2 in self._sorted_independence()
            )
        )
