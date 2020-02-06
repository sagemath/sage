"""
Linear Extensions of Directed Acyclic Graphs

A linear extension of a directed acyclic graph is a total (linear) ordering on
the vertices that is compatible with the graph in the following sense:
if there is a path from x to y in the graph, the x appears before y in the
linear extension.

The algorithm implemented in this module is from "Generating Linear Extensions
Fast" by Preusse and Ruskey, which can be found at
http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.52.3057 .  The algorithm
generates the extensions in constant amortized time (CAT) -- a constant amount
of time per extension generated, or linear in the number of extensions
generated.

EXAMPLES:

Here we generate the 5 linear extensions of the following directed
acyclic graph::

    sage: from sage.graphs.linearextensions import LinearExtensions
    sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
    sage: D.is_directed_acyclic()
    True
    sage: LinearExtensions(D).list()
    doctest:...: DeprecationWarning: LinearExtensions is deprecated; use FinitePoset.linear_extensions or DiGraph.topological_sort_generator instead
    See https://trac.sagemath.org/25864 for details.
    [[0, 1, 2, 3, 4],
     [0, 1, 2, 4, 3],
     [0, 2, 1, 3, 4],
     [0, 2, 1, 4, 3],
     [0, 2, 4, 1, 3]]

Notice how all of the total orders are compatible with the ordering
induced from the graph.

We can also get at the linear extensions directly from the graph.  From
the graph, the linear extensions are known as topological sorts ::

    sage: list(D.topological_sort_generator())
    [[0, 1, 2, 3, 4],
     [0, 2, 1, 3, 4],
     [0, 2, 1, 4, 3],
     [0, 2, 4, 1, 3],
     [0, 1, 2, 4, 3]]

"""
#*****************************************************************************
#      Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.combinat import CombinatorialClass

class LinearExtensionsOld(CombinatorialClass):
    def __init__(self, dag):
        r"""
        Creates an object representing the class of all linear extensions
        of the directed acyclic graph \code{dag}.

        EXAMPLES::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: l = LinearExtensions(D)
            doctest:...: DeprecationWarning: LinearExtensions is deprecated; use FinitePoset.linear_extensions or DiGraph.topological_sort_generator instead
            See https://trac.sagemath.org/25864 for details.

            sage: l == loads(dumps(l))
            True

        TESTS::

            sage: list(LinearExtensions(DiGraph({ })))
            [[]]

            sage: LinearExtensions(DiGraph({ 0:[1], 1:[0] }))
            Traceback (most recent call last):
            ...
            ValueError: The graph is not directed acyclic

        """
        from sage.combinat.posets.posets import Poset
        self._dag = Poset(dag) # this returns a copy

    def _repr_(self):
        """
        TESTS::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: LinearExtensions(D)
            doctest:...: DeprecationWarning: LinearExtensions is deprecated; use FinitePoset.linear_extensions or DiGraph.topological_sort_generator instead
            See https://trac.sagemath.org/25864 for details.
            Linear extensions of Finite poset containing 5 elements

        """
        return "Linear extensions of %s"%self._dag

    def list(self):
        """
        Returns a list of the linear extensions of the directed acyclic graph.

        EXAMPLES::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: LinearExtensions(D).list()
            doctest:...: DeprecationWarning: LinearExtensions is deprecated; use FinitePoset.linear_extensions or DiGraph.topological_sort_generator instead
            See https://trac.sagemath.org/25864 for details.
            [[0, 1, 2, 3, 4],
             [0, 1, 2, 4, 3],
             [0, 2, 1, 3, 4],
             [0, 2, 1, 4, 3],
             [0, 2, 4, 1, 3]]

        TESTS::

            sage: D = DiGraph({ "a":["b","c"], "b":["d"], "c":["d","e"] })
            sage: LinearExtensions(D).list()
            [['a', 'b', 'c', 'd', 'e'],
             ['a', 'b', 'c', 'e', 'd'],
             ['a', 'c', 'b', 'd', 'e'],
             ['a', 'c', 'b', 'e', 'd'],
             ['a', 'c', 'e', 'b', 'd']]

            sage: D = DiGraph({ 4:[3,2], 3:[1], 2:[1,0] })
            sage: LinearExtensions(D).list()
            [[4, 2, 0, 3, 1],
             [4, 2, 3, 0, 1],
             [4, 2, 3, 1, 0],
             [4, 3, 2, 0, 1],
             [4, 3, 2, 1, 0]]

        """
        from sage.combinat.combinat_cython import linear_extension_iterator
        elts = list(self._dag)
        return sorted([[elts[i] for i in e] for e in linear_extension_iterator(self._dag._hasse_diagram)])

def LinearExtensions(dag):
    r"""
    ``LinearExtensions`` is deprecated; use
    :meth:`sage.combinat.posets.FinitePoset.linear_extensions` or :meth:`sage.graphs.digraph.DiGraph.topological_sort_generator` instead.

    EXAMPLES::

        sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
        sage: Poset(D).linear_extensions().list()
        [[0, 1, 2, 3, 4],
         [0, 2, 1, 3, 4],
         [0, 2, 1, 4, 3],
         [0, 2, 4, 1, 3],
         [0, 1, 2, 4, 3]]

        sage: D.topological_sort_generator().list()
        [[0, 1, 2, 3, 4],
         [0, 2, 1, 3, 4],
         [0, 2, 1, 4, 3],
         [0, 2, 4, 1, 3],
         [0, 1, 2, 4, 3]]

        sage: D = DiGraph({ "a":["b","c"], "b":["d"], "c":["d","e"] })
        sage: Poset(D).linear_extensions().list()
        [['a', 'b', 'c', 'd', 'e'],
         ['a', 'c', 'b', 'd', 'e'],
         ['a', 'c', 'b', 'e', 'd'],
         ['a', 'c', 'e', 'b', 'd'],
         ['a', 'b', 'c', 'e', 'd']]

        sage: D.topological_sort_generator().list()
        [['a', 'b', 'c', 'd', 'e'],
         ['a', 'c', 'b', 'd', 'e'],
         ['a', 'c', 'b', 'e', 'd'],
         ['a', 'c', 'e', 'b', 'd'],
         ['a', 'b', 'c', 'e', 'd']]

    TESTS::

        sage: from sage.graphs.linearextensions import LinearExtensions
        sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
        sage: LinearExtensions(D).list()
        doctest:...: DeprecationWarning: LinearExtensions is deprecated; use FinitePoset.linear_extensions or DiGraph.topological_sort_generator instead
        See https://trac.sagemath.org/25864 for details.
        [[0, 1, 2, 3, 4],
         [0, 1, 2, 4, 3],
         [0, 2, 1, 3, 4],
         [0, 2, 1, 4, 3],
         [0, 2, 4, 1, 3]]

    """
    from sage.misc.superseded import deprecation
    deprecation(25864, "LinearExtensions is deprecated; use FinitePoset.linear_extensions or DiGraph.topological_sort_generator instead")

    return LinearExtensionsOld(dag)
