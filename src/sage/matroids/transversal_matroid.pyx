r"""
Transversal matroids

A transversal matroid arises from a ground set `E` and a collection `A` of sets
over the ground set. This can be modeled as a bipartite graph `B`, where the
vertices the left are ground set elements, the vertices on the right are the
sets, and edges represent containment. Then a set `X` from the ground set is
independent if and only if `X` has a matching in `B`.

Construction
============

To construct a transversal matroid, first import
:class:`~sage.matroids.transversal_matroid.TransversalMatroid` from
:mod:`sage.matroids.transversal_matroid`.
The input should be a set system, formatted as an iterable of iterables::

    sage: from sage.matroids.transversal_matroid import TransversalMatroid
    sage: sets = [[3, 4, 5, 6, 7, 8]] * 3
    sage: M = TransversalMatroid(sets); M
    Transversal matroid of rank 3 on 6 elements, with 3 sets
    sage: M.groundset()
    frozenset({3, 4, 5, 6, 7, 8})
    sage: M.is_isomorphic(matroids.Uniform(3, 6))
    True
    sage: M = TransversalMatroid([[0, 1], [1, 2, 3], [3, 4, 5]],
    ....: set_labels=['1', '2', '3'])
    sage: M.graph().vertices()
    ['1', '2', '3', 0, 1, 2, 3, 4, 5]

AUTHORS:

- Zachary Gershkoff (2017-08-07): initial version

REFERENCES:

- [Bonin2017]_
"""
#*****************************************************************************
#       Copyright (C) 2017 Zachary Gershkoff <zgersh2@lsu.edu>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function, absolute_import

from sage.matroids.matroid cimport Matroid
from sage.matroids.basis_exchange_matroid cimport BasisExchangeMatroid
from sage.matroids.minor_matroid import MinorMatroid
from sage.matroids.utilities import newlabel

from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph
from sage.graphs.bipartite_graph import BipartiteGraph

from cpython.object cimport Py_EQ, Py_NE
from copy import copy, deepcopy
from collections import Counter

import networkx as nx

cdef class TransversalMatroid(BasisExchangeMatroid):
    r"""
    The Transversal Matroid class.

    INPUT:

    - ``sets`` -- an iterable of iterables of elements
    - ``groundset`` -- (optional) an iterable containing names of ground set
      elements
    - ``set_labels`` -- (optional) a list of labels in 1-1 correspondence with
      the iterables in ``sets``
    - ``matching`` -- (optional) a dictionary specifying a maximal matching
      between elements and set labels

    OUTPUT:

    An instance of ``TransversalMatroid``. The sets specified in ``sets``
    define the matroid. If ``matching`` is not specified, the constructor
    will determine a matching to use for basis exchange.

    EXAMPLES::

        sage: from sage.matroids.transversal_matroid import TransversalMatroid
        sage: sets = [[0, 1, 2, 3]] * 3
        sage: M = TransversalMatroid(sets)
        sage: M.full_rank()
        3
        sage: M.bases_count()
        4
        sage: sum(1 for b in M.bases())
        4

    ::

        sage: from sage.matroids.transversal_matroid import TransversalMatroid
        sage: M = TransversalMatroid(sets=[['a', 'c']], groundset=['a', 'c', 'd'])
        sage: M.loops()
        frozenset({'d'})
        sage: M.full_rank()
        1
    """
    def __init__(self, sets, groundset=None, set_labels=None, matching=None):
        """
        See class definition for full documentation.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[0, 1, 2, 3], [1, 2], [1, 3, 4]]
            sage: set_labels = [5, 6, 7]
            sage: M = TransversalMatroid(sets, set_labels=set_labels)
            sage: TestSuite(M).run()

        TESTS::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: M = TransversalMatroid([[0, 1], [1, 2], [2, 3], [3, 4]])
            sage: TestSuite(M).run()

            sage: M = TransversalMatroid([[],[],[]], groundset=range(3)); M
            Transversal matroid of rank 0 on 3 elements, with 3 sets
            sage: sets = [[0, 1, 2, 3, 4], [4, 5]]
            sage: M = TransversalMatroid(sets, groundset=[0, 1, 2, 3])
            Traceback (most recent call last):
            ...
            ValueError: ground set and sets do not match
            sage: M = TransversalMatroid(sets, set_labels=[1, 2, 3])
            Traceback (most recent call last):
            ...
            ValueError: set labels do not match sets
            sage: M = TransversalMatroid(sets, set_labels=[1, 2])
            Traceback (most recent call last):
            ...
            ValueError: set labels cannot be element labels

            sage: sets = [['s1', 's2'], ['s1', 's3']]
            sage: M = TransversalMatroid(sets); M
            Transversal matroid of rank 2 on 3 elements, with 2 sets
            sage: M.sets()
            [['s1', 's2'], ['s1', 's3']]
            sage: M.set_labels()
            ['s0', 's4']

        Testing the matching parameter::

            sage: M = TransversalMatroid([range(5)] * 4, set_labels='abcd')
            sage: M.full_rank()
            4
            sage: M.rank([0,1,2])
            3

        ::

            sage: sets = [[0, 1, 2, 3], [1, 2], [1, 3, 4]]
            sage: set_labels = [5, 6, 7]
            sage: m = {0:5, 1:6, 3:7}
            sage: M = TransversalMatroid(sets, set_labels=set_labels, matching=m)
            sage: len(M.bases())
            9
        """
        contents = set([e for subset in sets for e in subset])
        if groundset is None:
            groundset = contents
        elif not contents.issubset(groundset):
            raise ValueError("ground set and sets do not match")
        groundset = tuple(groundset)

        # keep the original list as input so we don't lose order between minors etc.
        self._sets_input = [list(s) for s in sets]
        self._sets = Counter([frozenset(s) for s in self._sets_input])

        # This might be redundant with self._idx from BasisExchangeMatroid
        # However, BasisExchangeMatroid has not been called yet
        element_int_map = {e: i for i, e in enumerate(groundset)}
        int_element_map = {i: e for i, e in enumerate(groundset)}

        # we need a matching and a corresponding graph
        if set_labels:
            if len(set_labels) != len(sets):
                raise ValueError("set labels do not match sets")
            if not contents.isdisjoint(set_labels):
                raise ValueError("set labels cannot be element labels")
            if matching:
                matching_temp = matching
            self._set_labels_input = list(set_labels)

        else:
            if matching:
                raise ValueError("set labels must be provided if matching is provided")
            self._set_labels_input = []
            increment = 0
            for i in range(len(sets)): # make sure "sn" is not in ground set
                label = 's' + str(i + increment)
                while label in groundset:
                    increment += 1
                    label = 's' + str(i + increment)
                self._set_labels_input.append(label)
        # for internal use with integers representing elements
        set_labels = ['s' + str(i) for i in range(len(sets))]

        if not matching:
            B = BipartiteGraph()
            for e in groundset:
                B.add_vertex(element_int_map[e], left=True)
            for i, s in enumerate(sets):
                new_vertex = set_labels[i]
                for e in s:
                    B.add_edge(new_vertex, element_int_map[e])
            matching_temp = {}
            for u, v, _ in B.matching():
                if u in range(len(groundset)):
                    matching_temp[int_element_map[u]] = v
                else:
                    matching_temp[int_element_map[v]] = u

        self._set_labels = list(set_labels)

        # determine the basis from the matching
        basis = frozenset(matching_temp)

        # This creates self._groundset attribute, among other things
        # It takes the actual ground set labels, not the translated ones
        BasisExchangeMatroid.__init__(self, groundset, basis)

        # matching_temp uses actual ground set labels
        # self._matching will use the translated ones
        if matching:
            self._matching = {self._idx[e]: self._set_labels[self._set_labels_input.index(matching_temp[e])]
                              for e in matching_temp}
        else:
            self._matching = {self._idx[e]: matching_temp[e] for e in matching_temp}

        # Build a DiGraph for doing basis exchange
        self._D = nx.DiGraph()
        # Make sure we get isolated vertices, corresponding to loops
        self._D.add_nodes_from(self._idx.itervalues())
        # Also get isolated vertices corresponding to empty sets
        self._D.add_nodes_from(self._set_labels)

        # For sets in the matching, orient them as starting from the collections
        matching_reversed = [(v,k) for k,v in self._matching.iteritems()]
        self._D.add_edges_from(matching_reversed)

        other_edges = []
        for i, s in enumerate(sets):
            for e in s:
                if e not in matching_temp or matching_temp[e] != set_labels[i]:
                    other_edges.append((self._idx[e], set_labels[i]))
        self._D.add_edges_from(other_edges)

    cdef bint _is_exchange_pair(self, long x, long y) except -1:
        r"""
        Check for `M`-alternating path from `x` to `y`.
        """
        return nx.has_path(self._D, y, x)

    cdef int _exchange(self, long x, long y) except -1:
        r"""
        Replace ``self.basis() with ``self.basis() - x + y``.
        Internal method, does no checks.
        """
        # update the internal matching
        sh = nx.shortest_path(self._D, y, x)
        del self._matching[x]
        for i in range(0, len(sh)-1, 2):
            self._matching[sh[i]] = sh[i+1]

        # update the graph to reflect this new matching
        sh_edges = [(sh[i], sh[i + 1]) for i in xrange(len(sh) - 1)]
        sh_edges_r = [(sh[i + 1], sh[i]) for i in xrange(len(sh) - 1)]
        self._D.remove_edges_from(sh_edges)
        self._D.add_edges_from(sh_edges_r)

        BasisExchangeMatroid._exchange(self, x, y)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[0, 1, 2, 3]] * 3
            sage: M = TransversalMatroid(sets); M
            Transversal matroid of rank 3 on 4 elements, with 3 sets
        """
        return "Transversal matroid of rank {} on {} elements, with {} sets".format(
                            self.rank(), self.size(), len(self._sets_input))

    cpdef list sets(self):
        r"""
        Return the sets of ``self``.

        A transversal matroid can be viewed as a ground set with a collection
        from its powerset. This is represented as a bipartite graph, where
        an edge represents containment.

        OUTPUT:

        A list of lists that correspond to the sets of the transversal matroid.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[0, 1, 2, 3], [1, 2], [3, 4]]
            sage: set_labels = [5,6,7]
            sage: M = TransversalMatroid(sets, set_labels=set_labels)
            sage: sorted(M.sets()) == sorted(sets)
            True
        """
        return copy(self._sets_input)

    def __richcmp__(left, right, op):
        r"""
        Compare two matroids.

        We take a very restricted view on equality: the objects need to be of
        the exact same type (so no subclassing) and the internal data need to
        be the same. For transversal matroids, in particular, the presentation
        as a bipartite graph must be the same.

        .. WARNING::

            This method is linked to ``__hash__``. If you override one, you
            MUST override the other!

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [['a', 'b', 'c'], ['c', 'd'], ['a','d'], ['e']]
            sage: M = TransversalMatroid(sets)
            sage: M1 = TransversalMatroid(sets)
            sage: sets2 = [['a', 'b', 'c'], ['c', 'd'], ['a', 'd', 'e'], ['e']]
            sage: M2 = TransversalMatroid(sets2)
            sage: M1 == M2
            False
            sage: M1.equals(M2)
            True
            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[0, 1, 2, 3], [1, 2], [1, 3, 4]]
            sage: M = TransversalMatroid(sets, groundset=range(5), set_labels=[5, 6, 7])
            sage: M2 = TransversalMatroid(sets)
            sage: M == M2
            True

        TESTS::

            sage: M = TransversalMatroid([range(5)] * 4, set_labels='abcd', matching={0:'a', 1:'c'})
            sage: N = TransversalMatroid([range(5)] * 4, set_labels='abcd')
            sage: M == N
            False
        """
        if op not in [Py_EQ, Py_NE]:
            return NotImplemented
        if not isinstance(left, TransversalMatroid) or not isinstance(right, TransversalMatroid):
            return NotImplemented
        if left.__class__ != right.__class__:   # since we have some subclasses, an extra test
            return NotImplemented
        if op == Py_EQ:
            res = True
        if op == Py_NE:
            res = False
        # res gets inverted if matroids are deemed different.

        cdef TransversalMatroid rhs, lhs
        rhs = <TransversalMatroid> right
        lhs = <TransversalMatroid> left

        # If either matroid is not valid, deem them different
        if not (lhs.is_valid() and rhs.is_valid()):
            return not res
        if (lhs._groundset == rhs._groundset and
            Counter([frozenset(s) for s in lhs._sets_input]) == Counter([frozenset(s) for s in rhs._sets_input])):
            return res
        else:
            return not res

    def __hash__(self):
        r"""
        Return the hash of ``self``.

        The hash is based on the ground set and an internal Counter of the
        collection of sets.

        This function is called when matroids are added to a set. It is very
        desirable to override it so it can distinguish matroids on the same
        groundset, which is a very typical use case!

        .. WARNING::

            This method is linked to ``__richcmp__`` (in Cython) and
            ``__eq__`` (in Python). If you override one, you must
            override the other!

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets1 = [[0, 1, 2, 3], [1, 2], [3, 4]]
            sage: M1 = TransversalMatroid(sets1)
            sage: M2 = TransversalMatroid(sets1, set_labels=[5, 6, 7])
            sage: sets3 = [['a', 'b', 'c'], ['c', 'd'], ['a', 'd', 'e']]
            sage: M3 = TransversalMatroid(sets3)
            sage: hash(M1) == hash(M2)
            True
            sage: hash(M1) == hash(M3)
            False
        """
        return hash((self._groundset, frozenset(self._sets.items())))

    cdef dict _translate_matching(self):
        """
        Return a Python dictionary that can be used as input in ``__init__()``.
        """
        cdef dict matching = {}
        for x in self._matching:
            matching[self._E[x]] = self._set_labels_input[self._set_labels.index(self._matching[x])]
        return matching

    def __reduce__(self):
        """
        Save the matroid for later reloading.

        OUTPUT:

        A tuple ``(unpickle, (version, data))``, where ``unpickle`` is the
        name of a function that, when called with ``(version, data)``,
        produces a matroid isomorphic to ``self``. ``version`` is an integer
        (currently 0) and ``data`` is a tuple ``(sets, E, name)`` where
        ``E`` is the groundset of the matroid, ``sets`` is the subsets of the
        transversal, and ``name`` is a custom name.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import *
            sage: sets = [range(6)] * 3
            sage: M = TransversalMatroid(sets)
            sage: M == loads(dumps(M))
            True
            sage: M.rename('U36')
            sage: loads(dumps(M))
            U36
        """
        from sage.matroids.unpickling import unpickle_transversal_matroid
        data = (self._sets_input, self._E, self._set_labels_input, self._translate_matching(),
                self.get_custom_name())
        version = 0
        return unpickle_transversal_matroid, (version, data)

    def graph(self):
        """
        Return a bipartite graph representing the transversal matroid.

        A transversal matroid can be represented as a set system, or as a
        bipartite graph with one color class corresponding to the ground set
        and the other to the sets of the set system. This method returns
        that bipartite graph.

        OUTPUT:

        A :class:`Graph`.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: edgedict = {5:[0,1,2,3], 6:[1,2], 7:[1,3,4]}
            sage: B = BipartiteGraph(edgedict)
            sage: M = TransversalMatroid(edgedict.values(), set_labels=edgedict.keys())
            sage: M.graph() == B
            True
            sage: M2 = TransversalMatroid(edgedict.values())
            sage: B2 = M2.graph()
            sage: B2 == B
            False
            sage: B2.is_isomorphic(B)
            True
        """
        # cast the internal networkx as a sage DiGraph
        D = DiGraph(self._D)
        # relabel the vertices, then return as a BipartiteGraph
        vertex_map = {i: e for i, e in enumerate(self._E)}
        for i, l in enumerate(self._set_labels):
            vertex_map[l] = self._set_labels_input[i]
        D.relabel(vertex_map)
        partition = [list(self._E), self._set_labels_input]

        return BipartiteGraph(D, partition=partition)

    cpdef _minor(self, contractions, deletions):
        """
        Return a minor.

        Deletions will yield a new transversal matroid. Contractions will
        have to be a :class:`~sage.matroids.minor_matroid.MinorMatroid`
        until Gammoids are implemented.

        INPUT:

        - ``contractions`` -- frozenset; an independent subset of the ground set
        - ``deletions`` -- frozenset; a coindependent subset of the ground set

        OUTPUT:

        If ``contractions`` is the empty set, or if ``contractions``
        consists of only coloops,  an instance of
        :class:`~sage.matroids.transversal_matroid.TransversalMatroid`.
        Otherwise, an instance of
        :class:`~sage.matroids.minor_matroid.MinorMatroid`.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[0,1,2,3], [1,2], [1,3,4]]
            sage: M1 = TransversalMatroid(sets)
            sage: N1 = M1.delete([2,3])
            sage: sets2 = [[0,1], [1], [4]]
            sage: M2 = TransversalMatroid(sets2)
            sage: N1.is_isomorphic(M2)
            True
            sage: M1._minor(deletions=set([3]), contractions=set([4]))
            M / {4}, where M is Transversal matroid of rank 3 on 4 elements, with 3 sets

        ::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [['a', 'c'], ['e']]
            sage: gs = ['a', 'c', 'd', 'e']
            sage: M = TransversalMatroid(sets, groundset=gs)
            sage: N = M.delete(['d','e']); N
            Transversal matroid of rank 1 on 2 elements, with 1 sets

        TESTS::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [['a', 'c'], ['e']]
            sage: gs = ['a', 'c', 'd', 'e']
            sage: M = TransversalMatroid(sets, groundset=gs)
            sage: M.contract(['e'])
            Transversal matroid of rank 1 on 3 elements, with 1 sets
            sage: M.contract(['d', 'e'])
            Transversal matroid of rank 1 on 2 elements, with 1 sets
            sage: M.contract(['a', 'e'])
            M / {'a', 'e'}, where M is Transversal matroid of rank 2 on 4 elements, with 2 sets
        """
        cdef TransversalMatroid N
        cdef list new_sets, new_set_labels, s, new_s
        cdef size_t i
        cdef frozenset contract = frozenset(contractions)
        cdef frozenset delete = frozenset(deletions)

        # if contractions are just coloops, we can just delete them
        if self._corank(contract) == 0:
            delete = delete.union(contract)
            contract = frozenset()

        if delete:
            new_sets = []
            new_set_labels = []
            for i in range(len(self._sets_input)):
                s = <list> self._sets_input[i]
                new_s = [e for e in s if e not in delete]
                if new_s:
                    # skip over empty buckets, and do bookkeeping with the labels
                    new_sets.append(new_s)
                    new_set_labels.append(self._set_labels_input[i])
            groundset = self._groundset.difference(delete)
            N = TransversalMatroid(new_sets, groundset, new_set_labels)
            # Check if what remains is just coloops
            return N._minor(contractions=contract, deletions=frozenset())
        else:
            N = self

        if contractions:
            # Until gammoids are implemented
            return MinorMatroid(N, contractions=contract, deletions=frozenset())
        else:
            return N

    def set_labels(self):
        """
        Return the labels used for the transversal sets.

        This method will return a list of the labels used of the non-ground
        set vertices of the bipartite graph used to represent the transversal
        matroid. This method does not set anything.

        OUTPUT:

        A list.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: M = TransversalMatroid([[0,1], [1,2,3], [3,4,7]])
            sage: M.set_labels()
            ['s0', 's1', 's2']
            sage: M.graph().vertices()
            ['s0', 's1', 's2', 0, 1, 2, 3, 4, 7]
        """
        return copy(self._set_labels_input)

    cpdef reduce_presentation(self):
        """
        Return an equal transversal matroid where the number of sets equals the rank.

        Every transversal matroid `M` has a presentation with `r(M)` sets, and if `M`
        has no coloops, then every presentation has `r(M)` nonempty sets. This method
        discards extra sets if `M` has coloops.

        OUTPUT:

        A ``TransversalMatroid`` instance with a reduced presentation.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[0,1], [2], [2]]
            sage: M = TransversalMatroid(sets); M
            Transversal matroid of rank 2 on 3 elements, with 3 sets
            sage: N = M.reduce_presentation(); N
            Transversal matroid of rank 2 on 3 elements, with 2 sets
            sage: N.equals(M)
            True
            sage: N == M
            False
            sage: sets = [[0,1], [], [], [2]]
            sage: M1 = TransversalMatroid(sets); M1
            Transversal matroid of rank 2 on 3 elements, with 4 sets
            sage: M1.reduce_presentation()
            Transversal matroid of rank 2 on 3 elements, with 2 sets

        ::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[0, 1, 2, 3]] * 3
            sage: M = TransversalMatroid(sets)
            sage: N = M.reduce_presentation()
            sage: M == N
            True

        TESTS::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[4], [1,3], [4], [0,1], [2,3], [1]]
            sage: M = TransversalMatroid(sets)
            sage: M1 = M.reduce_presentation(); M1
            Transversal matroid of rank 5 on 5 elements, with 5 sets
            sage: len(M1.graph().edges())
            5
        """
        cdef TransversalMatroid N
        cdef size_t i
        if len(self._sets_input) == self.full_rank():
            return self

        element_int_map = {e: i for i, e in enumerate(self._groundset)}
        coloops = self.coloops()
        cdef list coloops_to_delete = [e for e in coloops if self._D.degree(element_int_map[e]) > 1]
        N = <TransversalMatroid?> self._minor(contractions=set(), deletions=set(coloops_to_delete))
        cdef list sets = list(N._sets_input) # make a (shallow) copy as we will remove objects
        # reuse the old set labels
        # this does not respect containment
        labels = N.set_labels()
        # remove empty sets HERE
        if [] in sets:
            i = len(sets)
            while i > 0:
                i -= 1
                if not sets[i]:
                    del sets[i]
                    del labels[i]
        cdef set free_labels = set(self._set_labels_input).difference(labels)
        for c in coloops_to_delete:
            l = free_labels.pop()
            sets.append([c])
            labels.append(l)
        return TransversalMatroid(sets, groundset=self.groundset(), set_labels=labels)

    cpdef transversal_extension(self, element=None, newset=False, sets=[]):
        r"""
        Return a TransversalMatroid extended by an element.

        This will modify the presentation of the transversal matroid by
        adding a new element, and placing this element in the specified
        sets. It is also possible to use this method to create a new set
        that will have the new element as its only member, making it a coloop.

        INPUT:

        - ``element`` -- (optional) the name for the new element
        - ``newset`` -- (optional) if specified, the element will be
          given its own set
        - ``sets`` -- (default: ``None``) an iterable of labels representing the
          sets in the current presentation that the new element will belong to

        OUTPUT:

        A :class:`~sage.matroids.transversal_matroids.TransversalMatroid`
        with a ground set element added to specified sets. Note that the
        ``newset`` option will make the new element a coloop. If
        ``newset == True``, a name will be generated; otherwise the
        value of ``newset`` will be used.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: M = TransversalMatroid([['a', 'c']], groundset=['a', 'c'], set_labels=['b'])
            sage: M1 = M.transversal_extension(element='d', newset='e')
            sage: M2 = M.transversal_extension(element='d', newset=True)
            sage: M1.coloops()
            frozenset({'d'})
            sage: True in M2.graph().vertices()
            False
            sage: M1.is_isomorphic(M2)
            True
            sage: M3 = M.transversal_extension('d', sets=['b'])
            sage: M3.is_isomorphic(matroids.Uniform(1, 3))
            True
            sage: M4 = M.transversal_extension('d', sets=['a'])
            Traceback (most recent call last):
            ...
            ValueError: sets do not match presentation
            sage: M4 = M.transversal_extension('a', sets=['b'])
            Traceback (most recent call last):
            ...
            ValueError: cannot extend by element already in ground set
            sage: M2.transversal_extension(newset='b')
            Traceback (most recent call last):
            ...
            ValueError: newset is already a vertex in the presentation
            sage: M5 = M1.transversal_extension()
            sage: len(M5.loops())
            1
            sage: M2.transversal_extension(element='b')
            Transversal matroid of rank 2 on 4 elements, with 2 sets

        ::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[0, 1, 2, 3], [1, 2], [1, 3, 4]]
            sage: M = TransversalMatroid(sets, groundset=range(5), set_labels=[5, 6, 7])
            sage: N = M.delete(2)
            sage: M1 = N.transversal_extension(element=2, sets=[5, 6])
            sage: M1 == M
            True

        ::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[0, 1, 2, 3]] * 3
            sage: M = TransversalMatroid(sets, set_labels=[4, 5, 6])
            sage: N = M.transversal_extension(element='a', newset=True, sets=[4])
            sage: N.graph().degree('a')
            2

        TESTS::

            sage: N = TransversalMatroid(['abc', 'abd', 'cde']); N
            Transversal matroid of rank 3 on 5 elements, with 3 sets
            sage: Ne = N.transversal_extension(element='f', sets=['s2'])
        """
        cdef set parsed_sets = set(sets)
        if element is None:
            element = newlabel(self._groundset)
        elif element in self._groundset:
            raise ValueError("cannot extend by element already in ground set")
        cdef list labels = self._set_labels_input
        if not parsed_sets.issubset(labels):
            raise ValueError("sets do not match presentation")

        # check for conflicts with new labels
        labels_map = {l: l for l in labels}
        if element in labels:
            new_label = newlabel(self._groundset.union(labels).union([newset]))
            labels_map[element] = new_label

        # newset should not be a ground set element or existing set
        if newset in self._E or newset in self._set_labels_input:
            # keywords `True` and `False` give us problems here
            if not isinstance(newset, bool):
                raise ValueError("newset is already a vertex in the presentation")

        cdef list s, new_sets = []
        cdef size_t i
        for i in range(len(self._sets_input)):
            s = <list> self._sets_input[i]
            if labels[i] in parsed_sets:
                new_sets.append(s + [element])
            else:
                new_sets.append(s)

        if newset:
            if newset is True:
                newset = newlabel(self._groundset.union(labels))
            new_sets.append([element])
            labels = list(labels)  # Make a shallow copy since we mutate it
            labels.append(newset)

        groundset = self._groundset.union([element])
        return TransversalMatroid(new_sets, groundset, labels)

    def transversal_extensions(self, element=None, sets=[]):
        r"""
        Return an iterator of extensions based on the transversal presentation.

        This method will take an extension by adding an element to every possible
        sub-collection of the collection of desired sets. No checking is done
        for equal matroids. It is advised to make sure the presentation has as
        few sets as possible by using
        :meth:`reduce_presentation() <sage.matroids.transversal_matroid.TransversalMatroid.reduce_presentation>`

        INPUT:

        - ``element`` -- (optional) the name of the new element
        - ``sets`` -- (optional) a list containing names of sets in the matroid's
          presentation.

        OUTPUT:

        An iterator of ``TransversalMatroids``.

        If ``sets`` is not specified, every set will be used.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[3, 4, 5, 6]] * 3
            sage: M = TransversalMatroid(sets, set_labels=[0, 1, 2])
            sage: len([N for N in M.transversal_extensions()])
            8
            sage: len([N for N in M.transversal_extensions(sets=[0, 1])])
            4
            sage: set(sorted([N.groundset() for N in M.transversal_extensions(element=7)]))
            {frozenset({3, 4, 5, 6, 7})}
        """
        if element is None:
            element = newlabel(self.groundset())
        elif element in self._groundset:
            raise ValueError("cannot extend by element already in ground set")

        labels = self._set_labels_input
        if not sets:
            sets = labels
        elif not set(sets).issubset(labels):
            raise ValueError("sets do not match presentation")

        # Adapted from the Python documentation
        from itertools import chain, combinations
        powerset = chain.from_iterable(combinations(sets, r) for r in range(len(sets)+1))

        for collection in powerset:
            yield self.transversal_extension(element=element, sets=collection)

    cpdef is_valid(self, certificate=False):
        """
        Test whether the matching in memory is a valid maximal matching.

        The data for a transversal matroid is a set system, which is always valid,
        but it is possible for a user to provide invalid input with the ``matching``
        parameter. This checks that the matching provided is indeed a matching, fits in
        the set system, and is maximal.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import *
            sage: sets = [[0, 1, 2, 3], [1, 2], [1, 3, 4]]
            sage: set_labels = [5, 6, 7]
            sage: M = TransversalMatroid(sets, set_labels=set_labels)
            sage: M.is_valid()
            True
            sage: m = {0:5, 1:5, 3:7} # not a matching
            sage: TransversalMatroid(sets, set_labels=set_labels, matching=m).is_valid()
            False
            sage: m = {2:6, 3:7} # not maximal
            sage: TransversalMatroid(sets, set_labels=set_labels, matching=m).is_valid()
            False
            sage: m = {0:6, 1:5, 3:7} # not in the set system
            sage: TransversalMatroid(sets, set_labels=set_labels, matching=m).is_valid()
            False
        """
        # Check that self._matching is a matching, that is, every entry appears at most once.
        # Because of how python dictionaries work, each element appears once, but a set
        # can appear more than once.
        if len(self._matching) != len(set(self._matching.values())):
            return False

        # Check that every element in the matching is actually in the set that the matching
        # specifies.
        # Note that self._matching uses the internal set and element labels
        for e in self._matching:
            if e not in self._sets_input[self._set_labels.index(self._matching[e])]:
                return False

        # Check that the matching is maximal.
        B = self.graph()
        if len(B.matching()) != len(self._matching):
            return False

        return True

