r"""
Incidence structures (i.e. hypergraphs, i.e. set systems)

An incidence structure is specified by a list of points, blocks, or an incidence
matrix ([1]_, [2]_). :class:`IncidenceStructure` instances have the following methods:

{METHODS_OF_IncidenceStructure}

REFERENCES:

.. [1] Block designs and incidence structures from wikipedia,
  :wikipedia:`Block_design`
  :wikipedia:`Incidence_structure`

.. [2] \E. Assmus, J. Key, Designs and their codes, CUP, 1992.

AUTHORS:

- Peter Dobcsanyi and David Joyner (2007-2008)

  This is a significantly modified form of part of the module block_design.py
  (version 0.6) written by Peter Dobcsanyi peter@designtheory.org.

- Vincent Delecroix (2014): major rewrite

Methods
-------
"""
# **************************************************************************
#                              Copyright (C) 2007                          #
#                                                                          #
#                Peter Dobcsanyi       and         David Joyner            #
#           <peter@designtheory.org>          <wdjoyner@gmail.com>         #
#                                                                          #
#                                                                          #
#    Distributed under the terms of the GNU General Public License (GPL)   #
#    as published by the Free Software Foundation; either version 2 of     #
#    the License, or (at your option) any later version.                   #
#                    https://www.gnu.org/licenses/                          #
# **************************************************************************
from __future__ import annotations

from sage.misc.latex import latex
from sage.misc.lazy_import import lazy_import
from sage.rings.integer import Integer
from sage.sets.set import Set

lazy_import('sage.libs.gap.libgap', 'libgap')


class IncidenceStructure:
    r"""
    A base class for incidence structures (i.e. hypergraphs, i.e. set systems)

    An incidence structure (i.e. hypergraph, i.e. set system) can be defined
    from a collection of blocks (i.e. sets, i.e. edges), optionally with an
    explicit ground set (i.e. point set, i.e. vertex set). Alternatively they
    can be defined from a binary incidence matrix.

    INPUT:

    - ``points`` -- (i.e. ground set, i.e. vertex set) the underlying set. If
      ``points`` is an integer `v`, then the set is considered to be `\{0, ...,
      v-1\}`.

      .. NOTE::

          The following syntax, where ``points`` is omitted, automatically
          defines the ground set as the union of the blocks::

              sage: H = IncidenceStructure([['a','b','c'],['c','d','e']])
              sage: sorted(H.ground_set())
              ['a', 'b', 'c', 'd', 'e']

    - ``blocks`` -- (i.e. edges, i.e. sets) the blocks defining the incidence
      structure; can be any iterable

    - ``incidence_matrix`` -- a binary incidence matrix; each column represents
      a set

    - ``name`` -- string (such as "Fano plane")

    - ``check`` -- whether to check the input

    - ``copy`` -- (use with caution) if set to ``False`` then ``blocks`` must be
      a list of lists of integers. The list will not be copied but will be
      modified in place (each block is sorted, and the whole list is
      sorted). Your ``blocks`` object will become the
      :class:`IncidenceStructure` instance's internal data.

    EXAMPLES:

    An incidence structure can be constructed by giving the number of points and
    the list of blocks::

        sage: IncidenceStructure(7, [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
        Incidence structure with 7 points and 7 blocks

    Only providing the set of blocks is sufficient. In this case, the ground set
    is defined as the union of the blocks::

        sage: IncidenceStructure([[1,2,3],[2,3,4]])
        Incidence structure with 4 points and 2 blocks

    Or by its adjacency matrix (a `\{0,1\}`-matrix in which rows are indexed by
    points and columns by blocks)::

        sage: m = matrix([[0,1,0],[0,0,1],[1,0,1],[1,1,1]])                             # needs sage.modules
        sage: IncidenceStructure(m)                                                     # needs sage.modules
        Incidence structure with 4 points and 3 blocks

    The points can be any (hashable) object::

        sage: V = [(0,'a'),(0,'b'),(1,'a'),(1,'b')]
        sage: B = [(V[0],V[1],V[2]), (V[1],V[2]), (V[0],V[2])]
        sage: I = IncidenceStructure(V, B)
        sage: I.ground_set()
        [(0, 'a'), (0, 'b'), (1, 'a'), (1, 'b')]
        sage: I.blocks()
        [[(0, 'a'), (0, 'b'), (1, 'a')], [(0, 'a'), (1, 'a')], [(0, 'b'), (1, 'a')]]

    The order of the points and blocks does not matter as they are sorted on
    input (see :issue:`11333`)::

        sage: A = IncidenceStructure([0,1,2], [[0],[0,2]])
        sage: B = IncidenceStructure([1,0,2], [[0],[2,0]])
        sage: B == A
        True

        sage: C = BlockDesign(2, [[0], [1,0]])
        sage: D = BlockDesign(2, [[0,1], [0]])
        sage: C == D
        True

    If you care for speed, you can set ``copy`` to ``False``, but in that
    case, your input must be a list of lists and the ground set must be `{0,
    ..., v-1}`::

        sage: blocks = [[0,1],[2,0],[1,2]]  # a list of lists of integers
        sage: I = IncidenceStructure(3, blocks, copy=False)
        sage: I._blocks is blocks
        True
    """
    def __init__(self, points=None, blocks=None, incidence_matrix=None,
            name=None, check=True, copy=True):
        r"""
        TESTS::

            sage: IncidenceStructure(3, [[4]])
            Traceback (most recent call last):
            ...
            ValueError: Block [4] is not contained in the point set

            sage: IncidenceStructure(3, [[0,1],[0,2]], check=True)
            Incidence structure with 3 points and 2 blocks

            sage: IncidenceStructure(2, [[0,1,2,3,4,5]], check=False)
            Incidence structure with 2 points and 1 blocks

        We avoid to convert to integers when the points are not (but compare
        equal to integers because of coercion)::

            sage: # needs sage.rings.finite_rings
            sage: V = GF(5)
            sage: e0,e1,e2,e3,e4 = V
            sage: [e0,e1,e2,e3,e4] == list(range(5))  # coercion makes them equal
            True
            sage: blocks = [[e0,e1,e2],[e0,e1],[e2,e4]]
            sage: I = IncidenceStructure(V, blocks)
            sage: type(I.ground_set()[0])
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
            sage: type(I.blocks()[0][0])
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>

        TESTS::

            sage: IncidenceStructure([])
            Incidence structure with 0 points and 0 blocks
        """
        from sage.structure.element import Matrix

        # Reformatting input
        if isinstance(points, Matrix):
            assert incidence_matrix is None, "'incidence_matrix' cannot be defined when 'points' is a matrix"
            assert blocks is None, "'blocks' cannot be defined when 'points' is a matrix"
            incidence_matrix = points
            points = blocks = None
        elif points is not None and blocks is None:
            blocks = points
            points = set().union(*blocks)
        if points:
            assert incidence_matrix is None, "'incidence_matrix' cannot be defined when 'points' is defined"

        if incidence_matrix:
            from sage.matrix.constructor import matrix
            M = matrix(incidence_matrix)
            v = M.nrows()
            self._points = list(range(v))
            self._point_to_index = None
            self._blocks = sorted(M.nonzero_positions_in_column(i) for i in range(M.ncols()))

        else:
            if isinstance(points, (int, Integer)):
                self._points = list(range(points))
                self._point_to_index = None
            else:
                self._points = list(points)
                if self._points == list(range(len(points))) and all(isinstance(x, (int, Integer)) for x in self._points):
                    self._point_to_index = None
                else:
                    self._point_to_index = {e: i for i, e in enumerate(self._points)}

            if check:
                for block in blocks:
                    if any(x not in self._points for x in block):
                        raise ValueError("Block {} is not contained in the point set".format(block))
                    if len(block) != len(set(block)):
                        raise ValueError("Repeated element in block {}".format(block))

            if self._point_to_index:
                # translate everything to integers between 0 and v-1
                blocks = [sorted(self._point_to_index[e] for e in block) for block in blocks]
            elif copy:
                # create a new list made of sorted blocks
                blocks = [sorted(block) for block in blocks]
            else:
                # sort the data but avoid copying it
                for b in blocks:
                    b.sort()

            blocks.sort()
            self._blocks = blocks

        self._name = str(name) if name is not None else 'IncidenceStructure'
        self._classes = None
        self._canonical_label = None

    def __iter__(self):
        """
        Iterator over the blocks.

        EXAMPLES::

            sage: sts = designs.steiner_triple_system(9)
            sage: list(sts)
            [[0, 1, 5], [0, 2, 4], [0, 3, 6], [0, 7, 8], [1, 2, 3], [1, 4, 7],
            [1, 6, 8], [2, 5, 8], [2, 6, 7], [3, 4, 8], [3, 5, 7], [4, 5, 6]]

            sage: b = IncidenceStructure('ab', ['a','ab'])
            sage: it = iter(b)
            sage: next(it)
            ['a']
            sage: next(it)
            ['a', 'b']
        """
        if self._point_to_index is None:
            for b in self._blocks:
                yield b[:]
        else:
            for b in self._blocks:
                yield [self._points[i] for i in b]

    def __repr__(self):
        """
        A print method.

        EXAMPLES::

            sage: BD = IncidenceStructure(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD
            Incidence structure with 7 points and 7 blocks
        """
        return 'Incidence structure with {} points and {} blocks'.format(
                self.num_points(), self.num_blocks())

    __str__ = __repr__

    def __eq__(self, other):
        """
        Test whether the two incidence structures are equal.

        TESTS::

            sage: blocks = [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]]
            sage: BD1 = IncidenceStructure(7, blocks)
            sage: M = BD1.incidence_matrix()                                            # needs sage.modules
            sage: BD2 = IncidenceStructure(incidence_matrix=M)                          # needs sage.modules
            sage: BD1 == BD2                                                            # needs sage.modules
            True

            sage: e1 = frozenset([0,1])
            sage: e2 = frozenset([2])
            sage: sorted([e1,e2]) == [e1,e2]
            True
            sage: sorted([e2,e1]) == [e2,e1]
            True
            sage: I1 = IncidenceStructure([e1,e2], [[e1],[e1,e2]])
            sage: I2 = IncidenceStructure([e1,e2], [[e2,e1],[e1]])
            sage: I3 = IncidenceStructure([e2,e1], [[e1,e2],[e1]])
            sage: I1 == I2 and I2 == I1 and I1 == I3 and I3 == I1 and I2 == I3 and I3 == I2
            True
        """
        # We are extra careful in this method since we cannot assume that a
        # total order is defined on the point set.
        if not isinstance(other, IncidenceStructure):
            return False

        if self._points == other._points:
            return self._blocks == other._blocks

        if (self.num_points() != other.num_points() or
            self.num_blocks() != other.num_blocks()):
            return False

        p_to_i = self._point_to_index if self._point_to_index else list(range(self.num_points()))

        if any(p not in p_to_i for p in other.ground_set()):
            return False

        other_blocks = sorted(sorted(p_to_i[p] for p in b) for b in other.blocks())
        return self._blocks == other_blocks

    def __ne__(self, other):
        r"""
        Difference test.

        EXAMPLES::

            sage: BD1 = IncidenceStructure(7, [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: M = BD1.incidence_matrix()                                            # needs sage.modules
            sage: BD2 = IncidenceStructure(incidence_matrix=M)                          # needs sage.modules
            sage: BD1 != BD2                                                            # needs sage.modules
            False
        """
        return not self == other

    def __contains__(self, block):
        r"""
        Test if a block belongs to the incidence structure.

        INPUT:

        - ``block`` -- a block

        EXAMPLES::

            sage: [1,2,3,4] in IncidenceStructure([[1,2,3,4]])
            True
            sage: [1,2,4,3] in IncidenceStructure([[1,2,3,4]])
            True
            sage: [1,2,"3",4] in IncidenceStructure([[1,2,3,4]])
            False
            sage: [1,2,"3",4] in IncidenceStructure([[1,2,"3",4]])
            True

        More complicated examples::

            sage: str="I had a dream of a time when a 3-lines patch does not kill one hour"
            sage: sets = Subsets(str.split(), 4)
            sage: IS = IncidenceStructure(sets)  # a complete 4-uniform hypergraph
            sage: ["I", "dream", "of", "one"] in IS
            True
            sage: ["does", "patch", "kill", "dream"] in IS
            True
            sage: ["Am", "I", "finally", "done ?"] in IS
            False
            sage: IS = designs.ProjectiveGeometryDesign(3, 1, GF(2),                    # needs sage.combinat
            ....:                                       point_coordinates=False)
            sage: [3,8,7] in IS                                                         # needs sage.combinat
            True
            sage: [3,8,9] in IS                                                         # needs sage.combinat
            False
        """
        try:
            iter(block)
        except TypeError:
            return False

        # Relabel to 0,...,n-1 if necessary
        if self._point_to_index is not None:
            try:
                block = [self._point_to_index[x] for x in block]
            except KeyError:
                return False

        return sorted(block) in self._blocks

    def canonical_label(self):
        r"""
        Return a canonical label for the incidence structure.

        A canonical label is relabeling of the points into integers
        `\{0,...,n-1\}` such that isomorphic incidence structures are
        relabelled to equal objects.

        EXAMPLES::

            sage: # needs sage.schemes
            sage: fano1 = designs.balanced_incomplete_block_design(7,3)
            sage: fano2 = designs.projective_plane(2)
            sage: fano1 == fano2
            False
            sage: fano1.relabel(fano1.canonical_label())
            sage: fano2.relabel(fano2.canonical_label())
            sage: fano1 == fano2
            True
        """
        if self._canonical_label is None:
            from sage.graphs.graph import Graph
            g = Graph()
            n = self.num_points()
            g.add_edges((i+n,x) for i,b in enumerate(self._blocks) for x in b)
            canonical_label = g.canonical_label([list(range(n)),list(range(n,n+self.num_blocks()))],certificate=True)[1]
            canonical_label = [canonical_label[x] for x in range(n)]
            self._canonical_label = canonical_label

        return dict(zip(self._points,self._canonical_label))

    def is_isomorphic(self, other, certificate=False):
        r"""
        Return whether the two incidence structures are isomorphic.

        INPUT:

        - ``other`` -- an incidence structure

        - ``certificate`` -- boolean (default: ``False``); whether to return an
          isomorphism from ``self`` to ``other`` instead of a boolean answer

        EXAMPLES::

            sage: # needs sage.schemes
            sage: fano1 = designs.balanced_incomplete_block_design(7,3)
            sage: fano2 = designs.projective_plane(2)
            sage: fano1.is_isomorphic(fano2)
            True
            sage: fano1.is_isomorphic(fano2,certificate=True)
            {0: 0, 1: 1, 2: 2, 3: 6, 4: 4, 5: 3, 6: 5}

        TESTS::

            sage: # needs sage.symbolic
            sage: IS  = IncidenceStructure([["A",5,pi],["A",5,"Wouhou"],
            ....:                           ["A","Wouhou",(9,9)],[pi,12]])
            sage: IS2 = IS.copy()
            sage: IS2.relabel(IS2.canonical_label())
            sage: IS.is_isomorphic(IS2)
            True
            sage: canon = IS.is_isomorphic(IS2, certificate=True)
            sage: IS.relabel(canon)
            sage: IS==IS2
            True

            sage: IS2 = IncidenceStructure([[1,2]])
            sage: IS2.is_isomorphic(IS)                                                 # needs sage.symbolic
            False
            sage: IS2.is_isomorphic(IS, certificate=True)                               # needs sage.symbolic
            {}

        Checking whether two :class:`IncidenceStructure` are isomorphic
        incidentally computes their canonical label (if necessary). Thus,
        subsequent calls to :meth:`is_isomorphic` will be faster::

            sage: # needs sage.schemes
            sage: IS1 = designs.projective_plane(3)
            sage: IS2 = IS1.relabel(Permutations(IS1.ground_set()).random_element(),
            ....:                   inplace=False)
            sage: IS2 = IncidenceStructure(IS2.blocks())
            sage: IS1._canonical_label is None and IS2._canonical_label is None
            True
            sage: IS1.is_isomorphic(IS2)
            True
            sage: IS1._canonical_label is None or IS2._canonical_label is None
            False
        """
        if (self.num_points() != other.num_points() or
            self.num_blocks() != other.num_blocks() or
            sorted(self.block_sizes()) != sorted(other.block_sizes())):
            return {} if certificate else False

        A_canon = self.canonical_label()
        B_canon = other.canonical_label()

        A = self.relabel(A_canon,inplace=False)
        B = other.relabel(B_canon,inplace=False)

        if A == B:
            if certificate:
                B_canon_rev = {y:x for x,y in B_canon.items()}
                return {x:B_canon_rev[xint] for x,xint in A_canon.items()}
            else:
                return True
        else:
            return {} if certificate else False

    def isomorphic_substructures_iterator(self, H2,induced=False):
        r"""
        Iterate over all copies of ``H2`` contained in ``self``.

        A hypergraph `H_1` contains an isomorphic copy of a hypergraph `H_2` if
        there exists an injection `f:V(H_2)\mapsto V(H_1)` such that for any set
        `S_2\in E(H_2)` the set `S_1=f(S2)` belongs to `E(H_1)`.

        It is an *induced* copy if no other set of `E(H_1)` is contained in
        `f(V(H_2))`, i.e. `|E(H_2)|=\{S:S\in E(H_1)\text{ and }f(V(H_2))\}`.

        This function lists all such injections. In particular, the number of
        copies of `H` in itself is equal to *the size of its automorphism
        group*.

        See :mod:`~sage.combinat.designs.subhypergraph_search` for more information.

        INPUT:

        - ``H2`` -- an :class:`IncidenceStructure` object

        - ``induced`` -- boolean (default: ``False``); whether to require the copies to be
          induced

        EXAMPLES:

        How many distinct `C_5` in Petersen's graph ? ::

            sage: P = graphs.PetersenGraph()
            sage: C = graphs.CycleGraph(5)
            sage: IP = IncidenceStructure(P.edges(sort=True, labels=False))
            sage: IC = IncidenceStructure(C.edges(sort=True, labels=False))
            sage: sum(1 for _ in IP.isomorphic_substructures_iterator(IC))
            120

        As the automorphism group of `C_5` has size 10, the number of distinct
        unlabelled copies is 12. Let us check that all functions returned
        correspond to an actual `C_5` subgraph::

            sage: for f in IP.isomorphic_substructures_iterator(IC):
            ....:     assert all(P.has_edge(f[x],f[y]) for x,y in C.edges(sort=True, labels=False))

        The number of induced copies, in this case, is the same::

            sage: sum(1 for _ in IP.isomorphic_substructures_iterator(IC,induced=True))
            120

        They begin to differ if we make one vertex universal::

            sage: P.add_edges([(0,x) for x in P], loops=False)
            sage: IP = IncidenceStructure(P.edges(sort=True, labels=False))
            sage: IC = IncidenceStructure(C.edges(sort=True, labels=False))
            sage: sum(1 for _ in IP.isomorphic_substructures_iterator(IC))
            420
            sage: sum(1 for _ in IP.isomorphic_substructures_iterator(IC,induced=True))
            60

        The number of copies of `H` in itself is the size of its automorphism
        group::

            sage: H = designs.projective_plane(3)                                       # needs sage.schemes
            sage: sum(1 for _ in H.isomorphic_substructures_iterator(H))                # needs sage.schemes
            5616
            sage: H.automorphism_group().cardinality()                                  # needs sage.groups sage.schemes
            5616
        """
        from sage.combinat.designs.subhypergraph_search import SubHypergraphSearch
        return SubHypergraphSearch(self,H2,induced=induced)

    def copy(self):
        r"""
        Return a copy of the incidence structure.

        EXAMPLES::

            sage: IS = IncidenceStructure([[1,2,3,"e"]], name='Test')
            sage: IS
            Incidence structure with 4 points and 1 blocks
            sage: copy(IS)
            Incidence structure with 4 points and 1 blocks
            sage: [1, 2, 3, 'e'] in copy(IS)
            True
            sage: copy(IS)._name
            'Test'
        """
        IS = IncidenceStructure(self._blocks,
                                name=self._name,
                                check=False)
        IS.relabel(dict(zip(range(self.num_points()),self._points)))
        IS._canonical_label = None if self._canonical_label is None else self._canonical_label[:]

        return IS

    __copy__ = copy

    def induced_substructure(self, points):
        r"""
        Return the substructure induced by a set of points.

        The substructure induced in `\mathcal H` by a set `X\subseteq V(\mathcal
        H)` of points is the incidence structure `\mathcal H_X` defined on `X`
        whose sets are all `S\in \mathcal H` such that `S\subseteq X`.

        INPUT:

        - ``points`` -- set of points

        .. NOTE::

            This method goes over all sets of ``self`` before building a new
            :class:`IncidenceStructure` (which involves some relabelling and
            sorting). It probably should not be called in a performance-critical
            code.

        EXAMPLES:

        A Fano plane with one point removed::

            sage: F = designs.steiner_triple_system(7)
            sage: F.induced_substructure([0..5])
            Incidence structure with 6 points and 4 blocks

        TESTS::

            sage: F.induced_substructure([0..50])
            Traceback (most recent call last):
            ...
            ValueError: 7 is not a point of the incidence structure
            sage: F.relabel(dict(enumerate("abcdefg")))
            sage: F.induced_substructure("abc")
            Incidence structure with 3 points and ...
            sage: F.induced_substructure("Y")
            Traceback (most recent call last):
            ...
            ValueError: 'Y' is not a point of the incidence structure
        """
        # Checking the input
        if self._point_to_index is None:
            n = self.num_points()
            for x in points:
                x = int(x)
                if x < 0 or x >= n:
                    raise ValueError("{} is not a point of the incidence structure".format(x))
            int_points = points
        else:
            try:
                int_points = [self._point_to_index[x] for x in points]
            except KeyError as bad_pt:
                raise ValueError("{} is not a point of the incidence structure".format(bad_pt))

        int_points = set(int_points)
        return IncidenceStructure(points,
                                  [[self._points[x] for x in S]
                                   for S in self._blocks
                                   if int_points.issuperset(S)])

    def trace(self, points, min_size=1, multiset=True):
        r"""
        Return the trace of a set of points.

        Given an hypergraph `\mathcal H`, the *trace* of a set `X` of points in
        `\mathcal H` is the hypergraph whose blocks are all non-empty `S \cap X`
        where `S \in \mathcal H`.

        INPUT:

        - ``points`` -- set of points

        - ``min_size`` -- integer (default: 1); minimum size of the sets to
          keep. By default all empty sets are discarded, i.e. ``min_size=1``

        - ``multiset`` -- boolean (default: ``True``); whether to keep multiple
          copies of the same set

        .. NOTE::

            This method goes over all sets of ``self`` before building a new
            :class:`IncidenceStructure` (which involves some relabelling and
            sorting). It probably should not be called in a performance-critical
            code.

        EXAMPLES:

        A Baer subplane of order 2 (i.e. a Fano plane) in a projective plane of order 4::

            sage: # needs sage.schemes
            sage: P4 = designs.projective_plane(4)
            sage: F = designs.projective_plane(2)
            sage: for x in Subsets(P4.ground_set(),7):
            ....:     if P4.trace(x,min_size=2).is_isomorphic(F):
            ....:         break
            sage: subplane = P4.trace(x,min_size=2); subplane
            Incidence structure with 7 points and 7 blocks
            sage: subplane.is_isomorphic(F)
            True

        TESTS::

            sage: # needs sage.schemes
            sage: F.trace([0..50])
            Traceback (most recent call last):
            ...
            ValueError: 7 is not a point of the incidence structure
            sage: F.relabel(dict(enumerate("abcdefg")))
            sage: F.trace("abc")
            Incidence structure with 3 points and ...
            sage: F.trace("Y")
            Traceback (most recent call last):
            ...
            ValueError: 'Y' is not a point of the incidence structure
        """
        # Checking the input
        if self._point_to_index is None:
            n = self.num_points()
            int_points = frozenset(int(x) for x in points)
            for x in int_points:
                if x < 0 or x >= n:
                    raise ValueError("{} is not a point of the incidence structure".format(x))
        else:
            try:
                int_points = frozenset(self._point_to_index[x] for x in points)
            except KeyError as bad_pt:
                raise ValueError("{} is not a point of the incidence structure".format(bad_pt))

        blocks = [int_points.intersection(S) for S in self._blocks]
        if min_size:
            blocks = [S for S in blocks if len(S) >= min_size]
        if not multiset:
            blocks = set(blocks)
        IS = IncidenceStructure(blocks)
        IS.relabel({i:self._points[i] for i in int_points})
        return IS

    def ground_set(self):
        r"""
        Return the ground set (i.e the list of points).

        EXAMPLES::

            sage: IncidenceStructure(3, [[0,1],[0,2]]).ground_set()
            [0, 1, 2]
        """
        return self._points[:]

    def num_points(self):
        r"""
        Return the size of the ground set.

        EXAMPLES::

            sage: designs.DesarguesianProjectivePlaneDesign(2).num_points()
            7
            sage: B = IncidenceStructure(4, [[0,1],[0,2],[0,3],[1,2], [1,2,3]])
            sage: B.num_points()
            4
        """
        return len(self._points)

    def num_blocks(self):
        r"""
        Return the number of blocks.

        EXAMPLES::

            sage: designs.DesarguesianProjectivePlaneDesign(2).num_blocks()
            7
            sage: B = IncidenceStructure(4, [[0,1],[0,2],[0,3],[1,2], [1,2,3]])
            sage: B.num_blocks()
            5
        """
        return len(self._blocks)

    def blocks(self):
        """
        Return the list of blocks.

        EXAMPLES::

            sage: BD = IncidenceStructure(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.blocks()
            [[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]
        """
        if self._point_to_index is None:
            return [b[:] for b in self._blocks]
        else:
            return [[self._points[i] for i in b] for b in self._blocks]

    def block_sizes(self):
        r"""
        Return the set of block sizes.

        EXAMPLES::

            sage: BD = IncidenceStructure(8, [[0,1,3],[1,4,5,6],[1,2],[5,6,7]])
            sage: BD.block_sizes()
            [3, 2, 4, 3]
            sage: BD = IncidenceStructure(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.block_sizes()
            [3, 3, 3, 3, 3, 3, 3]
        """
        return [len(b) for b in self._blocks]

    def degree(self, p=None, subset=False):
        r"""
        Return the degree of a point ``p`` (or a set of points).

        The degree of a point (or set of points) is the number of blocks that
        contain it.

        INPUT:

        - ``p`` -- a point (or a set of points) of the incidence structure

        - ``subset`` -- boolean (default: ``False``); whether to interpret the
          argument as a set of point or as a point (default)

        EXAMPLES::

            sage: designs.steiner_triple_system(9).degree(3)
            4
            sage: designs.steiner_triple_system(9).degree({1,2},subset=True)
            1

        TESTS::

            sage: designs.steiner_triple_system(9).degree(subset=True)
            Traceback (most recent call last):
            ...
            ValueError: subset must be False when p is None
        """
        if p is None:
            if subset is True:
                raise ValueError("subset must be False when p is None")

        # degree of a point
        if not subset:
            if self._point_to_index:
                p = self._point_to_index.get(p,-1)
            else:
                p = p if (p >= 0 and p < len(self._points)) else -1
            return sum((p in b) for b in self._blocks) if p != -1 else 0

        # degree of a set
        else:
            if self._point_to_index:
                p = set(self._point_to_index.get(x,-1) for x in p)
            else:
                p = set(p) if all(x >= 0 and x < len(self._points) for x in p) else set([-1])

            return sum(p.issubset(b) for b in self._blocks) if -1 not in p else 0

    def degrees(self, size=None):
        r"""
        Return the degree of all sets of given size, or the degree of all points.

        The degree of a point (or set of point) is the number of blocks that
        contain it.

        INPUT:

        - ``size`` -- integer; return the degree of all subsets of points of
          cardinality ``size``. When ``size=None``, the function outputs the
          degree of all points.

          .. NOTE::

              When ``size=None`` the output is indexed by the points. When
              ``size=1`` it is indexed by tuples of size 1. This is the same
              information, stored slightly differently.

        OUTPUT: a dictionary whose values are degrees and keys are either:

        - the points of the incidence structure if ``size=None`` (default)

        - the subsets of size ``size`` of the points stored as tuples

        EXAMPLES::

            sage: IncidenceStructure([[1,2,3],[1,4]]).degrees(2)
            {(1, 2): 1, (1, 3): 1, (1, 4): 1, (2, 3): 1, (2, 4): 0, (3, 4): 0}

        In a Steiner triple system, all pairs have degree 1::

            sage: S13 = designs.steiner_triple_system(13)
            sage: all(v == 1 for v in S13.degrees(2).values())
            True
        """
        if size is None:
            d = [0]*self.num_points()
            for b in self._blocks:
                for x in b:
                    d[x] += 1
            return {p: d[i] for i, p in enumerate(self._points)}
        else:
            from itertools import combinations
            d = {t:0 for t in combinations(range(self.num_points()),size)}
            for b in self._blocks:
                for s in combinations(b,size):
                    d[s] += 1
            if self._point_to_index:
                return {tuple([self._points[x] for x in s]):v for s,v in d.items()}
            else:
                return d

    def rank(self):
        r"""
        Return the rank of the hypergraph (the maximum size of a block).

        EXAMPLES::

            sage: h = Hypergraph(8, [[0,1,3],[1,4,5,6],[1,2]])
            sage: h.rank()
            4
        """
        return max(len(b) for b in self._blocks)

    def is_regular(self, r=None) -> bool | int:
        r"""
        Test whether the incidence structure is `r`-regular.

        An incidence structure is said to be `r`-regular if all its points are
        incident with exactly `r` blocks.

        INPUT:

        - ``r`` -- integer

        OUTPUT:

        If ``r`` is defined, a boolean is returned. If ``r`` is set to ``None``
        (default), the method returns either ``False`` or the integer ``r`` such
        that the incidence structure is `r`-regular.

        .. WARNING::

            In case of `0`-regular incidence structure, beware that ``if not
            H.is_regular()`` is a satisfied condition.

        EXAMPLES::

            sage: designs.balanced_incomplete_block_design(7,3).is_regular()            # needs sage.schemes
            3
            sage: designs.balanced_incomplete_block_design(7,3).is_regular(r=3)         # needs sage.schemes
            True
            sage: designs.balanced_incomplete_block_design(7,3).is_regular(r=4)         # needs sage.schemes
            False

        TESTS::

            sage: IncidenceStructure([]).is_regular()
            Traceback (most recent call last):
            ...
            ValueError: This incidence structure has no points.
        """
        if self.num_points() == 0:
            raise ValueError("This incidence structure has no points.")
        count = [0] * self.num_points()
        for b in self._blocks:
            for x in b:
                count[x] += 1
        scount = set(count)
        if len(scount) != 1:
            return False
        if r is None:
            return scount.pop()
        return scount.pop() == r

    def is_uniform(self, k=None) -> bool | int:
        r"""
        Test whether the incidence structure is `k`-uniform

        An incidence structure is said to be `k`-uniform if all its blocks have
        size `k`.

        INPUT:

        - ``k`` -- integer

        OUTPUT:

        If ``k`` is defined, a boolean is returned. If ``k`` is set to ``None``
        (default), the method returns either ``False`` or the integer ``k`` such
        that the incidence structure is `k`-uniform.

        .. WARNING::

            In case of `0`-uniform incidence structure, beware that ``if not
            H.is_uniform()`` is a satisfied condition.

        EXAMPLES::

            sage: designs.balanced_incomplete_block_design(7,3).is_uniform()            # needs sage.schemes
            3
            sage: designs.balanced_incomplete_block_design(7,3).is_uniform(k=3)         # needs sage.schemes
            True
            sage: designs.balanced_incomplete_block_design(7,3).is_uniform(k=4)         # needs sage.schemes
            False

        TESTS::

            sage: IncidenceStructure([]).is_uniform()
            Traceback (most recent call last):
            ...
            ValueError: This incidence structure has no blocks.
        """
        if self.num_blocks() == 0:
            raise ValueError("This incidence structure has no blocks.")
        sizes = set(self.block_sizes())
        if len(sizes) != 1:
            return False
        if k is None:
            return sizes.pop()
        return sizes.pop() == k

    def is_connected(self) -> bool:
        r"""
        Test whether the design is connected.

        EXAMPLES::

            sage: IncidenceStructure(3, [[0,1],[0,2]]).is_connected()
            True
            sage: IncidenceStructure(4, [[0,1],[2,3]]).is_connected()
            False
        """
        from sage.sets.disjoint_set import DisjointSet
        D = DisjointSet(self.num_points())
        for B in self._blocks:
            x = B[0]
            for i in range(1, len(B)):
                D.union(x, B[i])
        return D.number_of_subsets() == 1

    def is_simple(self) -> bool:
        r"""
        Test whether this design is simple (i.e. no repeated block).

        EXAMPLES::

            sage: IncidenceStructure(3, [[0,1],[1,2],[0,2]]).is_simple()
            True
            sage: IncidenceStructure(3, [[0],[0]]).is_simple()
            False

            sage: V = [(0,'a'),(0,'b'),(1,'a'),(1,'b')]
            sage: B = [[V[0],V[1]], [V[1],V[2]]]
            sage: I = IncidenceStructure(V, B)
            sage: I.is_simple()
            True
            sage: I2 = IncidenceStructure(V, B*2)
            sage: I2.is_simple()
            False
        """
        B = self._blocks
        return all(B[i] != B[i + 1] for i in range(len(B) - 1))

    def _gap_(self):
        """
        Return the GAP string describing the design.

        EXAMPLES::

            sage: BD = IncidenceStructure(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD._gap_()
            'BlockDesign(7,[[1, 2, 3], [1, 4, 5], [1, 6, 7], [2, 4, 6], [2, 5, 7], [3, 4, 7], [3, 5, 6]])'
        """
        v = self.num_points()
        gB = [[x + 1 for x in b] for b in self._blocks]
        return "BlockDesign({},{})".format(v, gB)

    def _libgap_(self):
        """
        Return the design as a GAP record.

        EXAMPLES::

            sage: D = IncidenceStructure(4, [[0,2],[1,2,3],[2,3]])
            sage: D._libgap_()                # optional - gap_package_design
            rec( blocks := [ [ 1, 3 ], [ 2, 3, 4 ], [ 3, 4 ] ],
            isBlockDesign := true, v := 4 )
        """
        libgap.load_package("design")
        v = self.num_points()
        gB = [[x + 1 for x in b] for b in self._blocks]
        return libgap.BlockDesign(v, gB)

    def intersection_graph(self, sizes=None):
        r"""
        Return the intersection graph of the incidence structure.

        The vertices of this graph are the :meth:`blocks` of the incidence
        structure. Two of them are adjacent if the size of their intersection
        belongs to the set ``sizes``.

        INPUT:

        - ``sizes`` -- list/set of integers; for convenience, setting
          ``sizes`` to ``5`` has the same effect as ``sizes=[5]``. When set to
          ``None`` (default), behaves as ``sizes=PositiveIntegers()``.

        EXAMPLES:

        The intersection graph of a
        :func:`~sage.combinat.designs.bibd.balanced_incomplete_block_design` is
        a :meth:`strongly regular graph <Graph.is_strongly_regular>` (when it is
        not trivial)::

            sage: BIBD =  designs.balanced_incomplete_block_design(19,3)
            sage: G = BIBD.intersection_graph(1)
            sage: G.is_strongly_regular(parameters=True)
            (57, 24, 11, 9)
        """
        from sage.sets.positive_integers import PositiveIntegers
        from sage.graphs.graph import Graph
        from sage.sets.set import Set
        if sizes is None:
            sizes = PositiveIntegers()
        elif sizes in PositiveIntegers():
            sizes = (sizes,)
        V = [Set(v)  for v in self]
        return Graph([V, lambda x,y: len(x & y) in sizes], loops=False)

    def incidence_matrix(self):
        r"""
        Return the incidence matrix `A` of the design. A is a `(v \times b)`
        matrix defined by: ``A[i,j] = 1`` if ``i`` is in block ``B_j`` and 0
        otherwise.

        EXAMPLES::

            sage: BD = IncidenceStructure(7, [[0,1,2],[0,3,4],[0,5,6],[1,3,5],
            ....:                             [1,4,6],[2,3,6],[2,4,5]])
            sage: BD.block_sizes()
            [3, 3, 3, 3, 3, 3, 3]
            sage: BD.incidence_matrix()                                                 # needs sage.modules
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [1 0 0 0 0 1 1]
            [0 1 0 1 0 1 0]
            [0 1 0 0 1 0 1]
            [0 0 1 1 0 0 1]
            [0 0 1 0 1 1 0]

            sage: I = IncidenceStructure('abc', ('ab','abc','ac','c'))
            sage: I.incidence_matrix()                                                  # needs sage.modules
            [1 1 1 0]
            [1 1 0 0]
            [0 1 1 1]
        """
        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import ZZ
        A = matrix(ZZ, self.num_points(), self.num_blocks(), sparse=True)
        for j, b in enumerate(self._blocks):
            for i in b:
                A[i, j] = 1
        return A

    def incidence_graph(self,labels=False):
        r"""
        Return the incidence graph of the incidence structure.

        A point and a block are adjacent in this graph whenever they are
        incident.

        INPUT:

        - ``labels`` -- boolean; whether to return a graph whose vertices are
          integers, or labelled elements

          - ``labels is False`` -- default; in this case the first vertices
            of the graphs are the elements of :meth:`ground_set`, and appear
            in the same order. Similarly, the following vertices represent the
            elements of :meth:`blocks`, and appear in the same order.

          - ``labels is True``, the points keep their original labels, and the
            blocks are :func:`Set <Set>` objects.

            Note that the labelled incidence graph can be incorrect when
            blocks are repeated, and on some (rare) occasions when the
            elements of :meth:`ground_set` mix :func:`Set` and non-:func:`Set
            <Set>` objects.

        EXAMPLES::

            sage: BD = IncidenceStructure(7, [[0,1,2],[0,3,4],[0,5,6],[1,3,5],
            ....:                             [1,4,6],[2,3,6],[2,4,5]])
            sage: BD.incidence_graph()                                                  # needs sage.modules
            Bipartite graph on 14 vertices
            sage: A = BD.incidence_matrix()                                             # needs sage.modules
            sage: Graph(block_matrix([[A*0, A],                                         # needs sage.modules
            ....:                     [A.transpose(),A*0]])) == BD.incidence_graph()
            True

        TESTS:

        With ``labels = True``::

            sage: BD.incidence_graph(labels=True).has_edge(0,Set([0,1,2]))
            True
        """
        if labels:
            from sage.graphs.graph import Graph
            from sage.sets.set import Set
            G = Graph()
            G.add_vertices(self.ground_set())
            for b in self.blocks():
                b = Set(b)
                G.add_vertex(b)
                G.add_edges((b,x) for x in b)
            return G

        else:
            from sage.graphs.bipartite_graph import BipartiteGraph
            A = self.incidence_matrix()
            return BipartiteGraph(A)

    def is_berge_cyclic(self):
        r"""
        Check whether ``self`` is a Berge-Cyclic uniform hypergraph.

        A `k`-uniform Berge cycle (named after Claude Berge) of length `\ell`
        is a cyclic list of distinct `k`-sets `F_1,\ldots,F_\ell`, `\ell>1`,
        and distinct vertices `C = \{v_1,\ldots,v_\ell\}` such that for each
        `1\le i\le \ell`, `F_i` contains `v_i` and `v_{i+1}` (where `v_{l+1} =
        v_1`).

        A uniform hypergraph is Berge-cyclic if its incidence graph is cyclic.
        It is called "Berge-acyclic" otherwise.

        For more information, see [Fag1983]_ and :wikipedia:`Hypergraph`.

        EXAMPLES::

            sage: Hypergraph(5, [[1, 2, 3], [2, 3, 4]]).is_berge_cyclic()               # needs sage.modules
            True
            sage: Hypergraph(6, [[1, 2, 3], [3, 4, 5]]).is_berge_cyclic()               # needs sage.modules
            False

        TESTS::

            sage: Hypergraph(5, [[1, 2, 3], [2, 3]]).is_berge_cyclic()
            Traceback (most recent call last):
            ...
            TypeError: Berge cycles are defined for uniform hypergraphs only
        """
        if not self.is_uniform():
            raise TypeError("Berge cycles are defined for uniform hypergraphs only")

        return not self.incidence_graph().is_forest()

    def complement(self,uniform=False):
        r"""
        Return the complement of the incidence structure.

        Two different definitions of "complement" are made available, according
        to the value of ``uniform``.

        INPUT:

        - ``uniform`` -- boolean

          - if set to ``False`` (default), returns the incidence structure whose
            blocks are the complements of all blocks of the incidence structure.

          - If set to ``True`` and the incidence structure is `k`-uniform,
            returns the incidence structure whose blocks are all `k`-sets of the
            ground set that do not appear in ``self``.

        EXAMPLES:

        The complement of a
        :class:`~sage.combinat.designs.bibd.BalancedIncompleteBlockDesign` is
        also a `2`-design::

            sage: bibd = designs.balanced_incomplete_block_design(13,4)                 # needs sage.schemes
            sage: bibd.is_t_design(return_parameters=True)                              # needs sage.schemes
            (True, (2, 13, 4, 1))
            sage: bibd.complement().is_t_design(return_parameters=True)                 # needs sage.schemes
            (True, (2, 13, 9, 6))

        The "uniform" complement of a graph is a graph::

            sage: g = graphs.PetersenGraph()
            sage: G = IncidenceStructure(g.edges(sort=True, labels=False))
            sage: H = G.complement(uniform=True)
            sage: h = Graph(H.blocks())
            sage: g == h
            False
            sage: g == h.complement()
            True

        TESTS::

            sage: bibd.relabel({i:str(i) for i in bibd.ground_set()})                   # needs sage.schemes
            sage: bibd.complement().ground_set()                                        # needs sage.schemes
            ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']

            sage: I = IncidenceStructure('abc', ['ab','ac','bc'])
            sage: I.is_t_design(return_parameters=True)
            (True, (2, 3, 2, 1))
        """
        if uniform:
            k = self.is_uniform()
            if k is False:
                raise ValueError("The incidence structure is not uniform.")

            blocks = []
            num_blocks = self.num_blocks()
            i = 0
            from itertools import combinations
            for B in combinations(range(self.num_points()),k):
                B = list(B)
                while i < num_blocks and self._blocks[i] < B:
                    i += 1
                if i < num_blocks and self._blocks[i] == B:
                    i += 1
                    continue
                blocks.append(B)
            I = IncidenceStructure(blocks,copy=False)
        else:
            X = set(range(self.num_points()))
            I = IncidenceStructure([X.difference(B) for B in self._blocks])

        I.relabel({i:self._points[i] for i in range(self.num_points())})
        return I

    def relabel(self, perm=None, inplace=True):
        r"""
        Relabel the ground set.

        INPUT:

        - ``perm`` -- can be one of

            - a dictionary -- then each point ``p`` (which should be a key of
              ``d``) is relabeled to ``d[p]``

            - a list or a tuple of length ``n`` -- the first point returned by
              :meth:`ground_set` is relabeled to ``l[0]``, the second to
              ``l[1]``, ...

            - ``None`` -- the incidence structure is relabeled to be on
              `\{0,1,...,n-1\}` in the ordering given by :meth:`ground_set`

        - ``inplace`` -- boolean (default: ``False``); if ``True`` then return
          a relabeled graph and does not touch ``self``

        EXAMPLES::

            sage: # needs sage.schemes
            sage: TD = designs.transversal_design(5,5)
            sage: TD.relabel({i: chr(97+i) for i in range(25)})
            sage: TD.ground_set()
            ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
             'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y']
            sage: TD.blocks()[:3]
            [['a', 'f', 'k', 'p', 'u'], ['a', 'g', 'm', 's', 'y'], ['a', 'h', 'o', 'q', 'x']]

        Relabel to integer points::

            sage: TD.relabel()                                                          # needs sage.schemes
            sage: TD.blocks()[:3]                                                       # needs sage.schemes
            [[0, 5, 10, 15, 20], [0, 6, 12, 18, 24], [0, 7, 14, 16, 23]]

        TESTS:

        Check that the relabel is consistent on a fixed incidence structure::

            sage: I = IncidenceStructure([0,1,2,3,4],
            ....:               [[0,1,3],[0,2,4],[2,3,4],[0,1]])
            sage: I.relabel()
            sage: from itertools import permutations
            sage: for p in permutations([0,1,2,3,4]):
            ....:     J = I.relabel(p,inplace=False)
            ....:     if I == J: print(p)
            (0, 1, 2, 3, 4)
            (0, 1, 4, 3, 2)

        And one can also verify that we have exactly two automorphisms::

            sage: I.automorphism_group()                                                # needs sage.groups
            Permutation Group with generators [(2,4)]
        """
        if not inplace:
            from copy import copy
            G = copy(self)
            G.relabel(perm=perm, inplace=True)
            return G

        if perm is None:
            self._points = list(range(self.num_points()))
            self._point_to_index = None
            return

        if isinstance(perm, (list,tuple)):
            perm = dict(zip(self._points, perm))

        if not isinstance(perm, dict):
            raise ValueError("perm argument must be None, a list or a dictionary")

        if len(set(perm.values())) != len(perm):
            raise ValueError("two points are getting relabelled with the same name")

        self._points = [perm[x] for x in self._points]
        if self._points == list(range(self.num_points())):
            self._point_to_index = None
        else:
            self._point_to_index = {v: i for i, v in enumerate(self._points)}

    # __hash__ = None
    # This object is mutable because of .relabel()

    #####################
    # real computations #
    #####################

    def packing(self, solver=None, verbose=0, *, integrality_tolerance=1e-3):
        r"""
        Return a maximum packing.

        A maximum packing in a hypergraph is collection of disjoint sets/blocks
        of maximal cardinality. This problem is NP-complete in general, and in
        particular on 3-uniform hypergraphs. It is solved here with an Integer
        Linear Program.

        For more information, see the :wikipedia:`Packing_in_a_hypergraph`.

        INPUT:

        - ``solver`` -- (default: ``None``) specify a Mixed Integer Linear
          Programming (MILP) solver to be used. If set to ``None``, the default
          one is used. For more information on LP solvers and which default
          solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: 0); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- parameter for use with MILP solvers over
          an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        EXAMPLES::

            sage: P = IncidenceStructure([[1,2],[3,4],[2,3]]).packing()                 # needs sage.numerical.mip
            sage: sorted(sorted(b) for b in P)                                          # needs sage.numerical.mip
            [[1, 2], [3, 4]]
            sage: len(designs.steiner_triple_system(9).packing())                       # needs sage.numerical.mip
            3
        """
        from sage.numerical.mip import MixedIntegerLinearProgram

        # List of blocks containing a given point x
        d = [[] for _ in self._points]
        for i, B in enumerate(self._blocks):
            for x in B:
                d[x].append(i)

        p = MixedIntegerLinearProgram(solver=solver)
        b = p.new_variable(binary=True)
        for x, L in enumerate(d):  # Set of disjoint blocks
            p.add_constraint(p.sum([b[i] for i in L]) <= 1)

        # Maximum number of blocks
        p.set_objective(p.sum([b[i] for i in range(self.num_blocks())]))

        p.solve(log=verbose)

        values = p.get_values(b, convert=bool, tolerance=integrality_tolerance)
        return [[self._points[x] for x in self._blocks[i]]
                for i, v in values.items() if v]

    def is_t_design(self, t=None, v=None, k=None, l=None, return_parameters=False):
        r"""
        Test whether ``self`` is a `t-(v,k,l)` design.

        A `t-(v,k,\lambda)` (sometimes called `t`-design for short) is a block
        design in which:

        - the underlying set has cardinality `v`
        - the blocks have size `k`
        - each `t`-subset of points is covered by `\lambda` blocks

        INPUT:

        - ``t``, ``v``, ``k``, ``l`` -- integers; their value is set to
          ``None`` by default. The function tests whether the design is a
          `t-(v,k,l)` design using the provided values and guesses the
          others. Note that ``l`` cannot be specified if ``t`` is not.

        - ``return_parameters`` -- boolean; whether to return the parameters of
          the `t`-design. If set to ``True``, the function returns a pair
          ``(boolean_answer,(t,v,k,l))``.

        EXAMPLES::

            sage: fano_blocks = [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]]
            sage: BD = IncidenceStructure(7, fano_blocks)
            sage: BD.is_t_design()
            True
            sage: BD.is_t_design(return_parameters=True)
            (True, (2, 7, 3, 1))
            sage: BD.is_t_design(2, 7, 3, 1)
            True
            sage: BD.is_t_design(1, 7, 3, 3)
            True
            sage: BD.is_t_design(0, 7, 3, 7)
            True

            sage: BD.is_t_design(0,6,3,7) or BD.is_t_design(0,7,4,7) or BD.is_t_design(0,7,3,8)
            False

            sage: BD = designs.AffineGeometryDesign(3, 1, GF(2))                        # needs sage.combinat
            sage: BD.is_t_design(1)                                                     # needs sage.combinat
            True
            sage: BD.is_t_design(2)                                                     # needs sage.combinat
            True

        Steiner triple and quadruple systems are other names for `2-(v,3,1)` and
        `3-(v,4,1)` designs::

            sage: S3_9 = designs.steiner_triple_system(9)
            sage: S3_9.is_t_design(2,9,3,1)
            True

            sage: blocks = designs.steiner_quadruple_system(8)
            sage: S4_8 = IncidenceStructure(8, blocks)
            sage: S4_8.is_t_design(3,8,4,1)
            True

            sage: blocks = designs.steiner_quadruple_system(14)
            sage: S4_14 = IncidenceStructure(14, blocks)
            sage: S4_14.is_t_design(3,14,4,1)
            True

        Some examples of Witt designs that need the gap database::

            sage: # optional - gap_package_design
            sage: BD = designs.WittDesign(9)
            sage: BD.is_t_design(2,9,3,1)
            True
            sage: W12 = designs.WittDesign(12)
            sage: W12.is_t_design(5,12,6,1)
            True
            sage: W12.is_t_design(4)
            True

        Further examples::

            sage: D = IncidenceStructure(4,[[],[]])
            sage: D.is_t_design(return_parameters=True)
            (True,  (0, 4, 0, 2))

            sage: D = IncidenceStructure(4, [[0,1],[0,2],[0,3]])
            sage: D.is_t_design(return_parameters=True)
            (True, (0, 4, 2, 3))

            sage: D = IncidenceStructure(4, [[0],[1],[2],[3]])
            sage: D.is_t_design(return_parameters=True)
            (True, (1, 4, 1, 1))

            sage: D = IncidenceStructure(4,[[0,1],[2,3]])
            sage: D.is_t_design(return_parameters=True)
            (True, (1, 4, 2, 1))

            sage: D = IncidenceStructure(4, [list(range(4))])
            sage: D.is_t_design(return_parameters=True)
            (True, (4, 4, 4, 1))

        TESTS::

            sage: blocks = designs.steiner_quadruple_system(8)
            sage: S4_8 = IncidenceStructure(8, blocks)
            sage: R = list(range(15))
            sage: [(v,k,l) for v in R for k in R for l in R if S4_8.is_t_design(3,v,k,l)]
            [(8, 4, 1)]
            sage: [(v,k,l) for v in R for k in R for l in R if S4_8.is_t_design(2,v,k,l)]
            [(8, 4, 3)]
            sage: [(v,k,l) for v in R for k in R for l in R if S4_8.is_t_design(1,v,k,l)]
            [(8, 4, 7)]
            sage: [(v,k,l) for v in R for k in R for l in R if S4_8.is_t_design(0,v,k,l)]
            [(8, 4, 14)]

            sage: # needs sage.rings.finite_rings
            sage: A = designs.AffineGeometryDesign(3, 1, GF(2))
            sage: A.is_t_design(return_parameters=True)
            (True, (2, 8, 2, 1))
            sage: A = designs.AffineGeometryDesign(4, 2, GF(2))
            sage: A.is_t_design(return_parameters=True)
            (True, (3, 16, 4, 1))

            sage: I = IncidenceStructure(2, [])
            sage: I.is_t_design(return_parameters=True)
            (True, (0, 2, 0, 0))
            sage: I = IncidenceStructure(2, [[0],[0,1]])
            sage: I.is_t_design(return_parameters=True)
            (False, (0, 0, 0, 0))
        """
        from sage.arith.misc import binomial

        # Missing parameters ?
        if v is None:
            v = self.num_points()

        if k is None:
            k = len(self._blocks[0]) if self._blocks else 0

        if l is not None and t is None:
            raise ValueError("t must be set when l=None")

        b = self.num_blocks()

        # Trivial wrong answers
        if (any(len(block) != k for block in self._blocks) or # non k-uniform
            v != self.num_points()):
            return (False, (0,0,0,0)) if return_parameters else False

        # Trivial case t>k
        if (t is not None and t > k):
            if (l is None or l == 0):
                return (True, (t,v,k,0)) if return_parameters else True
            else:
                return (False, (0,0,0,0)) if return_parameters else False

        # Trivial case k=0
        if k == 0:
            if (l is None or l == 0):
                return (True, (0,v,k,b)) if return_parameters else True
            else:
                return (False, (0,0,0,0)) if return_parameters else False

        # Trivial case k=v (includes v=0)
        if k == v:
            if t is None:
                t = v
            if l is None or b == l:
                return (True, (t,v,k,b)) if return_parameters else True
            else:
                return (True, (0,0,0,0)) if return_parameters else False

        # Handbook of combinatorial design theorem II.4.8:
        #
        # a t-(v,k,l) is also a t'-(v,k,l')
        # for t' < t and l' = l* binomial(v-t',t-t') / binomial(k-t',t-t')
        #
        # We look for the largest t such that self is a t-design
        from itertools import combinations
        for tt in (range(1,k+1) if t is None else [t]):
            # is lambda an integer?
            if (b*binomial(k,tt)) % binomial(v,tt) != 0:
                tt -= 1
                break

            s = {}
            for block in self._blocks:
                for i in combinations(block,tt):
                    s[i] = s.get(i,0) + 1

            if len(set(s.values())) != 1:
                tt -= 1
                break

            ll = b*binomial(k,tt) // binomial(v,tt)

        if ((t is not None and t != tt) or
            (l is not None and l != ll)):
            return (False, (0,0,0,0)) if return_parameters else False
        else:
            if tt == 0:
                ll = b
            return (True, (tt,v,k,ll)) if return_parameters else True

    def is_generalized_quadrangle(self, verbose=False, parameters=False):
        r"""
        Test if the incidence structure is a generalized quadrangle.

        An incidence structure is a generalized quadrangle iff (see [BH2012]_,
        section 9.6):

        - two blocks intersect on at most one point.

        - For every point `p` not in a block `B`, there is a unique block `B'`
          intersecting both `\{p\}` and `B`

        It is a *regular* generalized quadrangle if furthermore:

        - it is `s+1`-:meth:`uniform <is_uniform>` for some positive integer `s`.

        - it is `t+1`-:meth:`regular <is_regular>` for some positive integer `t`.

        For more information, see the :wikipedia:`Generalized_quadrangle`.

        .. NOTE::

            Some references (e.g. [PT2009]_ or
            :wikipedia:`Generalized_quadrangle`) only allow *regular*
            generalized quadrangles. To use such a definition, see the
            ``parameters`` optional argument described below, or the methods
            :meth:`is_regular` and :meth:`is_uniform`.

        INPUT:

        - ``verbose`` -- boolean; whether to print an explanation when the
          instance is not a generalized quadrangle

        - ``parameters`` -- (boolean; ``False``); if set to ``True``, the
          function returns a pair ``(s,t)`` instead of ``True`` answers. In this
          case, `s` and `t` are the integers defined above if they exist (each
          can be set to ``False`` otherwise).

        EXAMPLES::

            sage: h = designs.CremonaRichmondConfiguration()                            # needs networkx
            sage: h.is_generalized_quadrangle()                                         # needs networkx
            True

        This is actually a *regular* generalized quadrangle::

            sage: h.is_generalized_quadrangle(parameters=True)                          # needs networkx
            (2, 2)

        TESTS::

            sage: H = IncidenceStructure((2*graphs.CompleteGraph(3)).edges(sort=True, labels=False))
            sage: H.is_generalized_quadrangle(verbose=True)                             # needs sage.modules
            Some point is at distance >3 from some block.
            False

            sage: G = graphs.CycleGraph(5)
            sage: B = list(G.subgraph_search_iterator(graphs.PathGraph(3),              # needs sage.modules
            ....:                                     return_graphs=False))
            sage: H = IncidenceStructure(B)                                             # needs sage.modules
            sage: H.is_generalized_quadrangle(verbose=True)                             # needs sage.modules
            Two blocks intersect on >1 points.
            False

            sage: hypergraphs.CompleteUniform(4,2).is_generalized_quadrangle(verbose=1)             # needs sage.modules
            Some point has two projections on some line.
            False
        """
        # The distance between a point and a line in the incidence graph is odd
        # and must be <= 3. Thus, the diameter is at most 4
        g = self.incidence_graph()
        if g.diameter() > 4:
            if verbose:
                print("Some point is at distance >3 from some block.")
            return False

        # There is a unique projection of a point on a line. Thus, the girth of
        # g is at least 7
        girth = g.girth()
        if girth == 4:
            if verbose:
                print("Two blocks intersect on >1 points.")
            return False
        elif girth == 6:
            if verbose:
                print("Some point has two projections on some line.")
            return False

        if parameters:
            s = self.is_uniform()
            t = self.is_regular()
            s = s-1 if (s is not False and s >= 2) else False
            t = t-1 if (t is not False and t >= 2) else False
            return (s,t)
        else:
            return True

    def dual(self, algorithm=None):
        """
        Return the dual of the incidence structure.

        INPUT:

        - ``algorithm`` -- whether to use Sage's implementation
          (``algorithm=None``, default) or use GAP's (``algorithm='gap'``)

          .. NOTE::

              The ``algorithm='gap'`` option requires GAP's Design package
              (included in the ``gap_packages`` Sage spkg).

        EXAMPLES:

        The dual of a projective plane is a projective plane::

            sage: PP = designs.DesarguesianProjectivePlaneDesign(4)                     # needs sage.rings.finite_rings
            sage: PP.dual().is_t_design(return_parameters=True)                         # needs sage.modules sage.rings.finite_rings
            (True, (2, 21, 5, 1))

        TESTS::

            sage: D = IncidenceStructure(4, [[0,2],[1,2,3],[2,3]]); D
            Incidence structure with 4 points and 3 blocks
            sage: D.dual()                                                              # needs sage.modules
            Incidence structure with 3 points and 4 blocks
            sage: print(D.dual(algorithm='gap'))            # optional - gap_package_design
            Incidence structure with 3 points and 4 blocks
            sage: blocks = [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]]
            sage: BD = IncidenceStructure(7, blocks, name='FanoPlane'); BD
            Incidence structure with 7 points and 7 blocks
            sage: print(BD.dual(algorithm='gap'))           # optional - gap_package_design
            Incidence structure with 7 points and 7 blocks
            sage: BD.dual()                                                             # needs sage.modules
            Incidence structure with 7 points and 7 blocks

        REFERENCE:

        - Leonard Soicher, :gap_package:`Design package manual <design/htm/CHAP003.htm>`
        """
        if algorithm == "gap":
            libgap.load_package("design")
            DD = libgap(self).DualBlockDesign()
            v = DD['v'].sage()
            gB = [[x - 1 for x in b] for b in DD['blocks'].sage()]
            return IncidenceStructure(list(range(v)), gB, name=None, check=False)
        else:
            return IncidenceStructure(
                          incidence_matrix=self.incidence_matrix().transpose(),
                          check=False)

    def automorphism_group(self):
        r"""
        Return the subgroup of the automorphism group of the incidence graph
        which respects the P B partition. It is (isomorphic to) the automorphism
        group of the block design, although the degrees differ.

        EXAMPLES::

            sage: # needs sage.groups sage.rings.finite_rings
            sage: P = designs.DesarguesianProjectivePlaneDesign(2); P
            (7,3,1)-Balanced Incomplete Block Design
            sage: G = P.automorphism_group()
            sage: G.is_isomorphic(PGL(3,2))
            True
            sage: G
            Permutation Group with generators [...]
            sage: G.cardinality()
            168

        A non self-dual example::

            sage: IS = IncidenceStructure(list(range(4)), [[0,1,2,3],[1,2,3]])
            sage: IS.automorphism_group().cardinality()                                 # needs sage.groups
            6
            sage: IS.dual().automorphism_group().cardinality()                          # needs sage.groups sage.modules
            1

        Examples with non-integer points::

            sage: I = IncidenceStructure('abc', ('ab','ac','bc'))
            sage: I.automorphism_group()                                                # needs sage.groups
            Permutation Group with generators [('b','c'), ('a','b')]
            sage: IncidenceStructure([[(1,2),(3,4)]]).automorphism_group()              # needs sage.groups
            Permutation Group with generators [((1,2),(3,4))]
        """
        from sage.graphs.graph import Graph
        from sage.groups.perm_gps.permgroup import PermutationGroup
        g = Graph()
        n = self.num_points()
        g.add_edges((i+n,x) for i,b in enumerate(self._blocks) for x in b)
        ag = g.automorphism_group(partition=[list(range(n)),
                                             list(range(n,n+self.num_blocks()))])

        if self._point_to_index:
            gens = [[tuple([self._points[i] for i in cycle if (not cycle or cycle[0] < n)])
                     for cycle in g.cycle_tuples()]
                    for g in ag.gens()]
        else:
            gens = [[tuple(cycle) for cycle in g.cycle_tuples() if (not cycle or cycle[0] < n)]
                    for g in ag.gens()]

        return PermutationGroup(gens, domain=self._points)

    def is_resolvable(self, certificate=False, solver=None, verbose=0, check=True,
                      *, integrality_tolerance=1e-3):
        r"""
        Test whether the hypergraph is resolvable.

        A hypergraph is said to be resolvable if its sets can be partitionned
        into classes, each of which is a partition of the ground set.

        .. NOTE::

            This problem is solved using an Integer Linear Program, and GLPK
            (the default LP solver) has been reported to be very slow on some
            instances. If you hit this wall, consider installing a more powerful
            MILP solver (CPLEX, Gurobi, ...).

        INPUT:

        - ``certificate`` -- boolean; whether to return the classes along with
          the binary answer (see examples below)

        - ``solver`` -- (default: ``None``) specify a Mixed Integer Linear
          Programming (MILP) solver to be used. If set to ``None``, the default
          one is used. For more information on MILP solvers and which default
          solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: 0); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``check`` -- boolean (default: ``True``); whether to check that
          output is correct before returning it. As this is expected to be
          useless, you may want to disable it whenever you want speed.

        - ``integrality_tolerance`` -- parameter for use with MILP solvers over
          an inexact base ring; see :meth:`MixedIntegerLinearProgram.get_values`

        EXAMPLES:

        Some resolvable designs::

            sage: TD = designs.transversal_design(2,2,resolvable=True)
            sage: TD.is_resolvable()
            True

            sage: AG = designs.AffineGeometryDesign(3,1,GF(2))                          # needs sage.combinat
            sage: AG.is_resolvable()                                                    # needs sage.combinat
            True

        Their classes::

            sage: b, cls = TD.is_resolvable(True)
            sage: b
            True
            sage: cls # random
            [[[0, 3], [1, 2]], [[1, 3], [0, 2]]]

            sage: # needs sage.combinat
            sage: b, cls = AG.is_resolvable(True)
            sage: b
            True
            sage: cls  # random
            [[[6, 7], [4, 5], [0, 1], [2, 3]],
             [[5, 7], [0, 4], [3, 6], [1, 2]],
             [[0, 2], [4, 7], [1, 3], [5, 6]],
             [[3, 4], [0, 7], [1, 5], [2, 6]],
             [[3, 7], [1, 6], [0, 5], [2, 4]],
             [[0, 6], [2, 7], [1, 4], [3, 5]],
             [[4, 6], [0, 3], [2, 5], [1, 7]]]

        A non-resolvable design::

            sage: Fano = designs.balanced_incomplete_block_design(7,3)                  # needs sage.schemes
            sage: Fano.is_resolvable()                                                  # needs sage.schemes
            False
            sage: Fano.is_resolvable(True)                                              # needs sage.schemes
            (False, [])

        TESTS::

            sage: # needs sage.combinat
            sage: _, cls1 = AG.is_resolvable(certificate=True)
            sage: _, cls2 = AG.is_resolvable(certificate=True)
            sage: cls1 is cls2
            False
        """
        if self._classes is None:
            degrees = set(self.degrees().values())
            if len(degrees) != 1:
                self._classes = False
            else:
                from sage.numerical.mip import MixedIntegerLinearProgram
                from sage.numerical.mip import MIPSolverException
                n_classes = degrees.pop()
                p = MixedIntegerLinearProgram(solver=solver)
                b = p.new_variable(binary=True)
                domain = list(range(self.num_points()))

                # Lists of blocks containing i for every i
                dual = [[] for _ in domain]
                for i,B in enumerate(self._blocks):
                    for x in B:
                        dual[x].append(i)

                # Each class is a partition
                for t in range(n_classes):
                    for x in domain:
                        p.add_constraint(p.sum(b[t,i] for i in dual[x]) == 1)

                # Each set appears exactly once
                for i in range(len(self._blocks)):
                    p.add_constraint(p.sum(b[t,i] for t in range(n_classes)) == 1)

                try:
                    p.solve(log=verbose)
                except MIPSolverException:
                    self._classes = False
                else:
                    # each class is stored as the list of indices of its blocks
                    self._classes = [[] for _ in range(n_classes)]
                    for (t, i), v in p.get_values(b, convert=bool, tolerance=integrality_tolerance).items():
                        if v:
                            self._classes[t].append(self._blocks[i])

        if check and self._classes is not False:
            assert sorted(id(c) for cls in self._classes for c in cls) == sorted(id(b) for b in self._blocks), "some set does not appear exactly once"
            domain = list(range(self.num_points()))
            for i,c in enumerate(self._classes):
                assert sorted(sum(c,[])) == domain, "class {} is not a partition".format(i)

        if self._classes is False:
            return (False, []) if certificate else False

        if certificate:
            if self._point_to_index is None:
                classes = [[block[:] for block in classs] for classs in self._classes]
            else:
                classes = [[[self._points[i] for i in block] for block in classs] for classs in self._classes]

            return (True, classes)

        else:
            return True

    def coloring(self, k=None, solver=None, verbose=0,
                 *, integrality_tolerance=1e-3):
        r"""
        Compute a (weak) `k`-coloring of the hypergraph.

        A weak coloring of a hypergraph `\mathcal H` is an assignment of colors
        to its vertices such that no set is monochromatic.

        INPUT:

        - ``k`` -- integer; compute a coloring with `k` colors if an integer is
          provided, otherwise returns an optimal coloring (i.e. with the minimum
          possible number of colors).

        - ``solver`` -- (default: ``None``) specify a Mixed Integer Linear
          Programming (MILP) solver to be used. If set to ``None``, the default
          one is used. For more information on MILP solvers and which default
          solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- nonnegative integer (default: `0`); set the level
          of verbosity you want from the linear program solver. Since the
          problem is `NP`-complete, its solving may take some time depending on
          the graph. A value of `0` means that there will be no message printed by
          the solver.

        - ``integrality_tolerance`` -- parameter for use with MILP solvers over
          an inexact base ring; see :meth:`MixedIntegerLinearProgram.get_values`

        EXAMPLES:

        The Fano plane has chromatic number 3::

            sage: len(designs.steiner_triple_system(7).coloring())                      # needs sage.numerical.mip
            3

        One admissible 3-coloring::

            sage: designs.steiner_triple_system(7).coloring()   # not tested            # needs sage.numerical.mip
            [[0, 2, 5, 1], [4, 3], [6]]

        The chromatic number of a graph is equal to the chromatic number of its
        2-uniform corresponding hypergraph::

            sage: g = graphs.PetersenGraph()
            sage: H = IncidenceStructure(g.edges(sort=True, labels=False))
            sage: len(g.coloring())
            3
            sage: len(H.coloring())                                                     # needs sage.numerical.mip
            3
        """
        if k is None:
            for k in range(self.num_points()+1):
                try:
                    return self.coloring(k)
                except ValueError:
                    pass

        if k == 0:
            if self.num_points():
                raise ValueError("Only empty hypergraphs are 0-chromatic")
            return []
        elif any(len(x) == 1 for x in self._blocks):
            raise RuntimeError("No coloring can be defined "
                               "when there is a set of size 1")
        elif k == 1:
            if any(self._blocks):
                raise ValueError("This hypergraph contains a set. "
                                 "It is not 1-chromatic")
            return [self.ground_set()]

        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
        p = MixedIntegerLinearProgram(solver=solver)
        b = p.new_variable(binary=True)

        for x in range(self.num_points()):
            p.add_constraint(p.sum(b[x,i] for i in range(k)) == 1)

        for s in self._blocks:
            for i in range(k):
                p.add_constraint(p.sum(b[x,i] for x in s) <= len(s)-1)

        try:
            p.solve(log=verbose)
        except MIPSolverException:
            raise ValueError("This hypergraph is not {}-colorable".format(k))

        col = [[] for _ in range(k)]

        for (x,i),v in p.get_values(b, convert=bool, tolerance=integrality_tolerance).items():
            if v:
                col[i].append(self._points[x])

        return col

    def edge_coloring(self):
        r"""
        Compute a proper edge-coloring.

        A proper edge-coloring is an assignment of colors to the sets of the
        incidence structure such that two sets with non-empty intersection
        receive different colors. The coloring returned minimizes the number of
        colors.

        OUTPUT: a partition of the sets into color classes

        EXAMPLES::

            sage: H = Hypergraph([{1,2,3},{2,3,4},{3,4,5},{4,5,6}]); H
            Incidence structure with 6 points and 4 blocks
            sage: C = H.edge_coloring()
            sage: C # random
            [[[3, 4, 5]], [[2, 3, 4]], [[4, 5, 6], [1, 2, 3]]]
            sage: Set(map(Set,sum(C,[]))) == Set(map(Set,H.blocks()))
            True
        """
        from sage.graphs.graph import Graph
        blocks = self.blocks()
        blocks_sets = [frozenset(b) for b in blocks]
        g = Graph([list(range(self.num_blocks())),
                   lambda x, y: len(blocks_sets[x] & blocks_sets[y])],
                  loops=False)
        return [[blocks[i] for i in C] for C in g.coloring(algorithm='MILP')]

    def _spring_layout(self):
        r"""
        Return a spring layout for the points.

        The layout is computed by creating a graph `G` on the points *and* sets
        of the incidence structure. Each set is then made adjacent in `G` with
        all points it contains before a spring layout is computed for this
        graph. The position of the points in the graph gives the position of the
        points in the final drawing.

        .. NOTE::

            This method also returns the position of the "fake" points,
            i.e. those representing the sets.

        EXAMPLES::

            sage: # needs sage.plot
            sage: H = Hypergraph([{1,2,3},{2,3,4},{3,4,5},{4,5,6}]); H
            Incidence structure with 6 points and 4 blocks
            sage: L = H._spring_layout()
            sage: L  # random
            {1: (0.238, -0.926),
             2: (0.672, -0.518),
             3: (0.449, -0.225),
             4: (0.782, 0.225),
             5: (0.558, 0.518),
             6: (0.992, 0.926),
             {3, 4, 5}: (0.504, 0.173),
             {2, 3, 4}: (0.727, -0.173),
             {4, 5, 6}: (0.838, 0.617),
             {1, 2, 3}: (0.393, -0.617)}
            sage: all(v in L for v in H.ground_set())
            True
            sage: all(v in L for v in map(Set, H.blocks()))
            True
        """
        from sage.graphs.graph import Graph

        g = Graph()
        for s in map(Set, self.blocks()):
            for x in s:
                g.add_edge((0, s), (1, x))

        _ = g.plot(iterations=50000,save_pos=True)

        # The values are rounded as TikZ does not like accuracy.
        return {k[1]: (round(x, 3), round(y, 3))
                for k, (x, y) in g.get_pos().items()}

    def _latex_(self):
        r"""
        Return a TikZ representation of the incidence structure.

        EXAMPLES::

            sage: H = Hypergraph([{1,2,3},{2,3,4},{3,4,5},{4,5,6}]); H
            Incidence structure with 6 points and 4 blocks
            sage: view(H) # not tested

        With sets of size 4::

            sage: g = graphs.Grid2dGraph(5,5)
            sage: C4 = graphs.CycleGraph(4)
            sage: sets = Set(map(Set, g.subgraph_search_iterator(C4,                    # needs sage.modules
            ....:                                                return_graphs=False)))
            sage: H = Hypergraph(sets)                                                  # needs sage.modules
            sage: view(H)                       # not tested                            # needs sage.modules sage.plot

        TESTS::

            # verify that :issue:`30976` is fixed
            sage: IS = IncidenceStructure([1,2,3], [[1,2], [2,3]])
            sage: if latex.has_file("tikz.sty"):          # optional - latex
            ....:     IS._latex_()
            ...UserWarning:
            The hypergraph is drawn as a set of closed curves...
            \begin{tikzpicture}...
            \draw... -- ...;
            \draw... -- ...;
             \draw node...;
             \draw node...;
             \draw node...;
            \end{tikzpicture}
        """
        from sage.functions.trig import arctan2

        from warnings import warn
        warn("\nThe hypergraph is drawn as a set of closed curves. The curve "
             "representing a set S goes **THROUGH** the points contained "
             "in S.\n A point which is encircled by a curve but is not located "
             "on its boundary is **NOT** included in the corresponding set.\n"
             "\n"
             "The colors are picked for readability and have no other meaning.")

        latex.add_package_to_preamble_if_available("tikz")

        if not latex.has_file("tikz.sty"):
            raise RuntimeError("You must have TikZ installed in order "
                               "to draw a hypergraph.")

        domain = self.ground_set()
        pos = self._spring_layout()
        tex = "\\begin{tikzpicture}[scale=3]\n"

        colors = ["black", "red", "green", "blue", "cyan", "magenta", "yellow","pink","brown"]
        colored_sets = [(s,i) for i,S in enumerate(self.edge_coloring()) for s in S]

        # Prints each set with its color
        for s,i in colored_sets:
            current_color = colors[i % len(colors)]

            if len(s) == 2:
                s = list(s)
                tex += ("\\draw[color="+str(current_color)+"," +
                        "line width=.1cm,opacity = .6] " +
                        str(pos[s[0]])+" -- "+str(pos[s[1]])+";\n")
                continue

            tex += ("\\draw[color="+str(current_color)+","
                    "line width=.1cm,opacity = .6,"
                    "line cap=round,"
                    "line join=round]"
                    "plot [smooth cycle,tension=1] coordinates {")

            # Reorders the vertices of s according to their angle with the
            # "center", i.e. the vertex representing the set s
            cx, cy = pos[Set(s)]
            s = [pos[_] for _ in s]
            s = sorted(s, key=lambda x_y: arctan2(x_y[0] - cx, x_y[1] - cy))

            for x in s:
                tex += str(x)+" "
            tex += "};\n"

        # Prints each vertex
        for v in domain:
            tex += "\\draw node[fill,circle,scale=.5,label={90:$"+latex(v)+"$}] at "+str(pos[v])+" {};\n"

        tex += "\\end{tikzpicture}"
        return tex

    def is_spread(self, spread):
        r"""
        Check whether the input is a spread for ``self``.

        A spread of an incidence structure `(P, B)` is a subset of `B` which
        forms a partition of `P`.

        INPUT:

        - ``spread`` -- iterable; defines the spread

        EXAMPLES::

            sage: E = IncidenceStructure([[1, 2, 3], [4, 5, 6], [1, 5, 6]])
            sage: E.is_spread([[1, 2, 3], [4, 5, 6]])
            True
            sage: E.is_spread([1, 2, 3, 4, 5, 6])
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable
            sage: E.is_spread([[1, 2, 3, 4], [5, 6]])
            False

        Order of blocks or of points within each block doesn't matter::

            sage: E = IncidenceStructure([[1, 2, 3], [4, 5, 6], [1, 5, 6]])
            sage: E.is_spread([[5, 6, 4], [3, 1, 2]])
            True

        TESTS::

            sage: E = IncidenceStructure([])
            sage: E.is_spread([])
            True
            sage: E = IncidenceStructure([[1]])
            sage: E.is_spread([])
            False
            sage: E.is_spread([[1]])
            True
            sage: E = IncidenceStructure([[1], [1]])
            sage: E.is_spread([[1]])
            True
        """

        points = set(self.ground_set())
        allBlocks = set(map(frozenset, self.blocks()))
        for block in spread:
            sblock = set(block)

            if sblock not in allBlocks:
                return False

            if not points.issuperset(sblock):
                return False

            points.difference_update(sblock)

        if points:
            return False

        return True


from sage.misc.rest_index_of_methods import gen_rest_table_index
__doc__ = __doc__.format(METHODS_OF_IncidenceStructure=gen_rest_table_index(IncidenceStructure))
