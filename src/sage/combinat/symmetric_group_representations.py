# sage.doctest: needs sage.combinat sage.modules sage.groups
r"""
Representations of the Symmetric Group

.. TODO::

    - construct the product of two irreducible representations.

    - implement Induction/Restriction of representations.

.. WARNING::

    This code uses a different convention than in Sagan's book "The Symmetric
    Group"
"""
# ****************************************************************************
#       Copyright (C) 2009 Franco Saliola <saliola@gmail.com>
#                     2024 Travis Scrimshaw <tcscrims at gmail.com>
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import lazy_import
from sage.misc.functional import sqrt
from sage.combinat.partition import Partition, Partitions
from sage.combinat.permutation import Permutation, Permutations, from_cycles
from sage.combinat.tableau import StandardTableaux, Tableau
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.finite_enumerated_set import FiniteEnumeratedSets

lazy_import("sage.combinat.yang_baxter_graph", "YangBaxterGraph_partition")
lazy_import("sage.groups.perm_gps.constructor", "PermutationGroupElement", as_='PermutationConstructor')
lazy_import("sage.symbolic.ring", "SR")


# #### Constructor function ################################################


def SymmetricGroupRepresentation(partition, implementation='specht',
        ring=None, cache_matrices=True):
    r"""
    The irreducible representation of the symmetric group corresponding to
    ``partition``.

    INPUT:

    - ``partition`` -- a partition of a positive integer

    - ``implementation`` -- string (default: ``'specht'``); one of:

      * ``'seminormal'`` -- for Young's seminormal representation
      * ``'orthogonal'`` -- for Young's orthogonal representation
      * ``'specht'`` -- for Specht's representation

    - ``ring`` -- the ring over which the representation is defined

    - ``cache_matrices`` -- boolean (default: ``True``); if ``True``, then any
      representation matrices that are computed are cached

    EXAMPLES:

    Young's orthogonal representation: the matrices are orthogonal.

    ::

        sage: orth = SymmetricGroupRepresentation([2,1], "orthogonal"); orth            # needs sage.symbolic
        Orthogonal representation of the symmetric group corresponding to [2, 1]
        sage: all(a*a.transpose() == a.parent().identity_matrix() for a in orth)        # needs sage.symbolic
        True

    ::

        sage: # needs sage.symbolic
        sage: orth = SymmetricGroupRepresentation([3,2], "orthogonal"); orth
        Orthogonal representation of the symmetric group corresponding to [3, 2]
        sage: orth([2,1,3,4,5])
        [ 1  0  0  0  0]
        [ 0  1  0  0  0]
        [ 0  0 -1  0  0]
        [ 0  0  0  1  0]
        [ 0  0  0  0 -1]
        sage: orth([1,3,2,4,5])
        [          1           0           0           0           0]
        [          0        -1/2 1/2*sqrt(3)           0           0]
        [          0 1/2*sqrt(3)         1/2           0           0]
        [          0           0           0        -1/2 1/2*sqrt(3)]
        [          0           0           0 1/2*sqrt(3)         1/2]
        sage: orth([1,2,4,3,5])
        [       -1/3 2/3*sqrt(2)           0           0           0]
        [2/3*sqrt(2)         1/3           0           0           0]
        [          0           0           1           0           0]
        [          0           0           0           1           0]
        [          0           0           0           0          -1]

    The Specht representation::

        sage: spc = SymmetricGroupRepresentation([3,2], "specht")
        sage: spc.scalar_product_matrix(Permutation([1,2,3,4,5]))
        [ 1  0  0  0  0]
        [ 0 -1  0  0  0]
        [ 0  0  1  0  0]
        [ 0  0  0  1  0]
        [-1  0  0  0 -1]
        sage: spc.scalar_product_matrix(Permutation([5,4,3,2,1]))
        [ 1 -1  0  1  0]
        [ 0  0  1  0 -1]
        [ 0  0  0 -1  1]
        [ 0  1 -1 -1  1]
        [-1  0  0  0 -1]
        sage: spc([5,4,3,2,1])
        [ 1 -1  0  1  0]
        [ 0  0 -1  0  1]
        [ 0  0  0 -1  1]
        [ 0  1 -1 -1  1]
        [ 0  1  0 -1  1]
        sage: spc.verify_representation()
        True

    By default, any representation matrices that are computed are cached::

        sage: spc = SymmetricGroupRepresentation([3,2], "specht")
        sage: spc([5,4,3,2,1])
        [ 1 -1  0  1  0]
        [ 0  0 -1  0  1]
        [ 0  0  0 -1  1]
        [ 0  1 -1 -1  1]
        [ 0  1  0 -1  1]
        sage: spc._cache__representation_matrix
        {(([5, 4, 3, 2, 1],), ()): [ 1 -1  0  1  0]
        [ 0  0 -1  0  1]
        [ 0  0  0 -1  1]
        [ 0  1 -1 -1  1]
        [ 0  1  0 -1  1]}

    This can be turned off with the keyword ``cache_matrices``::

        sage: spc = SymmetricGroupRepresentation([3,2], "specht", cache_matrices=False)
        sage: spc([5,4,3,2,1])
        [ 1 -1  0  1  0]
        [ 0  0 -1  0  1]
        [ 0  0  0 -1  1]
        [ 0  1 -1 -1  1]
        [ 0  1  0 -1  1]
        sage: hasattr(spc, '_cache__representation_matrix')
        False

    .. NOTE::

        The implementation is based on the paper [Las]_.

    REFERENCES:

    .. [Las] Alain Lascoux, 'Young representations of the symmetric group.'
       http://phalanstere.univ-mlv.fr/~al/ARTICLES/ProcCrac.ps.gz

    AUTHORS:

    - Franco Saliola (2009-04-23)
    """
    partition = Partition(partition)
    Rep = SymmetricGroupRepresentations(sum(partition), implementation=implementation,
                                        ring=ring, cache_matrices=cache_matrices)
    return Rep(partition)


def SymmetricGroupRepresentations(n, implementation='specht', ring=None,
        cache_matrices=True):
    r"""
    Irreducible representations of the symmetric group.

    INPUT:

    - ``n`` -- positive integer

    - ``implementation`` -- string (default: ``'specht'``); one of:

      * ``'seminormal'`` -- for Young's seminormal representation
      * ``'orthogonal'`` -- for Young's orthogonal representation
      * ``'specht'`` -- for Specht's representation

    - ``ring`` -- the ring over which the representation is defined

    - ``cache_matrices`` -- boolean (default: ``True``); if ``True``, then any
      representation matrices that are computed are cached

    EXAMPLES:

    Young's orthogonal representation: the matrices are orthogonal.

    ::

        sage: orth = SymmetricGroupRepresentations(3, "orthogonal"); orth               # needs sage.symbolic
        Orthogonal representations of the symmetric group of order 3! over Symbolic Ring
        sage: orth.list()                                                               # needs sage.symbolic
        [Orthogonal representation of the symmetric group corresponding to [3],
         Orthogonal representation of the symmetric group corresponding to [2, 1],
         Orthogonal representation of the symmetric group corresponding to [1, 1, 1]]
        sage: orth([2,1])([1,2,3])                                                      # needs sage.symbolic
        [1 0]
        [0 1]

    Young's seminormal representation.

    ::

        sage: snorm = SymmetricGroupRepresentations(3, "seminormal"); snorm
        Seminormal representations of the symmetric group of order 3! over Rational Field
        sage: sgn = snorm([1,1,1]); sgn
        Seminormal representation of the symmetric group corresponding to [1, 1, 1]
        sage: list(map(sgn, Permutations(3)))
        [[1], [-1], [-1], [1], [1], [-1]]

    The Specht Representation.

    ::

        sage: spc = SymmetricGroupRepresentations(5, "specht"); spc
        Specht representations of the symmetric group of order 5! over Integer Ring
        sage: spc([3,2])([5,4,3,2,1])
        [ 1 -1  0  1  0]
        [ 0  0 -1  0  1]
        [ 0  0  0 -1  1]
        [ 0  1 -1 -1  1]
        [ 0  1  0 -1  1]

    .. NOTE::

        The implementation is based on the paper [Las]_.

    AUTHORS:

    - Franco Saliola (2009-04-23)
    """
    if implementation == "seminormal":
        return YoungRepresentations_Seminormal(n, ring=ring, cache_matrices=cache_matrices)
    elif implementation == "orthogonal":
        return YoungRepresentations_Orthogonal(n, ring=ring, cache_matrices=cache_matrices)
    elif implementation == "specht":
        return SpechtRepresentations(n, ring=ring, cache_matrices=cache_matrices)
    else:
        raise NotImplementedError("only seminormal, orthogonal and specht are implemented")

# #### Generic classes for symmetric group representations #################


class SymmetricGroupRepresentation_generic_class(Element):
    r"""
    Generic methods for a representation of the symmetric group.
    """
    _default_ring = None

    def __init__(self, parent, partition):
        r"""
        An irreducible representation of the symmetric group corresponding
        to ``partition``.

        For more information, see the documentation for
        :func:`SymmetricGroupRepresentation`.

        EXAMPLES::

            sage: spc = SymmetricGroupRepresentation([3])
            sage: spc([3,2,1])
            [1]
            sage: spc == loads(dumps(spc))
            True

            sage: spc = SymmetricGroupRepresentation([3], cache_matrices=False)
            sage: spc([3,2,1])
            [1]
            sage: spc == loads(dumps(spc))
            True
        """
        self._partition = Partition(partition)
        self._n = parent._n
        self._ring = parent._ring
        if not parent._cache_matrices:
            self.representation_matrix = self._representation_matrix_uncached
        Element.__init__(self, parent)

    def __hash__(self):
        r"""
        TESTS::

            sage: spc1 = SymmetricGroupRepresentation([3], cache_matrices=True)
            sage: hash(spc1) ^^ hash((3,)) == hash(ZZ)
            True
        """
        return hash(self._ring) ^ hash(self._partition)

    def __eq__(self, other):
        r"""
        Test for equality.

        EXAMPLES::

            sage: spc1 = SymmetricGroupRepresentation([3], cache_matrices=True)
            sage: spc1([3,1,2])
            [1]
            sage: spc2 = loads(dumps(spc1))
            sage: spc1 == spc2
            True

        ::

            sage: spc3 = SymmetricGroupRepresentation([3], cache_matrices=False)
            sage: spc3([3,1,2])
            [1]
            sage: spc4 = loads(dumps(spc3))
            sage: spc3 == spc4
            True
            sage: spc1 == spc3
            True

        TESTS:

        The following tests against some bug that was fixed in :issue:`8611`::

            sage: spc = SymmetricGroupRepresentation([3])
            sage: spc.important_info = 'Sage rules'
            sage: spc == SymmetricGroupRepresentation([3])
            True
        """
        if not isinstance(other, type(other)):
            return False
        return (self._ring, self._partition) == (other._ring, other._partition)

    def __ne__(self, other):
        """
        Test for inequality.

        EXAMPLES::

            sage: spc1 = SymmetricGroupRepresentation([3], cache_matrices=True)
            sage: loads(dumps(spc1)) != spc1
            False
            sage: spc2 = SymmetricGroupRepresentation([2,1])
            sage: spc1 != spc2
            True
            sage: spc3 = SymmetricGroupRepresentation([3], cache_matrices=False)
            sage: spc1 != spc3
            False
        """
        return not (self == other)

    def __call__(self, permutation):
        r"""
        Return the image of ``permutation`` in the representation.

        EXAMPLES::

            sage: spc = SymmetricGroupRepresentation([2,1])
            sage: spc([1,3,2])
            [ 1  0]
            [ 1 -1]
        """
        return self.representation_matrix(Permutation(permutation))

    def __iter__(self):
        r"""
        Iterate over the matrices representing the elements of the
        symmetric group.

        EXAMPLES::

            sage: spc = SymmetricGroupRepresentation([1,1,1])
            sage: list(spc)
            [[1], [-1], [-1], [1], [1], [-1]]
        """
        for permutation in Permutations(self._n):
            yield self.representation_matrix(permutation)

    def verify_representation(self):
        r"""
        Verify the representation.

        This tests that the images of the simple transpositions are
        involutions and tests that the braid relations hold.

        EXAMPLES::

            sage: spc = SymmetricGroupRepresentation([1,1,1])
            sage: spc.verify_representation()
            True
            sage: spc = SymmetricGroupRepresentation([4,2,1])
            sage: spc.verify_representation()
            True
        """
        n = self._n
        transpositions = [from_cycles(n, ((i, i + 1),)) for i in range(1, n)]
        repn_matrices = [self.representation_matrix(t) for t in transpositions]
        for i, si in enumerate(repn_matrices):
            for j, sj in enumerate(repn_matrices):
                if i == j:
                    if si * sj != si.parent().identity_matrix():
                        return False, "si si != 1 for i = %s" % (i,)
                elif abs(i - j) > 1:
                    if si * sj != sj * si:
                        return False, "si sj != sj si for (i,j) =(%s,%s)" % (i, j)
                else:
                    if si * sj * si != sj * si * sj:
                        return False, "si sj si != sj si sj for (i,j) = (%s,%s)" % (i, j)
        return True

    def to_character(self):
        r"""
        Return the character of the representation.

        EXAMPLES:

        The trivial character::

            sage: rho = SymmetricGroupRepresentation([3])
            sage: chi = rho.to_character(); chi
            Character of Symmetric group of order 3! as a permutation group
            sage: chi.values()
            [1, 1, 1]
            sage: all(chi(g) == 1 for g in SymmetricGroup(3))
            True

        The sign character::

            sage: rho = SymmetricGroupRepresentation([1,1,1])
            sage: chi = rho.to_character(); chi
            Character of Symmetric group of order 3! as a permutation group
            sage: chi.values()
            [1, -1, 1]
            sage: all(chi(g) == g.sign() for g in SymmetricGroup(3))
            True

        The defining representation::

            sage: triv = SymmetricGroupRepresentation([4])
            sage: hook = SymmetricGroupRepresentation([3,1])
            sage: def_rep = lambda p : triv(p).block_sum(hook(p)).trace()
            sage: list(map(def_rep, Permutations(4)))
            [4, 2, 2, 1, 1, 2, 2, 0, 1, 0, 0, 1, 1, 0, 2, 1, 0, 0, 0, 1, 1, 2, 0, 0]
            sage: [p.to_matrix().trace() for p in Permutations(4)]
            [4, 2, 2, 1, 1, 2, 2, 0, 1, 0, 0, 1, 1, 0, 2, 1, 0, 0, 0, 1, 1, 2, 0, 0]
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        Sym = SymmetricGroup(sum(self._partition))
        values = [self(g).trace() for g in Sym.conjugacy_classes_representatives()]
        return Sym.character(values)


class SymmetricGroupRepresentations_class(UniqueRepresentation,Parent):
    r"""
    Generic methods for the CombinatorialClass of irreducible
    representations of the symmetric group.
    """

    def __init__(self, n, ring=None, cache_matrices=True):
        r"""
        Irreducible representations of the symmetric group.

        See the documentation for :func:`SymmetricGroupRepresentations`
        for more information.

        EXAMPLES::

            sage: snorm = SymmetricGroupRepresentations(3, "seminormal")
            sage: snorm == loads(dumps(snorm))
            True
        """
        self._n = n
        self._ring = ring if ring is not None else self._default_ring
        self._cache_matrices = cache_matrices
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _element_constructor_(self, partition):
        r"""
        Return the irreducible representation corresponding to ``partition``.

        EXAMPLES::

            sage: sp = SymmetricGroupRepresentations(3, "specht")
            sage: sp([1,1,1])
            Specht representation of the symmetric group corresponding to [1, 1, 1]

            sage: snorm = SymmetricGroupRepresentations(3, "seminormal")
            sage: snorm([2,1])
            Seminormal representation of the symmetric group corresponding to [2, 1]
        """
        if Partition(partition).size() != self._n:
            raise TypeError("not a partition of %s" % self._n)
        return self.element_class(self, partition)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: sp = SymmetricGroupRepresentations(4, "specht")
            sage: sp.cardinality()
            5
        """
        return Partitions(self._n).cardinality()

    def __iter__(self):
        r"""
        Iterate through all the irreducible representations of the
        symmetric group.

        EXAMPLES::

            sage: orth = SymmetricGroupRepresentations(3, "orthogonal")                 # needs sage.symbolic
            sage: for x in orth: print(x)                                               # needs sage.symbolic
            Orthogonal representation of the symmetric group corresponding to [3]
            Orthogonal representation of the symmetric group corresponding to [2, 1]
            Orthogonal representation of the symmetric group corresponding to [1, 1, 1]
        """
        for partition in Partitions(self._n):
            yield self.element_class(self, partition)

# #### Young's Seminormal Representation ###################################


class YoungRepresentation_generic(SymmetricGroupRepresentation_generic_class):
    r"""
    Generic methods for Young's representations of the symmetric group.
    """
    @lazy_attribute
    def _yang_baxter_graph(self):
        r"""
        Return the Yang-Baxter graph associated with the representation,
        with vertices labelled by the vector of contents of the partition.

        EXAMPLES::

            sage: orth = SymmetricGroupRepresentation([3,2], "orthogonal")              # needs sage.symbolic
            sage: orth._yang_baxter_graph                                               # needs sage.symbolic
            Yang-Baxter graph of [3, 2], with top vertex (0, -1, 2, 1, 0)
        """
        Y = YangBaxterGraph_partition(self._partition)
        n = self._n
        # relabel vertices with "vector of contents"
        Y.relabel_vertices(partition_to_vector_of_contents(self._partition,
                                                           reverse=True))
        # relabel edges with "differences"
        edge_relabel_dict = {}
        for u, v, op in Y.edges():
            i = op.position() + 1
            edge_relabel_dict[u, v] = (n - i, QQ((1, u[i] - u[i - 1])))
        Y.relabel_edges(edge_relabel_dict)
        return Y

    @lazy_attribute
    def _tableau_dict(self):
        r"""
        A dictionary pairing the vertices of the underlying Yang-Baxter
        graph with standard tableau.

        EXAMPLES::

            sage: orth = SymmetricGroupRepresentation([3,2], "orthogonal")              # needs sage.symbolic
            sage: orth._tableau_dict                                                    # needs sage.symbolic
            {(0, -1, 2, 1, 0): [[1, 2, 3], [4, 5]],
             (0, 2, -1, 1, 0): [[1, 2, 4], [3, 5]],
             (0, 2, 1, -1, 0): [[1, 3, 4], [2, 5]],
             (2, 0, -1, 1, 0): [[1, 2, 5], [3, 4]],
             (2, 0, 1, -1, 0): [[1, 3, 5], [2, 4]]}
        """
        # construct a dictionary pairing vertices with tableau
        t = StandardTableaux(self._partition).last()
        tableau_dict = {self._yang_baxter_graph.root(): t}
        for u, w, (i, _) in self._yang_baxter_graph._edges_in_bfs():
            # TODO: improve the following
            si = PermutationConstructor((i, i + 1))
            tableau_dict[w] = Tableau([[si(b) for b in row]
                                       for row in tableau_dict[u]])
        return tableau_dict

    @lazy_attribute
    def _word_dict(self):
        r"""
        A dictionary pairing the vertices of the underlying Yang-Baxter
        graph with words readings of standard tableau.

        EXAMPLES::

            sage: orth = SymmetricGroupRepresentation([3,2], "orthogonal")              # needs sage.symbolic
            sage: orth._word_dict                                                       # needs sage.symbolic
            {(0, -1, 2, 1, 0): (4, 5, 1, 2, 3),
             (0, 2, -1, 1, 0): (3, 5, 1, 2, 4),
             (0, 2, 1, -1, 0): (2, 5, 1, 3, 4),
             (2, 0, -1, 1, 0): (3, 4, 1, 2, 5),
             (2, 0, 1, -1, 0): (2, 4, 1, 3, 5)}
        """
        return {v: sum(reversed(t), ())
                for v, t in self._tableau_dict.items()}

    @cached_method
    def representation_matrix_for_simple_transposition(self, i):
        r"""
        Return the matrix representing the transposition that swaps ``i`` and
        ``i+1``.

        EXAMPLES::

            sage: orth = SymmetricGroupRepresentation([2,1], "orthogonal")              # needs sage.symbolic
            sage: orth.representation_matrix_for_simple_transposition(1)                # needs sage.symbolic
            [ 1  0]
            [ 0 -1]
            sage: orth.representation_matrix_for_simple_transposition(2)                # needs sage.symbolic
            [       -1/2 1/2*sqrt(3)]
            [1/2*sqrt(3)         1/2]

            sage: norm = SymmetricGroupRepresentation([2,1], "seminormal")
            sage: norm.representation_matrix_for_simple_transposition(1)
            [ 1  0]
            [ 0 -1]
            sage: norm.representation_matrix_for_simple_transposition(2)
            [-1/2  3/2]
            [ 1/2  1/2]
        """
        from copy import copy
        if not (1 <= i < sum(self._partition)):
            raise TypeError
        Y = self._yang_baxter_graph
        index_lookup = {b: a for a, b in enumerate(list(Y))}
        digraph = copy(Y._digraph)
        digraph.delete_edges((u, v) for (u, v, (j, beta)) in digraph.edges(sort=True)
                             if j != i)
        M = matrix(self._ring, digraph.num_verts())
        for g in digraph.connected_components_subgraphs():
            if g.num_verts() == 1:
                [v] = g.vertices(sort=True)
                w = self._word_dict[v]
                trivial = None
                for j, a in enumerate(w):
                    if a == i and w[j + 1] == i + 1:
                        trivial = True
                        break
                    elif a == i + 1:
                        trivial = False
                        break
                j = index_lookup[v]
                M[j, j] = 1 if trivial is True else -1
            else:
                [(u, v, (j, beta))] = g.edges(sort=True)
                iu = index_lookup[u]
                iv = index_lookup[v]
                M[iu, iu], M[iu, iv], M[iv, iu], M[iv, iv] = \
                    self._2x2_matrix_entries(beta)
        return M

    def _representation_matrix_uncached(self, permutation):
        r"""
        Return the matrix representing ``permutation``.

        EXAMPLES::

            sage: orth = SymmetricGroupRepresentation([2,1], "orthogonal")              # needs sage.symbolic
            sage: orth._representation_matrix_uncached(Permutation([2,1,3]))            # needs sage.symbolic
            [ 1  0]
            [ 0 -1]
            sage: orth._representation_matrix_uncached(Permutation([1,3,2]))            # needs sage.symbolic
            [       -1/2 1/2*sqrt(3)]
            [1/2*sqrt(3)         1/2]

        ::

            sage: norm = SymmetricGroupRepresentation([2,1], "seminormal")
            sage: p = PermutationGroupElement([2,1,3])
            sage: norm._representation_matrix_uncached(p)
            [ 1  0]
            [ 0 -1]
            sage: p = PermutationGroupElement([1,3,2])
            sage: norm._representation_matrix_uncached(p)
            [-1/2  3/2]
            [ 1/2  1/2]
        """
        m = self._yang_baxter_graph._digraph.num_verts()
        M = matrix(self._ring, m, m, 1)
        for i in Permutation(permutation).reduced_word():
            M *= self.representation_matrix_for_simple_transposition(i)
        return M

    @cached_method
    def representation_matrix(self, permutation):
        r"""
        Return the matrix representing ``permutation``.

        EXAMPLES::

            sage: orth = SymmetricGroupRepresentation([2,1], "orthogonal")              # needs sage.symbolic
            sage: orth.representation_matrix(Permutation([2,1,3]))                      # needs sage.symbolic
            [ 1  0]
            [ 0 -1]
            sage: orth.representation_matrix(Permutation([1,3,2]))                      # needs sage.symbolic
            [       -1/2 1/2*sqrt(3)]
            [1/2*sqrt(3)         1/2]

        ::

            sage: norm = SymmetricGroupRepresentation([2,1], "seminormal")
            sage: p = PermutationGroupElement([2,1,3])
            sage: norm.representation_matrix(p)
            [ 1  0]
            [ 0 -1]
            sage: p = PermutationGroupElement([1,3,2])
            sage: norm.representation_matrix(p)
            [-1/2  3/2]
            [ 1/2  1/2]
        """
        ret = self._representation_matrix_uncached(permutation)
        ret.set_immutable()
        return ret


class YoungRepresentation_Seminormal(YoungRepresentation_generic):
    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: SymmetricGroupRepresentation([2,1], "seminormal")
            Seminormal representation of the symmetric group corresponding to [2, 1]
        """
        return "Seminormal representation of the symmetric group corresponding to {}".format(self._partition)

    def _2x2_matrix_entries(self, beta):
        r"""
        Young's representations are constructed by combining
        `2 \times 2`-matrices that depend on ``beta``.

        For the seminormal representation, this is the following matrix::

            [  -beta     1+beta ]
            [ 1-beta      beta  ]

        EXAMPLES::

            sage: snorm = SymmetricGroupRepresentation([2,1], "seminormal")
            sage: snorm._2x2_matrix_entries(1/2)
            (-1/2, 3/2, 1/2, 1/2)
        """
        return (-beta, 1 + beta, 1 - beta, beta)


class YoungRepresentations_Seminormal(SymmetricGroupRepresentations_class):
    _default_ring = QQ

    Element = YoungRepresentation_Seminormal

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.symmetric_group_representations import YoungRepresentations_Seminormal
            sage: YoungRepresentations_Seminormal(3)
            Seminormal representations of the symmetric group of order 3! over Rational Field
        """
        return "Seminormal representations of the symmetric group of order %s! over %s" % (self._n, self._ring)

# #### Young's Orthogonal Representation ###################################


class YoungRepresentation_Orthogonal(YoungRepresentation_generic):
    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: SymmetricGroupRepresentation([2,1], "orthogonal")                     # needs sage.symbolic
            Orthogonal representation of the symmetric group corresponding to [2, 1]
        """
        return "Orthogonal representation of the symmetric group corresponding to {}".format(self._partition)

    def _2x2_matrix_entries(self, beta):
        r"""
        Young's representations are constructed by combining
        `2 \times 2`-matrices that depend on ``beta``.

        For the orthogonal representation, this is the following matrix::

            [     -beta       sqrt(1-beta^2) ]
            [ sqrt(1-beta^2)       beta      ]

        EXAMPLES::

            sage: orth = SymmetricGroupRepresentation([2,1], "orthogonal")              # needs sage.symbolic
            sage: orth._2x2_matrix_entries(1/2)                                         # needs sage.symbolic
            (-1/2, 1/2*sqrt(3), 1/2*sqrt(3), 1/2)
        """
        return (-beta, sqrt(1 - beta**2), sqrt(1 - beta**2), beta)


class YoungRepresentations_Orthogonal(SymmetricGroupRepresentations_class):
    _default_ring = SR

    Element = YoungRepresentation_Orthogonal

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.symmetric_group_representations import YoungRepresentations_Orthogonal
            sage: YoungRepresentations_Orthogonal(3)                                    # needs sage.symbolic
            Orthogonal representations of the symmetric group of order 3! over Symbolic Ring
        """
        return "Orthogonal representations of the symmetric group of order %s! over %s" % (self._n, self._ring)

# #### Specht Representation ###############################################


class SpechtRepresentation(SymmetricGroupRepresentation_generic_class):
    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: SymmetricGroupRepresentation([2,1], "specht")
            Specht representation of the symmetric group corresponding to [2, 1]
        """
        return "Specht representation of the symmetric group corresponding to {}".format(self._partition)

    _default_ring = ZZ

    @lazy_attribute
    def _yang_baxter_graph(self):
        r"""
        Construct and cache the underlying Yang-Baxter graph.

        EXAMPLES::

            sage: rho = SymmetricGroupRepresentation([3,2], 'specht')
            sage: rho._yang_baxter_graph
            Yang-Baxter graph of [3, 2], with top vertex (1, 0, 2, 1, 0)
        """
        return YangBaxterGraph_partition(self._partition)

    @lazy_attribute
    def _dual_vertices(self):
        r"""
        Return a list of the dual vertices of the vertices of the underlying
        Yang-Baxter graph.

        EXAMPLES::

            sage: spc = SymmetricGroupRepresentation([3,2], 'specht')
            sage: spc._dual_vertices
            [(3, 3, 0, 0, 0), (3, 0, 3, 0, 0), (3, 0, 0, 3, 0), (0, 3, 3, 0, 0), (0, 3, 0, 3, 0)]
        """
        top = self._yang_baxter_graph.root()
        exponents = tuple(i - x for i, x in enumerate(reversed(top)))[::-1]
        relabelling = self._yang_baxter_graph.vertex_relabelling_dict(exponents)
        return [relabelling[u] for u in self._yang_baxter_graph]

    @cached_method
    def scalar_product(self, u, v):
        r"""
        Return ``0`` if ``u+v`` is not a permutation, and the signature of the
        permutation otherwise.

        This is the scalar product of a vertex ``u`` of the underlying
        Yang-Baxter graph with the vertex ``v`` in the 'dual' Yang-Baxter
        graph.

        EXAMPLES::

            sage: spc = SymmetricGroupRepresentation([3,2], 'specht')
            sage: spc.scalar_product((1,0,2,1,0),(0,3,0,3,0))
            -1
            sage: spc.scalar_product((1,0,2,1,0),(3,0,0,3,0))
            0
        """
        uv = [a + v[i] + 1 for i, a in enumerate(u)]
        if uv not in Permutations():
            return 0
        else:
            return Permutation(uv).signature()

    def scalar_product_matrix(self, permutation=None):
        r"""
        Return the scalar product matrix corresponding to ``permutation``.

        The entries are given by the scalar products of ``u`` and
        ``permutation.action(v)``, where ``u`` is a vertex in the underlying
        Yang-Baxter graph and ``v`` is a vertex in the dual graph.

        EXAMPLES::

            sage: spc = SymmetricGroupRepresentation([3,1], 'specht')
            sage: spc.scalar_product_matrix()
            [ 1  0  0]
            [ 0 -1  0]
            [ 0  0  1]
        """
        if permutation is None:
            permutation = Permutation(range(1, 1 + self._n))
        Q = matrix(QQ, len(self._yang_baxter_graph))
        for i, v in enumerate(self._dual_vertices):
            for j, u in enumerate(self._yang_baxter_graph):
                Q[i, j] = self.scalar_product(tuple(permutation.action(v)), u)
        return Q

    @lazy_attribute
    def _scalar_product_matrix_inverse(self):
        r"""
        Compute and store the inverse of the scalar product matrix.

        EXAMPLES::

            sage: spc = SymmetricGroupRepresentation([3,1], 'specht')
            sage: spc._scalar_product_matrix_inverse
            [ 1  0  0]
            [ 0 -1  0]
            [ 0  0  1]
        """
        return self.scalar_product_matrix().inverse()

    @cached_method
    def representation_matrix(self, permutation):
        r"""
        Return the matrix representing the ``permutation`` in this
        irreducible representation.

        .. NOTE::

            This method caches the results.

        EXAMPLES::

            sage: spc = SymmetricGroupRepresentation([3,1], 'specht')
            sage: spc.representation_matrix(Permutation([2,1,3,4]))
            [ 0 -1  0]
            [-1  0  0]
            [ 0  0  1]
            sage: spc.representation_matrix(Permutation([3,2,1,4]))
            [0 0 1]
            [0 1 0]
            [1 0 0]
        """
        ret = self._representation_matrix_uncached(permutation)
        ret.set_immutable()
        return ret

    def _representation_matrix_uncached(self, permutation):
        r"""
        Return the matrix representing the ``permutation`` in this
        irreducible representation.

        EXAMPLES::

            sage: spc = SymmetricGroupRepresentation([3,1], 'specht')
            sage: spc._representation_matrix_uncached(Permutation([2,1,3,4]))
            [ 0 -1  0]
            [-1  0  0]
            [ 0  0  1]
            sage: spc._representation_matrix_uncached(Permutation([3,2,1,4]))
            [0 0 1]
            [0 1 0]
            [1 0 0]
        """
        R = self.scalar_product_matrix(permutation)
        return self._scalar_product_matrix_inverse * R


class SpechtRepresentations(SymmetricGroupRepresentations_class):
    _default_ring = ZZ

    Element = SpechtRepresentation

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: spc = SymmetricGroupRepresentations(4)
            sage: spc
            Specht representations of the symmetric group of order 4! over Integer Ring
        """
        return "Specht representations of the symmetric group of order %s! over %s" % (self._n, self._ring)

# ##### Miscellaneous functions ############################################


def partition_to_vector_of_contents(partition, reverse=False):
    r"""
    Return the "vector of contents" associated to ``partition``.

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_representations import partition_to_vector_of_contents
        sage: partition_to_vector_of_contents([3,2])
        (0, 1, 2, -1, 0)
    """
    v = []
    for i, p in enumerate(partition):
        v.extend(range(-i, -i + p))
    if reverse:
        return tuple(v)[::-1]
    return tuple(v)


# #### Garsia-Procesi modules ################################################

from sage.rings.quotient_ring import QuotientRing_generic
from sage.combinat.specht_module import SymmetricGroupRepresentation as SymmetricGroupRepresentation_mixin
class GarsiaProcesiModule(UniqueRepresentation, QuotientRing_generic, SymmetricGroupRepresentation_mixin):
    r"""
    A Garsia-Procesi module.

    Let `\lambda` be a partition of `n` and `R` be a commutative
    ring. The *Garsia-Procesi module* is defined by `R_{\lambda}
    := R[x_1, \ldots, x_n] / I_{\lambda}`, where

    .. MATH::

        I_{\lambda} := \langle e_r(x_{i_1}, \ldots, x_{i_k}) \mid
        \{i_1, \ldots, i_k\} \subseteq [n] \text{ and }
        k \geq r > k - d_k(\lambda) \rangle,

    with `e_r` being the `r`-the elementary symmetric function and
    `d_k(\lambda) = \lambda'_n + \cdots + \lambda'_{n+1-k}`, is the
    *Tanisaki ideal*.

    If we consider `R = \QQ`, then the Garsia-Procesi module has the
    following interpretation. Let `\mathcal{F}_n = GL_n / B` denote
    the (complex type A) flag variety. Consider the Springer fiber
    `F_{\lambda} \subseteq \mathcal{F}_n` associated to a nilpotent
    matrix with Jordan blocks sizes `\lambda`. Springer showed that
    the cohomology ring `H^*(F_{\lambda})` admits a graded `S_n`-action
    that agrees with the induced representation of the sign representation
    of the Young subgroup `S_{\lambda}`. From work of De Concini
    and Procesi, this `S_n`-representation is isomorphic to `R_{\lambda}`.
    Moreover, the graded Frobenius image is known to be a modified
    Hall-Littlewood polynomial.

    EXAMPLES::

        sage: SGA = SymmetricGroupAlgebra(QQ, 7)
        sage: GP421 = SGA.garsia_procesi_module([4, 2, 1])
        sage: GP421.dimension()
        105
        sage: v = GP421.an_element(); v
        -gp1 - gp2 - gp3 - gp4 - gp5 - gp6
        sage: SGA.an_element() * v
        -6*gp1 - 6*gp2 - 6*gp3 - 6*gp4 - 6*gp5 - 5*gp6

    We verify the result is a modified Hall-Littlewood polynomial by using
    the `Q'` Hall-Littlewood polynomials, replacing `q \mapsto q^{-1}` and
    multiplying by the smallest power of `q` so the coefficients are again
    polynomials::

        sage: GP421.graded_frobenius_image()
        q^4*s[4, 2, 1] + q^3*s[4, 3] + q^3*s[5, 1, 1] + (q^3+q^2)*s[5, 2]
         + (q^2+q)*s[6, 1] + s[7]
        sage: R.<q> = QQ[]
        sage: Sym = SymmetricFunctions(R)
        sage: s = Sym.s()
        sage: Qp = Sym.hall_littlewood(q).Qp()
        sage: mHL = s(Qp[4,2,1]); mHL
        s[4, 2, 1] + q*s[4, 3] + q*s[5, 1, 1] + (q^2+q)*s[5, 2]
         + (q^3+q^2)*s[6, 1] + q^4*s[7]
        sage: mHL.map_coefficients(lambda c: R(q^4*c(q^-1)))
        q^4*s[4, 2, 1] + q^3*s[4, 3] + q^3*s[5, 1, 1] + (q^3+q^2)*s[5, 2]
         + (q^2+q)*s[6, 1] + s[7]

    We show that the maximal degree component corresponds to the Yamanouchi
    words of content `\lambda`::

        sage: B = GP421.graded_decomposition(4).basis()
        sage: top_deg = [Word([i+1 for i in b.lift().lift().exponents()[0]]) for b in B]
        sage: yamanouchi = [P.to_packed_word() for P in OrderedSetPartitions(range(7), [4, 2, 1])
        ....:               if P.to_packed_word().reversal().is_yamanouchi()]
        sage: set(top_deg) == set(yamanouchi)
        True
    """
    @staticmethod
    def __classcall_private__(cls, SGA, shape):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: from sage.combinat.symmetric_group_representations import GarsiaProcesiModule
            sage: GP1 = GarsiaProcesiModule(SGA, [2, 2, 1])
            sage: GP2 = GarsiaProcesiModule(SGA, Partitions(5)([2, 2, 1]))
            sage: GP1 is GP2
            True
            sage: GarsiaProcesiModule(SGA, [3])
            Traceback (most recent call last):
            ...
            ValueError: [3] is not a partition of 5
        """
        shape = Partition(shape)
        if sum(shape) != SGA.n:
            raise ValueError(f"{shape} is not a partition of {SGA.n}")
        return super().__classcall__(cls, SGA, shape)

    def __init__(self, SGA, shape):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: GP = SGA.garsia_procesi_module([2, 2])
            sage: TestSuite(GP).run()

            sage: SGA = SymmetricGroupAlgebra(GF(2), 5)
            sage: GP = SGA.garsia_procesi_module([3, 1, 1])
            sage: TestSuite(GP).run()
        """
        self._shape = shape
        SymmetricGroupRepresentation_mixin.__init__(self, SGA)

        # Construct the Tanisaki ideal
        from sage.combinat.sf.sf import SymmetricFunctions
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from itertools import combinations
        n = SGA.n

        conj = list(shape.conjugate())
        conj += [0]*(n - len(conj))

        def p(k):
            return sum(conj[i] for i in range(n-k, n))

        BR = SGA.base_ring()
        R = PolynomialRing(BR, 'x', n)
        gens = R.gens()
        e = SymmetricFunctions(BR).e()
        I = R.ideal([e[d].expand(k)(*S)
                     for k in range(n+1) for d in range(k-p(k)+1, k+1)
                     for S in combinations(gens, k)])

        # Finalize the initialization
        names = tuple([f"gp{i}" for i in range(n)])
        from sage.categories.commutative_rings import CommutativeRings
        from sage.categories.algebras import Algebras
        cat = CommutativeRings().Quotients() & Algebras(SGA.base_ring()).Graded().WithBasis().FiniteDimensional()
        QuotientRing_generic.__init__(self, R, I, names=names, category=cat)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: SGA.garsia_procesi_module([2, 2])
            Garsia-Procesi module of shape [2, 2] over Rational Field
        """
        return "Garsia-Procesi module of shape {} over {}".format(self._shape, self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: GP = SGA.garsia_procesi_module([2, 2])
            sage: latex(GP)
            R_{{\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
            \lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-2}
            \lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-2}
            \end{array}$}
            }}^{\Bold{Q}}
        """
        from sage.misc.latex import latex
        return "R_{{{}}}^{{{}}}".format(latex(self._shape), latex(self.base_ring()))

    def _coerce_map_from_base_ring(self):
        r"""
        Disable the coercion from the base ring from the category.

        TESTS::

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: GP = SGA.garsia_procesi_module([2, 2])
            sage: GP._coerce_map_from_base_ring() is None
            True
        """
        return None  # don't need anything special

    @cached_method
    def get_order(self):
        """
        Return the order of the elements in the basis.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: GP = SGA.garsia_procesi_module([2, 2])
            sage: GP.get_order()
            (0, 1, 2, 3, 4, 5)
        """
        return tuple(self.basis().keys())

    @cached_method
    def basis(self):
        r"""
        Return a basis of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: GP = SGA.garsia_procesi_module([2, 2])
            sage: GP.basis()
            Family (gp2*gp3, gp1*gp3, gp3, gp2, gp1, 1)
        """
        from sage.sets.family import Family
        B = self.defining_ideal().normal_basis()
        return Family([self.retract(b) for b in B])

    @cached_method
    def one_basis(self):
        r"""
        Return the index of the basis element `1`.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: GP = SGA.garsia_procesi_module([2, 2])
            sage: GP.one_basis()
            5
        """
        B = self.defining_ideal().normal_basis()
        for i, b in enumerate(B):
            if b.is_one():
                return ZZ(i)

    @cached_method
    def dimension(self):
        r"""
        Return the dimension of ``self``.

        The graded Frobenius character of the Garsia-Procesi module
        `R_{\lambda}` is given by the modified Hall-Littlewood polynomial
        `\widetilde{H}_{\lambda'}(x; q)`.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: Sym = SymmetricFunctions(QQ)
            sage: s = Sym.s()
            sage: Qp = Sym.hall_littlewood(1).Qp()
            sage: for la in Partitions(5):
            ....:     print(SGA.garsia_procesi_module(la).dimension(),
            ....:           sum(c * StandardTableaux(la).cardinality()
            ....:               for la, c in s(Qp[la])))
            1 1
            5 5
            10 10
            20 20
            30 30
            60 60
            120 120
        """
        return self.defining_ideal().vector_space_dimension()

    @cached_method
    def graded_frobenius_image(self):
        r"""
        Return the graded Frobenius image of ``self``.

        The graded Frobenius image is the sum of the :meth:`frobenius_image`
        of each graded component, which is known to result in the modified
        Hall-Littlewood polynomial `\widetilde{H}_{\lambda}(x; q)`.

        EXAMPLES:

        We verify that the result is the modified Hall-Littlewood polynomial
        for `n = 5`::

            sage: R.<q> = QQ[]
            sage: Sym = SymmetricFunctions(R)
            sage: s = Sym.s()
            sage: Qp = Sym.hall_littlewood(q).Qp()
            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: for la in Partitions(5):
            ....:     f = SGA.garsia_procesi_module(la).graded_frobenius_image()
            ....:     d = f[la].degree()
            ....:     assert f.map_coefficients(lambda c: R(c(~q)*q^d)) == s(Qp[la])
        """
        from sage.combinat.sf.sf import SymmetricFunctions
        R = QQ['q']
        q = R.gen()
        Sym = SymmetricFunctions(R)
        p = Sym.p()
        s = Sym.s()
        G = self._semigroup
        CCR = [(elt, elt.cycle_type()) for elt in G.conjugacy_classes_representatives()]
        B = self.basis()
        return s(p._from_dict({la: coeff / la.centralizer_size() for elt, la in CCR
                               if (coeff := sum(q**b.degree() * (elt * b).lift().monomial_coefficient(b.lift())
                                                for b in B))},
                              remove_zeros=False))

    @cached_method
    def graded_character(self):
        r"""
        Return the graded character of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: GP = SGA.garsia_procesi_module([2, 2, 1])
            sage: gchi = GP.graded_character(); gchi
            (5*q^4 + 11*q^3 + 9*q^2 + 4*q + 1, -q^4 + q^3 + 3*q^2 + 2*q + 1,
             q^4 - q^3 + q^2 + 1, -q^4 - q^3 + q + 1, -q^4 + q^3 - q + 1,
             q^4 - q^3 - q^2 + 1, q^3 - q^2 - q + 1)
            sage: R.<q> = QQ[]
            sage: gchi == sum(q^d * D.character()
            ....:             for d, D in GP.graded_decomposition().items())
            True
        """
        q = QQ['q'].gen()
        G = self._semigroup
        B = self.basis()
        from sage.modules.free_module_element import vector
        return vector([sum(q**b.degree() * (g * b).lift().monomial_coefficient(b.lift()) for b in B)
                       for g in G.conjugacy_classes_representatives()],
                      immutable=True)

    @lazy_attribute
    def _graded_decomposition(self):
        """
        Construct the (internal) dictionary that encodes the graded
        decomposition of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(2), 5)
            sage: GP32 = SGA.garsia_procesi_module([3, 2])
            sage: GP32._graded_decomposition
            {0: Subrepresentation with basis {0} of Garsia-Procesi ...,
             1: Subrepresentation with basis {0, 1, 2, 3} of Garsia-Procesi ...,
             2: Subrepresentation with basis {0, 1, 2, 3, 4} of Garsia-Procesi ...}
        """
        d = {}
        for b in self.basis():
            deg = b.degree()
            if deg not in d:
                d[deg] = [b]
            else:
                d[deg].append(b)
        return {deg: self.subrepresentation(gens, is_closed=True)
                for deg, gens in sorted(d.items())}

    def graded_decomposition(self, k=None):
        r"""
        Return the decomposition of ``self`` as a direct sum of
        representations given by a fixed grading.

        INPUT:

        - ``k`` -- (optional) integer; if given, return the `k`-th graded part

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(2), 5)
            sage: GP32 = SGA.garsia_procesi_module([3, 2])
            sage: decomp = GP32.graded_decomposition(); decomp
            {0: Subrepresentation with basis {0} of Garsia-Procesi ...,
             1: Subrepresentation with basis {0, 1, 2, 3} of Garsia-Procesi ...,
             2: Subrepresentation with basis {0, 1, 2, 3, 4} of Garsia-Procesi ...}
            sage: decomp[2] is GP32.graded_decomposition(2)
            True
            sage: GP32.graded_decomposition(10)
            Subrepresentation with basis {} of Garsia-Procesi module
             of shape [3, 2] over Finite Field of size 2
        """
        if k is None:
            # make a copy since mutable
            return dict(self._graded_decomposition)
        if k < 0 or k not in self._graded_decomposition:
            return self.subrepresentation([], is_closed=True)
        return self._graded_decomposition[k]

    def graded_representation_matrix(self, elt, q=None):
        r"""
        Return the matrix corresponding to the left action of the symmetric
        group (algebra) element ``elt`` on ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(3), 3)
            sage: GP = SGA.garsia_procesi_module([1, 1, 1])
            sage: elt = SGA.an_element(); elt
            [1, 2, 3] + 2*[1, 3, 2] + [3, 1, 2]
            sage: X = GP.graded_representation_matrix(elt); X
            [  0   0   0   0   0   0]
            [  0 q^2   0   0   0   0]
            [  0 q^2   0   0   0   0]
            [  0   0   0   q   0   0]
            [  0   0   0   q   0   0]
            [  0   0   0   0   0   1]
            sage: X.parent()
            Full MatrixSpace of 6 by 6 dense matrices over
             Univariate Polynomial Ring in q over Finite Field of size 3
            sage: R.<q> = GF(3)[]
            sage: t = R.quotient([q^2+2*q+1]).gen()
            sage: GP.graded_representation_matrix(elt, t)
            [       0        0        0        0        0        0]
            [       0 qbar + 2        0        0        0        0]
            [       0 qbar + 2        0        0        0        0]
            [       0        0        0     qbar        0        0]
            [       0        0        0     qbar        0        0]
            [       0        0        0        0        0        1]
        """
        if q is None:
            q = self.base_ring()['q'].gen()
        R = q.parent()
        return matrix(R, [q**b.degree() * (elt * b).to_vector().change_ring(R)
                          for b in self.basis()])

    def graded_brauer_character(self):
        r"""
        Return the graded Brauer character of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(2), 5)
            sage: GP311 = SGA.garsia_procesi_module([3, 1, 1])
            sage: GP311.graded_brauer_character()
            (6*q^3 + 9*q^2 + 4*q + 1, q + 1, q^3 - q^2 - q + 1)
        """
        q = QQ['q'].gen()
        return sum(q**d * SM.brauer_character() for d, SM in self._graded_decomposition.items())

    class Element(QuotientRing_generic.Element):
        def _acted_upon_(self, scalar, self_on_left=True):
            r"""
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: SGA = SymmetricGroupAlgebra(GF(3), 4)
                sage: GP22 = SGA.garsia_procesi_module([2, 2])
                sage: x = SGA.an_element(); x
                [1, 2, 3, 4] + 2*[1, 2, 4, 3] + [4, 1, 2, 3]
                sage: v = GP22.an_element(); v
                -gp1 - gp2 - gp3
                sage: g = SGA.group().an_element(); g
                [4, 1, 2, 3]
                sage: g * v  # indirect doctest
                gp3
                sage: x * v  # indirect doctest
                gp3
                sage: 2 * v  # indirect doctest
                gp1 + gp2 + gp3
            """
            P = self.parent()
            if scalar in P.base_ring():
                return super()._acted_upon_(scalar, self_on_left)
            if scalar in P._semigroup:
                gens = P.ambient().gens()
                return P.retract(self.lift().subs({g: gens[scalar(i+1)-1] for i, g in enumerate(gens)}))
            if not self_on_left and scalar in P._semigroup_algebra:
                scalar = P._semigroup_algebra(scalar)
                gens = P.ambient().gens()
                return P.sum(c * P.retract(self.lift().subs({g: gens[sigma(i+1)-1] for i, g in enumerate(gens)}))
                             for sigma, c in scalar.monomial_coefficients(copy=False).items())
            return super()._acted_upon_(scalar, self_on_left)

        def to_vector(self, order=None):
            r"""
            Return ``self`` as a (dense) free module vector.

            EXAMPLES::

                sage: SGA = SymmetricGroupAlgebra(GF(3), 4)
                sage: GP22 = SGA.garsia_procesi_module([2, 2])
                sage: v = GP22.an_element(); v
                -gp1 - gp2 - gp3
                sage: v.to_vector()
                (0, 0, 2, 2, 2, 0)
            """
            P = self.parent()
            B = P.basis()
            FM = P._dense_free_module()
            f = self.lift()
            return FM([f.monomial_coefficient(b.lift()) for b in B])

        _vector_ = to_vector

        def monomial_coefficients(self, copy=None):
            r"""
            Return the monomial coefficients of ``self``.

            EXAMPLES::

                sage: SGA = SymmetricGroupAlgebra(GF(3), 4)
                sage: GP31 = SGA.garsia_procesi_module([3, 1])
                sage: v = GP31.an_element(); v
                -gp1 - gp2 - gp3
                sage: v.monomial_coefficients()
                {0: 2, 1: 2, 2: 2, 3: 0}
            """
            B = self.parent().basis()
            f = self.lift()
            return {i: f.monomial_coefficient(b.lift()) for i, b in enumerate(B)}

        def degree(self):
            r"""
            Return the degree of ``self``.

            EXAMPLES::

                sage: SGA = SymmetricGroupAlgebra(GF(3), 4)
                sage: GP22 = SGA.garsia_procesi_module([2, 2])
                sage: for b in GP22.basis():
                ....:     print(b, b.degree())
                gp2*gp3 2
                gp1*gp3 2
                gp3 1
                gp2 1
                gp1 1
                1 0
                sage: v = sum(GP22.basis())
                sage: v.degree()
                2
            """
            return self.lift().degree()

        def homogeneous_degree(self):
            r"""
            Return the (homogeneous) degree of ``self`` if homogeneous
            otherwise raise an error.

            EXAMPLES::

                sage: SGA = SymmetricGroupAlgebra(GF(2), 4)
                sage: GP31 = SGA.garsia_procesi_module([3, 1])
                sage: for b in GP31.basis():
                ....:     print(b, b.homogeneous_degree())
                gp3 1
                gp2 1
                gp1 1
                1 0
                sage: v = sum(GP31.basis()); v
                gp1 + gp2 + gp3 + 1
                sage: v.homogeneous_degree()
                Traceback (most recent call last):
                ...
                ValueError: element is not homogeneous

            TESTS::

                sage: SGA = SymmetricGroupAlgebra(GF(3), 4)
                sage: GP4 = SGA.garsia_procesi_module([4])
                sage: GP4.zero().homogeneous_degree()
                Traceback (most recent call last):
                ...
                ValueError: the zero element does not have a well-defined degree
            """
            if not self:
                raise ValueError("the zero element does not have a well-defined degree")
            f = self.lift()
            if not f.is_homogeneous():
                raise ValueError("element is not homogeneous")
            return f.degree()
