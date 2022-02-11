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

from sage.symbolic.ring import SR
from sage.misc.functional import sqrt
from sage.combinat.partition import Partition, Partitions
from sage.combinat.permutation import Permutation, Permutations, from_cycles
from sage.combinat.tableau import StandardTableaux, Tableau
from sage.combinat.yang_baxter_graph import YangBaxterGraph_partition
from sage.groups.perm_gps.constructor import PermutationGroupElement as PermutationConstructor
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.finite_enumerated_set import FiniteEnumeratedSets

# #### Constructor function ################################################


def SymmetricGroupRepresentation(partition, implementation="specht",
        ring=None, cache_matrices=True):
    r"""
    The irreducible representation of the symmetric group corresponding to
    ``partition``.

    INPUT:

    - ``partition`` -- a partition of a positive integer

    - ``implementation`` -- string (default: ``"specht"``), one of:

      * ``"seminormal"`` - for Young's seminormal representation
      * ``"orthogonal"`` - for Young's orthogonal representation
      * ``"specht"`` - for Specht's representation

    - ``ring`` -- the ring over which the representation is defined

    - ``cache_matrices`` -- boolean (default: ``True``) if ``True``, then any
      representation matrices that are computed are cached

    EXAMPLES:

    Young's orthogonal representation: the matrices are orthogonal.

    ::

        sage: orth = SymmetricGroupRepresentation([2,1], "orthogonal"); orth
        Orthogonal representation of the symmetric group corresponding to [2, 1]
        sage: all(a*a.transpose() == a.parent().identity_matrix() for a in orth)
        True

    ::

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

def SymmetricGroupRepresentations(n, implementation="specht", ring=None,
        cache_matrices=True):
    r"""
    Irreducible representations of the symmetric group.

    INPUT:

    - ``n`` -- positive integer

    - ``implementation`` -- string (default: ``"specht"``), one of:

      * ``"seminormal"`` - for Young's seminormal representation
      * ``"orthogonal"`` - for Young's orthogonal representation
      * ``"specht"`` - for Specht's representation

    - ``ring`` -- the ring over which the representation is defined

    - ``cache_matrices`` -- boolean (default: ``True``) if ``True``, then any
      representation matrices that are computed are cached

    EXAMPLES:

    Young's orthogonal representation: the matrices are orthogonal.

    ::

        sage: orth = SymmetricGroupRepresentations(3, "orthogonal"); orth
        Orthogonal representations of the symmetric group of order 3! over Symbolic Ring
        sage: orth.list()
        [Orthogonal representation of the symmetric group corresponding to [3],
         Orthogonal representation of the symmetric group corresponding to [2, 1],
         Orthogonal representation of the symmetric group corresponding to [1, 1, 1]]
        sage: orth([2,1])([1,2,3])
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

        The following tests against some bug that was fixed in :trac:`8611`::

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
        Verify the representation: tests that the images of the simple
        transpositions are involutions and tests that the braid relations
        hold.

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

            sage: orth = SymmetricGroupRepresentations(3, "orthogonal")
            sage: for x in orth: print(x)
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

            sage: orth = SymmetricGroupRepresentation([3,2], "orthogonal")
            sage: orth._yang_baxter_graph
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

            sage: orth = SymmetricGroupRepresentation([3,2], "orthogonal")
            sage: orth._tableau_dict
            {(0, -1, 2, 1, 0): [[1, 2, 3], [4, 5]],
             (0, 2, -1, 1, 0): [[1, 2, 4], [3, 5]],
             (0, 2, 1, -1, 0): [[1, 3, 4], [2, 5]],
             (2, 0, -1, 1, 0): [[1, 2, 5], [3, 4]],
             (2, 0, 1, -1, 0): [[1, 3, 5], [2, 4]]}
        """
        # construct a dictionary pairing vertices with tableau
        t = StandardTableaux(self._partition).last()
        tableau_dict = {self._yang_baxter_graph.root(): t}
        for (u, w, (i, beta)) in self._yang_baxter_graph._edges_in_bfs():
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

            sage: orth = SymmetricGroupRepresentation([3,2], "orthogonal")
            sage: orth._word_dict
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

            sage: orth = SymmetricGroupRepresentation([2,1], "orthogonal")
            sage: orth.representation_matrix_for_simple_transposition(1)
            [ 1  0]
            [ 0 -1]
            sage: orth.representation_matrix_for_simple_transposition(2)
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
        if not(1 <= i < sum(self._partition)):
            raise TypeError
        Y = self._yang_baxter_graph
        index_lookup = {b: a for a, b in enumerate(list(Y))}
        digraph = copy(Y._digraph)
        digraph.delete_edges((u, v) for (u, v, (j, beta)) in digraph.edges()
                             if j != i)
        M = matrix(self._ring, digraph.num_verts())
        for g in digraph.connected_components_subgraphs():
            if g.num_verts() == 1:
                [v] = g.vertices()
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
                [(u, v, (j, beta))] = g.edges()
                iu = index_lookup[u]
                iv = index_lookup[v]
                M[iu, iu], M[iu, iv], M[iv, iu], M[iv, iv] = \
                    self._2x2_matrix_entries(beta)
        return M

    def _representation_matrix_uncached(self, permutation):
        r"""
        Return the matrix representing ``permutation``.

        EXAMPLES::

            sage: orth = SymmetricGroupRepresentation([2,1], "orthogonal")
            sage: orth._representation_matrix_uncached(Permutation([2,1,3]))
            [ 1  0]
            [ 0 -1]
            sage: orth._representation_matrix_uncached(Permutation([1,3,2]))
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

            sage: orth = SymmetricGroupRepresentation([2,1], "orthogonal")
            sage: orth.representation_matrix(Permutation([2,1,3]))
            [ 1  0]
            [ 0 -1]
            sage: orth.representation_matrix(Permutation([1,3,2]))
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

            sage: SymmetricGroupRepresentation([2,1], "orthogonal")
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

            sage: orth = SymmetricGroupRepresentation([2,1], "orthogonal")
            sage: orth._2x2_matrix_entries(1/2)
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
            sage: YoungRepresentations_Orthogonal(3)
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

