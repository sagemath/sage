"""
Composition species
"""
#*****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
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
#*****************************************************************************
from .species import GenericCombinatorialSpecies
from .structure import GenericSpeciesStructure
from .partition_species import PartitionSpecies
from sage.structure.unique_representation import UniqueRepresentation


class CompositionSpeciesStructure(GenericSpeciesStructure):
    def __init__(self, parent, labels, pi, f, gs):
        """
        TESTS::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: L = E(C)
            sage: a = L.structures(['a','b','c']).random_element()                      # needs sage.libs.flint
            sage: a == loads(dumps(a))                                                  # needs sage.libs.flint
            True
        """
        self._partition = pi
        GenericSpeciesStructure.__init__(self, parent, labels, [f, gs])

    def __repr__(self):
        """
        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: L = E(C)
            sage: L.structures(['a','b','c'])[0]                                        # needs sage.libs.flint
            F-structure: {{'a', 'b', 'c'}}; G-structures: (('a', 'b', 'c'),)
        """
        f, gs = self._list
        return "F-structure: %s; G-structures: %s" % (repr(f), repr(gs))

    def transport(self, perm):
        """
        EXAMPLES::

            sage: p = PermutationGroupElement((2,3))                                    # needs sage.groups
            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: L = E(C)
            sage: S = L.structures(['a','b','c']).list()                                # needs sage.libs.flint
            sage: a = S[2]; a                                                           # needs sage.libs.flint
            F-structure: {{'a', 'c'}, {'b'}}; G-structures: (('a', 'c'), ('b'))
            sage: a.transport(p)                                                        # needs sage.groups sage.libs.flint
            F-structure: {{'a', 'b'}, {'c'}}; G-structures: (('a', 'c'), ('b'))
        """
        f, gs = self._list
        pi = self._partition.transport(perm)
        f = f.change_labels(pi._list)
        _ = [g.change_labels(part) for g, part in zip(gs, pi)]  # TODO: BUG HERE ?
        return self.__class__(self, self._labels, pi, f, gs)

    def change_labels(self, labels):
        """
        Return a relabelled structure.

        INPUT:

        - ``labels`` -- list of labels

        OUTPUT:

        A structure with the `i`-th label of ``self`` replaced with the `i`-th
        label of the list.

        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: L = E(C)
            sage: S = L.structures(['a','b','c']).list()                                # needs sage.libs.flint
            sage: a = S[2]; a                                                           # needs sage.libs.flint
            F-structure: {{'a', 'c'}, {'b'}}; G-structures: (('a', 'c'), ('b'))
            sage: a.change_labels([1,2,3])                                              # needs sage.libs.flint
            F-structure: {{1, 3}, {2}}; G-structures: [(1, 3), (2)]
        """
        f, gs = self._list
        pi = self._partition.change_labels(labels)
        f = f.change_labels(list(pi))
        g = [g.change_labels(part) for g, part in zip(gs, pi)]
        return self.__class__(self, labels, pi, f, g)


class CompositionSpecies(GenericCombinatorialSpecies, UniqueRepresentation):
    def __init__(self, F, G, min=None, max=None, weight=None):
        """
        Return the composition of two species.

        EXAMPLES::

            sage: E = species.SetSpecies()
            sage: C = species.CycleSpecies()
            sage: S = E(C)
            sage: S.generating_series()[:5]
            [1, 1, 1, 1, 1]
            sage: E(C) is S
            True

        TESTS::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: L = E(C)
            sage: c = L.generating_series()[:3]
            sage: L._check()  #False due to isomorphism types not being implemented     # needs sage.libs.flint
            False
            sage: L == loads(dumps(L))
            True
        """
        self._F = F
        self._G = G
        self._name = "Composition of (%s) and (%s)" % (F, G)
        self._state_info = [F, G]
        GenericCombinatorialSpecies.__init__(self, min=None, max=None, weight=None)

    _default_structure_class = CompositionSpeciesStructure

    def _structures(self, structure_class, labels):
        """
        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: L = E(C)
            sage: L.structures(['a','b','c']).list()                                    # needs sage.libs.flint
            [F-structure: {{'a', 'b', 'c'}}; G-structures: (('a', 'b', 'c'),),
             F-structure: {{'a', 'b', 'c'}}; G-structures: (('a', 'c', 'b'),),
             F-structure: {{'a', 'c'}, {'b'}}; G-structures: (('a', 'c'), ('b')),
             F-structure: {{'a', 'b'}, {'c'}}; G-structures: (('a', 'b'), ('c')),
             F-structure: {{'b', 'c'}, {'a'}}; G-structures: (('b', 'c'), ('a')),
             F-structure: {{'a'}, {'b'}, {'c'}}; G-structures: (('a'), ('b'), ('c'))]

        TESTS::

            sage: # needs sage.libs.flint
            sage: a = _[2]
            sage: f, gs = a._list
            sage: f
            {{'a', 'c'}, {'b'}}
            sage: f.parent()
            Set species
            sage: f._list
            [1, 2]
            sage: f._labels
            [{'a', 'c'}, {'b'}]
            sage: [g.parent() for g in gs]
            [Cyclic permutation species, Cyclic permutation species]
            sage: [g._labels for g in gs]
            [['a', 'c'], ['b']]
            sage: [g._list for g in gs]
            [[1, 2], [1]]
        """
        from itertools import product
        P = PartitionSpecies()
        for pi in P.structures(labels):
            # The labels of the G-structures will be just be the things
            # in labels
            gs = product(*[self._G.structures(part.label_subset()) for part in pi])

            # The labels of the F-structure will be set objects
            fs = self._F.structures(list(pi))
            for f, gg in product(fs, gs):
                yield structure_class(self, labels, pi, f, gg)

    def _isotypes(self, structure_class, labels):
        """
        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: L = E(C)
            sage: L.isotypes(['a','b','c']).list()                                      # needs sage.modules
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _gs(self, series_ring, base_ring):
        """
        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: L = E(C)
            sage: L.generating_series()[:5]
            [1, 1, 1, 1, 1]
        """
        return self._F.generating_series(base_ring)(self._G.generating_series(base_ring))

    def _itgs(self, series_ring, base_ring):
        """
        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: L = E(C)
            sage: L.isotype_generating_series()[:10]                                    # needs sage.modules
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]
        """
        cis = self.cycle_index_series(base_ring)
        return cis.isotype_generating_series()

    def _cis(self, series_ring, base_ring):
        """
        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: L = E(C)
            sage: L.cycle_index_series()[:5]                                            # needs sage.modules
            [p[],
             p[1],
             p[1, 1] + p[2],
             p[1, 1, 1] + p[2, 1] + p[3],
             p[1, 1, 1, 1] + p[2, 1, 1] + p[2, 2] + p[3, 1] + p[4]]

        Here we (indirectly) check that the cycle index series for
        permutations weighted by the number of cycles is correctly
        computed.

        ::

            sage: t = QQ['t'].gen()
            sage: E = species.SetSpecies()
            sage: C = species.CycleSpecies(weight=t)
            sage: S = E(C)
            sage: S.isotype_generating_series()[:5]  #indirect                          # needs sage.modules
            [1, t, t^2 + t, t^3 + t^2 + t, t^4 + t^3 + 2*t^2 + t]

        We do the same thing with set partitions weighted by the number of
        blocks.

        ::

            sage: t = QQ['t'].gen()
            sage: E = species.SetSpecies()
            sage: E_t = species.SetSpecies(min=1,weight=t)
            sage: Par = E(E_t)
            sage: Par.isotype_generating_series()[:5]                                   # needs sage.modules
            [1, t, t^2 + t, t^3 + t^2 + t, t^4 + t^3 + 2*t^2 + t]
        """
        f_cis = self._F.cycle_index_series(base_ring)
        g_cis = self._G.cycle_index_series(base_ring)
        return f_cis(g_cis)

    def weight_ring(self):
        """
        Return the weight ring for this species. This is determined by
        asking Sage's coercion model what the result is when you multiply
        (and add) elements of the weight rings for each of the operands.

        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: L = E(C)
            sage: L.weight_ring()
            Rational Field
        """
        return self._common_parent([self._F.weight_ring(), self._G.weight_ring()])


CompositionSpecies_class = CompositionSpecies
