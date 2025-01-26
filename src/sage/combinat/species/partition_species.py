# sage.doctest: needs sage.libs.flint
"""
Partition species
"""
# ****************************************************************************
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
# ****************************************************************************

from .species import GenericCombinatorialSpecies
from sage.arith.misc import factorial
from .subset_species import SubsetSpeciesStructure
from .set_species import SetSpecies
from .structure import GenericSpeciesStructure
from sage.combinat.species.misc import accept_size
from functools import reduce


class PartitionSpeciesStructure(GenericSpeciesStructure):
    def __init__(self, parent, labels, list):
        """
        EXAMPLES::

            sage: from sage.combinat.species.partition_species import PartitionSpeciesStructure
            sage: P = species.PartitionSpecies()
            sage: s = PartitionSpeciesStructure(P, ['a','b','c'], [[1,2],[3]]); s
            {{'a', 'b'}, {'c'}}
            sage: s == loads(dumps(s))
            True
        """
        list = [SubsetSpeciesStructure(parent, labels, block) if not isinstance(block, SubsetSpeciesStructure) else block for block in list]
        list.sort(key=lambda block:(-len(block), block))
        GenericSpeciesStructure.__init__(self, parent, labels, list)

    def __repr__(self):
        """
        EXAMPLES::

            sage: S = species.PartitionSpecies()
            sage: a = S.structures(["a","b","c"])[0]; a
            {{'a', 'b', 'c'}}
        """
        s = GenericSpeciesStructure.__repr__(self)
        return "{"+s[1:-1]+"}"

    def canonical_label(self):
        """
        EXAMPLES::

            sage: P = species.PartitionSpecies()
            sage: S = P.structures(["a", "b", "c"])
            sage: [s.canonical_label() for s in S]
            [{{'a', 'b', 'c'}},
             {{'a', 'b'}, {'c'}},
             {{'a', 'b'}, {'c'}},
             {{'a', 'b'}, {'c'}},
             {{'a'}, {'b'}, {'c'}}]
        """
        P = self.parent()
        p = [len(block) for block in self._list]
        return P._canonical_rep_from_partition(self.__class__, self._labels, p)

    def transport(self, perm):
        """
        Return the transport of this set partition along the permutation
        perm. For set partitions, this is the direct product of the
        automorphism groups for each of the blocks.

        EXAMPLES::

            sage: p = PermutationGroupElement((2,3))
            sage: from sage.combinat.species.partition_species import PartitionSpeciesStructure
            sage: a = PartitionSpeciesStructure(None, [2,3,4], [[1,2],[3]]); a
            {{2, 3}, {4}}
            sage: a.transport(p)
            {{2, 4}, {3}}
        """
        l = [block.transport(perm)._list for block in self._list]
        l.sort(key=lambda block:(-len(block), block))
        return PartitionSpeciesStructure(self.parent(), self._labels, l)

    def automorphism_group(self):
        """
        Return the group of permutations whose action on this set
        partition leave it fixed.

        EXAMPLES::

            sage: p = PermutationGroupElement((2,3))
            sage: from sage.combinat.species.partition_species import PartitionSpeciesStructure
            sage: a = PartitionSpeciesStructure(None, [2,3,4], [[1,2],[3]]); a
            {{2, 3}, {4}}
            sage: a.automorphism_group()
            Permutation Group with generators [(1,2)]
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        return reduce(lambda a,b: a.direct_product(b, maps=False),
                      [SymmetricGroup(block._list) for block in self._list])

    def change_labels(self, labels):
        """
        Return a relabelled structure.

        INPUT:

        - ``labels`` -- list of labels

        OUTPUT:

        A structure with the `i`-th label of ``self`` replaced with the `i`-th
        label of the list.

        EXAMPLES::

            sage: p = PermutationGroupElement((2,3))
            sage: from sage.combinat.species.partition_species import PartitionSpeciesStructure
            sage: a = PartitionSpeciesStructure(None, [2,3,4], [[1,2],[3]]); a
            {{2, 3}, {4}}
            sage: a.change_labels([1,2,3])
            {{1, 2}, {3}}
        """
        return PartitionSpeciesStructure(self.parent(), labels, [block.change_labels(labels) for block in self._list])


class PartitionSpecies(GenericCombinatorialSpecies):
    @staticmethod
    @accept_size
    def __classcall__(cls, *args, **kwds):
        """
        EXAMPLES::

            sage: P = species.PartitionSpecies(); P
            Partition species
        """
        return super().__classcall__(cls, *args, **kwds)

    def __init__(self, min=None, max=None, weight=None):
        """
        Return the species of partitions.

        EXAMPLES::

            sage: P = species.PartitionSpecies()
            sage: P.generating_series()[0:5]
            [1, 1, 1, 5/6, 5/8]
            sage: P.isotype_generating_series()[0:5]
            [1, 1, 2, 3, 5]

            sage: P = species.PartitionSpecies()
            sage: P._check()
            True
            sage: P == loads(dumps(P))
            True
        """
        GenericCombinatorialSpecies.__init__(self, min=min, max=max, weight=weight)
        self._name = "Partition species"

    _default_structure_class = PartitionSpeciesStructure

    def _structures(self, structure_class, labels):
        """
        EXAMPLES::

            sage: P = species.PartitionSpecies()
            sage: P.structures([1,2,3]).list()
            [{{1, 2, 3}}, {{1, 3}, {2}}, {{1, 2}, {3}}, {{2, 3}, {1}}, {{1}, {2}, {3}}]
        """
        from sage.combinat.restricted_growth import RestrictedGrowthArrays
        n = len(labels)

        if n == 0:
            yield structure_class(self, labels, [])
            return

        u = list(range(n, 0, -1))
        s0 = u.pop()

        # Reconstruct the set partitions from
        # restricted growth arrays
        for a in RestrictedGrowthArrays(n):
            m = a.pop(0)
            r = [[] for _ in range(m)]
            i = n
            for i, z in enumerate(u):
                r[a[i]].append(z)
            r[0].append(s0)

            for sp in r:
                sp.reverse()

            r.sort(key=len, reverse=True)

            yield structure_class(self, labels, r)

    def _isotypes(self, structure_class, labels):
        """
        EXAMPLES::

            sage: P = species.PartitionSpecies()
            sage: P.isotypes([1,2,3,4]).list()
            [{{1, 2, 3, 4}},
             {{1, 2, 3}, {4}},
             {{1, 2}, {3, 4}},
             {{1, 2}, {3}, {4}},
             {{1}, {2}, {3}, {4}}]
        """
        from sage.combinat.partition import Partitions
        for p in Partitions(len(labels)):
            yield self._canonical_rep_from_partition(structure_class, labels, p)

    def _canonical_rep_from_partition(self, structure_class, labels, p):
        """
        Return the canonical representative corresponding to the partition
        p.

        EXAMPLES::

            sage: P = species.PartitionSpecies()
            sage: P._canonical_rep_from_partition(P._default_structure_class,[1,2,3],[2,1])
            {{1, 2}, {3}}
        """
        breaks = [sum(p[:i]) for i in range(len(p) + 1)]
        return structure_class(self, labels, [list(range(breaks[i]+1, breaks[i+1]+1)) for i in range(len(p))])

    def _gs_callable(self, base_ring, n):
        r"""
        EXAMPLES::

            sage: P = species.PartitionSpecies()
            sage: g = P.generating_series()
            sage: [g.coefficient(i) for i in range(5)]
            [1, 1, 1, 5/6, 5/8]
        """
        from sage.combinat.combinat import bell_number
        return self._weight * base_ring(bell_number(n) / factorial(n))

    def _itgs_callable(self, base_ring, n):
        r"""
        The isomorphism type generating series is given by
        `\frac{1}{1-x}`.

        EXAMPLES::

            sage: P = species.PartitionSpecies()
            sage: g = P.isotype_generating_series()
            sage: [g.coefficient(i) for i in range(10)]
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]
        """
        from sage.combinat.partition import number_of_partitions
        return self._weight*base_ring(number_of_partitions(n))

    def _cis(self, series_ring, base_ring):
        r"""
        The cycle index series for the species of partitions is given by.

        .. MATH::

             exp \sum_{n \ge 1} \frac{1}{n} \left( exp \left( \sum_{k \ge 1} \frac{x_{kn}}{k} \right) -1 \right).

        EXAMPLES::

            sage: P = species.PartitionSpecies()
            sage: g = P.cycle_index_series()                                            # needs sage.modules
            sage: g[0:5]                                                                # needs sage.modules
            [p[],
             p[1],
             p[1, 1] + p[2],
             5/6*p[1, 1, 1] + 3/2*p[2, 1] + 2/3*p[3],
             5/8*p[1, 1, 1, 1] + 7/4*p[2, 1, 1] + 7/8*p[2, 2] + p[3, 1] + 3/4*p[4]]
        """
        ciset = SetSpecies().cycle_index_series(base_ring)
        res = ciset(ciset - 1)
        if self.is_weighted():
            res *= self._weight
        return res


#Backward compatibility
PartitionSpecies_class = PartitionSpecies
