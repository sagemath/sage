"""
Cycle Species
"""

#*****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.arith.misc import divisors, euler_phi
from sage.combinat.species.misc import accept_size
from sage.combinat.species.species import GenericCombinatorialSpecies
from sage.combinat.species.structure import GenericSpeciesStructure
from sage.structure.unique_representation import UniqueRepresentation


class CycleSpeciesStructure(GenericSpeciesStructure):
    def __repr__(self):
        """
        EXAMPLES::

            sage: S = species.CycleSpecies()
            sage: S.structures(["a","b","c"])[0]
            ('a', 'b', 'c')
        """
        s = GenericSpeciesStructure.__repr__(self)
        return "("+s[1:-1]+")"

    def canonical_label(self):
        """
        EXAMPLES::

            sage: P = species.CycleSpecies()
            sage: P.structures(["a","b","c"]).random_element().canonical_label()
            ('a', 'b', 'c')
        """
        n = len(self._labels)
        return CycleSpeciesStructure(self.parent(), self._labels, range(1, n+1))

    def permutation_group_element(self):
        """
        Return this cycle as a permutation group element.

        EXAMPLES::

            sage: F = species.CycleSpecies()
            sage: a = F.structures(["a", "b", "c"])[0]; a
            ('a', 'b', 'c')
            sage: a.permutation_group_element()                                         # needs sage.groups
            (1,2,3)
        """
        from sage.groups.perm_gps.constructor import PermutationGroupElement
        return PermutationGroupElement(tuple(self._list))

    def transport(self, perm):
        """
        Return the transport of this structure along the permutation
        perm.

        EXAMPLES::

            sage: F = species.CycleSpecies()
            sage: a = F.structures(["a", "b", "c"])[0]; a
            ('a', 'b', 'c')
            sage: p = PermutationGroupElement((1,2))                                    # needs sage.groups
            sage: a.transport(p)                                                        # needs sage.groups
            ('a', 'c', 'b')
        """
        p = self.permutation_group_element()
        p = perm*p*~perm
        new_list = [1]
        for i in range(len(self._list)-1):
            new_list.append( p(new_list[-1]) )
        return CycleSpeciesStructure(self.parent(), self._labels, new_list)

    def automorphism_group(self):
        """
        Return the group of permutations whose action on this structure
        leave it fixed.

        EXAMPLES::

            sage: P = species.CycleSpecies()
            sage: a = P.structures([1, 2, 3, 4])[0]; a
            (1, 2, 3, 4)
            sage: a.automorphism_group()                                                # needs sage.groups
            Permutation Group with generators [(1,2,3,4)]

        ::

            sage: [a.transport(perm) for perm in a.automorphism_group()]                # needs sage.groups
            [(1, 2, 3, 4), (1, 2, 3, 4), (1, 2, 3, 4), (1, 2, 3, 4)]
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        from sage.groups.perm_gps.permgroup import PermutationGroup
        S = SymmetricGroup(len(self._labels))
        p = self.permutation_group_element()
        return PermutationGroup(S.centralizer(p).gens())


class CycleSpecies(GenericCombinatorialSpecies, UniqueRepresentation):
    @staticmethod
    @accept_size
    def __classcall__(cls, *args, **kwds):
        r"""
        EXAMPLES::

            sage: C = species.CycleSpecies(); C
            Cyclic permutation species
        """
        return super().__classcall__(cls, *args, **kwds)

    def __init__(self, min=None, max=None, weight=None):
        """
        Return the species of cycles.

        EXAMPLES::

            sage: C = species.CycleSpecies(); C
            Cyclic permutation species
            sage: C.structures([1,2,3,4]).list()
            [(1, 2, 3, 4),
             (1, 2, 4, 3),
             (1, 3, 2, 4),
             (1, 3, 4, 2),
             (1, 4, 2, 3),
             (1, 4, 3, 2)]

        TESTS:

        We check to verify that the caching of species is actually
        working.

        ::

            sage: species.CycleSpecies() is species.CycleSpecies()
            True

            sage: P = species.CycleSpecies()
            sage: c = P.generating_series()[:3]
            sage: P._check()
            True
            sage: P == loads(dumps(P))
            True
        """
        GenericCombinatorialSpecies.__init__(self, min=min, max=max, weight=weight)
        self._name = "Cyclic permutation species"

    _default_structure_class = CycleSpeciesStructure

    def _structures(self, structure_class, labels):
        """
        EXAMPLES::

            sage: P = species.CycleSpecies()
            sage: P.structures([1,2,3]).list()
            [(1, 2, 3), (1, 3, 2)]
        """
        from sage.combinat.permutation import CyclicPermutations
        for c in CyclicPermutations(range(1, len(labels)+1)):
            yield structure_class(self, labels, c)

    def _isotypes(self, structure_class, labels):
        """
        EXAMPLES::

            sage: P = species.CycleSpecies()
            sage: P.isotypes([1,2,3]).list()
            [(1, 2, 3)]
        """
        if len(labels) != 0:
            yield structure_class(self, labels, range(1, len(labels)+1))

    def _gs_callable(self, base_ring, n):
        r"""
        The generating series for cyclic permutations is
        `-\log(1-x) = \sum_{n=1}^\infty x^n/n`.

        EXAMPLES::

            sage: P = species.CycleSpecies()
            sage: g = P.generating_series()
            sage: g[0:10]
            [0, 1, 1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 1/9]

        TESTS::

            sage: P = species.CycleSpecies()
            sage: g = P.generating_series(RR)
            sage: g[0:3]
            [0.000000000000000, 1.00000000000000, 0.500000000000000]
        """
        if n:
            return self._weight * base_ring.one() / n
        return base_ring.zero()

    def _order(self):
        """
        Return the order of the generating series.

        EXAMPLES::

            sage: P = species.CycleSpecies()
            sage: P._order()
            1
        """
        return 1

    def _itgs_list(self, base_ring, n):
        """
        The isomorphism type generating series for cyclic permutations is
        given by `x/(1-x)`.

        EXAMPLES::

            sage: P = species.CycleSpecies()
            sage: g = P.isotype_generating_series()
            sage: g[0:5]
            [0, 1, 1, 1, 1]

        TESTS::

            sage: P = species.CycleSpecies()
            sage: g = P.isotype_generating_series(RR)
            sage: g[0:3]
            [0.000000000000000, 1.00000000000000, 1.00000000000000]
        """
        if n:
            return self._weight * base_ring.one()
        return base_ring.zero()

    def _cis_callable(self, base_ring, n):
        r"""
        The cycle index series of the species of cyclic permutations is
        given by

        .. MATH::

             -\sum_{k=1}^\infty \phi(k)/k * log(1 - x_k)


        which is equal to

        .. MATH::

             \sum_{n=1}^\infty \frac{1}{n} * \sum_{k|n} \phi(k) * x_k^{n/k}

        .

        EXAMPLES::

            sage: P = species.CycleSpecies()
            sage: cis = P.cycle_index_series()                                          # needs sage.modules
            sage: cis[0:7]                                                              # needs sage.modules
            [0,
             p[1],
             1/2*p[1, 1] + 1/2*p[2],
             1/3*p[1, 1, 1] + 2/3*p[3],
             1/4*p[1, 1, 1, 1] + 1/4*p[2, 2] + 1/2*p[4],
             1/5*p[1, 1, 1, 1, 1] + 4/5*p[5],
             1/6*p[1, 1, 1, 1, 1, 1] + 1/6*p[2, 2, 2] + 1/3*p[3, 3] + 1/3*p[6]]
        """
        from sage.combinat.sf.sf import SymmetricFunctions
        p = SymmetricFunctions(base_ring).power()

        zero = base_ring.zero()

        if not n:
            return zero
        res = zero
        for k in divisors(n):
            res += euler_phi(k)*p([k])**(n//k)
        res /= n
        return self._weight * res


#Backward compatibility
CycleSpecies_class = CycleSpecies
