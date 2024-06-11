"""
Linear-order Species
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
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.species.misc import accept_size


class LinearOrderSpeciesStructure(GenericSpeciesStructure):
    def canonical_label(self):
        """
        EXAMPLES::

            sage: P = species.LinearOrderSpecies()
            sage: s = P.structures(["a", "b", "c"]).random_element()
            sage: s.canonical_label()
            ['a', 'b', 'c']
        """
        return self.__class__(self.parent(), self._labels, range(1, len(self._labels)+1))

    def transport(self, perm):
        """
        Return the transport of this structure along the permutation
        perm.

        EXAMPLES::

            sage: F = species.LinearOrderSpecies()
            sage: a = F.structures(["a", "b", "c"])[0]; a
            ['a', 'b', 'c']
            sage: p = PermutationGroupElement((1,2))                                    # needs sage.groups
            sage: a.transport(p)                                                        # needs sage.groups
            ['b', 'a', 'c']
        """
        return LinearOrderSpeciesStructure(self.parent(), self._labels, [perm(i) for i in self._list])

    def automorphism_group(self):
        """
        Return the group of permutations whose action on this structure
        leave it fixed. For the species of linear orders, there is no
        non-trivial automorphism.

        EXAMPLES::

            sage: F = species.LinearOrderSpecies()
            sage: a = F.structures(["a", "b", "c"])[0]; a
            ['a', 'b', 'c']
            sage: a.automorphism_group()                                                # needs sage.groups
            Symmetric group of order 1! as a permutation group
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        return SymmetricGroup(1)


class LinearOrderSpecies(GenericCombinatorialSpecies, UniqueRepresentation):
    @staticmethod
    @accept_size
    def __classcall__(cls, *args, **kwds):
        r"""
        EXAMPLES::

            sage: L = species.LinearOrderSpecies(); L
            Linear order species
        """
        return super().__classcall__(cls, *args, **kwds)

    def __init__(self, min=None, max=None, weight=None):
        """
        Return the species of linear orders.

        EXAMPLES::

            sage: L = species.LinearOrderSpecies()
            sage: L.generating_series()[0:5]
            [1, 1, 1, 1, 1]

            sage: L = species.LinearOrderSpecies()
            sage: L._check()
            True
            sage: L == loads(dumps(L))
            True
        """
        GenericCombinatorialSpecies.__init__(self, min=min, max=max, weight=None)
        self._name = "Linear order species"

    _default_structure_class = LinearOrderSpeciesStructure

    def _structures(self, structure_class, labels):
        """
        EXAMPLES::

            sage: L = species.LinearOrderSpecies()
            sage: L.structures([1,2,3]).list()
            [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
        """
        from sage.combinat.permutation import Permutations
        for p in Permutations(len(labels)):
            yield structure_class(self, labels, p._list)

    def _isotypes(self, structure_class, labels):
        """
        EXAMPLES::

            sage: L = species.LinearOrderSpecies()
            sage: L.isotypes([1,2,3]).list()
            [[1, 2, 3]]
        """
        yield structure_class(self, labels, range(1, len(labels)+1))

    def _gs_list(self, base_ring, n):
        r"""
        The generating series for the species of linear orders is
        `\frac{1}{1-x}`.

        EXAMPLES::

            sage: L = species.LinearOrderSpecies()
            sage: g = L.generating_series()
            sage: g[0:10]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        return base_ring.one()

    def _itgs_list(self, base_ring, n):
        r"""
        The isomorphism type generating series is given by
        `\frac{1}{1-x}`.

        EXAMPLES::

            sage: L = species.LinearOrderSpecies()
            sage: g = L.isotype_generating_series()
            sage: g[0:10]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        return base_ring.one()

    def _cis_callable(self, base_ring, n):
        """
        EXAMPLES::

            sage: L = species.LinearOrderSpecies()
            sage: g = L.cycle_index_series()                                            # needs sage.modules
            sage: g[0:5]                                                                # needs sage.modules
            [p[], p[1], p[1, 1], p[1, 1, 1], p[1, 1, 1, 1]]
        """
        from sage.combinat.sf.sf import SymmetricFunctions
        p = SymmetricFunctions(base_ring).power()
        return p([1]*n)


#Backward compatibility
LinearOrderSpecies_class = LinearOrderSpecies
