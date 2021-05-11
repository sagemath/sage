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
from .generating_series import _integers_from
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
        Returns the transport of this structure along the permutation
        perm.

        EXAMPLES::

            sage: F = species.LinearOrderSpecies()
            sage: a = F.structures(["a", "b", "c"])[0]; a
            ['a', 'b', 'c']
            sage: p = PermutationGroupElement((1,2))
            sage: a.transport(p)
            ['b', 'a', 'c']
        """
        return LinearOrderSpeciesStructure(self.parent(), self._labels, [perm(i) for i in self._list])

    def automorphism_group(self):
        """
        Returns the group of permutations whose action on this structure
        leave it fixed. For the species of linear orders, there is no
        non-trivial automorphism.

        EXAMPLES::

            sage: F = species.LinearOrderSpecies()
            sage: a = F.structures(["a", "b", "c"])[0]; a
            ['a', 'b', 'c']
            sage: a.automorphism_group()
            Symmetric group of order 1! as a permutation group
        """
        from sage.groups.all import SymmetricGroup
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
        return super(LinearOrderSpecies, cls).__classcall__(cls, *args, **kwds)

    def __init__(self, min=None, max=None, weight=None):
        """
        Returns the species of linear orders.

        EXAMPLES::

            sage: L = species.LinearOrderSpecies()
            sage: L.generating_series().coefficients(5)
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

    def _gs_list(self, base_ring):
        r"""
        The generating series for the species of linear orders is
        `\frac{1}{1-x}`.

        EXAMPLES::

            sage: L = species.LinearOrderSpecies()
            sage: g = L.generating_series()
            sage: g.coefficients(10)
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        return [base_ring(1)]

    def _itgs_list(self, base_ring):
        r"""
        The isomorphism type generating series is given by
        `\frac{1}{1-x}`.

        EXAMPLES::

            sage: L = species.LinearOrderSpecies()
            sage: g = L.isotype_generating_series()
            sage: g.coefficients(10)
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        return [base_ring(1)]


    def _cis_iterator(self, base_ring):
        """
        EXAMPLES::

            sage: L = species.LinearOrderSpecies()
            sage: g = L.cycle_index_series()
            sage: g.coefficients(5)
            [p[], p[1], p[1, 1], p[1, 1, 1], p[1, 1, 1, 1]]
        """
        from sage.combinat.sf.sf import SymmetricFunctions
        p = SymmetricFunctions(base_ring).power()
        for n in _integers_from(0):
            yield p([1]*n)

#Backward compatibility
LinearOrderSpecies_class = LinearOrderSpecies
