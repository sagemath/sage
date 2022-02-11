"""
Set Species
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
from .generating_series import _integers_from
from sage.combinat.species.structure import GenericSpeciesStructure
from sage.combinat.species.misc import accept_size
from sage.structure.unique_representation import UniqueRepresentation
from sage.arith.misc import factorial


class SetSpeciesStructure(GenericSpeciesStructure):
    def __repr__(self):
        """
        EXAMPLES::

            sage: S = species.SetSpecies()
            sage: a = S.structures(["a","b","c"]).random_element(); a
            {'a', 'b', 'c'}
        """
        s = GenericSpeciesStructure.__repr__(self)
        return "{"+s[1:-1]+"}"

    def canonical_label(self):
        """
        EXAMPLES::

            sage: S = species.SetSpecies()
            sage: a = S.structures(["a","b","c"]).random_element(); a
            {'a', 'b', 'c'}
            sage: a.canonical_label()
            {'a', 'b', 'c'}
        """
        rng = list(range(1, len(self._labels) + 1))
        return SetSpeciesStructure(self.parent(), self._labels, rng)

    def transport(self, perm):
        """
        Returns the transport of this set along the permutation perm.

        EXAMPLES::

            sage: F = species.SetSpecies()
            sage: a = F.structures(["a", "b", "c"]).random_element(); a
            {'a', 'b', 'c'}
            sage: p = PermutationGroupElement((1,2))
            sage: a.transport(p)
            {'a', 'b', 'c'}
        """
        return self

    def automorphism_group(self):
        """
        Returns the group of permutations whose action on this set leave it
        fixed. For the species of sets, there is only one isomorphism
        class, so every permutation is in its automorphism group.

        EXAMPLES::

            sage: F = species.SetSpecies()
            sage: a = F.structures(["a", "b", "c"]).random_element(); a
            {'a', 'b', 'c'}
            sage: a.automorphism_group()
            Symmetric group of order 3! as a permutation group
        """
        from sage.groups.all import SymmetricGroup
        return SymmetricGroup(max(1,len(self._labels)))

class SetSpecies(GenericCombinatorialSpecies, UniqueRepresentation):
    @staticmethod
    @accept_size
    def __classcall__(cls, *args, **kwds):
        """
        EXAMPLES::

            sage: E = species.SetSpecies(); E
            Set species
        """
        return super(SetSpecies, cls).__classcall__(cls, *args, **kwds)

    def __init__(self, min=None, max=None, weight=None):
        """
        Returns the species of sets.

        EXAMPLES::

            sage: E = species.SetSpecies()
            sage: E.structures([1,2,3]).list()
            [{1, 2, 3}]
            sage: E.isotype_generating_series().coefficients(4)
            [1, 1, 1, 1]

            sage: S = species.SetSpecies()
            sage: c = S.generating_series().coefficients(3)
            sage: S._check()
            True
            sage: S == loads(dumps(S))
            True
        """
        GenericCombinatorialSpecies.__init__(self, min=min, max=max, weight=weight)
        self._name = "Set species"

    _default_structure_class = SetSpeciesStructure

    def _structures(self, structure_class, labels):
        """
        EXAMPLES::

            sage: S = species.SetSpecies()
            sage: S.structures([1,2,3]).list()
            [{1, 2, 3}]
        """
        n = len(labels)
        yield structure_class(self, labels, range(1,n+1))

    _isotypes = _structures

    def _gs_iterator(self, base_ring):
        r"""
        The generating series for the species of sets is given by
        `e^x`.

        EXAMPLES::

            sage: S = species.SetSpecies()
            sage: g = S.generating_series()
            sage: g.coefficients(10)
            [1, 1, 1/2, 1/6, 1/24, 1/120, 1/720, 1/5040, 1/40320, 1/362880]
            sage: [g.count(i) for i in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        for n in _integers_from(0):
            yield base_ring(self._weight / factorial(n))

    def _itgs_list(self, base_ring):
        r"""
        The isomorphism type generating series for the species of sets is
        `\frac{1}{1-x}`.

        EXAMPLES::

            sage: S = species.SetSpecies()
            sage: g = S.isotype_generating_series()
            sage: g.coefficients(10)
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: [g.count(i) for i in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        return [base_ring(self._weight)]

    def _cis(self, series_ring, base_ring):
        r"""
        The cycle index series for the species of sets is given by
        `exp\( \sum_{n=1}{\infty} = \frac{x_n}{n} \)`.

        EXAMPLES::

            sage: S = species.SetSpecies()
            sage: g = S.cycle_index_series()
            sage: g.coefficients(5)
            [p[],
             p[1],
             1/2*p[1, 1] + 1/2*p[2],
             1/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3],
             1/24*p[1, 1, 1, 1] + 1/4*p[2, 1, 1] + 1/8*p[2, 2] + 1/3*p[3, 1] + 1/4*p[4]]
        """
        from .generating_series import ExponentialCycleIndexSeries
        res = ExponentialCycleIndexSeries(base_ring)

        if self.is_weighted():
            res *= self._weight

        return res


#Backward compatibility
SetSpecies_class = SetSpecies
