"""
Characteristic Species
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
from sage.arith.misc import factorial
from .structure import GenericSpeciesStructure
from .set_species import SetSpecies
from sage.structure.unique_representation import UniqueRepresentation


class CharacteristicSpeciesStructure(GenericSpeciesStructure):
    def __repr__(self):
        """
        EXAMPLES::

            sage: F = species.CharacteristicSpecies(3)
            sage: a = F.structures([1, 2, 3]).random_element(); a
            {1, 2, 3}
            sage: F = species.SingletonSpecies()
            sage: F.structures([1]).list()
            [1]
            sage: F = species.EmptySetSpecies()
            sage: F.structures([]).list()
            [{}]
        """
        s = GenericSpeciesStructure.__repr__(self)
        if self.parent()._n == 1:
            return s[1:-1]
        else:
            return "{" + s[1:-1] + "}"

    def canonical_label(self):
        """
        EXAMPLES::

            sage: F = species.CharacteristicSpecies(3)
            sage: a = F.structures(["a", "b", "c"]).random_element(); a
            {'a', 'b', 'c'}
            sage: a.canonical_label()
            {'a', 'b', 'c'}
        """
        P = self.parent()
        rng = list(range(1, P._n + 1))
        return CharacteristicSpeciesStructure(P, self._labels, rng)

    def transport(self, perm):
        """
        Return the transport of this structure along the permutation ``perm``.

        EXAMPLES::

            sage: F = species.CharacteristicSpecies(3)
            sage: a = F.structures(["a", "b", "c"]).random_element(); a
            {'a', 'b', 'c'}
            sage: p = PermutationGroupElement((1,2))                                    # needs sage.groups
            sage: a.transport(p)                                                        # needs sage.groups
            {'a', 'b', 'c'}
        """
        return self

    def automorphism_group(self):
        """
        Return the group of permutations whose action on this structure
        leave it fixed. For the characteristic species, there is only one
        structure, so every permutation is in its automorphism group.

        EXAMPLES::

            sage: F = species.CharacteristicSpecies(3)
            sage: a = F.structures(["a", "b", "c"]).random_element(); a
            {'a', 'b', 'c'}
            sage: a.automorphism_group()                                                # needs sage.groups
            Symmetric group of order 3! as a permutation group
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        return SymmetricGroup(len(self._labels))


class CharacteristicSpecies(GenericCombinatorialSpecies, UniqueRepresentation):
    def __init__(self, n, min=None, max=None, weight=None):
        """
        Return the characteristic species of order `n`.

        This species has exactly one structure on a set of size `n`
        and no structures on sets of any other size.

        EXAMPLES::

            sage: X = species.CharacteristicSpecies(1)
            sage: X.structures([1]).list()
            [1]
            sage: X.structures([1,2]).list()
            []
            sage: X.generating_series()[0:4]
            [0, 1, 0, 0]
            sage: X.isotype_generating_series()[0:4]
            [0, 1, 0, 0]
            sage: X.cycle_index_series()[0:4]                                           # needs sage.modules
            [0, p[1], 0, 0]

            sage: F = species.CharacteristicSpecies(3)
            sage: c = F.generating_series()[0:4]
            sage: F._check()
            True
            sage: F == loads(dumps(F))
            True

        TESTS::

            sage: S1 = species.CharacteristicSpecies(1)
            sage: S2 = species.CharacteristicSpecies(1)
            sage: S3 = species.CharacteristicSpecies(2)
            sage: S4 = species.CharacteristicSpecies(2, weight=2)
            sage: S1 is S2
            True
            sage: S1 == S3
            False
        """
        self._n = n
        self._name = "Characteristic species of order %s" % n
        self._state_info = [n]
        GenericCombinatorialSpecies.__init__(self, min=min, max=max, weight=weight)

    _default_structure_class = CharacteristicSpeciesStructure

    def _structures(self, structure_class, labels):
        """
        EXAMPLES::

            sage: F = species.CharacteristicSpecies(2)
            sage: l = [1, 2, 3]
            sage: F.structures(l).list()
            []
            sage: F = species.CharacteristicSpecies(3)
            sage: F.structures(l).list()
            [{1, 2, 3}]
        """
        if len(labels) == self._n:
            yield structure_class(self, labels, range(1,self._n+1))

    _isotypes = _structures

    def _gs_term(self, base_ring):
        """
        EXAMPLES::

            sage: F = species.CharacteristicSpecies(2)
            sage: F.generating_series()[0:5]
            [0, 0, 1/2, 0, 0]
            sage: F.generating_series().count(2)
            1
        """
        return base_ring(self._weight) / base_ring(factorial(self._n))

    def _order(self):
        """
        Return the order of the generating series.

        EXAMPLES::

            sage: F = species.CharacteristicSpecies(2)
            sage: F._order()
            2
        """
        return self._n

    def _itgs_term(self, base_ring):
        """
        EXAMPLES::

            sage: F = species.CharacteristicSpecies(2)
            sage: F.isotype_generating_series()[0:5]
            [0, 0, 1, 0, 0]

        Here we test out weighting each structure by q.

        ::

            sage: R.<q> = ZZ[]
            sage: Fq = species.CharacteristicSpecies(2, weight=q)
            sage: Fq.isotype_generating_series()[0:5]
            [0, 0, q, 0, 0]
        """
        return base_ring(self._weight)

    def _cis_term(self, base_ring):
        """
        EXAMPLES::

            sage: F = species.CharacteristicSpecies(2)
            sage: g = F.cycle_index_series()                                            # needs sage.modules
            sage: g[0:5]                                                                # needs sage.modules
            [0, 0, 1/2*p[1, 1] + 1/2*p[2], 0, 0]
        """
        cis = SetSpecies(weight=self._weight).cycle_index_series(base_ring)
        return cis.coefficient(self._n)

    def _equation(self, var_mapping):
        """
        Return the right hand side of an algebraic equation satisfied by
        this species. This is a utility function called by the
        algebraic_equation_system method.

        EXAMPLES::

            sage: C = species.CharacteristicSpecies(2)
            sage: Qz = QQ['z']
            sage: R.<node0> = Qz[]
            sage: var_mapping = {'z':Qz.gen(), 'node0':R.gen()}
            sage: C._equation(var_mapping)
            z^2
        """
        return var_mapping['z']**(self._n)


#Backward compatibility
CharacteristicSpecies_class = CharacteristicSpecies


class EmptySetSpecies(CharacteristicSpecies):
    def __init__(self, min=None, max=None, weight=None):
        """
        Return the empty set species.

        This species has exactly one structure on the empty set. It is
        the same (and is implemented) as ``CharacteristicSpecies(0)``.

        EXAMPLES::

            sage: X = species.EmptySetSpecies()
            sage: X.structures([]).list()
            [{}]
            sage: X.structures([1,2]).list()
            []
            sage: X.generating_series()[0:4]
            [1, 0, 0, 0]
            sage: X.isotype_generating_series()[0:4]
            [1, 0, 0, 0]
            sage: X.cycle_index_series()[0:4]                                           # needs sage.modules
            [p[], 0, 0, 0]

        TESTS::

            sage: E1 = species.EmptySetSpecies()
            sage: E2 = species.EmptySetSpecies()
            sage: E1 is E2
            True

            sage: E = species.EmptySetSpecies()
            sage: E._check()
            True
            sage: E == loads(dumps(E))
            True
        """
        CharacteristicSpecies_class.__init__(self, 0, min=min, max=max, weight=weight)
        self._name = "Empty set species"
        self._state_info = []


#Backward compatibility
EmptySetSpecies_class = EmptySetSpecies._cached_constructor = EmptySetSpecies


class SingletonSpecies(CharacteristicSpecies):
    def __init__(self, min=None, max=None, weight=None):
        """
        Return the species of singletons.

        This species has exactly one structure on a set of size `1`. It
        is the same (and is implemented) as ``CharacteristicSpecies(1)``.

        EXAMPLES::

            sage: X = species.SingletonSpecies()
            sage: X.structures([1]).list()
            [1]
            sage: X.structures([1,2]).list()
            []
            sage: X.generating_series()[0:4]
            [0, 1, 0, 0]
            sage: X.isotype_generating_series()[0:4]
            [0, 1, 0, 0]
            sage: X.cycle_index_series()[0:4]                                           # needs sage.modules
            [0, p[1], 0, 0]

        TESTS::

            sage: S1 = species.SingletonSpecies()
            sage: S2 = species.SingletonSpecies()
            sage: S1 is S2
            True

            sage: S = species.SingletonSpecies()
            sage: S._check()
            True
            sage: S == loads(dumps(S))
            True
        """
        CharacteristicSpecies_class.__init__(self, 1, min=min, max=max, weight=weight)
        self._name = "Singleton species"
        self._state_info = []


#Backward compatibility
SingletonSpecies_class = SingletonSpecies
