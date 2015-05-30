"""
Species structures

We will illustrate the use of the structure classes using the
"balls and bars" model for integer compositions. An integer
composition of 6 such as [2, 1, 3] can be represented in this model
as 'oooooo' where the 6 o's correspond to the balls and the 2 's
correspond to the bars. If BB is our species for this model, the it
satisfies the following recursive definition:

BB = o + o\*BB + o\*|\*BB

Here we define this species using the default structures::

    sage: ball = species.SingletonSpecies(); o = var('o')
    sage: bar = species.EmptySetSpecies()
    sage: BB = CombinatorialSpecies()
    sage: BB.define(ball + ball*BB + ball*bar*BB)
    sage: BB.isotypes([o]*3).list()
    [o*(o*o), o*((o*{})*o), (o*{})*(o*o), (o*{})*((o*{})*o)]

If we ignore the parentheses, we can read off that the integer
compositions are [3], [2, 1], [1, 2], and [1, 1, 1].
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
from sage.combinat.combinat import CombinatorialClass, CombinatorialObject
from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from copy import copy

class SpeciesStructure(SageObject):
    def parent(self):
        """
        Returns the species that this structure is associated with.

        EXAMPLES::

            sage: L = species.LinearOrderSpecies()
            sage: a,b = L.structures([1,2])
            sage: a.parent()
            Linear order species
        """
        try:
            return self._parent
        except AttributeError:
            raise NotImplementedError

class GenericSpeciesStructure(CombinatorialObject, SpeciesStructure):
    def __init__(self, parent, labels, list):
        """
        EXAMPLES::

            sage: from sage.combinat.species.structure import GenericSpeciesStructure
            sage: a = GenericSpeciesStructure(None, [2,3,4], [1,2,3])
            sage: a
            [2, 3, 4]
            sage: a.parent() is None
            True
            sage: a == loads(dumps(a))
            True
        """
        self._parent = parent
        self._labels = labels
        CombinatorialObject.__init__(self, list)

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.species.structure import GenericSpeciesStructure
            sage: a = GenericSpeciesStructure(None, [2,3,4], [1,2,3])
            sage: a
            [2, 3, 4]
        """
        return repr(self.labels())

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: T = species.BinaryTreeSpecies()
            sage: t = T.structures([1,2,3])[0]; t
            1*(2*3)
            sage: t[0], t[1][0]
            (1, 2)
            sage: t[0] == t[1][0]
            False
        """
        if not isinstance(self, type(other)):
            return False
        return self._list == other._list and self.labels() == other.labels()

    def labels(self):
        """
        EXAMPLES::

            sage: P = species.SubsetSpecies()
            sage: S = P.structures(["a", "b", "c"])
            sage: [s.labels() for s in S]
            [[], ['a'], ['b'], ['c'], ['a', 'b'], ['a', 'c'], ['b', 'c'], ['a', 'b', 'c']]
        """
        return [self._relabel(i) for i in self._list]

    def change_labels(self, labels):
        """
        EXAMPLES::

            sage: P = species.SubsetSpecies()
            sage: S = P.structures(["a", "b", "c"])
            sage: [s.change_labels([1,2,3]) for s in S]
            [{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}]
        """
        c = copy(self)
        c._labels = labels
        return c

    def _relabel(self, i):
        """
        EXAMPLES::

            sage: from sage.combinat.species.structure import GenericSpeciesStructure
            sage: a = GenericSpeciesStructure(None, [2,3,4], [1,2,3])
            sage: a._relabel(1)
            2
            sage: a._relabel([1,2,3])
            [1, 2, 3]
        """
        if isinstance(i, (int, Integer)):
            return self._labels[i-1]
        else:
            return i

    def is_isomorphic(self, x):
        """
        EXAMPLES::

            sage: S = species.SetSpecies()
            sage: a = S.structures([1,2,3]).random_element(); a
            {1, 2, 3}
            sage: b = S.structures(['a','b','c']).random_element(); b
            {'a', 'b', 'c'}
            sage: a.is_isomorphic(b)
            True
        """
        if self.__class__ != x.__class__:
            return False
        if self.parent() != x.parent():
            return False

        #We don't care about the labels for isomorphism testing
        if self.canonical_label()._list == x.canonical_label()._list:
            return True
        else:
            return False

class SpeciesStructureWrapper(GenericSpeciesStructure):
    def __init__(self, parent, s, **options):
        """
        EXAMPLES::

            sage: E = species.SetSpecies(); B = E+E
            sage: s = B.structures([1,2,3]).random_element()
            sage: s.parent()
            Sum of (Set species) and (Set species)
            sage: s == loads(dumps(s))
            True
        """
        self._parent = parent
        self._s = s
        self._options = options
        GenericSpeciesStructure.__init__(self, parent, s._labels, s._list)

    def __getattr__(self, attr):
        """
        EXAMPLES::

            sage: E = species.SetSpecies(); B = E+E
            sage: s = B.structures([1,2,3]).random_element()
            sage: s
            {1, 2, 3}
        """
        if attr == "_s":
            return None

        return getattr(self._s, attr)

    def __repr__(self):
        """
        Returns the repr of the object which this one wraps.

        EXAMPLES::

            sage: E = species.SetSpecies()
            sage: s = (E+E).structures([1,2,3]).random_element(); s
            {1, 2, 3}
        """
        return repr(self._s)

    def transport(self, perm):
        """
        EXAMPLES::

            sage: P = species.PartitionSpecies()
            sage: s = (P+P).structures([1,2,3]).random_element(); s
            {{1, 3}, {2}}
            sage: s.transport(PermutationGroupElement((2,3)))
            {{1, 2}, {3}}
        """
        return self.__class__(self._parent, self._s.transport(perm), **self._options)

    def canonical_label(self):
        """
        EXAMPLES::

            sage: P = species.PartitionSpecies()
            sage: s = (P+P).structures([1,2,3]).random_element(); s
            {{1, 3}, {2}}
            sage: s.canonical_label()
            {{1, 2}, {3}}
        """
        return self.__class__(self._parent, self._s.canonical_label(), **self._options)

    def labels(self):
        """
        EXAMPLES::

            sage: P = species.PartitionSpecies()
            sage: s = (P+P).structures([1,2,3]).random_element(); s
            {{1, 3}, {2}}
            sage: s.labels()
            [{1, 3}, {2}]
            sage: type(_)
            <type 'list'>
        """
        return self._s.labels()


##############################################################


class SpeciesWrapper(CombinatorialClass):
    def __init__(self, species, labels, iterator, generating_series, name, structure_class):
        """
        EXAMPLES::

            sage: from sage.combinat.species.structure import SpeciesWrapper
            sage: F = species.SetSpecies()
            sage: S = SpeciesWrapper(F, [1,2,3], "_structures", "generating_series", 'Structures', None)
            sage: S
            Structures for Set species with labels [1, 2, 3]
            sage: S.list()
            [{1, 2, 3}]
            sage: S.cardinality()
            1
        """
        self._species = species
        self._labels = labels
        self._iterator = iterator
        self._generating_series = generating_series
        self._name = "%s for %s with labels %s"%(name, species, labels)
        self._structure_class = structure_class if structure_class is not None else species._default_structure_class

    def __iter__(self):
        """
        EXAMPLES::

            sage: F = species.SetSpecies()
            sage: F.structures([1,2,3]).list()
            [{1, 2, 3}]
        """
        #If the min and max are set, then we want to make sure
        #that the iterator respects those bounds.
        if (self._species._min is not None and
            len(self._labels) < self._species._min):
            return iter([])

        if (self._species._max is not None and
            len(self._labels) >= self._species._max):
            return iter([])

        #We check to see if the
        try:
            if self.cardinality() == 0:
                return iter([])
        except RuntimeError:
            raise NotImplementedError

        return getattr(self._species, self._iterator)(self._structure_class, self._labels)

    def cardinality(self):
        """
        EXAMPLES::

            sage: F = species.SetSpecies()
            sage: F.structures([1,2,3]).cardinality()
            1
        """
        return getattr(self._species, self._generating_series)().count(len(self._labels))

class StructuresWrapper(SpeciesWrapper):
    def __init__(self, species, labels, structure_class):
        """
        EXAMPLES::

            sage: F = species.SetSpecies()
            sage: S = F.structures([1,2,3])
            sage: S == loads(dumps(S))
            True
        """
        SpeciesWrapper.__init__(self, species, labels,
                                "_structures",
                                "generating_series",
                                "Structures",
                                structure_class)

class IsotypesWrapper(SpeciesWrapper):
    def __init__(self, species, labels, structure_class):
        """
        EXAMPLES::

            sage: F = species.SetSpecies()
            sage: S = F.isotypes([1,2,3])
            sage: S == loads(dumps(S))
            True
        """
        SpeciesWrapper.__init__(self, species, labels,
                                "_isotypes",
                                "isotype_generating_series",
                                "Isomorphism types",
                                structure_class)


class SimpleStructuresWrapper(SpeciesWrapper):
    def __init__(self, species, labels, structure_class):
        """
        EXAMPLES::

            sage: F = species.SetSpecies()
            sage: S = F.structures([1,2,3])
            sage: S == loads(dumps(S))
            True
        """
        SpeciesWrapper.__init__(self, species, labels,
                                "_simple_structures_selector",
                                "generating_series",
                                "Simple structures",
                                structure_class)

class SimpleIsotypesWrapper(SpeciesWrapper):
    def __init__(self, species, labels, structure_class):
        """
        EXAMPLES::

            sage: F = species.SetSpecies()
            sage: S = F.structures([1,2,3])
            sage: S == loads(dumps(S))
            True
        """
        SpeciesWrapper.__init__(self, species, labels,
                                "_simple_isotypes_selector",
                                "isotype_generating_series",
                                "Simple isomorphism types",
                                structure_class)
