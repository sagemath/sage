"""
Recursive Species
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
from sage.combinat.species.species import GenericCombinatorialSpecies
from sage.combinat.species.structure import SpeciesStructureWrapper
from sage.rings.rational_field import QQ


class CombinatorialSpeciesStructure(SpeciesStructureWrapper):
    pass


class CombinatorialSpecies(GenericCombinatorialSpecies):
    def __init__(self):
        """
        EXAMPLES::

            sage: F = CombinatorialSpecies()
            sage: loads(dumps(F))
            Combinatorial species

        ::

            sage: X = species.SingletonSpecies()
            sage: E = species.EmptySetSpecies()
            sage: L = CombinatorialSpecies()
            sage: L.define(E+X*L)
            sage: L.generating_series().coefficients(4)
            [1, 1, 1, 1]
            sage: LL = loads(dumps(L))
            sage: LL.generating_series().coefficients(4)
            [1, 1, 1, 1]
        """
        self._generating_series = {}
        self._isotype_generating_series = {}
        self._cycle_index_series = {}
        self._min = None
        self._max = None
        self._weight = 1
        GenericCombinatorialSpecies.__init__(self, min=None, max=None, weight=None)

    _default_structure_class = CombinatorialSpeciesStructure

    def __hash__(self):
        """
        EXAMPLES::

            sage: hash(CombinatorialSpecies) #random
            53751280

        ::

            sage: X = species.SingletonSpecies()
            sage: E = species.EmptySetSpecies()
            sage: L = CombinatorialSpecies()
            sage: L.define(E+X*L)
            sage: hash(L) #random
            -826511807095108317
        """
        try:
            return hash(('CombinatorialSpecies', id(self._reference)))
        except AttributeError:
            return hash('CombinatorialSpecies')

    def __eq__(self, other):
        """
        TESTS::

            sage: A = species.CombinatorialSpecies()
            sage: B = species.CombinatorialSpecies()
            sage: A == B
            False
            sage: X = species.SingletonSpecies()
            sage: A.define(X+A*A)
            sage: B.define(X+B*B)
            sage: A == B
            True

            sage: C = species.CombinatorialSpecies()
            sage: E = species.EmptySetSpecies()
            sage: C.define(E+X*C*C)
            sage: A == C
            False
        """
        if not isinstance(other, CombinatorialSpecies):
            return False
        if not hasattr(self, "_reference"):
            return False
        if hasattr(self, '_computing_eq'):
            return True

        self._computing_eq = True
        res = self._unique_info() == other._unique_info()
        del self._computing_eq
        return res

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: A = species.CombinatorialSpecies()
            sage: B = species.CombinatorialSpecies()
            sage: A != B
            True
            sage: X = species.SingletonSpecies()
            sage: A.define(X+A*A)
            sage: B.define(X+B*B)
            sage: A != B
            False

            sage: C = species.CombinatorialSpecies()
            sage: E = species.EmptySetSpecies()
            sage: C.define(E+X*C*C)
            sage: A != C
            True
        """
        return not (self == other)

    def _unique_info(self):
        """
        Return a tuple which should uniquely identify the species.

        EXAMPLES::

            sage: F = CombinatorialSpecies()
            sage: F._unique_info()
            (<class 'sage.combinat.species.recursive_species.CombinatorialSpecies'>,)

        ::

            sage: X = species.SingletonSpecies()
            sage: E = species.EmptySetSpecies()
            sage: L = CombinatorialSpecies()
            sage: L.define(E+X*L)
            sage: L._unique_info()
            (<class 'sage.combinat.species.recursive_species.CombinatorialSpecies'>,
             <class 'sage.combinat.species.sum_species.SumSpecies'>,
             None,
             None,
             1,
             Empty set species,
             Product of (Singleton species) and (Combinatorial species))
        """
        if hasattr(self, "_reference"):
            return (self.__class__,) + self._reference._unique_info()
        else:
            return (self.__class__,)

    def __getstate__(self):
        """
        EXAMPLES::

            sage: X = species.SingletonSpecies()
            sage: E = species.EmptySetSpecies()
            sage: L = CombinatorialSpecies()
            sage: L.define(E+X*L)
            sage: L.__getstate__()
            {'reference': Sum of (Empty set species) and (Product of (Singleton species) and (Combinatorial species))}
        """
        state = {}
        if hasattr(self, '_reference'):
            state['reference'] = self._reference
        return state

    def __setstate__(self, state):
        """
        EXAMPLES::

            sage: X = species.SingletonSpecies()
            sage: E = species.EmptySetSpecies()
            sage: L = CombinatorialSpecies()
            sage: L.define(E+X*L)
            sage: state = L.__getstate__(); state
            {'reference': Sum of (Empty set species) and (Product of (Singleton species) and (Combinatorial species))}
            sage: L._reference = None
            sage: L.__setstate__(state)
            sage: L._reference
            Sum of (Empty set species) and (Product of (Singleton species) and (Combinatorial species))
        """
        CombinatorialSpecies.__init__(self)
        if 'reference' in state:
            self.define(state['reference'])

    def _structures(self, structure_class, labels):
        """
        EXAMPLES::

            sage: F = CombinatorialSpecies()
            sage: list(F._structures(F._default_structure_class, [1,2,3]))
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if not hasattr(self, "_reference"):
            raise NotImplementedError
        for s in self._reference.structures(labels):
            yield structure_class(self, s)

    def _isotypes(self, structure_class, labels):
        """
        EXAMPLES::

            sage: F = CombinatorialSpecies()
            sage: list(F._isotypes(F._default_structure_class, [1,2,3]))
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if not hasattr(self, "_reference"):
            raise NotImplementedError
        for s in self._reference.isotypes(labels):
            yield structure_class(self, s)

    def _gs(self, series_ring, base_ring):
        """
        EXAMPLES::

            sage: F = CombinatorialSpecies()
            sage: F.generating_series()
            Uninitialized lazy power series
        """
        if base_ring not in self._generating_series:
            self._generating_series[base_ring] = series_ring()

        res = self._generating_series[base_ring]
        if hasattr(self, "_reference") and not hasattr(res, "_reference"):
            res._reference = None
            res.define(self._reference.generating_series(base_ring))
        return res

    def _itgs(self, series_ring, base_ring):
        """
        EXAMPLES::

            sage: F = CombinatorialSpecies()
            sage: F.isotype_generating_series()
            Uninitialized lazy power series
        """
        if base_ring not in self._isotype_generating_series:
            self._isotype_generating_series[base_ring] = series_ring()

        res = self._isotype_generating_series[base_ring]
        if hasattr(self, "_reference") and not hasattr(res, "_reference"):
            res._reference = None
            res.define(self._reference.isotype_generating_series(base_ring))
        return res

    def _cis(self, series_ring, base_ring):
        """
        EXAMPLES::

            sage: F = CombinatorialSpecies()
            sage: F.cycle_index_series()
            Uninitialized lazy power series
        """
        if base_ring not in self._cycle_index_series:
            self._cycle_index_series[base_ring] = series_ring()

        res = self._cycle_index_series[base_ring]
        if hasattr(self, "_reference") and not hasattr(res, "_reference"):
            res._reference = None
            res.define(self._reference.cycle_index_series(base_ring))
        return res

    def weight_ring(self):
        """
        EXAMPLES::

            sage: F = species.CombinatorialSpecies()
            sage: F.weight_ring()
            Rational Field

        ::

            sage: X = species.SingletonSpecies()
            sage: E = species.EmptySetSpecies()
            sage: L = CombinatorialSpecies()
            sage: L.define(E+X*L)
            sage: L.weight_ring()
            Rational Field
        """
        if not hasattr(self, "_reference"):
            return QQ

        if hasattr(self, "_weight_ring_been_called"):
            return QQ
        else:
            self._weight_ring_been_called = True
            res = self._reference.weight_ring()
            del self._weight_ring_been_called
            return res

    def define(self, x):
        """
        Define ``self`` to be equal to the combinatorial species ``x``.

        This is
        used to define combinatorial species recursively. All of the real
        work is done by calling the .set() method for each of the series
        associated to self.

        EXAMPLES: The species of linear orders L can be recursively defined
        by `L = 1 + X*L` where 1 represents the empty set species
        and X represents the singleton species.

        ::

            sage: X = species.SingletonSpecies()
            sage: E = species.EmptySetSpecies()
            sage: L = CombinatorialSpecies()
            sage: L.define(E+X*L)
            sage: L.generating_series().coefficients(4)
            [1, 1, 1, 1]
            sage: L.structures([1,2,3]).cardinality()
            6
            sage: L.structures([1,2,3]).list()
            [1*(2*(3*{})),
             1*(3*(2*{})),
             2*(1*(3*{})),
             2*(3*(1*{})),
             3*(1*(2*{})),
             3*(2*(1*{}))]

        ::

            sage: L = species.LinearOrderSpecies()
            sage: L.generating_series().coefficients(4)
            [1, 1, 1, 1]
            sage: L.structures([1,2,3]).cardinality()
            6
            sage: L.structures([1,2,3]).list()
            [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]

        TESTS::

            sage: A = CombinatorialSpecies()
            sage: A.define(E+X*A*A)
            sage: A.generating_series().coefficients(6)
            [1, 1, 2, 5, 14, 42]
            sage: A.generating_series().counts(6)
            [1, 1, 4, 30, 336, 5040]
            sage: len(A.structures([1,2,3,4]).list())
            336
            sage: A.isotype_generating_series().coefficients(6)
            [1, 1, 2, 5, 14, 42]
            sage: len(A.isotypes([1,2,3,4]).list())
            14

        ::

            sage: A = CombinatorialSpecies()
            sage: A.define(X+A*A)
            sage: A.generating_series().coefficients(6)
            [0, 1, 1, 2, 5, 14]
            sage: A.generating_series().counts(6)
            [0, 1, 2, 12, 120, 1680]
            sage: len(A.structures([1,2,3]).list())
            12
            sage: A.isotype_generating_series().coefficients(6)
            [0, 1, 1, 2, 5, 14]
            sage: len(A.isotypes([1,2,3,4]).list())
            5

        ::

            sage: X2 = X*X
            sage: X5 = X2*X2*X
            sage: A = CombinatorialSpecies()
            sage: B = CombinatorialSpecies()
            sage: C = CombinatorialSpecies()
            sage: A.define(X5+B*B)
            sage: B.define(X5+C*C)
            sage: C.define(X2+C*C+A*A)
            sage: A.generating_series().coefficients(Integer(10))
            [0, 0, 0, 0, 0, 1, 0, 0, 1, 2]
            sage: A.generating_series().coefficients(15)
            [0, 0, 0, 0, 0, 1, 0, 0, 1, 2, 5, 4, 14, 10, 48]
            sage: B.generating_series().coefficients(15)
            [0, 0, 0, 0, 1, 1, 2, 0, 5, 0, 14, 0, 44, 0, 138]
            sage: C.generating_series().coefficients(15)
            [0, 0, 1, 0, 1, 0, 2, 0, 5, 0, 15, 0, 44, 2, 142]

        ::

            sage: F = CombinatorialSpecies()
            sage: F.define(E+X+(X*F+X*X*F))
            sage: F.generating_series().counts(10)
            [1, 2, 6, 30, 192, 1560, 15120, 171360, 2217600, 32296320]
            sage: F.generating_series().coefficients(10)
            [1, 2, 3, 5, 8, 13, 21, 34, 55, 89]
            sage: F.isotype_generating_series().coefficients(10)
            [1, 2, 3, 5, 8, 13, 21, 34, 55, 89]
        """
        if not isinstance(x, GenericCombinatorialSpecies):
            raise TypeError("x must be a combinatorial species")

        if self.__class__ is not CombinatorialSpecies:
            raise TypeError("only undefined combinatorial species can be set")

        self._reference = x


    def _add_to_digraph(self, d):
        """
        Adds this species as a vertex to the digraph d along with any
        'children' of this species.

        Note that to avoid infinite recursion, we just return if this
        species already occurs in the digraph d.

        EXAMPLES::

            sage: d = DiGraph(multiedges=True)
            sage: X = species.SingletonSpecies()
            sage: B = species.CombinatorialSpecies()
            sage: B.define(X+B*B)
            sage: B._add_to_digraph(d); d
            Multi-digraph on 4 vertices

        TESTS::

            sage: C = species.CombinatorialSpecies()
            sage: C._add_to_digraph(d)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if self in d:
            return
        try:
            d.add_edge(self, self._reference)
            self._reference._add_to_digraph(d)
        except AttributeError:
            raise NotImplementedError

    def _equation(self, var_mapping):
        """
        Returns the right hand side of an algebraic equation satisfied by
        this species. This is a utility function called by the
        algebraic_equation_system method.

        EXAMPLES::

            sage: C = species.CombinatorialSpecies()
            sage: C.algebraic_equation_system()
            Traceback (most recent call last):
            ...
            NotImplementedError

        ::

            sage: B = species.BinaryTreeSpecies()
            sage: B.algebraic_equation_system()
            [-node3^2 + node1, -node1 + node3 + (-z)]
        """
        try:
            return var_mapping[self._reference]
        except AttributeError:
            raise NotImplementedError
