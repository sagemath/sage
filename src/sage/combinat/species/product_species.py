"""
Sum species
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
from .subset_species import SubsetSpecies
from sage.structure.unique_representation import UniqueRepresentation


class ProductSpeciesStructure(GenericSpeciesStructure):
    def __init__(self, parent, labels, subset, left, right):
        """
        TESTS::

            sage: S = species.SetSpecies()
            sage: F = S * S
            sage: a = F.structures(['a','b','c']).random_element()
            sage: a == loads(dumps(a))
            True
        """
        self._subset = subset
        GenericSpeciesStructure.__init__(self, parent, labels, [left, right])

    def __repr__(self):
        """
        Return the string representation of this object.

        EXAMPLES::

            sage: S = species.SetSpecies()
            sage: (S*S).structures(['a','b','c'])[0]
            {}*{'a', 'b', 'c'}
            sage: (S*S*S).structures(['a','b','c'])[13]
            ({'c'}*{'a'})*{'b'}
        """
        left, right = map(repr, self._list)
        if "*" in left:
            left = "(%s)" % left
        if "*" in right:
            right = "(%s)" % right
        return "%s*%s" % (left, right)

    def transport(self, perm):
        """
        EXAMPLES::

            sage: p = PermutationGroupElement((2,3))
            sage: S = species.SetSpecies()
            sage: F = S * S
            sage: a = F.structures(['a','b','c'])[4]; a
            {'a', 'b'}*{'c'}
            sage: a.transport(p)
            {'a', 'c'}*{'b'}
        """
        left, right = self._list
        new_subset = self._subset.transport(perm)
        left_labels = new_subset.label_subset()
        right_labels = new_subset.complement().label_subset()

        return self.__class__(self.parent(), self._labels,
                              new_subset,
                              left.change_labels(left_labels),
                              right.change_labels(right_labels))

    def canonical_label(self):
        """
        EXAMPLES::

            sage: S = species.SetSpecies()
            sage: F = S * S
            sage: S = F.structures(['a','b','c']).list(); S
            [{}*{'a', 'b', 'c'},
             {'a'}*{'b', 'c'},
             {'b'}*{'a', 'c'},
             {'c'}*{'a', 'b'},
             {'a', 'b'}*{'c'},
             {'a', 'c'}*{'b'},
             {'b', 'c'}*{'a'},
             {'a', 'b', 'c'}*{}]

        ::

            sage: F.isotypes(['a','b','c']).cardinality()
            4
            sage: [s.canonical_label() for s in S]
            [{}*{'a', 'b', 'c'},
             {'a'}*{'b', 'c'},
             {'a'}*{'b', 'c'},
             {'a'}*{'b', 'c'},
             {'a', 'b'}*{'c'},
             {'a', 'b'}*{'c'},
             {'a', 'b'}*{'c'},
             {'a', 'b', 'c'}*{}]
        """
        left, right = self._list
        new_subset = self._subset.canonical_label()
        left_labels = new_subset.label_subset()
        right_labels = new_subset.complement().label_subset()

        return self.__class__(self.parent(), self._labels,
                              new_subset,
                              left.canonical_label().change_labels(left_labels),
                              right.canonical_label().change_labels(right_labels))

    def change_labels(self, labels):
        """
        Return a relabelled structure.

        INPUT:

        - ``labels``, a list of labels.

        OUTPUT:

        A structure with the i-th label of self replaced with the i-th
        label of the list.

        EXAMPLES::

            sage: S = species.SetSpecies()
            sage: F = S * S
            sage: a = F.structures(['a','b','c'])[0]; a
            {}*{'a', 'b', 'c'}
            sage: a.change_labels([1,2,3])
            {}*{1, 2, 3}
        """
        left, right = self._list
        new_subset = self._subset.change_labels(labels)
        left_labels = new_subset.label_subset()
        right_labels = new_subset.complement().label_subset()
        return self.__class__(self.parent(), labels,
                              new_subset,
                              left.change_labels(left_labels),
                              right.change_labels(right_labels))

    def automorphism_group(self):
        """
        EXAMPLES::

            sage: p = PermutationGroupElement((2,3))
            sage: S = species.SetSpecies()
            sage: F = S * S
            sage: a = F.structures([1,2,3,4])[1]; a
            {1}*{2, 3, 4}
            sage: a.automorphism_group()
            Permutation Group with generators [(2,3), (2,3,4)]

        ::

            sage: [a.transport(g) for g in a.automorphism_group()]
            [{1}*{2, 3, 4},
             {1}*{2, 3, 4},
             {1}*{2, 3, 4},
             {1}*{2, 3, 4},
             {1}*{2, 3, 4},
             {1}*{2, 3, 4}]

        ::

            sage: a = F.structures([1,2,3,4])[8]; a
            {2, 3}*{1, 4}
            sage: [a.transport(g) for g in a.automorphism_group()]
            [{2, 3}*{1, 4}, {2, 3}*{1, 4}, {2, 3}*{1, 4}, {2, 3}*{1, 4}]
        """
        from sage.groups.all import PermutationGroupElement, PermutationGroup
        from sage.combinat.species.misc import change_support

        left, right = self._list

        #Get the supports for each of the sides
        l_support = self._subset._list
        r_support = self._subset.complement()._list

        #Get the automorphism group for the left object and
        #make it have the correct support. Do the same to the
        #right side.
        l_aut = change_support(left.automorphism_group(), l_support)
        r_aut = change_support(right.automorphism_group(), r_support)

        identity = PermutationGroupElement([])

        gens = l_aut.gens() + r_aut.gens()
        gens = [g for g in gens if g != identity]
        gens = sorted(set(gens)) if gens else [[]]
        return PermutationGroup(gens)


class ProductSpecies(GenericCombinatorialSpecies, UniqueRepresentation):
    def __init__(self, F, G, min=None, max=None, weight=None):
        """
        EXAMPLES::

            sage: X = species.SingletonSpecies()
            sage: A = X*X
            sage: A.generating_series().coefficients(4)
            [0, 0, 1, 0]

            sage: P = species.PermutationSpecies()
            sage: F = P * P; F
            Product of (Permutation species) and (Permutation species)
            sage: F == loads(dumps(F))
            True
            sage: F._check()
            True

        TESTS::

            sage: X = species.SingletonSpecies()
            sage: X*X is X*X
            True
        """
        self._F = F
        self._G = G
        self._state_info = [F, G]
        GenericCombinatorialSpecies.__init__(self, min=None, max=None, weight=weight)


    _default_structure_class = ProductSpeciesStructure

    def left_factor(self):
        """
        Returns the left factor of this product.

        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: X = species.SingletonSpecies()
            sage: F = P*X
            sage: F.left_factor()
            Permutation species
        """
        return self._F

    def right_factor(self):
        """
        Returns the right factor of this product.

        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: X = species.SingletonSpecies()
            sage: F = P*X
            sage: F.right_factor()
            Singleton species
        """
        return self._G

    def _name(self):
        """
        Note that we use a function to return the name of this species
        because we can't do it in the __init__ method due to it
        requiring that self.left_factor() and self.right_factor()
        already be unpickled.

        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: F = P * P
            sage: F._name()
            'Product of (Permutation species) and (Permutation species)'
        """
        return "Product of (%s) and (%s)"%(self.left_factor(), self.right_factor())

    def _structures(self, structure_class, labels):
        """
        EXAMPLES::

            sage: S = species.SetSpecies()
            sage: F = S * S
            sage: F.structures([1,2]).list()
            [{}*{1, 2}, {1}*{2}, {2}*{1}, {1, 2}*{}]
        """
        return self._times_gen(structure_class, "structures", labels)

    def _isotypes(self, structure_class, labels):
        """
        EXAMPLES::

            sage: S = species.SetSpecies()
            sage: F = S * S
            sage: F.isotypes([1,2,3]).list()
            [{}*{1, 2, 3}, {1}*{2, 3}, {1, 2}*{3}, {1, 2, 3}*{}]
        """
        return self._times_gen(structure_class, "isotypes", labels)

    def _times_gen(self, structure_class, attr, labels):
        """
        EXAMPLES::

            sage: S = species.SetSpecies()
            sage: F = S * S
            sage: list(F._times_gen(F._default_structure_class, 'structures',[1,2]))
            [{}*{1, 2}, {1}*{2}, {2}*{1}, {1, 2}*{}]
        """
        c = lambda F,n: F.generating_series().coefficient(n)
        S = SubsetSpecies()

        for u in getattr(S, attr)(labels):
            vl = u.complement().label_subset()
            ul = u.label_subset()
            if c(self.left_factor(), len(ul)) == 0 or c(self.right_factor(), len(vl)) == 0:
                continue
            for x in getattr(self.left_factor(), attr)(ul):
                for y in getattr(self.right_factor(), attr)(vl):
                    yield structure_class(self, labels, u, x, y)

    def _gs(self, series_ring, base_ring):
        """
        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: F = P * P
            sage: F.generating_series().coefficients(5)
            [1, 2, 3, 4, 5]
        """
        res = (self.left_factor().generating_series(base_ring) *
               self.right_factor().generating_series(base_ring))
        if self.is_weighted():
            res = self._weight * res
        return res

    def _itgs(self, series_ring, base_ring):
        """
        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: F = P * P
            sage: F.isotype_generating_series().coefficients(5)
            [1, 2, 5, 10, 20]
        """
        res =  (self.left_factor().isotype_generating_series(base_ring) *
                self.right_factor().isotype_generating_series(base_ring))
        if self.is_weighted():
            res = self._weight * res
        return res

    def _cis(self, series_ring, base_ring):
        """
        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: F = P * P
            sage: F.cycle_index_series().coefficients(5)
            [p[],
             2*p[1],
             3*p[1, 1] + 2*p[2],
             4*p[1, 1, 1] + 4*p[2, 1] + 2*p[3],
             5*p[1, 1, 1, 1] + 6*p[2, 1, 1] + 3*p[2, 2] + 4*p[3, 1] + 2*p[4]]
        """
        res =  (self.left_factor().cycle_index_series(base_ring) *
                self.right_factor().cycle_index_series(base_ring))
        if self.is_weighted():
            res = self._weight * res
        return res

    def weight_ring(self):
        """
        Returns the weight ring for this species. This is determined by
        asking Sage's coercion model what the result is when you multiply
        (and add) elements of the weight rings for each of the operands.

        EXAMPLES::

            sage: S = species.SetSpecies()
            sage: C = S*S
            sage: C.weight_ring()
            Rational Field

        ::

            sage: S = species.SetSpecies(weight=QQ['t'].gen())
            sage: C = S*S
            sage: C.weight_ring()
            Univariate Polynomial Ring in t over Rational Field

        ::

            sage: S = species.SetSpecies()
            sage: C = (S*S).weighted(QQ['t'].gen())
            sage: C.weight_ring()
            Univariate Polynomial Ring in t over Rational Field
        """
        return self._common_parent([self.left_factor().weight_ring(),
                                    self.right_factor().weight_ring(),
                                    self._weight.parent()])

    def _equation(self, var_mapping):
        """
        Returns the right hand side of an algebraic equation satisfied by
        this species. This is a utility function called by the
        algebraic_equation_system method.

        EXAMPLES::

            sage: X = species.SingletonSpecies()
            sage: S = X * X
            sage: S.algebraic_equation_system()
            [node0 + (-z^2)]
        """
        from sage.misc.misc_c import prod
        return prod(var_mapping[operand] for operand in self._state_info)


#Backward compatibility
ProductSpecies_class = ProductSpecies
