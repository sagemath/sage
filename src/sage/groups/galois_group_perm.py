# sage_setup: distribution = sagemath-gap
r"""
Galois groups of field extensions as permutation groups
"""
from sage.groups.galois_group import _GaloisMixin, _SubGaloisMixin
from sage.groups.perm_gps.permgroup import PermutationGroup_generic, PermutationGroup_subgroup
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.structure.category_object import normalize_names


class GaloisGroup_perm(_GaloisMixin, PermutationGroup_generic):
    r"""
    The group of automorphisms of a Galois closure of a given field.

    INPUT:

    - ``field`` -- a field, separable over its base

    - ``names`` -- string or tuple of length 1, giving a variable name for
      the splitting field

    - ``gc_numbering`` -- boolean, whether to express permutations in terms of
      the roots of the defining polynomial of the splitting field (versus the
      defining polynomial of the original extension); the default value may
      vary based on the type of field
    """
    @abstract_method
    def transitive_number(self, algorithm=None, recompute=False):
        """
        Return the transitive number (as in the GAP and Magma databases of transitive groups)
        for the action on the roots of the defining polynomial of the top field.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)                                    # needs sage.rings.number_field
            sage: G = K.galois_group()                                                  # needs sage.rings.number_field
            sage: G.transitive_number()                                                 # needs sage.rings.number_field
            2
        """

    @lazy_attribute
    def _gens(self):
        """
        The generators of this Galois group as permutations of the roots.

        It is important that this be computed lazily, since it is
        often possible to compute other attributes (such as the order
        or transitive number) more cheaply.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5 - 2)                                          # needs sage.rings.number_field
            sage: G = K.galois_group(gc_numbering=False)                                # needs sage.rings.number_field
            sage: G._gens                                                               # needs sage.rings.number_field
            [(1,2,3,5), (1,4,3,2,5)]
        """
        return NotImplemented

    def __init__(self, field, algorithm=None, names=None, gc_numbering=False):
        r"""
        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)                                    # needs sage.rings.number_field
            sage: G = K.galois_group()                                                  # needs sage.rings.number_field
            sage: TestSuite(G).run()                                                    # needs sage.rings.number_field
        """
        self._field = field
        self._default_algorithm = algorithm
        self._gc_numbering = gc_numbering
        if names is None:
            # add a c for Galois closure
            names = field.variable_name() + 'c'
        self._gc_names = normalize_names(1, names)
        # We do only the parts of the initialization of PermutationGroup_generic
        # that don't depend on _gens
        from sage.categories.permutation_groups import PermutationGroups
        category = PermutationGroups().FinitelyGenerated().Finite()
        # Note that we DON'T call the __init__ method for PermutationGroup_generic
        # Instead, the relevant attributes are computed lazily
        super(PermutationGroup_generic, self).__init__(category=category)

    @lazy_attribute
    def _deg(self):
        r"""
        The number of moved points in the permutation representation.

        This will be the degree of the original number field if ``_gc_numbering``
        is ``False``, or the degree of the Galois closure otherwise.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5 - 2)
            sage: G = K.galois_group(gc_numbering=False); G
            Galois group 5T3 (5:4) with order 20 of x^5 - 2
            sage: G._deg
            5
            sage: G = K.galois_group(gc_numbering=True); G._deg
            20
        """
        if self._gc_numbering:
            return self.order()

        try:
            return self._field.degree()
        except NotImplementedError:  # relative number fields don't support degree
            return self._field.relative_degree()

    @lazy_attribute
    def _domain(self):
        r"""
        The integers labeling the roots on which this Galois group acts.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5 - 2)
            sage: G = K.galois_group(gc_numbering=False); G
            Galois group 5T3 (5:4) with order 20 of x^5 - 2
            sage: G._domain
            {1, 2, 3, 4, 5}
            sage: G = K.galois_group(gc_numbering=True); G._domain
            {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20}
        """
        return FiniteEnumeratedSet(range(1, self._deg + 1))

    @lazy_attribute
    def _domain_to_gap(self) -> dict:
        r"""
        The dictionary implementing the identity (used by PermutationGroup_generic).

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5 - 2)                                          # needs sage.rings.number_field
            sage: G = K.galois_group(gc_numbering=False)                                # needs sage.rings.number_field
            sage: G._domain_to_gap[5]                                                   # needs sage.rings.number_field
            5
        """
        return {key: i + 1 for i, key in enumerate(self._domain)}

    @lazy_attribute
    def _domain_from_gap(self) -> dict:
        r"""
        The dictionary implementing the identity (used by PermutationGroup_generic).

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5 - 2)                                          # needs sage.rings.number_field
            sage: G = K.galois_group(gc_numbering=True)                                 # needs sage.rings.number_field
            sage: G._domain_from_gap[20]                                                # needs sage.rings.number_field
            20
        """
        return {i + 1: key for i, key in enumerate(self._domain)}

    def ngens(self) -> int:
        r"""
        Return the number of generators of this Galois group.

        EXAMPLES::

            sage: QuadraticField(-23, 'a').galois_group().ngens()                       # needs sage.rings.number_field
            1
        """
        return len(self._gens)


class GaloisSubgroup_perm(PermutationGroup_subgroup, _SubGaloisMixin):
    """
    Subgroups of Galois groups (implemented as permutation groups), specified
    by giving a list of generators.

    Unlike ambient Galois groups, where we use a lazy ``_gens`` attribute in order
    to enable creation without determining a list of generators,
    we require that generators for a subgroup be specified during initialization,
    as specified in the ``__init__`` method of permutation subgroups.
    """
    pass


GaloisGroup_perm.Subgroup = GaloisSubgroup_perm
