r"""
Galois groups of field extensions

AUTHORS:

- David Roe (2019): initial version
"""

from sage.groups.perm_gps.permgroup import PermutationGroup_generic
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.category_object import normalize_names

def _alg_key(self, algorithm=None, recompute=False):
    r"""
    Return a key for use in cached_method calls.

    If recompute is false, will cache using ``None`` as the key, so no recomputation will be done.

    If recompute is true, will cache by algorithm, yielding a recomputation for each different algorithm.

    EXAMPLES::

        sage: from sage.groups.galois_group import _alg_key
        sage: R.<x> = ZZ[]
        sage: K.<a> = NumberField(x^3 + 2*x + 2)
        sage: G = K.galois_group()
        sage: _alg_key(G, algorithm="pari", recompute=True)
        'pari'
    """
    if recompute:
        algorithm = self._get_algorithm(algorithm)
        return algorithm

class GaloisGroup(PermutationGroup_generic):
    r"""
    The group of automorphisms of a Galois closure of a given field.

    INPUT:

    - ``field`` -- a field, separable over its base

    - ``names`` -- a string or tuple of length 1, giving a variable name for the splitting field

    - ``gc_numbering`` -- boolean, whether to express permutations in terms of the
        roots of the defining polynomial of the splitting field (versus the defining polynomial
        of the original extension).  The default value may vary based on the type of field.
    """
    # Subclasses should implement the following methods and lazy attributes

    # methods (taking algorithm and recompute as arguments):
    # * transitive_number
    # * order
    # * _element_constructor_ -- for creating elements

    # lazy_attributes
    # * _gcdata -- a pair, the Galois closure and an embedding of the top field into it
    # * _gens -- the list of generators of this group, as elements.  This is not computed during __init__ for speed
    # * _elts -- the list of all elements of this group.

    # * Element (for coercion)

    def __init__(self, field, algorithm=None, names=None, gc_numbering=False):
        r"""
        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)
            sage: G = K.galois_group()
            sage: TestSuite(G).run()
        """
        self._field = field
        self._default_algorithm = algorithm
        self._base = field.base_field()
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

    def _repr_(self):
        """
        String representation of this Galois group

        EXAMPLES::

            sage: from sage.groups.galois_group import GaloisGroup
            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)
            sage: G = K.galois_group()
            sage: GaloisGroup._repr_(G)
            'Galois group of x^3 + 2*x + 2'
        """
        f = self._field.defining_polynomial()
        return "Galois group of %s" % f

    def _get_algorithm(self, algorithm):
        r"""
        Allows overriding the default algorithm specified at object creation.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)
            sage: G = K.galois_group()
            sage: G._get_algorithm(None)
            'pari'
            sage: G._get_algorithm('magma')
            'magma'
        """
        return self._default_algorithm if algorithm is None else algorithm

    def top_field(self):
        r"""
        Return the larger of the two fields in the extension defining this Galois group.

        Note that this field may not be Galois.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)
            sage: L = K.galois_closure('b')
            sage: GK = K.galois_group()
            sage: GK.top_field() is K
            True
            sage: GL = L.galois_group()
            sage: GL.top_field() is L
            True
        """
        return self._field

    def transitive_label(self):
        r"""
        Return the transitive label for the action of this Galois group on the roots of
        the defining polynomial of the field extension.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^8 - x^5 + x^4 - x^3 + 1)
            sage: G = K.galois_group()
            sage: G.transitive_label()
            '8T44'
        """
        try:
            return "%sT%s" % (self._field.degree(), self.transitive_number())
        except NotImplementedError: # relative number fields don't support degree
            return "%sT%s" % (self._field.relative_degree(), self.transitive_number())

    def is_galois(self):
        r"""
        Return whether the top field is Galois over its base.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^8 - x^5 + x^4 - x^3 + 1)
            sage: G = K.galois_group()
            sage: from sage.groups.galois_group import GaloisGroup
            sage: GaloisGroup.is_galois(G)
            False
        """
        return self.order() == self._field.degree()

    @lazy_attribute
    def _galois_closure(self):
        r"""
        The Galois closure of the top field.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)
            sage: G = K.galois_group(names='b')
            sage: G._galois_closure
            Number Field in b with defining polynomial x^6 + 12*x^4 + 36*x^2 + 140
        """
        return self._gcdata[0]

    def splitting_field(self):
        r"""
        The Galois closure of the top field.

        EXAMPLES::

            sage: K = NumberField(x^3 - x + 1, 'a')
            sage: K.galois_group(names='b').splitting_field()
            Number Field in b with defining polynomial x^6 - 6*x^4 + 9*x^2 + 23
            sage: L = QuadraticField(-23, 'c'); L.galois_group().splitting_field() is L
            True
        """
        return self._galois_closure

    @lazy_attribute
    def _gc_map(self):
        r"""
        The inclusion of the top field into the Galois closure.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)
            sage: G = K.galois_group(names='b')
            sage: G._gc_map
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 + 2*x + 2
              To:   Number Field in b with defining polynomial x^6 + 12*x^4 + 36*x^2 + 140
              Defn: a |--> 1/36*b^4 + 5/18*b^2 - 1/2*b + 4/9
        """
        return self._gcdata[1]

    @lazy_attribute
    def _deg(self):
        r"""
        The number of moved points in the permutation representation.

        This will be the degree of the original number field if `_gc_numbering``
        is ``False``, or the degree of the Galois closure otherwise.

        EXAMPES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5-2)
            sage: G = K.galois_group(gc_numbering=False); G
            Galois group 5T3 (5:4) with order 20 of x^5 - 2
            sage: G._deg
            5
            sage: G = K.galois_group(gc_numbering=True); G._deg
            20
        """
        if self._gc_numbering:
            return self.order()
        else:
            try:
                return self._field.degree()
            except NotImplementedError: # relative number fields don't support degree
                return self._field.relative_degree()

    @lazy_attribute
    def _domain(self):
        r"""
        The integers labeling the roots on which this Galois group acts.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5-2)
            sage: G = K.galois_group(gc_numbering=False); G
            Galois group 5T3 (5:4) with order 20 of x^5 - 2
            sage: G._domain
            {1, 2, 3, 4, 5}
            sage: G = K.galois_group(gc_numbering=True); G._domain
            {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20}
        """
        return FiniteEnumeratedSet(range(1, self._deg+1))

    @lazy_attribute
    def _domain_to_gap(self):
        r"""
        Dictionary implementing the identity (used by PermutationGroup_generic).

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5-2)
            sage: G = K.galois_group(gc_numbering=False)
            sage: G._domain_to_gap[5]
            5
        """
        return dict((key, i+1) for i, key in enumerate(self._domain))

    @lazy_attribute
    def _domain_from_gap(self):
        r"""
        Dictionary implementing the identity (used by PermutationGroup_generic).

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5-2)
            sage: G = K.galois_group(gc_numbering=True)
            sage: G._domain_from_gap[20]
            20
        """
        return dict((i+1, key) for i, key in enumerate(self._domain))

    def ngens(self):
        r"""
        Number of generators of this Galois group

        EXAMPLES::

            sage: QuadraticField(-23, 'a').galois_group().ngens()
            1
        """
        return len(self._gens)
