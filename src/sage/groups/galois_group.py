r"""
Galois groups of field extensions.

We don't necessarily require extensions to be normal, but we do require them to be separable.
When an extension is not normal, the Galois group refers to
the automorphism group of the normal closure.

AUTHORS:

- David Roe (2019): initial version
"""

from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic, PermutationGroup_subgroup
from sage.groups.abelian_gps.abelian_group import AbelianGroup_class, AbelianGroup_subgroup
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.structure.category_object import normalize_names
from sage.rings.integer_ring import ZZ

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

class _GMixin:
    r"""
    This class provides some methods for Galois groups to be used for both permutation groups
    and abelian groups, subgroups and full Galois groups.

    It is just intended to provide common functionality between various different Galois group classes.
    """
    @lazy_attribute
    def _default_algorithm(self):
        """
        A string, the default algorithm used for computing the Galois group

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)
            sage: G = K.galois_group()
            sage: G._default_algorithm
            'pari'
        """
        return NotImplemented

    @lazy_attribute
    def _gcdata(self):
        """
        A pair:

        - the Galois closure of the top field in the ambient Galois group;

        - an embedding of the top field into the Galois closure.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 - 2)
            sage: G = K.galois_group()
            sage: G._gcdata
            (Number Field in ac with defining polynomial x^6 + 108,
             Ring morphism:
               From: Number Field in a with defining polynomial x^3 - 2
               To:   Number Field in ac with defining polynomial x^6 + 108
               Defn: a |--> -1/36*ac^4 - 1/2*ac)
        """
        return NotImplemented

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

class _GaloisMixin(_GMixin):
    """
    This class provides methods for Galois groups, allowing concrete instances
    to inherit from both permutation group and abelian group classes.
    """
    @lazy_attribute
    def _field(self):
        """
        The top field, ie the field whose Galois closure elements of this group act upon.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)
            sage: G = K.galois_group()
            sage: G._field
            Number Field in a with defining polynomial x^3 + 2*x + 2
        """
        return NotImplemented

    def _repr_(self):
        """
        String representation of this Galois group

        EXAMPLES::

            sage: from sage.groups.galois_group import GaloisGroup_perm
            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)
            sage: G = K.galois_group()
            sage: GaloisGroup_perm._repr_(G)
            'Galois group of x^3 + 2*x + 2'
        """
        f = self._field.defining_polynomial()
        return "Galois group of %s" % f

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

    @lazy_attribute
    def _field_degree(self):
        """
        Degree of the top field over its base.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)
            sage: L.<b> = K.extension(x^2 + 3*a^2 + 8)
            sage: GK = K.galois_group()
            sage: GL = L.galois_group()
            doctest:warning
            ...
            DeprecationWarning: Use .absolute_field().galois_group() if you want the Galois group of the absolute field
            See https://trac.sagemath.org/28782 for details.
            sage: GK._field_degree
            3

        Despite the fact that `L` is a relative number field, the Galois group
        is computed for the corresponding absolute extension of the rationals.

        This behavior may change in the future::

            sage: GL._field_degree
            6
            sage: GL.transitive_label()
            '6T2'
            sage: GL
            Galois group 6T2 ([3]2) with order 6 of x^2 + 3*a^2 + 8
        """
        try:
            return self._field.degree()
        except NotImplementedError: # relative number fields don't support degree
            return self._field.absolute_degree()

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
        return "%sT%s" % (self._field_degree, self.transitive_number())

    def is_galois(self):
        r"""
        Return whether the top field is Galois over its base.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^8 - x^5 + x^4 - x^3 + 1)
            sage: G = K.galois_group()
            sage: from sage.groups.galois_group import GaloisGroup_perm
            sage: GaloisGroup_perm.is_galois(G)
            False
        """
        return self.order() == self._field_degree

class _SubGaloisMixin(_GMixin):
    """
    This class provides methods for subgroups of Galois groups, allowing concrete instances
    to inherit from both permutation group and abelian group classes.
    """
    @lazy_attribute
    def _ambient_group(self):
        """
        The ambient Galois group of which this is a subgroup.

        EXAMPLES::

            sage: L.<a> = NumberField(x^4 + 1)
            sage: G = L.galois_group()
            sage: H = G.decomposition_group(L.primes_above(3)[0])
            sage: H._ambient_group is G
            True
        """
        return NotImplemented

    @abstract_method(optional=True)
    def fixed_field(self, name=None, polred=None, threshold=None):
        """
        Return the fixed field of this subgroup (as a subfield of the Galois closure).

        INPUT:

        - ``name`` -- a variable name for the new field.

        - ``polred`` -- whether to optimize the generator of the newly created field
            for a simpler polynomial, using pari's polredbest.
            Defaults to ``True`` when the degree of the fixed field is at most 8.

        - ``threshold`` -- positive number; polred only performed if the cost is at most this threshold

        EXAMPLES::

            sage: k.<a> = GF(3^12)
            sage: g = k.galois_group()([8])
            sage: k0, embed = g.fixed_field()
            sage: k0.cardinality()
            81
        """

    @lazy_attribute
    def _gcdata(self):
        """
        The Galois closure data is just that of the ambient group.

        EXAMPLES::

            sage: L.<a> = NumberField(x^4 + 1)
            sage: G = L.galois_group()
            sage: H = G.decomposition_group(L.primes_above(3)[0])
            sage: H.splitting_field() # indirect doctest
            Number Field in a with defining polynomial x^4 + 1
        """
        return self._ambient_group._gcdata

class GaloisGroup_perm(_GaloisMixin, PermutationGroup_generic):
    r"""
    The group of automorphisms of a Galois closure of a given field.

    INPUT:

    - ``field`` -- a field, separable over its base

    - ``names`` -- a string or tuple of length 1, giving a variable name for the splitting field

    - ``gc_numbering`` -- boolean, whether to express permutations in terms of the
        roots of the defining polynomial of the splitting field (versus the defining polynomial
        of the original extension).  The default value may vary based on the type of field.
    """
    @abstract_method
    def transitive_number(self, algorithm=None, recompute=False):
        """
        The transitive number (as in the GAP and Magma databases of transitive groups)
        for the action on the roots of the defining polynomial of the top field.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)
            sage: G = K.galois_group()
            sage: G.transitive_number()
            2
        """

    @lazy_attribute
    def _gens(self):
        """
        The generators of this Galois group as permutations of the roots.  It's important that this
        be computed lazily, since it's often possible to compute other attributes (such as the order
        or transitive number) more cheaply.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5-2)
            sage: G = K.galois_group(gc_numbering=False)
            sage: G._gens
            [(1,2,3,5), (1,4,3,2,5)]
        """
        return NotImplemented

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

    @lazy_attribute
    def _deg(self):
        r"""
        The number of moved points in the permutation representation.

        This will be the degree of the original number field if `_gc_numbering``
        is ``False``, or the degree of the Galois closure otherwise.

        EXAMPLES::

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

class GaloisGroup_ab(_GaloisMixin, AbelianGroup_class):
    r"""
    Abelian Galois groups
    """
    def __init__(self, field, generator_orders, algorithm=None, gen_names='sigma'):
        r"""
        Initialize this Galois group.

        TESTS::

            sage: TestSuite(GF(9).galois_group()).run()
        """
        self._field = field
        self._default_algorithm = algorithm
        AbelianGroup_class.__init__(self, generator_orders, gen_names)

    def is_galois(self):
        r"""
        Abelian extensions are Galois.

        For compatibility with Galois groups of number fields.

        EXAMPLES::

            sage: GF(9).galois_group().is_galois()
            True
        """
        return True

    @lazy_attribute
    def _gcdata(self):
        r"""
        Return the Galois closure (ie, the finite field itself) together with the identity

        EXAMPLES::

            sage: GF(3^2).galois_group()._gcdata
            (Finite Field in z2 of size 3^2,
             Identity endomorphism of Finite Field in z2 of size 3^2)
        """
        k = self._field
        return k, k.Hom(k).identity()

    @cached_method
    def permutation_group(self):
        r"""
        Return a permutation group giving the action on the roots of a defining polynomial.

        This is the regular representation for the abelian group, which is not necessarily the smallest degree permutation representation.

        EXAMPLES::

            sage: GF(3^10).galois_group().permutation_group()
            Permutation Group with generators [(1,2,3,4,5,6,7,8,9,10)]
        """
        return PermutationGroup(gap_group=self._gap_().RegularActionHomomorphism().Image())

    @cached_method(key=_alg_key)
    def transitive_number(self, algorithm=None, recompute=False):
        r"""
        Return the transitive number for the action on the roots of the defining polynomial.

        For abelian groups, there is only one transitive action up to isomorphism
        (left multiplication of the group on itself), so we identify that action.

        EXAMPLES::

            sage: from sage.groups.galois_group import GaloisGroup_ab
            sage: Gtest = GaloisGroup_ab(field=None, generator_orders=(2,2,4))
            sage: Gtest.transitive_number()
            2
        """
        return ZZ(self.permutation_group()._gap_().TransitiveIdentification())

class GaloisGroup_cyc(GaloisGroup_ab):
    r"""
    Cyclic Galois groups
    """
    def transitive_number(self, algorithm=None, recompute=False):
        r"""
        Return the transitive number for the action on the roots of the defining polynomial.

        EXAMPLES::

            sage: GF(2^8).galois_group().transitive_number()
            1
            sage: GF(3^32).galois_group().transitive_number()
            33
            sage: GF(2^60).galois_group().transitive_number()
            Traceback (most recent call last):
            ...
            NotImplementedError: transitive database only computed up to degree 47
        """
        d = self.order()
        if d > 47:
            raise NotImplementedError("transitive database only computed up to degree 47")
        elif d == 32:
            # I don't know why this case is special, but you can check this in Magma (GAP only goes up to 22)
            return ZZ(33)
        else:
            return ZZ(1)

    def signature(self):
        r"""
        Return 1 if contained in the alternating group, -1 otherwise.

        EXAMPLES::

            sage: GF(3^2).galois_group().signature()
            -1
            sage: GF(3^3).galois_group().signature()
            1
        """
        return ZZ(1) if (self._field.degree() % 2) else ZZ(-1)

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

class GaloisSubgroup_ab(AbelianGroup_subgroup, _SubGaloisMixin):
    """
    Subgroups of abelian Galois groups.
    """
    pass

GaloisGroup_perm.Subgroup = GaloisSubgroup_perm
GaloisGroup_ab.Subgroup = GaloisSubgroup_ab
