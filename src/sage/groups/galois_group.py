r"""
Galois groups of field extensions.

We don't necessarily require extensions to be normal, but we do require them to be separable.
When an extension is not normal, the Galois group refers to
the automorphism group of the normal closure.

AUTHORS:

- David Roe (2019): initial version
"""

from sage.groups.abelian_gps.abelian_group import AbelianGroup_class, AbelianGroup_subgroup
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import lazy_import
from sage.rings.integer_ring import ZZ

lazy_import('sage.groups.galois_group_perm', ['GaloisGroup_perm', 'GaloisSubgroup_perm'])
lazy_import('sage.groups.perm_gps.permgroup', 'PermutationGroup')


def _alg_key(self, algorithm=None, recompute=False):
    r"""
    Return a key for use in cached_method calls.

    If recompute is false, will cache using ``None`` as the key, so no
    recomputation will be done.

    If recompute is true, will cache by algorithm, yielding a recomputation
    for each different algorithm.

    EXAMPLES::

        sage: from sage.groups.galois_group import _alg_key
        sage: R.<x> = ZZ[]
        sage: K.<a> = NumberField(x^3 + 2*x + 2)                                        # needs sage.rings.number_field
        sage: G = K.galois_group()                                                      # needs sage.rings.number_field
        sage: _alg_key(G, algorithm='pari', recompute=True)                             # needs sage.rings.number_field
        'pari'
    """
    if recompute:
        algorithm = self._get_algorithm(algorithm)
        return algorithm


class _GMixin:
    r"""
    This class provides some methods for Galois groups to be used for both
    permutation groups and abelian groups, subgroups and full Galois groups.

    It is just intended to provide common functionality between various
    different Galois group classes.
    """
    @lazy_attribute
    def _default_algorithm(self):
        """
        A string, the default algorithm used for computing the Galois group.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)                                    # needs sage.rings.number_field
            sage: G = K.galois_group()                                                  # needs sage.rings.number_field
            sage: G._default_algorithm                                                  # needs sage.rings.number_field
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
            sage: K.<a> = NumberField(x^3 - 2)                                          # needs sage.rings.number_field
            sage: G = K.galois_group()                                                  # needs sage.rings.number_field
            sage: G._gcdata                                                             # needs sage.rings.number_field
            (Number Field in ac with defining polynomial x^6 + 108,
             Ring morphism:
               From: Number Field in a with defining polynomial x^3 - 2
               To:   Number Field in ac with defining polynomial x^6 + 108
               Defn: a |--> -1/36*ac^4 - 1/2*ac)
        """
        return NotImplemented

    def _get_algorithm(self, algorithm):
        r"""
        Allow overriding the default algorithm specified at object creation.

        EXAMPLES::

            sage: # needs sage.rings.number_field
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
            sage: K.<a> = NumberField(x^3 + 2*x + 2)                                    # needs sage.rings.number_field
            sage: G = K.galois_group(names='b')                                         # needs sage.rings.number_field
            sage: G._galois_closure                                                     # needs sage.rings.number_field
            Number Field in b with defining polynomial x^6 + 12*x^4 + 36*x^2 + 140
        """
        return self._gcdata[0]

    def splitting_field(self):
        r"""
        The Galois closure of the top field.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
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
            sage: K.<a> = NumberField(x^3 + 2*x + 2)                                    # needs sage.rings.number_field
            sage: G = K.galois_group(names='b')                                         # needs sage.rings.number_field
            sage: G._gc_map                                                             # needs sage.rings.number_field
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
            sage: K.<a> = NumberField(x^3 + 2*x + 2)                                    # needs sage.rings.number_field
            sage: G = K.galois_group()                                                  # needs sage.rings.number_field
            sage: G._field                                                              # needs sage.rings.number_field
            Number Field in a with defining polynomial x^3 + 2*x + 2
        """
        return NotImplemented

    def _repr_(self):
        """
        String representation of this Galois group.

        EXAMPLES::

            sage: from sage.groups.galois_group import GaloisGroup_perm
            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)                                    # needs sage.rings.number_field
            sage: G = K.galois_group()                                                  # needs sage.rings.number_field
            sage: GaloisGroup_perm._repr_(G)                                            # needs sage.rings.number_field
            'Galois group of x^3 + 2*x + 2'
        """
        f = self._field.defining_polynomial()
        return "Galois group of %s" % f

    def top_field(self):
        r"""
        Return the larger of the two fields in the extension defining this Galois group.

        Note that this field may not be Galois.

        EXAMPLES::

            sage: # needs sage.rings.number_field
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

            sage: # needs sage.rings.number_field
            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^3 + 2*x + 2)
            sage: L.<b> = K.extension(x^2 + 3*a^2 + 8)
            sage: GK = K.galois_group()
            sage: GL = L.galois_group()
            doctest:warning
            ...
            DeprecationWarning: Use .absolute_field().galois_group() if you want the Galois group of the absolute field
            See https://github.com/sagemath/sage/issues/28782 for details.
            sage: GK._field_degree
            3

        Despite the fact that `L` is a relative number field, the Galois group
        is computed for the corresponding absolute extension of the rationals.

        This behavior may change in the future::

            sage: GL._field_degree                                                      # needs sage.rings.number_field
            6
            sage: GL.transitive_label()                                                 # needs sage.rings.number_field
            '6T2'
            sage: GL                                                                    # needs sage.rings.number_field
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
            sage: K.<a> = NumberField(x^8 - x^5 + x^4 - x^3 + 1)                        # needs sage.rings.number_field
            sage: G = K.galois_group()                                                  # needs sage.rings.number_field
            sage: G.transitive_label()                                                  # needs sage.rings.number_field
            '8T44'
        """
        return "%sT%s" % (self._field_degree, self.transitive_number())

    def is_galois(self):
        r"""
        Return whether the top field is Galois over its base.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^8 - x^5 + x^4 - x^3 + 1)                        # needs sage.rings.number_field
            sage: G = K.galois_group()                                                  # needs sage.rings.number_field
            sage: from sage.groups.galois_group import GaloisGroup_perm
            sage: GaloisGroup_perm.is_galois(G)                                         # needs sage.rings.number_field
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

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
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

        - ``name`` -- a variable name for the new field

        - ``polred`` -- whether to optimize the generator of the newly created field
            for a simpler polynomial, using Pari's :pari:`polredbest`;
            defaults to ``True`` when the degree of the fixed field is at most 8

        - ``threshold`` -- positive number; polred only performed if the cost
          is at most this threshold

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
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

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: L.<a> = NumberField(x^4 + 1)
            sage: G = L.galois_group()
            sage: H = G.decomposition_group(L.primes_above(3)[0])
            sage: H.splitting_field()  # indirect doctest
            Number Field in a with defining polynomial x^4 + 1
        """
        return self._ambient_group._gcdata


class GaloisGroup_ab(_GaloisMixin, AbelianGroup_class):
    r"""
    Abelian Galois groups
    """
    def __init__(self, field, generator_orders, algorithm=None, gen_names='sigma'):
        r"""
        Initialize this Galois group.

        TESTS::

            sage: TestSuite(GF(9).galois_group()).run()                                 # needs sage.rings.finite_rings
        """
        self._field = field
        self._default_algorithm = algorithm
        AbelianGroup_class.__init__(self, generator_orders, gen_names)

    def is_galois(self):
        r"""
        Abelian extensions are Galois.

        For compatibility with Galois groups of number fields.

        EXAMPLES::

            sage: GF(9).galois_group().is_galois()                                      # needs sage.rings.finite_rings
            True
        """
        return True

    @lazy_attribute
    def _gcdata(self):
        r"""
        Return the Galois closure (i.e., the finite field itself) together with the identity.

        EXAMPLES::

            sage: GF(3^2).galois_group()._gcdata                                        # needs sage.rings.finite_rings
            (Finite Field in z2 of size 3^2,
             Identity endomorphism of Finite Field in z2 of size 3^2)
        """
        k = self._field
        return k, k.Hom(k).identity()

    @cached_method
    def permutation_group(self):
        r"""
        Return a permutation group giving the action on the roots of a defining polynomial.

        This is the regular representation for the abelian group, which is
        not necessarily the smallest degree permutation representation.

        EXAMPLES::

            sage: GF(3^10).galois_group().permutation_group()                           # needs sage.libs.gap sage.rings.finite_rings
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
            sage: Gtest.transitive_number()                                             # needs sage.libs.gap
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

            sage: GF(2^8).galois_group().transitive_number()                            # needs sage.rings.finite_rings
            1
            sage: GF(3^32).galois_group().transitive_number()                           # needs sage.rings.finite_rings
            33
            sage: GF(2^60).galois_group().transitive_number()                           # needs sage.rings.finite_rings
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

            sage: GF(3^2).galois_group().signature()                                    # needs sage.rings.finite_rings
            -1
            sage: GF(3^3).galois_group().signature()                                    # needs sage.rings.finite_rings
            1
        """
        return ZZ(1) if (self._field.degree() % 2) else ZZ(-1)


class GaloisSubgroup_ab(AbelianGroup_subgroup, _SubGaloisMixin):
    """
    Subgroups of abelian Galois groups.
    """
    pass


GaloisGroup_ab.Subgroup = GaloisSubgroup_ab
