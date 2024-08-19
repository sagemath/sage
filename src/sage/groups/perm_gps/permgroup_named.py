r"""
"Named" Permutation groups (such as the symmetric group, S_n)

You can construct the following permutation groups:

-  :class:`SymmetricGroup`, `S_n` of order `n!` (`n` can also be a list `X` of distinct
   positive integers, in which case it returns `S_X`)

-  :class:`AlternatingGroup`, `A_n` of order `n!/2` (n can also be a list `X`
   of distinct positive integers, in which case it returns `A_X`)

-  :class:`DihedralGroup`, `D_n` of order `2n`

-  :class:`GeneralDihedralGroup`, `Dih(G)`, where `G` is an abelian group

-  :class:`CyclicPermutationGroup`, `C_n` of order `n`

-  :class:`DiCyclicGroup`, nonabelian groups of order `4m` with a unique element of order 2

-  :class:`TransitiveGroup`, `n`-th` transitive group of degree `d`
   from the GAP tables of transitive groups

-  ``TransitiveGroups(d)``, ``TransitiveGroups()``, set of all of the above

-  :class:`PrimitiveGroup`, `n`-th` primitive group of degree `d`
   from the GAP tables of primitive groups

-  ``PrimitiveGroups(d)``, ``PrimitiveGroups()``, set of all of the above

-  ``MathieuGroup(degree)``, Mathieu group of degree 9, 10, 11, 12, 21, 22, 23, or 24.

-  :class:`KleinFourGroup`, subgroup of `S_4` of order `4` which is not `C_2 \times C_2`

-  :class:`QuaternionGroup`, non-abelian group of order `8`, `\{\pm 1, \pm I, \pm J, \pm K\}`

-  :class:`SplitMetacyclicGroup`, nonabelian groups of order `p^m` with cyclic
   subgroups of index p

-  :class:`SemidihedralGroup`, nonabelian 2-groups with cyclic subgroups of index 2

-  ``PGL(n,q)``, projective general linear group of `n\times n` matrices over
   the finite field `\GF(q)`

-  ``PSL(n,q)``, projective special linear group of `n\times n` matrices over
   the finite field `\GF(q)`

-  ``PSp(2n,q)``, projective symplectic linear group of `2n\times 2n` matrices
   over the finite field `\GF(q)`

-  ``PSU(n,q)``, projective special unitary group of `n \times n` matrices having
   coefficients in the finite field `\GF(q^2)` that respect a
   fixed nondegenerate sesquilinear form, of determinant 1.

-  ``PGU(n,q)``, projective general unitary group of `n\times n` matrices having
   coefficients in the finite field `\GF(q^2)` that respect a
   fixed nondegenerate sesquilinear form, modulo the centre.

-  ``SuzukiGroup(q)``, Suzuki group over `\GF(q)`, `^2 B_2(2^{2k+1}) = Sz(2^{2k+1})`.

-  :class:`ComplexReflectionGroup`, the complex reflection group `G(m, p, n)` or
   the exceptional complex reflection group `G_m`

-  :class:`SmallPermutationGroup`, a permutation realization of a group specified by its GAP id.

AUTHOR:

- David Joyner (2007-06): split from permgp.py (suggested by Nick Alexander)

REFERENCES:

- Cameron, P., Permutation Groups. New York: Cambridge University Press, 1999.
- Wielandt, H., Finite Permutation Groups. New York: Academic Press, 1964.
- Dixon, J. and Mortimer, B., Permutation Groups, Springer-Verlag, Berlin/New York, 1996.

.. NOTE::

    Though Suzuki groups are okay, Ree groups should *not* be wrapped as
    permutation groups - the construction is too slow - unless (for
    small values or the parameter) they are made using explicit generators.
"""
# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                          David Joyner <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from pathlib import Path

from sage.arith.misc import factor, valuation
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.groups.abelian_gps.abelian_group import AbelianGroup
from sage.groups.perm_gps.permgroup import PermutationGroup_generic
from sage.groups.perm_gps.permgroup_element import SymmetricGroupElement
from sage.libs.gap.libgap import libgap
from sage.misc.cachefunc import cached_method, weak_cached_function
from sage.misc.functional import is_even
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.rings.integer import Integer
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.primes import Primes
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp
from sage.structure.unique_representation import CachedRepresentation


class PermutationGroup_unique(CachedRepresentation, PermutationGroup_generic):
    """
    .. TODO::

        Fix the broken hash. ::

            sage: G = SymmetricGroup(6)
            sage: G3 = G.subgroup([G((1,2,3,4,5,6)),G((1,2))])
            sage: hash(G) == hash(G3)  # todo: Should be True!
            False

    TESTS::

        sage: G = SymmetricGroup(6)
        sage: G3 = G.subgroup([G((1,2,3,4,5,6)),G((1,2))])
        sage: G == G3
        True
    """
    @weak_cached_function
    def __classcall__(cls, *args, **kwds):
        """
        This makes sure that domain is a FiniteEnumeratedSet before it gets passed
        on to the __init__ method.

        EXAMPLES::

            sage: SymmetricGroup(['a','b']).domain()  # indirect doctest
            {'a', 'b'}
        """
        domain = kwds.pop('domain', None)
        if domain is not None:
            if domain not in FiniteEnumeratedSets():
                domain = FiniteEnumeratedSet(domain)
            kwds['domain'] = domain
        return super().__classcall__(cls, *args, **kwds)


class PermutationGroup_symalt(PermutationGroup_unique):
    """
    This is a class used to factor out some of the commonality
    in the :class:`SymmetricGroup` and :class:`AlternatingGroup` classes.
    """

    @staticmethod
    def __classcall__(cls, domain):
        """
        Normalize the input of the constructor into a set.

        INPUT:

        - ``n`` -- integer or list or tuple thereof

        Calls the constructor with a tuple representing the set.

        EXAMPLES::

            sage: S1 = SymmetricGroup(4)
            sage: S2 = SymmetricGroup([1,2,3,4])
            sage: S3 = SymmetricGroup((1,2,3,4))
            sage: S1 is S2
            True
            sage: S1 is S3
            True

        TESTS::

            sage: SymmetricGroup(0)
            Symmetric group of order 0! as a permutation group
            sage: SymmetricGroup(1)
            Symmetric group of order 1! as a permutation group
            sage: SymmetricGroup(-1)
            Traceback (most recent call last):
            ...
            ValueError: domain (=-1) must be an integer >= 0 or a list
        """
        if domain not in FiniteEnumeratedSets():
            if not isinstance(domain, (tuple, list, range)):
                try:
                    domain = Integer(domain)
                except TypeError:
                    raise TypeError("domain (={}) must be an integer >= 0 or a finite set (but domain has type {})".format(domain, type(domain)))

                if domain < 0:
                    raise ValueError("domain (={}) must be an integer >= 0 or a list".format(domain))
                domain = list(range(1, domain + 1))
            v = FiniteEnumeratedSet(domain)
        else:
            v = domain

        return super().__classcall__(cls, domain=v)


class SymmetricGroup(PermutationGroup_symalt):
    r"""
    The full symmetric group of order `n!`, as a permutation group.

    If `n` is a list or tuple of positive integers then it returns the
    symmetric group of the associated set.

    INPUT:

    - ``n`` -- positive integer, or list or tuple thereof

    .. NOTE::

        This group is also available via ``groups.permutation.Symmetric()``.

    EXAMPLES::

        sage: G = SymmetricGroup(8)
        sage: G.order()
        40320
        sage: G
        Symmetric group of order 8! as a permutation group
        sage: G.degree()
        8
        sage: S8 = SymmetricGroup(8)
        sage: G = SymmetricGroup([1,2,4,5])
        sage: G
        Symmetric group of order 4! as a permutation group
        sage: G.domain()
        {1, 2, 4, 5}
        sage: G = SymmetricGroup(4)
        sage: G
        Symmetric group of order 4! as a permutation group
        sage: G.domain()
        {1, 2, 3, 4}
        sage: G.category()
        Join of Category of finite enumerated permutation groups and
        Category of finite Weyl groups and
        Category of well generated finite irreducible complex reflection groups

    TESTS::

        sage: groups.permutation.Symmetric(4)
        Symmetric group of order 4! as a permutation group
        sage: groups.permutation.Symmetric(1).gens()
        ()

    Check for :issue:`36204`::

        sage: h = SymmetricGroup(1).hom(SymmetricGroup(2))
    """
    def __init__(self, domain=None):
        """
        Initialize ``self``.

        TESTS::

            sage: TestSuite(SymmetricGroup(0)).run()                                    # needs sage.rings.number_field
            sage: TestSuite(SymmetricGroup(1)).run()                                    # needs sage.rings.number_field
            sage: TestSuite(SymmetricGroup(3)).run()                                    # needs sage.rings.number_field
        """
        from sage.categories.finite_weyl_groups import FiniteWeylGroups
        from sage.categories.finite_permutation_groups import FinitePermutationGroups
        from sage.categories.category import Category

        # Note that we skip the call to the superclass initializer in order to
        # avoid infinite recursion since SymmetricGroup is called by
        # PermutationGroupElement
        cat = Category.join([FinitePermutationGroups(),
                             FiniteWeylGroups().Irreducible()])
        super(PermutationGroup_generic, self).__init__(category=cat)

        self._domain = domain
        self._deg = len(self._domain)
        self._domain_to_gap = {key: i+1 for i, key in enumerate(self._domain)}
        self._domain_from_gap = {i+1: key for i, key in enumerate(self._domain)}

        # Create the generators for the symmetric group
        if self._deg <= 1:
            self._gens = ()
        else:
            gens = [tuple(self._domain)]
            if self._deg > 2:
                gens.append(tuple(self._domain[:2]))
            self._gens = tuple([self.element_class(g, self, check=False)
                                for g in gens])

    def _gap_init_(self, gap=None):
        """
        Return the string used to create this group in GAP.

        EXAMPLES::

            sage: S = SymmetricGroup(3)
            sage: S._gap_init_()
            'SymmetricGroup(3)'
            sage: S = SymmetricGroup(['a', 'b', 'c'])
            sage: S._gap_init_()
            'SymmetricGroup(3)'
        """
        return 'SymmetricGroup({})'.format(self.degree())

    @cached_method
    def index_set(self):
        """
        Return the index set for the descents of the symmetric group ``self``.

        EXAMPLES::

            sage: S8 = SymmetricGroup(8)
            sage: S8.index_set()
            (1, 2, 3, 4, 5, 6, 7)

            sage: S = SymmetricGroup([3,1,4,5])
            sage: S.index_set()
            (3, 1, 4)
        """
        return tuple(self.domain()[:-1])

    def __richcmp__(self, x, op):
        """
        Fast comparison for SymmetricGroups.

        EXAMPLES::

            sage: S8 = SymmetricGroup(8)
            sage: S3 = SymmetricGroup(3)
            sage: S8 > S3
            True
        """
        if isinstance(x, SymmetricGroup):
            return richcmp((self._deg, self._domain), (x._deg, x._domain), op)
        return super().__richcmp__(x, op)

    def _repr_(self):
        """
        EXAMPLES::

            sage: A = SymmetricGroup([2,3,7]); A
            Symmetric group of order 3! as a permutation group
        """
        return "Symmetric group of order {}! as a permutation group".format(self.degree())

    def _coerce_map_from_(self, G):
        """
        Return if there is a coercion map from ``G`` into ``self``.

        EXAMPLES::

            sage: J3 = groups.misc.Cactus(3)                                            # needs sage.rings.number_field
            sage: S5 = SymmetricGroup(5)
            sage: S5.coerce_map_from(J3)                                                # needs sage.rings.number_field
            Conversion via _from_cactus_group_element map:
              From: Cactus Group with 3 fruit
              To:   Symmetric group of order 5! as a permutation group
            sage: S2 = SymmetricGroup(2)
            sage: S2._coerce_map_from_(J3) is None                                      # needs sage.rings.number_field
            True
        """
        from sage.groups.cactus_group import CactusGroup
        if isinstance(G, CactusGroup) and G._n <= self._deg:
            return self._from_cactus_group_element
        return super()._coerce_map_from_(G)

    def _from_cactus_group_element(self, x):
        """
        Return an element of ``self`` from a cactus group element.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: J3 = groups.misc.Cactus(3)
            sage: s12,s13,s23 = J3.gens()
            sage: elt = s12*s23*s13
            sage: S5 = SymmetricGroup(5)
            sage: S5._from_cactus_group_element(elt)
            (2,3)
        """
        return self(x.to_permutation())

    def cartan_type(self):
        r"""
        Return the Cartan type of ``self``.

        The symmetric group `S_n` is a Coxeter group of type `A_{n-1}`.

        EXAMPLES::

            sage: A = SymmetricGroup([2,3,7]); A.cartan_type()
            ['A', 2]

            sage: A = SymmetricGroup([]); A.cartan_type()
            ['A', 0]
        """
        from sage.combinat.root_system.cartan_type import CartanType
        return CartanType(['A', max(self.degree() - 1, 0)])

    def coxeter_matrix(self):
        r"""
        Return the Coxeter matrix of ``self``.

        EXAMPLES::

            sage: A = SymmetricGroup([2,3,7,'a']); A.coxeter_matrix()                   # needs sage.graphs
            [1 3 2]
            [3 1 3]
            [2 3 1]
        """
        return self.cartan_type().coxeter_matrix()

    def simple_reflection(self, i):
        r"""
        For `i` in the index set of ``self``, return the
        elementary transposition `s_i = (i,i+1)`.

        EXAMPLES::

            sage: A = SymmetricGroup(5)
            sage: A.simple_reflection(3)
            (3,4)

            sage: A = SymmetricGroup([2,3,7])
            sage: A.simple_reflections()
            Finite family {2: (2,3), 3: (3,7)}
        """
        return self([(i, self._domain[self._domain.index(i)+1])], check=False)

    def reflections(self):
        """
        Return the list of all reflections in ``self``.

        EXAMPLES::

            sage: A = SymmetricGroup(3)
            sage: A.reflections()
            [(1,2), (1,3), (2,3)]
        """
        from itertools import combinations
        dom = self._domain
        return [self([(i, j)], check=False) for i, j in combinations(dom, 2)]

    def young_subgroup(self, comp):
        """
        Return the Young subgroup associated with the composition ``comp``.

        EXAMPLES::

            sage: S = SymmetricGroup(8)
            sage: c = Composition([2,2,2,2])
            sage: S.young_subgroup(c)
            Subgroup generated by [(1,2), (3,4), (5,6), (7,8)] of
             (Symmetric group of order 8! as a permutation group)

            sage: S = SymmetricGroup(['a','b','c'])
            sage: S.young_subgroup([2,1])
            Subgroup generated by [('a','b')] of
             (Symmetric group of order 3! as a permutation group)

            sage: Y = S.young_subgroup([2,2,2,2,2])
            Traceback (most recent call last):
            ...
            ValueError: the composition is not of expected size
        """
        if sum(comp) != self.degree():
            raise ValueError('the composition is not of expected size')

        domain = self._domain
        gens = []
        pos = 0
        for c in comp:
            gens.extend(self((domain[pos + i], domain[pos + i + 1]))
                        for i in range(c - 1))
            pos += c

        return self.subgroup(gens)

    def major_index(self, parameter=None):
        r"""
        Return the *major index generating polynomial* of ``self``.

        This is a gadget counting the elements of ``self`` by major index.

        INPUT:

        - ``parameter`` -- an element of a ring; the result is
          more explicit with a formal variable (default:
          element ``q`` of Univariate Polynomial Ring in ``q`` over
          Integer Ring)

        .. MATH::

            P(q) = \sum_{g\in S_n} q^{ \operatorname{major\ index}(g) }

        EXAMPLES::

            sage: S4 = SymmetricGroup(4)
            sage: S4.major_index()                                                      # needs sage.combinat
            q^6 + 3*q^5 + 5*q^4 + 6*q^3 + 5*q^2 + 3*q + 1
            sage: K.<t> = QQ[]
            sage: S4.major_index(t)                                                     # needs sage.combinat
            t^6 + 3*t^5 + 5*t^4 + 6*t^3 + 5*t^2 + 3*t + 1
        """
        from sage.combinat.q_analogues import q_factorial
        return q_factorial(self.degree(), parameter)

    def conjugacy_classes_representatives(self):
        r"""
        Return a complete list of representatives of conjugacy classes in
        a permutation group `G`.

        Let `S_n` be the symmetric group on `n` letters. The conjugacy
        classes are indexed by partitions `\lambda` of `n`. The ordering
        of the conjugacy classes is reverse lexicographic order of
        the partitions.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: G.conjugacy_classes_representatives()                                 # needs sage.combinat
            [(), (1,2), (1,2)(3,4), (1,2,3), (1,2,3)(4,5),
             (1,2,3,4), (1,2,3,4,5)]

        ::

            sage: S = SymmetricGroup(['a','b','c'])
            sage: S.conjugacy_classes_representatives()                                 # needs sage.combinat
            [(), ('a','b'), ('a','b','c')]

        TESTS:

        Check some border cases::

            sage: S = SymmetricGroup(0)
            sage: S.conjugacy_classes_representatives()                                 # needs sage.combinat
            [()]
            sage: S = SymmetricGroup(1)
            sage: S.conjugacy_classes_representatives()                                 # needs sage.combinat
            [()]
        """
        from sage.combinat.partition import Partitions_n
        from sage.groups.perm_gps.symgp_conjugacy_class import default_representative
        n = len(self.domain())
        return [default_representative(la, self)
                for la in reversed(Partitions_n(n))]

    def conjugacy_classes_iterator(self):
        """
        Iterate over the conjugacy classes of ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: list(G.conjugacy_classes_iterator()) == G.conjugacy_classes()         # needs sage.combinat
            True
        """
        from sage.combinat.partition import Partitions_n
        from sage.groups.perm_gps.symgp_conjugacy_class import SymmetricGroupConjugacyClass
        P = Partitions_n(len(self.domain()))
        for la in reversed(P):
            yield SymmetricGroupConjugacyClass(self, la)

    def conjugacy_classes(self):
        """
        Return a list of the conjugacy classes of ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: G.conjugacy_classes()                                                 # needs sage.combinat
            [Conjugacy class of cycle type [1, 1, 1, 1, 1] in
                 Symmetric group of order 5! as a permutation group,
             Conjugacy class of cycle type [2, 1, 1, 1] in
                 Symmetric group of order 5! as a permutation group,
             Conjugacy class of cycle type [2, 2, 1] in
                 Symmetric group of order 5! as a permutation group,
             Conjugacy class of cycle type [3, 1, 1] in
                 Symmetric group of order 5! as a permutation group,
             Conjugacy class of cycle type [3, 2] in
                 Symmetric group of order 5! as a permutation group,
             Conjugacy class of cycle type [4, 1] in
                 Symmetric group of order 5! as a permutation group,
             Conjugacy class of cycle type [5] in
                 Symmetric group of order 5! as a permutation group]
        """
        return list(self.conjugacy_classes_iterator())

    def conjugacy_class(self, g):
        r"""
        Return the conjugacy class of ``g`` inside the symmetric
        group ``self``.

        INPUT:

        - ``g`` -- a partition or an element of the symmetric group ``self``

        OUTPUT: a conjugacy class of a symmetric group

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: g = G((1,2,3,4))
            sage: G.conjugacy_class(g)                                                  # needs sage.combinat
            Conjugacy class of cycle type [4, 1] in
             Symmetric group of order 5! as a permutation group
        """
        from sage.groups.perm_gps.symgp_conjugacy_class import SymmetricGroupConjugacyClass
        return SymmetricGroupConjugacyClass(self, g)

    def algebra(self, base_ring, category=None):
        r"""
        Return the symmetric group algebra associated to ``self``.

        INPUT:

        - ``base_ring`` -- a ring
        - ``category`` -- a category (default: the category of ``self``)

        If ``self`` is the symmetric group on `1,\ldots,n`, then this
        is special cased to take advantage of the features in
        :class:`SymmetricGroupAlgebra`. Otherwise the usual group
        algebra is returned.

        EXAMPLES::

            sage: S4 = SymmetricGroup(4)
            sage: S4.algebra(QQ)                                                        # needs sage.combinat
            Symmetric group algebra of order 4 over Rational Field

            sage: S3 = SymmetricGroup([1,2,3])
            sage: A = S3.algebra(QQ); A                                                 # needs sage.combinat
            Symmetric group algebra of order 3 over Rational Field
            sage: a = S3.an_element(); a
            (2,3)
            sage: A(a)                                                                  # needs sage.combinat
            (2,3)

        We illustrate the choice of the category::

            sage: A.category()                                                          # needs sage.combinat
            Join of Category of Coxeter group algebras over Rational Field
                and Category of finite group algebras over Rational Field
                and Category of finite dimensional cellular algebras with basis
                     over Rational Field
            sage: A = S3.algebra(QQ, category=Semigroups())                             # needs sage.combinat
            sage: A.category()                                                          # needs sage.combinat
            Category of finite dimensional unital cellular semigroup algebras
             over Rational Field

        In the following case, a usual group algebra is returned::

            sage: S = SymmetricGroup([2,3,5])
            sage: S.algebra(QQ)                                                         # needs sage.combinat
            Algebra of Symmetric group of order 3! as a permutation group over Rational Field
            sage: a = S.an_element(); a
            (3,5)
            sage: S.algebra(QQ)(a)                                                      # needs sage.combinat
            (3,5)
        """
        from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
        domain = self.domain()
        if list(domain) == list(range(1, len(domain) + 1)):
            return SymmetricGroupAlgebra(base_ring, self, category=category)
        else:
            return super().algebra(base_ring)

    Element = SymmetricGroupElement


class AlternatingGroup(PermutationGroup_symalt):
    def __init__(self, domain=None):
        """
        The alternating group of order `n!/2`, as a permutation group.

        INPUT:

        - ``n`` -- positive integer, or list or tuple thereof

        .. NOTE::

            This group is also available via ``groups.permutation.Alternating()``.

        EXAMPLES::

            sage: G = AlternatingGroup(6)
            sage: G.order()
            360
            sage: G
            Alternating group of order 6!/2 as a permutation group
            sage: G.category()
            Category of finite enumerated permutation groups
            sage: TestSuite(G).run()  # long time

            sage: G = AlternatingGroup([1,2,4,5])
            sage: G
            Alternating group of order 4!/2 as a permutation group
            sage: G.domain()
            {1, 2, 4, 5}
            sage: G.category()
            Category of finite enumerated permutation groups
            sage: TestSuite(G).run()

        TESTS::

            sage: groups.permutation.Alternating(6)
            Alternating group of order 6!/2 as a permutation group

        Check that alternating groups have their generators in the correct
        order (:issue:`37590`)::

            sage: A = AlternatingGroup(6)
            sage: A.gens()
            ((1,2,3,4,5), (4,5,6))
            sage: A.gap().GeneratorsOfGroup()
            [ (1,2,3,4,5), (4,5,6) ]
            sage: h = A.hom(codomain=A, im_gens=A.gens())
        """
        PermutationGroup_symalt.__init__(self, gap_group='AlternatingGroup(%s)' % len(domain), domain=domain)

    def _repr_(self):
        """
        EXAMPLES::

            sage: A = AlternatingGroup([2,3,7]); A
            Alternating group of order 3!/2 as a permutation group
        """
        return "Alternating group of order %s!/2 as a permutation group" % self.degree()

    def _gap_init_(self, gap=None):
        """
        Return the string used to create this group in GAP.

        EXAMPLES::

            sage: A = AlternatingGroup(3)
            sage: A._gap_init_()
            'AlternatingGroup(3)'
            sage: A = AlternatingGroup(['a', 'b', 'c'])
            sage: A._gap_init_()
            'AlternatingGroup(3)'
        """
        return 'AlternatingGroup(%s)' % self.degree()


class CyclicPermutationGroup(PermutationGroup_unique):
    def __init__(self, n):
        """
        A cyclic group of order `n`, as a permutation group.

        INPUT:

        - ``n`` -- positive integer

        .. NOTE::

            This group is also available via ``groups.permutation.Cyclic()``.

        EXAMPLES::

            sage: G = CyclicPermutationGroup(8)
            sage: G.order()
            8
            sage: G
            Cyclic group of order 8 as a permutation group
            sage: G.category()
            Category of finite enumerated permutation groups
            sage: TestSuite(G).run()
            sage: C = CyclicPermutationGroup(10)
            sage: C.is_abelian()
            True
            sage: C = CyclicPermutationGroup(10)
            sage: C.as_AbelianGroup()
            Multiplicative Abelian group isomorphic to C2 x C5

        TESTS::

            sage: groups.permutation.Cyclic(6)
            Cyclic group of order 6 as a permutation group
        """
        n = Integer(n)
        if n < 1:
            raise ValueError("n (=%s) must be >= 1" % n)
        gens = tuple(range(1, n+1))
        PermutationGroup_generic.__init__(self, [gens], n)

    def _repr_(self):
        """
        EXAMPLES::

            sage: CyclicPermutationGroup(8)
            Cyclic group of order 8 as a permutation group
        """
        return "Cyclic group of order %s as a permutation group" % self.order()

    def is_commutative(self):
        """
        Return ``True`` if this group is commutative.

        EXAMPLES::

            sage: C = CyclicPermutationGroup(8)
            sage: C.is_commutative()
            True
        """
        return True

    def is_abelian(self):
        """
        Return ``True`` if this group is abelian.

        EXAMPLES::

            sage: C = CyclicPermutationGroup(8)
            sage: C.is_abelian()
            True
        """
        return True

    def as_AbelianGroup(self):
        """
        Return the corresponding :class:`AbelianGroup` instance.

        EXAMPLES::

            sage: C = CyclicPermutationGroup(8)
            sage: C.as_AbelianGroup()
            Multiplicative Abelian group isomorphic to C8
        """
        n = self.order()
        a = list(factor(n))
        invs = [x[0]**x[1] for x in a]
        G = AbelianGroup(len(a), invs)
        return G


class DiCyclicGroup(PermutationGroup_unique):
    r"""
    The dicyclic group of order `4n`, for `n\geq 2`.

    INPUT:

    - ``n`` -- positive integer, two or greater

    OUTPUT:

    This is a nonabelian group similar in some respects to the
    dihedral group of the same order, but with far fewer
    elements of order 2 (it has just one).  The permutation
    representation constructed here is based on the presentation

    .. MATH::

        \langle a, x\mid a^{2n}=1, x^{2}=a^{n}, x^{-1}ax=a^{-1}\rangle

    For `n=2` this is the group of quaternions
    (`{\pm 1, \pm I,\pm J, \pm K}`), which is the nonabelian
    group of order 8 that is not the dihedral group `D_4`,
    the symmetries of a square.  For `n=3` this is the nonabelian
    group of order 12 that is not the dihedral group `D_6`
    nor the alternating group `A_4`.  This group of order 12 is
    also the semi-direct product of `C_2` by `C_4`,
    `C_3\rtimes C_4`.  [Con]_


    When the order of the group is a
    power of 2 it is known as a "generalized quaternion group."

    IMPLEMENTATION:

    The presentation above means every element can be written as
    `a^{i}x^{j}` with `0\leq i<2n`, `j=0,1`.  We code `a^i` as the symbol
    `i+1` and code `a^{i}x` as the symbol `2n+i+1`.  The two generators
    are then represented using a left regular representation.

    .. NOTE::

        This group is also available via ``groups.permutation.DiCyclic()``.

    EXAMPLES:

    A dicyclic group of order 384, with a large power of 2 as a divisor::

        sage: n = 3*2^5
        sage: G = DiCyclicGroup(n)
        sage: G.order()
        384
        sage: a = G.gen(0)
        sage: x = G.gen(1)
        sage: a^(2*n)
        ()
        sage: a^n==x^2
        True
        sage: x^-1*a*x==a^-1
        True

    A large generalized quaternion group (order is a power of 2)::

        sage: n = 2^10
        sage: G = DiCyclicGroup(n)
        sage: G.order()
        4096
        sage: a = G.gen(0)
        sage: x = G.gen(1)
        sage: a^(2*n)
        ()
        sage: a^n==x^2
        True
        sage: x^-1*a*x==a^-1
        True

    Just like the dihedral group, the dicyclic group has
    an element whose order is half the order of the group.
    Unlike the dihedral group, the dicyclic group has only
    one element of order 2.  Like the dihedral groups of
    even order, the center of the dicyclic group is a
    subgroup of order 2 (thus has the unique element of
    order 2 as its non-identity element). ::

        sage: G = DiCyclicGroup(3*5*4)
        sage: G.order()
        240
        sage: two = [g for g in G if g.order()==2]; two
        [(1,5)(2,6)(3,7)(4,8)(9,13)(10,14)(11,15)(12,16)]
        sage: G.center().order()
        2

    For small orders, we check this is really a group
    we do not have in Sage otherwise. ::

        sage: G = DiCyclicGroup(2)
        sage: H = DihedralGroup(4)
        sage: G.is_isomorphic(H)
        False
        sage: G = DiCyclicGroup(3)
        sage: H = DihedralGroup(6)
        sage: K = AlternatingGroup(6)
        sage: G.is_isomorphic(H) or G.is_isomorphic(K)
        False

    TESTS::

        sage: groups.permutation.DiCyclic(6)
        Dicyclic group of order 24 as a permutation group

    AUTHOR:

    - Rob Beezer (2009-10-18)
    """
    def __init__(self, n):
        r"""
        The dicyclic group of order `4n`, as a permutation group.

        INPUT:

        - ``n`` -- a positive integer; 2 or greater

        EXAMPLES::

            sage: G = DiCyclicGroup(3*8)
            sage: G.order()
            96
            sage: TestSuite(G).run()
        """
        n = Integer(n)
        if n < 2:
            raise ValueError("n (=%s) must be 2 or greater" % n)

        # Certainly 2^2 is part of the first factor of the order
        #   r is maximum power of 2 in the order
        #   m is the rest, the odd part
        order = 4*n
        factored = order.factor()
        r = factored[0][0]**factored[0][1]
        m = order//r
        halfr, fourthr = r//2, r//4

        # Representation of  a
        # Two cycles of length halfr
        a = [tuple(range(1, halfr+1)), tuple(range(halfr+1, r+1))]
        # With an odd part, a cycle of length m will give the right order for a
        if m > 1:
            a.append(tuple(range(r + 1, r + m + 1)))

        # Representation of  x
        # Four-cycles that will conjugate the generator  a  properly
        x = [(i+1, (-i) % halfr + halfr + 1, (fourthr+i) % halfr + 1, (-fourthr-i) % halfr + halfr + 1)
             for i in range(0, fourthr)]
        # With an odd part, transpositions will conjugate the m-cycle to create inverse
        if m > 1:
            x += [(r+i+1, r+m-i) for i in range(0, (m-1)//2)]

        PermutationGroup_generic.__init__(self, gens=[a, x])

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: DiCyclicGroup(12)
            Dicyclic group of order 48 as a permutation group
        """
        return "Dicyclic group of order %s as a permutation group" % self.order()

    def is_commutative(self):
        r"""
        Return ``True`` if this group is commutative.

        EXAMPLES::

            sage: D = DiCyclicGroup(12)
            sage: D.is_commutative()
            False
        """
        return False

    def is_abelian(self):
        r"""
        Return ``True`` if this group is abelian.

        EXAMPLES::

            sage: D = DiCyclicGroup(12)
            sage: D.is_abelian()
            False
        """
        return False


class KleinFourGroup(PermutationGroup_unique):
    def __init__(self):
        r"""
        The Klein 4 Group, which has order `4` and exponent `2`, viewed
        as a subgroup of `S_4`.

        OUTPUT: the Klein 4 group of order 4, as a permutation group of degree 4

        .. NOTE::

            This group is also available via ``groups.permutation.KleinFour()``.

        EXAMPLES::

            sage: G = KleinFourGroup(); G
            The Klein 4 group of order 4, as a permutation group
            sage: sorted(G)
            [(), (3,4), (1,2), (1,2)(3,4)]

        TESTS::

            sage: G.category()
            Category of finite enumerated permutation groups
            sage: TestSuite(G).run()

            sage: groups.permutation.KleinFour()
            The Klein 4 group of order 4, as a permutation group

        AUTHOR:
            -- Bobby Moretti (2006-10)
        """
        gens = [(1, 2), (3, 4)]
        PermutationGroup_generic.__init__(self, gens)

    def _repr_(self):
        """
        EXAMPLES::

            sage: G = KleinFourGroup(); G
            The Klein 4 group of order 4, as a permutation group
        """
        return 'The Klein 4 group of order 4, as a permutation group'


class JankoGroup(PermutationGroup_unique):
    def __init__(self, n):
        r"""
        Janko Groups `J1, J2`, and `J3`.
        (Note that `J4` is too big to be treated here.)

        INPUT:

        - ``n`` -- integer among `\{1,2,3\}`

        EXAMPLES::

            sage: G = groups.permutation.Janko(1); G                            # optional - gap_package_atlasrep internet
            Janko group J1 of order 175560 as a permutation group

        TESTS::

            sage: G.category()                                                  # optional - gap_package_atlasrep internet
            Category of finite enumerated permutation groups
            sage: TestSuite(G).run(skip=["_test_enumerated_set_contains",       # optional - gap_package_atlasrep internet
            ....:                        "_test_enumerated_set_iter_list"])
        """
        if n not in [1, 2, 3]:
            raise ValueError("n must belong to {1,2,3}")
        self._n = n
        libgap.load_package("atlasrep")
        id = 'AtlasGroup("J%s")' % n
        PermutationGroup_generic.__init__(self, gap_group=id)

    def _repr_(self):
        """
        EXAMPLES::

            sage: G = groups.permutation.Janko(1); G                            # optional - gap_package_atlasrep internet
            Janko group J1 of order 175560 as a permutation group
        """
        return "Janko group J%s of order %s as a permutation group" % (self._n, self.order())


class SuzukiSporadicGroup(PermutationGroup_unique):
    def __init__(self):
        r"""
        Suzuki Sporadic Group.

        EXAMPLES::

            sage: G = groups.permutation.SuzukiSporadic(); G                    # optional - gap_package_atlasrep internet
            Sporadic Suzuki group acting on 1782 points

        TESTS::

            sage: G.category() # optional - gap_package_atlasrep internet
            Category of finite enumerated permutation groups
            sage: TestSuite(G).run(skip=["_test_enumerated_set_contains",       # optional - gap_package_atlasrep internet
            ....:                        "_test_enumerated_set_iter_list"])
        """
        libgap.load_package("atlasrep")
        PermutationGroup_generic.__init__(self, gap_group='AtlasGroup("Suz")')

    def _repr_(self):
        """
        EXAMPLES::

            sage: G = groups.permutation.SuzukiSporadic(); G                    # optional - gap_package_atlasrep internet
            Sporadic Suzuki group acting on 1782 points
        """
        return "Sporadic Suzuki group acting on 1782 points"


class QuaternionGroup(DiCyclicGroup):
    r"""
    The quaternion group of order 8.

    OUTPUT:

    The quaternion group of order 8, as a permutation group.
    See the :class:`DiCyclicGroup` class for a generalization of this
    construction.

    .. NOTE::

        This group is also available via ``groups.permutation.Quaternion()``.

    EXAMPLES:

    The quaternion group is one of two non-abelian groups of order 8,
    the other being the dihedral group `D_4`.  One way to describe this
    group is with three generators, `I, J, K`, so the whole group is
    then given as the set `\{\pm 1, \pm I, \pm J, \pm K\}` with relations
    such as `I^2=J^2=K^2=-1`, `IJ=K` and `JI=-K`.

    The examples below illustrate how to use this group in a similar
    manner, by testing some of these relations.  The representation used
    here is the left-regular representation. ::

        sage: Q = QuaternionGroup()
        sage: I = Q.gen(0)
        sage: J = Q.gen(1)
        sage: K = I*J
        sage: [I,J,K]
        [(1,2,3,4)(5,6,7,8), (1,5,3,7)(2,8,4,6), (1,8,3,6)(2,7,4,5)]
        sage: neg_one = I^2; neg_one
        (1,3)(2,4)(5,7)(6,8)
        sage: J^2 == neg_one and K^2 == neg_one
        True
        sage: J*I == neg_one*K
        True
        sage: Q.center().order() == 2
        True
        sage: neg_one in Q.center()
        True

    TESTS::

        sage: groups.permutation.Quaternion()
        Quaternion group of order 8 as a permutation group

    AUTHOR:

    - Rob Beezer (2009-10-09)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: Q = QuaternionGroup()
            sage: TestSuite(Q).run()
        """
        DiCyclicGroup.__init__(self, 2)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: Q = QuaternionGroup(); Q
            Quaternion group of order 8 as a permutation group
        """
        return "Quaternion group of order 8 as a permutation group"


class GeneralDihedralGroup(PermutationGroup_generic):
    r"""
    The Generalized Dihedral Group generated by the abelian group with
    direct factors in the input list.

    INPUT:

    - ``factors`` -- list of the sizes of the cyclic factors of the
      abelian group being dihedralized (this will be sorted once
      entered)

    OUTPUT:

    For a given abelian group (noting that each finite abelian group
    can be represented as the direct product of cyclic groups), the
    General Dihedral Group it generates is simply the semi-direct
    product of the given group with `C_2`, where the nonidentity
    element of `C_2` acts on the abelian group by turning each element
    into its inverse. In this implementation, each input abelian group
    will be standardized so as to act on a minimal amount of letters.
    This will be done by breaking the direct factors into products of
    p-groups, before this new set of factors is ordered from smallest
    to largest for complete standardization. Note that the generalized
    dihedral group corresponding to a cyclic group, `C_n`, is simply
    the dihedral group `D_n`.

    EXAMPLES:

    As is noted in [TW1980]_, `Dih(C_3 \times C_3)` has the presentation

    .. MATH::

        \langle a, b, c\mid a^{3}, b^{3}, c^{2}, ab = ba, ac = ca^{-1}, bc = cb^{-1} \rangle

    Note also the fact, verified by [TW1980]_, that the dihedralization of
    `C_3 \times C_3` is the only nonabelian group of order 18
    with no element of order 6. ::

        sage: G = GeneralDihedralGroup([3,3])
        sage: G
        Generalized dihedral group generated by C3 x C3
        sage: G.order()
        18
        sage: G.gens()
        ((4,5,6), (2,3)(5,6), (1,2,3))
        sage: a = G.gens()[2]; b = G.gens()[0]; c = G.gens()[1]
        sage: a.order() == 3, b.order() == 3, c.order() == 2
        (True, True, True)
        sage: a*b == b*a, a*c == c*a.inverse(), b*c == c*b.inverse()
        (True, True, True)
        sage: G.subgroup([a,b,c]) == G
        True
        sage: G.is_abelian()
        False
        sage: all(x.order() != 6 for x in G)
        True

    If all of the direct factors are `C_2`, then the action turning
    each element into its inverse is trivial, and the semi-direct
    product becomes a direct product. ::

        sage: G = GeneralDihedralGroup([2,2,2])
        sage: G.order()
        16
        sage: G.gens()
        ((7,8), (5,6), (3,4), (1,2))
        sage: G.is_abelian()
        True
        sage: H = KleinFourGroup()
        sage: G.is_isomorphic(H.direct_product(H)[0])
        True

    If two nonidentical input lists generate isomorphic abelian
    groups, then they will generate identical groups (with each direct
    factor broken up into its prime factors), but they will still have
    distinct descriptions. Note that if `\gcd(n,m)=1`, then `C_n \times
    C_m \cong C_{nm}`, while the general dihedral groups
    generated by isomorphic abelian groups should be themselves
    isomorphic. ::

        sage: G = GeneralDihedralGroup([6,34,46,14])
        sage: H = GeneralDihedralGroup([7,17,3,46,2,2,2])
        sage: G == H, G.gens() == H.gens()
        (True, True)
        sage: [x.order() for x in G.gens()]
        [23, 17, 7, 2, 3, 2, 2, 2, 2]
        sage: G
        Generalized dihedral group generated by C6 x C34 x C46 x C14
        sage: H
        Generalized dihedral group generated by C7 x C17 x C3 x C46 x C2 x C2 x C2

    A cyclic input yields a Classical Dihedral Group. ::

        sage: G = GeneralDihedralGroup([6])
        sage: D = DihedralGroup(6)
        sage: G.is_isomorphic(D)
        True

    A Generalized Dihedral Group will always have size twice the
    underlying group, be solvable (as it has an abelian subgroup with
    index 2), and, unless the underlying group is of the form
    `{C_2}^n`, be nonabelian (by the structure theorem of finite
    abelian groups and the fact that a semi-direct product is a
    direct product only when the underlying action is trivial). ::

        sage: G = GeneralDihedralGroup([6,18,33,60])
        sage: (6*18*33*60)*2
        427680
        sage: G.order()
        427680
        sage: G.is_solvable()
        True
        sage: G.is_abelian()
        False

    TESTS::

        sage: G = GeneralDihedralGroup("foobar")
        Traceback (most recent call last):
        ...
        TypeError: factors of abelian group must be a list, not foobar

        sage: GeneralDihedralGroup([])
        Traceback (most recent call last):
        ...
        ValueError: there must be at least one direct factor in the abelian group being dihedralized

        sage: GeneralDihedralGroup([3, 1.5])
        Traceback (most recent call last):
        ...
        TypeError: the input list must consist of Integers

        sage: GeneralDihedralGroup([4, -8])
        Traceback (most recent call last):
        ...
        ValueError: all direct factors must be greater than 1

    AUTHOR:

    - Kevin Halasz (2012-7-12)
    """
    def __init__(self, factors):
        r"""
        Init method of class <GeneralDihedralGroup>. See the docstring
        for :class:`<GeneralDihedralGroup>`.

        EXAMPLES::

            sage: G = GeneralDihedralGroup([5,5,5])
            sage: G.order()
            250
            sage: TestSuite(G).run() # long time
        """
        if not isinstance(factors, list):
            msg = "factors of abelian group must be a list, not {}"
            raise TypeError(msg.format(factors))

        if len(factors) < 1:
            raise ValueError('there must be at least one direct factor in the abelian group being dihedralized')

        if not all(isinstance(x, Integer) for x in factors):
            raise TypeError('the input list must consist of Integers')

        if not all(x >= 2 for x in factors):
            s = 'all direct factors must be greater than 1'
            raise ValueError(s)

        self.factors = factors
        # To get uniform outputs for isomorphic inputs, we break
        # each inputted cyclic group into a direct product of cyclic
        # p-groups
        simplified = sorted([term[0]**term[1] for a in factors for term in a.factor()])

        gens = []
        # genx is an element of order two that turns each of the
        # generators of the abelian group into its inverse via
        # conjugation
        genx = []
        jumppoint = Integer(1)
        for a in simplified:
            # create one of the generators for the abelian group
            gens.append([tuple(range(jumppoint, jumppoint + a))])
            # make contribution to the generator that dihedralizes the
            # abelian group
            genx.extend((jumppoint + i, jumppoint + a - i)
                        for i in range(1, (a//2) + 1) if i != a - i)
            jumppoint += a
        # If all of the direct factors are C2, then the action turning
        # each element into its inverse is trivial, and the
        # semi-direct product becomes a direct product, so we simply
        # tack on another disjoint transposition
        if all(x == 2 for x in simplified):
            genx.append((jumppoint, jumppoint+1))
        gens.append(genx)
        PermutationGroup_generic.__init__(self, gens=gens)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: G = GeneralDihedralGroup([2,4,8])
            sage: G
            Generalized dihedral group generated by C2 x C4 x C8
        """
        grouplist = [f'C{n}' for n in self.factors]
        return 'Generalized dihedral group generated by ' + ' x '.join(grouplist)


class DihedralGroup(PermutationGroup_unique):
    def __init__(self, n):
        r"""
        The Dihedral group of order `2n` for any integer `n\geq 1`.

        INPUT:

        - ``n`` -- positive integer

        OUTPUT: the dihedral group of order `2n`, as a permutation group

        .. NOTE::

          This group is also available via ``groups.permutation.Dihedral()``.

        EXAMPLES::

            sage: DihedralGroup(1)
            Dihedral group of order 2 as a permutation group

            sage: DihedralGroup(2)
            Dihedral group of order 4 as a permutation group
            sage: DihedralGroup(2).gens()
            ((3,4), (1,2))

            sage: DihedralGroup(5).gens()
            ((1,2,3,4,5), (1,5)(2,4))
            sage: sorted(DihedralGroup(5))
            [(), (2,5)(3,4), (1,2)(3,5), (1,2,3,4,5), (1,3)(4,5), (1,3,5,2,4),
             (1,4)(2,3), (1,4,2,5,3), (1,5,4,3,2), (1,5)(2,4)]

            sage: G = DihedralGroup(6)
            sage: G.order()
            12
            sage: G = DihedralGroup(5)
            sage: G.order()
            10
            sage: G
            Dihedral group of order 10 as a permutation group
            sage: G.gens()
            ((1,2,3,4,5), (1,5)(2,4))

            sage: DihedralGroup(0)
            Traceback (most recent call last):
            ...
            ValueError: n must be positive

        TESTS::

            sage: TestSuite(G).run()
            sage: G.category()
            Category of finite enumerated permutation groups
            sage: TestSuite(G).run()

            sage: groups.permutation.Dihedral(6)
            Dihedral group of order 12 as a permutation group
        """
        n = Integer(n)
        if n <= 0:
            raise ValueError("n must be positive")

        # the first generator generates the cyclic subgroup of D_n, <(1...n)> in
        # cycle notation
        gen0 = range(1, n + 1)

        if n < 1:
            raise ValueError("n (=%s) must be >= 1" % n)

        # D_1 is a subgroup of S_2, we need the cyclic group of order 2
        if n == 1:
            gens = CyclicPermutationGroup(2).gens()
        elif n == 2:
            gens = ((1, 2), (3, 4))
        else:
            gen1 = tuple((i, n - i + 1) for i in range(1, n // 2 + 1))
            gens = (tuple(gen0), gen1)

        PermutationGroup_generic.__init__(self, gens)

    def _repr_(self):
        """
        EXAMPLES::

            sage: G = DihedralGroup(5); G
            Dihedral group of order 10 as a permutation group
        """
        return f"Dihedral group of order {self.order()} as a permutation group"


class SplitMetacyclicGroup(PermutationGroup_unique):
    def __init__(self, p, m):
        r"""
        The split metacyclic group of order `p^m`.

        INPUT:

        - ``p`` -- a prime number that is the prime underlying this
          p-group

        - ``m`` -- positive integer such that the order of this
          group is the `p^m`. Be aware that, for even `p`, `m` must be
          greater than 3, while for odd `p`, `m` must be greater than
          2.

        OUTPUT:

        The split metacyclic group of order `p^m`. This family of
        groups has presentation

        .. MATH::

            \langle x, y\mid x^{p^{m-1}}, y^{p}, y^{-1}xy=x^{1+p^{m-2}} \rangle

        This family is notable because, for odd `p`, these are the
        only `p`-groups with a cyclic subgroup of index `p`, a
        result proven in [Gor1980]_. It is also shown in
        [Gor1980]_ that this is one of four families containing
        nonabelian 2-groups with a cyclic subgroup of index 2
        (with the others being the dicyclic groups, the dihedral
        groups, and the semidihedral groups).

        EXAMPLES:

        Using the last relation in the group's presentation,
        one can see that the elements of the form `y^{i}x`,
        `0 \leq i \leq p-1` all have order `p^{m-1}`, as it can be
        shown that their `p` th powers are all `x^{p^{m-2}+p}`,
        an element with order `p^{m-2}`. Manipulation of the same
        relation shows that none of these elements are powers of
        any other. Thus, there are `p` cyclic maximal subgroups in
        each split metacyclic group. It is also proven in
        [Gor1980]_ that this family has commutator subgroup
        of order `p`, and the Frattini subgroup is equal to the
        center, with this group being cyclic of order `p^{m-2}`.
        These characteristics are necessary to identify these
        groups in the case that `p=2`, although the possession of
        a cyclic maximal subgroup in a non-abelian `p`-group is
        enough for odd `p` given the group's order. ::

            sage: G = SplitMetacyclicGroup(2,8)
            sage: G.order() == 2**8
            True
            sage: G.is_abelian()
            False
            sage: len([H for H in G.subgroups() if H.order() == 2^7 and H.is_cyclic()])
            2
            sage: G.commutator().order()
            2
            sage: G.frattini_subgroup() == G.center()
            True
            sage: G.center().order() == 2^6
            True
            sage: G.center().is_cyclic()
            True

            sage: G = SplitMetacyclicGroup(3,3)
            sage: len([H for H in G.subgroups() if H.order() == 3^2 and H.is_cyclic()])
            3
            sage: G.commutator().order()
            3
            sage: G.frattini_subgroup() == G.center()
            True
            sage: G.center().order()
            3

        TESTS::

            sage: G = SplitMetacyclicGroup(3,2.5)
            Traceback (most recent call last):
            ...
            TypeError: both p and m must be integers

            sage: G = SplitMetacyclicGroup(4,3)
            Traceback (most recent call last):
            ...
            ValueError: p must be prime, 4 is not prime

            sage: G = SplitMetacyclicGroup(2,2)
            Traceback (most recent call last):
            ...
            ValueError: if prime is 2, the exponent must be greater than 3, not 2

            sage: G = SplitMetacyclicGroup(11,2)
            Traceback (most recent call last):
            ...
            ValueError: if prime is odd, the exponent must be greater than 2, not 2

        AUTHOR:

        - Kevin Halasz (2012-8-7)
        """
        if not isinstance(p, Integer) or not isinstance(m, Integer):
            raise TypeError('both p and m must be integers')

        if p not in Primes():
            raise ValueError('p must be prime, %s is not prime' % p)

        if p == 2 and m <= 3:
            raise ValueError('if prime is 2, the exponent must be greater than 3, not %s' % m)

        if p % 2 == 1 and m <= 2:
            raise ValueError('if prime is odd, the exponent must be greater than 2, not %s' % m)

        self.p = p
        self.m = m

        # x is created with list, rather than cycle, notation
        x = list(range(2, p ** (m - 1) + 1))
        x.append(1)
        # y is also created with list notation, where the list
        # used is equivalent to the cycle notation representation of
        # x^(1+p^(m-2)). This technique is inspired by exercise 5.30
        # Judson's "Abstract Algebra" (abstract.pugetsound.edu).
        y = [1]
        point = 1
        for i in range(p**(m-1)-1):
            next = (point + 1 + p**(m-2)) % (p**(m-1))
            if next == 0:
                next = p**(m-1)
            y.append(next)
            point = next
        PermutationGroup_unique.__init__(self, gens=[x, y])

    def _repr_(self):
        """
        EXAMPLES::

            sage: G = SplitMetacyclicGroup(7,4);G
            The split metacyclic group of order 7 ^ 4
        """
        return 'The split metacyclic group of order %s ^ %s' % (self.p, self.m)


class SemidihedralGroup(PermutationGroup_unique):
    def __init__(self, m):
        r"""
        The semidihedral group of order `2^m`.

        INPUT:

        - ``m`` -- positive integer; the power of 2 that is the
          group's order

        OUTPUT:

        The semidihedral group of order `2^m`. These groups can be
        thought of as a semidirect product of `C_{2^{m-1}}` with
        `C_2`, where the nontrivial element of `C_2` is sent to
        the element of the automorphism group of `C_{2^{m-1}}` that
        sends elements to their `-1+2^{m-2}` th power. Thus, the
        group has the presentation:

        .. MATH::

            \langle x, y\mid x^{2^{m-1}}, y^{2}, y^{-1}xy=x^{-1+2^{m-2}} \rangle

        This family is notable because it is made up of non-abelian
        2-groups that all contain cyclic subgroups of index 2. It
        is one of only four such families.

        EXAMPLES:

        In [Gor1980]_ it is shown that the semidihedral groups
        have center of order 2. It is also shown that they have a
        Frattini subgroup equal to their commutator, which is a
        cyclic subgroup of order `2^{m-2}`. ::

            sage: G = SemidihedralGroup(12)
            sage: G.order() == 2^12
            True
            sage: G.commutator() == G.frattini_subgroup()
            True
            sage: G.commutator().order() == 2^10
            True
            sage: G.commutator().is_cyclic()
            True
            sage: G.center().order()
            2

            sage: G = SemidihedralGroup(4)
            sage: len([H for H in G.subgroups() if H.is_cyclic() and H.order() == 8])
            1
            sage: G.gens()
            ((2,4)(3,7)(6,8), (1,2,3,4,5,6,7,8))
            sage: x = G.gens()[1]; y = G.gens()[0]
            sage: x.order() == 2^3; y.order() == 2
            True
            True
            sage: y*x*y == x^(-1+2^2)
            True

        TESTS::

            sage: G = SemidihedralGroup(4.4)
            Traceback (most recent call last):
            ...
            TypeError: m must be an integer, not 4.40000000000000

            sage: G = SemidihedralGroup(-5)
            Traceback (most recent call last):
            ...
            ValueError: the exponent must be greater than 3, not -5

        AUTHOR:

        - Kevin Halasz (2012-8-7)
        """
        if not isinstance(m, Integer):
            raise TypeError('m must be an integer, not %s' % m)

        if m <= 3:
            raise ValueError('the exponent must be greater than 3, not %s' % m)

        self.m = m

        # x is created with list, rather than cycle, notation
        x = list(range(2, 2 ** (m - 1) + 1))
        x.append(1)
        # y is also created with list notation, where the list
        # used is equivalent to the cycle notation representation of
        # x^(1+p^(m-2)). This technique is inspired by exercise 5.30
        # Judson's "Abstract Algebra" (abstract.pugetsound.edu).
        y = [1]
        k = 1
        for i in range(2**(m-1)-1):
            next = (k - 1 + 2**(m-2)) % (2**(m-1))
            if next == 0:
                next = 2**(m-1)
            y.append(next)
            k = next
        PermutationGroup_unique.__init__(self, gens=[x, y])

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: G = SemidihedralGroup(6); G
            The semidihedral group of order 64
        """
        return 'The semidihedral group of order %s' % (2**self.m)


class MathieuGroup(PermutationGroup_unique):
    def __init__(self, n):
        """
        The Mathieu group of degree `n`.

        INPUT:

        - ``n`` -- positive integer in {9, 10, 11, 12, 21, 22, 23, 24}

        OUTPUT: the Mathieu group of degree `n`, as a permutation group

        .. NOTE::

            This group is also available via ``groups.permutation.Mathieu()``.

        EXAMPLES::

            sage: G = MathieuGroup(12); G
            Mathieu group of degree 12 and order 95040 as a permutation group

        TESTS::

            sage: G.category()
            Category of finite enumerated permutation groups
            sage: TestSuite(G).run(skip=["_test_enumerated_set_contains", "_test_enumerated_set_iter_list"])

            sage: groups.permutation.Mathieu(9)
            Mathieu group of degree 9 and order 72 as a permutation group

        Note: this is a fairly big group, so we skip some tests that
        currently require to list all the elements of the group,
        because there is no proper iterator yet.
        """
        n = Integer(n)
        self._n = n
        if n not in [9, 10, 11, 12, 21, 22, 23, 24]:
            raise ValueError("argument must belong to {9, 10, 11, 12, 21, 22, 23, 24}")
        id = 'MathieuGroup(%s)' % n
        PermutationGroup_generic.__init__(self, gap_group=id)

    def _repr_(self):
        """
        EXAMPLES::

            sage: G = MathieuGroup(12); G
            Mathieu group of degree 12 and order 95040 as a permutation group
        """
        return "Mathieu group of degree %s and order %s as a permutation group" % (self._n, self.order())


class TransitiveGroup(PermutationGroup_unique):
    def __init__(self, d, n):
        """
        The transitive group from the GAP tables of transitive groups.

        INPUT:

        - ``d`` -- nonnegative integer; the degree
        - ``n`` -- positive integer; the index of the group in the GAP database,
          starting at 1

        OUTPUT: the `n`-th transitive group of degree `d`

        .. NOTE::

            This group is also available via ``groups.permutation.Transitive()``.

        EXAMPLES::

            sage: TransitiveGroup(0,1)
            Transitive group number 1 of degree 0
            sage: TransitiveGroup(1,1)
            Transitive group number 1 of degree 1
            sage: G = TransitiveGroup(5, 2); G
            Transitive group number 2 of degree 5
            sage: G.gens()
            ((1,2,3,4,5), (1,4)(2,3))

            sage: G.category()
            Category of finite enumerated permutation groups

        .. WARNING:: this follows GAP's naming convention of indexing
          the transitive groups starting from ``1``::

            sage: TransitiveGroup(5,0)
            Traceback (most recent call last):
            ...
            ValueError: index n must be in {1,..,5}

        .. WARNING:: only transitive groups of "small" degree are
          available in GAP's database::

            sage: TransitiveGroup(32,1)
            Traceback (most recent call last):
            ...
            NotImplementedError: only the transitive groups of degree at most 31
            are available in GAP's database

        TESTS::


            sage: groups.permutation.Transitive(1, 1)
            Transitive group number 1 of degree 1

            sage: TestSuite(TransitiveGroup(0,1)).run()
            sage: TestSuite(TransitiveGroup(1,1)).run()
            sage: TestSuite(TransitiveGroup(5,2)).run()

            sage: TransitiveGroup(1,5)
            Traceback (most recent call last):
            ...
            ValueError: index n must be in {1,..,1}
        """
        self._d = d = Integer(d)
        self._n = n = Integer(n)
        if d < 0:
            raise ValueError("degree d must not be negative")
        max_n = TransitiveGroups(d).cardinality()
        if n > max_n or n <= 0:
            raise ValueError("index n must be in {1,..,%s}" % max_n)
        if d <= 1:
            PermutationGroup_generic.__init__(self, gens=[()], domain=list(range(1, d+1)))
        else:
            gap_group = libgap.TransitiveGroup(d, n)
            PermutationGroup_generic.__init__(self, gap_group=gap_group)

    def transitive_number(self):
        """
        Return the index of this group in the GAP database, starting at 1.

        EXAMPLES::

            sage: TransitiveGroup(8, 44).transitive_number()
            44
        """
        return self._n

    def degree(self):
        """
        Return the degree of this permutation group.

        EXAMPLES::

            sage: TransitiveGroup(8, 44).degree()
            8
        """
        return self._d

    def _repr_(self):
        """
        EXAMPLES::

            sage: G = TransitiveGroup(1,1); G
            Transitive group number 1 of degree 1
        """
        return "Transitive group number %s of degree %s" % (self._n, self._d)


def TransitiveGroups(d=None):
    """
    INPUT:

    - ``d`` -- integer (optional)

    Return the set of all transitive groups of a given degree
    ``d`` up to isomorphisms. If ``d`` is not specified, it returns the set of all
    transitive groups up to isomorphisms.

    EXAMPLES::

        sage: TransitiveGroups(3)
        Transitive Groups of degree 3
        sage: TransitiveGroups(7)
        Transitive Groups of degree 7
        sage: TransitiveGroups(8)
        Transitive Groups of degree 8

        sage: TransitiveGroups()
        Transitive Groups

    .. WARNING:: in practice, the database currently only contains
      transitive groups up to degree 31::

        sage: TransitiveGroups(32).cardinality()
        Traceback (most recent call last):
        ...
        NotImplementedError: only the transitive groups of degree at most 31
        are available in GAP's database
    """
    if d is None:
        return TransitiveGroupsAll()

    d = Integer(d)
    if d < 0:
        raise ValueError("a transitive group acts on a nonnegative integer number of positions")
    return TransitiveGroupsOfDegree(d)


class TransitiveGroupsAll(DisjointUnionEnumeratedSets):
    """
    The infinite set of all transitive groups up to isomorphisms.

    EXAMPLES::

        sage: L = TransitiveGroups(); L
        Transitive Groups
        sage: L.category()
        Category of facade infinite enumerated sets
        sage: L.cardinality()
        +Infinity

        sage: p = L.__iter__()
        sage: (next(p), next(p), next(p), next(p), next(p), next(p), next(p), next(p))
        (Transitive group number 1 of degree 0, Transitive group number 1 of degree 1,
         Transitive group number 1 of degree 2, Transitive group number 1 of degree 3,
         Transitive group number 2 of degree 3, Transitive group number 1 of degree 4,
         Transitive group number 2 of degree 4, Transitive group number 3 of degree 4)
    """
    def __init__(self):
        """
        TESTS::

            sage: S = TransitiveGroups()
            sage: S.category()
            Category of facade infinite enumerated sets
            sage: TestSuite(TransitiveGroups()).run()
        """
        DisjointUnionEnumeratedSets.__init__(self,
                                             Family(NonNegativeIntegers(),
                                                    TransitiveGroups))

    # We override the __call__ as the elements are not instances of Element
    __call__ = DisjointUnionEnumeratedSets._element_constructor_facade

    def _repr_(self):
        """
        TESTS::

            sage: TransitiveGroups()  # indirect doctest
            Transitive Groups
        """
        return "Transitive Groups"

    def __contains__(self, G):
        r"""
        EXAMPLES::

            sage: TransitiveGroup(5,2) in TransitiveGroups()
            True
            sage: TransitiveGroup(6,5) in TransitiveGroups()
            True
            sage: 1 in TransitiveGroups()
            False
        """
        return isinstance(G, TransitiveGroup)


class TransitiveGroupsOfDegree(CachedRepresentation, Parent):
    r"""
    The set of all transitive groups of a given (small) degree up
    to isomorphism.

    EXAMPLES::

        sage: S = TransitiveGroups(4); S
        Transitive Groups of degree 4
        sage: list(S)
        [Transitive group number 1 of degree 4,
         Transitive group number 2 of degree 4,
         Transitive group number 3 of degree 4,
         Transitive group number 4 of degree 4,
         Transitive group number 5 of degree 4]

        sage: TransitiveGroups(5).an_element()
        Transitive group number 1 of degree 5

    We write the cardinality of all transitive groups of degree 5::

        sage: for G in TransitiveGroups(5):
        ....:     print(G.cardinality())
        5
        10
        20
        60
        120

    TESTS::

        sage: TestSuite(TransitiveGroups(3)).run()
    """
    def __init__(self, n):
        """
        TESTS::

            sage: S = TransitiveGroups(4)
            sage: S.category()
            Category of finite enumerated sets
        """
        self._degree = n
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: TransitiveGroups(6)
            Transitive Groups of degree 6
        """
        return "Transitive Groups of degree %s" % (self._degree,)

    def __contains__(self, G):
        r"""
        EXAMPLES::

            sage: TransitiveGroup(6,5) in TransitiveGroups(4)
            False
            sage: TransitiveGroup(4,3) in TransitiveGroups(4)
            True
            sage: 1 in TransitiveGroups(4)
            False
        """
        return isinstance(G, TransitiveGroup) and G._d == self._degree

    def __getitem__(self, n):
        r"""
        Return the `n`-th transitive group of a given degree.

        INPUT:

        - ``n`` -- positive integer

        EXAMPLES::

            sage: TransitiveGroups(5)[3]
            Transitive group number 3 of degree 5

        .. WARNING::

            This follows GAP's naming convention of indexing
            the transitive groups starting from ``1``::

                sage: TransitiveGroups(5)[0]
                Traceback (most recent call last):
                ...
                ValueError: index n must be in {1,..,5}
        """
        return TransitiveGroup(self._degree, n)

    def __iter__(self):
        """
        EXAMPLES::

            sage: list(TransitiveGroups(5))  # indirect doctest
            [Transitive group number 1 of degree 5,
             Transitive group number 2 of degree 5,
             Transitive group number 3 of degree 5,
             Transitive group number 4 of degree 5,
             Transitive group number 5 of degree 5]
        """
        for n in range(1, self.cardinality() + 1):
            yield self[n]

    def __call__(self, G):
        r"""
        Convert ``G`` to an element of ``self`` if possible.

        We override ``Parent.__call__`` as this is not a proper parent.

        EXAMPLES::

            sage: T2 = TransitiveGroups(2)
            sage: all(T2(H) is H for H in T2)
            True

            sage: T2(SymmetricGroup(2))
            Transitive group number 1 of degree 2

            sage: T2(SymmetricGroup(3))
            Traceback (most recent call last):
            ...
            ValueError: Symmetric group of order 3! as a permutation group is not a transitive group of degree 2
        """
        if G in self:
            return G
        for H in self:
            if H.is_isomorphic(G):
                return H
        raise ValueError("{} is not a transitive group of degree {}".format(G, self._degree))

    @cached_method
    def cardinality(self):
        r"""
        Return the cardinality of ``self``, that is the number of
        transitive groups of a given degree.

        EXAMPLES::

            sage: TransitiveGroups(0).cardinality()
            1
            sage: TransitiveGroups(2).cardinality()
            1
            sage: TransitiveGroups(7).cardinality()
            7
            sage: TransitiveGroups(12).cardinality()
            301
            sage: [TransitiveGroups(i).cardinality() for i in range(11)]
            [1, 1, 1, 2, 5, 5, 16, 7, 50, 34, 45]

        .. WARNING::

            GAP comes with a database containing all transitive groups
            up to degree 31::

                sage: TransitiveGroups(32).cardinality()
                Traceback (most recent call last):
                ...
                NotImplementedError: only the transitive groups of degree at most 31
                are available in GAP's database

        TESTS::

            sage: type(TransitiveGroups(12).cardinality())
            <class 'sage.rings.integer.Integer'>
            sage: type(TransitiveGroups(0).cardinality())
            <class 'sage.rings.integer.Integer'>
        """
        # gap.NrTransitiveGroups(0) fails, so Sage needs to handle this

        # While we are at it, and since Sage also handles the
        # transitive group of degree 1, we may as well handle 1
        if self._degree <= 1:
            return Integer(1)
        else:
            try:
                return Integer(libgap.NrTransitiveGroups(libgap(self._degree)))
            except RuntimeError:
                from sage.misc.verbose import verbose
                verbose("Error: TransitiveGroups should come with GAP.", level=0)
            except TypeError:
                raise NotImplementedError("only the transitive groups of degree at most 31 are available in GAP's database")


class PrimitiveGroup(PermutationGroup_unique):
    """
    The primitive group from the GAP tables of primitive groups.

    INPUT:

    - ``d`` -- nonnegative integer; the degree of the group

    - ``n`` -- positive integer; the index of the group in the GAP
      database, starting at 1

    OUTPUT: the ``n``-th primitive group of degree ``d``

    EXAMPLES::

        sage: PrimitiveGroup(0,1)
        Trivial group
        sage: PrimitiveGroup(1,1)
        Trivial group
        sage: G = PrimitiveGroup(5, 2); G
        D(2*5)
        sage: G.gens()
        ((2,4)(3,5), (1,2,3,5,4))
        sage: G.category()
        Category of finite enumerated permutation groups

    .. WARNING::

        this follows GAP's naming convention of indexing the primitive
        groups starting from ``1``::

            sage: PrimitiveGroup(5,0)
            Traceback (most recent call last):
            ...
            ValueError: index n must be in {1,..,5}

    Only primitive groups of "small" degree are available in GAP's
    database::

        sage: PrimitiveGroup(2^12,1)
        Traceback (most recent call last):
        ...
        GAPError: Error, Primitive groups of degree 4096 are not known!
    """

    def __init__(self, d, n):
        """
        The Python constructor.

        INPUT/OUTPUT: see :class:`PrimitiveGroup`

        TESTS::

            sage: TestSuite(PrimitiveGroup(0,1)).run()
            sage: TestSuite(PrimitiveGroup(1,1)).run()
            sage: TestSuite(PrimitiveGroup(5,2)).run()
            sage: PrimitiveGroup(6,5)
            Traceback (most recent call last):
            ...
            ValueError: index n must be in {1,..,4}
        """
        d = Integer(d)
        n = Integer(n)
        if d < 0:
            raise ValueError("degree d must not be negative")
        max_n = PrimitiveGroups(d).cardinality()
        if n > max_n or n <= 0:
            raise ValueError("index n must be in {1,..,%s}" % max_n)

        if d <= 1:
            PermutationGroup_generic.__init__(self, gens=[()], domain=list(range(1, d+1)))
            self._pretty_name = "Trivial group"
        else:
            gap_group = libgap.PrimitiveGroup(d, n)
            self._pretty_name = str(gap_group)
            PermutationGroup_generic.__init__(self, gap_group=gap_group)

        self._d = d
        self._n = n

    def _repr_(self):
        """
        Return a string representation.

        EXAMPLES::

            sage: G = PrimitiveGroup(5,1); G
            C(5)
        """
        return self._pretty_name

    def group_primitive_id(self):
        """
        Return the index of this group in the GAP database of primitive groups.

        OUTPUT: a positive integer, following GAP's conventions

        EXAMPLES::

            sage: G = PrimitiveGroup(5,2); G.group_primitive_id()
            2
        """
        return self._n


def PrimitiveGroups(d=None):
    """
    Return the set of all primitive groups of a given degree ``d``.

    INPUT:

    - ``d`` -- integer (optional)

    OUTPUT:

    The set of all primitive groups of a given degree ``d`` up to
    isomorphisms using GAP. If ``d`` is not specified, it returns the
    set of all primitive groups up to isomorphisms stored in GAP.

    EXAMPLES::

        sage: PrimitiveGroups(3)
        Primitive Groups of degree 3
        sage: PrimitiveGroups(7)
        Primitive Groups of degree 7
        sage: PrimitiveGroups(8)
        Primitive Groups of degree 8
        sage: PrimitiveGroups()
        Primitive Groups

    The database is currently limited::

         sage: PrimitiveGroups(2^12).cardinality()
         Traceback (most recent call last):
         ...
         GAPError: Error, Primitive groups of degree 4096 are not known!

    .. TODO::

        This enumeration helper could be extended based on
        ``PrimitiveGroupsIterator`` in GAP.  This method allows to
        enumerate groups with specified properties such as transitivity,
        solvability, ..., without creating all groups.
    """
    if d is None:
        return PrimitiveGroupsAll()
    else:
        d = Integer(d)
        if d < 0:
            raise ValueError("a primitive group acts on a nonnegative integer number of positions")
        return PrimitiveGroupsOfDegree(d)


class PrimitiveGroupsAll(DisjointUnionEnumeratedSets):
    """
    The infinite set of all primitive groups up to isomorphisms.

    EXAMPLES::

        sage: L = PrimitiveGroups(); L
        Primitive Groups
        sage: L.category()
        Category of facade infinite enumerated sets
        sage: L.cardinality()
        +Infinity

        sage: p = L.__iter__()
        sage: (next(p), next(p), next(p), next(p),
        ....:  next(p), next(p), next(p), next(p))
        (Trivial group, Trivial group, S(2), A(3), S(3), A(4), S(4), C(5))

    TESTS:

    The following test is broken, see :issue:`22576`::

        sage: TestSuite(PrimitiveGroups()).run()  # known bug, long time
    """
    def __init__(self):
        """
        TESTS::

            sage: S = PrimitiveGroups()
            sage: S.category()
            Category of facade infinite enumerated sets
        """
        DisjointUnionEnumeratedSets.__init__(self,
                                             Family(NonNegativeIntegers(),
                                                    PrimitiveGroups))

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT: string

        TESTS::

            sage: PrimitiveGroups()  # indirect doctest
            Primitive Groups
        """
        return "Primitive Groups"

    def __contains__(self, G):
        r"""
        Test whether `G` is in ``self``.

        INPUT:

        - ``G`` -- anything

        OUTPUT: boolean

        EXAMPLES::

            sage: PrimitiveGroup(5,2) in PrimitiveGroups()
            True
            sage: PrimitiveGroup(6,4) in PrimitiveGroups()
            True
            sage: 1 in PrimitiveGroups()
            False
        """
        return isinstance(G, PrimitiveGroup)


class PrimitiveGroupsOfDegree(CachedRepresentation, Parent):
    """
    The set of all primitive groups of a given degree up to isomorphisms.

    EXAMPLES::

        sage: S = PrimitiveGroups(5); S
        Primitive Groups of degree 5
        sage: S.list()
        [C(5), D(2*5), AGL(1, 5), A(5), S(5)]
        sage: S.an_element()
        C(5)

    We write the cardinality of all primitive groups of degree 5::

        sage: for G in PrimitiveGroups(5):
        ....:     print(G.cardinality())
        5
        10
        20
        60
        120

    TESTS::

        sage: TestSuite(PrimitiveGroups(3)).run()
    """
    def __init__(self, n):
        """
        TESTS::

            sage: S = PrimitiveGroups(4)
            sage: S.category()
            Category of finite enumerated sets
        """
        self._degree = n
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT: string

        TESTS::

            sage: PrimitiveGroups(6)
            Primitive Groups of degree 6
        """
        return "Primitive Groups of degree %s" % (self._degree)

    def __contains__(self, G):
        r"""
        Test whether `G` is in ``self``.

        INPUT:

        - ``G`` -- anything

        OUTPUT: boolean

        EXAMPLES::

            sage: PrimitiveGroup(6,4) in PrimitiveGroups(4)
            False
            sage: PrimitiveGroup(4,2) in PrimitiveGroups(4)
            True
            sage: 1 in PrimitiveGroups(4)
            False
        """
        return isinstance(G, PrimitiveGroup) and G._d == self._degree

    def __getitem__(self, n):
        r"""
        Return the `n`-th primitive group of a given degree.

        INPUT:

        - ``n`` -- positive integer

        EXAMPLES::

            sage: PrimitiveGroups(5)[3]
            AGL(1, 5)

        .. WARNING::

            this follows GAP's naming convention of indexing the
            primitive groups starting from ``1``::

                sage: PrimitiveGroups(5)[0]
                Traceback (most recent call last):
                ...
                ValueError: index n must be in {1,..,5}
        """
        return PrimitiveGroup(self._degree, n)

    def __iter__(self):
        """
        EXAMPLES::

            sage: list(PrimitiveGroups(5))  # indirect doctest
            [C(5), D(2*5), AGL(1, 5), A(5), S(5)]
        """
        for n in range(1, self.cardinality() + 1):
            yield self[n]

    @cached_method
    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        OUTPUT: integer; the number of primitive groups of a given degree
        up to isomorphism

        EXAMPLES::

            sage: PrimitiveGroups(0).cardinality()
            1
            sage: PrimitiveGroups(2).cardinality()
            1
            sage: PrimitiveGroups(7).cardinality()
            7
            sage: PrimitiveGroups(12).cardinality()
            6
            sage: [PrimitiveGroups(i).cardinality() for i in range(11)]
            [1, 1, 1, 2, 2, 5, 4, 7, 7, 11, 9]

        TESTS::

            sage: type(PrimitiveGroups(12).cardinality())
            <class 'sage.rings.integer.Integer'>
            sage: type(PrimitiveGroups(0).cardinality())
            <class 'sage.rings.integer.Integer'>

        Check for :issue:`31774`::

            sage: PrimitiveGroups(2500).cardinality()
            34
            sage: PrimitiveGroups(2^12).cardinality()
            Traceback (most recent call last):
            ...
            GAPError: Error, Primitive groups of degree 4096 are not known!
        """
        if self._degree <= 1:
            # gap.NrPrimitiveGroups(0) fails, so Sage needs to handle this
            # While we are at it, and since Sage also handles the
            # primitive group of degree 1, we may as well handle 1
            return Integer(1)
        else:
            return Integer(libgap.NrPrimitiveGroups(self._degree))


class PermutationGroup_plg(PermutationGroup_unique):
    def base_ring(self):
        """
        EXAMPLES::

            sage: G = PGL(2,3)
            sage: G.base_ring()
            Finite Field of size 3

            sage: G = PSL(2,3)
            sage: G.base_ring()
            Finite Field of size 3
        """
        return self._base_ring

    def matrix_degree(self):
        """
        EXAMPLES::

            sage: G = PSL(2,3)
            sage: G.matrix_degree()
            2
        """
        return self._n


class PGL(PermutationGroup_plg):
    def __init__(self, n, q, name='a'):
        r"""
        The projective general linear groups over `\GF(q)`.

        INPUT:

        - ``n`` -- positive integer; the degree
        - ``q`` -- prime power; the size of the ground field
        - ``name`` -- (default: ``'a'``) variable name of indeterminate of finite field `\GF(q)`

        OUTPUT: PGL(`n`, `q`)

        .. NOTE::

            This group is also available via ``groups.permutation.PGL()``.

        EXAMPLES::

            sage: G = PGL(2,3); G
            Permutation Group with generators [(3,4), (1,2,4)]
            sage: print(G)
            The projective general linear group of degree 2 over Finite Field of size 3
            sage: G.base_ring()
            Finite Field of size 3
            sage: G.order()
            24

            sage: G = PGL(2, 9, 'b'); G                                                 # needs sage.rings.finite_rings
            Permutation Group with generators [(3,10,9,8,4,7,6,5), (1,2,4)(5,6,8)(7,9,10)]
            sage: G.base_ring()
            Finite Field in b of size 3^2

            sage: G.category()
            Category of finite enumerated permutation groups
            sage: TestSuite(G).run()  # long time

        TESTS::

            sage: groups.permutation.PGL(2, 3)
            Permutation Group with generators [(3,4), (1,2,4)]
        """
        id = 'Group([()])' if n == 1 else 'PGL(%s,%s)' % (n, q)
        PermutationGroup_generic.__init__(self, gap_group=id)
        self._q = q
        self._base_ring = GF(q, name=name)
        self._n = n

    def __str__(self):
        """
        EXAMPLES::

            sage: G = PGL(2,3); G
            Permutation Group with generators [(3,4), (1,2,4)]
            sage: print(G)
            The projective general linear group of degree 2 over Finite Field of size 3
        """
        return "The projective general linear group of degree %s over %s" % (self._n, self.base_ring())


class PSL(PermutationGroup_plg):
    def __init__(self, n, q, name='a'):
        r"""
        The projective special linear groups over `\GF(q)`.

        INPUT:

        - ``n`` -- positive integer; the degree
        - ``q`` -- either a prime power (the size of the ground field) or a finite field
        - ``name`` -- (default: ``'a'``) variable name of indeterminate of finite field `\GF(q)`

        OUTPUT: the group PSL(`n`, `q`)

        .. NOTE::

            This group is also available via ``groups.permutation.PSL()``.

        EXAMPLES::

            sage: G = PSL(2,3); G
            Permutation Group with generators [(2,3,4), (1,2)(3,4)]
            sage: G.order()
            12
            sage: G.base_ring()
            Finite Field of size 3
            sage: print(G)
            The projective special linear group of degree 2 over Finite Field of size 3

        We create two groups over nontrivial finite fields::

            sage: G = PSL(2, 4, 'b'); G                                                 # needs sage.rings.finite_rings
            Permutation Group with generators [(3,4,5), (1,2,3)]
            sage: G.base_ring()
            Finite Field in b of size 2^2
            sage: G = PSL(2, 8); G                                                      # needs sage.rings.finite_rings
            Permutation Group with generators [(3,8,6,4,9,7,5), (1,2,3)(4,7,5)(6,9,8)]
            sage: G.base_ring()
            Finite Field in a of size 2^3

            sage: G.category()
            Category of finite enumerated permutation groups
            sage: TestSuite(G).run()  # long time

        TESTS::

            sage: groups.permutation.PSL(2, 3)
            Permutation Group with generators [(2,3,4), (1,2)(3,4)]

        Check that :issue:`7424` is handled::

            sage: PSL(2, GF(7,'x'))
            Permutation Group with generators [(3,7,5)(4,8,6), (1,2,6)(3,4,8)]
            sage: PSL(2, GF(3))
            Permutation Group with generators [(2,3,4), (1,2)(3,4)]
            sage: PSL(2, QQ)
            Traceback (most recent call last):
            ...
            ValueError: q must be a prime power or a finite field
        """
        from sage.categories.finite_fields import FiniteFields
        if q in FiniteFields():
            name = q.gen()
            q = q.cardinality()
        if q not in NonNegativeIntegers():
            raise ValueError('q must be a prime power or a finite field')
        if n == 1:
            id = 'Group([()])'
        else:
            id = 'PSL(%s,%s)' % (n, q)
        PermutationGroup_generic.__init__(self, gap_group=id)
        self._q = q
        self._base_ring = GF(q, name=name)
        self._n = n

    def __str__(self):
        """
        EXAMPLES::

            sage: G = PSL(2,3)
            sage: print(G)
            The projective special linear group of degree 2 over Finite Field of size 3
        """
        return "The projective special linear group of degree %s over %s" % (self._n, self.base_ring())

    def ramification_module_decomposition_hurwitz_curve(self):
        r"""
        Helps compute the decomposition of the ramification module
        for the Hurwitz curves X (over `\CC` say) with automorphism group
        G = PSL(2,q), q a "Hurwitz prime" (ie, p is `\pm 1 \pmod 7`).
        Using this computation and Borne's formula helps determine the
        G-module structure of the RR spaces of equivariant
        divisors can be determined explicitly.

        The output is a list of integer multiplicities: [m1,...,mn],
        where n is the number of conj classes of G=PSL(2,p) and mi is the
        multiplicity of pi_i in the ramification module of a
        Hurwitz curve with automorphism group G.
        Here IrrRepns(G) = [pi_1,...,pi_n] (in the order listed in the
        output of self.character_table()).

        REFERENCE:

        David Joyner, Amy Ksir, Roger Vogeler,
        "Group representations on Riemann-Roch spaces of some
        Hurwitz curves," preprint, 2006.

        EXAMPLES::

            sage: G = PSL(2,13)
            sage: G.ramification_module_decomposition_hurwitz_curve()  # random, optional - gap_packages
            [0, 7, 7, 12, 12, 12, 13, 15, 14]

        This means, for example, that the trivial representation does not
        occur in the ramification module of a Hurwitz curve with automorphism
        group PSL(2,13), since the trivial representation is listed first
        and that entry has multiplicity 0. The "randomness" is due to the
        fact that GAP randomly orders the conjugacy classes of the same order
        in the list of all conjugacy classes. Similarly, there is some
        randomness to the ordering of the characters.

        If you try to use this function on a group PSL(2,q) where q is
        not a (smallish) "Hurwitz prime", an error message will be printed.
        """
        from sage.env import SAGE_EXTCODE

        if self.matrix_degree() != 2:
            raise ValueError("degree must be 2")

        F = self.base_ring()
        q = F.order()
        libgap.Read(str(Path(SAGE_EXTCODE) / 'gap' / 'joyner' /
                        'hurwitz_crv_rr_sp.gap'))
        mults = libgap.eval(f"ram_module_hurwitz({q})")
        return mults.sage()

    def ramification_module_decomposition_modular_curve(self):
        r"""
        Helps compute the decomposition of the ramification module
        for the modular curve X(p) (over CC say) with automorphism group G = PSL(2,q),
        q a prime > 5. Using this computation and Borne's formula helps determine the
        G-module structure of the RR spaces of equivariant
        divisors can be determined explicitly.

        The output is a list of integer multiplicities: [m1,...,mn],
        where n is the number of conj classes of G=PSL(2,p) and mi is the
        multiplicity of pi_i in the ramification module of a
        modular curve with automorphism group G.
        Here IrrRepns(G) = [pi_1,...,pi_n] (in the order listed in the
        output of self.character_table()).

        REFERENCE: D. Joyner and A. Ksir, 'Modular representations
                   on some Riemann-Roch spaces of modular curves
                   `X(N)`', Computational Aspects of Algebraic Curves,
                   (Editor: T. Shaska) Lecture Notes in Computing, WorldScientific,
                   2005.)

        EXAMPLES::

            sage: G = PSL(2,7)
            sage: G.ramification_module_decomposition_modular_curve()  # random, optional - gap_packages
            [0, 4, 3, 6, 7, 8]

        This means, for example, that the trivial representation does not
        occur in the ramification module of X(7), since the trivial representation
        is listed first and that entry has multiplicity 0. The "randomness" is due to the
        fact that GAP randomly orders the conjugacy classes of the same order
        in the list of all conjugacy classes. Similarly, there is some
        randomness to the ordering of the characters.
        """
        from sage.env import SAGE_EXTCODE
        if self.matrix_degree() != 2:
            raise ValueError("degree must be 2")
        F = self.base_ring()
        q = F.order()
        libgap.Read(str(Path(SAGE_EXTCODE) / 'gap' / 'joyner' /
                        'modular_crv_rr_sp.gap'))
        mults = libgap.eval(f"ram_module_X({q})")
        return mults.sage()


class PSp(PermutationGroup_plg):
    def __init__(self, n, q, name='a'):
        r"""
        The projective symplectic linear groups over `\GF(q)`.

        INPUT:

        - ``n`` -- positive integer; the degree
        - ``q`` -- prime power; the size of the ground field
        - ``name`` -- (default: ``'a'``) variable name of indeterminate of finite field `\GF(q)`

        OUTPUT: PSp(`n`, `q`)

        .. NOTE::

            This group is also available via ``groups.permutation.PSp()``.

        EXAMPLES::

            sage: G = PSp(2,3); G
            Permutation Group with generators [(2,3,4), (1,2)(3,4)]
            sage: G.order()
            12
            sage: G = PSp(4,3); G
            Permutation Group with generators [(3,4)(6,7)(9,10)(12,13)(17,20)(18,21)(19,22)(23,32)(24,33)(25,34)(26,38)(27,39)(28,40)(29,35)(30,36)(31,37), (1,5,14,17,27,22,19,36,3)(2,6,32)(4,7,23,20,37,13,16,26,40)(8,24,29,30,39,10,33,11,34)(9,15,35)(12,25,38)(21,28,31)]
            sage: G.order()
            25920
            sage: print(G)
            The projective symplectic linear group of degree 4 over Finite Field of size 3
            sage: G.base_ring()
            Finite Field of size 3

            sage: G = PSp(2, 8, name='alpha'); G                                        # needs sage.rings.finite_rings
            Permutation Group with generators [(3,8,6,4,9,7,5), (1,2,3)(4,7,5)(6,9,8)]
            sage: G.base_ring()
            Finite Field in alpha of size 2^3

        TESTS::

            sage: groups.permutation.PSp(2, 3)
            Permutation Group with generators [(2,3,4), (1,2)(3,4)]
        """
        if n % 2:
            raise TypeError("the degree n must be even")
        else:
            id = 'PSp(%s,%s)' % (n, q)
        PermutationGroup_generic.__init__(self, gap_group=id)
        self._q = q
        self._base_ring = GF(q, name=name)
        self._n = n

    def __str__(self):
        """
        EXAMPLES::

            sage: G = PSp(4,3)
            sage: print(G)
            The projective symplectic linear group of degree 4 over Finite Field of size 3
        """
        return "The projective symplectic linear group of degree %s over %s" % (self._n, self.base_ring())


PSP = PSp


class PermutationGroup_pug(PermutationGroup_plg):
    def field_of_definition(self):
        """
        EXAMPLES::

            sage: PSU(2,3).field_of_definition()                                        # needs sage.rings.finite_rings
            Finite Field in a of size 3^2
        """
        return self._field_of_definition


class PSU(PermutationGroup_pug):
    def __init__(self, n, q, name='a'):
        r"""
        The projective special unitary groups over `\GF(q)`.

        INPUT:

        - ``n`` -- positive integer; the degree
        - ``q`` -- prime power; the size of the ground field
        - ``name`` -- (default: ``'a'``) variable name of indeterminate of finite
          field `\GF(q)`

        OUTPUT: PSU(`n`,`q`)

        .. NOTE::

            This group is also available via ``groups.permutation.PSU()``.

        EXAMPLES::

            sage: PSU(2,3)                                                              # needs sage.rings.finite_rings
            The projective special unitary group of degree 2 over Finite Field of size 3

            sage: G = PSU(2, 8, name='alpha'); G                                        # needs sage.rings.finite_rings
            The projective special unitary group of degree 2 over Finite Field in alpha of size 2^3
            sage: G.base_ring()                                                         # needs sage.rings.finite_rings
            Finite Field in alpha of size 2^3

        TESTS::

            sage: groups.permutation.PSU(2, 3)                                          # needs sage.rings.finite_rings
            The projective special unitary group of degree 2 over Finite Field of size 3
        """
        id = 'PSU(%s,%s)' % (n, q)
        PermutationGroup_generic.__init__(self, gap_group=id)
        self._q = q
        self._base_ring = GF(q, name=name)
        self._field_of_definition = GF(q**2, name)
        self._n = n

    def _repr_(self):
        """
        EXAMPLES::

            sage: PSU(2,3)                                                              # needs sage.rings.finite_rings
            The projective special unitary group of degree 2 over Finite Field of size 3
        """
        return "The projective special unitary group of degree %s over %s" % (self._n, self.base_ring())


class PGU(PermutationGroup_pug):
    def __init__(self, n, q, name='a'):
        r"""
        The projective general unitary groups over `\GF(q)`.

        INPUT:

        - ``n`` -- positive integer; the degree
        - ``q`` -- prime power; the size of the ground field
        - ``name`` -- (default: ``'a'``) variable name of indeterminate of finite field `\GF(q)`

        OUTPUT: PGU(`n`, `q`)

        .. NOTE::

            This group is also available via ``groups.permutation.PGU()``.

        EXAMPLES::

            sage: PGU(2,3)                                                              # needs sage.rings.finite_rings
            The projective general unitary group of degree 2 over Finite Field of size 3

            sage: G = PGU(2, 8, name='alpha'); G                                        # needs sage.rings.finite_rings
            The projective general unitary group of degree 2
             over Finite Field in alpha of size 2^3
            sage: G.base_ring()                                                         # needs sage.rings.finite_rings
            Finite Field in alpha of size 2^3

        TESTS::

            sage: groups.permutation.PGU(2, 3)                                          # needs sage.rings.finite_rings
            The projective general unitary group of degree 2 over Finite Field of size 3
        """
        id = 'PGU(%s,%s)' % (n, q)
        PermutationGroup_generic.__init__(self, gap_group=id)
        self._q = q
        self._base_ring = GF(q, name=name)
        self._field_of_definition = GF(q**2, name)
        self._n = n

    def _repr_(self):
        """
        EXAMPLES::

            sage: PGU(2,3)                                                              # needs sage.rings.finite_rings
            The projective general unitary group of degree 2 over Finite Field of size 3
        """
        return "The projective general unitary group of degree %s over %s" % (self._n, self.base_ring())


class SuzukiGroup(PermutationGroup_unique):
    def __init__(self, q, name='a'):
        r"""
        The Suzuki group over `\GF(q)`,
        `^2 B_2(2^{2k+1}) = Sz(2^{2k+1})`.

        A wrapper for the GAP function ``SuzukiGroup``.

        INPUT:

        - ``q`` -- `2^n`, an odd power of 2; the size of the ground
          field. (Strictly speaking, `n` should be greater than 1, or
          else this group is not simple.)
        - ``name`` -- (default: ``'a'``) variable name of indeterminate of
          finite field `\GF(q)`

        OUTPUT: a Suzuki group

        .. NOTE::

            This group is also available via ``groups.permutation.Suzuki()``.

        EXAMPLES::

            sage: SuzukiGroup(8)                                                        # needs sage.rings.finite_rings
            Permutation Group with generators [(1,28,10,44)(3,50,11,42)(4,43,53,64)(5,9,39,52)(6,36,63,13)(7,51,60,57)(8,33,37,16)(12,24,55,29)(14,30,48,47)(15,19,61,54)(17,59,22,62)(18,23,34,31)(20,38,49,25)(21,26,45,58)(27,32,41,65)(35,46,40,56),
             (1,2)(3,10)(4,42)(5,18)(6,50)(7,26)(8,58)(9,34)(12,28)(13,45)(14,44)(15,23)(16,31)(17,21)(19,39)(20,38)(22,25)(24,61)(27,60)(29,65)(30,55)(32,33)(35,52)(36,49)(37,59)(40,54)(41,62)(43,53)(46,48)(47,56)(51,63)(57,64)]
            sage: print(SuzukiGroup(8))                                                 # needs sage.rings.finite_rings
            The Suzuki group over Finite Field in a of size 2^3

            sage: # needs sage.rings.finite_rings
            sage: G = SuzukiGroup(32, name='alpha')
            sage: G.order()
            32537600
            sage: G.order().factor()
            2^10 * 5^2 * 31 * 41
            sage: G.base_ring()
            Finite Field in alpha of size 2^5

        TESTS::

            sage: groups.permutation.Suzuki(8)                                          # needs sage.rings.finite_rings
            Permutation Group with generators [(1,28,10,44)(3,50,11,42)(4,43,53,64)(5,9,39,52)(6,36,63,13)(7,51,60,57)(8,33,37,16)(12,24,55,29)(14,30,48,47)(15,19,61,54)(17,59,22,62)(18,23,34,31)(20,38,49,25)(21,26,45,58)(27,32,41,65)(35,46,40,56),
             (1,2)(3,10)(4,42)(5,18)(6,50)(7,26)(8,58)(9,34)(12,28)(13,45)(14,44)(15,23)(16,31)(17,21)(19,39)(20,38)(22,25)(24,61)(27,60)(29,65)(30,55)(32,33)(35,52)(36,49)(37,59)(40,54)(41,62)(43,53)(46,48)(47,56)(51,63)(57,64)]

        REFERENCES:

        -  :wikipedia:`Group_of_Lie_type\#Suzuki-Ree_groups`
        """
        q = Integer(q)
        t = valuation(q, 2)
        if 2**t != q or is_even(t):
            raise ValueError("the ground field size %s must be an odd power of 2" % q)
        id = 'SuzukiGroup(IsPermGroup,%s)' % q
        PermutationGroup_generic.__init__(self, gap_group=id)
        self._q = q
        self._base_ring = GF(q, name=name)

    def base_ring(self):
        """
        EXAMPLES::

            sage: G = SuzukiGroup(32, name='alpha')                                     # needs sage.rings.finite_rings
            sage: G.base_ring()                                                         # needs sage.rings.finite_rings
            Finite Field in alpha of size 2^5
        """
        return self._base_ring

    def __str__(self):
        """
        EXAMPLES::

            sage: G = SuzukiGroup(32, name='alpha')                                     # needs sage.rings.finite_rings
            sage: print(G)                                                              # needs sage.rings.finite_rings
            The Suzuki group over Finite Field in alpha of size 2^5
        """
        return "The Suzuki group over %s" % self.base_ring()


class ComplexReflectionGroup(PermutationGroup_unique):
    r"""
    A finite complex reflection group as a permutation group.

    We can realize `G(m,1,n)` as `m` copies of the symmetric group
    `S_n` with `s_i` for `1 \leq i < n` acting as the usual adjacent
    transposition on each copy of `S_n`. We construct the cycle
    `s_n = (n, 2n, \ldots, mn)`.

    We construct `G(m,p,n)` as a subgroup of `G(m,1,n)` by
    `s_i \mapsto s_i` for all `1 \leq i < n`,

    .. MATH::

        s_n \mapsto s_n^{-1} s_{n-1} s_n,
        \qquad\qquad
        s_{n+1} \mapsto s_n^p.

    Note that if `p = m`, then `s_{n+1} = 1`, in which case we
    do not consider it as a generator.

    The exceptional complex reflection groups `G_m` (in the Shephard-Todd
    classification) are not yet implemented.

    INPUT:

    One of the following:

    - ``m``, ``p``, ``n`` -- positive integers to construct `G(m,p,n)`

    - ``m`` -- integer such that `4 \leq m \leq 37` to construct
      an exceptional complex reflection `G_m`

    .. NOTE::

        This group is also available via
        ``groups.permutation.ComplexReflection()``.

    .. NOTE::

        The convention for the index set is for `G(m, 1, n)` to
        have the complex reflection of order `m` correspond to
        `s_n`; i.e., `s_n^m = 1` and `s_i^2 = 1` for all `i < m`.

    EXAMPLES::

        sage: G = groups.permutation.ComplexReflection(3, 1, 5)
        sage: G.order()
        29160
        sage: G
        Complex reflection group G(3, 1, 5) as a permutation group
        sage: G.category()
        Join of Category of finite enumerated permutation groups
            and Category of finite complex reflection groups

        sage: G = groups.permutation.ComplexReflection(3, 3, 4)
        sage: G.cardinality()
        648
        sage: s1, s2, s3, s4 = G.simple_reflections()
        sage: s4*s2*s4 == s2*s4*s2
        True
        sage: (s4*s3*s2)^2 == (s2*s4*s3)^2
        True

        sage: G = groups.permutation.ComplexReflection(6, 2, 3)
        sage: G.cardinality()
        648
        sage: s1, s2, s3, s4 = G.simple_reflections()
        sage: s3^2 == G.one()
        True
        sage: s4^3 == G.one()
        True
        sage: s4 * s3 * s2 == s3 * s2 * s4
        True
        sage: (s3*s2*s1)^2 == (s1*s3*s2)^2
        True
        sage: s3 * s1 * s3 == s1 * s3 * s1
        True
        sage: s4 * s3 * (s2*s3)^(2-1) == s2 * s4
        True

        sage: G = groups.permutation.ComplexReflection(4, 2, 5)
        sage: G.cardinality()
        61440

        sage: G = groups.permutation.ComplexReflection(4)
        Traceback (most recent call last):
        ...
        NotImplementedError: exceptional complex reflection groups are not yet implemented

    REFERENCES:

    - :wikipedia:`Complex_reflection_group`
    """
    def __init__(self, m, p=None, n=None):
        """
        Initialize ``self``.

        TESTS::

            sage: G = groups.permutation.ComplexReflection(3, 2, 5)
            Traceback (most recent call last):
            ...
            ValueError: p (=2) must divide m (=3)

            sage: G = groups.permutation.ComplexReflection(3, 1, 5)
            sage: TestSuite(G).run()  # long time

            sage: G = groups.permutation.ComplexReflection(3, 3, 4)
            sage: TestSuite(G).run()

            sage: G = groups.permutation.ComplexReflection(4, 2, 5)
            sage: TestSuite(G).run()  # long time

            sage: groups.permutation.ComplexReflection(4, 1, 3).cardinality()
            384
            sage: CP = ColoredPermutations(4, 3)
            sage: CP.cardinality()
            384
        """
        if p is None:
            raise NotImplementedError("exceptional complex reflection groups"
                                      " are not yet implemented")
        self._m = Integer(m)
        self._p = Integer(p)
        self._n = Integer(n)
        if min(self._m, self._p, self._n) < 1:
            raise ValueError("m (=%s) p (=%s) n (=%s) must all be >= 1" % (self._m, self._p, self._n))
        if self._m % self._p != 0:
            raise ValueError("p (=%s) must divide m (=%s)" % (self._p, self._m))

        from sage.categories.finite_complex_reflection_groups import FiniteComplexReflectionGroups
        from sage.categories.finite_permutation_groups import FinitePermutationGroups
        cat = FinitePermutationGroups() & FiniteComplexReflectionGroups()

        gens = [[(i+k, i+1+k) for k in range(0, self._m * self._n, self._n)]
                for i in range(1, self._n)]
        if self._p == 1:
            gens.append([tuple(range(self._n, self._m*self._n + 1, self._n))])
        else:
            from sage.groups.perm_gps.constructor import PermutationGroupElement
            sn = PermutationGroupElement([tuple(range(self._n, self._m*self._n + 1, self._n))])
            if self._n > 1:
                snm = PermutationGroupElement(gens[-1])
                gens.append(sn**(self._m-1) * snm * sn)
            if self._p != self._m:
                gens.append(sn**self._p)

        PermutationGroup_generic.__init__(self, gens, category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: groups.permutation.ComplexReflection(3, 1, 5)
            Complex reflection group G(3, 1, 5) as a permutation group
        """
        return "Complex reflection group G(%s, %s, %s) as a permutation group" % (self._m, self._p, self._n)

    def _latex_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.ComplexReflection(3, 1, 5)
            sage: latex(G)
            G(3,1,5)
        """
        return "G(%s,%s,%s)" % (self._m, self._p, self._n)

    @cached_method
    def index_set(self):
        r"""
        Return the index set of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.ComplexReflection(4, 1, 3)
            sage: G.index_set()
            (1, 2, 3)

            sage: G = groups.permutation.ComplexReflection(1, 1, 3)
            sage: G.index_set()
            (1, 2)

            sage: G = groups.permutation.ComplexReflection(4, 2, 3)
            sage: G.index_set()
            (1, 2, 3, 4)

            sage: G = groups.permutation.ComplexReflection(4, 4, 3)
            sage: G.index_set()
            (1, 2, 3)

            sage: G = groups.permutation.ComplexReflection(10, 2, 1)
            sage: G.index_set()
            (1,)

            sage: G = groups.permutation.ComplexReflection(10, 10, 1)
            sage: G.index_set()
            ()
        """
        n = self._n
        if self._m != 1 and self._n > 1:
            n += 1
        if self._p != 1 and self._p != self._m:
            n += 1
        return tuple(range(1, n))

    def simple_reflection(self, i):
        r"""
        Return the ``i``-th simple reflection of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.ComplexReflection(3, 1, 4)
            sage: G.simple_reflection(2)
            (2,3)(6,7)(10,11)
            sage: G.simple_reflection(4)
            (4,8,12)

            sage: G = groups.permutation.ComplexReflection(1, 1, 4)
            sage: G.simple_reflections()
            Finite family {1: (1,2), 2: (2,3), 3: (3,4)}

            sage: G = groups.permutation.ComplexReflection(10, 2, 1)
            sage: G.simple_reflections()
            Finite family {1: (1,3,5,7,9)(2,4,6,8,10)}

            sage: G = groups.permutation.ComplexReflection(10, 10, 1)
            sage: G.simple_reflections()
            Finite family {}
        """
        if i not in self.index_set():
            raise ValueError("not an index of a simple reflection")

        if i < self._n:
            return self([(i+k, i+1+k) for k in range(0, self._m * self._n, self._n)])
        if self._p == 1:
            return self([tuple(range(self._n, self._m*self._n + 1, self._n))])

        from sage.groups.perm_gps.constructor import PermutationGroupElement
        sn = PermutationGroupElement([tuple(range(self._n, self._m*self._n + 1, self._n))])
        if i == self._n + 1 or self._n == 1:
            return self(sn**self._p)

        snm = PermutationGroupElement([(self._n-1+k, self._n+k)
                                       for k in range(0, self._m * self._n, self._n)])
        return self(sn**(self._m-1) * snm * sn)

    def degrees(self):
        r"""
        Return the degrees of ``self``.

        The degrees of a complex reflection group are the degrees of
        the fundamental invariants of the ring of polynomial invariants.

        If `m = 1`, then we are in the special case of the symmetric group
        and the degrees are `(2, 3, \ldots, n, n+1)`. Otherwise the degrees
        are `(m, 2m, \ldots, (n-1)m, nm/p)`.

        EXAMPLES::

            sage: C = groups.permutation.ComplexReflection(4, 1, 3)
            sage: C.degrees()
            (4, 8, 12)
            sage: G = groups.permutation.ComplexReflection(4, 2, 3)
            sage: G.degrees()
            (4, 6, 8)
            sage: Gp = groups.permutation.ComplexReflection(4, 4, 3)
            sage: Gp.degrees()
            (3, 4, 8)
            sage: S = groups.permutation.ComplexReflection(1, 1, 3)
            sage: S.degrees()
            (2, 3)

        Check that the product of the degrees is equal to the
        cardinality::

            sage: prod(C.degrees()) == C.cardinality()
            True
            sage: prod(G.degrees()) == G.cardinality()
            True
            sage: prod(Gp.degrees()) == Gp.cardinality()
            True
            sage: prod(S.degrees()) == S.cardinality()
            True
        """
        # For the usual symmetric group (self._m=1) we need to start at 2
        start = 2 if self._m == 1 else 1
        ret = [self._m * i for i in range(start, self._n)]
        ret.append(self._m * self._n // self._p)
        return tuple(sorted(ret))

    def codegrees(self):
        r"""
        Return the codegrees of ``self``.

        Let `G` be a complex reflection group. The codegrees
        `d_1^* \leq d_2^* \leq \cdots \leq d_{\ell}^*` of `G` can be
        defined by:

        .. MATH::

            \prod_{i=1}^{\ell} (q - d_i^* - 1)
            = \sum_{g \in G} \det(g) q^{\dim(V^g)},

        where `V` is the natural complex vector space that `G` acts on
        and `\ell` is the :meth:`rank`.

        If `m = 1`, then we are in the special case of the symmetric group
        and the codegrees are `(n-2, n-3, \ldots 1, 0)`. Otherwise the
        codegrees are `((n-1)m, (n-2)m, \ldots, m, 0)`.

        EXAMPLES::

            sage: C = groups.permutation.ComplexReflection(4, 1, 3)
            sage: C.codegrees()
            (8, 4, 0)
            sage: G = groups.permutation.ComplexReflection(3, 3, 4)
            sage: G.codegrees()
            (6, 5, 3, 0)
            sage: S = groups.permutation.ComplexReflection(1, 1, 3)
            sage: S.codegrees()
            (1, 0)
        """
        # Special case for the usual symmetric group
        if self._m == 1:
            return tuple(reversed(range(self._n-1)))
        if self._p < self._m:
            return tuple([self._m * i for i in reversed(range(self._n))])
        ret = [self._m * i for i in reversed(range(self._n-1))]
        ret.append((self._n-1)*self._m - self._n)
        return tuple(sorted(ret, reverse=True))


class SmallPermutationGroup(PermutationGroup_generic):
    r"""
    A GAP SmallGroup, returned as a permutation group.

    GAP contains a library SGL of small groups, each identified by
    its GAP SmallGroup id. (MAGMA uses the same identifiers).
    The GAP SmallGroup id is a pair ``[n,k]`` consisting of
    ``n``, the order of the group, and ``k``, an index determining
    the group specifically. This class can construct the group as a
    permutation group from this data.

    INPUT:

    - ``order`` -- the order of the group

    - ``gap_id`` -- the numerical index in the GAP id of the group

    Generators may be obtained through the :meth:`gens` method.
    These could change for a particular group in later releases
    of GAP. In many instances the degree of the constructed group
    ``SmallPermutationGroup(n,k)`` will be a permutation group on
    `n` letters, but this will not always be true.

    EXAMPLES::

        sage: G = SmallPermutationGroup(12,4); G
        Group of order 12 and GAP Id 4 as a permutation group
        sage: G.gens()
        ((1,2)(3,5)(4,10)(6,8)(7,12)(9,11),
         (1,3)(2,5)(4,7)(6,9)(8,11)(10,12),
         (1,4,8)(2,6,10)(3,7,11)(5,9,12))
        sage: G.character_table()                                                       # needs sage.rings.number_field
        [ 1  1  1  1  1  1]
        [ 1 -1 -1  1  1 -1]
        [ 1 -1  1  1 -1  1]
        [ 1  1 -1  1 -1 -1]
        [ 2  0 -2 -1  0  1]
        [ 2  0  2 -1  0 -1]
        sage: def numgps(n): return ZZ(libgap.NumberSmallGroups(n))
        sage: all(SmallPermutationGroup(n,k).id() == [n,k]
        ....:     for n in [1..64] for k in [1..numgps(n)])
        True
        sage: H = SmallPermutationGroup(6,1)
        sage: H.is_abelian()
        False
        sage: [H.centralizer(g) for g in H.conjugacy_classes_representatives()]
        [Subgroup generated by [(1,2)(3,6)(4,5), (1,3,5)(2,4,6)] of
          (Group of order 6 and GAP Id 1 as a permutation group),
         Subgroup generated by [(1,2)(3,6)(4,5)] of
          (Group of order 6 and GAP Id 1 as a permutation group),
         Subgroup generated by [(1,3,5)(2,4,6), (1,5,3)(2,6,4)] of
          (Group of order 6 and GAP Id 1 as a permutation group)]
    """

    def __init__(self, order, gap_id):
        """
        Initialize ``self``.

        TESTS::

            sage: TestSuite(SmallPermutationGroup(60,5)).run()
        """
        self._n = order
        self._gap_id = gap_id
        self._gap_small_group = libgap.SmallGroup(order, gap_id)
        gap_permutation_group = self._gap_small_group.IsomorphismPermGroup().Image(self._gap_small_group)
        PermutationGroup_generic.__init__(self, gap_group=gap_permutation_group)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: G = SmallPermutationGroup(12,4); G
            Group of order 12 and GAP Id 4 as a permutation group
        """
        return "Group of order %s and GAP Id %s as a permutation group" % (self._n, self._gap_id)

    def order(self):
        """
        Return the order of the group corresponding to ``self``.

        EXAMPLES::

            sage: [SmallPermutationGroup(21,k).order() for k in [1,2]]
            [21, 21]
        """
        return self._n

    def gap_small_group(self):
        r"""
        Return the GAP small group object corresponding to ``self``.

        GAP realizes some small groups as ``PermutationGroup``, others as ``PcGroups``
        (polycyclic groups). The :class:`SmallPermutationGroup` class always
        returns a ``PermutationGroup``, but in the process of creating this group
        a GAP ``SmallGroup`` is generated. This method returns that group.

        EXAMPLES::

            sage: SmallPermutationGroup(168,41).gap_small_group()
            <pc group of size 168 with 5 generators>
            sage: SmallPermutationGroup(168,42).gap_small_group()
            Group([ (3,4)(5,6), (1,2,3)(4,5,7) ])
        """
        return self._gap_small_group
