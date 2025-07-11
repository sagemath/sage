# sage_setup: distribution = sagemath-coxeter3
# sage.doctest: optional - coxeter3
"""
Coxeter Groups implemented with Coxeter3
"""
# ****************************************************************************
#       Copyright (C) 2009-2013 Mike Hansen <mhansen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.libs.coxeter3.coxeter import get_CoxGroup, CoxGroupElement
from sage.misc.cachefunc import cached_method

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.richcmp import richcmp
from sage.categories.coxeter_groups import CoxeterGroups
from sage.structure.parent import Parent

from sage.combinat.root_system.coxeter_matrix import CoxeterMatrix

from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class CoxeterGroup(UniqueRepresentation, Parent):
    @staticmethod
    def __classcall__(cls, cartan_type, *args, **options):
        """
        TESTS::

            sage: from sage.libs.coxeter3.coxeter_group import CoxeterGroup
            sage: CoxeterGroup(['B',2])
            Coxeter group of type ['B', 2] implemented by Coxeter3
            sage: CoxeterGroup(CartanType(['B', 3]).relabel({1: 3, 2: 2, 3: 1}))
            Coxeter group of type ['B', 3] relabelled by {1: 3, 2: 2, 3: 1} implemented by Coxeter3
        """
        from sage.combinat.root_system.cartan_type import CartanType
        ct = CartanType(cartan_type)
        return super().__classcall__(cls, ct, *args, **options)

    def __init__(self, cartan_type):
        """
        TESTS::

            sage: from sage.libs.coxeter3.coxeter_group import CoxeterGroup
            sage: CoxeterGroup(['A',2])
            Coxeter group of type ['A', 2] implemented by Coxeter3

        As degrees and codegrees are not implemented, they are skipped in the
        testsuite::

            sage: to_skip = ['_test_degrees', '_test_codegrees']
            sage: TestSuite(CoxeterGroup(['A',2])).run(skip=to_skip)
        """
        category = CoxeterGroups()
        if cartan_type.is_finite():
            category = category.Finite()
        Parent.__init__(self, category=category)
        self._coxgroup = get_CoxGroup(cartan_type)
        self._cartan_type = cartan_type

    def _repr_(self):
        """
        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3'); W      # indirect doctest
            Coxeter group of type ['A', 3] implemented by Coxeter3
            sage: W = CoxeterGroup(['A', 3, 1], implementation='coxeter3'); W
            Coxeter group of type ['A', 3, 1] implemented by Coxeter3
        """
        return "Coxeter group of type %s implemented by Coxeter3" % (self.cartan_type())

    def __iter__(self):
        """
        EXAMPLES::

            sage: W = CoxeterGroup(['A', 2], implementation='coxeter3')
            sage: list(W)
            [[], [1], [2], [1, 2], [2, 1], [1, 2, 1]]
        """
        for x in self._coxgroup:
            yield CoxeterGroup.Element(self, x)

    def cartan_type(self):
        """
        Return the Cartan type for this Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: W.cartan_type()
            ['A', 3]
        """
        return self._cartan_type

    def index_set(self):
        """
        Return the index set for the generators of this Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: W.index_set()
            (1, 2, 3)
            sage: C = CoxeterGroup(['A', 3,1], implementation='coxeter3')
            sage: C.index_set()
            (0, 1, 2, 3)
        """
        return self.cartan_type().index_set()

    def bruhat_interval(self, u, v):
        """
        Return the Bruhat interval between ``u`` and ``v``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: W.bruhat_interval([1],[3,1,2,3])
            [[1], [1, 2], [1, 3], [1, 2, 3], [1, 3, 2], [1, 2, 3, 2]]
        """
        u, v = self(u), self(v)
        return self._coxgroup.bruhat_interval(u.value, v.value)

    def cardinality(self):
        """
        Return the cardinality of this Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: W.cardinality()
            24
        """
        return self._coxgroup.order()

    def one(self):
        """
        Return the identity element of this Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: W.one()
            []
        """
        return self.element_class(self, [])

    def simple_reflections(self):
        """
        Return the family of generators for this Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: s = W.simple_reflections()
            sage: s[2]*s[1]*s[2]
            [1, 2, 1]
        """
        from sage.sets.family import Family
        return Family(self.index_set(), lambda i: self.element_class(self, [i]))

    gens = simple_reflections

    def from_reduced_word(self, w):
        """
        Return an element of ``self`` from its (reduced) word.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: W.from_reduced_word([1, 3])
            [1, 3]
            sage: W.from_reduced_word([3, 1])
            [1, 3]
        """
        return self.element_class(self, w)

    def rank(self):
        """
        Return the rank of this Coxeter group, that is, the number of generators.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: W.rank()
            3
        """
        return self._coxgroup.rank()

    def is_finite(self):
        """
        Return ``True`` if this is a finite Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: W.is_finite()
            True
        """
        return self._coxgroup.is_finite()

    def length(self, x):
        """
        Return the length of an element ``x`` in this Coxeter group.
        This is just the length of a reduced word for ``x``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: W.length(W([1,2]))
            2
            sage: W.length(W([1,1]))
            0
        """
        return x.length()

    @cached_method
    def coxeter_matrix(self):
        """
        Return the Coxeter matrix for this Coxeter group.

        The columns and rows are ordered according to the result of
        :meth:`index_set`.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: m = W.coxeter_matrix(); m
            [1 3 2]
            [3 1 3]
            [2 3 1]
            sage: m.index_set() == W.index_set()
            True
        """
        return CoxeterMatrix(self._coxgroup.coxeter_matrix(), self.index_set())

    def root_system(self):
        """
        Return the root system associated with this Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: R = W.root_system(); R
            Root system of type ['A', 3]
            sage: alpha = R.root_space().basis()
            sage: alpha[2] + alpha[3]
            alpha[2] + alpha[3]
        """
        return self.cartan_type().root_system()

    def _an_element_(self):
        """
        Return an element of this Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: W._an_element_()
            []
        """
        return self.element_class(self, [])

    def m(self, i, j):
        r"""
        This is deprecated, use ``self.coxeter_matrix()[i,j]`` instead.

        TESTS::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: W.m(1, 1)
            doctest:warning...:
            DeprecationWarning: the .m(i, j) method has been deprecated; use .coxeter_matrix()[i,j] instead.
            See https://github.com/sagemath/sage/issues/30237 for details.
            1
        """
        from sage.misc.superseded import deprecation
        deprecation(30237, "the .m(i, j) method has been deprecated; use .coxeter_matrix()[i,j] instead.")
        return self.coxeter_matrix()[i, j]

    def kazhdan_lusztig_polynomial(self, u, v, constant_term_one=True):
        r"""
        Return the Kazhdan-Lusztig polynomial `P_{u,v}`.

        INPUT:

        - ``u``, ``v`` -- elements of the underlying Coxeter group
        - ``constant_term_one`` -- boolean (default: ``True``); ``True`` uses
          the constant equals one convention, ``False`` uses the Leclerc-Thibon
          convention

        .. SEEALSO::

            - :class:`~sage.combinat.kazhdan_lusztig.KazhdanLusztigPolynomial`
            - :meth:`parabolic_kazhdan_lusztig_polynomial`

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: W.kazhdan_lusztig_polynomial([], [1,2, 1])
            1
            sage: W.kazhdan_lusztig_polynomial([1],[3,2])
            0
            sage: W = CoxeterGroup(['A',3],implementation='coxeter3')
            sage: W.kazhdan_lusztig_polynomial([2],[2,1,3,2])
            q + 1

        .. NOTE::

            Coxeter3, as well as Sage's native implementation in
            :class:`~sage.combinat.kazhdan_lusztig.KazhdanLusztigPolynomial`
            use the convention under which Kazhdan-Lusztig
            polynomials give the change of basis from the `(C_w)_{w\in W}`
            basis to the `(T_w)_{w\in W}` of the Hecke algebra of `W` with
            parameters `q` and `q^{-1}`:

                .. MATH:: C_w = \sum_u  P_{u,w} T_u.

            In particular, `P_{u,u}=1`::

                sage: all(W.kazhdan_lusztig_polynomial(u,u) == 1 for u in W)
                True

            This convention differs from Theorem 2.7 in [LT1998]_ by:

            .. MATH::

                {}^{LT} P_{y,w}(q) = q^{\ell(w)-\ell(y)} P_{y,w}(q^{-2})

            To access the Leclerc-Thibon convention use::

                sage: W = CoxeterGroup(['A',3],implementation='coxeter3')
                sage: W.kazhdan_lusztig_polynomial([2],[2,1,3,2],constant_term_one=False)
                q^3 + q

        TESTS:

        We check that Coxeter3 and Sage's implementation give the same results::

            sage: C = CoxeterGroup(['B', 3], implementation='coxeter3')
            sage: W = WeylGroup("B3",prefix='s')
            sage: [s1,s2,s3] = W.simple_reflections()
            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: KL = KazhdanLusztigPolynomial(W,q)
            sage: all(KL.P(1,w) == C.kazhdan_lusztig_polynomial([],w.reduced_word()) for w in W)  # long (15s)
            True
        """
        u, v = self(u), self(v)
        p = u.value.kazhdan_lusztig_polynomial(v.value)
        if constant_term_one:
            return p
        ZZq = PolynomialRing(ZZ, 'q', sparse=True)
        # This is the same as q**len_diff * p(q**(-2))
        len_diff = v.length() - u.length()
        d = {-2 * deg + len_diff: coeff for deg, coeff in enumerate(p)
             if coeff != 0}
        return ZZq(d)

    def parabolic_kazhdan_lusztig_polynomial(self, u, v, J, constant_term_one=True):
        r"""
        Return the parabolic Kazhdan-Lusztig polynomial `P_{u,v}^{-,J}`.

        INPUT:

        - ``u``, ``v`` -- minimal length coset representatives of `W/W_J` for this Coxeter group `W`
        - ``J`` -- a subset of the index set of ``self`` specifying the parabolic subgroup

        This method implements the parabolic Kazhdan-Lusztig polynomials
        `P^{-,J}_{u,v}` of [Deo1987b]_, which are defined as
        `P^{-,J}_{u,v} = \sum_{z\in W_J} (-1)^{\ell(z)} P_{yz,w}(q)`
        with the conventions in Sage.
        As for :meth:`kazhdan_lusztig_polynomial` the convention
        differs from Theorem 2.7 in [LT1998]_ by:

        .. MATH::

            {}^{LT} P_{y,w}^{-,J}(q) = q^{\ell(w)-\ell(y)} P_{y,w}^{-,J}(q^{-2})

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], implementation='coxeter3')
            sage: W.parabolic_kazhdan_lusztig_polynomial([],[3,2],[1,3])
            0
            sage: W.parabolic_kazhdan_lusztig_polynomial([2],[2,1,3,2],[1,3])
            q

            sage: C = CoxeterGroup(['A',3,1], implementation='coxeter3')
            sage: C.parabolic_kazhdan_lusztig_polynomial([],[1],[0])
            1
            sage: C.parabolic_kazhdan_lusztig_polynomial([],[1,2,1],[0])
            1
            sage: C.parabolic_kazhdan_lusztig_polynomial([],[0,1,0,1,2,1],[0])
            q
            sage: w=[1, 2, 1, 3, 0, 2, 1, 0, 3, 0, 2]
            sage: v=[1, 2, 1, 3, 0, 1, 2, 1, 0, 3, 0, 2, 1, 0, 3, 0, 2]
            sage: C.parabolic_kazhdan_lusztig_polynomial(w,v,[1,3])
            q^2 + q
            sage: C.parabolic_kazhdan_lusztig_polynomial(w,v,[1,3],constant_term_one=False)
            q^4 + q^2

        TESTS::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
            sage: type(W.parabolic_kazhdan_lusztig_polynomial([2],[],[1]))
            <class 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>
        """
        u = self(u)
        v = self(v)
        if any(d in J for d in u.descents()) or any(d in J for d in v.descents()):
            raise ValueError("u and v have to be minimal coset representatives")
        J_set = set(J)
        WOI = self.weak_order_ideal(lambda x: J_set.issuperset(x.descents()))
        if constant_term_one:
            P = PolynomialRing(ZZ, 'q')
            return P.sum((-1)**(z.length()) * self.kazhdan_lusztig_polynomial(u * z, v)
                         for z in WOI if (u * z).bruhat_le(v))
        P = PolynomialRing(ZZ, 'q', sparse=True)
        return P.sum((-1)**(z.length()) * self.kazhdan_lusztig_polynomial(u * z, v, constant_term_one=False).shift(z.length())
                     for z in WOI if (u * z).bruhat_le(v))

    class Element(ElementWrapper):
        wrapped_class = CoxGroupElement

        def __init__(self, parent, x):
            """
            TESTS::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
                sage: W([2,1,2])
                [1, 2, 1]

            Check that :issue:`32266` is fixed::

                sage: A3 = CoxeterGroup('A3', implementation='coxeter3')
                sage: s1,s2,s3 = A3.simple_reflections()
                sage: s1*s3
                [1, 3]
                sage: s3*s1
                [1, 3]
                sage: s3*s1 == s1*s3
                True
            """
            if not isinstance(x, CoxGroupElement):
                x = CoxGroupElement(parent._coxgroup, x).reduced()
            x = x.normal_form()
            ElementWrapper.__init__(self, parent, x)

        def __iter__(self):
            """
            Return an iterator for the elements in the reduced word.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
                sage: w = W([1,2,1])
                sage: list(iter(w))
                [1, 2, 1]
            """
            return iter(self.value)

        def coatoms(self):
            """
            Return the coatoms (or co-covers) of this element in the Bruhat order.

            EXAMPLES::

                sage: W = CoxeterGroup(['B', 3], implementation='coxeter3')
                sage: w = W([1,2,3])
                sage: w.coatoms()
                [[2, 3], [3, 1], [1, 2]]
            """
            W = self.parent()
            return [W(w) for w in self.value.coatoms()]

        def _richcmp_(self, other, op):
            """
            Return lexicographic comparison of ``self`` and ``other``.

            EXAMPLES::

                sage: W = CoxeterGroup(['B', 3], implementation='coxeter3')
                sage: w = W([1,2,3])
                sage: v = W([3,1,2])
                sage: v < w
                False
                sage: w < v
                True

            Some tests for equality::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
                sage: W([1,2,1]) == W([2,1,2])
                True
                sage: W([1,2,1]) == W([2,1])
                False
            """
            return richcmp(list(self), list(other), op)

        def reduced_word(self):
            """
            Return the reduced word of ``self``.

            EXAMPLES::

                sage: W = CoxeterGroup(['B', 3], implementation='coxeter3')
                sage: w = W([1,2,3])
                sage: w.reduced_word()
                [1, 2, 3]
            """
            return list(self)

        def __invert__(self):
            """
            Return the inverse of this Coxeter group element.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
                sage: w = W([1,2,3])
                sage: ~w
                [3, 2, 1]
            """
            return self.__class__(self.parent(), ~self.value)

        inverse = __invert__

        def __getitem__(self, i):
            """
            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
                sage: w0 = W([1,2,1])
                sage: w0[0]
                1
                sage: w0[1]
                2
            """
            # Allow the error message to be raised by the underlying element
            return self.value[i]

        def _mul_(self, y):
            """
            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
                sage: s = W.gens()
                sage: s[1]._mul_(s[1])
                []
                sage: s[1]*s[2]*s[1]
                [1, 2, 1]
                sage: s[2]*s[1]*s[2]
                [1, 2, 1]
            """
            return self.__class__(self.parent(), self.value * y.value)

        def __len__(self):
            """
            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
                sage: w = W([1,2,1])
                sage: w.length()
                3
                sage: len(w)
                3
            """
            return len(self.value)

        length = __len__

        def bruhat_le(self, v):
            r"""
            Return whether ``self`` `\le` ``v`` in Bruhat order.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
                sage: W([]).bruhat_le([1,2,1])
                True
            """
            v = self.parent()(v)
            return self.value.bruhat_le(v.value)

        def poincare_polynomial(self):
            """
            Return the Poincaré polynomial associated with this element.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 2], implementation='coxeter3')
                sage: W.long_element().poincare_polynomial()
                t^3 + 2*t^2 + 2*t + 1
                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
                sage: W([2,1,3,2]).poincare_polynomial()
                t^4 + 4*t^3 + 5*t^2 + 3*t + 1
                sage: W([1,2,3,2,1]).poincare_polynomial()
                t^5 + 4*t^4 + 6*t^3 + 5*t^2 + 3*t + 1
                sage: rw = sage.combinat.permutation.from_reduced_word
                sage: p = [w.poincare_polynomial() for w in W]
                sage: [rw(w.reduced_word()) for i,w in enumerate(W) if p[i] != p[i].reverse()]
                [[3, 4, 1, 2], [4, 2, 3, 1]]
            """
            return self.value.poincare_polynomial()

        def has_right_descent(self, i) -> bool:
            """
            Return whether ``i`` is a right descent of this element.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 4], implementation='coxeter3')
                sage: W([1,2]).has_right_descent(1)
                False
                sage: W([1,2]).has_right_descent(2)
                True
            """
            return i in self.value.right_descents()

        def has_left_descent(self, i) -> bool:
            """
            Return ``True`` if ``i`` is a left descent of this element.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 4], implementation='coxeter3')
                sage: W([1,2]).has_left_descent(1)
                True
                sage: W([1,2]).has_left_descent(2)
                False
            """
            return i in self.value.left_descents()

        def action(self, v):
            """
            Return the action of this Coxeter group element on the root space.

            INPUT:

            - ``v`` -- an element of the root space associated with the Coxeter group for ``self``

            EXAMPLES::

                sage: W = CoxeterGroup(['B', 3], implementation='coxeter3')
                sage: R = W.root_system().root_space()
                sage: v = R.an_element(); v
                2*alpha[1] + 2*alpha[2] + 3*alpha[3]
                sage: w = W([1,2,3])
                sage: w.action(v)
                -alpha[1] + alpha[2] + alpha[3]
            """
            # TODO: Find a better way to do this
            W = self.parent().root_system().root_space().weyl_group()
            w = W.from_reduced_word(list(self))
            return w.action(v)

        def action_on_rational_function(self, f):
            r"""
            Return the natural action of this Coxeter group element on a
            polynomial considered as an element of `S(\mathfrak{h}^*)`.

            .. NOTE::

               Note that the number of variables in the polynomial
               ring must correspond to the rank of this Coxeter
               group. The ordering of the variables is assumed to
               coincide with the result of :meth:`index_set`.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')
                sage: S = PolynomialRing(QQ, 'x,y,z').fraction_field()
                sage: x,y,z = S.gens()
                sage: W([1]).action_on_rational_function(x+y+z)
                (x^2*y + x*z + 1)/x
                sage: W([2]).action_on_rational_function(x+y+z)
                (x*y^2 + y^2*z + 1)/y
                sage: W([3]).action_on_rational_function(x+y+z)
                (y*z^2 + x*z + 1)/z
            """
            Q = f.parent()
            Q_gens = Q.gens()
            W = self.parent()
            R = W.root_system().root_space()
            alpha = R.basis()
            n = W.rank()

            if Q.ngens() != n:
                raise ValueError("the number of generators for the polynomial "
                                 "ring must be the same as the rank of the "
                                 "root system")

            basis_elements = [alpha[i] for i in W.index_set()]
            basis_to_order = {s: i for i, s in enumerate(W.index_set())}

            results = []
            for poly in [f.numerator(), f.denominator()]:
                result = 0
                exponents = poly.exponents()

                for exponent in exponents:
                    # Construct something in the root lattice from the exponent vector
                    exponent = sum(e * b for e, b in zip(exponent, basis_elements))
                    exponent = self.action(exponent)

                    monomial = 1
                    for s, c in exponent.monomial_coefficients().items():
                        monomial *= Q_gens[basis_to_order[s]]**int(c)

                    result += monomial

                results.append(result)

            numerator, denominator = results
            return numerator / denominator
