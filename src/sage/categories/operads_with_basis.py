r"""
OperadsWithBasis
"""
# ****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.operads import Operads
from sage.categories.category_types import Category_over_base_ring
from sage.categories.graded_modules_with_basis import GradedModulesWithBasis


class OperadsWithBasis(Category_over_base_ring):
    """
    The category of operads with a distinguished basis

    EXAMPLES::

        sage: OperadsWithBasis(QQ)
        Category of operads with basis over Rational Field
        sage: OperadsWithBasis(QQ).super_categories()
        [Category of graded vector spaces with basis over Rational Field,
        Category of operads over Rational Field]

    TESTS::

        sage: C = OperadsWithBasis(QQ)
        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: OperadsWithBasis(QQ).super_categories()
            [Category of graded vector spaces with basis over Rational Field,
            Category of operads over Rational Field]
        """
        R = self.base_ring()
        return [GradedModulesWithBasis(R), Operads(R)]

    def example(self):
        """
        Return an example of operad with basis.

        Here, the associative operad.

        EXAMPLES::

            sage: OperadsWithBasis(QQ).example()
            An example of an operad with basis: the Associative operad
            over Rational Field
        """
        from sage.categories.examples.operads_with_basis import Example
        return Example(self.base_ring())

    class ParentMethods:

        @cached_method
        def one(self, letter):
            """
            Return the one of the operad.

            INPUT:

            ``letter`` -- the chosen labelling.

            EXAMPLES::

                sage: A = OperadsWithBasis(QQ).example()
                sage: A.one('x')
                B[word: x]
            """
            if self.one_basis is not NotImplemented:
                return self.basis()[self.one_basis(letter)]
            return NotImplemented

        @abstract_method(optional=True)
        def one_basis(self, letter):
            """
            When the one of an operad with basis is an element of
            this basis, this optional method can return the index of
            this element.

            EXAMPLES::

                sage: A = OperadsWithBasis(QQ).example()
                sage: A.one_basis("a")
                word: a
                sage: A.one("a")
                B[word: a]
            """

        @abstract_method(optional=True)
        def composition_on_basis(self, i, j, k):
            """
            The composition of the operad on the basis (optional)

            INPUT:

            - `i`, `j` -- the indices of two elements of the basis of ``self``

            - `k` -- the composition label

            Return the composition of the two corresponding basis
            elements at the chosen label.

            If implemented, :meth:`composition` is defined from it by
            bilinearity.

            EXAMPLES::

                sage: A = OperadsWithBasis(QQ).example()
                sage: Word = A.basis().keys()
                sage: A.composition_on_basis(Word("abc"),Word("de"),"b")
                B[word: adec]
            """

        def composition_on_basis_truncated(self, i, j, k, N):
            """
            The composition of the operad on the basis (optional)
            up to degree `N` (included) only.

            INPUT:

            - `i`, `j` -- the indices of two elements of the basis of ``self``

            - `k` -- the composition label

            - `N` -- an integer, the order of truncation

            Return the composition of the two corresponding basis
            elements at the chosen label, up to order `N`.

            If implemented, :meth:`composition_truncated` is defined
            from it by bilinearity.

            EXAMPLES::

                sage: A = OperadsWithBasis(QQ).example()
                sage: Word = A.basis().keys()
                sage: A.composition_on_basis_truncated(Word("abc"),Word("de"),"b",4)
                B[word: adec]
                sage: A.composition_on_basis_truncated(Word("abc"),Word("de"),"b",3)
                0
            """
            if self.degree_on_basis(i) + self.degree_on_basis(j) - 1 <= N:
                return self.composition_on_basis(i, j, k)
            return self.zero()

        def composition_on_basis_with_numbers(self, i, j, k):
            """
            This is variant of composition where one assumes that
            labels are integers from `1` to `n`.

            The basis indices must have a method `map_labels`.

            EXAMPLES::

                sage: from sage.combinat.abstract_tree import from_hexacode
                sage: A = PreLieOperad(QQ)
                sage: RT = A.basis().keys()
                sage: x = from_hexacode('200',RT).canonical_labelling()
                sage: y = from_hexacode('10',RT).canonical_labelling()
                sage: A.composition_on_basis_with_numbers(x,y,3)
                B[1[2[], 3[4[]]]]
            """
            # the operad must define map_labels !
            shifted_j = j.map_labels(lambda z: z + k - 1)

            def shift(l):
                if l > k:
                    return l + self.degree_on_basis(j) - 1
                return l
            shifted_i = i.map_labels(shift)
            return self.composition_on_basis(shifted_i, shifted_j, k)

        @lazy_attribute
        def composition(self):
            """
            The composition of the operad, as per ``Operads.ParentMethods.composition``

            By default, this is implemented from
            :meth:`.composition_on_basis`, if available.

            EXAMPLES::

                sage: A = OperadsWithBasis(QQ).example()
                sage: Word = A.basis().keys()
                sage: A.composition(A(Word("abc")),A(Word("de")),"b")
                B[word: adec]
            """
            if self.composition_on_basis is not NotImplemented:
                return self._module_morphism(self._module_morphism(self.composition_on_basis, position=0, codomain=self), position=1)
            return NotImplemented

        @lazy_attribute
        def composition_truncated(self):
            """
            The composition of the operad, as per ``Operads.ParentMethods.composition``
            up to some degree only.

            By default, this is implemented from
            :meth:`.composition_on_basis_truncated`, if available.

            EXAMPLES::

                sage: A = OperadsWithBasis(QQ).example()
                sage: Word = A.basis().keys()
                sage: A.composition_truncated(A(Word("abc")),A(Word("de")),"b",4)
                B[word: adec]
                sage: A.composition_truncated(A(Word("abc")),A(Word("de")),"b",3)
                0
            """
            if self.composition_on_basis_truncated is not NotImplemented:
                return self._module_morphism(self._module_morphism(self.composition_on_basis_truncated, position=0, codomain=self), position=1)
            return NotImplemented

        @lazy_attribute
        def composition_with_numbers(self):
            """
            This is a variant of composition, where one assumes that
            the objects are labelled by integers from `1` to `n`.

            The result is labelled in the same way.

            The basis indices must have a method `map_labels`.

            EXAMPLES::

                sage: from sage.combinat.abstract_tree import from_hexacode
                sage: A = PreLieOperad(QQ)
                sage: RT = A.basis().keys()
                sage: x = from_hexacode('200',RT).canonical_labelling()
                sage: y = from_hexacode('10',RT).canonical_labelling()
                sage: A.composition_with_numbers(A(x),A(y),2)
                B[1[4[], 2[3[]]]]
            """
            if self.composition_on_basis_with_numbers is not NotImplemented:
                return self._module_morphism(self._module_morphism(self.composition_on_basis_with_numbers, position=0, codomain=self), position=1)
            return NotImplemented

        @abstract_method(optional=True)
        def operad_morphism_on_basis(self, i, codomain):
            """
            Morphism from the operad ``self`` to the operad `P`
            defined on the basis of ``self`` (optional)

            INPUT:

            - `i` -- the index of an element of the basis of ``self``

            - `P` -- an operad with the method for ``self``

            OUTPUT:

            an element of the operad `P`

            EXAMPLES::

                sage: PL = PreLieOperad(QQ)
                sage: PLT = PL.basis().keys()
                sage: t0 = PLT([],'a')
                sage: PL.operad_morphism_on_basis(t0,PL)
                B[a[]]
            """

        def operad_morphism(self, codomain, *args, **opts):
            """
            Morphism from the given operad, as per ``Operads.ParentMethods.operad_morphism``

            By default, this is implemented from
            :meth:`.operad_morphism_on_basis`, if available.

            EXAMPLES::

                sage: PL = PreLieOperad(QQ)
                sage: DO = DendriformOperad(QQ)
                sage: PLT = PL.basis().keys()
                sage: t0 = PL(PLT([],'a'))
                sage: PL2PL = PL.operad_morphism(PL)
                sage: PL2PL(t0)
                B[a[]]
                sage: PL2DO = PL.operad_morphism(DO)
                sage: PL2DO(t0)
                B[a[., .]]
                sage: PL2DO(t0 < t0)
                -B[a[., a[., .]]] + B[a[a[., .], .]]
            """
            if self.operad_morphism_on_basis is not NotImplemented:
                return self._module_morphism(
                    lambda x: self.operad_morphism_on_basis(x, codomain,
                                                            *args, **opts),
                    position=0, codomain=codomain)
            raise NotImplementedError

        @lazy_attribute
        def labelling(self):
            """
            This ensures that elements are labelled by integers
            between `1` and their size, in a bijective way.

            EXAMPLES::

                sage: DO = DendriformOperad(QQ)
                sage: t = 2*DO(LabelledBinaryTree([None, None],label='@'))
                sage: DO.labelling(t)
                2*B[1[., .]]
            """
            if self.labelling_on_basis is not NotImplemented:
                return self._module_morphism(
                    lambda x: self.labelling_on_basis(x), position=0,
                    codomain=self)
            raise NotImplementedError

        @lazy_attribute
        def unlabelling(self):
            """
            This ensures that basis elements are not labelled.

            EXAMPLES::

                sage: DO = DendriformOperad(QQ)
                sage: t = 2*DO(LabelledBinaryTree([None, None],label='1'))
                sage: DO.unlabelling(t)
                2*B[[., .]]
            """
            if self.unlabelling_on_basis is not NotImplemented:
                return self._module_morphism(
                    lambda x: self.unlabelling_on_basis(x), position=0, codomain=self)
            raise NotImplementedError

        def suspension(self, s, q):
            r"""
            Return the suspension of `s` with argument `q`.

            INPUT:

            - `q` -- an element of the base ring of `s`

            For the element `s = \sum_n s_n`, this is defined as
            `\sum_\n q^(n-1) s_n`.

            EXAMPLES::

                sage: PL = PreLieOperad(QQ)
                sage: o = PL.one('o')
                sage: PL.suspension(o,4)
                B[o[]]
                sage: t = (o<o)+((o<o)<o)
                sage: PL.suspension(t,2)
                2*B[o[o[]]] + 4*B[o[o[o[]]]] + 4*B[o[o[], o[]]]
            """
            def susp(i, c):
                return (i, q ** (self.degree_on_basis(i) - 1) * c)

            return s.map_item(susp)

        def group_product_with_numbers(self, s, t, N):
            """
            Return the group product of `s` by `t` up to order `N`.

            INPUT:

            - `s` and `t` are elements of the operad, assumed
              to be labelled by consecutive integers
            - `N` -- an integer, the truncation order

            EXAMPLES::

                sage: PL = PreLieOperad(QQ)
                sage: x = PL.one(1)
                sage: PL.group_product_with_numbers(x, x, 2)
                B[1[]]
            """
            dom = s.parent()
            resu = dom.zero()
            for i in range(1, N + 1):
                si = s.homogeneous_component(i)
                t_short = t.truncate(1 + N - i)
                for j in range(1, i + 1):
                    si = dom.composition_with_numbers(si, t_short, i + 1 - j).truncate(N)
                resu += si
            return resu

        def group_product_without_numbers(self, s, t, N):
            """
            Return the group product of `s` by `t` up to order `N`.

            INPUT:

            - `s` and `t` are elements of the operad, assumed to be unlabelled.

            - `N` -- an integer, the truncation order

            EXAMPLES::

                sage: PL = PreLieOperad(QQ)
                sage: x = PL.unlabelling(PL.one(1))
                sage: PL.group_product_without_numbers(x,x,2)
                B[[]]
            """
            return self.unlabelling(self.group_product_with_numbers(self.labelling(s), self.labelling(t), N))

        def inverse_with_numbers(self, s, N):
            r"""
            Return the inverse of `s` up to order `N`.

            This takes place in the group associated with the
            composition of the operad.

            INPUT:

            - `s` -- an element of the operad, in which the coefficient of
              the unit must be non-zero

            - `N` -- an integer, the truncation order

            EXAMPLES::

                sage: PL = PreLieOperad(QQ)
                sage: x = PL.one(1)
                sage: PL.inverse_with_numbers(x,1)
                B[1[]]

                sage: x = 4*PL.one(1)
                sage: PL.inverse_with_numbers(x,1)
                1/4*B[1[]]

                sage: x = PL.one(1)+(PL.one(1)<PL.one(2))
                sage: inv = PL.inverse_with_numbers(x,3); inv
                B[1[]] - B[2[1[]]]
                sage: PL.group_product_with_numbers(x,inv,3)
                B[1[]]

            TESTS::

                sage: x = (PL.one(1)<PL.one(2))
                sage: PL.inverse_with_numbers(x,3)
                Traceback (most recent call last):
                ...
                ValueError: this series is not invertible
            """
            if s.homogeneous_component(1) == 0:
                raise ValueError("this series is not invertible")
            coeff_o = s.homogeneous_component(1).coefficients()[0]
            dom = s.parent()
            t = s / coeff_o
            inv = dom.one(1) / coeff_o
            for j in range(2, N + 1):
                inv += dom.one(1) / coeff_o - dom.group_product_with_numbers(t, inv, j)
            return inv

        def right_division_with_numbers(self, s, t, N):
            r"""
            Compute the product of `s` by `t^{-1}` up to order `N`.

            INPUT:

            - `s` and `t` are elements of the operad, assumed
              to be labelled by consecutive integers.
            - `N` -- an integer, the truncation order

            EXAMPLES::

                sage: PL = PreLieOperad(QQ)
                sage: x = PL.one(1)
                sage: PL.right_division_with_numbers(x, x, 1)
                B[1[]]
            """
            dom = s.parent()
            inv = dom.one(1)
            for j in range(2, N + 1):
                inv += s.truncate(j) - dom.group_product_with_numbers(inv, t, j)
            return inv

    class ElementMethods:

        def eval(self, *args, **kwds):
            """
            Evaluate the coefficients of the element ``self`` of some operad.

            EXAMPLES::

                sage: w = PreLieOperad(PolynomialRing(ZZ,'x')).one('A')
                sage: x = w.parent().base_ring().gen()
                sage: ww = (1+x)*w ; ww
                (x+1)*B[A[]]
                sage: ww.eval(x=1)
                2*B[A[]]
            """
            return self.map_coefficients(lambda g: g(*args, **kwds))

        def canonical_labelling(self):
            """
            Return the canonical labelling of ``self``

            EXAMPLES::

                sage: w = PreLieOperad(PolynomialRing(ZZ,'x')).one('A')
                sage: w.canonical_labelling()
                B[1[]]
            """
            Op = self.parent()
            return Op.labelling(Op.unlabelling(self))
