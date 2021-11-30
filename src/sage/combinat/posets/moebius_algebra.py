# -*- coding: utf-8 -*-
r"""
Möbius Algebras
"""
# ****************************************************************************
#       Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>,
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.bindable_class import BindableClass
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.algebras import Algebras
from sage.categories.realizations import Realizations, Category_realization_of_parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.free_module import CombinatorialFreeModule
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.integer_ring import ZZ


class BasisAbstract(CombinatorialFreeModule, BindableClass):
    """
    Abstract base class for a basis.
    """
    def __getitem__(self, x):
        """
        Return the basis element indexed by ``x``.

        INPUT:

        - ``x`` -- an element of the lattice

        EXAMPLES::

            sage: L = posets.BooleanLattice(4)
            sage: E = L.moebius_algebra(QQ).E()
            sage: E[5]
            E[5]
            sage: C = L.quantum_moebius_algebra().C()
            sage: C[5]
            C[5]
        """
        L = self.realization_of()._lattice
        return self.monomial(L(x))


class MoebiusAlgebra(Parent, UniqueRepresentation):
    r"""
    The Möbius algebra of a lattice.

    Let `L` be a lattice. The *Möbius algebra* `M_L` was originally
    constructed by Solomon [Solomon67]_ and has a natural basis
    `\{ E_x \mid x \in L \}` with multiplication given by
    `E_x \cdot E_y = E_{x \vee y}`. Moreover this has a basis given by
    orthogonal idempotents `\{ I_x \mid x \in L \}` (so
    `I_x I_y = \delta_{xy} I_x` where `\delta` is the Kronecker delta)
    related to the natural basis by

    .. MATH::

        I_x = \sum_{x \leq y} \mu_L(x, y) E_y,

    where `\mu_L` is the Möbius function of `L`.

    .. NOTE::

        We use the join `\vee` for our multiplication, whereas [Greene73]_
        and [Etienne98]_ define the Möbius algebra using the meet `\wedge`.
        This is done for compatibility with :class:`QuantumMoebiusAlgebra`.

    REFERENCES:

    .. [Solomon67] Louis Solomon.
       *The Burnside Algebra of a Finite Group*.
       Journal of Combinatorial Theory, **2**, 1967.
       :doi:`10.1016/S0021-9800(67)80064-4`.

    .. [Greene73] Curtis Greene.
       *On the Möbius algebra of a partially ordered set*.
       Advances in Mathematics, **10**, 1973.
       :doi:`10.1016/0001-8708(73)90106-0`.

    .. [Etienne98] Gwihen Etienne.
       *On the Möbius algebra of geometric lattices*.
       European Journal of Combinatorics, **19**, 1998.
       :doi:`10.1006/eujc.1998.0227`.
    """
    def __init__(self, R, L):
        """
        Initialize ``self``.

        TESTS::

            sage: L = posets.BooleanLattice(4)
            sage: M = L.moebius_algebra(QQ)
            sage: TestSuite(M).run()
        """
        cat = Algebras(R).Commutative().WithBasis()
        if L in FiniteEnumeratedSets():
            cat = cat.FiniteDimensional()
        self._lattice = L
        self._category = cat
        Parent.__init__(self, base=R, category=self._category.WithRealizations())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L = posets.BooleanLattice(4)
            sage: L.moebius_algebra(QQ)
            Moebius algebra of Finite lattice containing 16 elements over Rational Field
        """
        return "Moebius algebra of {} over {}".format(self._lattice, self.base_ring())

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the `B`-basis).

        EXAMPLES::

            sage: L = posets.BooleanLattice(4)
            sage: M = L.moebius_algebra(QQ)
            sage: M.a_realization()
            Moebius algebra of Finite lattice containing 16 elements
             over Rational Field in the natural basis
        """
        return self.E()

    def lattice(self):
        """
        Return the defining lattice of ``self``.

        EXAMPLES::

            sage: L = posets.BooleanLattice(4)
            sage: M = L.moebius_algebra(QQ)
            sage: M.lattice()
            Finite lattice containing 16 elements
            sage: M.lattice() == L
            True
        """
        return self._lattice

    class E(BasisAbstract):
        r"""
        The natural basis of a Möbius algebra.

        Let `E_x` and `E_y` be basis elements of `M_L` for some lattice `L`.
        Multiplication is given by `E_x E_y = E_{x \vee y}`.
        """
        def __init__(self, M, prefix='E'):
            """
            Initialize ``self``.

            TESTS::

                sage: L = posets.BooleanLattice(4)
                sage: M = L.moebius_algebra(QQ)
                sage: TestSuite(M.E()).run()
            """
            self._basis_name = "natural"
            CombinatorialFreeModule.__init__(self, M.base_ring(),
                                             tuple(M._lattice),
                                             prefix=prefix,
                                             category=MoebiusAlgebraBases(M))

        @cached_method
        def _to_idempotent_basis(self, x):
            """
            Convert the element indexed by ``x`` to the idempotent basis.

            EXAMPLES::

                sage: M = posets.BooleanLattice(4).moebius_algebra(QQ)
                sage: E = M.E()
                sage: all(E(E._to_idempotent_basis(x)) == E.monomial(x)
                ....:     for x in E.basis().keys())
                True
            """
            M = self.realization_of()
            I = M.idempotent()
            return I.sum_of_monomials(M._lattice.order_filter([x]))

        def product_on_basis(self, x, y):
            """
            Return the product of basis elements indexed by ``x`` and ``y``.

            EXAMPLES::

                sage: L = posets.BooleanLattice(4)
                sage: E = L.moebius_algebra(QQ).E()
                sage: E.product_on_basis(5, 14)
                E[15]
                sage: E.product_on_basis(2, 8)
                E[10]

            TESTS::

                sage: M = posets.BooleanLattice(4).moebius_algebra(QQ)
                sage: E = M.E()
                sage: I = M.I()
                sage: all(I(x)*I(y) == I(x*y) for x in E.basis() for y in E.basis())
                True
            """
            return self.monomial(self.realization_of()._lattice.join(x, y))

        @cached_method
        def one(self):
            """
            Return the element ``1`` of ``self``.

            EXAMPLES::

                sage: L = posets.BooleanLattice(4)
                sage: E = L.moebius_algebra(QQ).E()
                sage: E.one()
                E[0]
            """
            elts = self.realization_of()._lattice.minimal_elements()
            return self.sum_of_monomials(elts)

    natural = E

    class I(BasisAbstract):
        r"""
        The (orthogonal) idempotent basis of a Möbius algebra.

        Let `I_x` and `I_y` be basis elements of `M_L` for some lattice `L`.
        Multiplication is given by `I_x I_y = \delta_{xy} I_x` where
        `\delta_{xy}` is the Kronecker delta.
        """
        def __init__(self, M, prefix='I'):
            """
            Initialize ``self``.

            TESTS::

                sage: L = posets.BooleanLattice(4)
                sage: M = L.moebius_algebra(QQ)
                sage: TestSuite(M.I()).run()

            Check that the transition maps can be pickled::

                sage: L = posets.BooleanLattice(4)
                sage: M = L.moebius_algebra(QQ)
                sage: E = M.E()
                sage: I = M.I()
                sage: phi = E.coerce_map_from(I)
                sage: loads(dumps(phi))
                Generic morphism:
                ...
            """
            self._basis_name = "idempotent"
            CombinatorialFreeModule.__init__(self, M.base_ring(),
                                             tuple(M._lattice),
                                             prefix=prefix,
                                             category=MoebiusAlgebraBases(M))

            # Change of basis:
            E = M.E()
            self.module_morphism(self._to_natural_basis,
                                 codomain=E, category=self.category(),
                                 triangular='lower', unitriangular=True,
                                 key=M._lattice._element_to_vertex
                                 ).register_as_coercion()

            E.module_morphism(E._to_idempotent_basis,
                              codomain=self, category=self.category(),
                              triangular='lower', unitriangular=True,
                              key=M._lattice._element_to_vertex
                              ).register_as_coercion()

        @cached_method
        def _to_natural_basis(self, x):
            """
            Convert the element indexed by ``x`` to the natural basis.

            EXAMPLES::

                sage: M = posets.BooleanLattice(4).moebius_algebra(QQ)
                sage: I = M.I()
                sage: all(I(I._to_natural_basis(x)) == I.monomial(x)
                ....:     for x in I.basis().keys())
                True
            """
            M = self.realization_of()
            N = M.natural()
            moebius = M._lattice.moebius_function
            return N.sum_of_terms((y, moebius(x, y))
                                  for y in M._lattice.order_filter([x]))

        def product_on_basis(self, x, y):
            """
            Return the product of basis elements indexed by ``x`` and ``y``.

            EXAMPLES::

                sage: L = posets.BooleanLattice(4)
                sage: I = L.moebius_algebra(QQ).I()
                sage: I.product_on_basis(5, 14)
                0
                sage: I.product_on_basis(2, 2)
                I[2]

            TESTS::

                sage: M = posets.BooleanLattice(4).moebius_algebra(QQ)
                sage: E = M.E()
                sage: I = M.I()
                sage: all(E(x)*E(y) == E(x*y) for x in I.basis() for y in I.basis())
                True
            """
            if x == y:
                return self.monomial(x)
            return self.zero()

        @cached_method
        def one(self):
            """
            Return the element ``1`` of ``self``.

            EXAMPLES::

                sage: L = posets.BooleanLattice(4)
                sage: I = L.moebius_algebra(QQ).I()
                sage: I.one()
                I[0] + I[1] + I[2] + I[3] + I[4] + I[5] + I[6] + I[7] + I[8]
                 + I[9] + I[10] + I[11] + I[12] + I[13] + I[14] + I[15]
            """
            return self.sum_of_monomials(self.realization_of()._lattice)

        def __getitem__(self, x):
            """
            Return the basis element indexed by ``x``.

            INPUT:

            - ``x`` -- an element of the lattice

            EXAMPLES::

                sage: L = posets.BooleanLattice(4)
                sage: I = L.moebius_algebra(QQ).I()
                sage: I[5]
                I[5]
            """
            L = self.realization_of()._lattice
            return self.monomial(L(x))

    idempotent = I


class QuantumMoebiusAlgebra(Parent, UniqueRepresentation):
    r"""
    The quantum Möbius algebra of a lattice.

    Let `L` be a lattice, and we define the *quantum Möbius algebra* `M_L(q)`
    as the algebra with basis `\{ E_x \mid x \in L \}` with
    multiplication given by

    .. MATH::

        E_x E_y = \sum_{z \geq a \geq x \vee y} \mu_L(a, z)
        q^{\operatorname{crk} a} E_z,

    where `\mu_L` is the Möbius function of `L` and `\operatorname{crk}`
    is the corank function (i.e., `\operatorname{crk} a =
    \operatorname{rank} L - \operatorname{rank}` a). At `q = 1`, this
    reduces to the multiplication formula originally given by Solomon.
    """
    def __init__(self, L, q=None):
        """
        Initialize ``self``.

        TESTS::

            sage: L = posets.BooleanLattice(4)
            sage: M = L.quantum_moebius_algebra()
            sage: TestSuite(M).run() # long time

            sage: from sage.combinat.posets.moebius_algebra import QuantumMoebiusAlgebra
            sage: L = posets.Crown(2)
            sage: QuantumMoebiusAlgebra(L)
            Traceback (most recent call last):
            ...
            ValueError: L must be a lattice
        """
        if not L.is_lattice():
            raise ValueError("L must be a lattice")
        if q is None:
            q = LaurentPolynomialRing(ZZ, 'q').gen()
        self._q = q
        R = q.parent()
        cat = Algebras(R).WithBasis()
        if L in FiniteEnumeratedSets():
            cat = cat.Commutative().FiniteDimensional()
        self._lattice = L
        self._category = cat
        Parent.__init__(self, base=R, category=self._category.WithRealizations())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L = posets.BooleanLattice(4)
            sage: L.quantum_moebius_algebra()
            Quantum Moebius algebra of Finite lattice containing 16 elements
             with q=q over Univariate Laurent Polynomial Ring in q over Integer Ring
        """
        txt = "Quantum Moebius algebra of {} with q={} over {}"
        return txt.format(self._lattice, self._q, self.base_ring())

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the `B`-basis).

        EXAMPLES::

            sage: L = posets.BooleanLattice(4)
            sage: M = L.quantum_moebius_algebra()
            sage: M.a_realization()
            Quantum Moebius algebra of Finite lattice containing 16 elements
             with q=q over Univariate Laurent Polynomial Ring in q
             over Integer Ring in the natural basis
        """
        return self.E()

    def lattice(self):
        """
        Return the defining lattice of ``self``.

        EXAMPLES::

            sage: L = posets.BooleanLattice(4)
            sage: M = L.quantum_moebius_algebra()
            sage: M.lattice()
            Finite lattice containing 16 elements
            sage: M.lattice() == L
            True
        """
        return self._lattice

    class E(BasisAbstract):
        r"""
        The natural basis of a quantum Möbius algebra.

        Let `E_x` and `E_y` be basis elements of `M_L` for some lattice `L`.
        Multiplication is given by

        .. MATH::

            E_x E_y = \sum_{z \geq a \geq x \vee y} \mu_L(a, z)
            q^{\operatorname{crk} a} E_z,

        where `\mu_L` is the Möbius function of `L` and `\operatorname{crk}`
        is the corank function (i.e., `\operatorname{crk} a =
        \operatorname{rank} L - \operatorname{rank}` a).
        """
        def __init__(self, M, prefix='E'):
            """
            Initialize ``self``.

            TESTS::

                sage: L = posets.BooleanLattice(4)
                sage: M = L.quantum_moebius_algebra()
                sage: TestSuite(M.E()).run() # long time
            """
            self._basis_name = "natural"
            CombinatorialFreeModule.__init__(self, M.base_ring(),
                                             tuple(M._lattice),
                                             prefix=prefix,
                                             category=MoebiusAlgebraBases(M))

        def product_on_basis(self, x, y):
            """
            Return the product of basis elements indexed by ``x`` and ``y``.

            EXAMPLES::

                sage: L = posets.BooleanLattice(4)
                sage: E = L.quantum_moebius_algebra().E()
                sage: E.product_on_basis(5, 14)
                E[15]
                sage: E.product_on_basis(2, 8)
                q^2*E[10] + (q-q^2)*E[11] + (q-q^2)*E[14] + (1-2*q+q^2)*E[15]
            """
            L = self.realization_of()._lattice
            q = self.realization_of()._q
            moebius = L.moebius_function
            rank = L.rank_function()
            R = L.rank()
            j = L.join(x, y)
            return self.sum_of_terms((z, moebius(a, z) * q**(R - rank(a)))
                                     for z in L.order_filter([j])
                                     for a in L.closed_interval(j, z))

        @cached_method
        def one(self):
            """
            Return the element ``1`` of ``self``.

            EXAMPLES::

                sage: L = posets.BooleanLattice(4)
                sage: E = L.quantum_moebius_algebra().E()
                sage: all(E.one() * b == b for b in E.basis())
                True
            """
            L = self.realization_of()._lattice
            q = self.realization_of()._q
            moebius = L.moebius_function
            rank = L.rank_function()
            R = L.rank()
            return self.sum_of_terms((x, moebius(y, x) * q**(rank(y) - R))
                                     for x in L for y in L.order_ideal([x]))

    natural = E

    class C(BasisAbstract):
        r"""
        The characteristic basis of a quantum Möbius algebra.

        The characteristic basis `\{ C_x \mid x \in L \}` of `M_L`
        for some lattice `L` is defined by

        .. MATH::

            C_x = \sum_{a \geq x} P(F^x; q) E_a,

        where `F^x = \{ y \in L \mid y \geq x \}` is the principal order
        filter of `x` and `P(F^x; q)` is the characteristic polynomial
        of the (sub)poset `F^x`.
        """
        def __init__(self, M, prefix='C'):
            """
            Initialize ``self``.

            TESTS::

                sage: L = posets.BooleanLattice(3)
                sage: M = L.quantum_moebius_algebra()
                sage: TestSuite(M.C()).run() # long time
            """
            self._basis_name = "characteristic"
            CombinatorialFreeModule.__init__(self, M.base_ring(),
                                             tuple(M._lattice),
                                             prefix=prefix,
                                             category=MoebiusAlgebraBases(M))

            # Change of basis:
            E = M.E()
            phi = self.module_morphism(self._to_natural_basis,
                                       codomain=E, category=self.category(),
                                       triangular='lower', unitriangular=True,
                                       key=M._lattice._element_to_vertex)

            phi.register_as_coercion()
            (~phi).register_as_coercion()

        @cached_method
        def _to_natural_basis(self, x):
            """
            Convert the element indexed by ``x`` to the natural basis.

            EXAMPLES::

                sage: M = posets.BooleanLattice(4).quantum_moebius_algebra()
                sage: C = M.C()
                sage: all(C(C._to_natural_basis(x)) == C.monomial(x)
                ....:     for x in C.basis().keys())
                True
            """
            M = self.realization_of()
            N = M.natural()
            q = M._q
            L = M._lattice

            def poly(x, y):
                return L.subposet(L.closed_interval(x, y)).characteristic_polynomial()
            return N.sum_of_terms((y, poly(x, y)(q=q))
                                  for y in L.order_filter([x]))

    characteristic_basis = C

    class KL(BasisAbstract):
        r"""
        The Kazhdan-Lusztig basis of a quantum Möbius algebra.

        The Kazhdan-Lusztig basis `\{ B_x \mid x \in L \}` of `M_L`
        for some lattice `L` is defined by

        .. MATH::

            B_x = \sum_{y \geq x} P_{x,y}(q) E_a,

        where `P_{x,y}(q)` is the Kazhdan-Lusztig polynomial of `L`,
        following the definition given in [EPW14]_.

        EXAMPLES:

        We construct some examples of Proposition 4.5 of [EPW14]_::

            sage: M = posets.BooleanLattice(4).quantum_moebius_algebra()
            sage: KL = M.KL()
            sage: KL[4] * KL[5]
            (q^2+q^3)*KL[5] + (q+2*q^2+q^3)*KL[7] + (q+2*q^2+q^3)*KL[13]
             + (1+3*q+3*q^2+q^3)*KL[15]
            sage: KL[4] * KL[15]
            (1+3*q+3*q^2+q^3)*KL[15]
            sage: KL[4] * KL[10]
            (q+3*q^2+3*q^3+q^4)*KL[14] + (1+4*q+6*q^2+4*q^3+q^4)*KL[15]
        """
        def __init__(self, M, prefix='KL'):
            """
            Initialize ``self``.

            TESTS::

                sage: L = posets.BooleanLattice(4)
                sage: M = L.quantum_moebius_algebra()
                sage: TestSuite(M.KL()).run() # long time
            """
            self._basis_name = "Kazhdan-Lusztig"
            CombinatorialFreeModule.__init__(self, M.base_ring(),
                                             tuple(M._lattice),
                                             prefix=prefix,
                                             category=MoebiusAlgebraBases(M))

            # Change of basis:
            E = M.E()
            phi = self.module_morphism(self._to_natural_basis,
                                       codomain=E, category=self.category(),
                                       triangular='lower', unitriangular=True,
                                       key=M._lattice._element_to_vertex)

            phi.register_as_coercion()
            (~phi).register_as_coercion()

        @cached_method
        def _to_natural_basis(self, x):
            """
            Convert the element indexed by ``x`` to the natural basis.

            EXAMPLES::

                sage: M = posets.BooleanLattice(4).quantum_moebius_algebra()
                sage: KL = M.KL()
                sage: all(KL(KL._to_natural_basis(x)) == KL.monomial(x) # long time
                ....:     for x in KL.basis().keys())
                True
            """
            M = self.realization_of()
            L = M._lattice
            E = M.E()
            q = M._q
            rank = L.rank_function()
            return E.sum_of_terms((y, q**(rank(y) - rank(x)) *
                                   L.kazhdan_lusztig_polynomial(x, y)(q=q**-2))
                                  for y in L.order_filter([x]))

    kazhdan_lusztig = KL


class MoebiusAlgebraBases(Category_realization_of_parent):
    r"""
    The category of bases of a Möbius algebra.

    INPUT:

    - ``base`` -- a Möbius algebra

    TESTS::

        sage: from sage.combinat.posets.moebius_algebra import MoebiusAlgebraBases
        sage: M = posets.BooleanLattice(4).moebius_algebra(QQ)
        sage: bases = MoebiusAlgebraBases(M)
        sage: M.E() in bases
        True
    """
    def _repr_(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.posets.moebius_algebra import MoebiusAlgebraBases
            sage: M = posets.BooleanLattice(4).moebius_algebra(QQ)
            sage: MoebiusAlgebraBases(M)
            Category of bases of Moebius algebra of Finite lattice
             containing 16 elements over Rational Field
        """
        return "Category of bases of {}".format(self.base())

    def super_categories(self):
        r"""
        The super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.posets.moebius_algebra import MoebiusAlgebraBases
            sage: M = posets.BooleanLattice(4).moebius_algebra(QQ)
            sage: bases = MoebiusAlgebraBases(M)
            sage: bases.super_categories()
            [Category of finite dimensional commutative algebras with basis over Rational Field,
             Category of realizations of Moebius algebra of Finite lattice
                containing 16 elements over Rational Field]
        """
        return [self.base()._category, Realizations(self.base())]

    class ParentMethods:
        def _repr_(self):
            """
            Text representation of this basis of a Möbius algebra.

            EXAMPLES::

                sage: M = posets.BooleanLattice(4).moebius_algebra(QQ)
                sage: M.E()
                Moebius algebra of Finite lattice containing 16 elements
                 over Rational Field in the natural basis
                sage: M.I()
                Moebius algebra of Finite lattice containing 16 elements
                 over Rational Field in the idempotent basis
            """
            return "{} in the {} basis".format(self.realization_of(), self._basis_name)

        def product_on_basis(self, x, y):
            """
            Return the product of basis elements indexed by ``x`` and ``y``.

            EXAMPLES::

                sage: L = posets.BooleanLattice(4)
                sage: C = L.quantum_moebius_algebra().C()
                sage: C.product_on_basis(5, 14)
                q^3*C[15]
                sage: C.product_on_basis(2, 8)
                q^4*C[10]
            """
            R = self.realization_of().a_realization()
            return self(R(self.monomial(x)) * R(self.monomial(y)))

        @cached_method
        def one(self):
            """
            Return the element ``1`` of ``self``.

            EXAMPLES::

                sage: L = posets.BooleanLattice(4)
                sage: C = L.quantum_moebius_algebra().C()
                sage: all(C.one() * b == b for b in C.basis())
                True
            """
            R = self.realization_of().a_realization()
            return self(R.one())

    class ElementMethods:
        pass
