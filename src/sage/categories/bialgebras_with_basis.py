r"""
Bialgebras with basis
"""
# ****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#  Copyright (C) 2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.tensor import tensor


class BialgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    r"""
    The category of bialgebras with a distinguished basis.

    EXAMPLES::

        sage: C = BialgebrasWithBasis(QQ); C
        Category of bialgebras with basis over Rational Field

        sage: sorted(C.super_categories(), key=str)
        [Category of algebras with basis over Rational Field,
         Category of bialgebras over Rational Field,
         Category of coalgebras with basis over Rational Field]

    TESTS::

        sage: TestSuite(BialgebrasWithBasis(ZZ)).run()
    """
    class ParentMethods:

        def convolution_product(self, *maps):
            r"""
            Return the convolution product (a map) of the given maps.

            Let `A` and `B` be bialgebras over a commutative ring `R`.
            Given maps `f_i : A \to B` for `1 \leq i < n`, define the
            convolution product

            .. MATH::

                (f_1 * f_2 * \cdots * f_n) := \mu^{(n-1)} \circ (f_1 \otimes
                f_2 \otimes \cdots \otimes f_n) \circ \Delta^{(n-1)},

            where `\Delta^{(k)} := \bigl(\Delta \otimes
            \mathrm{Id}^{\otimes(k-1)}\bigr) \circ \Delta^{(k-1)}`,
            with `\Delta^{(1)} = \Delta` (the ordinary coproduct in `A`) and
            `\Delta^{(0)} = \mathrm{Id}`; and with `\mu^{(k)} := \mu \circ
            \bigl(\mu^{(k-1)} \otimes \mathrm{Id})` and `\mu^{(1)} = \mu`
            (the ordinary product in `B`). See [Swe1969]_.

            (In the literature, one finds, e.g., `\Delta^{(2)}` for what we
            denote above as `\Delta^{(1)}`. See [KMN2012]_.)

            INPUT:

            - ``maps`` -- any number `n \geq 0` of linear maps `f_1, f_2,
              \ldots, f_n` on ``self``; or a single ``list`` or ``tuple``
              of such maps

            OUTPUT:

            - the new map `f_1 * f_2 * \cdots * f_2` representing their
              convolution product

            .. SEEALSO::

                :meth:`sage.categories.bialgebras.ElementMethods.convolution_product`

            AUTHORS:

            - Aaron Lauve - 12 June 2015 - Sage Days 65

            .. TODO::

                Remove dependency on ``modules_with_basis`` methods.

            EXAMPLES:

            We construct some maps: the identity, the antipode and
            projection onto the homogeneous component of degree 2::

                sage: Id = lambda x: x
                sage: Antipode = lambda x: x.antipode()
                sage: Proj2 = lambda x: x.parent().sum_of_terms([(m, c) for (m, c) in x if m.size() == 2])

            Compute the convolution product of the identity with itself and
            with the projection ``Proj2`` on the Hopf algebra of
            non-commutative symmetric functions::

                sage: R = NonCommutativeSymmetricFunctions(QQ).ribbon()
                sage: T = R.convolution_product([Id, Id])
                sage: [T(R(comp)) for comp in Compositions(3)]
                [4*R[1, 1, 1] + R[1, 2] + R[2, 1],
                 2*R[1, 1, 1] + 4*R[1, 2] + 2*R[2, 1] + 2*R[3],
                 2*R[1, 1, 1] + 2*R[1, 2] + 4*R[2, 1] + 2*R[3],
                 R[1, 2] + R[2, 1] + 4*R[3]]
                sage: T = R.convolution_product(Proj2, Id)
                sage: [T(R([i])) for i in range(1, 5)]
                [0, R[2], R[2, 1] + R[3], R[2, 2] + R[4]]

            Compute the convolution product of no maps on the Hopf algebra of
            symmetric functions in non-commuting variables. This is the
            composition of the counit with the unit::

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).m()
                sage: T = m.convolution_product()
                sage: [T(m(lam)) for lam in SetPartitions(0).list() + SetPartitions(2).list()]
                [m{}, 0, 0]

            Compute the convolution product of the projection ``Proj2`` with
            the identity on the Hopf algebra of symmetric functions in
            non-commuting variables::

                sage: T = m.convolution_product(Proj2, Id)
                sage: [T(m(lam)) for lam in SetPartitions(3)]
                [0,
                 m{{1, 2}, {3}} + m{{1, 2, 3}},
                 m{{1, 2}, {3}} + m{{1, 2, 3}},
                 m{{1, 2}, {3}} + m{{1, 2, 3}},
                 3*m{{1}, {2}, {3}} + 3*m{{1}, {2, 3}} + 3*m{{1, 3}, {2}}]

            Compute the convolution product of the antipode with itself and the
            identity map on group algebra of the symmetric group::

                sage: G = SymmetricGroup(3)
                sage: QG = GroupAlgebra(G, QQ)
                sage: x = QG.sum_of_terms([(p,p.number_of_peaks() + p.number_of_inversions()) for p in Permutations(3)]); x
                2*[1, 3, 2] + [2, 1, 3] + 3*[2, 3, 1] + 2*[3, 1, 2] + 3*[3, 2, 1]
                sage: T = QG.convolution_product(Antipode, Antipode, Id)
                sage: T(x)
                2*[1, 3, 2] + [2, 1, 3] + 2*[2, 3, 1] + 3*[3, 1, 2] + 3*[3, 2, 1]
            """
            onbasis = lambda x: self.term(x).convolution_product(*maps)
            return self.module_morphism(on_basis=onbasis, codomain=self)

    class ElementMethods:

        def adams_operator(self, n):
            r"""
            Compute the `n`-th convolution power of the identity morphism
            `\mathrm{Id}` on ``self``.

            INPUT:

            - ``n`` -- a nonnegative integer

            OUTPUT:

            - the image of ``self`` under the convolution power `\mathrm{Id}^{*n}`

            .. NOTE::

                In the literature, this is also called a Hopf power or
                Sweedler power, cf. [AL2015]_.

            .. SEEALSO::

                :meth:`sage.categories.bialgebras.ElementMethods.convolution_product`

            .. TODO::

                Remove dependency on ``modules_with_basis`` methods.

            EXAMPLES::

                sage: h = SymmetricFunctions(QQ).h()
                sage: h[5].adams_operator(2)
                2*h[3, 2] + 2*h[4, 1] + 2*h[5]
                sage: h[5].plethysm(2*h[1])
                2*h[3, 2] + 2*h[4, 1] + 2*h[5]
                sage: h([]).adams_operator(0)
                h[]
                sage: h([]).adams_operator(1)
                h[]
                sage: h[3,2].adams_operator(0)
                0
                sage: h[3,2].adams_operator(1)
                h[3, 2]

            ::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: S[4].adams_operator(5)
                5*S[1, 1, 1, 1] + 10*S[1, 1, 2] + 10*S[1, 2, 1] + 10*S[1, 3] + 10*S[2, 1, 1] + 10*S[2, 2] + 10*S[3, 1] + 5*S[4]


            ::

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).m()
                sage: m[[1,3],[2]].adams_operator(-2)
                3*m{{1}, {2, 3}} + 3*m{{1, 2}, {3}} + 6*m{{1, 2, 3}} - 2*m{{1, 3}, {2}}
            """
            if n < 0:
                if hasattr(self, 'antipode'):
                    T = lambda x: x.antipode()
                    n = abs(n)
                else:
                    raise ValueError("antipode not defined; cannot take negative convolution powers: {} < 0".format(n))
            else:
                T = lambda x: x
            return self.convolution_product([T] * n)

        def convolution_product(self, *maps):
            r"""
            Return the image of ``self`` under the convolution product (map) of
            the maps.

            Let `A` and `B` be bialgebras over a commutative ring `R`.
            Given maps `f_i : A \to B` for `1 \leq i < n`, define the
            convolution product

            .. MATH::

                (f_1 * f_2 * \cdots * f_n) := \mu^{(n-1)} \circ (f_1 \otimes
                f_2 \otimes \cdots \otimes f_n) \circ \Delta^{(n-1)},

            where `\Delta^{(k)} := \bigl(\Delta \otimes
            \mathrm{Id}^{\otimes(k-1)}\bigr) \circ \Delta^{(k-1)}`,
            with `\Delta^{(1)} = \Delta` (the ordinary coproduct in `A`) and
            `\Delta^{(0)} = \mathrm{Id}`; and with `\mu^{(k)} := \mu \circ
            \bigl(\mu^{(k-1)} \otimes \mathrm{Id})` and `\mu^{(1)} = \mu`
            (the ordinary product in `B`). See [Swe1969]_.

            (In the literature, one finds, e.g., `\Delta^{(2)}` for what we
            denote above as `\Delta^{(1)}`. See [KMN2012]_.)

            INPUT:

            - ``maps`` -- any number `n \geq 0` of linear maps `f_1, f_2,
              \ldots, f_n` on ``self.parent()``; or a single ``list`` or
              ``tuple`` of such maps

            OUTPUT:

            - the convolution product of ``maps`` applied to ``self``

            AUTHORS:

            - Amy Pang - 12 June 2015 - Sage Days 65

            .. TODO::

                Remove dependency on ``modules_with_basis`` methods.

            EXAMPLES:

            We compute convolution products of the identity and antipode maps
            on Schur functions::

                sage: Id = lambda x: x
                sage: Antipode = lambda x: x.antipode()
                sage: s = SymmetricFunctions(QQ).schur()
                sage: s[3].convolution_product(Id, Id)
                2*s[2, 1] + 4*s[3]
                sage: s[3,2].convolution_product(Id) == s[3,2]
                True

            The method accepts multiple arguments, or a single argument
            consisting of a list of maps::

                sage: s[3,2].convolution_product(Id, Id)
                2*s[2, 1, 1, 1] + 6*s[2, 2, 1] + 6*s[3, 1, 1] + 12*s[3, 2] + 6*s[4, 1] + 2*s[5]
                sage: s[3,2].convolution_product([Id, Id])
                2*s[2, 1, 1, 1] + 6*s[2, 2, 1] + 6*s[3, 1, 1] + 12*s[3, 2] + 6*s[4, 1] + 2*s[5]

            We test the defining property of the antipode morphism; namely,
            that the antipode is the inverse of the identity map in the
            convolution algebra whose identity element is the composition of
            the counit and unit::

                sage: s[3,2].convolution_product() == s[3,2].convolution_product(Antipode, Id) == s[3,2].convolution_product(Id, Antipode)
                True

            ::

                sage: Psi = NonCommutativeSymmetricFunctions(QQ).Psi()
                sage: Psi[2,1].convolution_product(Id, Id, Id)
                3*Psi[1, 2] + 6*Psi[2, 1]
                sage: (Psi[5,1] - Psi[1,5]).convolution_product(Id, Id, Id)
                -3*Psi[1, 5] + 3*Psi[5, 1]

            ::

                sage: G = SymmetricGroup(3)
                sage: QG = GroupAlgebra(G,QQ)
                sage: x = QG.sum_of_terms([(p,p.length()) for p in Permutations(3)]); x
                [1, 3, 2] + [2, 1, 3] + 2*[2, 3, 1] + 2*[3, 1, 2] + 3*[3, 2, 1]
                sage: x.convolution_product(Id, Id)
                5*[1, 2, 3] + 2*[2, 3, 1] + 2*[3, 1, 2]
                sage: x.convolution_product(Id, Id, Id)
                4*[1, 2, 3] + [1, 3, 2] + [2, 1, 3] + 3*[3, 2, 1]
                sage: x.convolution_product([Id]*6)
                9*[1, 2, 3]

            TESTS::

                sage: Id = lambda x: x
                sage: Antipode = lambda x: x.antipode()

            ::

                sage: h = SymmetricFunctions(QQ).h()
                sage: h[5].convolution_product([Id, Id])
                2*h[3, 2] + 2*h[4, 1] + 2*h[5]
                sage: h.one().convolution_product([Id, Antipode])
                h[]
                sage: h[3,2].convolution_product([Id, Antipode])
                0
                sage: h.one().convolution_product([Id, Antipode]) == h.one().convolution_product()
                True

            ::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: S[4].convolution_product([Id]*5)
                5*S[1, 1, 1, 1] + 10*S[1, 1, 2] + 10*S[1, 2, 1] + 10*S[1, 3]
                 + 10*S[2, 1, 1] + 10*S[2, 2] + 10*S[3, 1] + 5*S[4]

            ::

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).m()
                sage: m[[1,3],[2]].convolution_product([Antipode, Antipode])
                3*m{{1}, {2, 3}} + 3*m{{1, 2}, {3}} + 6*m{{1, 2, 3}} - 2*m{{1, 3}, {2}}
                sage: m[[]].convolution_product([])
                m{}
                sage: m[[1,3],[2]].convolution_product([])
                0

            ::

                sage: QS = SymmetricGroupAlgebra(QQ, 5)
                sage: x = QS.sum_of_terms(zip(Permutations(5)[3:6],[1,2,3])); x
                [1, 2, 4, 5, 3] + 2*[1, 2, 5, 3, 4] + 3*[1, 2, 5, 4, 3]
                sage: x.convolution_product([Antipode, Id])
                6*[1, 2, 3, 4, 5]
                sage: x.convolution_product(Id, Antipode, Antipode, Antipode)
                3*[1, 2, 3, 4, 5] + [1, 2, 4, 5, 3] + 2*[1, 2, 5, 3, 4]

            ::

                sage: G = SymmetricGroup(3)
                sage: QG = GroupAlgebra(G,QQ)
                sage: x = QG.sum_of_terms([(p,p.length()) for p in Permutations(3)]); x
                [1, 3, 2] + [2, 1, 3] + 2*[2, 3, 1] + 2*[3, 1, 2] + 3*[3, 2, 1]
                sage: x.convolution_product(Antipode, Id)
                9*[1, 2, 3]
                sage: x.convolution_product([Id, Antipode, Antipode, Antipode])
                5*[1, 2, 3] + 2*[2, 3, 1] + 2*[3, 1, 2]

            ::

                sage: s[3,2].counit().parent() == s[3,2].convolution_product().parent()
                False
            """
            # Be flexible on how the maps are entered: accept a list/tuple of
            # maps as well as multiple arguments
            if len(maps) == 1 and isinstance(maps[0], (list, tuple)):
                T = tuple(maps[0])
            else:
                T = maps

            H = self.parent()

            n = len(T)
            if n == 0:
                return H.one() * self.counit()
            if n == 1:
                return T[0](self)

            # We apply the maps T_i and products concurrently with coproducts, as this
            # seems to be faster than applying a composition of maps, e.g., (H.nfold_product) * tensor(T) * (H.nfold_coproduct).

            out = tensor((H.one(), self))
            HH = tensor((H, H))

            for mor in T[:-1]:
                # ALGORITHM:
                # `split_convolve` moves terms of the form x # y to x*Ti(y1) # y2 in Sweedler notation.
                def split_convolve(x_y):
                    x, y = x_y
                    return (((xy1, y2), c * d)
                            for ((y1, y2), d) in H.term(y).coproduct()
                            for (xy1, c) in H.term(x) * mor(H.term(y1)))
                out = HH.module_morphism(on_basis=lambda t: HH.sum_of_terms(split_convolve(t)), codomain=HH)(out)

            # Apply final map `T_n` to last term, `y`, and multiply.
            return HH.module_morphism(on_basis=lambda xy: H.term(xy[0]) * T[-1](H.term(xy[1])), codomain=H)(out)
