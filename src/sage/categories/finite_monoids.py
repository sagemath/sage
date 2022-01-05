r"""
Finite monoids
"""
# ****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom


class FiniteMonoids(CategoryWithAxiom):
    """
    The category of finite (multiplicative) :class:`monoids <Monoids>`.

    A finite monoid is a :class:`finite sets <FiniteSets>` endowed
    with an associative unital binary operation `*`.

    EXAMPLES::

        sage: FiniteMonoids()
        Category of finite monoids
        sage: FiniteMonoids().super_categories()
        [Category of monoids, Category of finite semigroups]

    TESTS::

        sage: TestSuite(FiniteMonoids()).run()
    """
    class ParentMethods:

        def nerve(self):
            r"""
            The nerve (classifying space) of this monoid.

            OUTPUT: the nerve $BG$ (if $G$ denotes this monoid), as a
            simplicial set.  The $k$-dimensional simplices of this
            object are indexed by products of $k$ elements in the
            monoid:

            .. MATH::

                a_1 * a_2 * \cdots * a_k

            The 0th face of this is obtained by deleting $a_1$, and
            the $k$-th face is obtained by deleting $a_k$. The other
            faces are obtained by multiplying elements: the 1st face
            is

            .. MATH::

               (a1 * a_2) * \cdots * a_k

            and so on. See :wikipedia:`Nerve_(category_theory)`, which
            describes the construction of the nerve as a simplicial
            set.

            A simplex in this simplicial set will be degenerate if in
            the corresponding product of $k$ elements, one of those
            elements is the identity. So we only need to keep track of
            the products of non-identity elements. Similarly, if a
            product `a_{i-1} a_i` is the identity element, then the
            corresponding face of the simplex will be a degenerate
            simplex.

            EXAMPLES:

            The nerve (classifying space) of the cyclic group of order
            2 is infinite-dimensional real projective space. ::

                sage: Sigma2 = groups.permutation.Cyclic(2)
                sage: BSigma2 = Sigma2.nerve()
                sage: BSigma2.cohomology(4, base_ring=GF(2))
                Vector space of dimension 1 over Finite Field of size 2

            The `k`-simplices of the nerve are named after the chains
            of `k` non-unit elements to be multiplied. The group
            `\Sigma_2` has two elements, written ``()`` (the identity
            element) and ``(1,2)`` in Sage. So the 1-cells and 2-cells
            in `B\Sigma_2` are::

                sage: BSigma2.n_cells(1)
                [(1,2)]
                sage: BSigma2.n_cells(2)
                [(1,2) * (1,2)]

            Another construction of the group, with different names
            for its elements::

                sage: C2 = groups.misc.MultiplicativeAbelian([2])
                sage: BC2 = C2.nerve()
                sage: BC2.n_cells(0)
                [1]
                sage: BC2.n_cells(1)
                [f]
                sage: BC2.n_cells(2)
                [f * f]

            With mod `p` coefficients, `B \Sigma_p` should have its
            first nonvanishing homology group in dimension `p`::

                sage: Sigma3 = groups.permutation.Symmetric(3)
                sage: BSigma3 = Sigma3.nerve()
                sage: BSigma3.homology(range(4), base_ring=GF(3))
                {0: Vector space of dimension 0 over Finite Field of size 3,
                1: Vector space of dimension 0 over Finite Field of size 3,
                2: Vector space of dimension 0 over Finite Field of size 3,
                3: Vector space of dimension 1 over Finite Field of size 3}

            Note that we can construct the `n`-skeleton for
            `B\Sigma_2` for relatively large values of `n`, while for
            `B\Sigma_3`, the complexes get large pretty quickly::

                sage: Sigma2.nerve().n_skeleton(14)
                Simplicial set with 15 non-degenerate simplices

                sage: BSigma3 = Sigma3.nerve()
                sage: BSigma3.n_skeleton(3)
                Simplicial set with 156 non-degenerate simplices
                sage: BSigma3.n_skeleton(4)
                Simplicial set with 781 non-degenerate simplices

            Finally, note that the classifying space of the order `p`
            cyclic group is smaller than that of the symmetric group
            on `p` letters, and its first homology group appears
            earlier::

                sage: C3 = groups.misc.MultiplicativeAbelian([3])
                sage: list(C3)
                [1, f, f^2]
                sage: BC3 = C3.nerve()
                sage: BC3.n_cells(1)
                [f, f^2]
                sage: BC3.n_cells(2)
                [f * f, f * f^2, f^2 * f, f^2 * f^2]
                sage: len(BSigma3.n_cells(2))
                25
                sage: len(BC3.n_cells(3))
                8
                sage: len(BSigma3.n_cells(3))
                125

                sage: BC3.homology(range(5), base_ring=GF(3))
                {0: Vector space of dimension 0 over Finite Field of size 3,
                 1: Vector space of dimension 1 over Finite Field of size 3,
                 2: Vector space of dimension 1 over Finite Field of size 3,
                 3: Vector space of dimension 1 over Finite Field of size 3,
                 4: Vector space of dimension 1 over Finite Field of size 3}

                sage: BC5 = groups.permutation.Cyclic(5).nerve()
                sage: BC5.homology(range(5), base_ring=GF(5))
                {0: Vector space of dimension 0 over Finite Field of size 5,
                1: Vector space of dimension 1 over Finite Field of size 5,
                2: Vector space of dimension 1 over Finite Field of size 5,
                3: Vector space of dimension 1 over Finite Field of size 5,
                4: Vector space of dimension 1 over Finite Field of size 5}
            """
            from sage.topology.simplicial_set_examples import Nerve
            return Nerve(self)

        def rhodes_radical_congruence(self, base_ring=None):
            r"""
            Return the Rhodes radical congruence of the semigroup.

            The Rhodes radical congruence is the congruence induced on S by the
            map `S \rightarrow kS \rightarrow kS / rad kS` with k a field.

            INPUT:

            - ``base_ring`` (default: `\QQ`) a field

            OUTPUT:

            - A list of couples (m, n) with `m \neq n` in the lexicographic
              order for the enumeration of the monoid ``self``.

            EXAMPLES::

                sage: M = Monoids().Finite().example()
                sage: M.rhodes_radical_congruence()
                [(0, 6), (2, 8), (4, 10)]
                sage: from sage.monoids.hecke_monoid import HeckeMonoid
                sage: H3 = HeckeMonoid(SymmetricGroup(3))
                sage: H3.repr_element_method(style="reduced")
                sage: H3.rhodes_radical_congruence()
                [([1, 2], [2, 1]), ([1, 2], [1, 2, 1]), ([2, 1], [1, 2, 1])]

            By Maschke's theorem, every group algebra over `\QQ`
            is semisimple hence the Rhodes radical of a group must be trivial::

                sage: SymmetricGroup(3).rhodes_radical_congruence()
                []
                sage: DihedralGroup(10).rhodes_radical_congruence()
                []

            REFERENCES:

            - [Rho69]_
            """
            from sage.rings.rational_field import QQ
            if base_ring is None:
                base_ring = QQ
            kS = self.algebra(base_ring)
            kSrad = kS.radical()
            res = []
            for m in self:
                for n in self:
                    if (m == n) or ((n, m) in res):
                        continue
                    try:
                        kSrad.retract(kS(m) - kS(n))
                    except ValueError:
                        pass
                    else:
                        res.append((m, n))
            return res

    class ElementMethods:
        def pseudo_order(self):
            r"""
            Returns the pair `[k, j]` with `k` minimal and `0\leq j <k` such
            that ``self^k == self^j``.

            Note that `j` is uniquely determined.

            EXAMPLES::

                sage: M = FiniteMonoids().example(); M
                An example of a finite multiplicative monoid: the integers modulo 12

                sage: x = M(2)
                sage: [ x^i for i in range(7) ]
                [1, 2, 4, 8, 4, 8, 4]
                sage: x.pseudo_order()
                [4, 2]

                sage: x = M(3)
                sage: [ x^i for i in range(7) ]
                [1, 3, 9, 3, 9, 3, 9]
                sage: x.pseudo_order()
                [3, 1]

                sage: x = M(4)
                sage: [ x^i for i in range(7) ]
                [1, 4, 4, 4, 4, 4, 4]
                sage: x.pseudo_order()
                [2, 1]

                sage: x = M(5)
                sage: [ x^i for i in range(7) ]
                [1, 5, 1, 5, 1, 5, 1]
                sage: x.pseudo_order()
                [2, 0]

            TODO: more appropriate name? see, for example, Jean-Eric Pin's
            lecture notes on semigroups.
            """
            self_powers = {self.parent().one(): 0}
            k = 1
            self_power_k = self
            while self_power_k not in self_powers:
                self_powers[self_power_k] = k
                k += 1
                self_power_k = self_power_k * self
            return [k, self_powers[self_power_k]]
