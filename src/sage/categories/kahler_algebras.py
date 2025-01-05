r"""
Kähler Algebras

AUTHORS:

- Shriya M
"""

from sage.categories.category_types import Category_over_base_ring
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.finite_dimensional_algebras_with_basis import FiniteDimensionalAlgebrasWithBasis
from sage.categories.filtered_modules_with_basis import FilteredModulesWithBasis
from sage.misc.abstract_method import abstract_method
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.misc.cachefunc import cached_method


# ****************************************************************************
#       Copyright (C) 2024 Shriya M <25shriya at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
class KahlerAlgebras(Category_over_base_ring):
    r"""
    The category of graded algebras satisfying the Kähler package.
    A finite-dimensional graded algebra `\bigoplus_{k=1}^{r}A^k` satisfies
    the *Kähler package* if the following properties hold:

    - Poincaré duality: There exists a perfect `\mathbb{Z}`-bilinear pairing
      given by

      .. MATH::

          A^k \times A^{r-k} \longrightarrow \mathbb{Z} \\
            (a,b) \mapsto \text{deg}(a \cdot b).

    - Hard-Lefschetz Theorem: The graded algebra contains *Lefschetz elements*
      `\omega \in A^{1}_{\mathbb{R}}` such that multiplication by `\omega` is
      an injection from `A^k_{\mathbb{R}} \longrightarrow A^{k+1}_{\mathbb{R}}`
      for all `k < \frac{r}{2}`.

    - Hodge-Riemann-Minikowski Relations: Every Lefchetz element `\omega`,
      define quadratic forms on `A^{k}_{\mathbb{R}}` given by

      .. MATH::

          a \mapsto (-1)^k \text{deg}(a \cdot \omega^{r-2k} \cdot a)

      This quadratic form becomes positive definite upon restriction to the
      kernel of the following map

      .. MATH::

          A^k_\mathbb{R} \longrightarrow A^{r-k+1}_\mathbb{R} \\
                a \mapsto a \cdot \omega^{r-2k+1}.

    REFERENCES:

    - [ANR2023]_

    EXAMPLES::

        sage: from sage.categories.kahler_algebras import KahlerAlgebras

        sage: C = KahlerAlgebras(QQ); C
        Category of kahler algebras over Rational Field
        sage: sorted(C.super_categories(), key=str)
        [Category of finite dimensional algebras with basis over Rational Field]

    TESTS::

        sage: C = KahlerAlgebras(QQ)
        sage: TestSuite(C).run()
    """
    def super_categories(self):
        return [GradedAlgebrasWithBasis(self.base_ring()).FiniteDimensional()]

    class ParentMethods:
        @abstract_method
        def poincare_pairing(a,b):
            r"""
            Return the Poincaré pairing of two elements of the Kähler algebra.

                EXAMPLES::

                    sage: ch = matroids.catalog.Fano().chow_ring(QQ, True, 'fy')
                    sage: Ba, Bb, Bc, Bd, Be, Bf, Bg, Babf, Bace, Badg, Bbcd, Bbeg, Bcfg, Bdef, Babcdefg = ch.gens()[8:]
                    sage: u = ch(-Babf^2 + Bcfg^2 - 8/7*Bc*Babcdefg + 1/2*Bd*Babcdefg - Bf*Babcdefg - Bg*Babcdefg); u
                    -Babf^2 + Bcfg^2 - 8/7*Bc*Babcdefg + 1/2*Bd*Babcdefg - Bf*Babcdefg - Bg*Babcdefg
                    sage: v = ch(Bg - 2/37*Babf + Badg + Bbeg + Bcfg + Babcdefg); v
                    Bg - 2/37*Babf + Badg + Bbeg + Bcfg + Babcdefg
                    sage: ch.poincare_pairing(v, u)
                    3
            """
            pass

        @abstract_method
        def lefschetz_element():
            pass

        def hodge_riemann_relations(self, k):
            r"""
            Return the quadratic form for the corresponding k (< r/2) for the
            Kähler algebra.

            EXAMPLES::

                sage: ch = matroids.Uniform(4,6).chow_ring(QQ, False)
                sage: ch.hodge_riemann_relations(1)
                Quadratic form in 36 variables over Rational Field with coefficients:
                [ 3 -1 -1 3 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 3 ]
                [ * 3 -1 3 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 3 ]
                [ * * 3 3 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 3 ]
                [ * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * 3 -1 3 -1 3 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 3 ]
                [ * * * * * 3 3 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 3 ]
                [ * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * 3 3 3 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 ]
                [ * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * 3 -1 3 -1 -1 3 -1 -1 3 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 ]
                [ * * * * * * * * * * * 3 3 -1 3 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 ]
                [ * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * 3 3 3 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 -1 -1 -1 -1 3 ]
                [ * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * 3 3 3 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 -1 -1 -1 3 ]
                [ * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * 3 -1 3 -1 3 -1 -1 -1 -1 3 -1 -1 -1 3 -1 3 ]
                [ * * * * * * * * * * * * * * * * * * * * * 3 3 -1 -1 3 -1 -1 3 -1 -1 -1 3 -1 -1 3 ]
                [ * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * 3 3 3 -1 3 -1 -1 -1 3 -1 -1 -1 3 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * 3 3 3 3 -1 -1 -1 -1 3 3 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 3 3 3 3 3 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 -1 ]
                [ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 3 ]
                sage: ch.hodge_riemann_relations(3)
                Traceback (most recent call last):
                ...
                ValueError: k must be less than r < 2
            """
            r = self.top_degree()
            if k > (r/2):
                raise ValueError("k must be less than r < 2")
            basis_k = []
            lefschetz_el = self.lefschetz_element()
            for b in self.basis():
                if b.homogeneous_degree() == k:
                    basis_k.append(b)
            coeff = []
            for i,el in enumerate(basis_k):
                for j in range(i, len(basis_k)):
                    coeff.append((el * (lefschetz_el ** (r-(2*k)) * basis_k[j])).degree())
            return QuadraticForm(self.base_ring(), len(basis_k), coeff)