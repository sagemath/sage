r"""
Kähler Algebras

AUTHORS:

- Shriya M
"""
# ****************************************************************************
#       Copyright (C) 2024 Shriya M <25shriya at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.category_types import Category_over_base_ring
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.finite_dimensional_algebras_with_basis import FiniteDimensionalAlgebrasWithBasis
from sage.categories.filtered_modules_with_basis import FilteredModulesWithBasis
from sage.misc.abstract_method import abstract_method
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.misc.cachefunc import cached_method


class KahlerAlgebras(Category_over_base_ring):
    r"""
    The category of graded algebras satisfying the Kähler package.
    A finite-dimensional graded algebra `\bigoplus_{k=1}^{r}A^k` satisfies
    the *Kähler package* if the following properties hold:

    - Poincaré duality: There exists a perfect `\ZZ`-bilinear pairing
      given by

      .. MATH::

          A^k \times A^{r-k} \longrightarrow \ZZ \\
            (a,b) \mapsto \deg(a \cdot b).

    - Hard-Lefschetz Theorem: The graded algebra contains *Lefschetz elements*
      `\omega \in A^{1}_{\RR}` such that multiplication by `\omega` is
      an injection from `A^k_{\RR} \longrightarrow A^{k+1}_{\RR}`
      for all `k < \frac{r}{2}`.

    - Hodge-Riemann-Minikowski Relations: Every Lefchetz element `\omega`,
      define quadratic forms on `A^{k}_{\RR}` given by

      .. MATH::

          a \mapsto (-1)^k \deg(a \cdot \omega^{r-2k} \cdot a)

      This quadratic form becomes positive definite upon restriction to the
      kernel of the following map

      .. MATH::

          A^k_\RR \longrightarrow A^{r-k+1}_\RR \\
                a \mapsto a \cdot \omega^{r-2k+1}.

    REFERENCES:

    - [ANR2023]_

    TESTS::

        sage: from sage.categories.kahler_algebras import KahlerAlgebras

        sage: C = KahlerAlgebras(QQ)
        sage: TestSuite(C).run()
    """
    def super_categories(self):
        r"""
        Return the super categories of ``self``.

        EXAMPLES::

            sage: from sage.categories.kahler_algebras import KahlerAlgebras

            sage: C = KahlerAlgebras(QQ); C
            Category of kahler algebras over Rational Field
            sage: sorted(C.super_categories(), key=str)
            [Category of finite dimensional graded algebras with basis over
            Rational Field]
        """
        return [GradedAlgebrasWithBasis(self.base_ring()).FiniteDimensional()]

    class ParentMethods:
        @abstract_method
        def poincare_pairing(self, a, b):
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

        @abstract_method
        def lefschetz_element(self):
            r"""
            Return one Lefschetz element of the given Kähler algebra.

            EXAMPLES::

                sage: U46 = matroids.Uniform(4, 6)
                sage: C = U46.chow_ring(QQ, False)
                sage: w = C.lefschetz_element(); w
                -2*A01 - 2*A02 - 2*A03 - 2*A04 - 2*A05 - 2*A12 - 2*A13 - 2*A14
                - 2*A15 - 2*A23 - 2*A24 - 2*A25 - 2*A34 - 2*A35 - 2*A45 - 6*A012
                - 6*A013 - 6*A014 - 6*A015 - 6*A023 - 6*A024 - 6*A025 - 6*A034
                - 6*A035 - 6*A045 - 6*A123 - 6*A124 - 6*A125 - 6*A134 - 6*A135
                - 6*A145 - 6*A234 - 6*A235 - 6*A245 - 6*A345 - 30*A012345
                sage: basis_deg = {}
                sage: for b in C.basis():
                ....:     deg = b.homogeneous_degree()
                ....:     if deg not in basis_deg:
                ....:         basis_deg[deg] = []
                ....:     basis_deg[deg].append(b)
                sage: m = max(basis_deg); m
                3
                sage: len(basis_deg[1]) == len(basis_deg[2])
                True
                sage: matrix([(w*b).to_vector() for b in basis_deg[1]]).rank()
                36
                sage: len(basis_deg[2])
                36
            """

        def hodge_riemann_relations(self, k):
            r"""
            Return the quadratic form for the corresponding ``k``
            (`< \frac{r}{2}`) for the Kähler algebra, where `r` is the top degree.

            EXAMPLES::

                sage: ch = matroids.Uniform(4, 6).chow_ring(QQ, False)
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
                ValueError: k must be less than r/2 < 2
            """
            r = self.top_degree()
            if k > (r/2):
                raise ValueError("k must be less than r/2 < 2")
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