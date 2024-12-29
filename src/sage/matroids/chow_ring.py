r"""
Chow rings of matroids

AUTHORS:

- Shriya M
"""

from sage.matroids.chow_ring_ideal import ChowRingIdeal_nonaug, AugmentedChowRingIdeal_fy, AugmentedChowRingIdeal_atom_free
from sage.rings.quotient_ring import QuotientRing_generic
from sage.categories.kahler_algebras import KahlerAlgebras
from sage.categories.commutative_rings import CommutativeRings


class ChowRing(QuotientRing_generic):
    r"""
    The Chow ring of a matroid.

    The *Chow ring of the matroid* `M` is defined as the quotient ring

    .. MATH::

        A^*(M)_R := R[x_{F_1}, \ldots, x_{F_k}] / (I_M + J_M),

    where `(I_M + J_M)` is the :class:`Chow ring ideal
    <sage.matroids.chow_ring_ideal.ChowRingIdeal_nonaug>` of matroid `M`.

    The *augmented Chow ring of matroid* `M` has two different presentations
    as quotient rings:

    The *Feitchner-Yuzvinsky presentation* is the quotient ring

    .. MATH::

        A(M)_R := R[y_{e_1}, \ldots, y_{e_n}, x_{F_1}, \ldots, x_{F_k}] / I_{FY}(M),

    where `I_{FY}(M)` is the :class:`Feitchner-Yuzvinsky augmented Chow ring
    ideal <sage.matroids.chow_ring_ideal.AugmentedChowRingIdeal_fy>`
    of matroid `M`.

    The *atom-free presentation* is the quotient ring

    .. MATH::

        A(M)_R := R[x_{F_1}, \ldots, x_{F_k}] / I_{af}(M),

    where `I_{af}(M)` is the :class:`atom-free augmented Chow ring ideal
    <sage.matroids.chow_ring_ideal.AugmentedChowRingIdeal_atom_free>`
    of matroid `M`.

    .. SEEALSO::

        :mod:`sage.matroids.chow_ring_ideal`

    INPUT:

    - ``M`` -- matroid
    - ``R`` -- commutative ring
    - ``augmented`` -- boolean; when ``True``, this is the augmented
      Chow ring and if ``False``, this is the non-augmented Chow ring
    - ``presentation`` -- string (default: ``None``); one of the following
      (ignored if ``augmented=False``)

      * ``"fy"`` - the Feitchner-Yuzvinsky presentation
      * ``"atom-free"`` - the atom-free presentation

    REFERENCES:

    - [FY2004]_
    - [AHK2015]_

    EXAMPLES::

        sage: M1 = matroids.catalog.P8pp()
        sage: ch = M1.chow_ring(QQ, False)
        sage: ch
        Chow ring of P8'': Matroid of rank 4 on 8 elements with 8 nonspanning circuits
        over Rational Field
    """
    def __init__(self, R, M, augmented, presentation=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: ch = matroids.Wheel(3).chow_ring(QQ, False)
            sage: TestSuite(ch).run()
        """
        self._matroid = M
        self._augmented = augmented
        self._presentation = presentation
        if augmented is True:
            if presentation == 'fy':
                self._ideal = AugmentedChowRingIdeal_fy(M, R)
            elif presentation == 'atom-free':
                self._ideal = AugmentedChowRingIdeal_atom_free(M, R)
        else:
            self._ideal = ChowRingIdeal_nonaug(M, R)
        C = CommutativeRings().Quotients() & KahlerAlgebras(R).FiniteDimensional()
        QuotientRing_generic.__init__(self, R=self._ideal.ring(),
                                      I=self._ideal,
                                      names=self._ideal.ring().variable_names(),
                                      category=C)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: M1 = matroids.catalog.Fano()
            sage: ch = M1.chow_ring(QQ, False)
            sage: ch
            Chow ring of Fano: Binary matroid of rank 3 on 7 elements, type (3, 0)
            over Rational Field
        """
        output = "Chow ring of {}".format(self._matroid)
        if self._augmented is True:
            output = "Augmented " + output
            if self._presentation == 'fy':
                output += " in Feitchner-Yuzvinsky presentation"
            elif self._presentation == 'atom-free':
                output += " in atom-free presentation"
        return output + " over " + repr(self.base_ring())

    def _latex_(self):
        r"""
        Return the LaTeX output of the polynomial ring and Chow ring ideal.

        EXAMPLES::

            sage: M1 = matroids.Uniform(2,5)
            sage: ch = M1.chow_ring(QQ, False)
            sage: ch._latex_()
            'A(\\begin{array}{l}\n\\text{\\texttt{U(2,{ }5):{ }Matroid{ }of{ }rank{ }2{ }on{ }5{ }elements{ }with{ }circuit{-}closures}}\\\\\n\\text{\\texttt{{\\char`\\{}2:{ }{\\char`\\{}{\\char`\\{}0,{ }1,{ }2,{ }3,{ }4{\\char`\\}}{\\char`\\}}{\\char`\\}}}}\n\\end{array})_{\\Bold{Q}}'
        """
        from sage.misc.latex import latex
        base = "A({})_{{{}}}"
        if self._augmented:
            base += "^*"
        return base.format(latex(self._matroid), latex(self.base_ring()))

    def matroid(self):
        r"""
        Return the matroid of ``self``.

        EXAMPLES::

            sage: ch = matroids.Uniform(3,6).chow_ring(QQ, True, 'fy')
            sage: ch.matroid()
            U(3, 6): Matroid of rank 3 on 6 elements with circuit-closures
            {3: {{0, 1, 2, 3, 4, 5}}}
        """
        return self._matroid

    def _coerce_map_from_base_ring(self):
        r"""
        Disable the coercion from the base ring from the category.

        TESTS::

            sage: ch = matroids.Wheel(3).chow_ring(QQ, False)
            sage: ch._coerce_map_from_base_ring() is None
            True
        """
        return None  # don't need anything special

    def basis(self):
        r"""
        Return the monomial basis of the given Chow ring.

        EXAMPLES::

            sage: ch = matroids.Uniform(3, 6).chow_ring(QQ, True, 'fy')
            sage: ch.basis()
            Family (1, B1, B1*B012345, B0, B0*B012345, B01, B01^2, B2,
            B2*B012345, B02, B02^2, B12, B12^2, B3, B3*B012345, B03, B03^2,
            B13, B13^2, B23, B23^2, B4, B4*B012345, B04, B04^2, B14, B14^2,
            B24, B24^2, B34, B34^2, B5, B5*B012345, B05, B05^2, B15, B15^2,
            B25, B25^2, B35, B35^2, B45, B45^2, B012345, B012345^2, B012345^3)
            sage: set(ch.defining_ideal().normal_basis()) == set(ch.basis())
            True
            sage: ch = matroids.catalog.Fano().chow_ring(QQ, False)
            sage: ch.basis()
            Family (1, Abcd, Aace, Aabf, Adef, Aadg, Abeg, Acfg, Aabcdefg,
            Aabcdefg^2)
            sage: set(ch.defining_ideal().normal_basis()) == set(ch.basis())
            True
            sage: ch = matroids.Wheel(3).chow_ring(QQ, True, 'atom-free')
            sage: ch.basis()
            Family (1, A0, A0*A012345, A2, A2*A012345, A3, A3*A012345, A23,
            A23^2, A1, A1*A012345, A013, A013^2, A4, A4*A012345, A04, A04^2,
            A124, A124^2, A5, A5*A012345, A025, A025^2, A15, A15^2, A345,
            A345^2, A012345, A012345^2, A012345^3)
            sage: set(ch.defining_ideal().normal_basis()) == set(ch.basis())
            True
        """
        from sage.sets.family import Family
        monomial_basis = self._ideal.normal_basis()
        return Family([self.element_class(self, mon, reduce=False) for mon in monomial_basis])

    def flats_generator(self):
        r"""
        Return the corresponding generators of flats of the Chow ring.

        EXAMPLES::

            sage: ch = matroids.catalog.NonFano().chow_ring(ZZ, True, 'atom-free')
            sage: ch.flats_generator()
            {frozenset({'a'}): Aa,
            frozenset({'b'}): Ab,
            frozenset({'c'}): Ac,
            frozenset({'d'}): Ad,
            frozenset({'e'}): Ae,
            frozenset({'f'}): Af,
            frozenset({'g'}): Ag,
            frozenset({'a', 'b', 'f'}): Aabf,
            frozenset({'a', 'c', 'e'}): Aace,
            frozenset({'a', 'd', 'g'}): Aadg,
            frozenset({'b', 'c', 'd'}): Abcd,
            frozenset({'b', 'e', 'g'}): Abeg,
            frozenset({'c', 'f', 'g'}): Acfg,
            frozenset({'d', 'e'}): Ade,
            frozenset({'d', 'f'}): Adf,
            frozenset({'e', 'f'}): Aef,
            frozenset({'a', 'b', 'c', 'd', 'e', 'f', 'g'}): Aabcdefg}
        """
        flats = [X for i in range(1, self._matroid.rank() + 1)
                 for X in self._matroid.flats(i)]
        gens = self.gens()
        if not (self._augmented and self._presentation == 'fy'):
            return dict(zip(flats, gens))
        flats_gen = {}
        E = list(self.matroid().groundset())
        for i,F in enumerate(flats):
            flats_gen[F] = gens[len(E) + i]
        return flats_gen

    def lefschetz_element(self):
        r"""
        Return one Lefschetz element of the given Chow ring.

        EXAMPLES::

            sage: ch = matroids.catalog.P8pp().chow_ring(QQ, False)
            sage: ch.lefschetz_element()
            -2*Aab - 2*Aac - 2*Aad - 2*Aae - 2*Aaf - 2*Aag - 2*Aah - 2*Abc
            - 2*Abd - 2*Abe - 2*Abf - 2*Abg - 2*Abh - 2*Acd - 2*Ace - 2*Acf
            - 2*Acg - 2*Ach - 2*Ade - 2*Adf - 2*Adg - 2*Adh - 2*Aef - 2*Aeg
            - 2*Aeh - 2*Afg - 2*Afh - 2*Agh - 6*Aabc - 6*Aabd - 6*Aabe
            - 12*Aabfh - 6*Aabg - 6*Aacd - 12*Aacef - 12*Aacgh - 12*Aadeg
            - 6*Aadf - 6*Aadh - 6*Aaeh - 6*Aafg - 6*Abcd - 12*Abceg
            - 6*Abcf - 6*Abch - 12*Abdeh - 12*Abdfg - 6*Abef - 6*Abgh
            - 6*Acde - 12*Acdfh - 6*Acdg - 6*Aceh - 6*Acfg - 6*Adef
            - 6*Adgh - 6*Aefg - 6*Aefh - 6*Aegh - 6*Afgh - 56*Aabcdefgh

        The following example finds the Lefschetz element of the Chow ring
        of the uniform matroid of rank 4 on 5 elements (non-augmented).
        It is then multiplied with the elements of FY-monomial bases of
        different degrees::

            sage: ch = matroids.Uniform(4,5).chow_ring(QQ, False)
            sage: basis_deg = {}
            sage: for b in ch.basis():
            ....:     deg = b.homogeneous_degree()
            ....:     if deg not in basis_deg:
            ....:         basis_deg[deg] = []
            ....:     basis_deg[deg].append(b)
            ....:
            sage: basis_deg
            {0: [1], 1: [A02, A12, A01, A012, A03, A13, A013, A23, A023,
                A123, A04, A14, A014, A24, A024, A124, A34, A034, A134, A234,
                A01234], 2: [A02*A01234, A12*A01234, A01*A01234, A012^2,
                A03*A01234, A13*A01234, A013^2, A23*A01234, A023^2, A123^2,
                A04*A01234, A14*A01234, A014^2, A24*A01234, A024^2, A124^2,
                A34*A01234, A034^2, A134^2, A234^2, A01234^2], 3: [A01234^3]}
            sage: g_eq_maps = {}
            sage: lefschetz_el = ch.lefschetz_element(); lefschetz_el
            -2*A01 - 2*A02 - 2*A03 - 2*A04 - 2*A12 - 2*A13 - 2*A14 - 2*A23
            - 2*A24 - 2*A34 - 6*A012 - 6*A013 - 6*A014 - 6*A023 - 6*A024
            - 6*A034 - 6*A123 - 6*A124 - 6*A134 - 6*A234 - 20*A01234
            sage: for deg in basis_deg:
            ....:     if deg not in g_eq_maps:
            ....:         g_eq_maps[deg] = []
            ....:     g_eq_maps[deg].extend([i*lefschetz_el for i in basis_deg[deg]])
            ....:
            sage: g_eq_maps
            {0: [-2*A01 - 2*A02 - 2*A03 - 2*A04 - 2*A12 - 2*A13 - 2*A14
                - 2*A23 - 2*A24 - 2*A34 - 6*A012 - 6*A013 - 6*A014 - 6*A023
                - 6*A024 - 6*A034 - 6*A123 - 6*A124 - 6*A134 - 6*A234
                - 20*A01234], 1: [2*A012^2 + 2*A023^2 + 2*A024^2
                - 10*A02*A01234 + 2*A01234^2, 2*A012^2 + 2*A123^2 + 2*A124^2
                - 10*A12*A01234 + 2*A01234^2, 2*A012^2 + 2*A013^2 + 2*A014^2
                - 10*A01*A01234 + 2*A01234^2, -6*A012^2 + 2*A01*A01234
                + 2*A02*A01234 + 2*A12*A01234, 2*A013^2 + 2*A023^2 + 2*A034^2
                - 10*A03*A01234 + 2*A01234^2, 2*A013^2 + 2*A123^2 + 2*A134^2
                - 10*A13*A01234 + 2*A01234^2, -6*A013^2 + 2*A01*A01234
                + 2*A03*A01234 + 2*A13*A01234, 2*A023^2 + 2*A123^2 + 2*A234^2
                - 10*A23*A01234 + 2*A01234^2, -6*A023^2 + 2*A02*A01234
                + 2*A03*A01234 + 2*A23*A01234, -6*A123^2 + 2*A12*A01234
                + 2*A13*A01234 + 2*A23*A01234, 2*A014^2 + 2*A024^2 + 2*A034^2
                - 10*A04*A01234 + 2*A01234^2, 2*A014^2 + 2*A124^2 + 2*A134^2
                - 10*A14*A01234 + 2*A01234^2, -6*A014^2 + 2*A01*A01234
                + 2*A04*A01234 + 2*A14*A01234, 2*A024^2 + 2*A124^2 + 2*A234^2
                - 10*A24*A01234 + 2*A01234^2, -6*A024^2 + 2*A02*A01234
                + 2*A04*A01234 + 2*A24*A01234, -6*A124^2 + 2*A12*A01234
                + 2*A14*A01234 + 2*A24*A01234, 2*A034^2 + 2*A134^2 + 2*A234^2
                - 10*A34*A01234 + 2*A01234^2, -6*A034^2 + 2*A03*A01234
                + 2*A04*A01234 + 2*A34*A01234, -6*A134^2 + 2*A13*A01234
                + 2*A14*A01234 + 2*A34*A01234, -6*A234^2 + 2*A23*A01234
                + 2*A24*A01234 + 2*A34*A01234, -2*A01*A01234 - 2*A02*A01234
                - 2*A03*A01234 - 2*A04*A01234 - 2*A12*A01234 - 2*A13*A01234
                - 2*A14*A01234 - 2*A23*A01234 - 2*A24*A01234 - 2*A34*A01234
                - 20*A01234^2], 2: [2*A01234^3, 2*A01234^3, 2*A01234^3,
                6*A01234^3, 2*A01234^3, 2*A01234^3, 6*A01234^3, 2*A01234^3,
                6*A01234^3, 6*A01234^3, 2*A01234^3, 2*A01234^3, 6*A01234^3,
                2*A01234^3, 6*A01234^3, 6*A01234^3, 2*A01234^3, 6*A01234^3,
                6*A01234^3, 6*A01234^3, -20*A01234^3], 3: [0]}

        TESTS::

            sage: U46 = matroids.Uniform(4,6)
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
        w = sum(len(F) * (len(self.matroid().groundset()) - len(F)) * gen for F, gen in self.flats_generator().items())
        return w

    def poincare_pairing(self, el1, el2):
        r"""
        Return the Poincaré pairing of any two elements of the
        Chow ring.

        EXAMPLES::

            sage: ch = matroids.Wheel(3).chow_ring(QQ, True, 'atom-free')
            sage: A0, A1, A2, A3, A4, A5, A013, A025, A04, A124, A15, A23, A345, A012345 = ch.gens()
            sage: u = ch(-1/6*A2*A012345 + 41/48*A012345^2); u
            -1/6*A2*A012345 + 41/48*A012345^2
            sage: v = ch(-A345^2 - 1/4*A345); v
            -A345^2 - 1/4*A345
            sage: ch.poincare_pairing(v, u, ch.matroid().rank())
            3
        """
        r = self._top_degree()
        hom_components1 = el1.lift().homogeneous_components()
        hom_components2 = el2.lift().homogeneous_components()
        new_el = self.base_ring().zero()
        for i in hom_components1:
            for j in hom_components2:
                if r - i not in hom_components2:
                    continue
                new_el += hom_components1[i] * hom_components2[j]
        return new_el.degree()

    class Element(QuotientRing_generic.Element):
        def to_vector(self, order=None):
            r"""
            Return ``self`` as a (dense) free module vector.

            EXAMPLES::

                sage: ch = matroids.Uniform(3, 6).chow_ring(QQ, False)
                sage: v = ch.an_element(); v
                -A01 - A02 - A03 - A04 - A05 - A012345
                sage: v.to_vector()
                (0, -1, -1, 0, -1, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0)
            """
            P = self.parent()
            B = P.basis()
            FM = P._dense_free_module()
            f = self.lift()
            return FM([f.monomial_coefficient(b.lift()) for b in B])

        _vector_ = to_vector

        def monomial_coefficients(self, copy=None):
            r"""
            Return the monomial coefficients of ``self``.

            EXAMPLES::

                sage: ch = matroids.catalog.NonFano().chow_ring(QQ, True, 'atom-free')
                sage: v = ch.an_element(); v
                Aa
                sage: v.monomial_coefficients()
                {0: 0, 1: 1, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0,
                 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0,
                 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0,
                 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0,
                 34: 0, 35: 0}
            """
            B = self.parent().basis()
            f = self.lift()
            return {i: f.monomial_coefficient(b.lift()) for i, b in enumerate(B)}

        def degree(self):
            r"""
            Return the degree of ``self``.

            EXAMPLES::

                sage: ch = matroids.Uniform(3, 6).chow_ring(QQ, False)
                sage: for b in ch.basis():
                ....:     print(b, b.degree())
                1 0
                A01 1
                A02 1
                A12 1
                A03 1
                A13 1
                A23 1
                A04 1
                A14 1
                A24 1
                A34 1
                A05 1
                A15 1
                A25 1
                A35 1
                A45 1
                A012345 1
                A012345^2 2
                sage: v = sum(ch.basis())
                sage: v.degree()
                2
            """
            return self.lift().degree()

        def homogeneous_degree(self):
            r"""
            Return the (homogeneous) degree of ``self`` if homogeneous
            otherwise raise an error.

            EXAMPLES::

                sage: ch = matroids.catalog.Fano().chow_ring(QQ, True, 'fy')
                sage: for b in ch.basis():
                ....:     print(b, b.homogeneous_degree())
                1 0
                Ba 1
                Ba*Babcdefg 2
                Bb 1
                Bb*Babcdefg 2
                Bc 1
                Bc*Babcdefg 2
                Bd 1
                Bd*Babcdefg 2
                Bbcd 1
                Bbcd^2 2
                Be 1
                Be*Babcdefg 2
                Bace 1
                Bace^2 2
                Bf 1
                Bf*Babcdefg 2
                Babf 1
                Babf^2 2
                Bdef 1
                Bdef^2 2
                Bg 1
                Bg*Babcdefg 2
                Badg 1
                Badg^2 2
                Bbeg 1
                Bbeg^2 2
                Bcfg 1
                Bcfg^2 2
                Babcdefg 1
                Babcdefg^2 2
                Babcdefg^3 3
                sage: v = sum(ch.basis()); v
                Babcdefg^3 + Babf^2 + Bace^2 + Badg^2 + Bbcd^2 + Bbeg^2 +
                Bcfg^2 + Bdef^2 + Ba*Babcdefg + Bb*Babcdefg + Bc*Babcdefg +
                Bd*Babcdefg + Be*Babcdefg + Bf*Babcdefg + Bg*Babcdefg +
                Babcdefg^2 + Ba + Bb + Bc + Bd + Be + Bf + Bg + Babf + Bace +
                Badg + Bbcd + Bbeg + Bcfg + Bdef + Babcdefg + 1
                sage: v.homogeneous_degree()
                Traceback (most recent call last):
                ...
                ValueError: element is not homogeneous

            TESTS::

                sage: ch = matroids.Wheel(3).chow_ring(QQ, True, 'atom-free')
                sage: ch.zero().homogeneous_degree()
                Traceback (most recent call last):
                ...
                ValueError: the zero element does not have a well-defined degree
            """
            if not self:
                raise ValueError("the zero element does not have a well-defined degree")
            f = self.lift()
            if not f.is_homogeneous():
                raise ValueError("element is not homogeneous")
            return f.degree()
