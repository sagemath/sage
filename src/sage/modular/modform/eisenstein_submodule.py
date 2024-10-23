# sage.doctest: needs sage.libs.pari
"""
The Eisenstein subspace
"""

from sage.arith.functions import lcm
from sage.arith.misc import euler_phi
from sage.categories.objects import Objects
from sage.matrix.constructor import Matrix
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.rings.integer import Integer
from sage.structure.sequence import Sequence

lazy_import('sage.rings.number_field.number_field', 'CyclotomicField')

from . import eis_series
from . import element
from . import submodule


class EisensteinSubmodule(submodule.ModularFormsSubmodule):
    """
    The Eisenstein submodule of an ambient space of modular forms.
    """
    def __init__(self, ambient_space):
        """
        Return the Eisenstein submodule of the given space.

        EXAMPLES::

            sage: E = ModularForms(23,4).eisenstein_subspace()  # indirect doctest
            sage: E
            Eisenstein subspace of dimension 2 of Modular Forms space of dimension 7
             for Congruence Subgroup Gamma0(23) of weight 4 over Rational Field
            sage: E == loads(dumps(E))
            True
        """
        from sage.misc.verbose import verbose
        verbose('creating eisenstein submodule of %s' % ambient_space)
        d = ambient_space._dim_eisenstein()
        V = ambient_space.module()
        n = V.dimension()
        self._start_position = int(n - d)
        S = V.submodule([V.gen(i) for i in range(n-d,n)], check=False,
                        already_echelonized=True)
        submodule.ModularFormsSubmodule.__init__(self, ambient_space, S)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: E = ModularForms(23,4).eisenstein_subspace() # indirect doctest
            sage: E._repr_()
            'Eisenstein subspace of dimension 2 of Modular Forms space of dimension 7 for Congruence Subgroup Gamma0(23) of weight 4 over Rational Field'
        """
        return "Eisenstein subspace of dimension %s of %s" % (self.dimension(), self.ambient_module())

    def eisenstein_submodule(self):
        """
        Return the Eisenstein submodule of ``self``.
        (Yes, this is just self.)

        EXAMPLES::

            sage: E = ModularForms(23,4).eisenstein_subspace()
            sage: E == E.eisenstein_submodule()
            True
        """
        return self

    @cached_method
    def modular_symbols(self, sign=0):
        r"""
        Return the corresponding space of modular symbols with given sign. This
        will fail in weight 1.

        .. warning::

           If sign != 0, then the space of modular symbols will, in general,
           only correspond to a *subspace* of this space of modular forms.
           This can be the case for both sign +1 or -1.

        EXAMPLES::

            sage: E = ModularForms(11,2).eisenstein_submodule()
            sage: M = E.modular_symbols(); M
            Modular Symbols subspace of dimension 1 of Modular Symbols space
            of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: M.sign()
            0

            sage: M = E.modular_symbols(sign=-1); M
            Modular Symbols subspace of dimension 0 of Modular Symbols space of
            dimension 1 for Gamma_0(11) of weight 2 with sign -1 over Rational Field

            sage: E = ModularForms(1,12).eisenstein_submodule()
            sage: E.modular_symbols()
            Modular Symbols subspace of dimension 1 of Modular Symbols space of
            dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field

            sage: eps = DirichletGroup(13).0
            sage: E = EisensteinForms(eps^2, 2)
            sage: E.modular_symbols()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of
             dimension 4 and level 13, weight 2, character [zeta6], sign 0,
             over Cyclotomic Field of order 6 and degree 2

            sage: E = EisensteinForms(eps, 1); E
            Eisenstein subspace of dimension 1 of Modular Forms space of character
             [zeta12] and weight 1 over Cyclotomic Field of order 12 and degree 4
            sage: E.modular_symbols()
            Traceback (most recent call last):
            ...
            ValueError: the weight must be at least 2
        """
        A = self.ambient_module()
        return A.modular_symbols(sign).eisenstein_submodule()


class EisensteinSubmodule_params(EisensteinSubmodule):

    @cached_method
    def parameters(self):
        r"""
        Return a list of parameters for each Eisenstein series
        spanning ``self``. That is, for each such series, return a triple
        of the form (`\psi`, `\chi`, level), where `\psi` and `\chi`
        are the characters defining the Eisenstein series, and level
        is the smallest level at which this series occurs.

        EXAMPLES::

            sage: ModularForms(24,2).eisenstein_submodule().parameters()
            [(Dirichlet character modulo 24 of conductor 1 mapping 7 |--> 1, 13 |--> 1, 17 |--> 1,
              Dirichlet character modulo 24 of conductor 1 mapping 7 |--> 1, 13 |--> 1, 17 |--> 1, 2),
              ...
              Dirichlet character modulo 24 of conductor 1 mapping 7 |--> 1, 13 |--> 1, 17 |--> 1, 24)]
            sage: EisensteinForms(12,6).parameters()[-1]
            (Dirichlet character modulo 12 of conductor 1 mapping 7 |--> 1, 5 |--> 1,
             Dirichlet character modulo 12 of conductor 1 mapping 7 |--> 1, 5 |--> 1, 12)

            sage: pars = ModularForms(DirichletGroup(24).0,3).eisenstein_submodule().parameters()
            sage: [(x[0].values_on_gens(),x[1].values_on_gens(),x[2]) for x in pars]
            [((1, 1, 1), (-1, 1, 1), 1),
            ((1, 1, 1), (-1, 1, 1), 2),
            ((1, 1, 1), (-1, 1, 1), 3),
            ((1, 1, 1), (-1, 1, 1), 6),
            ((-1, 1, 1), (1, 1, 1), 1),
            ((-1, 1, 1), (1, 1, 1), 2),
            ((-1, 1, 1), (1, 1, 1), 3),
            ((-1, 1, 1), (1, 1, 1), 6)]
            sage: EisensteinForms(DirichletGroup(24).0,1).parameters()
            [(Dirichlet character modulo 24 of conductor 1 mapping 7 |--> 1, 13 |--> 1, 17 |--> 1,
              Dirichlet character modulo 24 of conductor 4 mapping 7 |--> -1, 13 |--> 1, 17 |--> 1, 1),
             (Dirichlet character modulo 24 of conductor 1 mapping 7 |--> 1, 13 |--> 1, 17 |--> 1,
              Dirichlet character modulo 24 of conductor 4 mapping 7 |--> -1, 13 |--> 1, 17 |--> 1, 2),
             (Dirichlet character modulo 24 of conductor 1 mapping 7 |--> 1, 13 |--> 1, 17 |--> 1,
              Dirichlet character modulo 24 of conductor 4 mapping 7 |--> -1, 13 |--> 1, 17 |--> 1, 3),
             (Dirichlet character modulo 24 of conductor 1 mapping 7 |--> 1, 13 |--> 1, 17 |--> 1,
              Dirichlet character modulo 24 of conductor 4 mapping 7 |--> -1, 13 |--> 1, 17 |--> 1, 6)]
        """
        char = self._parameters_character()
        if char is None:
            return eis_series.compute_eisenstein_params(self.level(), self.weight())
        else:
            return eis_series.compute_eisenstein_params(char, self.weight())

    def new_submodule(self, p=None):
        r"""
        Return the new submodule of ``self``.

        EXAMPLES::

            sage: e = EisensteinForms(Gamma0(225), 2).new_submodule(); e
            Modular Forms subspace of dimension 3 of Modular Forms space of dimension 42
             for Congruence Subgroup Gamma0(225) of weight 2 over Rational Field
            sage: e.basis()
            [
            q + O(q^6),
            q^2 + O(q^6),
            q^4 + O(q^6)
            ]
        """

        if p is not None:
            raise NotImplementedError
        return self.submodule([self(x) for x in self._compute_q_expansion_basis(self.sturm_bound(), new=True)], check=False)

    def _parameters_character(self):
        """
        Return the character defining ``self``.

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(33).1,5)._parameters_character()
            Dirichlet character modulo 33 of conductor 11 mapping 23 |--> 1, 13 |--> zeta10
        """
        return self.character()

    def change_ring(self, base_ring):
        """
        Return ``self`` as a module over ``base_ring``.

        EXAMPLES::

            sage: E = EisensteinForms(12,2) ; E
            Eisenstein subspace of dimension 5 of Modular Forms space of dimension 5
             for Congruence Subgroup Gamma0(12) of weight 2 over Rational Field
            sage: E.basis()
            [
            1 + O(q^6),
            q + 6*q^5 + O(q^6),
            q^2 + O(q^6),
            q^3 + O(q^6),
            q^4 + O(q^6)
            ]
            sage: E.change_ring(GF(5))
            Eisenstein subspace of dimension 5 of Modular Forms space of dimension 5
             for Congruence Subgroup Gamma0(12) of weight 2 over Finite Field of size 5
            sage: E.change_ring(GF(5)).basis()
            [
            1 + O(q^6),
            q + q^5 + O(q^6),
            q^2 + O(q^6),
            q^3 + O(q^6),
            q^4 + O(q^6)
            ]
        """
        if base_ring == self.base_ring():
            return self
        A = self.ambient_module()
        B = A.change_ring(base_ring)
        return B.eisenstein_submodule()

    def eisenstein_series(self):
        """
        Return the Eisenstein series that span this space (over the
        algebraic closure).

        EXAMPLES::

            sage: EisensteinForms(11,2).eisenstein_series()
            [
            5/12 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + O(q^6)
            ]
            sage: EisensteinForms(1,4).eisenstein_series()
            [
            1/240 + q + 9*q^2 + 28*q^3 + 73*q^4 + 126*q^5 + O(q^6)
            ]
            sage: EisensteinForms(1,24).eisenstein_series()
            [
            236364091/131040 + q + 8388609*q^2 + 94143178828*q^3
                + 70368752566273*q^4 + 11920928955078126*q^5 + O(q^6)
            ]
            sage: EisensteinForms(5,4).eisenstein_series()
            [
            1/240 + q + 9*q^2 + 28*q^3 + 73*q^4 + 126*q^5 + O(q^6),
            1/240 + q^5 + O(q^6)
            ]
            sage: EisensteinForms(13,2).eisenstein_series()
            [
            1/2 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + O(q^6)
            ]

            sage: E = EisensteinForms(Gamma1(7),2)
            sage: E.set_precision(4)
            sage: E.eisenstein_series()
            [
            1/4 + q + 3*q^2 + 4*q^3 + O(q^4),
            1/7*zeta6 - 3/7 + q + (-2*zeta6 + 1)*q^2 + (3*zeta6 - 2)*q^3 + O(q^4),
            q + (-zeta6 + 2)*q^2 + (zeta6 + 2)*q^3 + O(q^4),
            -1/7*zeta6 - 2/7 + q + (2*zeta6 - 1)*q^2 + (-3*zeta6 + 1)*q^3 + O(q^4),
            q + (zeta6 + 1)*q^2 + (-zeta6 + 3)*q^3 + O(q^4)
            ]

            sage: eps = DirichletGroup(13).0^2
            sage: ModularForms(eps,2).eisenstein_series()
            [
            -7/13*zeta6 - 11/13 + q + (2*zeta6 + 1)*q^2 + (-3*zeta6 + 1)*q^3
                + (6*zeta6 - 3)*q^4 - 4*q^5 + O(q^6),
            q + (zeta6 + 2)*q^2 + (-zeta6 + 3)*q^3 + (3*zeta6 + 3)*q^4 + 4*q^5 + O(q^6)
            ]

            sage: M = ModularForms(19,3).eisenstein_subspace()
            sage: M.eisenstein_series()
            [
            ]

            sage: M = ModularForms(DirichletGroup(13).0, 1)
            sage: M.eisenstein_series()
            [
            -1/13*zeta12^3 + 6/13*zeta12^2 + 4/13*zeta12 + 2/13 + q + (zeta12 + 1)*q^2
                + zeta12^2*q^3 + (zeta12^2 + zeta12 + 1)*q^4 + (-zeta12^3 + 1)*q^5 + O(q^6)
            ]

            sage: M = ModularForms(GammaH(15, [4]), 4)
            sage: M.eisenstein_series()
            [
            1/240 + q + 9*q^2 + 28*q^3 + 73*q^4 + 126*q^5 + O(q^6),
            1/240 + q^3 + O(q^6),
            1/240 + q^5 + O(q^6),
            1/240 + O(q^6),
            1 + q - 7*q^2 - 26*q^3 + 57*q^4 + q^5 + O(q^6),
            1 + q^3 + O(q^6),
            q + 7*q^2 + 26*q^3 + 57*q^4 + 125*q^5 + O(q^6),
            q^3 + O(q^6)
            ]
        """
        P = self.parameters()
        E = Sequence([element.EisensteinSeries(self.change_ring(chi.base_ring()),
              None, t, chi, psi) for chi, psi, t in P],
              immutable=True, cr=True, universe=Objects())
        assert len(E) == self.dimension(), "bug in enumeration of Eisenstein series."
        return E

    def new_eisenstein_series(self):
        r"""
        Return a list of the Eisenstein series in this space that are new.

        EXAMPLES::

            sage: E = EisensteinForms(25, 4)
            sage: E.new_eisenstein_series()
            [q + 7*zeta4*q^2 - 26*zeta4*q^3 - 57*q^4 + O(q^6),
             q - 9*q^2 - 28*q^3 + 73*q^4 + O(q^6),
             q - 7*zeta4*q^2 + 26*zeta4*q^3 - 57*q^4 + O(q^6)]
        """

        return [x for x in self.eisenstein_series() if x.new_level() == self.level()]

    def _compute_q_expansion_basis(self, prec=None, new=False):
        """
        Compute a `q`-expansion basis for ``self`` to precision ``prec``.

        EXAMPLES::

            sage: EisensteinForms(22,4)._compute_q_expansion_basis(6)
            [1 + O(q^6),
            q + 28*q^3 - 8*q^4 + 126*q^5 + O(q^6),
            q^2 + 9*q^4 + O(q^6),
            O(q^6)]
            sage: EisensteinForms(22,4)._compute_q_expansion_basis(15)
            [1 + O(q^15),
             q + 28*q^3 - 8*q^4 + 126*q^5 + 344*q^7 - 72*q^8 + 757*q^9 - 224*q^12 + 2198*q^13 + O(q^15),
             q^2 + 9*q^4 + 28*q^6 + 73*q^8 + 126*q^10 + 252*q^12 + 344*q^14 + O(q^15),
             q^11 + O(q^15)]
        """
        if prec is None:
            prec = self.prec()
        else:
            prec = Integer(prec)

        if new:
            E = self.new_eisenstein_series()
        else:
            E = self.eisenstein_series()
        K = self.base_ring()
        V = K**prec
        G = []
        for e in E:
            f = e.q_expansion(prec)
            w = f.padded_list(prec)
            L = f.base_ring()
            if K.has_coerce_map_from(L):
                G.append(V(w))
            else:
                # restrict scalars from L to K
                r,d = cyclotomic_restriction(L,K)
                s = [r(x) for x in w]
                for i in range(d):
                    G.append(V([x[i] for x in s]))

        W = V.submodule(G, check=False)
        R = self._q_expansion_ring()
        X = [R(f.list(), prec) for f in W.basis()]
        if not new:
            return X + [R(0,prec)]*(self.dimension() - len(X))
        else:
            return X

    def _q_expansion(self, element, prec):
        """
        Compute a `q`-expansion for a given element of self, expressed
        as a vector of coefficients for the basis vectors of self,
        viewing ``self`` as a subspace of the corresponding space of
        modular forms.

        EXAMPLES::

            sage: E = EisensteinForms(17,4)
            sage: (11*E.0 + 3*E.1).q_expansion(20)
            11 + 3*q + 27*q^2 + 84*q^3 + 219*q^4 + 378*q^5 + 756*q^6 + 1032*q^7 + 1755*q^8 + 2271*q^9 + 3402*q^10 + 3996*q^11 + 6132*q^12 + 6594*q^13 + 9288*q^14 + 10584*q^15 + 14043*q^16 + 17379*q^17 + 20439*q^18 + 20580*q^19 + O(q^20)
            sage: E._q_expansion([0,0,0,0,11,3],20)
            11 + 3*q + 27*q^2 + 84*q^3 + 219*q^4 + 378*q^5 + 756*q^6 + 1032*q^7 + 1755*q^8 + 2271*q^9 + 3402*q^10 + 3996*q^11 + 6132*q^12 + 6594*q^13 + 9288*q^14 + 10584*q^15 + 14043*q^16 + 17379*q^17 + 20439*q^18 + 20580*q^19 + O(q^20)
        """
        B = self.q_expansion_basis(prec)
        f = self._q_expansion_zero()
        for i in range(self._start_position, len(element)):
            if element[i] != 0:
                f += element[i] * B[i - self._start_position]
        return f


class EisensteinSubmodule_g0_Q(EisensteinSubmodule_params):
    r"""
    Space of Eisenstein forms for `\Gamma_0(N)`.
    """
    def _pari_init_(self):
        """
        Conversion to Pari.

        EXAMPLES::

            sage: E = EisensteinForms(17,4)
            sage: pari.mfdim(E)
            2
            sage: pari.mfparams(E)
            [17, 4, 1, 3, t - 1]
        """
        from sage.libs.pari import pari
        return pari.mfinit([self.level(), self.weight()], 3)


class EisensteinSubmodule_gH_Q(EisensteinSubmodule_params):
    r"""
    Space of Eisenstein forms for `\Gamma_H(N)`.
    """
    def _parameters_character(self):
        """
        Return the character defining ``self``. Since ``self`` is
        a space of Eisenstein forms on GammaH(N) rather than a space with fixed
        character, we return the group GammaH(N) itself.

        EXAMPLES::

            sage: EisensteinForms(GammaH(9, [4]),4)._parameters_character()
            Congruence Subgroup Gamma_H(9) with H generated by [4]
        """
        return self.group()

    def _convert_matrix_from_modsyms_eis(self, A):
        r"""
        Given a matrix acting on the space of modular symbols corresponding to
        this space, calculate the matrix of the operator it induces on this
        space itself. Used for Hecke and diamond operators.

        This is a minor modification of the code used for cusp forms, which is
        required because modular symbols "don't see the constant term": the
        modular symbol method calculates the matrix of the operator with
        respect to the unique basis of the modular forms space for which the
        *non-constant* coefficients are in echelon form, and we need to modify
        this to get a matrix with respect to the basis we're actually using.

        EXAMPLES::

            sage: EisensteinForms(Gamma1(6), 3).hecke_matrix(3) # indirect doctest
            [ 1  0 72  0]
            [ 0  0 36 -9]
            [ 0  0  9  0]
            [ 0  1 -4 10]
        """
        from .cuspidal_submodule import _convert_matrix_from_modsyms
        symbs = self.modular_symbols(sign=0)
        d = self.rank()
        wrong_mat, pivs = _convert_matrix_from_modsyms(symbs, A)
        c = Matrix(self.base_ring(), d, [self.basis()[i][j+1] for i in range(d) for j in pivs])
        return c * wrong_mat * ~c

    def _compute_hecke_matrix(self, n, bound=None):
        r"""
        Calculate the matrix of the Hecke operator `T_n` acting on this
        space, via modular symbols.

        INPUT:

        - ``n`` -- positive integer

        - ``bound`` -- integer such that any element of this space with
          coefficients a_1, ..., a_b all zero must be the zero
          element. If this turns out not to be true, the code will
          increase the bound and try again. Setting bound = None is
          equivalent to setting bound = self.dimension().

        OUTPUT: matrix (over `\QQ`)

        ALGORITHM:

            This uses the usual pairing between modular symbols and
            modular forms, but in a slightly non-standard way. As for
            cusp forms, we can find a basis for this space made up of
            forms with `q`-expansions `c_m(f) = a_{i,j}(T_m)`, where
            `T_m` denotes the matrix of the Hecke operator on the
            corresponding modular symbols space. Then `c_m(T_n f) =
            a_{i,j}(T_n* T_m)`. But we can't find the constant terms
            by this method, so an extra step is required.

        EXAMPLES::

            sage: EisensteinForms(Gamma1(6), 3).hecke_matrix(3) # indirect doctest
            [ 1  0 72  0]
            [ 0  0 36 -9]
            [ 0  0  9  0]
            [ 0  1 -4 10]
        """
        symbs = self.modular_symbols(sign=0)
        T = symbs.hecke_matrix(n)
        return self._convert_matrix_from_modsyms_eis(T)

    def _compute_diamond_matrix(self, d):
        r"""
        Calculate the matrix of the diamond bracket operator <d> on this space,
        using modular symbols.

        EXAMPLES::

            sage: E = EisensteinForms(Gamma1(7), 3)
            sage: E._compute_diamond_matrix(3)
            [  27  126  294  770 2142 3528]
            [56/3   85  200  530 1445 2408]
            [11/3   14   22   66  233  392]
            [  -1   -3   -3  -11  -51  -87]
            [  -1   -4   -7  -20  -67 -112]
            [-1/3   -2   -6  -15  -34  -56]
        """
        symbs = self.modular_symbols(sign=0)
        T = symbs.diamond_bracket_matrix(d)
        return self._convert_matrix_from_modsyms_eis(T)


class EisensteinSubmodule_g1_Q(EisensteinSubmodule_gH_Q):
    r"""
    Space of Eisenstein forms for `\Gamma_1(N)`.
    """
    def _parameters_character(self):
        r"""
        Return the character defining ``self``.

        Since ``self`` is a space of Eisenstein
        forms on `\Gamma_1(N)`, all characters modulo the level are possible,
        so we return the level.

        EXAMPLES::

            sage: EisensteinForms(Gamma1(7),4)._parameters_character()
            7
        """
        return self.level()


class EisensteinSubmodule_eps(EisensteinSubmodule_params):
    """
    Space of Eisenstein forms with given Dirichlet character.

    EXAMPLES::

        sage: e = DirichletGroup(27,CyclotomicField(3)).0**2
        sage: M = ModularForms(e,2,prec=10).eisenstein_subspace()
        sage: M.dimension()
        6

        sage: M.eisenstein_series()
        [
        -1/3*zeta6 - 1/3 + q + (2*zeta6 - 1)*q^2 + q^3
             + (-2*zeta6 - 1)*q^4 + (-5*zeta6 + 1)*q^5 + O(q^6),
        -1/3*zeta6 - 1/3 + q^3 + O(q^6),
        q + (-2*zeta6 + 1)*q^2 + (-2*zeta6 - 1)*q^4 + (5*zeta6 - 1)*q^5 + O(q^6),
        q + (zeta6 + 1)*q^2 + 3*q^3 + (zeta6 + 2)*q^4 + (-zeta6 + 5)*q^5 + O(q^6),
        q^3 + O(q^6),
        q + (-zeta6 - 1)*q^2 + (zeta6 + 2)*q^4 + (zeta6 - 5)*q^5 + O(q^6)
        ]
        sage: M.eisenstein_subspace().T(2).matrix().fcp()
        (x + 2*zeta3 + 1) * (x + zeta3 + 2) * (x - zeta3 - 2)^2 * (x - 2*zeta3 - 1)^2
        sage: ModularSymbols(e,2).eisenstein_subspace().T(2).matrix().fcp()
        (x + 2*zeta3 + 1) * (x + zeta3 + 2) * (x - zeta3 - 2)^2 * (x - 2*zeta3 - 1)^2

        sage: M.basis()
        [
        1 - 3*zeta3*q^6 + (-2*zeta3 + 2)*q^9 + O(q^10),
        q + (5*zeta3 + 5)*q^7 + O(q^10),
        q^2 - 2*zeta3*q^8 + O(q^10),
        q^3 + (zeta3 + 2)*q^6 + 3*q^9 + O(q^10),
        q^4 - 2*zeta3*q^7 + O(q^10),
        q^5 + (zeta3 + 1)*q^8 + O(q^10)
        ]
    """
    # TODO
    # def _compute_q_expansion_basis(self, prec):
    #     B = EisensteinSubmodule_params._compute_q_expansion_basis(self, prec)
    #     raise NotImplementedError("must restrict scalars down correctly.")


def cyclotomic_restriction(L, K):
    r"""
    Given two cyclotomic fields `L` and `K`, compute the compositum
    `M` of `K` and `L`, and return a function `f` and the index `[M:K]`.

    The function `f` is a map that acts as follows (here `M =\QQ(\zeta_m)`):

        INPUT: element alpha in `L`
        OUTPUT: a polynomial `f(x)` in `K[x]` such that `f(\zeta_m) = \alpha`,
        where we view alpha as living in `M`. (Note that `\zeta_m` generates
        `M`, not `L`.)

    EXAMPLES::

        sage: L = CyclotomicField(12); N = CyclotomicField(33); M = CyclotomicField(132)
        sage: z, n = sage.modular.modform.eisenstein_submodule.cyclotomic_restriction(L, N)
        sage: n
        2

        sage: z(L.0)
        -zeta33^19*x
        sage: z(L.0)(M.0)
        zeta132^11

        sage: z(L.0^3 - L.0 + 1)
        (zeta33^19 + zeta33^8)*x + 1
        sage: z(L.0^3 - L.0 + 1)(M.0)
        zeta132^33 - zeta132^11 + 1
        sage: z(L.0^3 - L.0 + 1)(M.0) - M(L.0^3 - L.0 + 1)
        0
    """
    if not L.has_coerce_map_from(K):
        M = CyclotomicField(lcm(L.zeta_order(), K.zeta_order()))
        f = cyclotomic_restriction_tower(M,K)

        def g(x):
            r"""
            Function returned by cyclotomic restriction.

            EXAMPLES::

                sage: L = CyclotomicField(12)
                sage: N = CyclotomicField(33)
                sage: g, n = sage.modular.modform.eisenstein_submodule.cyclotomic_restriction(L,N)
                sage: g(L.0)
                -zeta33^19*x
            """
            return f(M(x))
        return g, euler_phi(M.zeta_order())//euler_phi(K.zeta_order())
    else:
        return cyclotomic_restriction_tower(L,K), \
               euler_phi(L.zeta_order())//euler_phi(K.zeta_order())


def cyclotomic_restriction_tower(L, K):
    r"""
    Suppose `L/K` is an extension of cyclotomic fields and `L=Q(\zeta_m)`.
    This function computes a map with the following property:

        INPUT: element alpha in `L`
        OUTPUT: a polynomial `f(x)` in `K[x]` such that `f(\zeta_m) = alpha`

    EXAMPLES::

        sage: L = CyclotomicField(12) ; K = CyclotomicField(6)
        sage: z = sage.modular.modform.eisenstein_submodule.cyclotomic_restriction_tower(L,K)
        sage: z(L.0)
        x
        sage: z(L.0^2+L.0)
        x + zeta6
    """
    if not L.has_coerce_map_from(K):
        raise ValueError("K must be contained in L")
    f = L.defining_polynomial()
    R = K['x']
    g = R(f)
    h_ls = [t[0] for t in g.factor() if t[0](L.gen(0)) == 0]
    if not h_ls:
        raise ValueError(r"K (= Q(\zeta_%s)) is not contained in L (= Q(\zeta_%s))" % (K._n(), L._n()))
    h = h_ls[0]

    def z(a):
        """
        Function returned by cyclotomic_restriction_tower.

        EXAMPLES::

            sage: L = CyclotomicField(121); K = CyclotomicField(11)
            sage: z = sage.modular.modform.eisenstein_submodule.cyclotomic_restriction_tower(L,K)
            sage: z(L.0)
            x
            sage: z(L.0^11)
            zeta11
        """
        return R(a.polynomial()) % h
    return z
