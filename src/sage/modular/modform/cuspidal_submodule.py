"""
The Cuspidal Subspace

EXAMPLES::

    sage: S = CuspForms(SL2Z,12); S
    Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for
    Modular Group SL(2,Z) of weight 12 over Rational Field
    sage: S.basis()
    [
    q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
    ]

    sage: S = CuspForms(Gamma0(33),2); S
    Cuspidal subspace of dimension 3 of Modular Forms space of dimension 6 for
    Congruence Subgroup Gamma0(33) of weight 2 over Rational Field
    sage: S.basis()
    [
    q - q^5 + O(q^6),
    q^2 - q^4 - q^5 + O(q^6),
    q^3 + O(q^6)
    ]

    sage: S = CuspForms(Gamma1(3),6); S
    Cuspidal subspace of dimension 1 of Modular Forms space of dimension 3 for
    Congruence Subgroup Gamma1(3) of weight 6 over Rational Field
    sage: S.basis()
    [
    q - 6*q^2 + 9*q^3 + 4*q^4 + 6*q^5 + O(q^6)
    ]
"""

#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################
from sage.rings.all     import QQ, Integer
from sage.misc.all      import cached_method
from sage.misc.verbose  import verbose
from sage.matrix.all    import Matrix, identity_matrix

from .submodule import ModularFormsSubmodule
from . import vm_basis
from . import weight1

class CuspidalSubmodule(ModularFormsSubmodule):
    """
    Base class for cuspidal submodules of ambient spaces of modular forms.
    """
    def __init__(self, ambient_space):
        """
        The cuspidal submodule of an ambient space of modular forms.

        EXAMPLES::

            sage: S = CuspForms(SL2Z,12); S
            Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for
            Modular Group SL(2,Z) of weight 12 over Rational Field
            sage: S.basis()
            [
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
            ]

            sage: S = CuspForms(Gamma0(33),2); S
            Cuspidal subspace of dimension 3 of Modular Forms space of dimension 6 for
            Congruence Subgroup Gamma0(33) of weight 2 over Rational Field
            sage: S.basis()
            [
            q - q^5 + O(q^6),
            q^2 - q^4 - q^5 + O(q^6),
            q^3 + O(q^6)
            ]

            sage: S = CuspForms(Gamma1(3),6); S
            Cuspidal subspace of dimension 1 of Modular Forms space of dimension 3 for
            Congruence Subgroup Gamma1(3) of weight 6 over Rational Field
            sage: S.basis()
            [
            q - 6*q^2 + 9*q^3 + 4*q^4 + 6*q^5 + O(q^6)
            ]
            sage: S == loads(dumps(S))
            True
        """
        from sage.misc.verbose import verbose
        verbose('creating cuspidal submodule of %s'%ambient_space)
        d = ambient_space._dim_cuspidal()
        V = ambient_space.module()
        G = [V.gen(i) for i in range(d)]
        S = V.submodule(G, check=False, already_echelonized=True)
        ModularFormsSubmodule.__init__(self, ambient_space, S)

    def _compute_q_expansion_basis(self, prec):
        r"""
        Compute a basis of q-expansions of self to the given precision. Not
        implemented in this abstract base class.

        EXAMPLES::

            sage: M = CuspForms(GammaH(11,[2]), 2)
            sage: sage.modular.modform.cuspidal_submodule.CuspidalSubmodule._compute_q_expansion_basis(M, 6)
            Traceback (most recent call last):
            ...
            NotImplementedError: q-expansion basis not implemented for "Cuspidal subspace of ..."
        """
        raise NotImplementedError('q-expansion basis not implemented for "%s"' % self)

    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES::

            sage: S = CuspForms(Gamma1(3),6); S._repr_()
            'Cuspidal subspace of dimension 1 of Modular Forms space of dimension 3 for Congruence Subgroup Gamma1(3) of weight 6 over Rational Field'
        """
        return "Cuspidal subspace of dimension %s of %s"%(self.dimension(), self.ambient_module())

    def is_cuspidal(self):
        """
        Return True since spaces of cusp forms are cuspidal.

        EXAMPLES::

            sage: CuspForms(4,10).is_cuspidal()
            True
        """
        return True

    @cached_method
    def modular_symbols(self, sign=0):
        """
        Return the corresponding space of modular symbols with the given sign.

        EXAMPLES::

            sage: S = ModularForms(11,2).cuspidal_submodule()
            sage: S.modular_symbols()
            Modular Symbols subspace of dimension 2 of Modular Symbols space
            of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field

            sage: S.modular_symbols(sign=-1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space
            of dimension 1 for Gamma_0(11) of weight 2 with sign -1 over Rational Field

            sage: M = S.modular_symbols(sign=1); M
            Modular Symbols subspace of dimension 1 of Modular Symbols space of
            dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
            sage: M.sign()
            1

            sage: S = ModularForms(1,12).cuspidal_submodule()
            sage: S.modular_symbols()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of
            dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field

            sage: eps = DirichletGroup(13).0
            sage: S = CuspForms(eps^2, 2)

            sage: S.modular_symbols(sign=0)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 and level 13, weight 2, character [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2

            sage: S.modular_symbols(sign=1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 and level 13, weight 2, character [zeta6], sign 1, over Cyclotomic Field of order 6 and degree 2

            sage: S.modular_symbols(sign=-1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 1 and level 13, weight 2, character [zeta6], sign -1, over Cyclotomic Field of order 6 and degree 2
        """
        A = self.ambient_module()
        return A.modular_symbols(sign).cuspidal_submodule()


    def change_ring(self, R):
        r"""
        Change the base ring of self to R, when this makes sense. This differs
        from :meth:`~sage.modular.modform.space.ModularFormsSpace.base_extend`
        in that there may not be a canonical map from self to the new space, as
        in the first example below. If this space has a character then this may
        fail when the character cannot be defined over R, as in the second
        example.

        EXAMPLES::

            sage: chi = DirichletGroup(109, CyclotomicField(3)).0
            sage: S9 = CuspForms(chi, 2, base_ring = CyclotomicField(9)); S9
            Cuspidal subspace of dimension 8 of Modular Forms space of dimension 10, character [zeta3 + 1] and weight 2 over Cyclotomic Field of order 9 and degree 6
            sage: S9.change_ring(CyclotomicField(3))
            Cuspidal subspace of dimension 8 of Modular Forms space of dimension 10, character [zeta3 + 1] and weight 2 over Cyclotomic Field of order 3 and degree 2
            sage: S9.change_ring(QQ)
            Traceback (most recent call last):
            ...
            ValueError: Space cannot be defined over Rational Field
        """
        return self.ambient_module().change_ring(R).cuspidal_submodule()

class CuspidalSubmodule_R(CuspidalSubmodule):
    """
    Cuspidal submodule over a non-minimal base ring.
    """
    def _compute_q_expansion_basis(self, prec):
        r"""
        EXAMPLES::

            sage: CuspForms(Gamma1(13), 2, base_ring=QuadraticField(-7, 'a')).q_expansion_basis() # indirect doctest
            [
            q - 4*q^3 - q^4 + 3*q^5 + O(q^6),
            q^2 - 2*q^3 - q^4 + 2*q^5 + O(q^6)
            ]
        """
        return ModularFormsSubmodule._compute_q_expansion_basis(self, prec)


class CuspidalSubmodule_modsym_qexp(CuspidalSubmodule):
    """
    Cuspidal submodule with q-expansions calculated via modular symbols.
    """
    def _compute_q_expansion_basis(self, prec=None):
        """
        Compute q-expansions of a basis for self (via modular symbols).

        EXAMPLES::

            sage: sage.modular.modform.cuspidal_submodule.CuspidalSubmodule_modsym_qexp(ModularForms(11,2))._compute_q_expansion_basis()
            [
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6)
            ]
        """
        if prec is None:
            prec = self.prec()
        else:
            prec = Integer(prec)
        if self.dimension() == 0:
            return []
        M = self.modular_symbols(sign = 1)
        return M.q_expansion_basis(prec)

    def _compute_hecke_matrix_prime(self, p):
        """
        Compute the matrix of a Hecke operator.

        EXAMPLES::

            sage: C=CuspForms(38, 2)
            sage: C._compute_hecke_matrix_prime(7)
            [-1  0  0  0]
            [ 0 -1  0  0]
            [-2 -2  1 -2]
            [ 2  2 -2  1]
        """
        A = self.modular_symbols(sign=1)
        T = A.hecke_matrix(p)
        return _convert_matrix_from_modsyms(A, T)[0]

    def hecke_polynomial(self, n, var='x'):
        r"""
        Return the characteristic polynomial of the Hecke operator T_n on this
        space. This is computed via modular symbols, and in particular is
        faster to compute than the matrix itself.

        EXAMPLES::

            sage: CuspForms(105, 2).hecke_polynomial(2, 'y')
            y^13 + 5*y^12 - 4*y^11 - 52*y^10 - 34*y^9 + 174*y^8 + 212*y^7 - 196*y^6 - 375*y^5 - 11*y^4 + 200*y^3 + 80*y^2

        Check that this gives the same answer as computing the Hecke matrix::

            sage: CuspForms(105, 2).hecke_matrix(2).charpoly(var='y')
            y^13 + 5*y^12 - 4*y^11 - 52*y^10 - 34*y^9 + 174*y^8 + 212*y^7 - 196*y^6 - 375*y^5 - 11*y^4 + 200*y^3 + 80*y^2

        Check that :trac:`21546` is fixed (this example used to take about 5 hours)::

            sage: CuspForms(1728, 2).hecke_polynomial(2) # long time (20 sec)
            x^253 + x^251 - 2*x^249
        """
        return self.modular_symbols(sign=1).hecke_polynomial(n, var)

    def new_submodule(self, p=None):
        r"""
        Return the new subspace of this space of cusp forms. This is computed
        using modular symbols.

        EXAMPLES::

            sage: CuspForms(55).new_submodule()
            Modular Forms subspace of dimension 3 of Modular Forms space of dimension 8 for Congruence Subgroup Gamma0(55) of weight 2 over Rational Field
        """
        symbs = self.modular_symbols(sign=1).new_subspace(p)
        bas = []
        for x in symbs.q_expansion_basis(self.sturm_bound()):
            bas.append(self(x))
        return self.submodule(bas, check=False)


class CuspidalSubmodule_level1_Q(CuspidalSubmodule):
    r"""
    Space of cusp forms of level 1 over `\QQ`.
    """
    def _compute_q_expansion_basis(self, prec=None):
        """
        Compute q-expansions of a basis for self.

        EXAMPLES::

            sage: sage.modular.modform.cuspidal_submodule.CuspidalSubmodule_level1_Q(ModularForms(1,12))._compute_q_expansion_basis()
            [
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
            ]
        """
        if prec is None:
            prec = self.prec()
        else:
            prec = Integer(prec)
        return vm_basis.victor_miller_basis(self.weight(), prec,
                                            cusp_only=True, var='q')

    def _pari_init_(self):
        """
        Conversion to Pari.

        EXAMPLES::

            sage: A = ModularForms(1,12).cuspidal_submodule()
            sage: pari.mfparams(A)
            [1, 12, 1, 1, t - 1]
            sage: pari.mfdim(A)
            1
        """
        from sage.libs.pari import pari
        return pari.mfinit([self.level(), self.weight()], 1)


class CuspidalSubmodule_wt1_eps(CuspidalSubmodule):
    r"""
    Space of cusp forms of weight 1 with specified character.
    """

    def _compute_q_expansion_basis(self, prec=None):
        r"""
        Compute q-expansion basis using Schaeffer's algorithm.

        EXAMPLES::

            sage: CuspForms(DirichletGroup(23, QQ).0, 1).q_echelon_basis() # indirect doctest
            [
            q - q^2 - q^3 + O(q^6)
            ]
        """
        if prec is None:
            prec = self.prec()
        else:
            prec = Integer(prec)
        chi = self.character()
        return [weight1.modular_ratio_to_prec(chi, f, prec) for f in
            weight1.hecke_stable_subspace(chi)]


class CuspidalSubmodule_wt1_gH(CuspidalSubmodule):
    r"""
    Space of cusp forms of weight 1 for a GammaH group.
    """

    def _compute_q_expansion_basis(self, prec=None):
        r"""
        Compute q-expansion basis using Schaeffer's algorithm.

        EXAMPLES::

            sage: CuspForms(GammaH(31, [7]), 1).q_expansion_basis() # indirect doctest
            [
            q - q^2 - q^5 + O(q^6)
            ]

        A more elaborate example (two Galois-conjugate characters each giving a
        2-dimensional space)::

            sage: CuspForms(GammaH(124, [85]), 1).q_expansion_basis()
            [
            q - q^4 - q^6 + O(q^7),
            q^2 + O(q^7),
            q^3 + O(q^7),
            q^5 - q^6 + O(q^7)
            ]
        """
        if prec is None:
            prec = self.prec()
        else:
            prec = Integer(prec)

        chars=self.group().characters_mod_H(sign=-1, galois_orbits=True)

        B = []
        dim = 0
        for c in chars:
            chi = c.minimize_base_ring()
            Bchi = [weight1.modular_ratio_to_prec(chi, f, prec)
                for f in weight1.hecke_stable_subspace(chi) ]
            if Bchi == []:
                continue
            if chi.base_ring() == QQ:
                B += [f.padded_list(prec) for f in Bchi]
                dim += len(Bchi)
            else:
                d = chi.base_ring().degree()
                dim += d * len(Bchi)
                for f in Bchi:
                    w = f.padded_list(prec)
                    for i in range(d):
                        B.append([x[i] for x in w])

        basis_mat = Matrix(QQ, B)
        if basis_mat.is_zero():
            return []
        # Daft thing: "transformation=True" parameter to echelonize
        # is ignored for rational matrices!
        big_mat = basis_mat.augment(identity_matrix(dim))
        big_mat.echelonize()
        c = big_mat.pivots()[-1]

        echelon_basis_mat = big_mat[:, :prec]

        R = self._q_expansion_ring()

        if c >= prec:
            verbose("Precision %s insufficient to determine basis" % prec, level=1)
        else:
            verbose("Minimal precision for basis: %s" % (c+1), level=1)
            t = big_mat[:, prec:]
            assert echelon_basis_mat == t * basis_mat
            self.__transformation_matrix = t
            self._char_basis = [R(f.list(), c+1) for f in basis_mat.rows()]

        return [R(f.list(), prec) for f in echelon_basis_mat.rows() if f != 0]

    def _transformation_matrix(self):
        r"""
        Return the transformation matrix between a basis of Galois orbits of
        forms with character, and the echelon-form basis.

        EXAMPLES::

            sage: CuspForms(GammaH(31, [7]), 1)._transformation_matrix()
            [1]
            sage: CuspForms(GammaH(124, [33]), 1)._transformation_matrix() # long time
            [ 1  1  0  0  0  0  1]
            [ 0  0  0  0  0  1  0]
            [ 1  0  1  1 -1 -1  1]
            [ 0  0  0  0  0  0  1]
            [ 0  1  0  0  0  0  0]
            [ 1  0  0  0 -1  0  1]
            [ 0  0  1  0  0 -1  0]
        """
        self.q_expansion_basis() # triggers iterative computation
        return self.__transformation_matrix

    def _compute_diamond_matrix(self, d):
        r"""
        EXAMPLES::

            sage: CuspForms(GammaH(31, [7]), 1).diamond_bracket_matrix(3)
            [-1]

            sage: C = CuspForms(GammaH(124, [33]), 1)   # long time
            sage: D = C.diamond_bracket_matrix(3); D    # long time
            [ 0  0  0 -1  1  0  0]
            [ 0 -1  0  0  0  0  0]
            [ 2  1  1 -2 -1 -2 -1]
            [ 0  0  0 -1  0  0  0]
            [-1  0  0  1  1  0  0]
            [ 2  0  0 -2 -1 -1  0]
            [ 0  2  1  0  0 -1  0]
            sage: t = C._transformation_matrix(); ~t * D * t   # long time
            [ 1  1  0  0  0  0  0]
            [-1  0  0  0  0  0  0]
            [ 0  0  1  1  0  0  0]
            [ 0  0 -1  0  0  0  0]
            [ 0  0  0  0 -1  0  0]
            [ 0  0  0  0  0 -1  0]
            [ 0  0  0  0  0  0 -1]
        """
        chars=self.group().characters_mod_H(sign=-1, galois_orbits=True)
        A = Matrix(QQ, 0, 0)
        for c in chars:
            chi = c.minimize_base_ring()
            dim = weight1.dimension_wt1_cusp_forms(chi)
            if chi.base_ring() == QQ:
                m = Matrix(QQ, 1, 1, [chi(d)])
            else:
                m = chi(d).matrix().transpose()
            for i in range(dim):
                A = A.block_sum(m)
        t = self._transformation_matrix()
        return t * A * ~t

    def _compute_hecke_matrix(self, n):
        r"""
        EXAMPLES::

            sage: CuspForms(GammaH(31, [7]), 1).hecke_matrix(7)
            [-1]
            sage: C = CuspForms(GammaH(124, [33]), 1) # long time
            sage: C.hecke_matrix(2) # long time
            [ 0  0 -1 -1  0  1  0]
            [ 1  0  0 -1 -1 -1  0]
            [ 0  0  0 -1  1  1 -1]
            [ 0  1  0 -1  0  0  0]
            [ 0  0 -1  0  0  1  1]
            [ 0  0  0 -1  0  0 -1]
            [ 0  0  0  0  0  1  0]
            sage: C.hecke_matrix(7) # long time
            [ 0  1  0 -1  0  0  1]
            [ 0 -1  0  0  0  0  0]
            [ 0  1 -1  0  0  0  1]
            [ 0  0  0 -1  0  0  0]
            [ 0  1  1  0  0 -1  0]
            [ 1  0 -1 -1 -1  0  1]
            [ 0  1  0  0  1  0  0]
            sage: C.hecke_matrix(23) == 0 # long time
            True

        """
        chars=self.group().characters_mod_H(sign=-1, galois_orbits=True)
        A = Matrix(QQ, 0, 0)
        for c in chars:
            chi = c.minimize_base_ring()
            d = weight1.dimension_wt1_cusp_forms(chi)
            e = chi.base_ring().degree()
            H = Matrix(QQ, d*e, d*e)
            from .constructor import CuspForms
            M = CuspForms(chi, 1).hecke_matrix(n)
            if e == 1:
                H = M
            else:
                for i in range(d):
                    for j in range(d):
                        H[e*i:e*(i+1), e*j:e*(j+1)] = M[i,j].matrix().transpose()
            A = A.block_sum(H)
        t = self._transformation_matrix()
        return t * A * ~t


class CuspidalSubmodule_g0_Q(CuspidalSubmodule_modsym_qexp):
    r"""
    Space of cusp forms for `\Gamma_0(N)` over `\QQ`.
    """
    def _pari_init_(self):
        """
        Conversion to Pari.

        EXAMPLES::

            sage: h = Newforms(37)[1]
            sage: MF = h.parent()
            sage: pari.mfparams(MF)
            [37, 2, 1, 1, t - 1]
            sage: pari.mfdim(MF)
            2
        """
        from sage.libs.pari import pari
        return pari.mfinit([self.level(), self.weight()], 1)


class CuspidalSubmodule_gH_Q(CuspidalSubmodule_modsym_qexp):
    r"""
    Space of cusp forms for `\Gamma_H(N)` over `\QQ`.
    """

    def _compute_hecke_matrix(self, n):
        r"""
        Compute the matrix of the Hecke operator T_n acting on this space.
        This is done directly using modular symbols, rather than using
        q-expansions as for spaces with fixed character.

        EXAMPLES::

            sage: CuspForms(Gamma1(8), 4)._compute_hecke_matrix(2)
            [  0 -16  32]
            [  1 -10  18]
            [  0  -4   8]

            sage: CuspForms(GammaH(15, [4]), 3)._compute_hecke_matrix(17)
            [ 18  22 -30 -60]
            [  4   0   6 -18]
            [  6   0   2 -20]
            [  6  12 -24 -20]
        """
        symbs = self.modular_symbols(sign=1)
        return _convert_matrix_from_modsyms(symbs, symbs.hecke_matrix(n))[0]

    def _compute_diamond_matrix(self, d):
        r"""
        EXAMPLES::

            sage: CuspForms(Gamma1(5), 6).diamond_bracket_matrix(3) # indirect doctest
            [ -1   0   0]
            [  3   5 -12]
            [  1   2  -5]
            sage: CuspForms(GammaH(15, [4]), 3)._compute_diamond_matrix(7)
            [ 2  3 -9 -3]
            [ 0  2 -3  0]
            [ 0  1 -2  0]
            [ 1  1 -3 -2]
        """
        symbs = self.modular_symbols(sign=1)
        return _convert_matrix_from_modsyms(symbs, symbs.diamond_bracket_matrix(d))[0]

class CuspidalSubmodule_g1_Q(CuspidalSubmodule_gH_Q):
    r"""
    Space of cusp forms for `\Gamma_1(N)` over `\QQ`.
    """

class CuspidalSubmodule_eps(CuspidalSubmodule_modsym_qexp):
    """
    Space of cusp forms with given Dirichlet character.

    EXAMPLES::

        sage: S = CuspForms(DirichletGroup(5).0,5); S
        Cuspidal subspace of dimension 1 of Modular Forms space of dimension 3, character [zeta4] and weight 5 over Cyclotomic Field of order 4 and degree 2

        sage: S.basis()
        [
        q + (-zeta4 - 1)*q^2 + (6*zeta4 - 6)*q^3 - 14*zeta4*q^4 + (15*zeta4 + 20)*q^5 + O(q^6)
        ]
        sage: f = S.0
        sage: f.qexp()
        q + (-zeta4 - 1)*q^2 + (6*zeta4 - 6)*q^3 - 14*zeta4*q^4 + (15*zeta4 + 20)*q^5 + O(q^6)
        sage: f.qexp(7)
        q + (-zeta4 - 1)*q^2 + (6*zeta4 - 6)*q^3 - 14*zeta4*q^4 + (15*zeta4 + 20)*q^5 + 12*q^6 + O(q^7)
        sage: f.qexp(3)
        q + (-zeta4 - 1)*q^2 + O(q^3)
        sage: f.qexp(2)
        q + O(q^2)
        sage: f.qexp(1)
        O(q^1)
    """
    pass

def _convert_matrix_from_modsyms(symbs, T):
    r"""
    Given a space of modular symbols and a matrix T acting on it, calculate the
    matrix of the corresponding operator on the echelon-form basis of the
    corresponding space of modular forms.

    The matrix T *must* commute with the Hecke operators! We use this when T is
    either a Hecke operator, or a diamond operator. This will *not work* for
    the Atkin-Lehner operators, for instance, when there are oldforms present.

    OUTPUT:
        A pair `(T_e, ps)` with `T_e` the converted matrix and `ps` a list
        of pivot elements of the echelon basis.

    EXAMPLES::

        sage: CuspForms(Gamma1(5), 6).diamond_bracket_matrix(3) # indirect doctest
        [ -1   0   0]
        [  3   5 -12]
        [  1   2  -5]
    """
    d = symbs.rank()

    # create a vector space of appropriate dimension to
    # contain our q-expansions
    A = symbs.base_ring()
    r = symbs.sturm_bound()
    X = A**r
    Y = X.zero_submodule()
    basis = []
    basis_images = []

    # we repeatedly use these matrices below, so we store them
    # once as lists to save time.
    hecke_matrix_ls = [ symbs.hecke_matrix(m).list() for m in range(1,r+1) ]
    hecke_image_ls = [ (T*symbs.hecke_matrix(m)).list() for m in range(1,r+1) ]

    # compute the q-expansions of some cusp forms and their
    # images under T_n
    for i in range(d**2):
        v = X([ hecke_matrix_ls[m][i] for m in range(r) ])
        Ynew = Y.span(Y.basis() + [v])
        if Ynew.rank() > Y.rank():
            basis.append(v)
            basis_images.append(X([ hecke_image_ls[m][i] for m in range(r) ]))
            Y = Ynew
            if len(basis) == d:
                break

    # now we can compute the matrix acting on the echelonized space of mod forms
    # need to pass A as base ring since otherwise there are problems when the
    # space has dimension 0
    bigmat = Matrix(A, basis).augment(Matrix(A, basis_images))
    bigmat.echelonize()
    pivs = bigmat.pivots()
    return bigmat.matrix_from_rows_and_columns(list(range(d)),
                                               [r + x for x in pivs]), pivs
