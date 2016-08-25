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
from __future__ import absolute_import

#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

from sage.rings.all import Integer
from sage.misc.all import verbose
from sage.matrix.all import Matrix

from .submodule import ModularFormsSubmodule
from . import vm_basis

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

        EXAMPLE::

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
        try:
            return self.__modular_symbols[sign]
        except AttributeError:
            self.__modular_symbols = {}
        except KeyError:
            pass
        A = self.ambient_module()
        S = A.modular_symbols(sign).cuspidal_submodule()
        self.__modular_symbols[sign] = S
        return S


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
        EXAMPLE::

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

    def new_submodule(self, p=None):
        r"""
        Return the new subspace of this space of cusp forms. This is computed
        using modular symbols.

        EXAMPLE::

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

class CuspidalSubmodule_g0_Q(CuspidalSubmodule_modsym_qexp):
    r"""
    Space of cusp forms for `\Gamma_0(N)` over `\QQ`.
    """

class CuspidalSubmodule_gH_Q(CuspidalSubmodule_modsym_qexp):
    r"""
    Space of cusp forms for `\Gamma_1(N)` over `\QQ`.
    """

    def _compute_hecke_matrix(self, n):
        r"""
        Compute the matrix of the Hecke operator T_n acting on this space.
        This is done directly using modular symbols, rather than using
        q-expansions as for spaces with fixed character.

        EXAMPLE::

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
        EXAMPLE::

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

    #def _repr_(self):
    #    A = self.ambient_module()
    #    return "Cuspidal subspace of dimension %s of Modular Forms space with character %s and weight %s over %s"%(self.dimension(), self.character(), self.weight(), self.base_ring())


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

    EXAMPLE::

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
    for i in xrange(d**2):
        v = X([ hecke_matrix_ls[m][i] for m in xrange(r) ])
        Ynew = Y.span(Y.basis() + [v])
        if Ynew.rank() > Y.rank():
            basis.append(v)
            basis_images.append(X([ hecke_image_ls[m][i] for m in xrange(r) ]))
            Y = Ynew
            if len(basis) == d: break

    # now we can compute the matrix acting on the echelonized space of mod forms
    # need to pass A as base ring since otherwise there are problems when the
    # space has dimension 0
    bigmat = Matrix(A, basis).augment(Matrix(A, basis_images))
    bigmat.echelonize()
    pivs = bigmat.pivots()
    return bigmat.matrix_from_rows_and_columns(range(d), [ r+x for x in pivs ]), pivs
