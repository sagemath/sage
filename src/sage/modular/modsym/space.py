# sage.doctest: needs sage.libs.flint sage.libs.pari
"""
Base class of the space of modular symbols

All the spaces of modular symbols derive from this class. This class is an
abstract base class.
"""

# ****************************************************************************
#       Sage: Open Source Mathematical Software
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

from sage.arith.misc import divisors, next_prime
from sage.categories.fields import Fields
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.modules.free_module import EchelonMatrixKey, FreeModule, VectorSpace
from sage.modules.free_module_element import FreeModuleElement
from sage.rings.fast_arith import prime_range
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.rings.infinity import infinity
from sage.rings.number_field.number_field_base import NumberField
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ
from sage.structure.all import Sequence, SageObject
from sage.structure.richcmp import (richcmp_method, richcmp,
                                    rich_to_bool, richcmp_not_equal)

from sage.modular.arithgroup.congroup_gamma0 import Gamma0_constructor as Gamma0, Gamma0_class # for Sturm bound given a character
from sage.modular.hecke.module import HeckeModule_free_module
from sage.modular.modsym.element import ModularSymbolsElement

lazy_import('sage.modular.modsym', 'hecke_operator')


def is_ModularSymbolsSpace(x):
    r"""
    Return ``True`` if ``x`` is a space of modular symbols.

    EXAMPLES::

        sage: M = ModularForms(3, 2)
        sage: sage.modular.modsym.space.is_ModularSymbolsSpace(M)
        doctest:warning...
        DeprecationWarning: The function is_ModularSymbolsSpace is deprecated; use 'isinstance(..., ModularForms)' instead.
        See https://github.com/sagemath/sage/issues/38035 for details.
        False
        sage: sage.modular.modsym.space.is_ModularSymbolsSpace(M.modular_symbols(sign=1))
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38035, "The function is_ModularSymbolsSpace is deprecated; use 'isinstance(..., ModularForms)' instead.")
    return isinstance(x, ModularSymbolsSpace)


@richcmp_method
class ModularSymbolsSpace(HeckeModule_free_module):
    r"""
    Base class for spaces of modular symbols.
    """

    Element = ModularSymbolsElement

    def __init__(self, group, weight, character, sign, base_ring, category=None):
        """
        Create a space of modular symbols.

        EXAMPLES::

            sage: M = ModularSymbols(22,6) ; M
            Modular Symbols space of dimension 30 for Gamma_0(22) of weight 6 with sign 0 over Rational Field
            sage: M == loads(dumps(M))
            True
        """
        self.__group = group
        self.__character = character
        self.__sign = sign
        HeckeModule_free_module.__init__(self, base_ring, group.level(), weight, category=category)

    def __richcmp__(self, other, op):
        """
        Compare ``self`` and ``other``.

        When spaces are in a common ambient space, we order
        lexicographically by the sequence of traces of Hecke operators
        `T_p`, for all primes `p`. In general we order
        first by the group, then the weight, then the character, then the
        sign then the base ring, then the dimension.

        EXAMPLES::

            sage: M = ModularSymbols(21,4) ; N = ModularSymbols(Gamma1(5),6)
            sage: M.cuspidal_submodule() > N
            True
            sage: M.cuspidal_submodule() == N
            False
        """
        if not isinstance(other, ModularSymbolsSpace):
            return NotImplemented

        lx = self.__group
        rx = other.__group
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx = self.weight()
        rx = other.weight()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx = self.__character
        rx = other.__character
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx = self.__sign
        rx = other.__sign
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx = self.base_ring()
        rx = other.base_ring()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        if self.is_ambient() or other.is_ambient():
            # if one is ambient they are equal iff they have same
            # dimension at this point, since all defining properties
            # are the same, so they are in the same ambient space.
            return richcmp(self.dimension(), other.dimension(), op)

        lx = self.ambient_hecke_module()
        rx = other.ambient_hecke_module()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx = self.dimension()
        rx = other.dimension()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx = EchelonMatrixKey(self.free_module())
        rx = EchelonMatrixKey(other.free_module())

        if lx == rx:
            return rich_to_bool(op, 0)

        # distinguish using Hecke operators, if possible.
        try:
            for p in prime_range(self.hecke_bound()):
                ap = self.hecke_matrix(p).trace()
                bp = other.hecke_matrix(p).trace()
                if ap != bp:
                    return richcmp_not_equal(ap, bp, op)
        except ArithmeticError:
            pass
        # fallback on subspace comparison
        return richcmp_not_equal(lx, rx, op)

    def _hecke_operator_class(self):
        """
        Return the class to be used for instantiating Hecke operators
        acting on ``self``.

        EXAMPLES::

            sage: ModularSymbols(81,2)._hecke_operator_class()
            <class 'sage.modular.modsym.hecke_operator.HeckeOperator'>
        """
        return hecke_operator.HeckeOperator

    def compact_system_of_eigenvalues(self, v, names='alpha', nz=None):
        r"""
        Return a compact system of eigenvalues `a_n` for
        `n\in v`.

        This should only be called on simple factors of modular
        symbols spaces.

        INPUT:

        - ``v`` -- list of positive integers
        - ``nz`` -- (default: ``None``) if given specifies a column index
          such that the dual module has that column nonzero

        OUTPUT:

        - ``E`` -- matrix such that E\*v is a vector with components
          the eigenvalues `a_n` for `n \in v`
        - ``v`` -- a vector over a number field

        EXAMPLES::

            sage: M = ModularSymbols(43,2,1)[2]; M
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1 over Rational Field
            sage: E, v = M.compact_system_of_eigenvalues(prime_range(10))
            sage: E
            [ 2/3 -4/3]
            [-2/3  4/3]
            [ 4/3  4/3]
            [-4/3 -4/3]
            sage: v
            (1, -3/4*alpha + 1/2)
            sage: E*v
            (alpha, -alpha, -alpha + 2, alpha - 2)

        TESTS:

        Verify that :issue:`12772` is fixed::

            sage: M = ModularSymbols(1,12,sign=1).cuspidal_subspace().new_subspace()
            sage: A = M.decomposition()[0]
            sage: A.compact_system_of_eigenvalues(prime_range(10))
            (
            [   -24]
            [   252]
            [  4830]
            [-16744], (1)
            )
        """
        if nz is None:
            nz = self._eigen_nonzero()
        M = self.ambient()
        try:
            E = M.hecke_images(nz, v) * self.dual_free_module().basis_matrix().transpose()
        except AttributeError:
            # TODO!!!
            raise NotImplementedError("ambient space must implement hecke_images but doesn't yet")
        v = self.dual_eigenvector(names=names, lift=False, nz=nz)
        return E, v

    def character(self):
        """
        Return the character associated to ``self``.

        EXAMPLES::

            sage: ModularSymbols(12,8).character()
            Dirichlet character modulo 12 of conductor 1 mapping 7 |--> 1, 5 |--> 1
            sage: ModularSymbols(DirichletGroup(25).0, 4).character()
            Dirichlet character modulo 25 of conductor 25 mapping 2 |--> zeta20
        """
        return self.__character

    def cuspidal_submodule(self):
        """
        Return the cuspidal submodule of ``self``.

        .. NOTE::

            This should be overridden by all derived classes.

        EXAMPLES::

            sage: sage.modular.modsym.space.ModularSymbolsSpace(Gamma0(11),2,DirichletGroup(11).gens()[0]**10,0,QQ).cuspidal_submodule()
            Traceback (most recent call last):
            ...
            NotImplementedError: computation of cuspidal submodule not yet implemented for this class
            sage: ModularSymbols(Gamma0(11),2).cuspidal_submodule()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
        """
        raise NotImplementedError("computation of cuspidal submodule not yet implemented for this class")

    def cuspidal_subspace(self):
        """
        Synonym for cuspidal_submodule.

        EXAMPLES::

            sage: m = ModularSymbols(Gamma1(3),12); m.dimension()
            8
            sage: m.cuspidal_subspace().new_subspace().dimension()
            2
        """
        return self.cuspidal_submodule()

    def new_subspace(self, p=None):
        """
        Synonym for new_submodule.

        EXAMPLES::

            sage: m = ModularSymbols(Gamma0(5),12); m.dimension()
            12
            sage: m.new_subspace().dimension()
            6
            sage: m = ModularSymbols(Gamma1(3),12); m.dimension()
            8
            sage: m.new_subspace().dimension()
            2
        """
        return self.new_submodule(p)

    def old_subspace(self, p=None):
        """
        Synonym for old_submodule.

        EXAMPLES::

            sage: m = ModularSymbols(Gamma1(3),12); m.dimension()
            8
            sage: m.old_subspace().dimension()
            6
        """
        return self.old_submodule(p)

    def eisenstein_subspace(self):
        """
        Synonym for eisenstein_submodule.

        EXAMPLES::

            sage: m = ModularSymbols(Gamma1(3),12); m.dimension()
            8
            sage: m.eisenstein_subspace().dimension()
            2
            sage: m.cuspidal_subspace().dimension()
            6
        """
        return self.eisenstein_submodule()

    def dimension_of_associated_cuspform_space(self):
        """
        Return the dimension of the corresponding space of cusp forms.

        The input space must be cuspidal, otherwise there is no
        corresponding space of cusp forms.

        EXAMPLES::

            sage: m = ModularSymbols(Gamma0(389),2).cuspidal_subspace(); m.dimension()
            64
            sage: m.dimension_of_associated_cuspform_space()
            32
            sage: m = ModularSymbols(Gamma0(389),2,sign=1).cuspidal_subspace(); m.dimension()
            32
            sage: m.dimension_of_associated_cuspform_space()
            32
        """
        if not self.is_cuspidal():
            raise ArithmeticError("space must be cuspidal")
        if self.sign() == 0:
            return self.dimension() // 2
        return self.dimension()

    def dual_star_involution_matrix(self):
        """
        Return the matrix of the dual star involution, which is induced by
        complex conjugation on the linear dual of modular symbols.

        .. NOTE::

            This should be overridden in all derived classes.

        EXAMPLES::

            sage: sage.modular.modsym.space.ModularSymbolsSpace(Gamma0(11),2,DirichletGroup(11).gens()[0]**10,0,QQ).dual_star_involution_matrix()
            Traceback (most recent call last):
            ...
            NotImplementedError: computation of dual star involution matrix not yet implemented for this class
            sage: ModularSymbols(Gamma0(11),2).dual_star_involution_matrix()
            [ 1  0  0]
            [ 0 -1  0]
            [ 0  1  1]
        """
        raise NotImplementedError("computation of dual star involution matrix not yet implemented for this class")

    def group(self):
        """
        Return the group of this modular symbols space.

        INPUT:

        - ``ModularSymbols self`` -- an arbitrary space of
          modular symbols

        OUTPUT:

        - ``CongruenceSubgroup`` -- the congruence subgroup
          that this is a space of modular symbols for

        ALGORITHM: The group is recorded when this space is created.

        EXAMPLES::

            sage: m = ModularSymbols(20)
            sage: m.group()
            Congruence Subgroup Gamma0(20)
        """
        return self.__group

    def is_ambient(self):
        """
        Return ``True`` if ``self`` is an ambient space of modular symbols.

        EXAMPLES::

            sage: ModularSymbols(21,4).is_ambient()
            True
            sage: ModularSymbols(21,4).cuspidal_submodule().is_ambient()
            False
        """
        from sage.modular.modsym.ambient import ModularSymbolsAmbient
        return isinstance(self, ModularSymbolsAmbient)

    def is_cuspidal(self):
        """
        Return ``True`` if ``self`` is a cuspidal space of modular symbols.

        .. NOTE::

            This should be overridden in all derived classes.

        EXAMPLES::

            sage: sage.modular.modsym.space.ModularSymbolsSpace(Gamma0(11),2,DirichletGroup(11).gens()[0]**10,0,QQ).is_cuspidal()
            Traceback (most recent call last):
            ...
            NotImplementedError: computation of cuspidal subspace not yet implemented for this class
            sage: ModularSymbols(Gamma0(11),2).is_cuspidal()
            False
        """
        raise NotImplementedError("computation of cuspidal subspace not yet implemented for this class")

    def is_simple(self):
        r"""
        Return whether this modular symbols space is simple as a module
        over the anemic Hecke algebra adjoin \*.

        EXAMPLES::

            sage: m = ModularSymbols(Gamma0(33),2,sign=1)
            sage: m.is_simple()
            False
            sage: o = m.old_subspace()
            sage: o.decomposition()
            [Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 6 for Gamma_0(33) of weight 2 with sign 1 over Rational Field,
             Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 6 for Gamma_0(33) of weight 2 with sign 1 over Rational Field]
            sage: C = ModularSymbols(1,14,0,GF(5)).cuspidal_submodule(); C
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 for Gamma_0(1) of weight 14 with sign 0 over Finite Field of size 5
            sage: C.is_simple()
            True
        """
        try:
            return self._is_simple
        except AttributeError:
            D = self.factorization()
            if len(D) == 0 or len(D) == 1 and D[0][1] == 1:
                self._is_simple = True
            else:
                self._is_simple = False
            return self._is_simple

    def multiplicity(self, S, check_simple=True):
        """
        Return the multiplicity of the simple modular symbols space S in
        ``self``. S must be a simple anemic Hecke module.

        ASSUMPTION: ``self`` is an anemic Hecke module with the same weight and
        group as S, and S is simple.

        EXAMPLES::

            sage: M = ModularSymbols(11,2,sign=1)
            sage: N1, N2 = M.decomposition()
            sage: N1.multiplicity(N2)
            0
            sage: M.multiplicity(N1)
            1
            sage: M.multiplicity(ModularSymbols(14,2))
            0
        """
        if self.level() % S.level() != 0 or S.weight() != self.weight():
            return 0
        if check_simple and not S.is_simple():
            raise ArithmeticError("S must be simple")
        A = self.ambient_hecke_module()
        B = A.submodule_generated_by_images(S)
        C = self.intersection(B)
        d = C.rank()
        n = S.rank()
        assert d % n == 0, "the dimension of intersection must be a multiple of dimension of simple space.  bug!"
        return d // n

    def ngens(self):
        """
        Return the number of generators of ``self``.

        INPUT:

        - ``ModularSymbols self`` -- arbitrary space of modular symbols

        OUTPUT:

        - ``int`` -- the number of generators, which is the same as the
          dimension of ``self``

        ALGORITHM: Call the dimension function.

        EXAMPLES::

            sage: m = ModularSymbols(33)
            sage: m.ngens()
            9
            sage: m.rank()
            9
            sage: ModularSymbols(100, weight=2, sign=1).ngens()
            18
        """
        return self.rank()

    #########################################################################
    #
    #  Computation of q-expansions
    #
    #########################################################################

    def default_prec(self):
        r"""
        Get the default precision for computation of `q`-expansion
        associated to the ambient space of this space of modular symbols
        (and all subspaces). Use ``set_default_prec`` to
        change the default precision.

        EXAMPLES::

            sage: M = ModularSymbols(15)
            sage: M.cuspidal_submodule().q_expansion_basis()
            [q - q^2 - q^3 - q^4 + q^5 + q^6 + O(q^8)]
            sage: M.set_default_prec(20)

        Notice that setting the default precision of the ambient space
        affects the subspaces.

        ::

            sage: M.cuspidal_submodule().q_expansion_basis()
            [q - q^2 - q^3 - q^4 + q^5 + q^6 + 3*q^8 + q^9 - q^10 - 4*q^11 + q^12 - 2*q^13 - q^15 - q^16 + 2*q^17 - q^18 + 4*q^19 + O(q^20)]
            sage: M.cuspidal_submodule().default_prec()
            20
        """
        if not self.is_ambient():
            return self.ambient_hecke_module().default_prec()
        try:
            return self.__default_prec
        except AttributeError:
            self.__default_prec = Integer(8)
            return self.__default_prec

    def set_default_prec(self, prec):
        """
        Set the default precision for computation of `q`-expansion
        associated to the ambient space of this space of modular symbols
        (and all subspaces).

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(13),2)
            sage: M.set_default_prec(5)
            sage: M.cuspidal_submodule().q_expansion_basis()
            [q - 4*q^3 - q^4 + O(q^5), q^2 - 2*q^3 - q^4 + O(q^5)]
        """
        if not self.is_ambient():
            return self.ambient_hecke_module().set_default_prec(prec)
        else:
            self.__default_prec = Integer(prec)

    def set_precision(self, prec):
        """
        Same as self.set_default_prec(prec).

        EXAMPLES::

            sage: M = ModularSymbols(17,2)
            sage: M.cuspidal_submodule().q_expansion_basis()
            [q - q^2 - q^4 - 2*q^5 + 4*q^7 + O(q^8)]
            sage: M.set_precision(10)
            sage: M.cuspidal_submodule().q_expansion_basis()
            [q - q^2 - q^4 - 2*q^5 + 4*q^7 + 3*q^8 - 3*q^9 + O(q^10)]
        """
        self.set_default_prec(prec)

    def q_expansion_basis(self, prec=None, algorithm='default'):
        r"""
        Return a basis of `q`-expansions (as power series) to precision ``prec``
        of the space of modular forms associated to ``self``.

        The `q`-expansions are defined over the same base ring as ``self``,
        and a put in echelon form.

        INPUT:

        - ``self`` -- a space of CUSPIDAL modular symbols

        - ``prec`` -- integer

        - ``algorithm`` -- string; one of

          - ``'default'`` -- (default) decide which algorithm to
            use based on heuristics

          - ``'hecke'`` -- compute basis by computing
            homomorphisms T - K, where T is the Hecke algebra

          - ``'eigen'`` -- compute basis using eigenvectors for
            the Hecke action and Atkin-Lehner-Li theory to patch them together

          - ``'all'`` -- compute using hecke_dual and eigen
            algorithms and verify that the results are the same


        The computed basis is *not* cached, though of course Hecke
        operators used in computing the basis are cached.

        EXAMPLES::

            sage: M = ModularSymbols(1, 12).cuspidal_submodule()
            sage: M.q_expansion_basis(8)
            [q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 - 16744*q^7 + O(q^8)]

        ::

            sage: M.q_expansion_basis(8, algorithm='eigen')
            [q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 - 16744*q^7 + O(q^8)]

        ::

            sage: M = ModularSymbols(1, 24).cuspidal_submodule()
            sage: M.q_expansion_basis(8, algorithm='eigen')
            [q + 195660*q^3 + 12080128*q^4 + 44656110*q^5 - 982499328*q^6 - 147247240*q^7 + O(q^8),
             q^2 - 48*q^3 + 1080*q^4 - 15040*q^5 + 143820*q^6 - 985824*q^7 + O(q^8)]

        ::

            sage: M = ModularSymbols(11, 2, sign=-1).cuspidal_submodule()
            sage: M.q_expansion_basis(8, algorithm='eigen')
            [q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 + O(q^8)]

        ::

            sage: M = ModularSymbols(Gamma1(13), 2, sign=1).cuspidal_submodule()
            sage: M.q_expansion_basis(8, algorithm='eigen')
            [q - 4*q^3 - q^4 + 3*q^5 + 6*q^6 + O(q^8),
             q^2 - 2*q^3 - q^4 + 2*q^5 + 2*q^6 + O(q^8)]

        ::

            sage: M = ModularSymbols(Gamma1(5), 3, sign=-1).cuspidal_submodule()
            sage: M.q_expansion_basis(8, algorithm='eigen')   # dimension is 0
            []

        ::

            sage: M = ModularSymbols(Gamma1(7), 3, sign=-1).cuspidal_submodule()
            sage: M.q_expansion_basis(8)
            [q - 3*q^2 + 5*q^4 - 7*q^7 + O(q^8)]

        ::

            sage: M = ModularSymbols(43, 2, sign=0).cuspidal_submodule()
            sage: M[0]
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 7 for Gamma_0(43) of weight 2 with sign 0 over Rational Field
            sage: M[0].q_expansion_basis()
            [q - 2*q^2 - 2*q^3 + 2*q^4 - 4*q^5 + 4*q^6 + O(q^8)]
            sage: M[1]
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 7 for Gamma_0(43) of weight 2 with sign 0 over Rational Field
            sage: M[1].q_expansion_basis()
            [q + 2*q^5 - 2*q^6 - 2*q^7 + O(q^8), q^2 - q^3 - q^5 + q^7 + O(q^8)]
        """
        if prec is None:
            prec = self.default_prec()
        else:
            prec = Integer(prec)

        if prec < 1:
            raise ValueError("prec (=%s) must be >= 1" % prec)

        if not self.is_cuspidal():
            raise ArithmeticError("space must be cuspidal")

        if self.sign() == 0:
            P = self.plus_submodule(compute_dual=True)
            return Sequence(P.q_expansion_basis(prec=prec, algorithm=algorithm), cr=True)

        if self.dimension() == 0:
            return Sequence([])

        if algorithm == 'default':
            algorithm = 'hecke'
        if algorithm == 'hecke':
            return Sequence(self._q_expansion_basis_hecke_dual(prec), cr=True)
        elif algorithm == 'eigen':
            return Sequence(self._q_expansion_basis_eigen(prec, 'alpha'), cr=True)
        elif algorithm == 'all':
            B1 = self._q_expansion_basis_hecke_dual(prec)
            B2 = self._q_expansion_basis_eigen(prec, 'alpha')
            if B1 != B2:
                raise RuntimeError("There is a bug in q_expansion_basis -- basis computed differently with two algorithms:\n%s\n%s\n" % (B1, B2,))
            return Sequence(B1, cr=True)
        else:
            raise ValueError("no algorithm '%s'" % algorithm)

    def q_expansion_module(self, prec=None, R=None):
        r"""
        Return a basis over R for the space spanned by the coefficient
        vectors of the `q`-expansions corresponding to ``self``.

        If R is not the base ring of ``self``, this returns the
        restriction of scalars down to R (for this, ``self`` must have
        base ring `\QQ` or a number field).

        INPUT:

        - ``self`` -- must be cuspidal

        - ``prec`` -- integer (default: ``self.default_prec()``)

        - ``R`` -- either `\ZZ`, `\QQ`, or the ``base_ring`` of ``self``
           (which is the default)

        OUTPUT: a free module over `R`

        .. TODO::

           extend to more general R (though that is fairly easy for the
           user to get by just doing base_extend or change_ring on the
           output of this function).

        Note that the prec needed to distinguish elements of the
        restricted-down-to-R basis may be bigger than ``self.hecke_bound()``,
        since one must use the Sturm bound for modular forms on `\Gamma_H(N)`.

        EXAMPLES WITH SIGN 1 and R=QQ:

        Basic example with sign 1::

            sage: M = ModularSymbols(11, sign=1).cuspidal_submodule()
            sage: M.q_expansion_module(5, QQ)
            Vector space of degree 5 and dimension 1 over Rational Field
            Basis matrix:
            [ 0  1 -2 -1  2]

        Same example with sign -1::

            sage: M = ModularSymbols(11, sign=-1).cuspidal_submodule()
            sage: M.q_expansion_module(5, QQ)
            Vector space of degree 5 and dimension 1 over Rational Field
            Basis matrix:
            [ 0  1 -2 -1  2]

        An example involving old forms::

            sage: M = ModularSymbols(22, sign=1).cuspidal_submodule()
            sage: M.q_expansion_module(5, QQ)
            Vector space of degree 5 and dimension 2 over Rational Field
            Basis matrix:
            [ 0  1  0 -1 -2]
            [ 0  0  1  0 -2]

        An example that (somewhat spuriously) is over a number field::

            sage: x = polygen(QQ)
            sage: k = NumberField(x^2+1, 'a')
            sage: M = ModularSymbols(11, base_ring=k, sign=1).cuspidal_submodule()
            sage: M.q_expansion_module(5, QQ)
            Vector space of degree 5 and dimension 1 over Rational Field
            Basis matrix:
            [ 0  1 -2 -1  2]

        An example that involves an eigenform with coefficients in a number
        field::

            sage: M = ModularSymbols(23, sign=1).cuspidal_submodule()
            sage: M.q_eigenform(4, 'gamma')
            q + gamma*q^2 + (-2*gamma - 1)*q^3 + O(q^4)
            sage: M.q_expansion_module(11, QQ)
            Vector space of degree 11 and dimension 2 over Rational Field
            Basis matrix:
            [ 0  1  0 -1 -1  0 -2  2 -1  2  2]
            [ 0  0  1 -2 -1  2  1  2 -2  0 -2]

        An example that is genuinely over a base field besides QQ.

        ::

            sage: eps = DirichletGroup(11).0
            sage: M = ModularSymbols(eps,3,sign=1).cuspidal_submodule(); M
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 and level 11, weight 3, character [zeta10], sign 1, over Cyclotomic Field of order 10 and degree 4
            sage: M.q_eigenform(4, 'beta')
            q + (-zeta10^3 + 2*zeta10^2 - 2*zeta10)*q^2 + (2*zeta10^3 - 3*zeta10^2 + 3*zeta10 - 2)*q^3 + O(q^4)
            sage: M.q_expansion_module(7, QQ)
            Vector space of degree 7 and dimension 4 over Rational Field
            Basis matrix:
            [  0   1   0   0   0 -40  64]
            [  0   0   1   0   0 -24  41]
            [  0   0   0   1   0 -12  21]
            [  0   0   0   0   1  -4   4]

        An example involving an eigenform rational over the base, but the
        base is not QQ.

        ::

            sage: k.<a> = NumberField(x^2-5)
            sage: M = ModularSymbols(23, base_ring=k, sign=1).cuspidal_submodule()
            sage: D = M.decomposition(); D
            [Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(23) of weight 2 with sign 1 over Number Field in a with defining polynomial x^2 - 5,
             Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(23) of weight 2 with sign 1 over Number Field in a with defining polynomial x^2 - 5]
            sage: M.q_expansion_module(8, QQ)
            Vector space of degree 8 and dimension 2 over Rational Field
            Basis matrix:
            [ 0  1  0 -1 -1  0 -2  2]
            [ 0  0  1 -2 -1  2  1  2]

        An example involving an eigenform not rational over the base and
        for which the base is not QQ.

        ::

            sage: eps = DirichletGroup(25).0^2
            sage: M = ModularSymbols(eps,2,sign=1).cuspidal_submodule(); M
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 and level 25, weight 2, character [zeta10], sign 1, over Cyclotomic Field of order 10 and degree 4
            sage: D = M.decomposition(); D
            [Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 and level 25, weight 2, character [zeta10], sign 1, over Cyclotomic Field of order 10 and degree 4]
            sage: D[0].q_eigenform(4, 'mu')
            q + mu*q^2 + ((zeta10^3 + zeta10 - 1)*mu + zeta10^2 - 1)*q^3 + O(q^4)
            sage: D[0].q_expansion_module(11, QQ)
            Vector space of degree 11 and dimension 8 over Rational Field
            Basis matrix:
            [  0   1   0   0   0   0   0   0 -20  -3   0]
            [  0   0   1   0   0   0   0   0 -16  -1   0]
            [  0   0   0   1   0   0   0   0 -11  -2   0]
            [  0   0   0   0   1   0   0   0  -8  -1   0]
            [  0   0   0   0   0   1   0   0  -5  -1   0]
            [  0   0   0   0   0   0   1   0  -3  -1   0]
            [  0   0   0   0   0   0   0   1  -2   0   0]
            [  0   0   0   0   0   0   0   0   0   0   1]
            sage: D[0].q_expansion_module(11)
            Vector space of degree 11 and dimension 2 over Cyclotomic Field of order 10 and degree 4
            Basis matrix:
            [                                  0                                   1                                   0                        zeta10^2 - 1                       -zeta10^2 - 1                -zeta10^3 - zeta10^2                   zeta10^2 - zeta10           2*zeta10^3 + 2*zeta10 - 1    zeta10^3 - zeta10^2 - zeta10 + 1        zeta10^3 - zeta10^2 + zeta10   -2*zeta10^3 + 2*zeta10^2 - zeta10]
            [                                  0                                   0                                   1               zeta10^3 + zeta10 - 1                         -zeta10 - 1                -zeta10^3 - zeta10^2 -2*zeta10^3 + zeta10^2 - zeta10 + 1                            zeta10^2                                   0                        zeta10^3 + 1  2*zeta10^3 - zeta10^2 + zeta10 - 1]

        EXAMPLES WITH SIGN 0 and R=QQ:

        .. TODO::  This doesn't work yet as it's not implemented!!

        ::

            sage: M = ModularSymbols(11,2).cuspidal_submodule() #not tested
            sage: M.q_expansion_module() #not tested
            ... boom ...

        EXAMPLES WITH SIGN 1 and R=ZZ (computes saturation)::

            sage: M = ModularSymbols(43,2, sign=1).cuspidal_submodule()
            sage: M.q_expansion_module(8, QQ)
            Vector space of degree 8 and dimension 3 over Rational Field
            Basis matrix:
            [   0    1    0    0    0    2   -2   -2]
            [   0    0    1    0 -1/2    1 -3/2    0]
            [   0    0    0    1 -1/2    2 -3/2   -1]
            sage: M.q_expansion_module(8, ZZ)
            Free module of degree 8 and rank 3 over Integer Ring
            Echelon basis matrix:
            [ 0  1  0  0  0  2 -2 -2]
            [ 0  0  1  1 -1  3 -3 -1]
            [ 0  0  0  2 -1  4 -3 -2]
        """
        if prec is None:
            prec = self.default_prec()
        if R == ZZ and self.base_ring() == QQ:
            return self._q_expansion_module_integral(prec)
        elif R == QQ:
            return self._q_expansion_module_rational(prec)
        elif R is None or R == self.base_ring():
            # names is never used in this case
            return self._q_expansion_module(prec)
        else:
            raise NotImplementedError("R must be ZZ, QQ, or the base ring of the modular symbols space.")

    def _q_eigenform_images(self, A, prec, names):
        """
        Return list of images in space corresponding to ``self`` of eigenform
        corresponding to A under the degeneracy maps.

        This is mainly a helper function for other internal functions.

        INPUT:

        - ``self`` -- space of modular symbols

        - ``A`` -- cuspidal simple space of level dividing the
          level of ``self`` and the same weight

        - ``prec`` -- positive integer

        EXAMPLES::

            sage: M = ModularSymbols(33,2,sign=1)
            sage: A = M.modular_symbols_of_level(11).cuspidal_submodule()
            sage: M._q_eigenform_images(A, 10, names='a')
            [q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 - 2*q^9 + O(q^10),
             q^3 - 2*q^6 - q^9 + 2*q^12 + q^15 + 2*q^18 - 2*q^21 - 2*q^27 + O(q^30)]
        """
        f = A.q_eigenform(prec, names)
        if A.level() == self.level():
            return [f]
        D = divisors(self.level() // A.level())
        q = f.parent().gen()
        return [f] + [f(q**d) for d in D if d > 1]

    def _q_expansion_module(self, prec, algorithm='hecke'):
        """
        Return module spanned by the `q`-expansions corresponding to ``self``.

        See the docstring for ``q_expansion_module`` (without
        underscore) for further details. Note that this will not work
        if ``algorithm=eigen`` and the sign is 0.

        EXAMPLES::

            sage: ModularSymbols(11, 2, base_ring=GF(4,'a')).cuspidal_submodule()._q_expansion_module(prec=4, algorithm='hecke')
            Vector space of degree 4 and dimension 1 over Finite Field in a of size 2^2
            Basis matrix:
            [0 1 0 1]
            sage: ModularSymbols(11, 2, base_ring=QuadraticField(-7,'b'), sign=1).cuspidal_submodule()._q_expansion_module(prec=4, algorithm='eigen')
            Vector space of degree 4 and dimension 1 over Number Field in b with defining polynomial x^2 + 7 with b = 2.645751311064591?*I
            Basis matrix:
            [ 0  1 -2 -1]
        """
        if not self.is_cuspidal():
            raise ValueError("self must be cuspidal")
        R = self.base_ring()
        if R not in Fields():
            if R == ZZ:
                return self._q_expansion_module_integral(prec)
            raise NotImplementedError("base ring must be a field (or ZZ).")

        V = R ** prec
        if algorithm == 'hecke' or algorithm == 'default':
            return V.span([f.padded_list(prec) for f in self.q_expansion_basis(prec, algorithm)])

        if algorithm != 'eigen':
            raise ValueError("unknown algorithm '%s'" % algorithm)

        def q_eigen_gens(d, f):
            r"""
            Temporary function for internal use.

            EXAMPLES::

                sage: ModularSymbols(11, 4, base_ring=QuadraticField(-7,'b'),sign=1).cuspidal_submodule()._q_expansion_module(prec=5, algorithm='eigen') # indirect doctest
                Vector space of degree 5 and dimension 2 over Number Field in b with defining polynomial x^2 + 7 with b = 2.645751311064591?*I
                Basis matrix:
                [ 0  1  0  3 -6]
                [ 0  0  1 -4  2]
            """
            X = f.padded_list(prec)
            if d == 1:
                # X is just a list of elements of R
                return [X]
            else:
                # X is a list of elements of a poly quotient ring
                return [[X[i][j] for i in range(prec)] for j in range(d)]

        if self.sign() == 0:
            X = self.plus_submodule(compute_dual=True)
        else:
            X = self

        B = [sum([q_eigen_gens(A.dimension(), f)
                  for f in self._q_eigenform_images(A, prec, 'zeta')], [])
             for A, _ in X.factorization()]

        return V.span(sum(B, []))

    def _q_expansion_module_rational(self, prec):
        r"""
        Return a vector space over `\QQ` for the space spanned by the
        `q`-expansions corresponding to ``self``.

        The base ring of ``self`` must be `\QQ` or a number field, and
        ``self`` must be cuspidal. The returned space is a
        `\QQ`-vector space, where the coordinates are the coefficients
        of `q`-expansions.

        INPUT:

        - ``prec`` -- integer; number of `q`-expansion terms to calculate

        EXAMPLES::

            sage: ModularSymbols(11, 4).cuspidal_submodule()._q_expansion_module_rational(5)
            Vector space of degree 5 and dimension 2 over Rational Field
            Basis matrix:
            [ 0  1  0  3 -6]
            [ 0  0  1 -4  2]
        """
        if not self.is_cuspidal():
            raise ValueError("self must be cuspidal")
        K = self.base_ring()
        if not isinstance(K, NumberField):
            raise TypeError("self must be over QQ or a number field.")
        n = K.degree()
        if n == 1:
            return self._q_expansion_module(prec)

        # Construct the vector space over QQ of dimension equal to
        # the degree of the base field times the dimension over C
        # of the space of cusp forms corresponding to self.
        def q_eigen_gens(d, f):
            r"""
            Temporary function for internal use.

            EXAMPLES::

                sage: ModularSymbols(13, 6).cuspidal_submodule()._q_expansion_module_rational(4) # indirect doctest
                Vector space of degree 4 and dimension 3 over Rational Field
                Basis matrix:
                [0 1 0 0]
                [0 0 1 0]
                [0 0 0 1]
            """
            # Return restricted down to QQ gens for cusp space corresponding
            # to the simple factor A.
            X = f.padded_list(prec)
            if d == 1:
                return [[X[i][j] for i in range(prec)] for j in range(n)]
            else:
                # This looks like it might be really slow -- though
                # perhaps it's nothing compared to the time taken by
                # whatever computed this in the first place.
                return [[(X[i].list())[j][k] for i in range(prec)]
                        for j in range(d) for k in range(n)]
        if self.sign() == 0:
            X = self.plus_submodule(compute_dual=True)
        else:
            X = self

        B = [sum([q_eigen_gens(A.dimension(), f)
                  for f in self._q_eigenform_images(A, prec, 'alpha')], [])
             for A, _ in X.factorization()]
        A = QQ**prec
        return A.span(sum(B, []))

    def _q_expansion_module_integral(self, prec):
        r"""
        Return module over `\ZZ` for the space spanned by
        the `q`-expansions corresponding to ``self``.

        The base ring of ``self`` must be `\QQ` or a number field, and
        ``self`` must be cuspidal. The returned space is a `\ZZ`-module,
        where the coordinates are the coefficients of `q`-expansions.

        EXAMPLES::

            sage: M = ModularSymbols(11, sign=1).cuspidal_submodule()
            sage: M._q_expansion_module_integral(5)
            Free module of degree 5 and rank 1 over Integer Ring
            Echelon basis matrix:
            [ 0  1 -2 -1  2]
            sage: V = M.complement().cuspidal_submodule()
            sage: V._q_expansion_module_integral(5)
            Free module of degree 5 and rank 0 over Integer Ring
            Echelon basis matrix:
            []
        """
        V = self.q_expansion_module(prec, QQ)
        return FreeModule(ZZ, V.degree()).span(V.basis()).saturation()

    def congruence_number(self, other, prec=None):
        r"""
        Given two cuspidal spaces of modular symbols, compute the
        congruence number, using ``prec`` terms of the `q`-expansions.

        The congruence number is defined as follows. If `V` is the
        submodule of integral cusp forms corresponding to ``self`` (saturated in
        `\ZZ[[q]]`, by definition) and `W` is the
        submodule corresponding to other, each computed to precision ``prec``,
        the congruence number is the index of `V+W` in its
        saturation in `\ZZ[[q]]`.

        If ``prec`` is not given it is set equal to the max of the
        ``hecke_bound`` function called on each space.

        EXAMPLES::

            sage: A, B = ModularSymbols(48, 2).cuspidal_submodule().decomposition()
            sage: A.congruence_number(B)
            2
        """
        if not self.is_cuspidal():
            raise ValueError("self must be cuspidal")
        if not other.is_cuspidal():
            raise ValueError("right must be cuspidal")
        if prec is None:
            prec = max(self.hecke_bound(), other.hecke_bound())
        prec = int(prec)

        V = self.q_expansion_module(prec, ZZ)
        W = other.q_expansion_module(prec, ZZ)
        K = V + W
        return K.index_in_saturation()

    #########################################################################
    #
    #  Computation of a basis using eigenforms
    #
    #########################################################################

    @cached_method
    def q_eigenform_character(self, names=None):
        """
        Return the Dirichlet character associated to the specific
        choice of `q`-eigenform attached to this simple cuspidal
        modular symbols space.

        INPUT:

        - ``names`` -- string; name of the variable

        OUTPUT:

        - a Dirichlet character taking values in the Hecke eigenvalue
          field, where the indeterminate of that field is determined
          by the given variable name.

        EXAMPLES::

            sage: f = ModularSymbols(Gamma1(13),2,sign=1).cuspidal_subspace().decomposition()[0]
            sage: eps = f.q_eigenform_character('a'); eps
            Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -a - 1
            sage: parent(eps)
            Group of Dirichlet characters modulo 13 with values in Number Field in a with defining polynomial x^2 + 3*x + 3
            sage: eps(3)
            a + 1

        The modular symbols space must be simple.::

            sage: ModularSymbols(Gamma1(17),2,sign=1).cuspidal_submodule().q_eigenform_character('a')
            Traceback (most recent call last):
            ...
            ArithmeticError: self must be simple

        If the character is specified when making the modular symbols
        space, then names need not be given and the returned character
        is just the character of the space.::

            sage: f = ModularSymbols(kronecker_character(19),2,sign=1).cuspidal_subspace().decomposition()[0]
            sage: f
            Modular Symbols subspace of dimension 8 of Modular Symbols space of dimension 10 and level 76, weight 2, character [-1, -1], sign 1, over Rational Field
            sage: f.q_eigenform_character()
            Dirichlet character modulo 76 of conductor 76 mapping 39 |--> -1, 21 |--> -1
            sage: f.q_eigenform_character() is f.character()
            True

        The input space need not be cuspidal::

            sage: M = ModularSymbols(Gamma1(13),2,sign=1).eisenstein_submodule()[0]
            sage: M.q_eigenform_character('a')
            Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -1

        The modular symbols space does not have to come from a decomposition::

            sage: ModularSymbols(Gamma1(16),2,sign=1).cuspidal_submodule().q_eigenform_character('a')
            Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> -a - 1
        """
        eps = self.character()
        if eps is not None:
            # easy case
            return eps

        v = self.dual_eigenvector(names=names)
        i = v.nonzero_positions()[0]
        K = v.base_ring()
        from sage.modular.dirichlet import DirichletGroup
        G = DirichletGroup(self.level(), K)
        M = self.ambient_module()
        # act on right since v is a in the dual
        b = [(M.diamond_bracket_matrix(u) * v)[i] / v[i]
             for u in G.unit_gens()]
        return G(b)

    def q_eigenform(self, prec, names=None):
        """
        Return the `q`-expansion to precision ``prec`` of a new eigenform
        associated to ``self``.

        Here ``self`` must be new, cuspidal, and simple.

        EXAMPLES::

            sage: ModularSymbols(2, 8)[1].q_eigenform(5, 'a')
            q - 8*q^2 + 12*q^3 + 64*q^4 + O(q^5)
            sage: ModularSymbols(2, 8)[0].q_eigenform(5,'a')
            Traceback (most recent call last):
            ...
            ArithmeticError: self must be cuspidal.
        """
        if self.dimension() > 1 and names is None:
            raise ValueError("please specify a name to use for the field of eigenvalues")

        if prec is None:
            prec = self.default_prec()
        try:
            f = self._q_expansion_dict[names]
        except (AttributeError, KeyError):
            self._q_expansion_dict = {}
            if not self.is_cuspidal():
                raise ArithmeticError("self must be cuspidal.")

            if not self.is_simple():
                if self.sign() == 0:
                    return self.plus_submodule(compute_dual=True).q_eigenform(prec, names)
                raise ArithmeticError("self must be simple.")
            a2 = self.eigenvalue(2, names)
            R = PowerSeriesRing(a2.parent(), "q")
            q = R.gen(0)
            f = (q + a2 * q**2).O(3)

        if f.prec() < prec:
            R = f.parent()
            ext = [self.eigenvalue(n, names) for n in range(f.prec(), prec)]
            f = R(f.padded_list(f.prec()) + ext)
            self._q_expansion_dict[names] = f.add_bigoh(prec)
            return self._q_expansion_dict[names]
        else:
            return f.O(prec)

    def _q_expansion_basis_eigen(self, prec, names):
        r"""
        Return a basis of eigenforms corresponding to this space, which must be
        new, cuspidal and simple.

        EXAMPLES::

            sage: ModularSymbols(17, 2,sign=1).cuspidal_submodule()._q_expansion_basis_eigen(2, "a")
            [q + O(q^2)]
        """
        if self.is_simple():
            # should we perhaps check at this point if self is new?
            f = self.q_eigenform(prec, names)
            R = PowerSeriesRing(self.base_ring(), 'q')
            B = [R([f[i][j] for i in range(prec)], prec)
                 for j in range(self.rank())]
            return B
        else:
            raise NotImplementedError

    #########################################################################
    #
    #  Computation of a basis using linear functionals on the Hecke algebra.
    #
    #########################################################################

    def q_expansion_cuspforms(self, prec=None):
        r"""
        Return a function f(i,j) such that each value f(i,j) is the
        `q`-expansion, to the given precision, of an element of the
        corresponding space `S` of cusp forms.

        Together these functions span `S`. Here `i,j` are integers
        with `0\leq i,j < d`, where `d` is the dimension of ``self``.

        For a reduced echelon basis, use the function
        ``q_expansion_basis`` instead.

        More precisely, this function returns the `q`-expansions
        obtained by taking the `ij` entry of the matrices of the
        Hecke operators `T_n` acting on the subspace of the linear
        dual of modular symbols corresponding to ``self``.

        EXAMPLES::

            sage: S = ModularSymbols(11,2, sign=1).cuspidal_submodule()
            sage: f = S.q_expansion_cuspforms(8)
            sage: f(0,0)
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 + O(q^8)

        ::

            sage: S = ModularSymbols(37,2).cuspidal_submodule()
            sage: f = S.q_expansion_cuspforms(8)
            sage: f(0,0)
            q + q^3 - 2*q^4 - q^7 + O(q^8)
            sage: f(3,3)
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + 6*q^6 - q^7 + O(q^8)
            sage: f(1,2)
            q^2 + 2*q^3 - 2*q^4 + q^5 - 3*q^6 + O(q^8)

        ::

            sage: S = ModularSymbols(Gamma1(13),2,sign=-1).cuspidal_submodule()
            sage: f = S.q_expansion_cuspforms(8)
            sage: f(0,0)
            q - 2*q^2 + q^4 - q^5 + 2*q^6 + O(q^8)
            sage: f(0,1)
            -q^2 + 2*q^3 + q^4 - 2*q^5 - 2*q^6 + O(q^8)

        ::

            sage: S = ModularSymbols(1,12,sign=-1).cuspidal_submodule()
            sage: f = S.q_expansion_cuspforms(8)
            sage: f(0,0)
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 - 16744*q^7 + O(q^8)
        """
        if prec is None:
            prec = self.default_prec()
        if not self.is_cuspidal():
            raise ArithmeticError("self must be cuspidal")
        T = [self.dual_hecke_matrix(n) for n in range(1, prec)]
        R = PowerSeriesRing(self.base_ring(), 'q')
        return lambda i, j: R([0] + [t[i, j] for t in T], prec)

    def _q_expansion_basis_hecke_dual(self, prec):
        r"""
        Compute a basis of `q`-expansions for the associated space of cusp forms
        to the given precision, by using linear functionals on the Hecke
        algebra as described in William Stein's book (Algorithm 3.26, page 56)

        EXAMPLES::

            sage: ModularSymbols(37, 2).cuspidal_submodule()._q_expansion_basis_hecke_dual(12)
            [q + q^3 - 2*q^4 - q^7 - 2*q^9 + 3*q^11 + O(q^12),
            q^2 + 2*q^3 - 2*q^4 + q^5 - 3*q^6 - 4*q^9 - 2*q^10 + 4*q^11 + O(q^12)]
            sage: ModularSymbols(37, 2).cuspidal_submodule()._q_expansion_basis_hecke_dual(2)
            [q + O(q^2)]
        """
        from sage.misc.verbose import verbose
        d = self.dimension_of_associated_cuspform_space()
        prec = Integer(prec)
        if prec < 1:
            raise ValueError("prec (=%s) must be >= 1" % prec)
        d = min(prec - 1, d)
        K = self.base_ring()

        A = VectorSpace(K, prec - 1)
        M = MatrixSpace(K, prec - 1, self.dimension())

        V = A.zero_submodule()
        i = self.dimension() - 1
        j = 0

        t = verbose('computing basis to precision %s' % prec)
        while V.dimension() < d and i >= 0:
            v = [self.dual_hecke_matrix(n).column(i) for n in range(1, prec)]
            t = verbose('iteration: %s' % j, t)
            X = M(v).transpose()
            V += X.row_space()
            t = verbose('addition of row space: %s' % j, t)
            i -= 1
            j += 1

        R = PowerSeriesRing(K, 'q')
        B = V.basis()
        if len(B) < d:
            B += [V(0)] * (d - len(B))
        return [R([0] + b.list(), prec) for b in B]

    #########################################################################
    #
    #  Decomposition of spaces
    #
    ##########################################################################

#    def factorization(self):
#        """
#        Return a list of pairs `(S,e)` where `S` is simple
#        spaces of modular symbols and self is isomorphic to the direct sum
#        of the `S^e` as a module over the *anemic* Hecke algebra
#        adjoin the star involution.
#
#        ASSUMPTION: self is a module over the anemic Hecke algebra.
#        """
#        try:
#            return self._factorization
#        except AttributeError:
#            raise NotImplementedError

    def hecke_module_of_level(self, level):
        r"""
        Alias for ``self.modular_symbols_of_level(level)``.

        EXAMPLES::

            sage: ModularSymbols(11, 2).hecke_module_of_level(22)
            Modular Symbols space of dimension 7 for Gamma_0(22) of weight 2 with sign 0 over Rational Field
        """
        return self.modular_symbols_of_level(Integer(level))

    def sign(self):
        r"""
        Return the sign of ``self``.

        For efficiency reasons, it is often useful to compute in the
        (largest) quotient of modular symbols where the \* involution acts
        as +1, or where it acts as -1.

        INPUT:

        - ``ModularSymbols self`` -- arbitrary space of modular symbols

        OUTPUT:

        - ``-1`` -- if this is factor of quotient where \* acts as -1,

        - ``+1`` -- if this is factor of quotient where \* acts as +1,

        - ``0`` -- if this is full space of modular symbols (no quotient)

        EXAMPLES::

            sage: m = ModularSymbols(33)
            sage: m.rank()
            9
            sage: m.sign()
            0
            sage: m = ModularSymbols(33, sign=0)
            sage: m.sign()
            0
            sage: m.rank()
            9
            sage: m = ModularSymbols(33, sign=-1)
            sage: m.sign()
            -1
            sage: m.rank()
            3
        """
        return self.__sign

    def simple_factors(self):
        """
        Return a list modular symbols spaces `S` where `S`
        is simple spaces of modular symbols (for the anemic Hecke algebra)
        and ``self`` is isomorphic to the direct sum of the `S` with
        some multiplicities, as a module over the *anemic* Hecke algebra.

        For the multiplicities use factorization() instead.

        ASSUMPTION: self is a module over the anemic Hecke algebra.

        EXAMPLES::

            sage: ModularSymbols(1,100,sign=-1).simple_factors()
            [Modular Symbols subspace of dimension 8 of Modular Symbols space of dimension 8 for Gamma_0(1) of weight 100 with sign -1 over Rational Field]
            sage: ModularSymbols(1,16,0,GF(5)).simple_factors()
            [Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(1) of weight 16 with sign 0 over Finite Field of size 5,
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(1) of weight 16 with sign 0 over Finite Field of size 5,
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(1) of weight 16 with sign 0 over Finite Field of size 5]
        """
        return [S for S, _ in self.factorization()]

    def star_eigenvalues(self):
        """
        Return the eigenvalues of the star involution acting on ``self``.

        EXAMPLES::

            sage: M = ModularSymbols(11)
            sage: D = M.decomposition()
            sage: M.star_eigenvalues()
            [1, -1]
            sage: D[0].star_eigenvalues()
            [1]
            sage: D[1].star_eigenvalues()
            [1, -1]
            sage: D[1].plus_submodule().star_eigenvalues()
            [1]
            sage: D[1].minus_submodule().star_eigenvalues()
            [-1]
        """
        try:
            return self.__star_eigenvalues
        except AttributeError:
            pass
        if self.sign() != 0:
            return [self.sign()]
        M = self.star_involution().matrix()
        R = self.base_ring()
        if M == R(1):
            self.__star_eigenvalues = [R(1)]
        elif M == R(-1):
            self.__star_eigenvalues = [R(-1)]
        else:
            self.__star_eigenvalues = [R(1), R(-1)]
        return self.__star_eigenvalues

    def star_decomposition(self):
        r"""
        Decompose ``self`` into subspaces which are eigenspaces for the star
        involution.

        EXAMPLES::

            sage: ModularSymbols(Gamma1(19), 2).cuspidal_submodule().star_decomposition()
            [Modular Symbols subspace of dimension 7 of Modular Symbols space of dimension 31 for Gamma_1(19) of weight 2 with sign 0 over Rational Field,
             Modular Symbols subspace of dimension 7 of Modular Symbols space of dimension 31 for Gamma_1(19) of weight 2 with sign 0 over Rational Field]
        """
        S = self.star_involution()
        return S.decomposition()

    def integral_structure(self):
        r"""
        Return the `\ZZ`-structure of this modular symbols
        spaces generated by all integral modular symbols.

        EXAMPLES::

            sage: M = ModularSymbols(11,4)
            sage: M.integral_structure()
            Free module of degree 6 and rank 6 over Integer Ring
            Echelon basis matrix:
            [    1     0     0     0     0     0]
            [    0  1/14   1/7  5/14   1/2 13/14]
            [    0     0   1/2     0     0   1/2]
            [    0     0     0     1     0     0]
            [    0     0     0     0     1     0]
            [    0     0     0     0     0     1]
            sage: M.cuspidal_submodule().integral_structure()
            Free module of degree 6 and rank 4 over Integer Ring
            Echelon basis matrix:
            [     0   1/14    1/7   5/14    1/2 -15/14]
            [     0      0    1/2      0      0   -1/2]
            [     0      0      0      1      0     -1]
            [     0      0      0      0      1     -1]
        """
        try:
            return self.__integral_structure
        except AttributeError:
            pass
        A = self.ambient_hecke_module()
        I = A.integral_structure()
        if I.echelonized_basis_matrix().is_one() and (self.free_module().denominator() == 1):
            J = I.submodule(self.free_module().basis(), check=False, already_echelonized=True)
        else:
            J = self.free_module().intersection(I)
        self.__integral_structure = J
        return J

    def intersection_number(self, M):
        """
        Given modular symbols spaces ``self`` and ``M`` in some common ambient
        space, returns the intersection number of these two spaces.

        This is the index in their saturation of the sum of their
        underlying integral structures.

        If ``self`` and ``M`` are of weight two and defined over QQ,
        and correspond to newforms f and g, then this number equals
        the order of the intersection of the modular abelian varieties
        attached to f and g.

        EXAMPLES::

            sage: m = ModularSymbols(389,2)
            sage: d = m.decomposition(2)
            sage: eis = d[0]
            sage: ell = d[1]
            sage: af = d[-1]
            sage: af.intersection_number(eis)
            97
            sage: af.intersection_number(ell)
            400
        """
        if not isinstance(M, ModularSymbolsSpace):
            raise TypeError("M must be a modular symbols space")
        if M.ambient() != self.ambient():
            raise ValueError("self and M must be in the same ambient space.")
        A = self.integral_structure()
        B = M.integral_structure()
        return (A + B).index_in_saturation()

    def integral_basis(self):
        r"""
        Return a basis for the `\ZZ`-submodule of this
        modular symbols space spanned by the generators.

        Modular symbols spaces for congruence subgroups have a
        `\ZZ`-structure. Computing this
        `\ZZ`-structure is expensive, so by default modular
        symbols spaces for congruence subgroups in Sage are defined over
        `\QQ`. This function returns a tuple of independent
        elements in this modular symbols space whose
        `\ZZ`-span is the corresponding space of modular
        symbols over `\ZZ`.

        EXAMPLES::

            sage: M = ModularSymbols(11)
            sage: M.basis()
            ((1,0), (1,8), (1,9))
            sage: M.integral_basis()
            ((1,0), (1,8), (1,9))
            sage: S = M.cuspidal_submodule()
            sage: S.basis()
            ((1,8), (1,9))
            sage: S.integral_basis()
            ((1,8), (1,9))

        ::

            sage: M = ModularSymbols(13,4)
            sage: M.basis()
            ([X^2,(0,1)], [X^2,(1,4)], [X^2,(1,5)], [X^2,(1,7)], [X^2,(1,9)], [X^2,(1,10)], [X^2,(1,11)], [X^2,(1,12)])
            sage: M.integral_basis()
            ([X^2,(0,1)], 1/28*[X^2,(1,4)] + 2/7*[X^2,(1,5)] + 3/28*[X^2,(1,7)] + 11/14*[X^2,(1,9)] + 2/7*[X^2,(1,10)] + 11/28*[X^2,(1,11)] + 3/28*[X^2,(1,12)], [X^2,(1,5)], 1/2*[X^2,(1,7)] + 1/2*[X^2,(1,9)], [X^2,(1,9)], [X^2,(1,10)], [X^2,(1,11)], [X^2,(1,12)])
            sage: S = M.cuspidal_submodule()
            sage: S.basis()
            ([X^2,(1,4)] - [X^2,(1,12)], [X^2,(1,5)] - [X^2,(1,12)], [X^2,(1,7)] - [X^2,(1,12)], [X^2,(1,9)] - [X^2,(1,12)], [X^2,(1,10)] - [X^2,(1,12)], [X^2,(1,11)] - [X^2,(1,12)])
            sage: S.integral_basis()
            (1/28*[X^2,(1,4)] + 2/7*[X^2,(1,5)] + 3/28*[X^2,(1,7)] + 11/14*[X^2,(1,9)] + 2/7*[X^2,(1,10)] + 11/28*[X^2,(1,11)] - 53/28*[X^2,(1,12)], [X^2,(1,5)] - [X^2,(1,12)], 1/2*[X^2,(1,7)] + 1/2*[X^2,(1,9)] - [X^2,(1,12)], [X^2,(1,9)] - [X^2,(1,12)], [X^2,(1,10)] - [X^2,(1,12)], [X^2,(1,11)] - [X^2,(1,12)])

        This function currently raises a NotImplementedError on modular
        symbols spaces with character of order bigger than `2`:

        EXAMPLES::

            sage: M = ModularSymbols(DirichletGroup(13).0^2, 2); M
            Modular Symbols space of dimension 4 and level 13, weight 2, character [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2
            sage: M.basis()
            ((1,0), (1,5), (1,10), (1,11))
            sage: M.integral_basis()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        try:
            return self.__integral_basis
        except AttributeError:
            pass
        B = self.integral_structure().basis()
        self.__integral_basis = tuple([self(b) for b in B])
        return self.__integral_basis

    def integral_hecke_matrix(self, n):
        r"""
        Return the matrix of the `n`-th Hecke operator acting on the integral
        structure on ``self`` (as returned by ``self.integral_structure()``).

        This is often (but not always) different from the matrix
        returned by ``self.hecke_matrix``, even if the latter has
        integral entries.

        EXAMPLES::

            sage: M = ModularSymbols(6,4)
            sage: M.hecke_matrix(3)
            [27  0  0  0  6 -6]
            [ 0  1 -4  4  8 10]
            [18  0  1  0  6 -6]
            [18  0  4 -3  6 -6]
            [ 0  0  0  0  9 18]
            [ 0  0  0  0 12 15]
            sage: M.integral_hecke_matrix(3)
            [ 27   0   0   0   6  -6]
            [  0   1  -8   8  12  14]
            [ 18   0   5  -4  14   8]
            [ 18   0   8  -7   2 -10]
            [  0   0   0   0   9  18]
            [  0   0   0   0  12  15]
        """
        n = int(n)
        try:
            return self.__integral_hecke_matrix[n]
        except AttributeError:
            self.__integral_hecke_matrix = {}
        except KeyError:
            pass
        A = self.ambient_hecke_module()
        T = A.hecke_matrix(n)
        S = T.restrict(self.integral_structure()).change_ring(ZZ)
        self.__integral_hecke_matrix[n] = S
        return S

    def sturm_bound(self):
        r"""
        Return the Sturm bound for this space of modular symbols.

        Type ``sturm_bound?`` for more details.

        EXAMPLES::

            sage: ModularSymbols(11,2).sturm_bound()
            2
            sage: ModularSymbols(389,2).sturm_bound()
            65
            sage: ModularSymbols(1,12).sturm_bound()
            1
            sage: ModularSymbols(1,36).sturm_bound()
            3
            sage: ModularSymbols(DirichletGroup(31).0^2).sturm_bound()
            6
            sage: ModularSymbols(Gamma1(31)).sturm_bound()
            160
        """
        # For Gamma_0(N), n = \frac{k}{12}[\SL_2(\Z):\Gamma_0(N)]
        try:
            return self.__sturm_bound
        except AttributeError:
            if self.character() is not None:
                self.__sturm_bound = Gamma0(self.level()).sturm_bound(self.weight())
            else:
                self.__sturm_bound = self.group().sturm_bound(self.weight())
        return self.__sturm_bound

    def plus_submodule(self, compute_dual=True):
        """
        Return the subspace of ``self`` on which the star involution acts as +1.

        INPUT:

        - ``compute_dual`` -- boolean (default: ``True``); also
          compute dual subspace. This is useful for many algorithms.

        OUTPUT: subspace of modular symbols

        EXAMPLES::

            sage: ModularSymbols(17,2)
            Modular Symbols space of dimension 3 for Gamma_0(17) of weight 2 with sign 0 over Rational Field
            sage: ModularSymbols(17,2).plus_submodule()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(17) of weight 2 with sign 0 over Rational Field
        """
        return self.sign_submodule(+1, compute_dual)

    def minus_submodule(self, compute_dual=True):
        """
        Return the subspace of ``self`` on which the star involution acts as -1.

        INPUT:

        - ``compute_dual`` -- boolean (default: ``True``); also
          compute dual subspace. This is useful for many algorithms.

        OUTPUT: subspace of modular symbols

        EXAMPLES::

            sage: ModularSymbols(14,4)
            Modular Symbols space of dimension 12 for Gamma_0(14) of weight 4 with sign 0 over Rational Field
            sage: ModularSymbols(14,4).minus_submodule()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 12 for Gamma_0(14) of weight 4 with sign 0 over Rational Field
        """
        return self.sign_submodule(-1, compute_dual)

    def _compute_sign_submodule(self, sign, compute_dual=True):
        r"""
        Compute the submodule of ``self`` which is an eigenspace for the star involution with the given sign.

        INPUT:

        - ``sign`` -- integer; 1 or -1

        - ``compute_dual`` -- boolean (default: ``True``); also compute the
          dual submodule (useful for some algorithms)

        OUTPUT: a submodule of ``self``

        EXAMPLES::

            sage: ModularSymbols(Gamma1(11), 3)._compute_sign_submodule(-1)
            Modular Symbols subspace of dimension 10 of Modular Symbols space of dimension 20 for Gamma_1(11) of weight 3 with sign 0 over Rational Field
        """
        A = self.ambient()
        S = A.sign_submodule(sign, compute_dual=compute_dual)
        V = S.free_module().intersection(self.free_module())
        if compute_dual:
            W = S.dual_free_module()
            Y = self.dual_free_module()
            D = W.intersection(Y)
            M = A.submodule(V, D, check=False)
        else:
            M = A.submodule(V, check=False)
        M._set_sign(sign)
        return M

    def _set_sign(self, sign):
        r"""
        Store the sign of this module (used by various initialisation
        routines).

        INPUT:

        - (integer) sign (must be -1, 0 or 1)

        OUTPUT: none

        EXAMPLES::

            sage: ModularSymbols(11, 2)._set_sign(123)
            Traceback (most recent call last):
            ...
            ValueError: sign (=123) must be -1, 0, or 1
            sage: M = ModularSymbols(11, 2); M
            Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: M._set_sign(-1); M
            Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign -1 over Rational Field
            sage: M._set_sign(0)
        """
        sign = int(sign)
        if sign not in [-1, 0, 1]:
            raise ValueError("sign (=%s) must be -1, 0, or 1" % sign)
        self.__sign = sign

    def sign_submodule(self, sign, compute_dual=True):
        """
        Return the subspace of ``self`` that is fixed under the star involution.

        INPUT:

        - ``sign`` -- integer (either -1, 0 or +1)

        - ``compute_dual`` -- boolean (default: ``True``); also
          compute dual subspace. This is useful for many algorithms.

        OUTPUT: subspace of modular symbols

        EXAMPLES::

            sage: M = ModularSymbols(29,2)
            sage: M.sign_submodule(1)
            Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 5 for Gamma_0(29) of weight 2 with sign 0 over Rational Field
            sage: M.sign_submodule(-1)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 5 for Gamma_0(29) of weight 2 with sign 0 over Rational Field
            sage: M.sign_submodule(-1).sign()
            -1
        """
        sign = int(sign)
        if sign not in [-1, 0, 1]:
            raise ValueError("sign must be -1, 0 or 1")
        if self.sign() == sign:  # an easy case
            return self
        if self.sign() == -sign:  # another easy case
            return self.zero_submodule()
        if sign == 0:
            # if sign is zero then self.sign() isn't 0 because
            # of the above checks.
            raise ArithmeticError("There is no sign 0 subspace of a space of modular symbols with nonzero sign.")
        try:
            return self.__plus_submodule[(sign, compute_dual)]
        except AttributeError:
            self.__plus_submodule = {}
        except KeyError:
            pass
        P = self._compute_sign_submodule(sign, compute_dual)
        P.__star_eigenvalue = sign
        self.__plus_submodule[(sign, compute_dual)] = P
        return P

    def star_involution(self):
        """
        Return the star involution on ``self``, which is induced by complex
        conjugation on modular symbols.

        Not implemented in this abstract base class.

        EXAMPLES::

            sage: M = ModularSymbols(11, 2); sage.modular.modsym.space.ModularSymbolsSpace.star_involution(M)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def abelian_variety(self):
        """
        Return the corresponding abelian variety.

        INPUT:

        - ``self`` -- modular symbols space of weight 2 for a
          congruence subgroup such as Gamma0, Gamma1 or GammaH

        EXAMPLES::

            sage: ModularSymbols(Gamma0(11)).cuspidal_submodule().abelian_variety()
            Abelian variety J0(11) of dimension 1
            sage: ModularSymbols(Gamma1(11)).cuspidal_submodule().abelian_variety()
            Abelian variety J1(11) of dimension 1
            sage: ModularSymbols(GammaH(11,[3])).cuspidal_submodule().abelian_variety()
            Abelian variety JH(11,[3]) of dimension 1

        The abelian variety command only works on cuspidal modular symbols
        spaces::

            sage: M = ModularSymbols(37)
            sage: M[0].abelian_variety()
            Traceback (most recent call last):
            ...
            ValueError: self must be cuspidal
            sage: M[1].abelian_variety()
            Abelian subvariety of dimension 1 of J0(37)
            sage: M[2].abelian_variety()
            Abelian subvariety of dimension 1 of J0(37)
        """
        try:
            return self.__modular_abelian_variety
        except AttributeError:
            if not self.is_cuspidal():
                raise ValueError("self must be cuspidal")
            from sage.modular.abvar.abvar import ModularAbelianVariety_modsym
            A = ModularAbelianVariety_modsym(self, check=False)
            self.__modular_abelian_variety = A
            return A

    def rational_period_mapping(self):
        r"""
        Return the rational period mapping associated to ``self``.

        This is a homomorphism to a vector space whose kernel is the
        same as the kernel of the period mapping associated to
        ``self``. For this to exist, ``self`` must be Hecke equivariant.

        Use :meth:`integral_period_mapping` to obtain a homomorphism to a
        `\ZZ`-module, normalized so the image of integral modular symbols is
        exactly `\ZZ^n`.

        EXAMPLES::

            sage: M = ModularSymbols(37)
            sage: A = M[1]; A
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field
            sage: r = A.rational_period_mapping(); r
            Rational period mapping associated to Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field
            sage: r(M.0)
            (0, 0)
            sage: r(M.1)
            (1, 0)
            sage: r.matrix()
            [ 0  0]
            [ 1  0]
            [ 0  1]
            [-1 -1]
            [ 0  0]
            sage: r.domain()
            Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field
            sage: r.codomain()
            Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
        """
        try:
            return self.__rational_period_mapping
        except AttributeError:
            pass
        V = self.dual_free_module()
        # the rational period mapping is just dotting with
        # each element of V.
        A = V.basis_matrix().transpose()
        A.set_immutable()
        R = RationalPeriodMapping(self, A)
        self.__rational_period_mapping = R
        return R

    def integral_period_mapping(self):
        r"""
        Return the integral period mapping associated to ``self``.

        This is a homomorphism to a vector space whose kernel is the same
        as the kernel of the period mapping associated to ``self``, normalized
        so the image of integral modular symbols is exactly `\ZZ^n`.

        EXAMPLES::

            sage: m = ModularSymbols(23).cuspidal_submodule()
            sage: i = m.integral_period_mapping()
            sage: i
            Integral period mapping associated to Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
            sage: i.matrix()
            [-1/11  1/11     0  3/11]
            [    1     0     0     0]
            [    0     1     0     0]
            [    0     0     1     0]
            [    0     0     0     1]
            sage: [i(b) for b in m.integral_structure().basis()]
            [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
            sage: [i(b) for b in m.ambient_module().basis()]
            [(-1/11, 1/11, 0, 3/11),
             (1, 0, 0, 0),
             (0, 1, 0, 0),
             (0, 0, 1, 0),
             (0, 0, 0, 1)]

        We compute the image of the winding element::

            sage: m = ModularSymbols(37,sign=1)
            sage: a = m[1]
            sage: f = a.integral_period_mapping()
            sage: e = m([0,oo])
            sage: f(e)
            (-2/3)

        The input space must be cuspidal::

            sage: m = ModularSymbols(37,2,sign=1)
            sage: m.integral_period_mapping()
            Traceback (most recent call last):
            ...
            ValueError: integral mapping only defined for cuspidal spaces
        """
        try:
            return self.__integral_period_mapping
        except AttributeError:
            pass
        if self.base_ring() != QQ:
            raise ValueError("integral mapping only defined for spaces over QQ")
        if not self.is_cuspidal():
            raise ValueError("integral mapping only defined for cuspidal spaces")
        D = self.dual_free_module().basis_matrix().transpose()
        I = self.ambient_module().cuspidal_submodule().integral_structure().basis_matrix()
        # image of cuspidal integral submodule
        C = I * D
        if not C.is_one():
            if not C.is_square():
                C = (ZZ**C.ncols()).span(C.rows()).basis_matrix()
            D = D * C**(-1)
        D.set_immutable()
        R = IntegralPeriodMapping(self, D)
        self.__integral_period_mapping = R
        return R

    @cached_method
    def modular_symbols_of_sign(self, sign, bound=None):
        """
        Return a space of modular symbols with the same defining
        properties (weight, level, etc.) and Hecke eigenvalues as this
        space except with given sign.

        INPUT:

        - ``self`` -- a cuspidal space of modular symbols

        - ``sign`` -- integer, one of -1, 0, or 1

        - ``bound`` -- integer (default: ``None``); if specified
          only use Hecke operators up to the given bound

        EXAMPLES::

            sage: S = ModularSymbols(Gamma0(11),2,sign=0).cuspidal_subspace()
            sage: S
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: S.modular_symbols_of_sign(-1)
            Modular Symbols space of dimension 1 for Gamma_0(11) of weight 2 with sign -1 over Rational Field

        ::

            sage: S = ModularSymbols(43,2,sign=1)[2]; S
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1 over Rational Field
            sage: S.modular_symbols_of_sign(-1)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(43) of weight 2 with sign -1 over Rational Field

        ::

            sage: S.modular_symbols_of_sign(0)
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 7 for Gamma_0(43) of weight 2 with sign 0 over Rational Field

        ::

            sage: S = ModularSymbols(389,sign=1)[3]; S
            Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 33 for Gamma_0(389) of weight 2 with sign 1 over Rational Field
            sage: S.modular_symbols_of_sign(-1)
            Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 32 for Gamma_0(389) of weight 2 with sign -1 over Rational Field
            sage: S.modular_symbols_of_sign(0)
            Modular Symbols subspace of dimension 6 of Modular Symbols space of dimension 65 for Gamma_0(389) of weight 2 with sign 0 over Rational Field

        ::

            sage: S = ModularSymbols(23,sign=1,weight=4)[2]; S
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 7 for Gamma_0(23) of weight 4 with sign 1 over Rational Field
            sage: S.modular_symbols_of_sign(1) is S
            True
            sage: S.modular_symbols_of_sign(0)
            Modular Symbols subspace of dimension 8 of Modular Symbols space of dimension 12 for Gamma_0(23) of weight 4 with sign 0 over Rational Field
            sage: S.modular_symbols_of_sign(-1)
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 4 with sign -1 over Rational Field
        """
        if sign == self.sign():
            return self
        if not self.is_cuspidal():
            raise ValueError("self must be cuspidal for modular symbols space with given sign to be defined.")
        d = self.dimension()
        if d == 0:
            return self
        if sign != 0:
            if self.sign() == 0:
                d = d // 2
        elif sign == 0:      # self has nonzero sign
            d = 2 * d
        B = self.ambient_module().modular_symbols_of_sign(sign)
        p = 2
        if bound is None:
            bound = self.hecke_bound()
        while B.dimension() > d and p <= bound:
            while self.level() % p == 0:
                p = next_prime(p)
            f = self.hecke_polynomial(p)
            g = prod(g for g, _ in f.factor())   # square free part
            t = B.hecke_operator(p)
            s = g(t)
            B = s.kernel()
            p = next_prime(p)
        return B

    #########################################################
    # Cuspidal torsion groups
    #########################################################

    def abvarquo_cuspidal_subgroup(self):
        """
        Compute the rational subgroup of the cuspidal subgroup (as an
        abstract abelian group) of the abelian variety quotient A of
        the relevant modular Jacobian attached to this modular symbols
        space.

        We assume that ``self`` is defined over QQ and has weight 2.  If
        the sign of ``self`` is not 0, then the power of 2 may be wrong.

        EXAMPLES::

            sage: D = ModularSymbols(66,2,sign=0).cuspidal_subspace().new_subspace().decomposition()
            sage: D[0].abvarquo_cuspidal_subgroup()
            Finitely generated module V/W over Integer Ring with invariants (3)
            sage: [A.abvarquo_cuspidal_subgroup().invariants() for A in D]
            [(3,), (2,), ()]
            sage: D = ModularSymbols(66,2,sign=1).cuspidal_subspace().new_subspace().decomposition()
            sage: [A.abvarquo_cuspidal_subgroup().invariants() for A in D]
            [(3,), (2,), ()]
            sage: D = ModularSymbols(66,2,sign=-1).cuspidal_subspace().new_subspace().decomposition()
            sage: [A.abvarquo_cuspidal_subgroup().invariants() for A in D]
            [(), (), ()]
        """
        try:
            return self.__abvarquo_cuspidal_subgroup
        except AttributeError:
            pass
        if self.base_ring() != QQ:
            raise ValueError("base ring must be QQ")
        if self.weight() != 2:
            raise NotImplementedError("only implemented when weight is 2")
        M = self.ambient_module()
        phi = self.integral_period_mapping()

        # Make a list of all the finite cusps.
        P = [c for c in M.cusps() if not c.is_infinity()]

        # Compute the images of the cusp classes (c)-(oo) in the
        # rational homology of the quotient modular abelian variety.
        ims = [phi(M([c, infinity])) for c in P]

        # Take the span of the ims over ZZ
        A = phi.codomain().span(ims, ZZ)

        # The cuspidal subgroup is then the quotient of that module +
        # H_1(A) by H_1(A)
        C = (A.ambient_module() + A) / A.ambient_module()

        self.__abvarquo_cuspidal_subgroup = C
        return C

    def abvarquo_rational_cuspidal_subgroup(self):
        r"""
        Compute the rational subgroup of the cuspidal subgroup (as an
        abstract abelian group) of the abelian variety quotient A of
        the relevant modular Jacobian attached to this modular symbols
        space.

        If C is the subgroup of A generated by differences of
        cusps, then C is equipped with an action of Gal(Qbar/Q), and
        this function computes the fixed subgroup, i.e., C(Q).

        We assume that ``self`` is defined over QQ and has weight 2.  If
        the sign of ``self`` is not 0, then the power of 2 may be wrong.

        EXAMPLES:

        First we consider the fairly straightforward level 37 case,
        where the torsion subgroup of the optimal quotients (which are
        all elliptic curves) are all cuspidal::

            sage: M = ModularSymbols(37).cuspidal_subspace().new_subspace()
            sage: D = M.decomposition()
            sage: [(A.abvarquo_rational_cuspidal_subgroup().invariants(), A.T(19)[0,0]) for A in D]
            [((), 0), ((3,), 2)]
            sage: [(E.torsion_subgroup().invariants(),E.ap(19)) for E in cremona_optimal_curves([37])]
            [((), 0), ((3,), 2)]

        Next we consider level 54, where the rational cuspidal
        subgroups of the quotients are also cuspidal::

            sage: M = ModularSymbols(54).cuspidal_subspace().new_subspace()
            sage: D = M.decomposition()
            sage: [A.abvarquo_rational_cuspidal_subgroup().invariants() for A in D]
            [(3,), (3,)]
            sage: [E.torsion_subgroup().invariants() for E in cremona_optimal_curves([54])]
            [(3,), (3,)]

        Level 66 is interesting, since not all torsion of the quotient
        is rational. In fact, for each elliptic curve quotient, the
        `\QQ`-rational subgroup of the image of the cuspidal subgroup
        in the quotient is a nontrivial subgroup of `E(\QQ)_{tor}`.
        Thus not all torsion in the quotient is cuspidal!::

            sage: M = ModularSymbols(66).cuspidal_subspace().new_subspace()
            sage: D = M.decomposition()
            sage: [(A.abvarquo_rational_cuspidal_subgroup().invariants(), A.T(19)[0,0]) for A in D]
            [((3,), -4), ((2,), 4), ((), 0)]
            sage: [(E.torsion_subgroup().invariants(),E.ap(19)) for E in cremona_optimal_curves([66])]
            [((6,), -4), ((4,), 4), ((10,), 0)]
            sage: [A.abelian_variety().rational_cuspidal_subgroup().invariants() for A in D]
            [[6], [4], [10]]

        In this example, the abelian varieties involved all having
        dimension bigger than 1 (unlike above).  We find that all torsion
        in the quotient in each of these cases is cuspidal::

            sage: M = ModularSymbols(125).cuspidal_subspace().new_subspace()
            sage: D = M.decomposition()
            sage: [A.abvarquo_rational_cuspidal_subgroup().invariants() for A in D]
            [(), (5,), (5,)]
            sage: [A.abelian_variety().rational_torsion_subgroup().multiple_of_order() for A in D]
            [1, 5, 5]
        """
        try:
            return self.__abvarquo_rational_cuspidal_subgroup
        except AttributeError:
            pass
        if self.base_ring() != QQ:
            raise ValueError("base ring must be QQ")
        if self.weight() != 2:
            raise NotImplementedError("only implemented when weight is 2")
        if not isinstance(self.group(), Gamma0_class):
            # todo -- do Gamma1 and GammaH, which are easy
            raise NotImplementedError("only implemented when group is Gamma0")
        N = self.level()
        if N.is_squarefree():
            return self.abvarquo_cuspidal_subgroup()

        M = self.ambient_module()
        phi = self.integral_period_mapping()

        # Make a list of all the finite cusps.
        P = [c for c in M.cusps() if not c.is_infinity()]

        # Define the vector space V, which we think of as
        # the vector space with basis (c)-(oo), where c runs
        # through the finite cusp *classes*.
        V = ZZ**len(P)   # vector space on (c)-(oo)

        # Compute the images of the cusp classes (c)-(oo) in the
        # rational homology of the quotient modular abelian variety.
        ims = [phi(M([c, infinity])) for c in P]

        # Take the span of the ims over ZZ
        A = phi.codomain().span(ims, ZZ)

        # The cuspidal subgroup is then the quotient of that module +
        # H_1(A) by H_1(A)
        C = (A.ambient_module() + A) / A.ambient_module()

        # Make fgp module version of V.
        D = V / V.zero_submodule()
        psi = D.hom([C(x) for x in ims])

        # The rational cuspidal subgroup is got by intersecting kernels
        # of tau - 1, for all automorphisms tau.
        G = Zmod(N).unit_gens()
        CQ = C
        for t in G:
            T = self._matrix_of_galois_action(t, P) - 1
            if not T:
                continue
            im_gens = [psi(psi.lift(g).lift() * T) for g in CQ.gens()]
            h = CQ.hom(im_gens)
            CQ = h.kernel()
            if CQ.cardinality() == 1:
                break  # done -- no point in wasting more time shrinking CQ

        self.__abvarquo_rational_cuspidal_subgroup = CQ
        return CQ

    def _matrix_of_galois_action(self, t, P):
        """
        Compute the matrix of the action of the element of the
        cyclotomic Galois group defined by t on the set of cusps in P
        (which is the set of finite cusps).

        This function is used internally by the (rational) cuspidal
        subgroup and quotient functions.

        INPUT:

        - ``t`` -- integer

        - ``P`` -- list of cusps

        EXAMPLES:

        We compute the matrix of the element of the Galois group
        associated to 5 and 31 for level 32.  In the first case the
        Galois action is trivial, and in the second it is
        nontrivial. ::

            sage: M = ModularSymbols(32)
            sage: P = [c for c in Gamma0(32).cusps() if not c.is_infinity()]
            sage: M._matrix_of_galois_action(5, P)
            [1 0 0 0 0 0 0]
            [0 1 0 0 0 0 0]
            [0 0 1 0 0 0 0]
            [0 0 0 1 0 0 0]
            [0 0 0 0 1 0 0]
            [0 0 0 0 0 1 0]
            [0 0 0 0 0 0 1]
            sage: z = M._matrix_of_galois_action(31, P); z
            [1 0 0 0 0 0 0]
            [0 1 0 0 0 0 0]
            [0 0 0 0 1 0 0]
            [0 0 0 0 0 0 1]
            [0 0 1 0 0 0 0]
            [0 0 0 0 0 1 0]
            [0 0 0 1 0 0 0]
            sage: z.charpoly().factor()
            (x + 1)^2 * (x - 1)^5
        """
        N = self.level()
        from sage.matrix.constructor import matrix
        A = matrix(ZZ, len(P))
        for i, c in enumerate(P):
            d = c.galois_action(t, N)
            for j, e in enumerate(P):
                if d.is_gamma0_equiv(e, N, False):
                    A[i, j] = 1
        A.set_immutable()
        return A


class PeriodMapping(SageObject):
    r"""
    Base class for representing a period mapping attached to a space of modular
    symbols.

    To be used via the derived classes :class:`RationalPeriodMapping` and
    :class:`IntegralPeriodMapping`.
    """
    def __init__(self, modsym, A):
        r"""
        Standard initialisation function.

        INPUT:

        - ``modsym`` -- a space of modular symbols

        - ``A`` -- matrix of the associated period map

        EXAMPLES::

            sage: ModularSymbols(2, 8).cuspidal_submodule().integral_period_mapping() # indirect doctest
            Integral period mapping associated to Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 for Gamma_0(2) of weight 8 with sign 0 over Rational Field
        """
        self.__modsym = modsym
        self.__domain = modsym.ambient_module()
        self.__A = A
        A.set_immutable()

    def modular_symbols_space(self):
        r"""
        Return the space of modular symbols to which this period mapping
        corresponds.

        EXAMPLES::

            sage: ModularSymbols(17, 2).rational_period_mapping().modular_symbols_space()
            Modular Symbols space of dimension 3 for Gamma_0(17) of weight 2 with sign 0 over Rational Field
        """
        return self.__modsym

    def __call__(self, x):
        r"""
        Evaluate this mapping at an element of the domain.

        EXAMPLES::

            sage: M = ModularSymbols(17, 2).cuspidal_submodule().integral_period_mapping()
            sage: M(vector([1,0,2]))
            (0, 9/4)
        """
        if isinstance(x, FreeModuleElement):
            v = x
        else:
            v = self.__domain(x).element()
        return v * self.__A

    def matrix(self):
        r"""
        Return the matrix of this period mapping.

        EXAMPLES::

            sage: ModularSymbols(11, 2).cuspidal_submodule().integral_period_mapping().matrix()
            [  0 1/5]
            [  1   0]
            [  0   1]
        """
        return self.__A

    def domain(self):
        r"""
        Return the domain of this mapping (which is the ambient space of the
        corresponding modular symbols space).

        EXAMPLES::

            sage: ModularSymbols(17, 2).cuspidal_submodule().integral_period_mapping().domain()
            Modular Symbols space of dimension 3 for Gamma_0(17) of weight 2 with sign 0 over Rational Field
        """
        return self.__domain

    def codomain(self):
        r"""
        Return the codomain of this mapping.

        EXAMPLES:

        Note that this presently returns the wrong answer, as a consequence of
        various bugs in the free module routines::

            sage: ModularSymbols(11, 2).cuspidal_submodule().integral_period_mapping().codomain()
            Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
        """
        return self.__A.row_module()


class RationalPeriodMapping(PeriodMapping):
    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: ModularSymbols(40,2).rational_period_mapping()._repr_()
            'Rational period mapping associated to Modular Symbols space of dimension 13 for Gamma_0(40) of weight 2 with sign 0 over Rational Field'
        """
        return "Rational period mapping associated to %s" % self.modular_symbols_space()


class IntegralPeriodMapping(PeriodMapping):
    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: ModularSymbols(40,2).cuspidal_submodule().integral_period_mapping()._repr_()
            'Integral period mapping associated to Modular Symbols subspace of dimension 6 of Modular Symbols space of dimension 13 for Gamma_0(40) of weight 2 with sign 0 over Rational Field'
        """
        return "Integral period mapping associated to %s" % self.modular_symbols_space()
