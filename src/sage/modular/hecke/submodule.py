# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.libs.flint sage.libs.pari
"""
Submodules of Hecke modules
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

import sage.arith.all as arith
from sage.misc.verbose import verbose
from sage.misc.cachefunc import cached_method
from sage.modules.free_module import FreeModule_generic
from sage.structure.richcmp import richcmp_method, richcmp_not_equal

from . import module


def is_HeckeSubmodule(x):
    r"""
    Return ``True`` if x is of type HeckeSubmodule.

    EXAMPLES::

        sage: sage.modular.hecke.submodule.is_HeckeSubmodule(ModularForms(1, 12))
        doctest:warning...
        DeprecationWarning: the function is_HeckeSubmodule is deprecated;
        use 'isinstance(..., HeckeSubmodule)' instead
        See https://github.com/sagemath/sage/issues/37895 for details.
        False
        sage: sage.modular.hecke.submodule.is_HeckeSubmodule(CuspForms(1, 12))
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(37895, "the function is_HeckeSubmodule is deprecated; use 'isinstance(..., HeckeSubmodule)' instead")
    return isinstance(x, HeckeSubmodule)


@richcmp_method
class HeckeSubmodule(module.HeckeModule_free_module):
    """
    Submodule of a Hecke module.
    """
    def __init__(self, ambient, submodule, dual_free_module=None, check=True):
        r"""
        Initialise a submodule of an ambient Hecke module.

        INPUT:

        - ``ambient`` -- an ambient Hecke module

        - ``submodule`` -- a free module over the base ring which is a submodule
          of the free module attached to the ambient Hecke module. This should
          be invariant under all Hecke operators.

        - ``dual_free_module`` -- the submodule of the dual of the ambient
          module corresponding to this submodule (or ``None``)

        - ``check`` -- whether or not to explicitly check that the submodule is
          Hecke equivariant

        EXAMPLES::

            sage: CuspForms(1,60) # indirect doctest
            Cuspidal subspace of dimension 5 of Modular Forms space of dimension 6 for Modular Group SL(2,Z) of weight 60 over Rational Field

            sage: M = ModularForms(4,10)
            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.submodule(M.basis()[:3]).free_module())
            sage: S
            Rank 3 submodule of a Hecke module of level 4

            sage: S == loads(dumps(S))
            True
        """
        from . import ambient_module
        if not isinstance(ambient, ambient_module.AmbientHeckeModule):
            raise TypeError("ambient must be an ambient Hecke module")
        if not isinstance(submodule, FreeModule_generic):
            raise TypeError("submodule must be a free module")
        if not submodule.is_submodule(ambient.free_module()):
            raise ValueError("submodule must be a submodule of the ambient free module")

        if check:
            if not ambient._is_hecke_equivariant_free_module(submodule):
                raise ValueError("The submodule must be invariant under all Hecke operators.")

        self.__ambient = ambient
        self.__submodule = submodule
        module.HeckeModule_free_module.__init__(self, ambient.base_ring(),
                                                ambient.level(),
                                                ambient.weight())
        if dual_free_module is not None:
            if not isinstance(dual_free_module, FreeModule_generic):
                raise TypeError("dual_free_module must be a free module")
            if dual_free_module.rank() != submodule.rank():
                raise ArithmeticError("dual_free_module must have the same rank as submodule")
            self.dual_free_module.set_cache(dual_free_module)

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: M = ModularForms(4,10)
            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.submodule(M.basis()[:3]).free_module())
            sage: S._repr_()
            'Rank 3 submodule of a Hecke module of level 4'
        """
        return "Rank %s submodule of a Hecke module of level %s" % (
            self.rank(), self.level())

    def __add__(self, other):
        r"""
        Sum of ``self`` and ``other`` (as submodules of a common ambient
        module).

        EXAMPLES::

            sage: M = ModularForms(4,10)
            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.submodule(M.basis()[:3]).free_module())
            sage: E = sage.modular.hecke.submodule.HeckeSubmodule(M, M.submodule(M.basis()[3:]).free_module())
            sage: S + E # indirect doctest
            Modular Forms subspace of dimension 6 of Modular Forms space of dimension 6 for Congruence Subgroup Gamma0(4) of weight 10 over Rational Field
        """
        if not isinstance(other, module.HeckeModule_free_module):
            raise TypeError("other (=%s) must be a Hecke module." % other)
        if self.ambient() != other.ambient():
            raise ArithmeticError("sum only defined for submodules of a common ambient space")
        if other.is_ambient():
            return other
        # Neither is ambient
        M = self.free_module() + other.free_module()
        return self.ambient().submodule(M, check=False)

    def _element_constructor_(self, x, check=True):
        """
        Coerce x into the ambient module and checks that x is in this
        submodule.

        EXAMPLES::

            sage: M = ModularSymbols(37)
            sage: S = M.cuspidal_submodule()
            sage: M([0,oo])
            -(1,0)
            sage: S([0,oo])
            Traceback (most recent call last):
            ...
            TypeError: x does not coerce to an element of this Hecke module
            sage: S([-1/23,0])
            (1,23)
        """
        z = self.ambient_hecke_module()(x).element()
        if check and z not in self.__submodule:
            raise TypeError("x does not coerce to an element of this Hecke module")
        return self.element_class(self, z)

    def __richcmp__(self, other, op):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: M = ModularSymbols(12,6)
            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.cuspidal_submodule().free_module())
            sage: T = sage.modular.hecke.submodule.HeckeSubmodule(M, M.new_submodule().free_module())
            sage: S
            Rank 14 submodule of a Hecke module of level 12
            sage: T
            Rank 0 submodule of a Hecke module of level 12
            sage: S > T
            True
            sage: T < S
            True
            sage: S == S
            True
        """
        if not isinstance(other, module.HeckeModule_free_module):
            return NotImplemented
        lx = self.ambient()
        rx = other.ambient()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)
        return self.free_module()._echelon_matrix_richcmp(other.free_module(), op)

    ################################
    # Semi-Private functions
    ################################
    def _compute_dual_hecke_matrix(self, n):
        """
        Compute the matrix for the `n`-th Hecke operator acting on
        the dual of ``self``.

        EXAMPLES::

            sage: M = ModularForms(4,10)
            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.submodule(M.basis()[:3]).free_module())
            sage: S._compute_dual_hecke_matrix(3)
            [    0     0     1]
            [    0  -156     0]
            [35568     0    72]
            sage: CuspForms(4,10).dual_hecke_matrix(3)
            [    0     0     1]
            [    0  -156     0]
            [35568     0    72]
        """
        A = self.ambient_hecke_module().dual_hecke_matrix(n)
        check = arith.gcd(self.level(), n) != 1
        return A.restrict(self.dual_free_module(), check=check)

    def _compute_hecke_matrix(self, n):
        r"""
        Compute the matrix of the `n`-th Hecke operator acting on this space, by
        calling the corresponding function for the ambient space and
        restricting. If `n` is not coprime to the level, we check that the
        restriction is well-defined.

        EXAMPLES::

            sage: R.<q> = QQ[[]]
            sage: M = ModularForms(2, 12)
            sage: f = M(q^2 - 24*q^4 + O(q^6))
            sage: A = M.submodule(M.free_module().span([f.element()]),check=False)
            sage: sage.modular.hecke.submodule.HeckeSubmodule._compute_hecke_matrix(A, 3)
            [252]
            sage: sage.modular.hecke.submodule.HeckeSubmodule._compute_hecke_matrix(A, 4)
            Traceback (most recent call last):
            ...
            ArithmeticError: subspace is not invariant under matrix
        """
        A = self.ambient_hecke_module().hecke_matrix(n)
        check = arith.gcd(self.level(), n) != 1
        return A.restrict(self.free_module(), check=check)

    def _compute_diamond_matrix(self, d):
        r"""
        EXAMPLES::

            sage: f = ModularSymbols(Gamma1(13),2,sign=1).cuspidal_subspace().decomposition()[0]
            sage: a = f.diamond_bracket_operator(2).matrix() # indirect doctest
            sage: a.charpoly()
            x^2 - x + 1
            sage: a^12
            [1 0]
            [0 1]
        """
        return self.ambient_hecke_module().diamond_bracket_matrix(d).restrict(self.free_module(), check=False)

    def _compute_atkin_lehner_matrix(self, d):
        """
        Compute the Atkin-Lehner matrix corresponding to the
        divisor d of the level of ``self``.

        EXAMPLES::

            sage: M = ModularSymbols(4,10)
            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.cuspidal_submodule().free_module())
            sage: S
            Rank 6 submodule of a Hecke module of level 4
            sage: S._compute_atkin_lehner_matrix(1)
            [1 0 0 0 0 0]
            [0 1 0 0 0 0]
            [0 0 1 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 1 0]
            [0 0 0 0 0 1]
        """
        A = self.ambient_hecke_module()._compute_atkin_lehner_matrix(d)
        return A.restrict(self.free_module(), check=True)

    def _set_dual_free_module(self, V):
        """
        Set the dual free module of ``self`` to V. Here V must be a vector
        space of the same dimension as ``self``, embedded in a space of
        the same dimension as the ambient space of ``self``.

        EXAMPLES::

            sage: M = ModularSymbols(4,10)
            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.cuspidal_submodule().free_module())
            sage: S._set_dual_free_module(M.cuspidal_submodule().dual_free_module())
            sage: S._set_dual_free_module(S)
        """
        if V.degree() != self.ambient_hecke_module().rank():
            raise ArithmeticError("The degree of V must equal the rank of the ambient space.")
        if V.rank() != self.rank():
            raise ArithmeticError("The rank of V must equal the rank of self.")
        self.dual_free_module.set_cache(V)

    ################################
    # Public functions
    ################################

    def ambient_hecke_module(self):
        r"""
        Return the ambient Hecke module of which this is a submodule.

        EXAMPLES::

            sage: CuspForms(2, 12).ambient_hecke_module()
            Modular Forms space of dimension 4 for Congruence Subgroup Gamma0(2) of weight 12 over Rational Field
        """
        return self.__ambient

    def ambient(self):
        r"""
        Synonym for ambient_hecke_module.

        EXAMPLES::

            sage: CuspForms(2, 12).ambient()
            Modular Forms space of dimension 4 for Congruence Subgroup Gamma0(2) of weight 12 over Rational Field
        """
        return self.__ambient

    @cached_method
    def complement(self, bound=None):
        """
        Return the largest Hecke-stable complement of this space.

        EXAMPLES::

            sage: M = ModularSymbols(15, 6).cuspidal_subspace()
            sage: M.complement()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 20 for Gamma_0(15) of weight 6 with sign 0 over Rational Field
            sage: E = EllipticCurve("128a")
            sage: ME = E.modular_symbol_space()
            sage: ME.complement()
            Modular Symbols subspace of dimension 17 of Modular Symbols space of dimension 18 for Gamma_0(128) of weight 2 with sign 1 over Rational Field
        """

        if self.dual_free_module.is_in_cache():
            D = self.dual_free_module()
            V = D.basis_matrix().right_kernel()
            return self.submodule(V, check=False)

        if self.is_ambient():
            return self.ambient_hecke_module().zero_submodule()

        if self.is_zero():
            return self.ambient_hecke_module()

        if self.is_full_hecke_module():
            anemic = False
        else:
            anemic = True

        # TODO: optimize in some cases by computing image of
        # complementary factor instead of kernel...?
        verbose("computing")
        N = self.level()
        A = self.ambient_hecke_module()
        V = A.free_module()
        p = 2
        if bound is None:
            bound = A.hecke_bound()
        while True:
            if anemic:
                while N % p == 0:
                    p = arith.next_prime(p)
            verbose("using T_%s" % p)
            f = self.hecke_polynomial(p)
            T = A.hecke_matrix(p)
            g = T.charpoly('x')
            V = T.kernel_on(V, poly=g // f, check=False)
            if V.rank() + self.rank() <= A.rank():
                break
            p = arith.next_prime(p)
            if p > bound:  # to avoid computing Hecke bound unless necessary
                break

        if V.rank() + self.rank() == A.rank():
            C = A.submodule(V, check=False)
            return C

        # first attempt to compute the complement failed, we now try
        # the following naive approach: decompose the ambient space,
        # decompose self, and sum the pieces of ambient that are not
        # subspaces of self
        verbose("falling back on naive algorithm")
        D = A.decomposition()
        C = A.zero_submodule()
        for X in D:
            if self.intersection(X).dimension() == 0:
                C = C + X
        if C.rank() + self.rank() == A.rank():
            return C

        # failed miserably
        raise RuntimeError("Computation of complementary space failed (cut down to rank %s, but should have cut down to rank %s)." % (
            V.rank(), A.rank() - self.rank()))

    def degeneracy_map(self, level, t=1):
        """
        The `t`-th degeneracy map from ``self`` to the space of ambient modular
        symbols of the given level. The level of ``self`` must be a divisor or
        multiple of level, and `t` must be a divisor of the quotient.

        INPUT:

        - ``level`` -- positive integer; the level of the codomain of the
          map

        - ``t`` -- integer; the parameter of the degeneracy map,
          i.e., the map is related to `f(q)` - `f(q^t)`

        OUTPUT:

        A linear function from ``self`` to the space of modular symbols of given
        level with the same weight, character, sign, etc., as this space.

        EXAMPLES::

            sage: D = ModularSymbols(10,4).cuspidal_submodule().decomposition(); D
            [
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 10 for Gamma_0(10) of weight 4 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 10 for Gamma_0(10) of weight 4 with sign 0 over Rational Field
            ]
            sage: d = D[1].degeneracy_map(5); d
            Hecke module morphism defined by the matrix
            [   0    0   -1    1]
            [   0  1/2  3/2   -2]
            [   0   -1    1    0]
            [   0 -3/4 -1/4    1]
            Domain: Modular Symbols subspace of dimension 4 of Modular Symbols space ...
            Codomain: Modular Symbols space of dimension 4 for Gamma_0(5) of weight ...

        ::

            sage: d.rank()
            2
            sage: d.kernel()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 10 for Gamma_0(10) of weight 4 with sign 0 over Rational Field
            sage: d.image()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 for Gamma_0(5) of weight 4 with sign 0 over Rational Field
        """
        d = self.ambient_hecke_module().degeneracy_map(level, t)
        return d.restrict_domain(self)

    @cached_method
    def dual_free_module(self, bound=None, anemic=True, use_star=True):
        r"""
        Compute embedded dual free module if possible.

        In general this will not be possible, e.g., if this space is
        not Hecke equivariant, possibly if it is not cuspidal, or if
        the characteristic is not 0. In all these cases we raise a
        :exc:`RuntimeError` exception.

        If ``use_star`` is ``True`` (which is the default), we also use the +/-
        eigenspaces for the star operator to find the dual free module of ``self``.
        If ``self`` does not have a star involution, ``use_star`` will automatically be
        set to ``False``.

        EXAMPLES::

            sage: M = ModularSymbols(11, 2)
            sage: M.dual_free_module()
            Vector space of dimension 3 over Rational Field
            sage: Mpc = M.plus_submodule().cuspidal_submodule()
            sage: Mcp = M.cuspidal_submodule().plus_submodule()
            sage: Mcp.dual_free_module() == Mpc.dual_free_module()
            True
            sage: Mpc.dual_free_module()
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [  1 5/2   5]

            sage: M = ModularSymbols(35,2).cuspidal_submodule()
            sage: M.dual_free_module(use_star=False)
            Vector space of degree 9 and dimension 6 over Rational Field
            Basis matrix:
            [   1    0    0    0   -1    0    0    4   -2]
            [   0    1    0    0    0    0    0 -1/2  1/2]
            [   0    0    1    0    0    0    0 -1/2  1/2]
            [   0    0    0    1   -1    0    0    1    0]
            [   0    0    0    0    0    1    0   -2    1]
            [   0    0    0    0    0    0    1   -2    1]

            sage: M = ModularSymbols(40,2)
            sage: Mmc = M.minus_submodule().cuspidal_submodule()
            sage: Mcm = M.cuspidal_submodule().minus_submodule()
            sage: Mcm.dual_free_module() == Mmc.dual_free_module()
            True
            sage: Mcm.dual_free_module()
            Vector space of degree 13 and dimension 3 over Rational Field
            Basis matrix:
            [ 0  1  0  0  0  0  1  0 -1 -1  1 -1  0]
            [ 0  0  1  0 -1  0 -1  0  1  0  0  0  0]
            [ 0  0  0  0  0  1  1  0 -1  0  0  0  0]

            sage: M = ModularSymbols(43).cuspidal_submodule()
            sage: S = M[0].plus_submodule() + M[1].minus_submodule()
            sage: S.dual_free_module(use_star=False)
            Traceback (most recent call last):
            ...
            RuntimeError: Computation of complementary space failed (cut down to rank 7, but should have cut down to rank 4).
            sage: S.dual_free_module().dimension() == S.dimension()
            True

        We test that :issue:`5080` is fixed::

            sage: EllipticCurve('128a').congruence_number()
            32
        """

        # if we know the complement we can read off the dual module
        if self.complement.is_in_cache():
            verbose('This module knows its complement already -- cheating in dual_free_module')
            C = self.complement()
            V = C.basis_matrix().right_kernel()
            return V

        verbose("computing dual")

        A = self.ambient_hecke_module()

        if self.dimension() == 0:
            return A.zero_submodule()

        if A.dimension() == self.dimension():
            return A.free_module()

        # ALGORITHM: Compute the char poly of each Hecke operator on
        # the submodule, then use it to cut out a submodule of the
        # dual.  If the dimension cuts down to the dimension of self
        # terminate with success.  If it stays larger beyond the Sturm
        # bound, raise a RuntimeError exception.

        # In the case that the sign of self is not 1, we need to use
        # the star involution as well as the Hecke operators in order
        # to find the dual of self.
        #
        # Note that one needs to comment out the line caching the
        # result of this computation below in order to get meaningful
        # timings.

        # If the star involution doesn't make sense for self, then we
        # can't use it.
        if not hasattr(self, 'star_eigenvalues'):
            use_star = False

        if use_star:
            # If the star involution has both + and - eigenspaces on self,
            # then we compute the dual on each eigenspace, then put them
            # together.
            if len(self.star_eigenvalues()) == 2:
                V = self.plus_submodule(compute_dual=False).dual_free_module() + \
                    self.minus_submodule(compute_dual=False).dual_free_module()
                return V

            # At this point, we know that self is an eigenspace for star.
            V = A.sign_submodule(self.sign()).dual_free_module()
        else:
            V = A.free_module()

        N = self.level()
        p = 2
        if bound is None:
            bound = A.hecke_bound()
        while True:
            if anemic:
                while N % p == 0:
                    p = arith.next_prime(p)
            verbose("using T_%s" % p)
            f = self.hecke_polynomial(p)
            T = A.dual_hecke_matrix(p)
            V = T.kernel_on(V, poly=f, check=False)
            if V.dimension() <= self.dimension():
                break
            p = arith.next_prime(p)
            if p > bound:
                break

        if V.rank() == self.rank():
            return V

        # Failed to reduce V to the appropriate dimension
        W = self.complement()
        V2 = W.basis_matrix().right_kernel()
        if V2.rank() == self.rank():
            return V2

        raise RuntimeError("Computation of embedded dual vector space failed "
                           "(cut down to rank %s, but should have cut down to rank %s)." % (V.rank(), self.rank()))

    def free_module(self):
        """
        Return the free module corresponding to ``self``.

        EXAMPLES::

            sage: M = ModularSymbols(33,2).cuspidal_subspace() ; M
            Modular Symbols subspace of dimension 6 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field
            sage: M.free_module()
            Vector space of degree 9 and dimension 6 over Rational Field
            Basis matrix:
            [ 0  1  0  0  0  0  0 -1  1]
            [ 0  0  1  0  0  0  0 -1  1]
            [ 0  0  0  1  0  0  0 -1  1]
            [ 0  0  0  0  1  0  0 -1  1]
            [ 0  0  0  0  0  1  0 -1  1]
            [ 0  0  0  0  0  0  1 -1  0]
        """
        return self.__submodule

    def module(self):
        r"""
        Alias for \code{self.free_module()}.

        EXAMPLES::

            sage: M = ModularSymbols(17,4).cuspidal_subspace()
            sage: M.free_module() is M.module()
            True
        """
        return self.free_module()

    def intersection(self, other):
        """
        Return the intersection of ``self`` and ``other``, which must both lie in
        a common ambient space of modular symbols.

        EXAMPLES::

            sage: M = ModularSymbols(43, sign=1)
            sage: A = M[0] + M[1]
            sage: B = M[1] + M[2]
            sage: A.dimension(), B.dimension()
            (2, 3)
            sage: C = A.intersection(B); C.dimension()
            1

        TESTS::

            sage: M = ModularSymbols(1,80)
            sage: M.plus_submodule().cuspidal_submodule().sign() # indirect doctest
            1
        """
        if self.ambient_hecke_module() != other.ambient_hecke_module():
            raise ArithmeticError("intersection only defined for subspaces of"
                                  " a common ambient modular symbols space")
        if other.is_ambient():
            return self
        if self.is_ambient():
            return other

        # Neither is ambient
        V = self.free_module().intersection(other.free_module())
        M = self.ambient_hecke_module().submodule(V, check=False)

        # if sign is nonzero, the intersection will be, too
        # this only makes sense for modular symbols spaces (and hence shouldn't really be in this file)
        try:
            if self.sign():
                M._set_sign(self.sign())
            elif other.sign():
                M._set_sign(other.sign())
        except AttributeError:
            pass

        return M

    def is_ambient(self):
        r"""
        Return ``True`` if ``self`` is an ambient space of modular symbols.

        EXAMPLES::

            sage: M = ModularSymbols(17,4)
            sage: M.cuspidal_subspace().is_ambient()
            False
            sage: A = M.ambient_hecke_module()
            sage: S = A.submodule(A.basis())
            sage: sage.modular.hecke.submodule.HeckeSubmodule.is_ambient(S)
            True
        """
        return self.free_module() == self.ambient_hecke_module().free_module()

    def is_new(self, p=None):
        """
        Return ``True`` if this Hecke module is `p`-new. If `p` is None,
        returns ``True`` if it is new.

        EXAMPLES::

            sage: M = ModularSymbols(1,16)
            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.cuspidal_submodule().free_module())
            sage: S.is_new()
            True
        """
        try:
            return self.__is_new[p]
        except AttributeError:
            self.__is_new = {}
        except KeyError:
            pass
        N = self.ambient_hecke_module().new_submodule(p)
        self.__is_new[p] = self.is_submodule(N)
        return self.__is_new[p]

    def is_old(self, p=None):
        """
        Return ``True`` if this Hecke module is `p`-old. If `p` is ``None``,
        returns ``True`` if it is old.

        EXAMPLES::

            sage: M = ModularSymbols(50,2)
            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.old_submodule().free_module())
            sage: S.is_old()
            True
            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.new_submodule().free_module())
            sage: S.is_old()
            False
        """
        try:
            return self.__is_old[p]
        except AttributeError:
            self.__is_old = {}
        except KeyError:
            pass
        O = self.ambient_hecke_module().old_submodule(p)
        self.__is_old[p] = self.is_submodule(O)
        return self.__is_old[p]

    def is_submodule(self, V):
        """
        Return ``True`` if and only if ``self`` is a submodule of V.

        EXAMPLES::

            sage: M = ModularSymbols(30,4)
            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.cuspidal_submodule().free_module())
            sage: S.is_submodule(M)
            True
            sage: SS = sage.modular.hecke.submodule.HeckeSubmodule(M, M.old_submodule().free_module())
            sage: S.is_submodule(SS)
            False
        """
        if not isinstance(V, module.HeckeModule_free_module):
            return False
        return self.ambient_hecke_module() == V.ambient_hecke_module() and \
            self.free_module().is_subspace(V.free_module())

    def linear_combination_of_basis(self, v):
        """
        Return the linear combination of the basis of ``self`` given
        by the entries of `v`.

        The result can be of different types, and is printed
        accordingly, depending on the type of submodule.

        EXAMPLES::

            sage: M = ModularForms(Gamma0(2),12)

            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.cuspidal_submodule().free_module())
            sage: S.basis()
            ((1, 0, 0, 0), (0, 1, 0, 0))
            sage: S.linear_combination_of_basis([3, 10])
            (3, 10, 0, 0)

            sage: S = M.cuspidal_submodule()
            sage: S.basis()
            [
            q + 252*q^3 - 2048*q^4 + 4830*q^5 + O(q^6),
            q^2 - 24*q^4 + O(q^6)
            ]
            sage: S.linear_combination_of_basis([3, 10])
            3*q + 10*q^2 + 756*q^3 - 6384*q^4 + 14490*q^5 + O(q^6)
        """
        x = self.free_module().linear_combination_of_basis(v)
        return self(x)

    def new_submodule(self, p=None):
        """
        Return the new or `p`-new submodule of this space of modular
        symbols.

        EXAMPLES::

            sage: M = ModularSymbols(20,4)
            sage: M.new_submodule()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 18 for Gamma_0(20) of weight 4 with sign 0 over Rational Field
            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.cuspidal_submodule().free_module())
            sage: S
            Rank 12 submodule of a Hecke module of level 20
            sage: S.new_submodule()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 18 for Gamma_0(20) of weight 4 with sign 0 over Rational Field
        """
        try:
            if self.__is_new[p]:
                return self
        except AttributeError:
            self.__is_new = {}
        except KeyError:
            pass

        if self.rank() == 0:
            self.__is_new[p] = True
            return self
        try:
            return self.__new_submodule[p]
        except AttributeError:
            self.__new_submodule = {}
        except KeyError:
            pass

        S = self.ambient_hecke_module().new_submodule(p)
        ns = S.intersection(self)
        if ns.rank() == self.rank():
            self.__is_new[p] = True
        ns.__is_new = {p: True}
        self.__new_submodule[p] = ns
        return ns

    def nonembedded_free_module(self):
        """
        Return the free module corresponding to ``self`` as an abstract
        free module, i.e. not as an embedded vector space.

        EXAMPLES::

            sage: M = ModularSymbols(12,6)
            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.cuspidal_submodule().free_module())
            sage: S
            Rank 14 submodule of a Hecke module of level 12
            sage: S.nonembedded_free_module()
            Vector space of dimension 14 over Rational Field
        """
        return self.free_module().nonembedded_free_module()

    def old_submodule(self, p=None):
        r"""
        Return the old or `p`-old submodule of this space of modular
        symbols.

        EXAMPLES: We compute the old and new submodules of
        `\mathbf{S}_2(\Gamma_0(33))`.

        ::

            sage: M = ModularSymbols(33); S = M.cuspidal_submodule(); S
            Modular Symbols subspace of dimension 6 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field
            sage: S.old_submodule()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field
            sage: S.new_submodule()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field
        """
        try:
            if self.__is_old[p]:
                return self
        except AttributeError:
            self.__is_old = {}
        except KeyError:
            pass

        if self.rank() == 0:
            self.__is_old[p] = True
            return self
        try:
            return self.__old_submodule[p]
        except AttributeError:
            self.__old_submodule = {}
        except KeyError:
            pass

        S = self.ambient_hecke_module().old_submodule(p)
        os = S.intersection(self)
        if os.rank() == self.rank():
            self.__is_old[p] = True
        os.__is_old = {p: True}
        self.__old_submodule[p] = os
        return os

    def rank(self):
        r"""
        Return the rank of ``self`` as a free module over the base ring.

        EXAMPLES::

            sage: ModularSymbols(6, 4).cuspidal_subspace().rank()
            2
            sage: ModularSymbols(6, 4).cuspidal_subspace().dimension()
            2
        """
        return self.__submodule.rank()

    def submodule(self, M, Mdual=None, check=True):
        """
        Construct a submodule of ``self`` from the free module M, which
        must be a subspace of ``self``.

        EXAMPLES::

            sage: M = ModularSymbols(18,4)
            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.cuspidal_submodule().free_module())
            sage: S[0]
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 18 for Gamma_0(18) of weight 4 with sign 0 over Rational Field
            sage: S.submodule(S[0].free_module())
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 18 for Gamma_0(18) of weight 4 with sign 0 over Rational Field
        """
        if not isinstance(M, FreeModule_generic):
            V = self.ambient_module().free_module()
            if isinstance(M, (list, tuple)):
                M = V.span([V(x.element()) for x in M])
            else:
                M = V.span(M)

        if check:
            if not M.is_submodule(self.free_module()):
                raise TypeError("M (=%s) must be a submodule of the free module (=%s) associated to this module." % (M, self.free_module()))

        return self.ambient().submodule(M, Mdual, check=check)

    def submodule_from_nonembedded_module(self, V, Vdual=None, check=True):
        """
        Construct a submodule of ``self`` from V. Here V should be a
        subspace of a vector space whose dimension is the same as that
        of ``self``.

        INPUT:

        - ``V`` -- submodule of ambient free module of the same
          rank as the rank of ``self``

        - ``check`` -- whether to check that V is Hecke
          equivariant

        OUTPUT: Hecke submodule of self

        EXAMPLES::

            sage: M = ModularSymbols(37,2)
            sage: S = sage.modular.hecke.submodule.HeckeSubmodule(M, M.cuspidal_submodule().free_module())
            sage: V = (QQ**4).subspace([[1,-1,0,1/2],[0,0,1,-1/2]])
            sage: S.submodule_from_nonembedded_module(V)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field
        """
        E = self.free_module()
        M_V = V.matrix()
        M_E = E.matrix()
        # We encode the operation of taking the linear combinations of
        # the basis of E given by the basis of V as a single matrix
        # multiplication, since matrix multiplication is (presumed to be)
        # so fast, and their are asymptotically fast algorithms.
        A = M_V * M_E
        V = A.row_space()
        if Vdual is not None:
            E = self.dual_free_module()
            M_Vdual = Vdual.matrix()
            M_E = E.matrix()
            A = M_Vdual * M_E
            Vdual = A.row_space()
        return self.ambient_hecke_module().submodule(V, Vdual, check=check)

    def hecke_bound(self):
        r"""
        Compute the Hecke bound for ``self``.

        This is a number `n` such that the `T_m` for `m \leq n`
        generate the Hecke algebra.

        EXAMPLES::

            sage: M = ModularSymbols(24,8)
            sage: M.hecke_bound()
            53
            sage: M.cuspidal_submodule().hecke_bound()
            32
            sage: M.eisenstein_submodule().hecke_bound()
            53
        """
        if self.is_cuspidal():
            return self.sturm_bound()
        else:
            return self.ambient_hecke_module().hecke_bound()
