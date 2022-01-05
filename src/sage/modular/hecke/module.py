"""
Hecke modules
"""

# ****************************************************************************
#       Copyright (C) 2004,2005,2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.all import ZZ, QQ, CommutativeRing
import sage.arith.all as arith
from sage.misc.verbose import verbose
import sage.modules.module
from sage.structure.all import Sequence
import sage.matrix.matrix_space as matrix_space

import sage.misc.prandom as random

from . import algebra
from . import element
from . import hecke_operator

from sage.modules.all import FreeModule


def is_HeckeModule(x):
    r"""
    Return ``True`` if ``x`` is a Hecke module.

    EXAMPLES::

        sage: from sage.modular.hecke.module import is_HeckeModule
        sage: is_HeckeModule(ModularForms(Gamma0(7), 4))
        True
        sage: is_HeckeModule(QQ^3)
        False
        sage: is_HeckeModule(J0(37).homology())
        True
    """
    return isinstance(x, HeckeModule_generic)


class HeckeModule_generic(sage.modules.module.Module):
    r"""
    A very general base class for Hecke modules.

    We define a Hecke module of weight `k` to be a module over a commutative
    ring equipped with an action of operators `T_m` for all positive integers `m`
    coprime to some integer `n`(the level), which satisfy `T_r T_s = T_{rs}` for
    `r,s` coprime, and for powers of a prime `p`, `T_{p^r} = T_{p} T_{p^{r-1}} -
    \varepsilon(p) p^{k-1} T_{p^{r-2}}`, where `\varepsilon(p)` is some
    endomorphism of the module which commutes with the `T_m`.

    We distinguish between *full* Hecke modules, which also have an action of
    operators `T_m` for `m` not assumed to be coprime to the level, and
    *anemic* Hecke modules, for which this does not hold.
    """

    Element = element.HeckeModuleElement

    def __init__(self, base_ring, level, category=None):
        r"""
        Create a Hecke module. Not intended to be called directly.

        EXAMPLES::

            sage: CuspForms(Gamma0(17),2) # indirect doctest
            Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(17) of weight 2 over Rational Field
            sage: ModularForms(3, 3).category()
            Category of Hecke modules over Rational Field
        """
        if not isinstance(base_ring, CommutativeRing):
            raise TypeError("base_ring must be commutative ring")

        from sage.categories.hecke_modules import HeckeModules
        default_category = HeckeModules(base_ring)
        if category is None:
            category = default_category
        else:
            assert category.is_subcategory(default_category), "%s is not a subcategory of %s" % (category, default_category)

        sage.modules.module.Module.__init__(self, base_ring, category=category)

        level = ZZ(level)
        if level <= 0:
            raise ValueError("level (=%s) must be positive" % level)
        self.__level = level
        self._hecke_matrices = {}
        self._diamond_matrices = {}

    def __setstate__(self, state):
        r"""
        Ensure that the category is initialized correctly on unpickling.

        EXAMPLES::

            sage: loads(dumps(ModularSymbols(11))).category() # indirect doctest
            Category of Hecke modules over Rational Field
        """
        if not self._is_category_initialized():
            from sage.categories.hecke_modules import HeckeModules
            self._init_category_(HeckeModules(state['_base']))
        sage.modules.module.Module.__setstate__(self, state)

    def __hash__(self):
        r"""
        The hash is determined by the base ring and the level.

        EXAMPLES::

            sage: MS = sage.modular.hecke.module.HeckeModule_generic(QQ,1)
            sage: hash(MS) == hash((MS.base_ring(), MS.level()))
            True
        """
        return hash((self.base_ring(), self.__level))

    def _compute_hecke_matrix_prime_power(self, p, r, **kwds):
        r"""
        Compute the Hecke matrix T_{p^r}, where `p` is prime and `r \ge 2`, assuming that
        `T_p` is known.

        This is carried out by recursion.

        All derived classes must override either this function or
        ``self.character()``.

        EXAMPLES::

            sage: M = ModularForms(SL2Z, 24)
            sage: M._compute_hecke_matrix_prime_power(3, 3)
            [                -4112503986561480              53074162446443642880                                 0]
            [                    2592937954080                 -1312130996155080                                 0]
            [                                0                                 0 834385168339943471891603972970040]
        """
        # convert input arguments to int's.
        p, r = (int(p), int(r))
        if not arith.is_prime(p):
            raise ArithmeticError("p must be a prime")
        # T_{p^r} := T_p * T_{p^{r-1}} - eps(p)p^{k-1} T_{p^{r-2}}.
        pow = p**(r - 1)
        if pow not in self._hecke_matrices:
            # The following will force computation of T_{p^s}
            # for all s<=r-1, except possibly s=0.
            self._hecke_matrices[pow] = self._compute_hecke_matrix(pow)
        if 1 not in self._hecke_matrices:
            self._hecke_matrices[1] = self._compute_hecke_matrix(1)
        Tp = self._hecke_matrices[p]
        Tpr1 = self._hecke_matrices[pow]
        eps = self.character()
        if eps is None:
            raise NotImplementedError("either character or _compute_hecke_matrix_prime_power must be overloaded in a derived class")
        k = self.weight()
        Tpr2 = self._hecke_matrices[pow // p]
        return Tp * Tpr1 - eps(p) * (p**(k - 1)) * Tpr2

    def _compute_hecke_matrix_general_product(self, F, **kwds):
        r"""
        Compute the matrix of a general Hecke operator acting on this space, by
        factorising n into prime powers and multiplying together the Hecke
        operators for each of these.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma0(3), 4)
            sage: M._compute_hecke_matrix_general_product(factor(10))
            [1134    0]
            [   0 1134]
        """
        prod = None
        for p, r in F:
            pow = int(p**r)
            if pow not in self._hecke_matrices:
                self._hecke_matrices[pow] = self._compute_hecke_matrix(pow)
            if prod is None:
                prod = self._hecke_matrices[pow]
            else:
                prod *= self._hecke_matrices[pow]
        return prod

    def _compute_dual_hecke_matrix(self, n):
        r"""
        Compute the matrix of the Hecke operator `T_n` acting on the dual of ``self``.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma0(3), 4)
            sage: M._compute_dual_hecke_matrix(10)
            [1134    0]
            [   0 1134]
        """
        return self.hecke_matrix(n).transpose()

    def _compute_hecke_matrix(self, n, **kwds):
        r"""
        Compute the matrix of the Hecke operator `T_n` acting on ``self``.

        EXAMPLES::

            sage: M = EisensteinForms(DirichletGroup(3).0, 3)
            sage: M._compute_hecke_matrix(16)
            [205   0]
            [  0 205]
        """
        n = int(n)
        if n < 1:
            raise ValueError("Hecke operator T_%s is not defined." % n)
        if n == 1:
            Mat = matrix_space.MatrixSpace(self.base_ring(), self.rank())
            return Mat(1)

        if arith.is_prime(n):
            return self._compute_hecke_matrix_prime(n, **kwds)

        F = arith.factor(n)
        if len(F) == 1:  # nontrivial prime power case
            return self._compute_hecke_matrix_prime_power(F[0][0], F[0][1], **kwds)

        else:
            return self._compute_hecke_matrix_general_product(F, **kwds)

    def _compute_hecke_matrix_prime(self, p, **kwds):
        """
        Compute and return the matrix of the p-th Hecke operator for p prime.

        Derived classes should overload this function, and they will inherit
        the machinery for calculating general Hecke operators.

        EXAMPLES::

            sage: M = EisensteinForms(DirichletGroup(3).0, 3)
            sage: sage.modular.hecke.module.HeckeModule_generic._compute_hecke_matrix_prime(M, 3)
            Traceback (most recent call last):
            ...
            NotImplementedError: All subclasses must implement _compute_hecke_matrix_prime
        """
        raise NotImplementedError("All subclasses must implement _compute_hecke_matrix_prime")

    def _compute_diamond_matrix(self, d):
        r"""
        Compute the matrix of the diamond bracket operator `\langle d \rangle` on this space,
        in cases where this is not self-evident (i.e. when this is not a space
        with fixed character).

        EXAMPLES::

            sage: M = EisensteinForms(Gamma1(5), 3)
            sage: sage.modular.hecke.module.HeckeModule_generic._compute_diamond_matrix(M, 2)
            Traceback (most recent call last):
            ...
            NotImplementedError: All subclasses without fixed character must implement _compute_diamond_matrix
        """
        raise NotImplementedError("All subclasses without fixed character must implement _compute_diamond_matrix")

    def _hecke_operator_class(self):
        """
        Return the class to be used for instantiating Hecke operators
        acting on ``self``.

        EXAMPLES::

            sage: sage.modular.hecke.module.HeckeModule_generic(QQ,1)._hecke_operator_class()
            <class 'sage.modular.hecke.hecke_operator.HeckeOperator'>
            sage: ModularSymbols(1,12)._hecke_operator_class()
            <class 'sage.modular.modsym.hecke_operator.HeckeOperator'>
        """
        return hecke_operator.HeckeOperator

    def _diamond_operator_class(self):
        r"""
        Return the class to be used for instantiating diamond bracket operators
        acting on ``self``.

        EXAMPLES::

            sage: sage.modular.hecke.module.HeckeModule_generic(QQ,1)._diamond_operator_class()
            <class 'sage.modular.hecke.hecke_operator.DiamondBracketOperator'>
            sage: ModularSymbols(1,12)._diamond_operator_class()
            <class 'sage.modular.hecke.hecke_operator.DiamondBracketOperator'>
        """
        return hecke_operator.DiamondBracketOperator

    def anemic_hecke_algebra(self):
        """
        Return the Hecke algebra associated to this Hecke module.

        EXAMPLES::

            sage: T = ModularSymbols(1,12).hecke_algebra()
            sage: A = ModularSymbols(1,12).anemic_hecke_algebra()
            sage: T == A
            False
            sage: A
            Anemic Hecke algebra acting on Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
            sage: A.is_anemic()
            True
        """
        return algebra.AnemicHeckeAlgebra(self)

    def character(self):
        r"""
        Return the character of this space.

        As this is an abstract base class, return ``None``.

        EXAMPLES::

            sage: sage.modular.hecke.module.HeckeModule_generic(QQ, 10).character() is None
            True
        """
        return None

    def dimension(self):
        r"""
        Synonym for :meth:`rank`.

        EXAMPLES::

            sage: M = sage.modular.hecke.module.HeckeModule_generic(QQ, 10).dimension()
            Traceback (most recent call last):
            ...
            NotImplementedError: Derived subclasses must implement rank
        """
        return self.rank()

    def hecke_algebra(self):
        """
        Return the Hecke algebra associated to this Hecke module.

        EXAMPLES::

            sage: T = ModularSymbols(Gamma1(5),3).hecke_algebra()
            sage: T
            Full Hecke algebra acting on Modular Symbols space of dimension 4 for Gamma_1(5) of weight 3 with sign 0 over Rational Field
            sage: T.is_anemic()
            False

        ::

            sage: M = ModularSymbols(37,sign=1)
            sage: E, A, B = M.decomposition()
            sage: A.hecke_algebra() == B.hecke_algebra()
            False
        """
        return algebra.HeckeAlgebra(self)

    def is_zero(self):
        """
        Return ``True`` if this Hecke module has dimension 0.

        EXAMPLES::

            sage: ModularSymbols(11).is_zero()
            False
            sage: ModularSymbols(11).old_submodule().is_zero()
            True
            sage: CuspForms(10).is_zero()
            True
            sage: CuspForms(1,12).is_zero()
            False
        """
        return self.dimension() == 0

    def is_full_hecke_module(self):
        """
        Return ``True`` if this space is invariant under all Hecke operators.

        Since ``self`` is guaranteed to be an anemic Hecke module, the
        significance of this function is that it also ensures
        invariance under Hecke operators of index that divide the level.

        EXAMPLES::

            sage: M = ModularSymbols(22); M.is_full_hecke_module()
            True
            sage: M.submodule(M.free_module().span([M.0.list()]), check=False).is_full_hecke_module()
            False
        """
        try:
            return self._is_full_hecke_module
        except AttributeError:
            pass

        # now compute whether invariant under Hecke operators of index
        # dividing the level
        verbose("Determining if Hecke module is full.")
        N = self.level()
        for p in arith.prime_divisors(N):
            if not self.is_hecke_invariant(p):
                self._is_full_hecke_module = False
                return False
        self._is_full_hecke_module = True
        return True

    def is_hecke_invariant(self, n):
        """
        Return ``True`` if ``self`` is invariant under the Hecke operator `T_n`.

        Since ``self`` is guaranteed to be an anemic Hecke module it is only
        interesting to call this function when `n` is not coprime
        to the level.

        EXAMPLES::

            sage: M = ModularSymbols(22).cuspidal_subspace()
            sage: M.is_hecke_invariant(2)
            True

        We use ``check=False`` to create a nasty "module" that is not invariant
        under `T_2`::

            sage: S = M.submodule(M.free_module().span([M.0.list()]), check=False); S
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 7 for Gamma_0(22) of weight 2 with sign 0 over Rational Field
            sage: S.is_hecke_invariant(2)
            False
            sage: [n for n in range(1,12) if S.is_hecke_invariant(n)]
            [1, 3, 5, 7, 9, 11]
        """
        if arith.gcd(n, self.level()) == 1:
            return True
        if self.is_ambient():
            return True
        try:
            self.hecke_operator(n).matrix()
        except ArithmeticError:
            return False
        return True

    def level(self):
        """
        Return the level of this modular symbols space.

        INPUT:

        -  ``ModularSymbols self`` - an arbitrary space of modular symbols

        OUTPUT:

        -  ``int`` - the level

        EXAMPLES::

            sage: m = ModularSymbols(20)
            sage: m.level()
            20
        """
        return self.__level

    def rank(self):
        r"""
        Return the rank of this module over its base ring.

        This raises a ``NotImplementedError``, since this is an
        abstract base class.

        EXAMPLES::

            sage: sage.modular.hecke.module.HeckeModule_generic(QQ, 10).rank()
            Traceback (most recent call last):
            ...
            NotImplementedError: Derived subclasses must implement rank
        """
        raise NotImplementedError("Derived subclasses must implement rank")

    def submodule(self, X):
        r"""
        Return the submodule of ``self`` corresponding to ``X``.

        As this is an abstract base class, this raises a
        ``NotImplementedError``.

        EXAMPLES::

            sage: sage.modular.hecke.module.HeckeModule_generic(QQ, 10).submodule(0)
            Traceback (most recent call last):
            ...
            NotImplementedError: Derived subclasses should implement submodule
        """
        raise NotImplementedError("Derived subclasses should implement submodule")


class HeckeModule_free_module(HeckeModule_generic):
    """
    A Hecke module modeled on a free module over a commutative ring.
    """
    def __init__(self, base_ring, level, weight, category=None):
        r"""
        Initialise a module.

        EXAMPLES::

            sage: M = sage.modular.hecke.module.HeckeModule_free_module(QQ, 12, -4); M
            <class 'sage.modular.hecke.module.HeckeModule_free_module_with_category'>
            sage: skipped = ["_test_additive_associativity",
            ....:    "_test_an_element", "_test_elements",
            ....:    "_test_elements_eq_reflexive",
            ....:    "_test_elements_eq_symmetric",
            ....:    "_test_elements_eq_transitive", "_test_elements_neq",
            ....:    "_test_pickling", "_test_some_elements",
            ....:    "_test_zero", "_test_eq"]
            sage: TestSuite(M).run(skip=skipped)

        .. NOTE:: Is this supposed to be an abstract parent without elements?
        """
        HeckeModule_generic.__init__(self, base_ring, level, category=category)
        self.__weight = weight

    def _repr_(self):
        r"""

        EXAMPLES::

            sage: M = sage.modular.hecke.module.HeckeModule_free_module(QQ, 12, -4); M
            <class 'sage.modular.hecke.module.HeckeModule_free_module_with_category'>

        .. TODO::

            Implement a nicer repr, or implement the methods required
            by :class:`ModulesWithBasis` to benefit from
            :meth:`ModulesWithBasis.ParentMethods._repr_`.
        """
        return repr(type(self))

    def __getitem__(self, n):
        r"""
        Return the nth term in the decomposition of ``self``.

        See the docstring for :meth:`decomposition` for further information.

        EXAMPLES::

            sage: ModularSymbols(22)[0]
            Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 7 for Gamma_0(22) of weight 2 with sign 0 over Rational Field
        """
        n = int(n)
        D = self.decomposition()
        if n < 0 or n >= len(D):
            raise IndexError("index (=%s) must be between 0 and %s" % (n, len(D) - 1))
        return D[n]

    def __hash__(self):
        r"""
        The hash is determined by the weight, the level and the base ring.

        EXAMPLES::

            sage: MS = ModularSymbols(22)
            sage: hash(MS) == hash((MS.weight(), MS.level(), MS.base_ring()))
            True

        """
        return hash((self.__weight, self.level(), self.base_ring()))

    def __len__(self):
        r"""
        Return the number of factors in the decomposition of ``self``.

        EXAMPLES::

            sage: len(ModularSymbols(22))
            2
        """
        return len(self.decomposition())

    def _eigen_nonzero(self):
        """
        Return the smallest integer `i` such that the `i`-th entries of
        the entries of a basis for the dual vector space are not all 0.

        EXAMPLES::

            sage: M = ModularSymbols(31,2)
            sage: M._eigen_nonzero()
            0
            sage: M.dual_free_module().basis()
            [
            (1, 0, 0, 0, 0),
            (0, 1, 0, 0, 0),
            (0, 0, 1, 0, 0),
            (0, 0, 0, 1, 0),
            (0, 0, 0, 0, 1)
            ]
            sage: M.cuspidal_submodule().minus_submodule()._eigen_nonzero()
            1
            sage: M.cuspidal_submodule().minus_submodule().dual_free_module().basis()
            [
            (0, 1, 0, 0, 0),
            (0, 0, 1, 0, 0)
            ]
        """
        try:
            return self.__eigen_nonzero
        except AttributeError:
            pass
        V = self.dual_free_module()
        B = V.basis()
        for i in range(V.degree()):
            for b in B:
                if b[i] != 0:
                    self.__eigen_nonzero = i
                    return i
        assert False, 'bug in _eigen_nonzero'

    def _eigen_nonzero_element(self, n=1):
        r"""
        Return `T_n(x)` where `x` is a sparse modular
        symbol such that the image of `x` is nonzero under the dual
        projection map associated to this space, and `T_n` is the
        `n^{th}` Hecke operator.

        Used in the :meth:`dual_eigenvector` and :meth:`eigenvalue` methods.

        EXAMPLES::

            sage: ModularSymbols(22)._eigen_nonzero_element(3)
            4*(1,0) + (2,21) - (11,1) + (11,2)
        """
        if self.rank() == 0:
            raise ArithmeticError("the rank of self must be positive")
        A = self.ambient_hecke_module()
        i = self._eigen_nonzero()
        return A._hecke_image_of_ith_basis_vector(n, i)

    def _hecke_image_of_ith_basis_vector(self, n, i):
        r"""
        Return `T_n(e_i)`, where `e_i` is the
        `i`th basis vector of the ambient space.

        EXAMPLES::

            sage: ModularSymbols(Gamma0(3))._hecke_image_of_ith_basis_vector(4, 0)
            7*(1,0)
            sage: ModularForms(Gamma0(3))._hecke_image_of_ith_basis_vector(4, 0)
            7 + 84*q + 252*q^2 + 84*q^3 + 588*q^4 + 504*q^5 + O(q^6)
        """
        T = self.hecke_operator(n)
        return T.apply_sparse(self.gen(i))

    def _element_eigenvalue(self, x, name='alpha'):
        r"""
        Return the dot product of ``self`` with the eigenvector returned by dual_eigenvector.

        EXAMPLES::

            sage: M = ModularSymbols(11)[0]
            sage: M._element_eigenvalue(M.0)
            1
        """
        if not element.is_HeckeModuleElement(x):
            raise TypeError("x must be a Hecke module element.")
        if x not in self.ambient_hecke_module():
            raise ArithmeticError("x must be in the ambient Hecke module.")
        v = self.dual_eigenvector(names=name)
        return v.dot_product(x.element())

    def _is_hecke_equivariant_free_module(self, submodule):
        """
        Return ``True`` if the given free submodule of the ambient free module
        is invariant under all Hecke operators.

        EXAMPLES::

            sage: M = ModularSymbols(11); V = M.free_module()
            sage: M._is_hecke_equivariant_free_module(V.span([V.0]))
            False
            sage: M._is_hecke_equivariant_free_module(V)
            True
            sage: M._is_hecke_equivariant_free_module(M.cuspidal_submodule().free_module())
            True

        We do the same as above, but with a modular forms space::

            sage: M = ModularForms(11); V = M.free_module()
            sage: M._is_hecke_equivariant_free_module(V.span([V.0 + V.1]))
            False
            sage: M._is_hecke_equivariant_free_module(V)
            True
            sage: M._is_hecke_equivariant_free_module(M.cuspidal_submodule().free_module())
            True
        """
        verbose("Determining if free module is Hecke equivariant.")
        bound = self.hecke_bound()
        for p in arith.primes(bound + 1):
            try:
                self.T(p).matrix().restrict(submodule, check=True)
            except ArithmeticError:
                return False
        return True

    def _set_factor_number(self, i):
        r"""
        For internal use. If this Hecke module was computed via a decomposition of another
        Hecke module, this method stores the index of this space in that decomposition.

        EXAMPLES::

            sage: ModularSymbols(Gamma0(3))[0].factor_number() # indirect doctest
            0
        """
        self.__factor_number = i

    def ambient(self):
        r"""
        Return the ambient module associated to this module.

        Synonym for :meth:`ambient_hecke_module`.

        EXAMPLES::

            sage: CuspForms(1, 12).ambient()
            Modular Forms space of dimension 2 for Modular Group SL(2,Z) of weight 12 over Rational Field
        """
        return self.ambient_hecke_module()

    def ambient_module(self):
        r"""
        Return the ambient module associated to this module.

        Synonym for :meth:`ambient_hecke_module`.

        EXAMPLES::

            sage: CuspForms(1, 12).ambient_module()
            Modular Forms space of dimension 2 for Modular Group SL(2,Z) of weight 12 over Rational Field
            sage: sage.modular.hecke.module.HeckeModule_free_module(QQ, 10, 3).ambient_module()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return self.ambient_hecke_module()

    def ambient_hecke_module(self):
        r"""
        Return the ambient module associated to this module.

        As this is an abstract base class, raise ``NotImplementedError``.

        EXAMPLES::

            sage: sage.modular.hecke.module.HeckeModule_free_module(QQ, 10, 3).ambient_hecke_module()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def atkin_lehner_operator(self, d=None):
        r"""
        Return the Atkin-Lehner operator `W_d` on this space, if defined, where
        `d` is a divisor of the level `N` such that `N/d` and `d` are coprime.
        If `d` is not given, we take `d = N`.  If `N/d` is not coprime to `d`,
        then we replace `d` with the unique integer having this property which
        has the same prime factors as `d`.

        .. NOTE::

            The operator `W_d` is given by the action of any matrix of the form

            .. math::

                W_d = \begin{pmatrix} dx & y \\ Nz & dw \end{pmatrix}

            with `\det W_d = d` and such that `x = 1 \bmod N/d`, `y = 1 \bmod
            d`, as in [AL1978]_. However, our definition of the weight `k`
            action differs from theirs by a power of the determinant, so our
            operator `W_d` is `d^{k/2 - 1}` times the operator of Atkin-Li. In
            particular, if `k = 2` our conventions are identical to Atkin and
            Li's.

            With Sage's conventions, the operator `W_d` satisfies

            .. math::

                W_d^2 = d^{k - 2} \langle x^{-1} \rangle

            where `x` is congruent to `d` modulo `N/d` and to `-1` modulo `d`.
            In particular, the operator is an involution in weight 2 and
            trivial character (but not in most other situations).

        EXAMPLES::

            sage: M = ModularSymbols(11)
            sage: w = M.atkin_lehner_operator()
            sage: w
            Hecke module morphism Atkin-Lehner operator W_11 defined by the matrix
            [-1  0  0]
            [ 0 -1  0]
            [ 0  0 -1]
            Domain: Modular Symbols space of dimension 3 for Gamma_0(11) of weight ...
            Codomain: Modular Symbols space of dimension 3 for Gamma_0(11) of weight ...
            sage: M = ModularSymbols(Gamma1(13))
            sage: w = M.atkin_lehner_operator()
            sage: w.fcp('x')
            (x - 1)^7 * (x + 1)^8

        ::

            sage: M = ModularSymbols(33)
            sage: S = M.cuspidal_submodule()
            sage: S.atkin_lehner_operator()
            Hecke module morphism Atkin-Lehner operator W_33 defined by the matrix
            [ 0 -1  0  1 -1  0]
            [ 0 -1  0  0  0  0]
            [ 0 -1  0  0 -1  1]
            [ 1 -1  0  0 -1  0]
            [ 0  0  0  0 -1  0]
            [ 0 -1  1  0 -1  0]
            Domain: Modular Symbols subspace of dimension 6 of Modular Symbols space ...
            Codomain: Modular Symbols subspace of dimension 6 of Modular Symbols space ...

        ::

            sage: S.atkin_lehner_operator(3)
            Hecke module morphism Atkin-Lehner operator W_3 defined by the matrix
            [ 0  1  0 -1  1  0]
            [ 0  1  0  0  0  0]
            [ 0  1  0  0  1 -1]
            [-1  1  0  0  1  0]
            [ 0  0  0  0  1  0]
            [ 0  1 -1  0  1  0]
            Domain: Modular Symbols subspace of dimension 6 of Modular Symbols space ...
            Codomain: Modular Symbols subspace of dimension 6 of Modular Symbols space ...

        ::

            sage: N = M.new_submodule()
            sage: N.atkin_lehner_operator()
            Hecke module morphism Atkin-Lehner operator W_33 defined by the matrix
            [  1 2/5 4/5]
            [  0  -1   0]
            [  0   0  -1]
            Domain: Modular Symbols subspace of dimension 3 of Modular Symbols space ...
            Codomain: Modular Symbols subspace of dimension 3 of Modular Symbols space ...
        """
        if d is None:
            d = self.level()
        d = int(d)
        if self.level() % d:
            raise ArithmeticError("d (=%s) must be a divisor of the level (=%s)" % (d, self.level()))

        N = self.level()
        for p, e in arith.factor(d):
            v = arith.valuation(N, p)
            if e < v:
                d *= p**(v - e)
        d = int(d)
        try:
            return self.__atkin_lehner_operator[d]
        except AttributeError:
            self.__atkin_lehner_operator = {}
        except KeyError:
            pass
        Wmat = self._compute_atkin_lehner_matrix(d)
        H = self.endomorphism_ring()
        W = H(Wmat, "Atkin-Lehner operator W_%s" % d)
        self.__atkin_lehner_operator[d] = W
        return W

    def basis(self):
        """
        Return a basis for ``self``.

        EXAMPLES::

            sage: m = ModularSymbols(43)
            sage: m.basis()
            ((1,0), (1,31), (1,32), (1,38), (1,39), (1,40), (1,41))
        """
        try:
            return self.__basis
        except AttributeError:
            self.__basis = self.gens()
        return self.__basis

    def basis_matrix(self):
        r"""
        Return the matrix of the basis vectors of ``self`` (as vectors in some
        ambient module)

        EXAMPLES::

            sage: CuspForms(1, 12).basis_matrix()
            [1 0]
        """
        return self.free_module().basis_matrix()

    def coordinate_vector(self, x):
        """
        Write ``x`` as a vector with respect to the basis given by
        self.basis().

        EXAMPLES::

            sage: S = ModularSymbols(11,2).cuspidal_submodule()
            sage: S.0
            (1,8)
            sage: S.basis()
            ((1,8), (1,9))
            sage: S.coordinate_vector(S.0)
            (1, 0)
        """
        return self.free_module().coordinate_vector(x.element())

    def decomposition(self, bound=None, anemic=True, height_guess=1,
                      sort_by_basis=False, proof=None):
        """
        Return the maximal decomposition of this Hecke module under the
        action of Hecke operators of index coprime to the level.

        This is the finest decomposition of ``self`` that we can obtain
        using factors obtained by taking kernels of Hecke operators.

        Each factor in the decomposition is a Hecke submodule obtained as
        the kernel of `f(T_n)^r` acting on self, where n is
        coprime to the level and `r=1`. If anemic is False, instead
        choose `r` so that `f(X)^r` exactly divides the
        characteristic polynomial.

        INPUT:

        -  ``anemic`` - bool (default: True), if True, use only
           Hecke operators of index coprime to the level.

        -  ``bound`` - int or None, (default: None). If None,
           use all Hecke operators up to the Sturm bound, and hence obtain the
           same result as one would obtain by using every element of the Hecke
           ring. If a fixed integer, decompose using only Hecke operators
           `T_p`, with `p` prime, up to bound.
        -  ``sort_by_basis`` - bool (default: ``False``); If True the resulting
           decomposition will be sorted as if it was free modules, ignoring the
           Hecke module structure. This will save a lot of time.

        OUTPUT:

        -  ``list`` - a list of subspaces of ``self``.

        EXAMPLES::

            sage: ModularSymbols(17,2).decomposition()
            [
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(17) of weight 2 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(17) of weight 2 with sign 0 over Rational Field
            ]
            sage: ModularSymbols(Gamma1(10),4).decomposition()
            [
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 18 for Gamma_1(10) of weight 4 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 18 for Gamma_1(10) of weight 4 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 18 for Gamma_1(10) of weight 4 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 18 for Gamma_1(10) of weight 4 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 18 for Gamma_1(10) of weight 4 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 18 for Gamma_1(10) of weight 4 with sign 0 over Rational Field
            ]
            sage: ModularSymbols(GammaH(12, [11])).decomposition()
            [
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 9 for Congruence Subgroup Gamma_H(12) with H generated by [11] of weight 2 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 9 for Congruence Subgroup Gamma_H(12) with H generated by [11] of weight 2 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 9 for Congruence Subgroup Gamma_H(12) with H generated by [11] of weight 2 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 9 for Congruence Subgroup Gamma_H(12) with H generated by [11] of weight 2 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 5 of Modular Symbols space of dimension 9 for Congruence Subgroup Gamma_H(12) with H generated by [11] of weight 2 with sign 0 over Rational Field
            ]

        TESTS::

            sage: M = ModularSymbols(1000,2,sign=1).new_subspace().cuspidal_subspace()
            sage: M.decomposition(3, sort_by_basis = True)
            [
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 154 for Gamma_0(1000) of weight 2 with sign 1 over Rational Field,
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 154 for Gamma_0(1000) of weight 2 with sign 1 over Rational Field,
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 154 for Gamma_0(1000) of weight 2 with sign 1 over Rational Field,
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 154 for Gamma_0(1000) of weight 2 with sign 1 over Rational Field,
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 154 for Gamma_0(1000) of weight 2 with sign 1 over Rational Field,
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 154 for Gamma_0(1000) of weight 2 with sign 1 over Rational Field
            ]
        """
        if not isinstance(anemic, bool):
            raise TypeError("anemic must be of type bool.")

        key = (bound, anemic)

        try:
            if self.__decomposition[key] is not None:
                return self.__decomposition[key]
        except AttributeError:
            self.__decomposition = {}
        except KeyError:
            pass
        if self.rank() == 0:
            self.__decomposition[key] = Sequence([], immutable=True, cr=True)
            return self.__decomposition[key]

        is_rational = self.base_ring() == QQ

        time = verbose("Decomposing %s" % self)
        T = self.ambient_hecke_module().hecke_algebra()
        if bound is None:
            bound = self.ambient_hecke_module().hecke_bound()
        D = Sequence([], cr=True)
        U = [self.free_module()]
        p = 2
        while U and p <= bound:
            verbose(mesg="p=%s" % p, t=time)
            if anemic:
                while arith.GCD(p, self.level()) != 1:
                    p = arith.next_prime(p)
            verbose("Decomposition using p=%s" % p)
            t = T.hecke_operator(p).matrix()
            Uprime = []
            for i in range(len(U)):
                is_diagonalizable = (not self.base_ring().characteristic() and
                                     self.level() % p)
                if is_rational:
                    X = t.decomposition_of_subspace(U[i], check_restrict=False,
                                                    algorithm='multimodular',
                                                    height_guess=height_guess, proof=proof)
                else:
                    X = t.decomposition_of_subspace(U[i], check_restrict=False,
                                                    is_diagonalizable=is_diagonalizable)
                for Xi in X:
                    W, is_irred = Xi
                    if is_irred:
                        A = self.submodule(W, check=False)
                        D.append(A)
                    else:
                        Uprime.append(W)
            # end for
            p = arith.next_prime(p)
            U = Uprime
        # end while
        for i in range(len(U)):
            A = self.submodule(U[i], check=False)
            D.append(A)
        for A in D:
            if anemic:
                A.__is_splittable_anemic = False
                A.__is_splittable = False
            else:
                A.__is_splittable = False
        self.__is_splittable = len(D) > 1
        if anemic:
            self.__is_splittable_anemic = len(D) > 1
        from sage.modules.free_module import EchelonMatrixKey
        D.sort(key=None if not sort_by_basis
               else lambda ss: EchelonMatrixKey(ss.free_module()))
        D.set_immutable()
        self.__decomposition[key] = D
        for i in range(len(D)):
            self.__decomposition[key][i]._set_factor_number(i)
        return self.__decomposition[key]

    def degree(self):
        r"""
        Return the degree of this Hecke module.

        This is the rank of the ambient free module.

        EXAMPLES::

            sage: CuspForms(1, 12).degree()
            2
        """
        return self.free_module().degree()

    def dual_eigenvector(self, names='alpha', lift=True, nz=None):
        """
        Return an eigenvector for the Hecke operators acting on the linear
        dual of this space.

        This eigenvector will have entries in an extension of the base
        ring of degree equal to the dimension of this space.

        .. WARNING::

           The input space must be simple.

        INPUT:

        -  ``name`` -- print name of generator for eigenvalue
           field.

        -  ``lift`` -- bool (default: ``True``)

        -  ``nz`` -- if not ``None``, then normalize vector so dot
           product with this basis vector of ambient space is 1.

        OUTPUT: A vector with entries possibly in an extension of the base
        ring. This vector is an eigenvector for all Hecke operators acting
        via their transpose.

        If lift = False, instead return an eigenvector in the subspace for
        the Hecke operators on the dual space. I.e., this is an eigenvector
        for the restrictions of Hecke operators to the dual space.

        .. NOTE::

           #. The answer is cached so subsequent calls always return
              the same vector. However, the algorithm is randomized,
              so calls during another session may yield a different
              eigenvector. This function is used mainly for computing
              systems of Hecke eigenvalues.

           #. One can also view a dual eigenvector as defining (via
              dot product) a functional phi from the ambient space of
              modular symbols to a field. This functional phi is an
              eigenvector for the dual action of Hecke operators on
              functionals.

        EXAMPLES::

            sage: SF = ModularSymbols(14).cuspidal_subspace().simple_factors()
            sage: sorted([u.dual_eigenvector() for u in SF])
            [(0, 1, 0, 0, 0), (1, 0, -3, 2, -1)]
        """
        # TODO -- optimize by computing the answer for i not None in terms
        # of the answer for a given i if known !!
        try:
            w, w_lift = self.__dual_eigenvector[(names, nz)]
            if lift:
                return w_lift
            else:
                return w
        except KeyError:
            pass
        except AttributeError:
            self.__dual_eigenvector = {}

        if not self.is_simple():
            raise ArithmeticError("self must be simple")

        # Find a Hecke operator that acts irreducibly on this space:
        p = 2
        t = self.dual_hecke_matrix(p)
        while True:
            f = t.charpoly('x')
            if f.is_irreducible():
                break
            p = arith.next_prime(p)
            t += random.choice([-2, -1, 1, 2]) * self.dual_hecke_matrix(p)

        # Write down the eigenvector.
        # Write f(x) = (x-alpha)*g(x), where alpha is a root
        # of f(x).
        n = f.degree()
        if n > 1:
            R = f.parent()
            K = R.base_ring().extension(f, names=names)
            alpha = K.gen()
            beta = ~alpha   # multiplicative inverse of alpha
            c = [-f[0] * beta]
            for i in range(1, n - 1):
                c.append((c[i - 1] - f[i]) * beta)
            c.append(K.one())
        else:
            K = self.base_ring()
            c = [1]

        # The entries of c are the coefficients of g (as stated in
        # William Stein's Ph.D. thesis, Section 3.5.3).  We compute
        # g(t)v for a some vector v, and get an eigenvector.
        V = FreeModule(K, n)
        t = t.change_ring(K)      # coerce t to be over K.
        for j in range(n):
            v = V.gen(j)
            I = t.iterates(v, n)  # iterates v, v*t, v*t^2, ...
            w = V(0)
            for i in range(n):
                w += c[i] * V(I.row(i).list())
            if w != 0:
                break

        # Now w is an eigenvector for the action of the Hecke
        # operators on the subspace.  We need an eigenvector
        # in the original space, so we take the linear combination
        # of the basis for the embedded dual vector space given
        # by the entries of w.
        Vdual = self.dual_free_module().change_ring(K)
        w_lift = Vdual.linear_combination_of_basis(w)

        # Finally rescale so the dot product of this vector and
        # the _eigen_nonzero_element is 1.
        if nz is not None:
            x = self.ambient().gen(nz)
        else:
            x = self._eigen_nonzero_element()
        alpha = w_lift.dot_product(x.element())
        beta = ~alpha
        w_lift = w_lift * beta
        w = w * beta

        self.__dual_eigenvector[(names, nz)] = (w, w_lift)
        if lift:
            return w_lift
        else:
            return w

    def dual_hecke_matrix(self, n):
        """
        Return the matrix of the `n^{th}` Hecke operator acting on the dual
        embedded representation of ``self``.

        EXAMPLES::

            sage: CuspForms(1, 24).dual_hecke_matrix(5)
            [     44656110        -15040]
            [-307849789440      28412910]
        """
        n = int(n)
        try:
            self._dual_hecke_matrices
        except AttributeError:
            self._dual_hecke_matrices = {}
        if n not in self._dual_hecke_matrices:
            T = self._compute_dual_hecke_matrix(n)
            self._dual_hecke_matrices[n] = T
        return self._dual_hecke_matrices[n]

    def eigenvalue(self, n, name='alpha'):
        r"""
        Assuming that ``self`` is a simple space, return the eigenvalue of the
        `n^{th}` Hecke operator on ``self``.

        INPUT:

        -  ``n`` - index of Hecke operator

        -  ``name`` - print representation of generator of
           eigenvalue field

        EXAMPLES::

            sage: A = ModularSymbols(125,sign=1).new_subspace()[0]
            sage: A.eigenvalue(7)
            -3
            sage: A.eigenvalue(3)
            -alpha - 2
            sage: A.eigenvalue(3,'w')
            -w - 2
            sage: A.eigenvalue(3,'z').charpoly('x')
            x^2 + 3*x + 1
            sage: A.hecke_polynomial(3)
            x^2 + 3*x + 1

        ::

            sage: M = ModularSymbols(Gamma1(17)).decomposition()[8].plus_submodule()
            sage: M.eigenvalue(2,'a')
            a
            sage: M.eigenvalue(4,'a')
            4/3*a^3 + 17/3*a^2 + 28/3*a + 8/3

        .. NOTE::

           #. In fact there are `d` systems of eigenvalues
              associated to self, where `d` is the rank of
              self. Each of the systems of eigenvalues is conjugate
              over the base field. This function chooses one of the
              systems and consistently returns eigenvalues from that
              system. Thus these are the coefficients `a_n` for
              `n\geq 1` of a modular eigenform attached to self.

           #. This function works even for Eisenstein subspaces,
              though it will not give the constant coefficient of one
              of the corresponding Eisenstein series (i.e., the
              generalized Bernoulli number).

        TESTS:

        This checks that :trac:`15201` is fixed::

            sage: M = ModularSymbols(5, 6, sign=1)
            sage: f = M.decomposition()[0]
            sage: f.eigenvalue(10)
            50
        """
        if not self.is_simple():
            raise ArithmeticError("self must be simple")
        n = int(n)
        try:
            return self.__eigenvalues[n][name]
        except AttributeError:
            self.__eigenvalues = {}
        except KeyError:
            pass
        if n <= 0:
            raise IndexError("n must be a positive integer")

        ev = self.__eigenvalues

        if n == 1 or arith.is_prime(n):
            Tn_e = self._eigen_nonzero_element(n)
            an = self._element_eigenvalue(Tn_e, name=name)
            _dict_set(ev, n, name, an)
            return an

        # Now use the Hecke eigenvalue recurrence, since arithmetic in
        # a field is faster than computing Heilbronn matrices for
        # non-prime n and doing some big sum (i.e., computing T_n(e)).
        # Also by computing using the recurrence on eigenvalues
        # we use information about divisors.
        F = arith.factor(n)
        prod = None
        for p, r in F:
            (p, r) = (int(p), int(r))
            pow = p**r
            if not (pow in ev and name in ev[pow]):
                # TODO: Optimization -- do something much more
                # intelligent in case character is not defined.  For
                # example, compute it using the diamond operators <d>
                eps = self.character()
                if eps is None:
                    Tn_e = self._eigen_nonzero_element(pow)
                    _dict_set(ev, pow, name, self._element_eigenvalue(Tn_e, name=name))
                else:
                    # a_{p^r} := a_p * a_{p^{r-1}} - eps(p)p^{k-1} a_{p^{r-2}}
                    ap = self.eigenvalue(p, name=name)
                    if r == 1:
                        apow = ap
                    else:
                        apr1 = self.eigenvalue(pow // p, name=name)
                        k = self.weight()
                        apr2 = self.eigenvalue(pow // (p * p), name=name)
                        apow = ap * apr1 - eps(p) * (p**(k - 1)) * apr2
                    _dict_set(ev, pow, name, apow)
            if prod is None:
                prod = ev[pow][name]
            else:
                prod *= ev[pow][name]
        _dict_set(ev, n, name, prod)
        return prod

    def factor_number(self):
        """
        If this Hecke module was computed via a decomposition of another
        Hecke module, this is the corresponding number. Otherwise return
        -1.

        EXAMPLES::

            sage: ModularSymbols(23)[0].factor_number()
            0
            sage: ModularSymbols(23).factor_number()
            -1
        """
        try:
            return self.__factor_number
        except AttributeError:
            return -1

    def gens(self):
        """
        Return a tuple of basis elements of ``self``.

        EXAMPLES::

            sage: ModularSymbols(23).gens()
            ((1,0), (1,17), (1,19), (1,20), (1,21))
        """
        return tuple(self(x) for x in self.free_module().gens())

    def gen(self, n):
        r"""
        Return the nth basis vector of the space.

        EXAMPLES::

            sage: ModularSymbols(23).gen(1)
            (1,17)
        """
        return self(self.free_module().gen(n))

    def hecke_matrix(self, n):
        """
        Return the matrix of the `n^{th}` Hecke operator acting on given
        basis.

        EXAMPLES::

            sage: C = CuspForms(1, 16)
            sage: C.hecke_matrix(3)
            [-3348]
        """
        n = int(n)
        if n <= 0:
            raise IndexError("n must be positive.")
        if n not in self._hecke_matrices:
            T = self._compute_hecke_matrix(n)
            T.set_immutable()
            self._hecke_matrices[n] = T
        return self._hecke_matrices[n]

    def hecke_operator(self, n):
        """
        Return the `n`-th Hecke operator `T_n`.

        INPUT:

        -  ``ModularSymbols self`` - Hecke equivariant space of
           modular symbols

        -  ``int n`` - an integer at least 1.

        EXAMPLES::

            sage: M = ModularSymbols(11,2)
            sage: T = M.hecke_operator(3) ; T
            Hecke operator T_3 on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: T.matrix()
            [ 4  0 -1]
            [ 0 -1  0]
            [ 0  0 -1]
            sage: T(M.0)
            4*(1,0) - (1,9)
            sage: S = M.cuspidal_submodule()
            sage: T = S.hecke_operator(3) ; T
            Hecke operator T_3 on Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: T.matrix()
            [-1  0]
            [ 0 -1]
            sage: T(S.0)
            -(1,8)
        """
        return self.hecke_algebra().hecke_operator(n)

    def diamond_bracket_matrix(self, d):
        r"""
        Return the matrix of the diamond bracket operator `\langle d \rangle` on ``self``.

        EXAMPLES::

            sage: M = ModularSymbols(DirichletGroup(5).0, 3)
            sage: M.diamond_bracket_matrix(3)
            [-zeta4      0]
            [     0 -zeta4]
            sage: ModularSymbols(Gamma1(5), 3).diamond_bracket_matrix(3)
            [ 0  1  0  0]
            [-1  0  0  0]
            [ 0  0  0  1]
            [ 0  0 -1  0]
        """
        d = int(d) % self.level()
        if d not in self._diamond_matrices:
            if self.character() is not None:
                D = matrix_space.MatrixSpace(self.base_ring(), self.rank())(self.character()(d))
            else:
                D = self._compute_diamond_matrix(d)
            D.set_immutable()
            self._diamond_matrices[d] = D
        return self._diamond_matrices[d]

    def diamond_bracket_operator(self, d):
        r"""
        Return the diamond bracket operator `\langle d \rangle` on self.

        EXAMPLES::

            sage: M = ModularSymbols(DirichletGroup(5).0, 3)
            sage: M.diamond_bracket_operator(3)
            Diamond bracket operator <3> on Modular Symbols space of dimension 2 and level 5, weight 3, character [zeta4], sign 0, over Cyclotomic Field of order 4 and degree 2
        """
        return self.hecke_algebra().diamond_bracket_operator(d)

    def T(self, n):
        r"""
        Return the `n^{th}` Hecke operator `T_n`.

        This function is a synonym for :meth:`hecke_operator`.

        EXAMPLES::

            sage: M = ModularSymbols(11,2)
            sage: M.T(3)
            Hecke operator T_3 on Modular Symbols ...
        """
        return self.hecke_operator(n)

    def hecke_polynomial(self, n, var='x'):
        """
        Return the characteristic polynomial of the `n^{th}` Hecke operator
        acting on this space.

        INPUT:

        -  ``n`` - integer

        OUTPUT: a polynomial

        EXAMPLES::

            sage: ModularSymbols(11,2).hecke_polynomial(3)
            x^3 - 2*x^2 - 7*x - 4
        """
        return self.hecke_operator(n).charpoly(var)

    def is_simple(self):
        r"""
        Return ``True`` if this space is simple as a module for the
        corresponding Hecke algebra.

        Raises ``NotImplementedError``, as this is an abstract base
        class.

        EXAMPLES::

            sage: sage.modular.hecke.module.HeckeModule_free_module(QQ, 10, 3).is_simple()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def is_splittable(self):
        """
        Return ``True`` if and only if only it is possible to split
        off a nontrivial generalized eigenspace of ``self`` as the
        kernel of some Hecke operator (not necessarily prime to the level).

        Note that the direct sum of several copies of the same simple
        module is not splittable in this sense.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma0(64)).cuspidal_subspace()
            sage: M.is_splittable()
            True
            sage: M.simple_factors()[0].is_splittable()
            False
        """
        if not hasattr(self, "__is_splittable"):
            self.decomposition(anemic=False)
        return self.__is_splittable

    def is_submodule(self, other):
        r"""
        Return ``True`` if ``self`` is a submodule of ``other``.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma0(64))
            sage: M[0].is_submodule(M)
            True
            sage: CuspForms(1, 24).is_submodule(ModularForms(1, 24))
            True
            sage: CuspForms(1, 12).is_submodule(CuspForms(3, 12))
            False
        """
        if not isinstance(other, HeckeModule_free_module):
            return False
        return (self.ambient_free_module() == other.ambient_free_module() and
                self.free_module().is_submodule(other.free_module()))

    def is_splittable_anemic(self):
        """
        Return ``True`` if and only if only it is possible to split off a
        nontrivial generalized eigenspace of ``self`` as the kernel of some
        Hecke operator of index coprime to the level.

        Note that the direct sum of several copies of the same simple
        module is not splittable in this sense.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma0(64)).cuspidal_subspace()
            sage: M.is_splittable_anemic()
            True
            sage: M.simple_factors()[0].is_splittable_anemic()
            False
        """
        if not hasattr(self, "__is_splittable_anemic"):
            self.decomposition(anemic=True)
        return self.__is_splittable_anemic

    def ngens(self):
        r"""
        Return the number of generators of ``self``.

        This is equal to the rank.

        EXAMPLES::

            sage: ModularForms(1, 12).ngens()
            2
        """
        return self.rank()

    def projection(self):
        r"""
        Return the projection map from the ambient space to ``self``.

        ALGORITHM: Let `B` be the matrix whose columns are obtained
        by concatenating together a basis for the factors of the ambient
        space. Then the projection matrix onto self is the submatrix of
        `B^{-1}` obtained from the rows corresponding to self,
        i.e., if the basis vectors for self appear as columns `n`
        through `m` of `B`, then the projection matrix is
        got from rows `n` through `m` of `B^{-1}`.
        This is because projection with respect to the B basis is just
        given by an `m-n+1` row slice `P` of a diagonal
        matrix D with 1's in the `n` through `m` positions,
        so projection with respect to the standard basis is given by
        `P\cdot B^{-1}`, which is just rows `n`
        through `m` of `B^{-1}`.

        EXAMPLES::

            sage: e = EllipticCurve('34a')
            sage: m = ModularSymbols(34); s = m.cuspidal_submodule()
            sage: d = s.decomposition(7)
            sage: d
            [
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 9 for Gamma_0(34) of weight 2 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 9 for Gamma_0(34) of weight 2 with sign 0 over Rational Field
            ]
            sage: a = d[0]; a
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 9 for Gamma_0(34) of weight 2 with sign 0 over Rational Field
            sage: pi = a.projection()
            sage: pi(m([0,oo]))
            -1/6*(2,7) + 1/6*(2,13) - 1/6*(2,31) + 1/6*(2,33)
            sage: M = ModularSymbols(53,sign=1)
            sage: S = M.cuspidal_subspace()[1] ; S
            Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 5 for Gamma_0(53) of weight 2 with sign 1 over Rational Field
            sage: p = S.projection()
            sage: S.basis()
            ((1,43) - (1,45), (1,47), (1,50))
            sage: [ p(x) for x in S.basis() ]
            [(1,43) - (1,45), (1,47), (1,50)]
            sage: all(p(x)==x for x in S.basis())
            True
        """

        # Compute the Hecke-stable projection map pi from the ambient
        # space of M to M.  Computing the projection map is the same
        # as writing the ambient space as a direct sum of M and its
        # Hecke-stable complement, which is the old subspace plus the
        # other new factors, then *inverting*.  With the projection
        # map in hand, we can compute Hecke operators directly on M
        # fairly quickly without having to compute them on the whole
        # ambient space.  Of course, computing this inverse is way too
        # much work to be useful in general (!).  (I sort of learned
        # this trick from Joe Wetherell, or at least he was aware of
        # it when I mentioned it to him in an airport once.  It's also
        # sort of like a trick Cremona uses in his book for elliptic
        # curves.)  It's not a very good trick though.

        try:
            return self.__projection
        except AttributeError:
            i = self.factor_number()
            if i == -1:
                raise NotImplementedError("Computation of projection only implemented "
                                          "for decomposition factors.")
            A = self.ambient_hecke_module()
            B = A.decomposition_matrix_inverse()
            i = A.decomposition().index(self)
            n = sum([A[j].rank() for j in range(i)])
            C = B.matrix_from_columns(range(n, n + self.rank()))
            H = A.Hom(self)
            pi = H(C, "Projection" % self)
            self.__projection = pi
            return self.__projection

    def system_of_eigenvalues(self, n, name='alpha'):
        r"""
        Assuming that ``self`` is a simple space of modular symbols, return the
        eigenvalues `[a_1, \ldots, a_nmax]` of the Hecke
        operators on self. See ``self.eigenvalue(n)`` for more
        details.

        INPUT:

        -  ``n`` - number of eigenvalues

        -  ``alpha`` - name of generate for eigenvalue field

        EXAMPLES:

        The outputs of the following tests are very unstable. The algorithms
        are randomized and depend on cached results. A slight change in the
        sequence of pseudo-random numbers or a modification in caching is
        likely to modify the results. We reset the random number generator and
        clear some caches for reproducibility::

            sage: set_random_seed(0)
            sage: ModularSymbols_clear_cache()

        We compute eigenvalues for newforms of level 62::

            sage: M = ModularSymbols(62,2,sign=-1)
            sage: S = M.cuspidal_submodule().new_submodule()
            sage: [[o.minpoly() for o in A.system_of_eigenvalues(3)] for A in S.decomposition()]
            [[x - 1, x - 1, x], [x - 1, x + 1, x^2 - 2*x - 2]]

        Next we define a function that does the above::

            sage: def b(N,k=2):
            ....:    t=cputime()
            ....:    S = ModularSymbols(N,k,sign=-1).cuspidal_submodule().new_submodule()
            ....:    for A in S.decomposition():
            ....:        print("{} {}".format(N, A.system_of_eigenvalues(5)))

        ::

            sage: b(63)
            63 [1, 1, 0, -1, 2]
            63 [1, alpha, 0, 1, -2*alpha]

        This example illustrates finding field over which the eigenvalues
        are defined::

            sage: M = ModularSymbols(23,2,sign=1).cuspidal_submodule().new_submodule()
            sage: v = M.system_of_eigenvalues(10); v
            [1, alpha, -2*alpha - 1, -alpha - 1, 2*alpha, alpha - 2, 2*alpha + 2, -2*alpha - 1, 2, -2*alpha + 2]
            sage: v[0].parent()
            Number Field in alpha with defining polynomial x^2 + x - 1

        This example illustrates setting the print name of the eigenvalue
        field.

        ::

            sage: A = ModularSymbols(125,sign=1).new_subspace()[0]
            sage: A.system_of_eigenvalues(10)
            [1, alpha, -alpha - 2, -alpha - 1, 0, -alpha - 1, -3, -2*alpha - 1, 3*alpha + 2, 0]
            sage: A.system_of_eigenvalues(10,'x')
            [1, x, -x - 2, -x - 1, 0, -x - 1, -3, -2*x - 1, 3*x + 2, 0]
        """
        return [self.eigenvalue(m, name=name) for m in range(1, n + 1)]

    def weight(self):
        """
        Return the weight of this Hecke module.

        INPUT:

        -  ``self`` - an arbitrary Hecke module

        OUTPUT:

        -  ``int`` - the weight

        EXAMPLES::

            sage: m = ModularSymbols(20, weight=2)
            sage: m.weight()
            2
        """
        return self.__weight

    def zero_submodule(self):
        """
        Return the zero submodule of ``self``.

        EXAMPLES::

            sage: ModularSymbols(11,4).zero_submodule()
            Modular Symbols subspace of dimension 0 of Modular Symbols space of dimension 6 for Gamma_0(11) of weight 4 with sign 0 over Rational Field
            sage: CuspForms(11,4).zero_submodule()
            Modular Forms subspace of dimension 0 of Modular Forms space of dimension 4 for Congruence Subgroup Gamma0(11) of weight 4 over Rational Field
        """
        return self.submodule(self.free_module().zero_submodule(), check=False)


def _dict_set(v, n, key, val):
    r"""
    Rough-and-ready implementation of a two-layer-deep dictionary.

    EXAMPLES::

        sage: from sage.modular.hecke.module import _dict_set
        sage: v = {}
        sage: _dict_set(v, 1, 2, 3)
        sage: v
        {1: {2: 3}}
        sage: _dict_set(v, 1, 3, 4); v
        {1: {2: 3, 3: 4}}
    """
    if n in v:
        v[n][key] = val
    else:
        v[n] = {key: val}
