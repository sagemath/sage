# sage.doctest: needs sage.libs.flint sage.libs.pari
"""
Ambient Hecke modules
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
from . import degenmap
from . import module
from . import submodule

from sage.modules.free_module import FreeModule, FreeModule_generic
from sage.rings.integer import Integer

import sage.arith.all as arith

import sage.matrix.matrix_space as matrix_space
from sage.matrix.constructor import matrix

from sage.modular.arithgroup.all import Gamma0  # for Sturm bound


def is_AmbientHeckeModule(x) -> bool:
    r"""
    Return ``True`` if ``x`` is of type ``AmbientHeckeModule``.

    EXAMPLES::

        sage: from sage.modular.hecke.ambient_module import is_AmbientHeckeModule
        sage: is_AmbientHeckeModule(ModularSymbols(6))
        doctest:warning...
        DeprecationWarning: the function is_AmbientHeckeModule is deprecated;
        use 'isinstance(..., AmbientHeckeModule)' instead
        See https://github.com/sagemath/sage/issues/37895 for details.
        True
        sage: is_AmbientHeckeModule(ModularSymbols(6).cuspidal_subspace())
        False
        sage: is_AmbientHeckeModule(ModularForms(11))
        True
        sage: is_AmbientHeckeModule(BrandtModule(2, 3))
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(37895, "the function is_AmbientHeckeModule is deprecated; use 'isinstance(..., AmbientHeckeModule)' instead")
    return isinstance(x, AmbientHeckeModule)


class AmbientHeckeModule(module.HeckeModule_free_module):
    """
    An ambient Hecke module, i.e. a Hecke module that is isomorphic as a module
    over its base ring `R` to the standard free module `R^k` for some `k`. This
    is the base class for ambient spaces of modular forms and modular symbols,
    and for Brandt modules.
    """
    def __init__(self, base_ring, rank, level, weight, category=None):
        r"""
        Create an ambient Hecke module.

        EXAMPLES::

            sage: ModularSymbols(6) # indirect doctest
            Modular Symbols space of dimension 3 for Gamma_0(6) of weight 2 with sign 0 over Rational Field
            sage: sage.modular.hecke.ambient_module.AmbientHeckeModule(QQ, 3, 2, 4)
            Generic ambient Hecke module of rank 3, level 2 and weight 4 over Rational Field
        """
        rank = Integer(rank)
        if rank < 0:
            raise ValueError("rank (=%s) must be nonnegative" % rank)
        self.__rank = rank
        module.HeckeModule_free_module.__init__(self, base_ring, level,
                                                weight, category=category)

    def rank(self):
        """
        Return the rank of this ambient Hecke module.

        OUTPUT: integer

        EXAMPLES::

            sage: M = sage.modular.hecke.ambient_module.AmbientHeckeModule(QQ, 3, 11, 2); M
            Generic ambient Hecke module of rank 3, level 11 and weight 2 over Rational Field
            sage: M.rank()
            3
        """
        return self.__rank

    def __add__(self, other):
        r"""
        Return the sum of ``self`` and ``other``.

        As ``self`` is an ambient space, this will only make sense if
        ``other`` is a subspace of ``self``, in which case the answer
        is ``self``.

        EXAMPLES::

            sage: M = ModularSymbols(23)
            sage: M + M is M
            True
            sage: M + 3
            Traceback (most recent call last):
            ...
            TypeError: other (=3) must be a Hecke module
        """
        if not isinstance(other, module.HeckeModule_free_module):
            raise TypeError("other (=%s) must be a Hecke module" % other)
        if other.ambient_hecke_module() == self:
            return self
        raise ArithmeticError("sum only defined for subspaces of a common ambient Hecke module")

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        This should be overridden by derived classes.

        EXAMPLES::

            sage: sage.modular.hecke.ambient_module.AmbientHeckeModule(QQ, 3, 2, 4)._repr_()
            'Generic ambient Hecke module of rank 3, level 2 and weight 4 over Rational Field'
        """
        return "Generic ambient Hecke module of rank %s, level %s and weight %s over %s" % (self.rank(), self.level(), self.weight(), self.base_ring())

    def _degeneracy_raising_matrix(self, codomain):
        """
        Matrix of the degeneracy map (with t = 1) from ``self`` to ``codomain``,
        whose level should be a multiple of the level of ``self``.

        EXAMPLES::

            sage: sage.modular.hecke.ambient_module.AmbientHeckeModule(QQ, 3, 2, 4)._degeneracy_raising_matrix(4)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _degeneracy_lowering_matrix(self, codomain, t):
        """
        Matrix of the degeneracy map of index t from ``self`` to ``codomain``,
        whose level should be a divisor of the level of ``self``.

        EXAMPLES::

            sage: sage.modular.hecke.ambient_module.AmbientHeckeModule(QQ, 3, 2, 4)._degeneracy_lowering_matrix(2, 2)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _hecke_image_of_ith_basis_element(self, n, i):
        """
        Return the image under the Hecke operator `T_n` of the `i`-th basis
        element.

        EXAMPLES::

            sage: sage.modular.hecke.ambient_module.AmbientHeckeModule(QQ, 3, 2, 4)._hecke_image_of_ith_basis_element(4, 2)
            Traceback (most recent call last):
            ...
            NotImplementedError: All subclasses must implement _compute_hecke_matrix_prime
        """
        return self.hecke_operator(n)(self.gen(i))

    def _set_dual_free_module(self, V):
        r"""
        Store the embedded dual module of this module.

        Since this module is an ambient module, this is not necessary.

        EXAMPLES::

            sage: ModularForms(11, 2)._set_dual_free_module(None)
        """
        pass  # setting dual free module of ambient space is not necessary

    def ambient_hecke_module(self):
        r"""
        Return the ambient space that contains this ambient space.

        This is, of course, just this space again.

        EXAMPLES::

            sage: M = ModularForms(11, 4); M.ambient_hecke_module() is M
            True
        """
        return self

    def complement(self):
        """
        Return the largest Hecke-stable complement of this space.

        EXAMPLES::

            sage: M = ModularSymbols(11,2,1); M
            Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
            sage: M.complement()
            Modular Symbols subspace of dimension 0 of Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
            sage: C = M.cuspidal_subspace(); C
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
            sage: C.complement()
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
        """
        return self.zero_submodule()

    def decomposition_matrix(self):
        r"""
        Return the matrix whose columns form a basis for the canonical
        sorted decomposition of ``self`` coming from the Hecke operators.

        If the simple factors are `D_0, \ldots, D_n`, then the
        first few columns are an echelonized basis for `D_0`, the
        next an echelonized basis for `D_1`, the next for
        `D_2`, etc.

        EXAMPLES::

            sage: S = ModularSymbols(37, 2)
            sage: S.decomposition_matrix()
            [   1    0    0    0 -1/3]
            [   0    1   -1    0  1/2]
            [   0    0    0    1 -1/2]
            [   0    1    1    1    0]
            [   0    0    0    0    1]
        """
        try:
            return self.__decomposition_matrix_cache
        except AttributeError:
            rows = []
            for A in self.decomposition():
                for x in A.basis():
                    rows.append(x.list())
            A = matrix_space.MatrixSpace(self.base_ring(), self.rank())(rows)
            self.__decomposition_matrix_cache = A
            return self.__decomposition_matrix_cache

    def decomposition_matrix_inverse(self):
        """
        Return the inverse of the decomposition matrix.

        This is the inverse of the matrix returned by
        :meth:`decomposition_matrix`.

        EXAMPLES::

            sage: S = ModularSymbols(37, 2)
            sage: t = S.decomposition_matrix_inverse(); t
            [   1    0    0    0  1/3]
            [   0  1/2 -1/2  1/2 -1/2]
            [   0 -1/2 -1/2  1/2    0]
            [   0    0    1    0  1/2]
            [   0    0    0    0    1]
            sage: t * S.decomposition_matrix() == 1
            True
        """
        try:
            return self.__decomposition_matrix_inverse_cache
        except AttributeError:
            self.__decomposition_matrix_inverse_cache = ~self.decomposition_matrix()
            return self.__decomposition_matrix_inverse_cache

    def degeneracy_map(self, codomain, t=1):
        """
        The `t`-th degeneracy map from ``self`` to the module ``codomain``.

        The level of the codomain must be a divisor or multiple of
        level, and `t` must be a divisor of the quotient.

        INPUT:

        - ``codomain`` -- a Hecke module, which should be of the same type as
          self, or a positive integer (in which case Sage will use
          :meth:`~hecke_module_of_level` to find the "natural" module of the
          corresponding level).
        - ``t`` -- integer; the parameter of the degeneracy map, i.e., the map
          is related to `f(q)` - `f(q^t)`

        OUTPUT: a morphism from ``self`` to ``codomain``

        EXAMPLES::

            sage: M = ModularSymbols(11,sign=1)
            sage: d1 = M.degeneracy_map(33); d1
            Hecke module morphism degeneracy map corresponding to f(q) |--> f(q)
            defined by the matrix
            [ 1  0  0  0 -2 -1]
            [ 0 -1  1  0  0  0]
            Domain:   Modular Symbols space of dimension 2 for Gamma_0(11) of weight ...
            Codomain: Modular Symbols space of dimension 6 for Gamma_0(33) of weight ...
            sage: M.degeneracy_map(33, 3).matrix()
            [ 3  2  0  2 -2  1]
            [ 0  0 -1  1  0  0]
            sage: M = ModularSymbols(33, sign=1)
            sage: d2 = M.degeneracy_map(11); d2.matrix()
            [ 1  0]
            [ 0 -2]
            [ 0  2]
            [ 0  1]
            [-1  0]
            [-1  0]
            sage: (d2*d1).matrix()
            [4 0]
            [0 4]

        ::

            sage: M = ModularSymbols(3, 12, sign=1)
            sage: M.degeneracy_map(1)
            Hecke module morphism degeneracy map corresponding to f(q) |--> f(q)
            defined by the matrix
            [1 0]
            [0 0]
            [0 1]
            [0 1]
            [0 1]
            Domain: Modular Symbols space of dimension 5 for Gamma_0(3) of weight ...
            Codomain: Modular Symbols space of dimension 2 for Gamma_0(1) of weight ...

        ::

            sage: S = M.cuspidal_submodule()
            sage: S.degeneracy_map(1)
            Hecke module morphism defined by the matrix
            [1 0]
            [0 0]
            [0 0]
            Domain: Modular Symbols subspace of dimension 3 of Modular Symbols space ...
            Codomain: Modular Symbols space of dimension 2 for Gamma_0(1) of weight ...

        ::

            sage: D = ModularSymbols(10,4).cuspidal_submodule().decomposition()
            sage: D
            [
             Modular Symbols subspace of dimension 2 of
              Modular Symbols space of dimension 10 for
               Gamma_0(10) of weight 4 with sign 0 over Rational Field,
             Modular Symbols subspace of dimension 4 of
              Modular Symbols space of dimension 10 for
               Gamma_0(10) of weight 4 with sign 0 over Rational Field
            ]
            sage: D[1].degeneracy_map(5)
            Hecke module morphism defined by the matrix
            [   0    0   -1    1]
            [   0  1/2  3/2   -2]
            [   0   -1    1    0]
            [   0 -3/4 -1/4    1]
            Domain:   Modular Symbols subspace of dimension 4 of Modular Symbols space ...
            Codomain: Modular Symbols space of dimension 4 for Gamma_0(5) of weight ...

        We check for a subtle caching bug that came up in work on :issue:`10453`::

            sage: loads(dumps(J0(33).decomposition()[0].modular_symbols()))
            Modular Symbols subspace of dimension 2 of
             Modular Symbols space of dimension 9 for
              Gamma_0(33) of weight 2 with sign 0 over Rational Field

        We check that certain absurd inputs are correctly caught::

            sage: chi = kronecker_character(7)
            sage: ModularSymbols(Gamma0(7), 4).degeneracy_map(ModularSymbols(chi, 4))
            Traceback (most recent call last):
            ...
            ValueError: the characters of the domain and codomain must match
        """
        if isinstance(codomain, AmbientHeckeModule):
            M = codomain
            level = int(M.level())
        else:
            level = int(codomain)
            M = None

        t = int(t)

        err = False
        if self.level() % level == 0:
            quo = self.level() // level
            if quo % t != 0:
                err = True
        elif level % self.level() == 0:
            quo = level // self.level()
            if quo % t != 0:
                err = True
        else:
            err = True
        if err:
            raise ValueError(("the level of self (=%s) must be a divisor or multiple of "
                              "level (=%s) and t (=%s) must be a divisor of the quotient") % (self.level(), level, t))

        eps = self.character()
        if eps is not None and level % eps.conductor() != 0:
            raise ArithmeticError("the conductor of the character of this space "
                                  "(=%s) must be divisible by the level (=%s)" % (eps.conductor(), level))

        if M is None:
            M = self.hecke_module_of_level(level)

        if eps is not None and M.character() is not None:
            if eps.primitive_character() != M.character().primitive_character():
                raise ValueError("the characters of the domain and codomain must match")

        key = (M.group(), t)
        # bad idea to use (M, t) as the key, because using complicated objects
        # like modular forms spaces as dictionary keys causes weird behaviour;
        # on the other hand, (M.level(), t) isn't enough information.
        try:
            self._degeneracy_maps
        except AttributeError:
            self._degeneracy_maps = {}

        if key in self._degeneracy_maps:
            return self._degeneracy_maps[key]

        if M.rank() == 0:
            A = matrix_space.MatrixSpace(self.base_ring(), self.rank(), 0)(0)

        elif self.level() % level == 0:  # lower the level
            A = self._degeneracy_lowering_matrix(M, t)

        elif level % self.level() == 0:  # raise the level
            A = self._degeneracy_raising_matrix(M, t)

        d = degenmap.DegeneracyMap(A, self, M, t)
        self._degeneracy_maps[key] = d
        return d

    def dual_free_module(self):
        r"""
        The free module dual to ``self``, as a submodule of the dual
        module of the ambient space. As this space is ambient anyway,
        this just returns ``self.free_module()``.

        EXAMPLES::

            sage: M = ModularForms(2,8); M.dual_free_module()
            Vector space of dimension 3 over Rational Field
            sage: M.dual_free_module() is M.free_module()
            True
        """
        return self.free_module()

    def fcp(self, n, var='x'):
        """
        Return the factorization of the characteristic polynomial of
        the Hecke operator `T_n` of index `n` acting on this space.

        INPUT:

        - ``self`` -- Hecke module invariant under the Hecke operator of index
           `n`

        - ``n`` --- a positive integer.

        - ``var`` --- variable of polynomial (default: ``'x'``)

        OUTPUT:

        - ``list`` -- list of the pairs `(g,e)`, where `g` is an
          irreducible factor of the characteristic polynomial of `T_n`, and `e`
          is its multiplicity.

        EXAMPLES::

            sage: m = ModularSymbols(23, 2, sign=1)
            sage: m.fcp(2)
            (x - 3) * (x^2 + x - 1)
            sage: m.hecke_operator(2).charpoly('x').factor()
            (x - 3) * (x^2 + x - 1)
        """
        n = int(n)
        if n <= 0:
            raise ArithmeticError("n (=%s) must be positive" % n)
        return self.hecke_operator(n).fcp(var)

    def free_module(self):
        """
        Return the free module underlying this ambient Hecke module (the
        forgetful functor from Hecke modules to modules over the base ring)

        EXAMPLES::

            sage: ModularForms(59, 2).free_module()
            Vector space of dimension 6 over Rational Field
        """
        try:
            return self.__free_module
        except AttributeError:
            M = FreeModule(self.base_ring(), self.rank())
            self.__free_module = M
            return M

    def hecke_bound(self):
        r"""
        Return an integer `B` such that the Hecke operators `T_n`, for `n\leq B`,
        generate the full Hecke algebra as a module over the base ring. Note
        that we include the `n` with `n` not coprime to the level.

        At present this returns an unproven guess for non-cuspidal spaces which
        appears to be valid for `M_k(\Gamma_0(N))`, where `k` and `N` are the
        weight and level of ``self``. (It is clearly valid for *cuspidal* spaces
        of any fixed character, as a consequence of the Sturm bound theorem.)
        It returns a hopelessly wrong answer for spaces of full level
        `\Gamma_1`.

        TODO: Get rid of this dreadful bit of code.

        EXAMPLES::

            sage: ModularSymbols(17, 4).hecke_bound()
            15
            sage: ModularSymbols(Gamma1(17), 4).hecke_bound()  # wrong!
            15
        """
        from sage.misc.verbose import verbose
        try:
            if self.is_cuspidal():
                return Gamma0(self.level()).sturm_bound(self.weight())
        except AttributeError:
            pass
        verbose("WARNING: ambient.py -- hecke_bound; returning unproven guess.")
        return Gamma0(self.level()).sturm_bound(self.weight()) + 2 * Gamma0(self.level()).dimension_eis(self.weight()) + 5

    def hecke_module_of_level(self, level):
        r"""
        Return the Hecke module corresponding to ``self`` at the given level, which
        should be either a divisor or a multiple of the level of ``self``.

        This raises :exc:`NotImplementedError`, and should be overridden in
        derived classes.

        EXAMPLES::

            sage: sage.modular.hecke.ambient_module.AmbientHeckeModule.hecke_module_of_level(ModularForms(2, 8),6)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def hecke_images(self, i, v):
        """
        Return images of the `i`-th standard basis vector under the
        Hecke operators `T_p` for all integers in `v`.

        INPUT:

        - ``i`` -- nonnegative integer

        - ``v`` -- list of positive integer

        OUTPUT: ``matrix`` -- whose rows are the Hecke images

        EXAMPLES::

            sage: M = ModularSymbols(DirichletGroup(13).0, 3)
            sage: M.T(2)(M.0).element()
            (zeta12 + 4, 0, -1, 1)
            sage: M.hecke_images(0, [1,2])
            [         1          0          0          0]
            [zeta12 + 4          0         -1          1]
        """
        try:
            return self._hecke_images(i, v)
        except (AttributeError, NotImplementedError):
            pass
        # Use slow generic algorithm
        x = self.gen(i)
        X = [self.hecke_operator(n).apply_sparse(x).element() for n in v]
        return matrix(self.base_ring(), X)

    def intersection(self, other):
        """
        Return the intersection of ``self`` and ``other``, which must both lie in
        a common ambient space of modular symbols.

        EXAMPLES::

            sage: M = ModularSymbols(43, sign=1)
            sage: A = M[0] + M[1]
            sage: B = M[1] + M[2]
            sage: A.rank(), B.rank()
            (2, 3)
            sage: C = A.intersection(B); C.rank()  # TODO
            1
        """
        if not isinstance(other, module.HeckeModule_free_module):
            raise TypeError("other (=%s) must be a Hecke module" % other)
        if self.ambient_hecke_module() != other.ambient_hecke_module():
            raise ArithmeticError("intersection only defined for subspaces of a common ambient Hecke module")
        return other  # since self is ambient, so the intersection must equal other.

    def is_ambient(self) -> bool:
        r"""
        Return ``True`` if and only if ``self`` is an ambient Hecke module.

        .. warning::

           ``self`` can only be ambient by being of type
           :class:`AmbientHeckeModule`.

           For example, decomposing a simple ambient space yields a
           single factor, and that factor is *not* considered an
           ambient space.

        EXAMPLES::

            sage: m = ModularSymbols(10)
            sage: m.is_ambient()
            True

        ::

            sage: a = m[0]  # the unique simple factor
            sage: a == m
            True
            sage: a.is_ambient()
            False
        """
        return True

    def is_full_hecke_module(self, compute=True) -> bool:
        """
        Return ``True`` if this space is invariant under the action of
        all Hecke operators, even those that divide the level. This is
        always true for ambient Hecke modules, so return ``True``.

        EXAMPLES::

            sage: ModularSymbols(11, 4).is_full_hecke_module()
            True
        """
        return True

    def is_new(self, p=None) -> bool:
        r"""
        Return ``True`` if this module is entirely new.

        EXAMPLES::

            sage: ModularSymbols(11, 4).is_new()
            False
            sage: ModularSymbols(1, 12).is_new()
            True
        """
        try:
            if p in self.__is_new:
                return self.__is_new[p]
        except AttributeError:
            pass
        AmbientHeckeModule.new_submodule(self, p)
        return self.__is_new[p]

    def is_old(self, p=None):
        r"""
        Return ``True`` if this module is entirely old.

        EXAMPLES::

            sage: ModularSymbols(22).is_old()
            True
            sage: ModularSymbols(3, 12).is_old()
            False
        """
        try:
            if p in self.__is_old:
                return self.__is_old[p]
        except AttributeError:
            pass
        self.old_submodule(p)
        return self.__is_old[p]

    def is_submodule(self, V):
        """
        Return ``True`` if and only if ``self`` is a submodule of ``V``.

        Since this is an ambient space, this returns ``True`` if and
        only if ``V`` is equal to ``self``.

        EXAMPLES::

            sage: ModularSymbols(1, 4).is_submodule(ModularSymbols(11,4))
            False
            sage: ModularSymbols(11, 4).is_submodule(ModularSymbols(11,4))
            True
        """
        if not isinstance(V, module.HeckeModule_free_module):
            raise TypeError("V must be a Hecke module")
        if not V.is_ambient():
            return False
        return V.ambient_hecke_module() == self

    def linear_combination_of_basis(self, v):
        r"""
        Given a list or vector of length equal to the dimension of ``self``,
        construct the appropriate linear combination of the basis vectors of
        ``self``.

        EXAMPLES::

            sage: ModularForms(3, 12).linear_combination_of_basis([1,0,0,0,1])
            2*q + 2049*q^2 + 177147*q^3 + 4196177*q^4 + 48830556*q^5 + O(q^6)
        """
        return self(v)

    def new_submodule(self, p=None):
        """
        Return the new or `p`-new submodule of ``self``.

        INPUT:

        - ``p`` -- (default: ``None``) if not ``None``, return only
          the `p`-new submodule

        OUTPUT: the new or `p`-new submodule of ``self``, i.e. the intersection of
        the kernel of the degeneracy lowering maps to level `N/p` (for the
        given prime `p`, or for all prime divisors of `N` if `p` is not given).

        If ``self`` is cuspidal, this is a Hecke-invariant complement of the
        corresponding old submodule, but this may break down on Eisenstein
        subspaces (see the amusing example in William Stein's book of a form
        which is new and old at the same time).

        EXAMPLES::

            sage: m = ModularSymbols(33); m.rank()
            9
            sage: m.new_submodule().rank()
            3
            sage: m.new_submodule(3).rank()
            4
            sage: m.new_submodule(11).rank()
            8
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

        # Construct the degeneracy map d.
        N = self.level()
        d = None
        eps = self.character()
        if eps is None:
            f = 1
        else:
            f = eps.conductor()
        if p is None:
            D = arith.prime_divisors(N)
        else:
            if N % p != 0:
                raise ValueError("p must divide the level")
            D = [p]
        for q in D:
            # Here we are only using degeneracy *lowering* maps, so it is fine
            # to be careless and pass an integer for the level. One needs to be
            # a bit more careful with degeneracy *raising* maps for the Gamma1
            # and GammaH cases.
            if ((N // q) % f) == 0:
                NN = N // q
                d1 = self.degeneracy_map(NN, 1).matrix()
                if d is None:
                    d = d1
                else:
                    d = d.augment(d1)
                d = d.augment(self.degeneracy_map(NN, q).matrix())

        if d is None or d == 0:
            self.__is_new[p] = True
            return self
        else:
            self.__is_new[p] = False
        ns = self.submodule(d.kernel(), check=False)
        ns.__is_new = {p: True}
        ns._is_full_hecke_module = True
        self.__new_submodule[p] = ns
        return ns

    def nonembedded_free_module(self):
        r"""
        Return the free module corresponding to ``self`` as an abstract free module
        (rather than as a submodule of an ambient free module).

        As this module is ambient anyway, this just returns
        ``self.free_module()``.

        EXAMPLES::

            sage: M = ModularSymbols(11, 2)
            sage: M.nonembedded_free_module() is M.free_module()
            True
        """
        return self.free_module()

    def old_submodule(self, p=None):
        """
        Return the old or `p`-old submodule of ``self``, i.e. the sum of the images
        of the degeneracy maps from level `N/p` (for the given prime `p`, or
        for all primes `p` dividing `N` if `p` is not given).

        INPUT:

        - ``p`` -- (default: ``None``) if not ``None``, return only the `p`-old
          submodule

        OUTPUT: the old or `p`-old submodule of ``self``

        EXAMPLES::

            sage: m = ModularSymbols(33); m.rank()
            9
            sage: m.old_submodule().rank()
            7
            sage: m.old_submodule(3).rank()
            6
            sage: m.new_submodule(11).rank()
            8

        ::

            sage: e = DirichletGroup(16)([-1, 1])
            sage: M = ModularSymbols(e, 3, sign=1); M
            Modular Symbols space of dimension 4 and level 16, weight 3,
             character [-1, 1], sign 1, over Rational Field
            sage: M.old_submodule()
            Modular Symbols subspace of dimension 3 of
             Modular Symbols space of dimension 4 and level 16, weight 3,
              character [-1, 1], sign 1, over Rational Field

        Illustrate that :issue:`10664` is fixed::

            sage: ModularSymbols(DirichletGroup(42)[7], 6, sign=1).old_subspace(3)
            Modular Symbols subspace of dimension 0 of
             Modular Symbols space of dimension 40 and level 42, weight 6,
              character [-1, -1], sign 1, over Rational Field
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

        # Construct the degeneracy map d.
        N = self.level()
        d = None

        eps = self.character()
        if eps is None:
            f = 1
        else:
            f = eps.conductor()

        if p is None:
            D = arith.prime_divisors(N)
        else:
            if N % p != 0:
                raise ValueError("p must divide the level")
            D = [p]

        for q in D:
            NN = N // q
            if NN % f == 0:
                M = self.hecke_module_of_level(NN)

                # Here it is vital to pass self as an argument to
                # degeneracy_map, because M and the level N don't uniquely
                # determine self (e.g. the degeneracy map from level 1 to level
                # N could go to Gamma0(N), Gamma1(N) or anything in between)
                d1 = M.degeneracy_map(self, 1).matrix()

                if d is None:
                    d = d1
                else:
                    d = d.stack(d1)
                d = d.stack(M.degeneracy_map(self, q).matrix())

        if d is None:
            os = self.zero_submodule()
        else:
            os = self.submodule(d.image(), check=False)

        self.__is_old[p] = (os == self)

        os.__is_old = {p: True}
        os._is_full_hecke_module = True
        self.__old_submodule[p] = os
        return os

    def submodule(self, M, Mdual=None, check=True):
        """
        Return the Hecke submodule of ``self`` generated by `M`, which may be a
        submodule of the free module of ``self``, or a list of elements of ``self``.

        EXAMPLES::

            sage: M = ModularForms(37, 2)
            sage: A = M.submodule([M.newforms()[0].element(), M.newforms()[1].element()]); A
            Modular Forms subspace of dimension 2 of
             Modular Forms space of dimension 3 for
              Congruence Subgroup Gamma0(37) of weight 2 over Rational Field
        """
        if check:
            if not isinstance(M, FreeModule_generic):
                V = self.free_module()
                if isinstance(M, (list, tuple)):
                    M = V.span([V(x.element()) for x in M])
                else:
                    M = V.span(M)
            if not M.is_submodule(self.free_module()):
                raise TypeError("M must be a submodule of the free module associated to this module")
            if M == self.free_module():
                return self
        return self._submodule_class()(self, M, Mdual, check=check)

    def _submodule_class(self):
        r"""
        The class of submodules of this module. This is a separate method so it
        can be overridden in derived classes.

        EXAMPLES::

            sage: sage.modular.hecke.ambient_module.AmbientHeckeModule._submodule_class(ModularForms(1, 24))
            <class 'sage.modular.hecke.submodule.HeckeSubmodule'>
            sage: ModularForms(1, 24)._submodule_class()
            <class 'sage.modular.modform.submodule.ModularFormsSubmodule'>
        """
        return submodule.HeckeSubmodule

    def submodule_from_nonembedded_module(self, V, Vdual=None, check=True):
        """
        Create a submodule of this module, from a submodule of an ambient free
        module of the same rank as the rank of ``self``.

        INPUT:

        - ``V`` -- submodule of ambient free module of the same rank as the
          rank of ``self``

        - ``Vdual`` -- used to pass in dual submodule (may be ``None``)

        - ``check`` -- whether to check that submodule is Hecke-equivariant

        OUTPUT: Hecke submodule of ``self``

        EXAMPLES::

            sage: V = QQ^8
            sage: ModularForms(24, 2).submodule_from_nonembedded_module(V.submodule([0]))
            Modular Forms subspace of dimension 0 of
             Modular Forms space of dimension 8 for
              Congruence Subgroup Gamma0(24) of weight 2 over Rational Field
        """
        return self.submodule(V, Vdual, check=check)

    def submodule_generated_by_images(self, M):
        """
        Return the submodule of this ambient modular symbols space
        generated by the images under all degeneracy maps of `M`.

        The space `M` must have the same weight, sign, and group or
        character as this ambient space.

        EXAMPLES::

            sage: ModularSymbols(6, 12).submodule_generated_by_images(ModularSymbols(1,12))
            Modular Symbols subspace of dimension 12 of
             Modular Symbols space of dimension 22 for
              Gamma_0(6) of weight 12 with sign 0 over Rational Field
        """
        S = self.zero_submodule()
        if self.level() % M.level() == 0:
            D = arith.divisors(self.level() // M.level())
        elif M.level() % self.level() == 0:
            D = arith.divisors(M.level() // self.level())
        else:
            D = []
        for t in D:
            d = M.degeneracy_map(self, t)
            if d.codomain() != self:
                raise ArithmeticError("incompatible spaces of modular symbols")
            S += d.image()

        if self.is_full_hecke_module(compute=False):
            S._is_full_hecke_module = True

        return S
