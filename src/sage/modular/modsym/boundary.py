# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.libs.flint sage.libs.pari
r"""
Space of boundary modular symbols

Used mainly for computing the cuspidal subspace of modular symbols. The space
of boundary symbols of sign 0 is isomorphic as a Hecke module to the dual of
the space of Eisenstein series, but this does not give a useful method of
computing Eisenstein series, since there is no easy way to extract the constant
terms.

We represent boundary modular symbols as a sum of Manin symbols of the form
`[P, u/v]`, where `u/v` is a cusp for our group `G`. The group of boundary
modular symbols naturally embeds into a vector space `B_k(G)` (see Stein,
section 8.4, or Merel, section 1.4, where this space is called `\CC[\Gamma
\backslash \QQ]_k`, for a definition), which is a finite dimensional `\QQ`
vector space, of dimension equal to the number of cusps for `G` (if `k` is
even), or the number of regular cusps (if `k` is odd). The embedding
takes `[P, u/v]` to `P(u,v)\cdot [(u,v)]`. We represent the basis vectors by
pairs `[(u,v)]` with u, v coprime. On `B_k(G)`, we have the relations

.. MATH::

     [\gamma \cdot (u,v)] = [(u,v)]

for all `\gamma \in G` and

.. MATH::

     [(\lambda u, \lambda v)] = \operatorname{sign}(\lambda)^k [(u,v)]


for all `\lambda \in \QQ^\times`.

It's possible for these relations to kill a class, i.e., for a pair `[(u,v)]`
to be 0. For example, when `N=4` and `k=3` then `(-1,-2)` is equivalent mod
`\Gamma_1(4)` to `(1,2)` since `2=-2 \bmod 4` and `1=-1 \bmod 2`. But since `k`
is odd, `[(-1,-2)]` is also equivalent to `-[(1,2)]`. Thus this symbol is
equivalent to its negative, hence 0 (notice that this wouldn't be the case in
characteristic 2). This happens for any irregular cusp when the weight is odd;
there are no irregular cusps on `\Gamma_1(N)` except when `N = 4`, but there
can be more on `\Gamma_H` groups. See also prop 2.30 of Stein's Ph.D. thesis.

In addition, in the case that our space is of sign `\sigma = 1` or `-1`, we
also have the relation `[(-u,v)] = \sigma \cdot [(u,v)]`. This relation can
also combine with the above to kill a cusp class - for instance, take (u,v) =
(1,3) for `\Gamma_1(5)`. Then since the cusp `\tfrac{1}{3}` is
`\Gamma_1(5)`-equivalent to the cusp `-\tfrac{1}{3}`, we have that `[(1,3)] =
[(-1,3)]`. Now, on the minus subspace, we also have that `[(-1,3)] = -[(1,3)]`,
which means this class must vanish. Notice that this cannot be used to show
that `[(1,0)]` or `[(0,1)]` is 0.

.. NOTE::

   Special care must be taken when working with the images of the cusps 0 and
   `\infty` in `B_k(G)`. For all cusps *except* 0 and `\infty`, multiplying the
   cusp by -1 corresponds to taking `[(u,v)]` to `[(-u,v)]` in `B_k(G)`.  This
   means that `[(u,v)]` is equivalent to `[(-u,v)]` whenever `\tfrac{u}{v}` is
   equivalent to `-\tfrac{u}{v}`, except in the case of 0 and `\infty`.  We
   have the following conditions for `[(1,0)]` and `[(0,1)]`:

   - `[(0,1)] = \sigma \cdot [(0,1)]`, so `[(0,1)]` is 0 exactly when `\sigma =
     -1`.

   - `[(1,0)] = \sigma \cdot [(-1,0)]` and `[(1,0)] = (-1)^k [(-1,0)]`, so
     `[(1,0)] = 0` whenever `\sigma \ne (-1)^k`.

.. NOTE::

   For all the spaces of boundary symbols below, no work is done to determine
   the cusps for G at creation time. Instead, cusps are added as they are
   discovered in the course of computation. As a result, the rank of a space
   can change as a computation proceeds.

REFERENCES:

- Merel, "Universal Fourier expansions of modular
  forms." Springer LNM 1585 (1994), pg. 59-95.

- Stein, "Modular Forms, a computational approach." AMS (2007).
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.misc.repr import repr_lincomb
from sage.structure.richcmp import richcmp_method, richcmp

import sage.modules.free_module as free_module
from sage.modules.free_module_element import FreeModuleElement

import sage.modular.arithgroup.all as arithgroup
import sage.modular.cusps as cusps
import sage.modular.dirichlet as dirichlet
import sage.modular.hecke.all as hecke
from sage.modular.modsym.manin_symbol import ManinSymbol

from sage.rings.rational_field import Q as QQ
from sage.categories.rings import Rings

from . import element


class BoundarySpaceElement(hecke.HeckeModuleElement):
    def __init__(self, parent, x):
        """
        Create a boundary symbol.

        INPUT:

        - ``parent`` -- BoundarySpace; a space of boundary
          modular symbols

        - ``x`` -- dictionary with integer keys and values in the
          base field of parent

        EXAMPLES::

            sage: B = ModularSymbols(Gamma0(32), sign=-1).boundary_space()
            sage: B(Cusp(1,8))
            [1/8]
            sage: B.0
            [1/8]
            sage: type(B.0)
            <class 'sage.modular.modsym.boundary.BoundarySpaceElement'>
        """
        self.__x = x
        self.__vec = parent.free_module()(x)
        hecke.HeckeModuleElement.__init__(self, parent, self.__vec)

    def coordinate_vector(self):
        r"""
        Return ``self`` as a vector on the `\QQ`-vector space with basis
        ``self.parent()._known_cusps()``.

        EXAMPLES::

            sage: B = ModularSymbols(18,4,sign=1).boundary_space()
            sage: x = B(Cusp(1/2)) ; x
            [1/2]
            sage: x.coordinate_vector()
            (1)
            sage: ((18/5)*x).coordinate_vector()
            (18/5)
            sage: B(Cusp(0))
            [0]
            sage: x.coordinate_vector()
            (1)
            sage: x = B(Cusp(1/2)) ; x
            [1/2]
            sage: x.coordinate_vector()
            (1, 0)
        """
        return self.__vec

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: ModularSymbols(Gamma0(11), 2).boundary_space()(Cusp(0))._repr_()
            '[0]'
            sage: (-6*ModularSymbols(Gamma0(11), 2).boundary_space()(Cusp(0)))._repr_()
            '-6*[0]'
        """
        return repr_lincomb([('[' + repr(self.parent()._known_gens[i]) + ']', c)
                             for i, c in sorted(self.__x.items())])

    # can't inherit arithmetic operations from HeckeModule, because basis
    # dimension might change!

    def _add_(self, other):
        """
        Return ``self + other``. Assumes that other is a BoundarySpaceElement.

        EXAMPLES::

            sage: B = ModularSymbols(Gamma1(16), 4).boundary_space()
            sage: x = B(Cusp(2/7)) ; y = B(Cusp(13/16))
            sage: x + y # indirect doctest
            [2/7] + [13/16]
            sage: x + x # indirect doctest
            2*[2/7]
        """
        z = dict(other.__x)
        for i, c in self.__x.items():
            if i in z:
                z[i] += c
            else:
                z[i] = c
        return BoundarySpaceElement(self.parent(), z)

    def _sub_(self, other):
        """
        Return ``self - other``. Assumes that other is a BoundarySpaceElement.

        EXAMPLES::

            sage: B = ModularSymbols(Gamma1(16), 4).boundary_space()
            sage: x = B(Cusp(2/7)) ; y = B(Cusp(13/16))
            sage: x - y # indirect doctest
            [2/7] - [13/16]
            sage: x - x # indirect doctest
            0
        """
        z = dict(self.__x)
        for i, c in other.__x.items():
            if i in z:
                z[i] -= c
            else:
                z[i] = -c
        return BoundarySpaceElement(self.parent(), z)

    def _rmul_(self, other):
        r"""
        Return ``self * other``.

        Assumes that other can be coerced into
        ``self.parent().base_ring()``.

        EXAMPLES::

            sage: B = ModularSymbols(Gamma1(16), 4).boundary_space()
            sage: x = B(Cusp(2/7))
            sage: x*5 # indirect doctest
            5*[2/7]
            sage: x*-3/5 # indirect doctest
            -3/5*[2/7]
        """
        x = {}
        for i, c in self.__x.items():
            x[i] = c * other
        return BoundarySpaceElement(self.parent(), x)

    def _lmul_(self, other):
        r"""
        Return ``other * self``.

        Assumes that other can be coerced into
        ``self.parent().base_ring()``.

        EXAMPLES::

            sage: B = ModularSymbols(Gamma1(16), 4).boundary_space()
            sage: x = B(Cusp(13/16))
            sage: 11*x # indirect doctest
            11*[13/16]
            sage: 1/3*x # indirect doctest
            1/3*[13/16]
        """
        x = {}
        for i, c in self.__x.items():
            x[i] = other * c
        return BoundarySpaceElement(self.parent(), x)

    def __neg__(self):
        """
        Return -self.

        EXAMPLES::

            sage: B = ModularSymbols(Gamma1(16), 4).boundary_space()
            sage: x = B(Cusp(2/7))
            sage: -x # indirect doctest
            -[2/7]
            sage: -x + x # indirect doctest
            0
        """
        return self * (-1)


@richcmp_method
class BoundarySpace(hecke.HeckeModule_generic):
    def __init__(self,
                 group=arithgroup.Gamma0(1),
                 weight=2,
                 sign=0,
                 base_ring=QQ,
                 character=None):
        """
        Space of boundary symbols for a congruence subgroup of SL_2(Z).

        This class is an abstract base class, so only derived classes
        should be instantiated.

        INPUT:

        - ``weight`` -- integer; the weight

        - ``group`` -- arithgroup.congroup_generic.CongruenceSubgroup, a
          congruence subgroup

        - ``sign`` -- integer; either -1, 0, or 1

        - ``base_ring`` -- rings.Ring (defaults to the rational numbers)

        EXAMPLES::

            sage: B = ModularSymbols(Gamma0(11),2).boundary_space()
            sage: isinstance(B, sage.modular.modsym.boundary.BoundarySpace)
            True
            sage: B == loads(dumps(B))
            True
        """
        weight = int(weight)
        if weight <= 1:
            raise ArithmeticError("weight must be at least 2")
        if not isinstance(group, arithgroup.CongruenceSubgroupBase):
            raise TypeError("group must be a congruence subgroup")
        sign = int(sign)
        if base_ring not in Rings().Commutative():
            raise TypeError("base_ring must be a commutative ring")
        if character is None and isinstance(group, arithgroup.Gamma0_class):
            character = dirichlet.TrivialCharacter(group.level(), base_ring)
        (self.__group, self.__weight, self.__character,
         self.__sign, self.__base_ring) = (group, weight,
                                           character, sign, base_ring)
        self._known_gens = []
        self._zero_cusps = []
        hecke.HeckeModule_generic.__init__(self, base_ring, group.level())

    def __richcmp__(self, other, op):
        """
        EXAMPLES::

            sage: B2 = ModularSymbols(11, 2).boundary_space()
            sage: B4 = ModularSymbols(11, 4).boundary_space()
            sage: B2 == B4
            False
            sage: B2 == ModularSymbols(17, 2).boundary_space()
            False
        """
        if type(self) is not type(other):
            return NotImplemented

        return richcmp((self.group(), self.weight(), self.character()),
                       (other.group(), other.weight(), other.character()),
                       op)

    def _known_cusps(self):
        """
        Return the list of non-vanishing cusps found so far.

        EXAMPLES::

            sage: B = ModularSymbols(Gamma1(12), 4).boundary_space()
            sage: B._known_cusps()
            []
            sage: ls = [ B(Cusp(i,10)) for i in range(10) ]
            sage: B._known_cusps()
            [0, 1/10, 1/5]
        """
        return list(self._known_gens)

    def is_ambient(self):
        """
        Return ``True`` if ``self`` is a space of boundary symbols
        associated to an ambient space of modular symbols.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(6), 4)
            sage: M.is_ambient()
            True
            sage: M.boundary_space().is_ambient()
            True
        """
        return True

    def group(self):
        """
        Return the congruence subgroup associated to this space of boundary
        modular symbols.

        EXAMPLES::

            sage: ModularSymbols(GammaH(14,[9]), 2).boundary_space().group()
            Congruence Subgroup Gamma_H(14) with H generated by [9]
        """
        return self.__group

    def weight(self):
        """
        Return the weight of this space of boundary modular symbols.

        EXAMPLES::

            sage: ModularSymbols(Gamma1(9), 5).boundary_space().weight()
            5
        """
        return self.__weight

    def character(self):
        """
        Return the Dirichlet character associated to this space of boundary
        modular symbols.

        EXAMPLES::

            sage: ModularSymbols(DirichletGroup(7).0, 6).boundary_space().character()
            Dirichlet character modulo 7 of conductor 7 mapping 3 |--> zeta6
        """
        return self.__character

    def sign(self):
        """
        Return the sign of the complex conjugation involution on this space
        of boundary modular symbols.

        EXAMPLES::

            sage: ModularSymbols(13,2,sign=-1).boundary_space().sign()
            -1
        """
        return self.__sign

    def gen(self, i=0):
        """
        Return the i-th generator of this space.

        EXAMPLES::

            sage: B = ModularSymbols(Gamma0(24), 4).boundary_space()
            sage: B.gen(0)
            Traceback (most recent call last):
            ...
            ValueError: only 0 generators known for Space of Boundary Modular Symbols for Congruence Subgroup Gamma0(24) of weight 4 over Rational Field
            sage: B(Cusp(1/3))
            [1/3]
            sage: B.gen(0)
            [1/3]
        """
        if i >= len(self._known_gens) or i < 0:
            raise ValueError("only %s generators known for %s" % (len(self._known_gens), self))
        return BoundarySpaceElement(self, {i: 1})

    def free_module(self):
        """
        Return the underlying free module for ``self``.

        EXAMPLES::

            sage: B = ModularSymbols(Gamma1(7), 5, sign=-1).boundary_space()
            sage: B.free_module()
            Sparse vector space of dimension 0 over Rational Field
            sage: x = B(Cusp(0)) ; y = B(Cusp(1/7)) ; B.free_module()
            Sparse vector space of dimension 1 over Rational Field
        """
        return free_module.FreeModule(self.__base_ring, len(self._known_gens), sparse=True)

    def rank(self):
        """
        The rank of the space generated by boundary symbols that have been
        found so far in the course of computing the boundary map.

        .. warning::

           This number may change as more elements are coerced into
           this space!! (This is an implementation detail that will
           likely change.)

        EXAMPLES::

            sage: M = ModularSymbols(Gamma0(72), 2) ; B = M.boundary_space()
            sage: B.rank()
            0
            sage: _ = [ B(x) for x in M.basis() ]
            sage: B.rank()
            16

        Test that :issue:`7837` is fixed::

            sage: ModularSymbols(Gamma1(4),7).boundary_map().codomain().dimension()
            2
            sage: ModularSymbols(Gamma1(4),7, sign=1).boundary_map().codomain().dimension()
            1
            sage: ModularSymbols(Gamma1(4),7, sign=-1).boundary_map().codomain().dimension()
            1
        """
        return len(self._known_gens)

    #####################################################################
    # Coercion
    #####################################################################

    def _coerce_in_manin_symbol(self, x):
        """
        Coerce the Manin symbol x into ``self``. (That is, return the image of
        x under the boundary map.)

        Assumes that x is associated to the same space of modular symbols
        as ``self``.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(5), 4) ; B = M.boundary_space()
            sage: [ B(x) for x in M.basis() ]
            [-[2/5], -[Infinity], -[1/3], -[-1/4], -[-1/4], -[-1/4]]
            sage: [ B._coerce_in_manin_symbol(x) for x in M.manin_symbols_basis() ]
            [-[2/5], -[Infinity], -[1/3], -[-1/4], -[-1/4], -[-1/4]]
        """
        i = x.i
        alpha, beta = x.endpoints(self.level())
        if self.weight() == 2:
            return self(alpha) - self(beta)
        if i == 0:
            return self(alpha)
        elif i == self.weight() - 2:
            return -self(beta)
        else:
            return self(0)

    def __call__(self, x):
        """
        Coerce x into a boundary symbol space.

        If x is a modular symbol (with the same group, weight, character,
        sign, and base field), this returns the image of that modular
        symbol under the boundary map.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma0(15), 2) ; B = M.boundary_space()
            sage: B(M.0)
            [Infinity] - [0]
            sage: B(Cusp(1))
            [0]
            sage: B(Cusp(oo))
            [Infinity]
            sage: B(7)
            Traceback (most recent call last):
            ...
            TypeError: Coercion of 7 (of type <class 'sage.rings.integer.Integer'>) into Space of Boundary Modular Symbols for Congruence Subgroup Gamma0(15) of weight 2 over Rational Field not (yet) defined.
        """
        from .ambient import ModularSymbolsAmbient
        if isinstance(x, int) and x == 0:
            return BoundarySpaceElement(self, {})

        elif isinstance(x, cusps.Cusp):
            return self._coerce_cusp(x)

        elif isinstance(x, ManinSymbol):
            return self._coerce_in_manin_symbol(x)

        elif isinstance(x, element.ModularSymbolsElement):
            M = x.parent()
            if not isinstance(M, ModularSymbolsAmbient):
                raise TypeError("x (=%s) must be an element of a space of modular symbols of type ModularSymbolsAmbient" % x)
            if M.level() != self.level():
                raise TypeError("x (=%s) must have level %s but has level %s" % (
                    x, self.level(), M.level()))
            S = x.manin_symbol_rep()
            if len(S) == 0:
                return self(0)
            return sum([c * self._coerce_in_manin_symbol(v) for c, v in S])

        elif isinstance(x, FreeModuleElement):
            y = dict(enumerate(x))
            return BoundarySpaceElement(self, y)

        raise TypeError("Coercion of %s (of type %s) into %s not (yet) defined." % (x, type(x), self))

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: sage.modular.modsym.boundary.BoundarySpace(Gamma0(3), 2)._repr_()
            'Space of Boundary Modular Symbols of weight 2 for Congruence Subgroup Gamma0(3) with sign 0 and character [1] over Rational Field'
        """
        return ("Space of Boundary Modular Symbols of weight %s for" +
                " %s with sign %s and character %s over %s") % (
                    self.weight(), self.group(), self.sign(),
                    self.character()._repr_short_(), self.base_ring())

    def _cusp_index(self, cusp):
        """
        Return the index of the (unique) cusp in self._known_cusps()
        equivalent to cusp, or one of the following special values:

        - ``-1`` if ``cusp`` is not equivalent to any cusp
          found so far.
        - ``-2`` if ``cusp`` is equivalent to a cusp that's known to vanish
          from the relations in this space.

        EXAMPLES::

            sage: B = ModularSymbols(Gamma0(21), 4).boundary_space()
            sage: B._cusp_index(Cusp(0))
            -1
            sage: _ = B(Cusp(oo))
            sage: _ = B(Cusp(0))
            sage: B._cusp_index(Cusp(0))
            1
        """
        g = self._known_gens
        for i in range(len(g)):
            if self._is_equiv(cusp, g[i]):
                return i
        if any(self._is_equiv(cusp, x) for x in self._zero_cusps):
            return -2
        return -1


class BoundarySpace_wtk_g0(BoundarySpace):
    def __init__(self, level, weight, sign, F):
        """
        Initialize a space of boundary symbols of weight k for Gamma_0(N)
        over base field F.

        INPUT:

        - ``level`` -- integer; the level

        - ``weight`` -- integer; weight = 2

        - ``sign`` -- integer; either -1, 0, or 1

        - ``F`` -- field

        EXAMPLES::

            sage: B = ModularSymbols(Gamma0(2), 5).boundary_space()
            sage: type(B)
            <class 'sage.modular.modsym.boundary.BoundarySpace_wtk_g0_with_category'>
            sage: B == loads(dumps(B))
            True
        """
        level = int(level)
        sign = int(sign)
        weight = int(weight)
        if sign not in [-1, 0, 1]:
            raise ArithmeticError("sign must be an int in [-1,0,1]")
        if level <= 0:
            raise ArithmeticError("level must be positive")
        BoundarySpace.__init__(self,
                               weight=weight,
                               group=arithgroup.Gamma0(level),
                               sign=sign,
                               base_ring=F)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: B = ModularSymbols(Gamma0(97), 3).boundary_space()
            sage: B._repr_()
            'Space of Boundary Modular Symbols for Congruence Subgroup Gamma0(97) of weight 3 over Rational Field'
        """
        return ("Space of Boundary Modular Symbols for %s of weight %s over %s"
                % (self.group(), self.weight(), self.base_ring()))

    def _coerce_cusp(self, c):
        """
        Coerce the cusp c into this boundary symbol space.

        EXAMPLES::

            sage: B = ModularSymbols(Gamma0(17), 6).boundary_space()
            sage: B._coerce_cusp(Cusp(0))
            [0]
            sage: B = ModularSymbols(Gamma0(17), 6, sign=-1).boundary_space()
            sage: B._coerce_cusp(Cusp(0))
            0
            sage: B = ModularSymbols(Gamma0(16), 4).boundary_space()
            sage: [ B(Cusp(i,4)) for i in range(4) ]
            [[0], [1/4], [1/2], [3/4]]
            sage: B = ModularSymbols(Gamma0(16), 4, sign=1).boundary_space()
            sage: [ B(Cusp(i,4)) for i in range(4) ]
            [[0], [1/4], [1/2], [1/4]]
            sage: B = ModularSymbols(Gamma0(16), 4, sign=-1).boundary_space()
            sage: [ B(Cusp(i,4)) for i in range(4) ]
            [0, [1/4], 0, -[1/4]]
        """
        if self.weight() % 2:
            return self(0)

        # see if we've already found this cusp
        i = self._cusp_index(c)
        if i != -1:
            if i == -2:
                return self(0)
            return BoundarySpaceElement(self, {i: 1})

        # see if we've already found -c
        sign = self.sign()
        if sign != 0:
            i2 = self._cusp_index(-c)
            if i2 != -1:
                if i2 == -2:
                    return self(0)
                else:
                    return BoundarySpaceElement(self, {i2: sign})

        # found a new cusp class
        g = self._known_gens
        g.append(c)

        # See if the new cusp is killed by sign relations. The
        # relevant relations (for cusps other than 0 and Infinity)
        # are:
        #
        #    [(u,v)] = (-1)^k [(-u,-v)]
        #    [(u,v)] = [gamma * (u,v)]
        #   [(-u,v)] = sign * [(u,v)]
        #
        # So since k is always even on Gamma0, we have that [(u,v)] =
        # 0 from the above relations exactly when (u,v) = gamma*(-u,v)
        # and the sign is -1.
        if sign == -1:
            # NOTE: this code looks wrong. One should do the
            # following:
            #
            #  - if c is 0, if the sign is -1, append & return 0
            #  - if c is Infinity, then if the sign
            #    is not equal to (-1)**self.weight(), then
            #    append & return 0
            #  - otherwise, if the sign is -1, and c is
            #    equivalent to -c, append & return 0.
            #
            # Interestingly, the code below does precisely that.
            # (It's important to recall that for Gamma0, odd weight
            # spaces are 0.)
            if self._is_equiv(c, -c):
                self._zero_cusps.append(c)
                del self._known_gens[-1]
                return self(0)

        return BoundarySpaceElement(self, {(len(g) - 1): 1})

    def _is_equiv(self, c1, c2):
        """
        Determine whether or not c1 and c2 are equivalent for ``self``.

        EXAMPLES::

            sage: B = ModularSymbols(Gamma0(24), 6).boundary_space()
            sage: B._is_equiv(Cusp(0), Cusp(oo))
            False
            sage: B._is_equiv(Cusp(0), Cusp(1))
            True
        """
        return c1.is_gamma0_equiv(c2, self.level())


class BoundarySpace_wtk_g1(BoundarySpace):
    def __init__(self, level, weight, sign, F):
        """
        Initialize a space of boundary modular symbols for Gamma1(N).

        INPUT:

        - ``level`` -- integer; the level

        - ``weight`` -- integer; the weight = 2

        - ``sign`` -- integer; either -1, 0, or 1

        - ``F`` -- base ring

        EXAMPLES::

            sage: from sage.modular.modsym.boundary import BoundarySpace_wtk_g1
            sage: B = BoundarySpace_wtk_g1(17, 2, 0, QQ) ; B
            Boundary Modular Symbols space for Gamma_1(17) of weight 2 over Rational Field
            sage: B == loads(dumps(B))
            True
        """
        level = int(level)
        sign = int(sign)
        if sign not in [-1, 0, 1]:
            raise ArithmeticError("sign must be an int in [-1,0,1]")
        if level <= 0:
            raise ArithmeticError("level must be positive")

        BoundarySpace.__init__(self,
                               weight=weight,
                               group=arithgroup.Gamma1(level),
                               sign=sign,
                               base_ring=F)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: ModularSymbols(Gamma1(5), 3, sign=1).boundary_space()._repr_()
            'Boundary Modular Symbols space for Gamma_1(5) of weight 3 over Rational Field'
        """
        return ("Boundary Modular Symbols space for Gamma_1(%s) of weight %s " +
                "over %s") % (self.level(), self.weight(), self.base_ring())

    def _is_equiv(self, c1, c2):
        """
        Return whether c1 and c2 are equivalent cusps for ``self``.

        EXAMPLES::

            sage: B = ModularSymbols(Gamma1(10), 4).boundary_space()
            sage: B._is_equiv(Cusp(0), Cusp(1/5))
            (False, 0)
            sage: B._is_equiv(Cusp(4/5), Cusp(1/5))
            (True, -1)
            sage: B._is_equiv(Cusp(-4/5), Cusp(1/5))
            (True, 1)
        """
        return c1.is_gamma1_equiv(c2, self.level())

    def _cusp_index(self, cusp):
        """
        Return a pair (i, t), where i is the index of the first cusp in
        self._known_cusps() which is equivalent to cusp, and t is 1 or -1
        as cusp is Gamma1-equivalent to plus or minus
        self._known_cusps()[i]. If cusp is not equivalent to any known
        cusp, return (-1, 0); if cusp is equivalent to a cusp that is known to
        vanish, return (-2, 0).

        EXAMPLES::

            sage: B = ModularSymbols(Gamma1(11),2).boundary_space()
            sage: B._cusp_index(Cusp(1/11))
            (-1, 0)
            sage: B._cusp_index(Cusp(10/11))
            (-1, 0)
            sage: B._coerce_cusp(Cusp(1/11))
            [1/11]
            sage: B._cusp_index(Cusp(1/11))
            (0, 1)
            sage: B._cusp_index(Cusp(10/11))
            (0, -1)
        """
        g = self._known_gens
        for i in range(len(g)):
            t, eps = self._is_equiv(cusp, g[i])
            if t:
                return i, eps
        if any(self._is_equiv(cusp, c)[0] for c in self._zero_cusps):
            return (-2, 0)
        return -1, 0

    def _coerce_cusp(self, c):
        """
        Coerce a cusp into this boundary symbol space.

        EXAMPLES::

            sage: B = ModularSymbols(Gamma1(4), 4).boundary_space()
            sage: B._coerce_cusp(Cusp(1/2))
            [1/2]
            sage: B._coerce_cusp(Cusp(1/4))
            [1/4]
            sage: B._coerce_cusp(Cusp(3/4))
            [1/4]
            sage: B = ModularSymbols(Gamma1(5), 3, sign=-1).boundary_space()
            sage: B._coerce_cusp(Cusp(0))
            0
            sage: B._coerce_cusp(Cusp(oo))
            [Infinity]
            sage: B = ModularSymbols(Gamma1(2), 3, sign=-1).boundary_space()
            sage: B._coerce_cusp(Cusp(0))
            0
            sage: B._coerce_cusp(Cusp(oo))
            0
            sage: B = ModularSymbols(Gamma1(7), 3).boundary_space()
            sage: [ B(Cusp(i,7)) for i in range(7) ]
            [[0], [1/7], [2/7], [3/7], -[3/7], -[2/7], -[1/7]]
            sage: B._is_equiv(Cusp(1,6), Cusp(5,6))
            (True, 1)
            sage: B._is_equiv(Cusp(1,6), Cusp(0))
            (True, -1)
            sage: B(Cusp(0))
            [0]
            sage: B = ModularSymbols(Gamma1(7), 3, sign=1).boundary_space()
            sage: [ B(Cusp(i,7)) for i in range(7) ]
            [[0], 0, 0, 0, 0, 0, 0]
            sage: B = ModularSymbols(Gamma1(7), 3, sign=-1).boundary_space()
            sage: [ B(Cusp(i,7)) for i in range(7) ]
            [0, [1/7], [2/7], [3/7], -[3/7], -[2/7], -[1/7]]
        """
        N = self.level()
        k = self.weight()
        sign = self.sign()
        i, eps = self._cusp_index(c)
        if i != -1:
            if i == -2:
                return self(0)
            else:
                return BoundarySpaceElement(self, {i: eps**k})

        if sign != 0:
            i2, eps = self._cusp_index(-c)
            if i2 != -1:
                if i2 == -2:
                    return self(0)
                else:
                    return BoundarySpaceElement(self, {i2: sign * (eps**k)})

        # found a new cusp class
        g = self._known_gens
        g.append(c)

        # Does cusp class vanish because of irregularity? For Gamma1 groups
        # this occurs if and only if N = 4 and c is equivalent to 1/2.

        if N == 4 and k % 2:
            if self._is_equiv(c, (1, 2))[0]:
                self._zero_cusps.append(c)
                del self._known_gens[-1]
                return self(0)

        # Does class vanish because of sign relations?  The relevant
        # relations are
        #
        #    [(u,v)] = (-1)^k [(-u,-v)]
        #    [(u,v)] = sign * [(-u,v)]
        #    [(u,v)] = eps * (-1)^k [(-u,v)]
        #
        # where, in the last line, (u,v) is Gamma1-equivalent to
        # (-u,v) or (u,-v) as eps is 1 or -1.
        #
        # Thus (other than for 0 and Infinity), we have that [(u,v)]
        # can only be killed by sign relations when:
        #
        #  - (u,v) is Gamma1-equivalent to (-u,v) or (u,-v), and
        #  - eps is 1 and sign is -1, or eps is -1 and sign is not
        #    (-1)^k.
        #
        if sign:
            if (c.is_infinity() and sign != (-1)**self.weight()) or \
               (c.is_zero() and sign == -1):
                self._zero_cusps.append(c)
                del self._known_gens[-1]
                return self(0)
            elif (not c.is_infinity() and not c.is_zero()):
                t, eps = self._is_equiv(c, -c)
                if t and ((eps == 1 and sign == -1) or
                          (eps == -1 and sign != (-1)**self.weight())):
                    self._zero_cusps.append(c)
                    del self._known_gens[-1]
                    return self(0)

        return BoundarySpaceElement(self, {(len(g) - 1): 1})


class BoundarySpace_wtk_gamma_h(BoundarySpace):
    def __init__(self, group, weight, sign, F):
        """
        Initialize a space of boundary modular symbols for GammaH(N).

        INPUT:

        - ``group`` -- congruence subgroup Gamma_H(N)

        - ``weight`` -- integer; the weight = 2

        - ``sign`` -- integer; either -1, 0, or 1

        - ``F`` -- base ring

        EXAMPLES::

            sage: from sage.modular.modsym.boundary import BoundarySpace_wtk_gamma_h
            sage: B = BoundarySpace_wtk_gamma_h(GammaH(13,[3]), 2, 0, QQ) ; B
            Boundary Modular Symbols space for Congruence Subgroup Gamma_H(13) with H generated by [3] of weight 2 over Rational Field
            sage: B == loads(dumps(B))
            True

        A test case from :issue:`6072`::

            sage: ModularSymbols(GammaH(8,[5]), 3).boundary_map()
            Hecke module morphism boundary map defined by the matrix
            [-1  0  0  0]
            [ 0 -1  0  0]
            [ 0  0 -1  0]
            [ 0  0  0 -1]
            Domain: Modular Symbols space of dimension 4 for Congruence Subgroup ...
            Codomain: Boundary Modular Symbols space for Congruence Subgroup Gamma_H(8) ...
        """
        sign = int(sign)
        if sign not in [-1, 0, 1]:
            raise ArithmeticError("sign must be an int in [-1,0,1]")

        BoundarySpace.__init__(self,
                               weight=weight,
                               group=group,
                               sign=sign,
                               base_ring=F)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: ModularSymbols(GammaH(7,[2]), 4).boundary_space()._repr_()
            'Boundary Modular Symbols space for Congruence Subgroup Gamma_H(7) with H generated by [2] of weight 4 over Rational Field'
        """
        return ("Boundary Modular Symbols space for %s of weight %s " +
                "over %s") % (self.group(), self.weight(), self.base_ring())

    def _is_equiv(self, c1, c2):
        """
        Return a pair of the form (b, t), where b is ``True`` if c1 and c2 are
        equivalent cusps for ``self``, and ``False`` otherwise, and t gives extra
        information about the equivalence between c1 and c2.

        EXAMPLES::

            sage: B = ModularSymbols(GammaH(7,[2]), 4).boundary_space()
            sage: B._is_equiv(Cusp(0), Cusp(1/7))
            (False, 0)
            sage: B._is_equiv(Cusp(2/7), Cusp(1/7))
            (True, 1)
            sage: B._is_equiv(Cusp(3/7), Cusp(1/7))
            (True, -1)
        """
        return c1.is_gamma_h_equiv(c2, self.group())

    def _cusp_index(self, cusp):
        """
        Return a pair (i, t), where i is the index of the first cusp in
        self._known_cusps() which is equivalent to cusp, and t is 1 or -1 as
        cusp is GammaH-equivalent to plus or minus self._known_cusps()[i]. If
        cusp is not equivalent to any known cusp, return (-1, 0); if cusp is
        equivalent to a cusp known to vanish, return (-2, 0).

        EXAMPLES::

            sage: M = ModularSymbols(GammaH(9,[4]), 3)
            sage: B = M.boundary_space()
            sage: B._cusp_index(Cusp(0))
            (-1, 0)
            sage: _ = [ B(x) for x in M.basis() ]
            sage: B._cusp_index(Cusp(0))
            (1, -1)
            sage: B._cusp_index(Cusp(5/6))
            (3, 1)
        """
        g = self._known_gens
        for i in range(len(g)):
            t, eps = self._is_equiv(cusp, g[i])
            if t:
                return i, eps
        if any(self._is_equiv(cusp, c)[0] for c in self._zero_cusps):
            return -2, 0
        return -1, 0

    def _coerce_cusp(self, c):
        """
        Coerce the cusp c into ``self``.

        EXAMPLES::

            sage: B = ModularSymbols(GammaH(10,[9]), 2).boundary_space()
            sage: B(Cusp(0))
            [0]
            sage: B(Cusp(1/3))
            [1/3]
            sage: B(Cusp(1/13))
            [1/3]
            sage: B = ModularSymbols(GammaH(25, [6]), 2).boundary_space()
            sage: B._coerce_cusp(Cusp(0))
            [0]

        ::

            sage: B = ModularSymbols(GammaH(11,[3]), 3).boundary_space()
            sage: [ B(Cusp(i,11)) for i in range(11) ]
            [[0],
            [1/11],
            -[1/11],
            [1/11],
            [1/11],
            [1/11],
            -[1/11],
            -[1/11],
            -[1/11],
            [1/11],
            -[1/11]]
            sage: B._is_equiv(Cusp(0), Cusp(1,11))
            (False, 0)
            sage: B._is_equiv(Cusp(oo), Cusp(1,11))
            (True, 1)
            sage: B = ModularSymbols(GammaH(11,[3]), 3, sign=1).boundary_space()
            sage: [ B(Cusp(i,11)) for i in range(11) ]
            [[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            sage: B = ModularSymbols(GammaH(11,[3]), 3, sign=-1).boundary_space()
            sage: [ B(Cusp(i,11)) for i in range(11) ]
            [0,
            [1/11],
            -[1/11],
            [1/11],
            [1/11],
            [1/11],
            -[1/11],
            -[1/11],
            -[1/11],
            [1/11],
            -[1/11]]

        Test that :issue:`6072` is fixed. ::

            sage: G = GammaH(8,[3])

            sage: B2 = G.modular_symbols(weight=2).boundary_space()
            sage: B2._coerce_cusp(Cusp(1/4))
            [1/4]

            sage: B3 = G.modular_symbols(weight=3).boundary_space()
            sage: B3._coerce_cusp(Cusp(1/4))
            0

            sage: G = GammaH(10, [9])
            sage: B3 = G.modular_symbols(weight=3).boundary_space()
            sage: [B3(x) for x in G.cusps()]
            [0, 0, 0, 0, 0, 0, 0, 0]
            sage: B3.rank()
            0

        Test that :issue:`25268` is fixed::

            sage: G = GammaH(20, [13])
            sage: A = ModularSymbols(G, 3)
            sage: S = A.cuspidal_submodule()
            sage: S.dimension()
            8
        """
        k = self.weight()
        sign = self.sign()
        i, eps = self._cusp_index(c)
        if i != -1:
            if i == -2:
                return self(0)
            else:
                return BoundarySpaceElement(self, {i: eps**k})

        if sign != 0:
            i2, eps = self._cusp_index(-c)
            if i2 != -1:
                if i2 == -2:
                    return self(0)
                else:
                    return BoundarySpaceElement(self, {i2: sign * (eps**k)})

        # found a new cusp class
        g = self._known_gens
        g.append(c)

        # Does the cusp class vanish because the cusp is irregular?

        if k % 2:
            if self.group().is_even() or not self.group().is_regular_cusp(c):
                self._zero_cusps.append(c)
                del self._known_gens[-1]
                return self(0)

        # Does class vanish because of sign relations?  The relevant
        # relations are
        #
        #    [(u,v)] = (-1)^k [(-u,-v)]
        #    [(u,v)] = sign * [(-u,v)]
        #    [(u,v)] = eps * (-1)^k [(-u,v)]
        #
        # where, in the last line, (u,v) is GammaH-equivalent to
        # (-u,v) or (u,-v) as eps is 1 or -1.
        #
        # Thus (other than for 0 and Infinity), we have that [(u,v)]
        # can only be killed by sign relations when:
        #
        #  - (u,v) is GammaH-equivalent to (-u,v) or (u,-v), and
        #  - eps is 1 and sign is -1, or eps is -1 and sign is not
        #    (-1)^k.
        #
        # (Notice that while this description looks identical to that
        # of Gamma1, it differs in that the condition of being GammaH
        # equivalent is weaker than that of being Gamma1 equivalent
        # when H is larger than {1}.)
        #
        if sign:
            if (c.is_infinity() and sign != (-1)**self.weight()) or \
               (c.is_zero() and sign == -1):
                self._zero_cusps.append(c)
                del self._known_gens[-1]
                return self(0)
            elif (not c.is_infinity() and not c.is_zero()):
                t, eps = self._is_equiv(c, -c)
                if t and ((eps == 1 and sign == -1) or
                          (eps == -1 and sign != (-1)**self.weight())):
                    self._zero_cusps.append(c)
                    del self._known_gens[-1]
                    return self(0)

        return BoundarySpaceElement(self, {(len(g) - 1): 1})


class BoundarySpace_wtk_eps(BoundarySpace):
    def __init__(self, eps, weight, sign=0):
        """
        Space of boundary modular symbols with given weight, character, and
        sign.

        INPUT:

        - ``eps`` -- dirichlet.DirichletCharacter, the
           "Nebentypus" character

        - ``weight`` -- integer; the weight = 2

        - ``sign`` -- integer; either -1, 0, or 1

        EXAMPLES::

            sage: B = ModularSymbols(DirichletGroup(6).0, 4).boundary_space() ; B
            Boundary Modular Symbols space of level 6, weight 4, character [-1] and dimension 0 over Rational Field
            sage: type(B)
            <class 'sage.modular.modsym.boundary.BoundarySpace_wtk_eps_with_category'>
            sage: B == loads(dumps(B))
            True
        """
        level = eps.modulus()
        sign = int(sign)
        self.__eps = eps
        if sign not in [-1, 0, 1]:
            raise ArithmeticError("sign must be an int in [-1,0,1]")
        if level <= 0:
            raise ArithmeticError("level must be positive")
        BoundarySpace.__init__(self,
                               weight=weight,
                               group=arithgroup.Gamma1(level),
                               sign=sign,
                               base_ring=eps.base_ring(),
                               character=eps)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: ModularSymbols(DirichletGroup(6).0, 4).boundary_space()._repr_()
            'Boundary Modular Symbols space of level 6, weight 4, character [-1] and dimension 0 over Rational Field'
        """
        return ("Boundary Modular Symbols space of level %s, weight %s, character %s " +
                "and dimension %s over %s") % (self.level(), self.weight(),
                                               self.character()._repr_short_(),
                                               self.rank(),
                                               self.base_ring())

    def _is_equiv(self, c1, c2):
        """
        Return a pair (b, t), where b is ``True`` if c1 and c2 are equivalent
        cusps for ``self``, and ``False`` otherwise, and t gives extra
        information about the equivalence of c1 and c2.

        EXAMPLES::

            sage: B = ModularSymbols(DirichletGroup(12).1, 3).boundary_space()
            sage: B._is_equiv(Cusp(0), Cusp(1/3))
            (False, None)
            sage: B._is_equiv(Cusp(2/3), Cusp(1/3))
            (True, 5)
            sage: B._is_equiv(Cusp(3/4), Cusp(1/4))
            (True, 7)
        """
        return c1.is_gamma0_equiv(c2, self.level(), transformation=True)

    def _cusp_index(self, cusp):
        """
        Return a pair (i, s), where i is the index of the first cusp in
        self._known_cusps() which is equivalent to cusp, and such that cusp is
        Gamma0-equivalent to self.character()(s) times self._known_cusps()[i].
        If cusp is not equivalent to any known cusp, return (-1, 0); if cusp is
        equivalent to a cusp known to vanish, return (-2, 0).

        EXAMPLES::

            sage: B = ModularSymbols(DirichletGroup(11).0**3, 5).boundary_space()
            sage: B._cusp_index(Cusp(0))
            (-1, 0)
            sage: B._coerce_cusp(Cusp(0))
            [0]
            sage: B._cusp_index(Cusp(0))
            (0, 1)
            sage: B._coerce_cusp(Cusp(1,11))
            [1/11]
            sage: B._cusp_index(Cusp(2,11))
            (1, -zeta10^2)
        """
        g = self._known_gens
        for i in range(len(g)):
            t, s = self._is_equiv(cusp, g[i])
            if t:
                return i, self.__eps(s)
        if any(self._is_equiv(cusp, c)[0] for c in self._zero_cusps):
            return (-2, 0)
        return -1, 0

    def _coerce_cusp(self, c):
        """
        Coerce the cusp c into ``self``.

        EXAMPLES::

            sage: B = ModularSymbols(DirichletGroup(13).0**3, 5, sign=0).boundary_space()
            sage: [ B(Cusp(i,13)) for i in range(13) ]
            [[0],
            [1/13],
            -zeta4*[1/13],
            [1/13],
            -[1/13],
            -zeta4*[1/13],
            -zeta4*[1/13],
            zeta4*[1/13],
            zeta4*[1/13],
            [1/13],
            -[1/13],
            zeta4*[1/13],
            -[1/13]]
            sage: B._is_equiv(Cusp(oo), Cusp(1,13))
            (True, 1)
            sage: B._is_equiv(Cusp(0), Cusp(1,13))
            (False, None)
            sage: B = ModularSymbols(DirichletGroup(13).0**3, 5, sign=1).boundary_space()
            sage: [ B(Cusp(i,13)) for i in range(13) ]
            [[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            sage: B._coerce_cusp(Cusp(oo))
            0
            sage: B = ModularSymbols(DirichletGroup(13).0**3, 5, sign=-1).boundary_space()
            sage: [ B(Cusp(i,13)) for i in range(13) ]
            [0,
            [1/13],
            -zeta4*[1/13],
            [1/13],
            -[1/13],
            -zeta4*[1/13],
            -zeta4*[1/13],
            zeta4*[1/13],
            zeta4*[1/13],
            [1/13],
            -[1/13],
            zeta4*[1/13],
            -[1/13]]
            sage: B = ModularSymbols(DirichletGroup(13).0**4, 5, sign=1).boundary_space()
            sage: B._coerce_cusp(Cusp(0))
            [0]
            sage: B = ModularSymbols(DirichletGroup(13).0**4, 5, sign=-1).boundary_space()
            sage: B._coerce_cusp(Cusp(0))
            0
        """
        sign = self.sign()
        i, eps = self._cusp_index(c)
        if i != -1:
            if i == -2:
                return self(0)
            else:
                return BoundarySpaceElement(self, {i: eps})

        if sign != 0:
            i2, eps = self._cusp_index(-c)
            if i2 != -1:
                if i2 == -2:
                    return self(0)
                else:
                    return BoundarySpaceElement(self, {i2: sign * eps})

        # found a new cusp class
        g = self._known_gens
        g.append(c)

        # Does cusp vanish because of stabiliser?
        s = arithgroup.Gamma0(self.level()).cusp_data(c)[0]
        if self.__eps(s[1, 1]) != 1:
            self._zero_cusps.append(c)
            del self._known_gens[-1]
            return self(0)

        # Does class vanish because of sign relations?  The relevant
        # relations are
        #
        #    [(u,v)] = (-1)^k [(-u,-v)]
        #    [(u,v)] = sign * [(-u,v)]
        #    [(u,v)] = eps(d) * [(-u,v)]
        #
        # where, in the last line, eps is the character defining
        # our space, and [a,b;c,d] takes (u,v) to (-u,v).
        #
        # Thus (other than for 0 and Infinity), we have that [(u,v)]
        # can only be killed by sign relations when the sign is not
        # equal to eps(d).
        #
        if sign:
            if c.is_zero():
                if sign == -1:
                    self._zero_cusps.append(c)
                    del self._known_gens[-1]
                    return self(0)
            elif c.is_infinity():
                if sign != (-1)**self.weight():
                    self._zero_cusps.append(c)
                    del self._known_gens[-1]
                    return self(0)
            else:
                t, s = self._is_equiv(c, -c)
                if t:
                    if sign != self.__eps(s):
                        self._zero_cusps.append(c)
                        del self._known_gens[-1]
                        return self(0)

        return BoundarySpaceElement(self, {(len(g) - 1): 1})
