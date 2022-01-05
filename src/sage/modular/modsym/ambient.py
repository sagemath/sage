# -*- coding: utf-8 -*-
r"""
Ambient spaces of modular symbols

This module defines the following classes.  There is an abstract base
class ``ModularSymbolsAmbient``, derived from
``space.ModularSymbolsSpace`` and ``hecke.AmbientHeckeModule``.  As
this is an abstract base class, only derived classes should be
instantiated.  There are five derived classes:

- ``ModularSymbolsAmbient_wtk_g0``, for modular symbols of general
  weight `k` for `\Gamma_0(N)`;

- ``ModularSymbolsAmbient_wt2_g0`` (derived from
  ``ModularSymbolsAmbient_wtk_g0``), for modular symbols of weight 2
  for `\Gamma_0(N)`;

- ``ModularSymbolsAmbient_wtk_g1``, for modular symbols of general
  weight `k` for `\Gamma_1(N)`;

- ``ModularSymbolsAmbient_wtk_gamma_h``, for modular symbols of
  general weight `k` for `\Gamma_H`, where `H` is a subgroup of
  `\ZZ/N\ZZ`;

- ``ModularSymbolsAmbient_wtk_eps``, for modular symbols of general
  weight `k` and character `\epsilon`.



EXAMPLES:

We compute a space of modular symbols modulo 2. The dimension is
different from that of the corresponding space in characteristic
0::

    sage: M = ModularSymbols(11,4,base_ring=GF(2)); M
    Modular Symbols space of dimension 7 for Gamma_0(11) of weight 4
    with sign 0 over Finite Field of size 2
    sage: M.basis()
    ([X*Y,(1,0)], [X*Y,(1,8)], [X*Y,(1,9)], [X^2,(0,1)], [X^2,(1,8)], [X^2,(1,9)], [X^2,(1,10)])
    sage: M0 =  ModularSymbols(11,4,base_ring=QQ); M0
    Modular Symbols space of dimension 6 for Gamma_0(11) of weight 4
    with sign 0 over Rational Field
    sage: M0.basis()
    ([X^2,(0,1)], [X^2,(1,6)], [X^2,(1,7)], [X^2,(1,8)], [X^2,(1,9)], [X^2,(1,10)])

The characteristic polynomial of the Hecke operator `T_2` has an extra
factor `x`.

::

    sage: M.T(2).matrix().fcp('x')
    (x + 1)^2 * x^5
    sage: M0.T(2).matrix().fcp('x')
    (x - 9)^2 * (x^2 - 2*x - 2)^2
"""

################################################################################
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
################################################################################
# Sage packages
from sage.misc.cachefunc import cached_method
import sage.misc.latex as latex
from sage.misc.verbose import verbose

import sage.matrix.matrix_space as matrix_space
from sage.modular.arithgroup.arithgroup_element import M2Z
import sage.modules.free_module_element as free_module_element
import sage.modules.free_module as free_module
import sage.modular.arithgroup.all as arithgroup
import sage.modular.dirichlet as dirichlet
import sage.modular.hecke.all as hecke
from sage.rings.all import Integer, QQ, ZZ, Ring
from sage.arith.all import is_prime, divisors, number_of_divisors, crt
import sage.rings.polynomial.multi_polynomial_element
import sage.structure.formal_sum as formal_sum
import sage.categories.all as cat
from sage.modular.cusps import Cusp
from sage.modular.modsym.apply import apply_to_monomial
from sage.modular.modsym.manin_symbol import ManinSymbol
from sage.modular.modsym.manin_symbol_list import (ManinSymbolList_gamma0,
                                                   ManinSymbolList_gamma1,
                                                   ManinSymbolList_gamma_h,
                                                   ManinSymbolList_character)


from . import boundary
from . import element
from . import heilbronn
from . import modular_symbols
from . import modsym
from . import p1list
from . import relation_matrix
from .space import ModularSymbolsSpace
from . import subspace


class ModularSymbolsAmbient(ModularSymbolsSpace, hecke.AmbientHeckeModule):
    r"""
    An ambient space of modular symbols for a congruence subgroup of
    `SL_2(\ZZ)`.

    This class is an abstract base class, so only derived classes
    should be instantiated.

    INPUT:

    - ``weight`` - an integer
    - ``group`` - a congruence subgroup.
    - ``sign`` - an integer, either -1, 0, or 1
    - ``base_ring`` - a commutative ring
    - ``custom_init`` - a function that is called with self as input
      before any computations are done using self; this could be used
      to set a custom modular symbols presentation.

    TESTS::

        sage: ModularSymbols(11,2) == ModularSymbols(11,2)
        True
        sage: ModularSymbols(11,2) == ModularSymbols(11,4)
        False
        sage: ModularSymbols(11,2) != ModularSymbols(11,2)
        False
        sage: ModularSymbols(11,2) != ModularSymbols(11,4)
        True
        sage: hash(ModularSymbols(11,2)) != hash(ModularSymbols(11,4))
        True
    """
    def __init__(self, group, weight, sign, base_ring,
                 character=None, custom_init=None, category=None):
        """
        Initialize a space of modular symbols.

        INPUT:

        -  ``weight`` - an integer

        -  ``group`` - a congruence subgroup.

        -  ``sign`` - an integer, either -1, 0, or 1

        -  ``base_ring`` - a commutative ring

        EXAMPLES::

            sage: ModularSymbols(2,2)
            Modular Symbols space of dimension 1 for Gamma_0(2) of weight 2 with sign 0 over Rational Field

        """
        weight = int(weight)
        if weight <= 1:
            raise ValueError("Weight (=%s) Modular symbols of weight <= 1 not defined."%weight)
        if not arithgroup.is_CongruenceSubgroup(group):
            raise TypeError("group must be a congruence subgroup")

        sign = int(sign)
        if not isinstance(base_ring, Ring) and base_ring.is_field():
            raise TypeError("base_ring must be a commutative ring")

        if character is None and arithgroup.is_Gamma0(group):
            character = dirichlet.TrivialCharacter(group.level(), base_ring)

        ModularSymbolsSpace.__init__(self, group, weight,
                                           character, sign, base_ring,
                                           category=category)

        if custom_init is not None:
            custom_init(self)

        try:
            formula = self._dimension_formula()
        except NotImplementedError:
            formula = None

        rank = self.rank()
        if formula is not None:
            assert rank == formula, \
                   "Computed dimension (=%s) of ambient space \"%s\" doesn't match dimension formula (=%s)!\n"%(rank, self, formula) + \
                   "ModularSymbolsAmbient: group = %s, weight = %s, sign = %s, base_ring = %s, character = %s"%(
                         group, weight, sign, base_ring, character)

        hecke.AmbientHeckeModule.__init__(self, base_ring, rank, group.level(), weight, category=category)

    def new_submodule(self, p=None):
        r"""
        Return the new or `p`-new submodule of this modular symbols ambient space.

        INPUT:


        -  ``p`` - (default: None); if not None, return only
           the `p`-new submodule.


        OUTPUT:

        The new or `p`-new submodule of this modular symbols ambient space.

        EXAMPLES::

            sage: ModularSymbols(100).new_submodule()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 31 for Gamma_0(100) of weight 2 with sign 0 over Rational Field
            sage: ModularSymbols(389).new_submodule()
            Modular Symbols space of dimension 65 for Gamma_0(389) of weight 2 with sign 0 over Rational Field
        """
        # Check for special cases where the answer is easy.
        # If not in one of those cases, use the generic code.
        if self.level().is_prime() and self.weight() == 2:
            return self
        return hecke.AmbientHeckeModule.new_submodule(self, p=p)

    def manin_symbols(self):
        """
        Return the list of Manin symbols for this modular symbols ambient space.

        EXAMPLES::

            sage: ModularSymbols(11,2).manin_symbols()
            Manin Symbol List of weight 2 for Gamma0(11)
        """
        raise NotImplementedError

    def manin_generators(self):
        """
        Return list of all Manin symbols for this space. These are the
        generators in the presentation of this space by Manin symbols.

        EXAMPLES::

            sage: M = ModularSymbols(2,2)
            sage: M.manin_generators()
            [(0,1), (1,0), (1,1)]

        ::

            sage: M = ModularSymbols(1,6)
            sage: M.manin_generators()
            [[Y^4,(0,0)], [X*Y^3,(0,0)], [X^2*Y^2,(0,0)], [X^3*Y,(0,0)], [X^4,(0,0)]]
        """
        return self._manin_generators

    def manin_basis(self):
        r"""
        Return a list of indices into the list of Manin generators (see
        ``self.manin_generators()``) such that those symbols
        form a basis for the quotient of the `\QQ`-vector
        space spanned by Manin symbols modulo the relations.

        EXAMPLES::

            sage: M = ModularSymbols(2,2)
            sage: M.manin_basis()
            [1]
            sage: [M.manin_generators()[i] for i in M.manin_basis()]
            [(1,0)]
            sage: M = ModularSymbols(6,2)
            sage: M.manin_basis()
            [1, 10, 11]
            sage: [M.manin_generators()[i] for i in M.manin_basis()]
            [(1,0), (3,1), (3,2)]
        """
        try:
            return self._manin_basis
        except AttributeError:
            self.compute_presentation()
        return self._manin_basis

    def p1list(self):
        """
        Return a P1list of the level of this modular symbol space.

        EXAMPLES::

            sage: ModularSymbols(11,2).p1list()
            The projective line over the integers modulo 11
        """
        try:
            return self.__p1list
        except AttributeError:
            self.__p1list = p1list.P1List(self.level())
        return self.__p1list

#     See the file relation_matrix.py
#
#     def relation_matrix(self):
#         raise NotImplementedError

    def compute_presentation(self):
        r"""
        Compute and cache the presentation of this space.

        EXAMPLES::

            sage: ModularSymbols(11,2).compute_presentation() # no output

        """
        B, basis, mod = relation_matrix.compute_presentation(
                self.manin_symbols(), self.sign(),
                self.base_ring())
        self._manin_generators = self.manin_symbols().manin_symbol_list()
        self._manin_basis = basis
        self._manin_gens_to_basis = B
        self._mod2term = mod

    def manin_gens_to_basis(self):
        r"""
        Return the matrix expressing the manin symbol generators in terms of the basis.

        EXAMPLES::

            sage: ModularSymbols(11,2).manin_gens_to_basis()
            [-1  0  0]
            [ 1  0  0]
            [ 0  0  0]
            [ 0  0  1]
            [ 0 -1  1]
            [ 0 -1  0]
            [ 0  0 -1]
            [ 0  0 -1]
            [ 0  1 -1]
            [ 0  1  0]
            [ 0  0  1]
            [ 0  0  0]
        """
        try:
            return self._manin_gens_to_basis
        except AttributeError:
            self.compute_presentation()
            return self._manin_gens_to_basis


    #####################################################################
    # Coercion
    #####################################################################
    def _element_constructor_(self, x, computed_with_hecke=False):
        r"""
        Coerce `x` into this modular symbols space. The result is
        either an element of self or a subspace of self.

        INPUT:

        The allowed input types for `x` are as follows:


        -  ``Vector`` - a vector of the same degree. This
           defines the corresponding linear combination of the basis of self.

        -  ``ManinSymbol`` - a Manin symbol of the same weight
           as the space

        -  ``ModularSymbolsElement`` - a modular symbol whose
           ambient parent is this space of modular symbols. (TODO: make more
           sophisticated)

        -  0 - the integer 0; results in the 0 modular symbol.

        -  3-tuple - Given a 3-tuple (i,u,v), returns the modular symbol
           element defined by the Manin symbol
           `[X^{i}\cdot Y^{k-2-i}, (u,v)]`, where k is the weight.
           Note that we must have `0\leq i \leq k-2`.

        -  2-tuple - Given a 2-tuple (u,v), returns the element defined by
           the Manin symbol `[X^0 \cdot Y^{2-k}, (u,v)]`.

        -  2-elements list - Given a list ``[alpha, beta]``,
           where `\alpha` and `\beta` are (coercible to)
           cusps, return the modular symbol `\{\alpha, \beta\}`. When
           the weight `k > 2` return
           `Y^{k-2} \{\alpha, \beta\}`.

        -  3-element list - Given a list ``[i, alpha, beta]``,
           where `i` is an integer, and `\alpha`,
           `\beta` are (coercible to) cusps, return the modular symbol
           `X^i Y^{k-2-i} \{\alpha, \beta\}`.

           If our list is ``[f, alpha, beta]``, where `f`
           is a homogeneous polynomial in two variables of degree k-2 with
           integer coefficients, and alpha and beta are cusps, return the
           corresponding sum of modular symbols as an element of self. So if
           `f = \sum_{i=0}^{k-2} a_i X^i Y^{k-2-i}`, return
           `\sum_{i=0}^{k-2} a_i * [ i, alpha, beta ]`.

        EXAMPLES::

            sage: M = ModularSymbols(37,2)

        M(0) is the 0 element of the space::

            sage: M(0)
            0
            sage: type(M(0))
            <class 'sage.modular.modsym.ambient.ModularSymbolsAmbient_wt2_g0_with_category.element_class'>

        From a vector of the correct dimension we construct the
        corresponding linear combination of the basis elements::

            sage: M.dimension()
            5
            sage: M.basis()
            ((1,0), (1,23), (1,32), (1,34), (1,35))
            sage: M(vector([1,2,3,4,5]))
            (1,0) + 2*(1,23) + 3*(1,32) + 4*(1,34) + 5*(1,35)
            sage: M(vector([1/2,2/3,3/4,4/5,5/6]))
            1/2*(1,0) + 2/3*(1,23) + 3/4*(1,32) + 4/5*(1,34) + 5/6*(1,35)

        Manin symbols can be converted to elements of the space::

            sage: S = M.manin_symbols()
            sage: S((0,2,3))
            (2,3)
            sage: M( S((0,2,3)) )
            (1,34) - (1,35)

        However, it is easier to use one of the following forms.
        Either a 3-tuple `(i,u,v)` or a 2-tuple `(u,v)` with `i=0`
        assumed::

            sage: M((0,2,3))
            (1,34) - (1,35)
            sage: M((2,3))
            (1,34) - (1,35)

        Or a 3-list `[i,\alpha,\beta]` where `i` is the degree and
        `\alpha` and `beta` are cusps, or a 2-tuple `[\alpha,\beta]`
        with `i=0` assumed::

            sage: M([0,Cusp(1/2),Cusp(0)])
            (1,35)
            sage: M([Cusp(1/2),Cusp(0)])
            (1,35)


        """
        if isinstance(x, free_module_element.FreeModuleElement):
            if x.degree() != self.dimension():
                raise TypeError("Incompatible degrees: x has degree %s\
                    but modular symbols space has dimension %s"%(
                    x.degree(), self.dimension()))
            #if x.parent().base_ring() != self.base_ring():
            #    raise TypeError, "Vector x is over %s, but modular symbols space is over %s."%(
            #        x.parent().base_ring(), self.base_ring())
            return self.element_class(self, x)

        elif isinstance(x, (ManinSymbol, element.ModularSymbolsElement)):
            return self.element(x)

        elif isinstance(x, modular_symbols.ModularSymbol):
            return self(x.manin_symbol_rep())

        elif isinstance(x, (int, Integer)) and x==0:
            return self.element_class(self, self.free_module()(0))

        elif isinstance(x, tuple):
            return self.manin_symbol(x)

        elif isinstance(x, formal_sum.FormalSum):
            return sum([c*self(y) for c, y in x], self(0))

        elif isinstance(x, list):
            if len(x) == 3 and sage.rings.polynomial.multi_polynomial_element.is_MPolynomial(x[0]):
                return self.modular_symbol_sum(x)
            else:
                return self.modular_symbol(x)

        raise TypeError("No coercion of %s into %s defined."%(x, self))


    def change_ring(self, R):
        r"""
        Change the base ring to R.

        EXAMPLES::

            sage: ModularSymbols(Gamma1(13), 2).change_ring(GF(17))
            Modular Symbols space of dimension 15 for Gamma_1(13) of weight 2 with sign 0 over Finite Field of size 17
            sage: M = ModularSymbols(DirichletGroup(5).0, 7); MM=M.change_ring(CyclotomicField(8)); MM
            Modular Symbols space of dimension 6 and level 5, weight 7, character [zeta8^2], sign 0, over Cyclotomic Field of order 8 and degree 4
            sage: MM.change_ring(CyclotomicField(4)) == M
            True
            sage: M.change_ring(QQ)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce zeta4 to a rational

        Similarly with :meth:`base_extend`::

            sage: M = ModularSymbols(DirichletGroup(5).0, 7); MM = M.base_extend(CyclotomicField(8)); MM
            Modular Symbols space of dimension 6 and level 5, weight 7, character [zeta8^2], sign 0, over Cyclotomic Field of order 8 and degree 4
            sage: MM.base_extend(CyclotomicField(4))
            Traceback (most recent call last):
            ...
            TypeError: Base extension of self (over 'Cyclotomic Field of order 8 and degree 4') to ring 'Cyclotomic Field of order 4 and degree 2' not defined.
        """
        if self.character() is None:
            return modsym.ModularSymbols(self.group(), self.weight(), self.sign(), R)
        else:
            return modsym.ModularSymbols(self.character(), self.weight(), self.sign(), R)

    def _action_on_modular_symbols(self, g):
        r"""
        Return the matrix of the action of a 2x2 matrix on this space.

        INPUT:

        `g` (list) -- `g=[a,b,c,d]` where `a,b,c,d` are integers
        defining a `2\times2` integer matrix.

        OUTPUT:

        (matrix) The matrix of the action of `g` on this Modular
        Symbol space, with respect to the standard basis.

        .. NOTE::

            Use ``_matrix_of_operator_on_modular_symbols`` for more general
            operators.

        EXAMPLES::

            sage: M = ModularSymbols(11,4,1)
            sage: M._action_on_modular_symbols([1,2,3,7])
            [0 0 1 0]
            [0 0 0 1]
            [0 1 0 0]
            [0 1 0 0]

        """
        if not isinstance(g, list):
            raise TypeError("g must be a list")
        if not len(g) == 4:
            raise TypeError("g must be a list of length 4")
        return self._matrix_of_operator_on_modular_symbols(self, [g])

    def manin_symbol(self, x, check=True):
        r"""
        Construct a Manin Symbol from the given data.

        INPUT:

        - ``x`` (list) -- either `[u,v]` or `[i,u,v]`, where `0\le
          i\le k-2` where `k` is the weight, and `u`,`v` are integers
          defining a valid element of `\mathbb{P}^1(N)`, where `N` is
          the level.

        OUTPUT:

        (ManinSymbol) the Manin Symbol associated to `[i;(u,v)]`, with
        `i=0` if not supplied, corresponding to the monomial symbol
        `[X^i*Y^{k-2-i}, (u,v)]`.

        EXAMPLES::

            sage: M = ModularSymbols(11,4,1)
            sage: M.manin_symbol([2,5,6])
            -2/3*[X^2,(1,6)] + 5/3*[X^2,(1,9)]
        """
        if check:
            if len(x) == 2:
                x = (0,x[0],x[1])
            if len(x) == 3:
                # Manin symbol of the form (i, u, v), which corresponds to [X^i*Y^(k-2-i), (u,v)].
                if x[0] < 0 or x[0] > self.weight()-2:
                    raise ValueError("The first entry of the tuple (=%s)\
                        must be an integer between 0 and k-2 (=%s)."%(
                        x, self.weight()-2))
            else:
                raise ValueError("x (=%s) must be of length 2 or 3"%x)
        # end check

        N = self.level()
        x = (x[0], x[1]%N, x[2]%N)
        try:
            return self.__manin_symbol[x]
        except AttributeError:
            self.__manin_symbol = {}
        except KeyError:
            pass
        y = self.manin_symbols()(x)
        z = self(y)
        self.__manin_symbol[x] = z
        return z

    def _modular_symbol_0_to_alpha(self, alpha, i=0):
        r"""
        Return the modular symbol `\{0,\alpha\}` in this space.

        INPUT:

        - ``alpha`` (rational or Infinity) -- a cusp

        - ``i`` (int, default 0) -- the degree of the symbol.

        OUTPUT:

        (ModularSymbol) The modular symbol `X^iY^{k-2-i}\{0,\alpha\}`.

        EXAMPLES::

            sage: M = ModularSymbols(11,4,1)
            sage: M._modular_symbol_0_to_alpha(Cusp(3/5))
            11*[X^2,(1,4)] + 40/3*[X^2,(1,6)] - 1/3*[X^2,(1,9)]
            sage: M._modular_symbol_0_to_alpha(Cusp(3/5),1)
            15/2*[X^2,(1,4)] + 20/3*[X^2,(1,6)] + 5/6*[X^2,(1,9)]
            sage: M._modular_symbol_0_to_alpha(Cusp(Infinity))
            2/3*[X^2,(1,6)] - 5/3*[X^2,(1,9)]
            sage: M._modular_symbol_0_to_alpha(Cusp(Infinity),1)
            0
        """
        if alpha.is_infinity():
            return self.manin_symbol((i,0,1), check=False)
        # v, c = arith.continued_fraction_list(alpha._rational_(), partial_convergents=True)
        cf = alpha._rational_().continued_fraction()
        v = list(cf)
        c = [(cf.p(k),cf.q(k)) for k in range(len(cf))]
        a = self(0)
        one = ZZ.one()
        two = ZZ(2)
        if self.weight() > two:
            R = ZZ['X']
            X = R.gen(0)
            ## need to add first two terms, which aren't necessarily
            ## zero in this case. we do the first here, and the
            ## second in the k=0 case below, so as to avoid code
            ## duplication
            a += self.manin_symbol((i,0,1), check=False)
            for k in range(0,len(c)):
                ## matrix entries associated to this partial sum
                if k == 0:
                    x = c[0][0]
                    y = -1
                    z = 1
                    w = 0
                else:
                    x = c[k][0]
                    y = c[k-1][0]
                    z = c[k][1]
                    w = c[k-1][1]
                    if k%2 == 0:
                        y = -y
                        w = -w

                ## two options here: write out the polynomial directly,
                ## and deal with all the separate cases, or create two
                ## polynomials and then exponentiate and multiply them.
                ## given how fast ntl/flint/etc are, the second may
                ## be faster.

                ## method 1: write out solution. this is currently
                ## incorrect, because it ends up doing 0^0 in the sum,
                ## so I'll fix it and do timings soon.
##                for s in range(0,self.weight()-two+1):
##                    coeff = sum([ binomial(i,t)*binomial(self.weight()-two-i,s-t)*
##                                  x**t * y**(i-t) * z**(s-t) *
##                                  w**(self.weight()-two-i-s+t) for t in range(0,s) ])
##                    m = coeff * self.manin_symbol((s, y, w), check=False)
##                    a += m

                ## method 2
                p1 = x*X+y
                p2 = z*X+w
                if i == 0:
                    p1 = R(one)
                if (self.weight()-2-i == 0):
                    p2 = R(one)
                poly = (p1**i) * (p2**(self.weight()-2-i))
                for s in range(0,self.weight()-1):  # k-2+1 = k-1
                    a += poly[s] * self.manin_symbol((s,z,w), check=False)
        else:
            for k in range(1,len(c)):
                u = c[k][1]
                v = c[k-1][1]
                if k % 2 == 0:
                    v = -v
                x = self.manin_symbol((i, u, v), check=False)
                a += x
        return a

    def modular_symbol(self, x, check=True):
        r"""
        Create a modular symbol in this space.

        INPUT:

        -  ``x`` (list) -- a list of either 2 or 3 entries:

            - 2 entries: `[\alpha, \beta]` where `\alpha` and `\beta`
              are cusps;

            - 3 entries: `[i, \alpha, \beta]` where `0\le i\le k-2`
              and `\alpha` and `\beta` are cusps;

        - ``check`` (bool, default True) -- flag that determines
          whether the input ``x`` needs processing: use check=False
          for efficiency if the input ``x`` is a list of length 3 whose
          first entry is an Integer, and whose second and third
          entries are Cusps (see examples).

        OUTPUT:

        (Modular Symbol) The modular symbol `Y^{k-2}\{\alpha,
        \beta\}`. or `X^i Y^{k-2-i}\{\alpha,\beta\}`.

        EXAMPLES::

            sage: set_modsym_print_mode('modular')
            sage: M = ModularSymbols(11)
            sage: M.modular_symbol([2/11, oo])
            -{-1/9, 0}
            sage: M.1
            {-1/8, 0}
            sage: M.modular_symbol([-1/8, 0])
            {-1/8, 0}
            sage: M.modular_symbol([0, -1/8, 0])
            {-1/8, 0}
            sage: M.modular_symbol([10, -1/8, 0])
            Traceback (most recent call last):
            ...
            ValueError: The first entry of the tuple (=[10, -1/8, 0]) must be an integer between 0 and k-2 (=0).

        ::

            sage: N = ModularSymbols(6,4)
            sage: set_modsym_print_mode('manin')
            sage: N([1,Cusp(-1/4),Cusp(0)])
            17/2*[X^2,(2,3)] - 9/2*[X^2,(2,5)] + 15/2*[X^2,(3,1)] - 15/2*[X^2,(3,2)]
            sage: N([1,Cusp(-1/2),Cusp(0)])
            1/2*[X^2,(2,3)] + 3/2*[X^2,(2,5)] + 3/2*[X^2,(3,1)] - 3/2*[X^2,(3,2)]

        Use check=False for efficiency if the input x is a list of length 3
        whose first entry is an Integer, and whose second and third entries
        are cusps::

            sage: M.modular_symbol([0, Cusp(2/11), Cusp(oo)], check=False)
            -(1,9)

        ::

            sage: set_modsym_print_mode()   # return to default.
        """

        if check:
            if len(x) == 2:
                x = [0,x[0],x[1]]
            elif len(x) == 3:
                if x[0] < 0 or x[0] > self.weight()-2:
                    raise ValueError("The first entry of the tuple (=%s)\
                        must be an integer between 0 and k-2 (=%s)."%(
                        x, self.weight()-2))
            else:
                raise ValueError("x (=%s) must be of length 2 or 3"%x)
            i = Integer(x[0])
            alpha = Cusp(x[1])
            beta = Cusp(x[2])
        else:
            i = x[0]
            alpha = x[1]
            beta = x[2]

        # Compute {0,beta} - {0,alpha}
        a = self._modular_symbol_0_to_alpha(alpha, i)
        b = self._modular_symbol_0_to_alpha(beta, i)
        return b - a

    def modular_symbol_sum(self, x, check=True):
        r"""
        Construct a modular symbol sum.

        INPUT:

        - ``x`` (list) -- `[f, \alpha, \beta]` where `f =
          \sum_{i=0}^{k-2} a_i X^i Y^{k-2-i}` is a homogeneous
          polynomial over `\ZZ` of degree `k` and `\alpha` and `\beta`
          are cusps.

        - ``check`` (bool, default True) -- if True check the validity
          of the input tuple ``x``

        OUTPUT:

        The sum `\sum_{i=0}^{k-2} a_i [ i, \alpha, \beta ]` as an
        element of this modular symbol space.

        EXAMPLES::

            sage: M = ModularSymbols(11,4)
            sage: R.<X,Y>=QQ[]
            sage: M.modular_symbol_sum([X*Y,Cusp(0),Cusp(Infinity)])
            -3/14*[X^2,(1,6)] + 1/14*[X^2,(1,7)] - 1/14*[X^2,(1,8)] + 1/2*[X^2,(1,9)] - 2/7*[X^2,(1,10)]
        """
        if check:
            if len(x) != 3:
                raise ValueError("%s must have length 3"%x)
            f = x[0]
            R = self.base_ring()['X','Y']
            X = R.gen(0)
            try:
                f = R(f)
            except TypeError:
                raise ValueError("f must be coercible to a polynomial \
                    over %s"%self.base_ring())
            if (not f.is_homogeneous()) or (f.degree() != self.weight()-2):
                raise ValueError("f must be a homogeneous polynomial of degree k-2")
            alpha = Cusp(x[1])
            beta = Cusp(x[2])
        else:
            f = x[0]
            R = self.base_ring()
            X = R.gen(0)
            alpha = x[1]
            beta = x[2]

        s = self(0)

        for term in f.monomials():
            deg = term.degree(X)
            a = self._modular_symbol_0_to_alpha(alpha, deg)
            b = self._modular_symbol_0_to_alpha(beta, deg)
            s += f.monomial_coefficient(term) * (b - a)

        return s

    def _compute_dual_hecke_matrix(self, n):
        r"""
        Return the matrix of the dual Hecke operator `T(n)`.

        INPUT:

        - ``n`` (int) -- a positive integer

        OUTPUT:

        (matrix) The matrix of the dual of `T(n)`.

        EXAMPLES::

            sage: M = ModularSymbols(11,4,1)
            sage: M._compute_dual_hecke_matrix(5)
            [  126     0     0     0]
            [    2    63    62    38]
            [ 26/3   -20   -27    -4]
            [-32/3    83    91    92]
        """
        return self.hecke_matrix(n).transpose()

    def _compute_hecke_matrix_prime(self, p, rows=None):
        """
        Return the matrix of the Hecke operator `T(p)`.

        INPUT:

        - ``p`` (int) -- a prime number.

        - ``rows`` (list or None (default)) -- if not None, a list of
          the rows which should be computed; otherwise the complete
          matrix will be computed,

        .. note::

           `p` does not have to be, prime despite the function name.

        OUTPUT:

        (matrix) The matrix of the Hecke operator `T(p)` on this
        space, with respect to its standard basis.

        ALGORITHM:

        Use Heilbronn-Cremona matrices if `p` is prime, else use
        Heilbronn-Merel matrices.

        EXAMPLES:

        We first compute some examples for Gamma0(N)::

            sage: m = ModularSymbols(2, weight=4)
            sage: m._compute_hecke_matrix_prime(2).charpoly('x')
            x^2 - 9*x + 8

        ::

            sage: m = ModularSymbols(1,weight=12)
            sage: m._compute_hecke_matrix_prime(2).charpoly('x')
            x^3 - 2001*x^2 - 97776*x - 1180224
            sage: m._compute_hecke_matrix_prime(13).charpoly('x')
            x^3 - 1792159238562*x^2 - 2070797989680255444*x - 598189440899986203208472

        ::

            sage: m = ModularSymbols(1,weight=12, sign=-1)
            sage: m._compute_hecke_matrix_prime(5)
            [4830]
            sage: m._compute_hecke_matrix_prime(23)
            [18643272]

        ::

            sage: m = ModularSymbols(3,4)
            sage: m._compute_hecke_matrix_prime(2).charpoly('x')
            x^2 - 18*x + 81

        ::

            sage: m = ModularSymbols(6,4)
            sage: m._compute_hecke_matrix_prime(2).charpoly('x')
            x^6 - 14*x^5 + 29*x^4 + 172*x^3 - 124*x^2 - 320*x + 256
            sage: m._compute_hecke_matrix_prime(3).charpoly('x')
            x^6 - 50*x^5 + 511*x^4 + 3012*x^3 - 801*x^2 - 9234*x + 6561

        ::

            sage: m = ModularSymbols(15,4, sign=-1)
            sage: m._compute_hecke_matrix_prime(3).charpoly('x')
            x^4 - 2*x^3 + 18*x^2 + 18*x - 243

        ::

            sage: m = ModularSymbols(6,4)
            sage: m._compute_hecke_matrix_prime(7).charpoly('x')
            x^6 - 1344*x^5 + 666240*x^4 - 140462080*x^3 + 8974602240*x^2 + 406424518656*x + 3584872677376

        ::

            sage: m = ModularSymbols(4,4)
            sage: m._compute_hecke_matrix_prime(3).charpoly('x')
            x^3 - 84*x^2 + 2352*x - 21952

        We now compute some examples for modular symbols on Gamma1(N)::

            sage: m = ModularSymbols(Gamma1(13),2, sign=-1)
            sage: m._compute_hecke_matrix_prime(2).charpoly('x')
            x^2 + 3*x + 3

        The following is an example with odd weight::

            sage: m = ModularSymbols(Gamma1(5),3)
            sage: m._compute_hecke_matrix_prime(2).charpoly('x')
            x^4 - 10*x^3 + 50*x^2 - 170*x + 289

        This example has composite conductor and weight2 dividing the
        conductor and nontrivial sign::

            sage: m = ModularSymbols(Gamma1(9),3, sign=1)
            sage: m._compute_hecke_matrix_prime(3).charpoly('x')
            x^6 + 3*x^4 - 19*x^3 + 24*x^2 - 9*x

        In some situations we do not need all the rows of the result, and can thereby save time::

            sage: m = ModularSymbols(1,weight=12)
            sage: m._compute_hecke_matrix_prime(2)
            [ -24    0    0]
            [   0  -24    0]
            [4860    0 2049]
            sage: m._compute_hecke_matrix_prime(2,rows=[0,1])
            [-24   0   0]
            [  0 -24   0]
            sage: m._compute_hecke_matrix_prime(2,rows=[1,2])
            [   0  -24    0]
            [4860    0 2049]
        """
        # note -- p doesn't have to be prime despite the function name
        p = int(Integer(p))   # go through Integer so p = 2.5 gives an error.
        if isinstance(rows, list):
            rows = tuple(rows)
        try:
            return self._hecke_matrices[(p, rows)]
        except AttributeError:
            self._hecke_matrices = {}
        except KeyError:
            pass
        tm = verbose("Computing Hecke operator T_%s" % p)

        if is_prime(p):
            H = heilbronn.HeilbronnCremona(p)
        else:
            H = heilbronn.HeilbronnMerel(p)

        B = self.manin_basis()
        if rows is not None:
            B = [B[i] for i in rows]
        mod2term = self._mod2term
        R = self.manin_gens_to_basis()
        K = self.base_ring()
        W = R.new_matrix(nrows=len(B), ncols=R.nrows())
        syms = self.manin_symbols()
        j = 0
        for i in B:
            for h in H:
                entries = syms.apply(i,h)
                for k, x in entries:
                    f, s = mod2term[k]
                    if s:
                        # W[j,f] = W[j,f] + s*K(x)
                        W.add_to_entry(j, f, s * K(x))
            j += 1
        tm = verbose("start matrix multiply",tm)
        if hasattr(W, '_matrix_times_matrix_dense'):
            Tp = W._matrix_times_matrix_dense(R)
            verbose("done matrix multiply and computing Hecke operator",tm)
        else:
            Tp = W * R
            tm = verbose("done matrix multiply",tm)
            Tp = Tp.dense_matrix()
            verbose("done making Hecke operator matrix dense",tm)
        self._hecke_matrices[(p,rows)] = Tp
        return Tp


    def __heilbronn_operator(self, M, H, t=1):
        r"""
        Return the matrix function to the space `M` defined by `H`, `t`.

        .. note::

           Users will instead use the simpler interface defined, for
           example, by ``hecke_matrix()`` (see examples).

        INPUT:


        -  ``M`` (ModularSymbols) -- codomain (a space of modular
           symbols);

        -  ``H`` (list) -- a list of matrices in `M_2(\ZZ)`;

        -  ``t`` (int, default 1) -- an integer.


        OUTPUT:

        (free module morphism) A function from the Modular Symbol
         space to the Modular Symbol space `M` defined by `t` and the
         matrices in `H`.

         EXAMPLES::

             sage: M = ModularSymbols(37,2)
             sage: M._ModularSymbolsAmbient__heilbronn_operator(M,HeilbronnCremona(3))
             Hecke module morphism Heilbronn operator(The Cremona-Heilbronn matrices of determinant 3,1) defined by the matrix
             [ 4  0  0  0 -1]
             [ 0 -1  2  2 -2]
             [ 0  2 -1  2  0]
             [ 0  0  0 -3  2]
             [ 0  0  0  0  1]
             Domain: Modular Symbols space of dimension 5 for Gamma_0(37) of weight ...
             Codomain: Modular Symbols space of dimension 5 for Gamma_0(37) of weight ...

             sage: M.hecke_matrix(3)
             [ 4  0  0  0 -1]
             [ 0 -1  2  2 -2]
             [ 0  2 -1  2  0]
             [ 0  0  0 -3  2]
             [ 0  0  0  0  1]


        """

        MS = matrix_space.MatrixSpace(self.base_ring(), self.dimension(), M.dimension())
        hom = self.Hom(M)
        if self.dimension() == 0 or M.dimension() == 0:
            A = MS(0)
            phi = hom(A, "Heilbronn operator(%s,%s)"%(H,t))
            return phi

        rows = []
        B = self.manin_basis()
        syms = self.manin_symbols()
        k = self.weight()
        for n in B:
            z = M(0)
            i, u, v = syms[n].tuple()
            # We apply each Heilbronn matrix to the
            #    Manin symbol [X^i*Y^(k-2-i), (u,v)]
            for h in H:
                # Apply h to the polynomial part
                (a,b,c,d) = tuple(h)
                # P gives the ordered coefficients of (a*X+b*Y)^i*(c*X+d*Y)^(j-i)
                P = apply_to_monomial(i, k-2, a,b,c,d)
                # Apply h to the (u,v) part of the Manin symbol
                (uu,vv) = (u*a+v*c, u*b+v*d)

                # For the generalized Heilbronn operator, we through away any
                # symbols for which the (u,v) part of the symbol doesn't have
                # both entries divisible by t.
                if t != 1:
                    if uu%t != 0 or vv%t != 0:
                        continue
                    uu = uu//t
                    vv = vv//t

                # Now coerce each Manin symbol
                #
                #         P[m]*[X^m*Y^(k-2-m), (uu,vv)],    for m=0,...,len(P)
                #
                # into the image space M and add that to z.
                # Note that we coerce in Manin symbols as tuples.
                for m in range(len(P)):
                    x = M((m,uu,vv))
                    z += x*P[m]

            rows.append(z.element())

        A = MS(rows)
        return hom(A, "Heilbronn operator(%s,%s)"%(H,t))

    def _repr_(self):
        r"""
        String representation of this Modular Symbols space.

        EXAMPLES::

            sage: m = ModularSymbols(1,weight=12)
            sage: m # indirect doctest
            Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
        """
        return "Modular Symbols space of dimension %s and weight %s for %s with sign %s and character %s over %s"%(
                self.dimension(), self.weight(), self.group(), self.sign(), self.character()._repr_short_(), self.base_ring())

    def _latex_(self):
        r"""
        Latex representation of this Modular Symbols space.

        EXAMPLES::

            sage: m = ModularSymbols(11,weight=12)
            sage: latex(m) # indirect doctest
            \mathrm{ModSym}_{12}(\Gamma_0(11),\left[1\right];\Bold{Q})

            sage: chi = DirichletGroup(7).0
            sage: m = ModularSymbols(chi)
            sage: latex(m)
            \mathrm{ModSym}_{2}(\Gamma_1(7),\left[\zeta_{6}\right];\Bold{Q}(\zeta_{6}))

        """
        return "\\mathrm{ModSym}_{%s}(%s,%s;%s)"%(self.weight(),
                                                     latex.latex(self.group()),
                                                     latex.latex(list(self.character().values_on_gens())),
                                                     latex.latex(self.base_ring()))

    def _matrix_of_operator_on_modular_symbols(self, codomain, R):
        r"""
        Return the matrix of a modular symbols operator.

        .. note::

           Users will usually instead use the simpler interface
           defined, for example, by ``hecke_matrix()`` (see examples),
           though this function allows one to compute much more
           general operators.

        INPUT:

        -  ``codomain`` - space of modular symbols

        - ``R`` (list) -- a list of lists `[a,b,c,d]` of length 4,
          which we view as elements of `GL_2(`QQ)`.


        OUTPUT:

         -- (matrix) The matrix of the operator

        .. MATH::

                            x \mapsto \sum_{g in R} g.x,


        where `g.x` is the formal linear fractional transformation on modular
        symbols, with respect to the standard basis.

        EXAMPLES::

            sage: M = ModularSymbols(37,2)
            sage: M._matrix_of_operator_on_modular_symbols(M,HeilbronnCremona(3))
            [ 4  0  0  0  0]
            [ 0 -3  1  1  0]
            [ 0  3  0  5 -2]
            [ 0 -3  1 -5  3]
            [ 0  0  2  3 -3]

        """
        rows = []
        for b in self.basis():
            v = formal_sum.FormalSum(0, check=False)
            for c, x in b.modular_symbol_rep():
                for g in R:
                    y = x.apply(g)
                    v += y*c
            w = codomain(v).element()
            rows.append(w)
        M = matrix_space.MatrixSpace(self.base_ring(), len(rows), codomain.degree(), sparse=False)
        return M(rows)

    def _compute_atkin_lehner_matrix(self, d):
        r"""
        Return the matrix of the Atkin-Lehner operator `W_d`.

        INPUT:

        - ``d`` (int) -- an integer that divides the level.

        OUTPUT:

        (matrix) The matrix of the operator `W_d` with respect to
        the standard basis.

        EXAMPLES: An example at level 29::

            sage: M = ModularSymbols((DirichletGroup(29,QQ).0), 2,1); M
            Modular Symbols space of dimension 4 and level 29, weight 2, character [-1], sign 1, over Rational Field
            sage: w = M._compute_atkin_lehner_matrix(29)
            sage: w^2 == 1
            True
            sage: w.fcp()
            (x - 1)^2 * (x + 1)^2

        This doesn't work since the character is not trivial or quadratic::

            sage: M = ModularSymbols((DirichletGroup(13).0), 2,1); M
            Modular Symbols space of dimension 0 and level 13, weight 2, character [zeta12], sign 1, over Cyclotomic Field of order 12 and degree 4
            sage: M._compute_atkin_lehner_matrix(13)
            Traceback (most recent call last):
            ...
            ValueError: Atkin-Lehner W_d only defined when d-primary part of character is trivial or quadratic

        Note that Atkin-Lehner does make sense on `\Gamma_1(13)`,
        but doesn't commute with the Hecke operators::

            sage: M = ModularSymbols(Gamma1(13),2)
            sage: w = M.atkin_lehner_operator(13).matrix()
            sage: t = M.T(2).matrix()
            sage: t*w == w*t
            False
            sage: t * w * ~t * ~w == M.diamond_bracket_matrix(2)
            True
            sage: w^2 == 1
            True

        For `\Gamma_1(N)` levels, when `d` is a proper factor of `N`, the
        square of the operator `W_d` is not a scalar any more::

            sage: M = ModularSymbols(Gamma1(10), 2)
            sage: w = M.atkin_lehner_operator(2).matrix()
            sage: w^2 == M.diamond_bracket_matrix(7)
            True

        In higher weights, the operator is defined, but its eigenvalues are no longer roots of unity::

            sage: M = ModularSymbols(Gamma1(13), 3)
            sage: w = M.atkin_lehner_operator(13).matrix()
            sage: w**2 == -13
            True

        TESTS:

        Check that signed spaces are handled gracefully::

            sage: M = ModularSymbols(Gamma1(18), 3, sign=1)
            sage: M.atkin_lehner_operator(2).matrix().fcp()
            (x^2 + 2)^3 * (x^4 - 2*x^2 + 4)^3
            sage: M.atkin_lehner_operator(9)
            Traceback (most recent call last):
            ...
            ValueError: Atkin-Lehner operator not defined on signed space (use sign=0)

        GammaH spaces work::

            sage: G = GammaH(25, [6])
            sage: ModularSymbols(G, sign=0).atkin_lehner_operator().fcp()
            (x - 1)^5 * (x + 1)^6
        """
        N = self.level()

        chi = self.character()
        if chi is not None:
            dec = [u for u in chi.decomposition() if chi.modulus().divides(d)]
            if not all((u**2).is_trivial() for u in dec):
                raise ValueError("Atkin-Lehner W_d only defined when d-primary part of character is trivial or quadratic")

        if self.sign() != 0:
            # AL operator problematic on signed spaces
            if self.diamond_bracket_matrix(crt(-1, 1, d, N/d)) != 1:
                raise ValueError("Atkin-Lehner operator not defined on signed space (use sign=0)")

        W = self.group().atkin_lehner_matrix(d).list()
        return self._action_on_modular_symbols(W)

    def boundary_map(self):
        r"""
        Return the boundary map to the corresponding space of boundary modular
        symbols.

        EXAMPLES::

            sage: ModularSymbols(20,2).boundary_map()
            Hecke module morphism boundary map defined by the matrix
            [ 1 -1  0  0  0  0]
            [ 0  1 -1  0  0  0]
            [ 0  1  0 -1  0  0]
            [ 0  0  0 -1  1  0]
            [ 0  1  0 -1  0  0]
            [ 0  0  1 -1  0  0]
            [ 0  1  0  0  0 -1]
            Domain: Modular Symbols space of dimension 7 for Gamma_0(20) of weight ...
            Codomain: Space of Boundary Modular Symbols for Congruence Subgroup Gamma0(20) ...
            sage: type(ModularSymbols(20,2).boundary_map())
            <class 'sage.modular.hecke.morphism.HeckeModuleMorphism_matrix'>
        """
        try:
            return self.__boundary_map
        except AttributeError:
            # compute boundary map
            B = self.boundary_space()
            I = [B(b) for b in self.basis()]
            W = matrix_space.MatrixSpace(self.base_ring(), len(I), B.rank(), sparse=True)

            # Note -- the underlying elements have degree the number of distinct
            # cusps known when the element was computed.  This isn't constant,
            # so we pad the elements.
            E = [x.element() for x in I]
            zero = self.base_ring()(0)
            n = int(B.dimension())
            E = sum([ list(x) + [zero]*(n - len(x)) for x in E ], [])

            A = W( E )
            H = cat.Hom(self, B)
            self.__boundary_map = H(A, "boundary map")
            return self.__boundary_map

    def cusps(self):
        r"""
        Return the set of cusps for this modular symbols space.

        EXAMPLES::

            sage: ModularSymbols(20,2).cusps()
            [Infinity, 0, -1/4, 1/5, -1/2, 1/10]
        """
        try:
            return self.__cusps
        except AttributeError:
            f = self.boundary_map()
            B = f.codomain()
            C = B._known_cusps()
            self.__cusps = C
            return C

    def boundary_space(self):
        r"""
        Return the subspace of boundary modular symbols of this modular symbols ambient space.

        EXAMPLES::

            sage: M = ModularSymbols(20, 2)
            sage: B = M.boundary_space(); B
            Space of Boundary Modular Symbols for Congruence Subgroup Gamma0(20) of weight 2 over Rational Field
            sage: M.cusps()
            [Infinity, 0, -1/4, 1/5, -1/2, 1/10]
            sage: M.dimension()
            7
            sage: B.dimension()
            6
        """
        raise NotImplementedError

    def cuspidal_submodule(self):
        """
        The cuspidal submodule of this modular symbols ambient space.

        EXAMPLES::

            sage: M = ModularSymbols(12,2,0,GF(5)) ; M
            Modular Symbols space of dimension 5 for Gamma_0(12) of weight 2 with sign 0 over Finite Field of size 5
            sage: M.cuspidal_submodule()
            Modular Symbols subspace of dimension 0 of Modular Symbols space of dimension 5 for Gamma_0(12) of weight 2 with sign 0 over Finite Field of size 5
            sage: ModularSymbols(1,24,-1).cuspidal_submodule()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 2 for Gamma_0(1) of weight 24 with sign -1 over Rational Field

        The cuspidal submodule of the cuspidal submodule is itself::

            sage: M = ModularSymbols(389)
            sage: S = M.cuspidal_submodule()
            sage: S.cuspidal_submodule() is S
            True
        """
        try:
            return self.__cuspidal_submodule
        except AttributeError:
            try:
                if self.__is_cuspidal:
                    return self
            except AttributeError:
                pass
            S = self.boundary_map().kernel()
            S._set_is_cuspidal(True)
            S._is_full_hecke_module = True
            ## We know the cuspidal subspace is stable, so
            ## if it's one-dimensional, it must be simple
            if S.dimension() == 1:
                S._is_simple = True
            if self.base_ring().characteristic() == 0:
                d = self._cuspidal_submodule_dimension_formula()
                if d is not None:
                    assert d == S.dimension(), "According to dimension formulas the cuspidal subspace of \"%s\" has dimension %s; however, computing it using modular symbols we obtained %s, so there is a bug (please report!)."%(self, d, S.dimension())
            self.__cuspidal_submodule = S
        return self.__cuspidal_submodule

    def _degeneracy_raising_matrix(self, M, t):
        r"""
        Return the matrix of the level-raising degeneracy map from self to M,
        of index t. This is calculated by composing the level-raising matrix
        for `t = 1` with a Hecke operator.

        INPUT:

        - ``M`` (int) -- a space of modular symbols whose level is an integer
          multiple of the level of self

        - ``t`` (int) -- a positive integer dividing the quotient of the two
          levels.

        OUTPUT:

        (matrix) The matrix of the degeneracy map of index `t` from this space
        of level `N` to the space `M` (of level a multiple of `N`). Here `t` is
        a divisor of the quotient.

        EXAMPLES::

            sage: A = ModularSymbols(11, 2); B = ModularSymbols(22, 2)
            sage: A._degeneracy_raising_matrix(B, 1)
            [ 1  0  0  0  0 -1 -1]
            [ 0  1  0 -3  1  1 -1]
            [ 0  1  1 -1 -1  0  0]
            sage: A._degeneracy_raising_matrix(B, 2)
            [ 2  0  0  0  1  0 -1]
            [ 0  0 -1  3 -1 -1  1]
            [ 0 -1 -1  1  0  1 -1]

        Check that :trac:`13198` is fixed::

            sage: M22 = ModularSymbols(Gamma1(22), sign=1)
            sage: M2 = ModularSymbols(Gamma1(2))
            sage: d1 = M2.degeneracy_map(M22,1)
            sage: d2 = M2.degeneracy_map(M22,11)
            sage: M22.hecke_matrix(17).restrict((d1.image() + d2.image()).free_module())
            [18  0]
            [ 0 18]
            sage: S = M22.cuspidal_submodule()
            sage: S.new_submodule().intersection(S.old_submodule()) == S.zero_submodule()
            True
        """
        if t == 1:
            return self._degeneracy_raising_matrix_1(M)
        else:
            # use Hecke operator and t=1 case.
            d1 = self.degeneracy_map(M, 1).matrix()
            T = M.hecke_matrix(t)
            return (~self.base_ring()(t)) * d1 * T

    def _degeneracy_raising_matrix_1(self, M):
        r"""
        Return the matrix of the degeneracy map to the given level
        (which must be a multiple of the level of self).

        .. note::

           Not implemented in the base class, only in the derived classes.

        EXAMPLES::

            sage: M = ModularSymbols(37,4)
            sage: M._degeneracy_raising_matrix_1(ModularSymbols(74, 4))
            20 x 58 dense matrix over Rational Field (use the '.str()' method to see the entries)
        """
        raise NotImplementedError

    def _degeneracy_lowering_matrix(self, M, t):
        r"""
        Return the matrix of the level-lowering degeneracy map from self to M.

        INPUT:

        - ``M`` -- a modular symbols space whose level divides the level of
          self

        - ``t`` (int) -- a positive integer dividing the quotient of the
          levels.

        OUTPUT:

        (matrix) The matrix of the degeneracy map from this space to the space
        `M` of index `t`, where `t` is a divisor of the quotient of the levels
        of self and `M`.

        EXAMPLES::

            sage: M = ModularSymbols(22,2)
            sage: M._degeneracy_lowering_matrix(ModularSymbols(11, 2), 2)
            [ 1  0  0]
            [ 0  1 -1]
            [ 0  0 -1]
            [ 0  1  0]
            [ 0  0  0]
            [-1  0  1]
            [-1  0  0]
        """
        # Use Proposition 2.6.15 in Merel's 1585 paper (or Prop 15 in
        # electronic version of that paper).
        H = heilbronn.HeilbronnMerel(t)
        return self.__heilbronn_operator(M,H,t).matrix()

    def rank(self):
        """
        Return the rank of this modular symbols ambient space.

        OUTPUT:

        (int) The rank of this space of modular symbols.

        EXAMPLES::

            sage: M = ModularSymbols(389)
            sage: M.rank()
            65

        ::

            sage: ModularSymbols(11,sign=0).rank()
            3
            sage: ModularSymbols(100,sign=0).rank()
            31
            sage: ModularSymbols(22,sign=1).rank()
            5
            sage: ModularSymbols(1,12).rank()
            3
            sage: ModularSymbols(3,4).rank()
            2
            sage: ModularSymbols(8,6,sign=-1).rank()
            3
        """
        try:
            return self.__rank
        except AttributeError:
            self.__rank = len(self.manin_basis())
        return self.__rank

    def eisenstein_submodule(self):
        """
        Return the Eisenstein submodule of this space of modular symbols.

        EXAMPLES::

            sage: ModularSymbols(20,2).eisenstein_submodule()
            Modular Symbols subspace of dimension 5 of Modular Symbols space of dimension 7 for Gamma_0(20) of weight 2 with sign 0 over Rational Field
        """
        try:
            return self.__eisenstein_submodule
        except AttributeError:
            self.__eisenstein_submodule = self.cuspidal_submodule().complement()
            return self.__eisenstein_submodule

    def element(self, x):
        """
        Creates and returns an element of self from a modular symbol, if
        possible.

        INPUT:


        -  ``x`` - an object of one of the following types:
           ModularSymbol, ManinSymbol.


        OUTPUT:

        ModularSymbol - a modular symbol with parent self.

        EXAMPLES::

            sage: M = ModularSymbols(11,4,1)
            sage: M.T(3)
            Hecke operator T_3 on Modular Symbols space of dimension 4 for Gamma_0(11) of weight 4 with sign 1 over Rational Field
            sage: M.T(3)(M.0)
            28*[X^2,(0,1)] + 2*[X^2,(1,4)] + 2/3*[X^2,(1,6)] - 8/3*[X^2,(1,9)]
            sage: M.T(3)(M.0).element()
            (28, 2, 2/3, -8/3)
        """
        if isinstance(x, ManinSymbol):
            if not x.parent().weight() == self.weight():
                raise ArithmeticError("incompatible weights: Manin symbol\
                    has weight %s, but modular symbols space has weight %s"%(
                    x.parent().weight(), self.weight()))
            t = self.manin_symbols().index(x.tuple())
            if isinstance(t, tuple):
                i, scalar = t
                v = self.manin_gens_to_basis().row(i) * scalar
            else:
                v = self.manin_gens_to_basis().row(t)
            return self.element_class(self, v)

        elif isinstance(x, element.ModularSymbolsElement):
            M = x.parent()
            if M.ambient_hecke_module() != self:
                # TODO -- sometimes do something more sophisticated here.
                raise TypeError("Modular symbol (%s) does not lie in this space."%x)
            return self(x.element())

        else:
            raise ValueError("Cannot create element of %s from %s."%(x,self))

    def dual_star_involution_matrix(self):
        """
        Return the matrix of the dual star involution, which is induced by
        complex conjugation on the linear dual of modular symbols.

        EXAMPLES::

            sage: ModularSymbols(20,2).dual_star_involution_matrix()
            [1 0 0 0 0 0 0]
            [0 1 0 0 0 0 0]
            [0 0 0 0 1 0 0]
            [0 0 0 1 0 0 0]
            [0 0 1 0 0 0 0]
            [0 0 0 0 0 1 0]
            [0 0 0 0 0 0 1]
        """
        try:
            return self.__dual_star_involution_matrix
        except AttributeError:
            pass
        self.__dual_star_involution_matrix = self.star_involution().matrix().transpose()
        return self.__dual_star_involution_matrix

    def factorization(self):
        r"""
        Return a list of pairs `(S,e)` where `S` is spaces
        of modular symbols and self is isomorphic to the direct sum of the
        `S^e` as a module over the *anemic* Hecke algebra adjoin
        the star involution. The cuspidal `S` are all simple, but
        the Eisenstein factors need not be simple.

        EXAMPLES::

            sage: ModularSymbols(Gamma0(22), 2).factorization()
            (Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field)^2 *
            (Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field)^2 *
            (Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 7 for Gamma_0(22) of weight 2 with sign 0 over Rational Field)

        ::

            sage: ModularSymbols(1,6,0,GF(2)).factorization()
            (Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 for Gamma_0(1) of weight 6 with sign 0 over Finite Field of size 2) *
            (Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 for Gamma_0(1) of weight 6 with sign 0 over Finite Field of size 2)

        ::

            sage: ModularSymbols(18,2).factorization()
            (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 7 for Gamma_0(18) of weight 2 with sign 0 over Rational Field) *
            (Modular Symbols subspace of dimension 5 of Modular Symbols space of dimension 7 for Gamma_0(18) of weight 2 with sign 0 over Rational Field)

        ::

            sage: M = ModularSymbols(DirichletGroup(38,CyclotomicField(3)).0^2,  2, +1); M
            Modular Symbols space of dimension 7 and level 38, weight 2, character [zeta3], sign 1, over Cyclotomic Field of order 3 and degree 2
            sage: M.factorization()                    # long time (about 8 seconds)
            (Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 7 and level 38, weight 2, character [zeta3], sign 1, over Cyclotomic Field of order 3 and degree 2) *
            (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 7 and level 38, weight 2, character [zeta3], sign 1, over Cyclotomic Field of order 3 and degree 2) *
            (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 7 and level 38, weight 2, character [zeta3], sign 1, over Cyclotomic Field of order 3 and degree 2) *
            (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 7 and level 38, weight 2, character [zeta3], sign 1, over Cyclotomic Field of order 3 and degree 2)
        """

##         EXAMPLES::

##             sage: M = ModularSymbols(Gamma0(22), 2); M
##             Modular Symbols space of dimension 7 for Gamma_0(22) of weight 2 with sign 0 over Rational Field
##             sage: M.factorization():
##             ...    print b.dimension(), b.level(), e
##             1 11 2
##             1 11 2
##             1 11 2
##             1 22 1

##         An example with sign 1::

##             sage: M = ModularSymbols(Gamma0(22), 2, sign=1); M
##             Modular Symbols space of dimension 5 for Gamma_0(22) of weight 2 with sign 1 over Rational Field
##             sage: for b, e in M.factorization():
##             ...    print b.dimension(), b.level(), e
##             1 11 2
##             1 11 2
##             1 22 1

##         An example for Gamma1::

##             sage: M = ModularSymbols(Gamma1(26), 2, sign=1); M
##             Modular Symbols space of dimension 33 for Gamma_1(26) of weight 2 with sign 1 over Rational Field
##             sage: for b, e in M.factorization():
##             ...    print b.dimension(), b.level(), e
##             1 13 2
##             1 13 2
##             1 13 2
##             2 13 2
##             2 13 2
##             2 13 2
##             2 13 2
##             2 13 2
##             1 26 1
##             1 26 1
##             1 26 1
##             2 26 1
##             2 26 1

##         An example with level divisible by a square::

##             sage: M = ModularSymbols(Gamma0(2*9),2); M
##             ???
##             sage: for b, e in M.factorization():
##             ...    print b.dimension(), b.level(), e
##             ???
        try:
            return self._factorization
        except AttributeError:
            pass

        try:
            if self._is_simple:
                return [self]
        except AttributeError:
            pass

        D = []

        # Treat the cuspidal and eisenstein parts separately.  The
        # cuspidal part is very straightforward because of
        # Atkin-Lehner-Li theory.  The eisenstein part is trickier,
        # because of E2 and that the new and old Eisenstein subspaces
        # can intersect (e.g., they do for M_2(Gamma_0(6))), even
        # in a way that involves forms other than E_2 (i.e., twists
        # of E2).

        # 1. Cuspidal part -- compute the factors and their multiplicities
        #                     using Atkin-Lehner-Li.

        # 2. Eisenstein part -- just call normal decomposition.

        # In the special case of weight 2 we have to do a bunch of
        # annoying extra work below to deal with the Eisenstein series E_2.

        ## If the characteristic of the base ring is 2,
        ## the star involution is the identity, so we
        ## want to avoid adding each cuspidal submodule
        ## twice.
        if self.base_ring().characteristic() == 2:
            skip_minus = True
        else:
            skip_minus = False

        # The cuspidal part
        # We only run through spaces of level a multiple of the conductor of the character, which
        # we compute below, or set to 1 in case of Gamma_H or Gamma_1
        chi = self.character()
        cond = 1 if chi is None   else   chi.conductor()
        # Now actually run through the divisor levels, taking only the ones with that are
        # a multiple of the conductor.
        for d in reversed(divisors(self.level())):
            if d % cond:
                continue
            n = number_of_divisors(self.level() // d)
            M = self.modular_symbols_of_level(d)
            N = M.new_submodule().cuspidal_submodule().decomposition()
            for A in N:
                if self.sign() == 0:
                    V = A.plus_submodule()
                    V._is_simple = True
                    D.append((V,n))
                    if skip_minus:
                        continue
                    V = A.minus_submodule()
                    V._is_simple = True
                    D.append((V,n))
                else:
                    A._is_simple = True
                    D.append((A,n))
        # The eisenstein part
        for E in self.eisenstein_submodule().decomposition(anemic=True):
            D.append((E,1))

        r = self.dimension()
        s = sum(A.rank() * mult for A, mult in D)
        D = sage.structure.all.Factorization(D, cr=True, sort=False)
        D.sort()
        assert r == s, "bug in factorization --  self has dimension %s, but sum of dimensions of factors is %s\n%s" % (r, s, D)
        self._factorization = D
        return self._factorization

    factor = factorization

    def is_cuspidal(self):
        r"""
        Return True if this space is cuspidal, else False.

        EXAMPLES::

            sage: M = ModularSymbols(20,2)
            sage: M.is_cuspidal()
            False
            sage: S = M.cuspidal_subspace()
            sage: S.is_cuspidal()
            True
            sage: S = M.eisenstein_subspace()
            sage: S.is_cuspidal()
            False
        """
        try:
            return self.__is_cuspidal
        except AttributeError:
            S = self.ambient_hecke_module().cuspidal_submodule()
            self.__is_cuspidal = (S.dimension() == self.dimension())
        return self.__is_cuspidal

    def is_eisenstein(self):
        r"""
        Return True if this space is Eisenstein, else False.

        EXAMPLES::

            sage: M = ModularSymbols(20,2)
            sage: M.is_eisenstein()
            False
            sage: S = M.eisenstein_submodule()
            sage: S.is_eisenstein()
            True
            sage: S = M.cuspidal_subspace()
            sage: S.is_eisenstein()
            False
        """
        try:
            return self.__is_eisenstein
        except AttributeError:
            S = self.ambient_hecke_module().eisenstein_submodule()
            self.__is_eisenstein = self.dimension()==S.dimension()
        return self.__is_eisenstein

    def manin_symbols_basis(self):
        """
        A list of Manin symbols that form a basis for the ambient space
        ``self``.

        OUTPUT:

        -  ``list`` - a list of 2-tuples (if the weight is 2)
           or 3-tuples, which represent the Manin symbols basis for self.


        EXAMPLES::

            sage: m = ModularSymbols(23)
            sage: m.manin_symbols_basis()
            [(1,0), (1,17), (1,19), (1,20), (1,21)]
            sage: m = ModularSymbols(6, weight=4, sign=-1)
            sage: m.manin_symbols_basis()
            [[X^2,(2,1)]]
        """
        s = self.manin_symbols()
        return [s.manin_symbol(i) for i in self.manin_basis()]

    def modular_symbols_of_level(self, G):
        """
        Return a space of modular symbols with the same parameters as
        this space, except the congruence subgroup is changed to `G`.

        INPUT:

        - ``G`` -- either a congruence subgroup or an integer to use
          as the level of such a group.  The given group must either
          contain or be contained in the group defining ``self``.

        TESTS::

            sage: M = ModularSymbols(11)
            sage: M.modular_symbols_of_level(22)
            Modular Symbols space of dimension 7 for Gamma_0(22) of weight 2 with sign 0 over Rational Field
            sage: M.modular_symbols_of_level(Gamma1(22))
            Modular Symbols space of dimension 31 for Gamma_1(22) of weight 2 with sign 0 over Rational Field

            sage: M = ModularSymbols(Gamma1(6))
            sage: M.modular_symbols_of_level(12)
            Modular Symbols space of dimension 9 for Gamma_1(12) of weight 2 with sign 0 over Rational Field
            sage: M.modular_symbols_of_level(Gamma0(3))
            Modular Symbols space of dimension 1 for Gamma_0(3) of weight 2 with sign 0 over Rational Field
            sage: M.modular_symbols_of_level(Gamma0(12))
            Traceback (most recent call last):
            ...
            ValueError: one subgroup must contain the other

            sage: M = ModularSymbols(Gamma1(30),4); M
            Modular Symbols space of dimension 144 for Gamma_1(30) of weight 4 with sign 0 over Rational Field
            sage: M.modular_symbols_of_level(22)
            Traceback (most recent call last):
            ...
            ValueError: one level must divide the other

            sage: M = ModularSymbols(GammaH(15,[7]),6)
            sage: M.modular_symbols_of_level(5)
            Modular Symbols space of dimension 4 for Gamma_0(5) of weight 6 with sign 0 over Rational Field
            sage: M.modular_symbols_of_level(30)
            Modular Symbols space of dimension 60 for Congruence Subgroup Gamma_H(30) with H generated by [7] of weight 6 with sign 0 over Rational Field
            sage: M.modular_symbols_of_level(73)
            Traceback (most recent call last):
            ...
            ValueError: one level must divide the other
        """
        if G in ZZ:
            G = self.group()._new_group_from_level(G)
        elif not (self.group().is_subgroup(G) or G.is_subgroup(self.group())):
            raise ValueError('one subgroup must contain the other')
        return modsym.ModularSymbols(G, self.weight(), self.sign(), self.base_ring())

    def modular_symbols_of_sign(self, sign):
        r"""
        Return a space of modular symbols with the same defining
        properties (weight, level, etc.) as this space except with given
        sign.

        INPUT:

        - ``sign`` (int) -- A sign (`+1`, `-1` or `0`).

        OUTPUT:

        (ModularSymbolsAmbient) A space of modular symbols with the
        same defining properties (weight, level, etc.) as this space
        except with given sign.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma0(11),2,sign=0)
            sage: M
            Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: M.modular_symbols_of_sign(-1)
            Modular Symbols space of dimension 1 for Gamma_0(11) of weight 2 with sign -1 over Rational Field
            sage: M = ModularSymbols(Gamma1(11),2,sign=0)
            sage: M.modular_symbols_of_sign(-1)
            Modular Symbols space of dimension 1 for Gamma_1(11) of weight 2 with sign -1 over Rational Field
        """
        if sign == self.sign():
            return self
        return modsym.ModularSymbols(self.group(), self.weight(), sign=sign, base_ring=self.base_ring())


    def modular_symbols_of_weight(self, k):
        r"""
        Return a space of modular symbols with the same defining
        properties (weight, sign, etc.) as this space except with weight
        `k`.

        INPUT:

        - ``k`` (int) -- A positive integer.

        OUTPUT:

        (ModularSymbolsAmbient) A space of modular symbols with the
        same defining properties (level, sign) as this space
        except with given weight.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(6),2,sign=0)
            sage: M.modular_symbols_of_weight(3)
            Modular Symbols space of dimension 4 for Gamma_1(6) of weight 3 with sign 0 over Rational Field
        """
        if k == self.weight():
            return self
        return modsym.ModularSymbols(self.group(), weight=k, sign=self.sign(), base_ring=self.base_ring())

    def _compute_sign_submodule(self, sign, compute_dual=True):
        r"""
        Return the subspace of self that is fixed under the star
        involution.

        INPUT:


        -  ``sign`` - int (either -1 or +1)

        -  ``compute_dual`` - bool (default: True) also
           compute dual subspace. This are useful for many algorithms.


        OUTPUT:

        A subspace of modular symbols

        EXAMPLES::

            sage: ModularSymbols(1,12,0,GF(5)).minus_submodule() ## indirect doctest
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Finite Field of size 5
        """
        S = self.star_involution().matrix() - self.base_ring()(sign)
        V = S.kernel()
        if compute_dual:
            Vdual = S.transpose().kernel()
            M = self.submodule(V, Vdual, check=False)
        else:
            M = self.submodule(V, check=False)
        M._set_sign(sign)
        return M

    def star_involution(self):
        r"""
        Return the star involution on this modular symbols space.

        OUTPUT:

        (matrix) The matrix of the star involution on this space,
        which is induced by complex conjugation on modular symbols,
        with respect to the standard basis.

        EXAMPLES::

            sage: ModularSymbols(20,2).star_involution()
            Hecke module morphism Star involution on Modular Symbols space of dimension 7 for Gamma_0(20) of weight 2 with sign 0 over Rational Field defined by the matrix
            [1 0 0 0 0 0 0]
            [0 1 0 0 0 0 0]
            [0 0 0 0 1 0 0]
            [0 0 0 1 0 0 0]
            [0 0 1 0 0 0 0]
            [0 0 0 0 0 1 0]
            [0 0 0 0 0 0 1]
            Domain: Modular Symbols space of dimension 7 for Gamma_0(20) of weight ...
            Codomain: Modular Symbols space of dimension 7 for Gamma_0(20) of weight ...
        """
        try:
            return self.__star_involution
        except AttributeError:
            pass
        S = self.__heilbronn_operator(self, [[-1,0, 0,1]], 1)
        S.name("Star involution on %s"%self)
        self.__star_involution = S
        return self.__star_involution

    def _compute_diamond_matrix(self, d):
        r"""
        Return the diamond bracket d operator on this modular symbols space.

        INPUT:

            - `d` -- integer

        OUTPUT:

            - ``matrix`` - the matrix of the diamond bracket operator
              on this space.

        EXAMPLES::

            sage: e = kronecker_character(7)
            sage: M = ModularSymbols(e,2,sign=1)
            sage: D = M.diamond_bracket_operator(5); D
            Diamond bracket operator <5> on Modular Symbols space ...
            sage: D.matrix() # indirect doctest
            [-1  0  0  0]
            [ 0 -1  0  0]
            [ 0  0 -1  0]
            [ 0  0  0 -1]
            sage: [M.diamond_bracket_operator(d).matrix()[0,0] for d in [0..6]]
            [0, 1, 0, 1, 0, -1, 0]
            sage: [e(d) for d in [0..6]]
            [0, 1, 0, 1, 0, -1, 0]

        We test that the sign issue at :trac:`8620` is fixed::

            sage: M = Newforms(Gamma1(13),names = 'a')[0].modular_symbols(sign=0)
            sage: M.diamond_bracket_operator(4).matrix()
            [-1  1  1 -1]
            [-1  0  1  0]
            [ 0  0  0 -1]
            [ 0  0  1 -1]

        We check that the result is correctly normalised for weight > 2::

            sage: ModularSymbols(Gamma1(13), 5).diamond_bracket_operator(6).charpoly().factor()
            (x^2 + 1)^8 * (x^4 - x^2 + 1)^10
        """
        return self.__heilbronn_operator(self, [[d,0, 0,d]], 1).matrix() * d**(2 - self.weight())

    def submodule(self, M, dual_free_module=None, check=True):
        r"""
        Return the submodule with given generators or free module `M`.

        INPUT:


        -  ``M`` - either a submodule of this ambient free module, or
           generators for a submodule;

        - ``dual_free_module`` (bool, default None) -- this may be
           useful to speed up certain calculations; it is the
           corresponding submodule of the ambient dual module;

        - ``check`` (bool, default True) -- if True, check that `M` is
           a submodule, i.e. is invariant under all Hecke operators.

        OUTPUT:

        A subspace of this modular symbol space.

        EXAMPLES::

            sage: M = ModularSymbols(11)
            sage: M.submodule([M.0])
            Traceback (most recent call last):
            ...
            ValueError: The submodule must be invariant under all Hecke operators.
            sage: M.eisenstein_submodule().basis()
            ((1,0) - 1/5*(1,9),)
            sage: M.basis()
            ((1,0), (1,8), (1,9))
            sage: M.submodule([M.0 - 1/5*M.2])
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field

        .. note::

           It would make more sense to only check that `M` is invariant
           under the Hecke operators with index coprime to the level.
           Unfortunately, I do not know a reasonable algorithm for
           determining whether a module is invariant under just the
           anemic Hecke algebra, since I do not know an analogue of
           the Sturm bound for the anemic Hecke algebra. - William
           Stein, 2007-07-27
        """
        if check:
            if not free_module.is_FreeModule(M):
                V = self.free_module()
                if not isinstance(M, (list,tuple)):
                    M = M.gens()
                M = V.span([V(x.element()) for x in M])
        return subspace.ModularSymbolsSubspace(self, M, dual_free_module=dual_free_module, check=check)

    def twisted_winding_element(self, i, eps):
        r"""
        Return the twisted winding element of given degree and character.

        INPUT:

        - ``i`` (int) -- an integer, `0\le i\le k-2` where `k` is the weight.

        - ``eps`` (character) -- a Dirichlet character

        OUTPUT:

        (modular symbol) The so-called 'twisted winding element':

        .. MATH::

                \sum_{a \in (\ZZ/m\ZZ)^\times} \varepsilon(a) * [ i, 0, a/m ].

        .. note::

           This will only work if the base ring of the modular symbol
           space contains the character values.

        EXAMPLES::

            sage: eps = DirichletGroup(5)[2]
            sage: K = eps.base_ring()
            sage: M = ModularSymbols(37,2,0,K)
            sage: M.twisted_winding_element(0,eps)
            2*(1,23) - 2*(1,32) + 2*(1,34)
        """
        if not dirichlet.is_DirichletCharacter(eps):
            raise TypeError("eps must be a Dirichlet character.")
        if (i < 0) or (i > self.weight() - 2):
            raise ValueError("i must be between 0 and k-2.")

        m = eps.modulus()
        return self.sum(eps(a) * self.modular_symbol([i, Cusp(0), Cusp(a / m)])
                        for a in m.coprime_integers(m))

    ######################################################################
    # Z-module of integral modular symbols.
    #######################################################################
    def integral_structure(self, algorithm='default'):
        r"""
        Return the `\ZZ`-structure of this modular symbols
        space, generated by all integral modular symbols.

        INPUT:


        -  ``algorithm`` - string (default: 'default' - choose
           heuristically)

           -  ``'pari'`` - use pari for the HNF computation

           -  ``'padic'`` - use p-adic algorithm (only good for
              dense case)


        ALGORITHM: It suffices to consider lattice generated by the free
        generating symbols `X^iY^{k-2-i}.(u,v)` after quotienting
        out by the `S` (and `I`) relations, since the
        quotient by these relations is the same over any ring.

        EXAMPLES: In weight 2 the rational basis is often integral.

        ::

            sage: M = ModularSymbols(11,2)
            sage: M.integral_structure()
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]

        This is rarely the case in higher weight::

            sage: M = ModularSymbols(6,4)
            sage: M.integral_structure()
            Free module of degree 6 and rank 6 over Integer Ring
            Echelon basis matrix:
            [  1   0   0   0   0   0]
            [  0   1   0   0   0   0]
            [  0   0 1/2 1/2 1/2 1/2]
            [  0   0   0   1   0   0]
            [  0   0   0   0   1   0]
            [  0   0   0   0   0   1]

        Here is an example involving `\Gamma_1(N)`.

        ::

            sage: M = ModularSymbols(Gamma1(5),6)
            sage: M.integral_structure()
            Free module of degree 10 and rank 10 over Integer Ring
            Echelon basis matrix:
            [    1     0     0     0     0     0     0     0     0     0]
            [    0     1     0     0     0     0     0     0     0     0]
            [    0     0  1/96  1/32 23/24     0  1/96     0  7/24 67/96]
            [    0     0     0  1/24 23/24     0     0  1/24   1/4 17/24]
            [    0     0     0     0     1     0     0     0     0     0]
            [    0     0     0     0     0   1/6     0  1/48 23/48   1/3]
            [    0     0     0     0     0     0  1/24  1/24 11/24 11/24]
            [    0     0     0     0     0     0     0  1/16  7/16   1/2]
            [    0     0     0     0     0     0     0     0   1/2   1/2]
            [    0     0     0     0     0     0     0     0     0     1]
        """
        if not self.base_ring() == QQ:
            raise NotImplementedError

        try:
            return self.__integral_structure
        except AttributeError:
            pass

        # The attribute _mod2term is set by self.compute_presentation().
        # It is a list of pairs (n, c), such that the ith element of the list
        # is equivalent to c times the n-th basis Manin symbol.
        G = set([i for i, _ in self._mod2term])

        # Now G is a set of integer i such that these integers gives
        # indices of Manin symbols that together generate the integral
        # structure.  We next obtain the corresponding list of elements
        # by passing to the quotient by the remaining relations
        # via the _manin_gens_to_basis attribute.

        # Next we take each element of X, which gives a linear combination
        # of the basis of the underlying vector space of self, and compute
        # the Z-module they span.

        G = sorted(G)
        B = self._manin_gens_to_basis.matrix_from_rows(list(G)).dense_matrix()
        B, d = B._clear_denom()
        if algorithm == 'default':
            # pari is much better in the weight 2 case when the input
            # matrix is extremely sparse; the p-adic algorithm is
            # terrible in the sparse case.
            if self.weight() == 2:
                algorithm = 'pari'
            else:
                algorithm = 'padic'
        if algorithm == 'pari':
            B = B.echelon_form(algorithm='pari', include_zero_rows=False)
        elif algorithm == 'padic':
            B = B.echelon_form(algorithm='padic', include_zero_rows=False)
        else:
            raise ValueError("unknown algorithm '%s'"%algorithm)
        W = B.row_module()
        if d != 1:
            W = W.scale(1/d)
        self.__integral_structure = W
        assert W.rank() == self.rank(), "there is a bug in computing integral structure on modular symbols"
        return self.__integral_structure

    ######################################################################
    # Eigenvalues
    #######################################################################
    def compact_newform_eigenvalues(self, v, names='alpha'):
        r"""
        Return compact systems of eigenvalues for each Galois conjugacy
        class of cuspidal newforms in this ambient space.

        INPUT:


        -  ``v`` - list of positive integers


        OUTPUT:


        -  ``list`` - of pairs (E, x), where E\*x is a vector
           with entries the eigenvalues `a_n` for
           `n \in v`.


        EXAMPLES::

            sage: M = ModularSymbols(43,2,1)
            sage: X = M.compact_newform_eigenvalues(prime_range(10))
            sage: X[0][0] * X[0][1]
            (-2, -2, -4, 0)
            sage: X[1][0] * X[1][1]
            (alpha1, -alpha1, -alpha1 + 2, alpha1 - 2)

        ::

            sage: M = ModularSymbols(DirichletGroup(24,QQ).1,2,sign=1)
            sage: M.compact_newform_eigenvalues(prime_range(10),'a')
            [(
            [-1/2 -1/2]
            [ 1/2 -1/2]
            [  -1    1]
            [  -2    0], (1, -2*a0 - 1)
            )]
            sage: a = M.compact_newform_eigenvalues([1..10],'a')[0]
            sage: a[0]*a[1]
            (1, a0, a0 + 1, -2*a0 - 2, -2*a0 - 2, -a0 - 2, -2, 2*a0 + 4, -1, 2*a0 + 4)
            sage: M = ModularSymbols(DirichletGroup(13).0^2,2,sign=1)
            sage: M.compact_newform_eigenvalues(prime_range(10),'a')
            [(
            [  -zeta6 - 1]
            [ 2*zeta6 - 2]
            [-2*zeta6 + 1]
            [           0], (1)
            )]
            sage: a = M.compact_newform_eigenvalues([1..10],'a')[0]
            sage: a[0]*a[1]
            (1, -zeta6 - 1, 2*zeta6 - 2, zeta6, -2*zeta6 + 1, -2*zeta6 + 4, 0, 2*zeta6 - 1, -zeta6, 3*zeta6 - 3)
        """
        if self.sign() == 0:
            raise ValueError("sign must be nonzero")
        v = list(v)

        # Get decomposition of this space
        D = self.cuspidal_submodule().new_subspace().decomposition()
        for A in D:
            # since sign is zero and we're on the new cuspidal subspace
            # each factor is definitely simple.
            A._is_simple = True
        B = [A.dual_free_module().basis_matrix().transpose() for A in D]

        # Normalize the names strings.
        names = ['%s%s'%(names,i) for i in range(len(B))]

        # Find an integer i such that the i-th columns of the basis for the
        # dual modules corresponding to the factors in D are all nonzero.
        nz = None
        for i in range(self.dimension()):
            # Decide if this i works, i.e., ith row of every element of B is nonzero.
            bad = False
            for C in B:
                if C.row(i) == 0:
                    # i is bad.
                    bad = True
                    continue
            if bad:
                continue
            # It turns out that i is not bad.
            nz = i
            break

        if nz is not None:
            R = self.hecke_images(nz, v)
            return [(R * m, D[i].dual_eigenvector(names=names[i],
                                                  lift=False, nz=nz))
                    for i, m in enumerate(B)]

        # No single i works, so we do something less uniform.
        ans = []
        cache = {}
        for i in range(len(D)):
            nz = D[i]._eigen_nonzero()
            if nz in cache:
                R = cache[nz]
            else:
                R = self.hecke_images(nz, v)
                cache[nz] = R
            ans.append((R * B[i], D[i].dual_eigenvector(names=names[i],
                                                        lift=False, nz=nz)))
        return ans

    def __pari__(self):
        """
        Return a PARI object corresponding to ``self``.

        TESTS::

            sage: ModularSymbols(Gamma1(5), 2).__pari__()
            Traceback (most recent call last):
            ...
            NotImplementedError: PARI modular symbols are only implemented for Gamma0(n)
        """
        raise NotImplementedError('PARI modular symbols are only implemented for Gamma0(n)')

    def _pari_pairing(self):
        r"""
        Return the matrix of the canonical pairing between this space and
        the corresponding space of modular symbols in PARI.

        Let `M` be an ambient space of modular symbols over a field `K`
        of characteristic 0.  The corresponding space `P` in PARI is
        canonically dual to `M`, giving rise to a `K`-bilinear map

        .. MATH::

            E\colon M \times P \to K.

        OUTPUT: The matrix of the bilinear map `E`.

        This is currently only implemented for spaces of modular
        symbols of trivial character.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma0(5), 6)
            sage: P = M.__pari__()
            sage: E = M._pari_pairing(); E
            [  0  -1   0   0]
            [  0   0   8 -27]
            [  8   0   0  13]
            [ 24   0   8  37]

        The duality is compatible with (for example) Hecke operators
        and the star involution::

            sage: M.hecke_matrix(5) * E == E * P.mshecke(5)
            True
            sage: M.star_involution().matrix() * E == E * P.msstar()
            True
        """
        if self.weight() == 2:
            return self._pari_tensor().inverse()
        from sage.matrix.constructor import matrix
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        P = self.__pari__()
        I = matrix.identity(self.rank()).__pari__()
        m = Integer(self.weight() - 2)
        R = PolynomialRing(QQ, 'x')
        x = R.gen()
        def ev(s):
            # The Manin symbol s = X^i (c, d) corresponds to the
            # modular symbol (dX - bY)^i (-cX + aY)^(m - i) {b/d, a/c}.
            # The code below computes the canonical pairing of this
            # modular symbol with the distinguished basis of P.
            i = s.i
            a, b, c, d = s.lift_to_sl2z()
            e = [R(p) for p in P.mseval(I, matrix(2, 2, [b, a, d, c]))]
            g = (d*x - b)**i * (-c*x + a)**(m - i)
            return [sum(f[j] * g[m - j] / m.binomial(j) for j in range(m + 1))
                    for f in e]
        return matrix([ev(s) for s in self.manin_symbols_basis()])

    def _pari_tensor(self):
        r"""
        Return a matrix expressing the duality between this space and the
        corresponding space of modular symbols in PARI.

        Let `M` be an ambient space of modular symbols over a field `K`
        of characteristic 0.  The corresponding space `P` in PARI is
        canonically dual to `M`, giving rise to an element

        .. MATH::

            T \in P \otimes_K M.

        OUTPUT: The matrix of the element `T \in P \otimes_K M`.
        This is the inverse of the matrix returned by
        :meth:`_pari_pairing`.

        This is currently only implemented for spaces of modular
        symbols of trivial character.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma0(37), 2)
            sage: P = M.__pari__()
            sage: T = M._pari_tensor(); T
            [ 1  0  0  0  0]
            [ 0  0 -1  0  0]
            [ 0  0  1 -1  0]
            [ 0 -1  0 -1  1]
            [ 0  0  0  1 -1]

        The duality is compatible with (for example) Hecke operators
        and the star involution::

            sage: T * M.hecke_matrix(3) == P.mshecke(3) * T
            True
            sage: T * M.star_involution().matrix() == P.msstar() * T
            True
        """
        if self.weight() != 2:
            return self._pari_pairing().inverse()
        from sage.matrix.constructor import matrix
        gens = self.__pari__().mspathgens()[0][:self.rank()].sage()
        # gens is a basis for the space of modular symbols of weight 2
        # (in the sense of Sage), and the distinguished basis of the
        # corresponding PARI space of modular symbols is dual to this.
        return matrix([self.modular_symbol(g).element() for g in gens])


class ModularSymbolsAmbient_wtk_g0(ModularSymbolsAmbient):
    r"""
    Modular symbols for `\Gamma_0(N)` of integer weight
    `k > 2` over the field `F`.

    For weight `2`, it is faster to use ``ModularSymbols_wt2_g0``.

    INPUT:


    -  ``N`` - int, the level

    -  ``k`` - integer weight = 2.

    -  ``sign`` - int, either -1, 0, or 1

    -  ``F`` - field


    EXAMPLES::

        sage: ModularSymbols(1,12)
        Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
        sage: ModularSymbols(1,12, sign=1).dimension()
        2
        sage: ModularSymbols(15,4, sign=-1).dimension()
        4
        sage: ModularSymbols(6,6).dimension()
        10
        sage: ModularSymbols(36,4).dimension()
        36
    """
    def __init__(self, N, k, sign, F, custom_init=None, category=None):
        r"""
        Initialize a space of modular symbols of weight `k` for
        `\Gamma_0(N)`, over `\QQ`.

        For weight `2`, it is faster to use
        ``ModularSymbols_wt2_g0``.

        INPUT:


        -  ``N`` - int, the level

        -  ``k`` - integer weight = 2.

        -  ``sign`` - int, either -1, 0, or 1

        -  ``F`` - field


        EXAMPLES::

            sage: ModularSymbols(1,12)
            Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
            sage: ModularSymbols(1,12, sign=1).dimension()
            2
            sage: ModularSymbols(15,4, sign=-1).dimension()
            4
            sage: ModularSymbols(6,6).dimension()
            10
            sage: ModularSymbols(36,4).dimension()
            36
        """
        N = int(N)
        k = int(k)
        sign = int(sign)
        if sign not in [-1, 0, 1]:
            raise TypeError("sign must be an int in [-1,0,1]")

        ModularSymbolsAmbient.__init__(self, weight=k, group=arithgroup.Gamma0(N),
                                       sign=sign, base_ring=F,
                                       custom_init=custom_init, category=category)


    def _dimension_formula(self):
        r"""
        Return the dimension of this space using the formula.

        EXAMPLES::

            sage: M = ModularSymbols(37,6)
            sage: M.dimension()
            32
            sage: M._dimension_formula()
            32
        """
        if self.base_ring().characteristic() == 0:
            k, sign = self.weight(), self.sign()
            if sign != 0:
                return None
            if k % 2:
                return 0
            elif k > 2:
                return 2*self.group().dimension_cusp_forms(k) + self.group().ncusps()
            else:
                return 2*self.group().dimension_cusp_forms(k) + self.group().ncusps() - 1
        else:
            raise NotImplementedError

    def _repr_(self):
        r"""
        Return the string representation of this Modular Symbols space.

        EXAMPLES::

            sage: M = ModularSymbols(37,6)
            sage: M # indirect doctest
            Modular Symbols space of dimension 32 for Gamma_0(37) of weight 6 with sign 0 over Rational Field
        """
        return ("Modular Symbols space of dimension %s for Gamma_0(%s) of weight %s with sign %s " + \
                "over %s")%(self.dimension(), self.level(),self.weight(), self.sign(),
                            self.base_ring())

    def _cuspidal_submodule_dimension_formula(self):
        r"""
        Return the dimension of the cuspidal subspace, using the formula.

        EXAMPLES::

            sage: M = ModularSymbols(37,4)
            sage: M.cuspidal_subspace().dimension()
            18
            sage: M._cuspidal_submodule_dimension_formula()
            18
        """
        if self.base_ring().characteristic() == 0:
            k, sign = self.weight(), self.sign()
            if sign == 0:
                m = 2
            else:
                m = 1
            return m * self.group().dimension_cusp_forms(k)
        else:
            raise NotImplementedError


    def _degeneracy_raising_matrix_1(self, M):
        r"""
        Return the matrix of the degeneracy map (with t = 1) to level
        `N`, where `N` is a multiple of the level.

        INPUT:

        - ``M`` -- A space of Gamma0 modular symbols of the same weight as
          self, with level an integer multiple of the level of self.

        OUTPUT:

        (matrix) The matrix of the degeneracy raising map to `M`.

        EXAMPLES::

            sage: M = ModularSymbols(37,4)
            sage: M._degeneracy_raising_matrix_1(ModularSymbols(74, 4))
            20 x 58 dense matrix over Rational Field (use the '.str()' method to see the entries)
            sage: M.dimension()
            20
            sage: ModularSymbols(74,4).dimension()
            58
        """
        level = int(M.level())
        N = self.level()

        # 1. Find coset representatives H for Gamma_0(M.level()) \ Gamma_0(self.level())
        #    (need to be careful in some small levels, cf. #13198)

        if arithgroup.is_Gamma0(M.group()):
            H = arithgroup.degeneracy_coset_representatives_gamma0(level, N, 1)
        elif arithgroup.is_Gamma1(M.group()):
            H = arithgroup.degeneracy_coset_representatives_gamma1(level, N, 1)
        else:
            raise NotImplementedError("Degeneracy raising maps not implemented for GammaH levels")

        # 2. The map is
        #        [P,pi(g)] |--> sum_{h in H} [P, pi(h*g)]
        #
        MS = matrix_space.MatrixSpace(self.base_ring(), self.dimension(), M.dimension())
        if self.dimension() == 0 or M.dimension() == 0:
            return MS(0)
        rows = []
        B = self.manin_basis()
        syms = self.manin_symbols()
        H = [M2Z(h) for h in H]
        for n in B:
            z = M(0)
            s = syms.manin_symbol(n)
            g = M2Z(list(s.lift_to_sl2z(N)))
            i = s.i
            # We apply each matrix in H according to the above formula
            for h in H:
                hg = h*g
                z += M((i, hg[1,0], hg[1,1]))
            rows.append(z.element())

        A = MS(rows)
        return A


    def _cuspidal_new_submodule_dimension_formula(self):
        r"""
        Return the dimension of the new cuspidal subspace, via the formula.

        EXAMPLES::

            sage: M = ModularSymbols(100,2)
            sage: M._cuspidal_new_submodule_dimension_formula()
            2
            sage: M.cuspidal_subspace().new_subspace().dimension()
            2
        """
        if self.base_ring().characteristic() == 0:
            k, sign = self.weight(), self.sign()
            if sign == 0:
                m = 2
            else:
                m = 1
            return m * self.group().dimension_new_cusp_forms(k)
        else:
            raise NotImplementedError

    def boundary_space(self):
        r"""
        Return the space of boundary modular symbols for this space.

        EXAMPLES::

            sage: M = ModularSymbols(100,2)
            sage: M.boundary_space()
            Space of Boundary Modular Symbols for Congruence Subgroup Gamma0(100) of weight 2 over Rational Field
        """
        try:
            return self.__boundary_space
        except AttributeError:
            pass
        self.__boundary_space = boundary.BoundarySpace_wtk_g0(
            self.level(), self.weight(), self.sign(), self.base_ring())
        return self.__boundary_space

    def manin_symbols(self):
        r"""
        Return the Manin symbol list of this modular symbol space.

        EXAMPLES::

            sage: M = ModularSymbols(100,4)
            sage: M.manin_symbols()
            Manin Symbol List of weight 4 for Gamma0(100)
            sage: len(M.manin_symbols())
            540
        """
        try:
            return self.__manin_symbols
        except AttributeError:
            self.__manin_symbols = ManinSymbolList_gamma0(
                level=self.level(), weight=self.weight())
        return self.__manin_symbols

    def _hecke_images(self, i, v):
        """
        Return matrix whose rows are the images of the `i`-th
        standard basis vector under the Hecke operators `T_p` for
        all integers in `v`.

        INPUT:


        -  ``i`` - nonnegative integer

        -  ``v`` - a list of positive integer


        OUTPUT:

        -  ``matrix`` - whose rows are the Hecke images

        EXAMPLES::

            sage: M = ModularSymbols(11,4,1)
            sage: mat = M._hecke_images(0,[1,2,3,4])
            sage: M.T(1)(M.0).element() == mat[0]
            True
            sage: M.T(2)(M.0).element() == mat[1]
            True
            sage: M.T(3)(M.0).element() == mat[2]
            True
            sage: M.T(4)(M.0).element() == mat[3]
            True

            sage: M = ModularSymbols(12,4)
            sage: mat = M._hecke_images(0,[1,2,3,4])
            sage: M.T(1)(M.0).element() == mat[0]
            True
            sage: M.T(2)(M.0).element() == mat[1]
            True
            sage: M.T(3)(M.0).element() == mat[2]
            True
            sage: M.T(4)(M.0).element() == mat[3]
            True
        """
        # Find basis vector for ambient space such that it is not in
        # the kernel of the dual space corresponding to self.
        c = self.manin_generators()[self.manin_basis()[i]]
        N = self.level()
        return heilbronn.hecke_images_gamma0_weight_k(c.u,c.v, c.i, N, self.weight(),
                                                      v, self.manin_gens_to_basis())

    @cached_method
    def __pari__(self):
        """
        Return a PARI object corresponding to ``self``.

        EXAMPLES::

            sage: ModularSymbols(Gamma0(1), 2).__pari__()
            [[[[Vecsmall([0, 1])], [0], 1, [Vecsmall([]), Vecsmall([])],
               Vecsmall([1]), Vecsmall([]), Vecsmall([])],
              0, 0, 0, Vecsmall([1]), 0, 0, [[1, 1; [0, 1; -1, 0], 1]],
              [[1, 1; [0, -1; 1, -1], 1; [-1, 1; -1, 0], 1]], 0,
              Vecsmall([0, 0, 1, 1, 2]), [[Vecsmall([1, 0]), Vecsmall([0, 1])]],
              0, 0, 0, [Vecsmall([1, 1]), [Vecsmall([0, 1]), 0], [Vecsmall([1, 0])]]],
             [0, [;], [[;], [;], 1, Vecsmall([])]],
             [[], Vecsmall([2, 0]), Vecsmall([0, 0]), 0, [[;], [;], 1, Vecsmall([])]]]

        .. NOTE::

            Spaces of modular symbols as implemented in PARI are
            canonically dual to those implemented in Sage.  See
            :meth:`ModularSymbolsAmbient._pari_pairing` and
            :meth:`ModularSymbolsAmbient._pari_tensor` for how to use
            this duality.
        """
        return self.level().__pari__().msinit(self.weight(), self.sign())


class ModularSymbolsAmbient_wt2_g0(ModularSymbolsAmbient_wtk_g0):
    r"""
    Modular symbols for `\Gamma_0(N)` of integer weight `2` over the field
    `F`.

    INPUT:

    -  ``N`` - int, the level

    -  ``sign`` - int, either -1, 0, or 1


    OUTPUT:

    The space of modular symbols of weight `2`, trivial character, level
    `N` and given sign.

    EXAMPLES::

        sage: ModularSymbols(Gamma0(12),2)
        Modular Symbols space of dimension 5 for Gamma_0(12) of weight 2 with sign 0 over Rational Field
    """
    def __init__(self, N, sign, F, custom_init=None, category=None):
        """
        Initialize a space of modular symbols. INPUT:

        INPUT:

        -  ``N`` - int, the level

        -  ``sign`` - int, either -1, 0, or 1


        OUTPUT:

        The space of modular symbols of weight 2, trivial character,
        level N and given sign.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma0(12),2)
        """
        ModularSymbolsAmbient_wtk_g0.__init__(self, N=N, k=2, sign=sign, F=F,
                                              custom_init=custom_init, category=category)

    def _dimension_formula(self):
        r"""
        Return the dimension of this space using the formula.

        EXAMPLES::

            sage: M = ModularSymbols(37,6)
            sage: M.dimension()
            32
            sage: M._dimension_formula()
            32
        """
        if self.base_ring().characteristic() == 0:
            if self.sign() != 0:
                return None
            return 2*self.group().dimension_cusp_forms(2) + self.group().ncusps() - 1
        else:
            raise NotImplementedError

    def _cuspidal_submodule_dimension_formula(self):
        r"""
        Return the dimension of the cuspidal subspace, using the formula.

        EXAMPLES::

            sage: M = ModularSymbols(37,4)
            sage: M.cuspidal_subspace().dimension()
            18
            sage: M._cuspidal_submodule_dimension_formula()
            18
        """
        if self.base_ring().characteristic() == 0:
            if self.sign() == 0:
                m = 2
            else:
                m = 1
            return m * self.group().dimension_cusp_forms(2)
        else:
            raise NotImplementedError

    def _cuspidal_new_submodule_dimension_formula(self):
        r"""
        Return the dimension of the new cuspidal subspace, via the formula.

        EXAMPLES::

            sage: M = ModularSymbols(100,2)
            sage: M._cuspidal_new_submodule_dimension_formula()
            2
            sage: M.cuspidal_subspace().new_subspace().dimension()
            2
        """
        if self.base_ring().characteristic() == 0:
            if self.sign() == 0:
                m = 2
            else:
                m = 1
            return m * self.group().dimension_new_cusp_forms(2)
        else:
            raise NotImplementedError


    def _compute_hecke_matrix_prime(self, p, rows=None):
        r"""
        Compute and return the matrix of the `p`-th Hecke operator.

        EXAMPLES::

            sage: m = ModularSymbols(37,2)
            sage: m._compute_hecke_matrix_prime(2).charpoly('x')
            x^5 + x^4 - 8*x^3 - 12*x^2
        """
        # note -- p doesn't have to be prime.
        if isinstance(rows, list):
            rows = tuple(rows)
        try:
            return self._hecke_matrices[(p, rows)]
        except AttributeError:
            self._hecke_matrices = {}
        except KeyError:
            pass
        tm = verbose("Computing Hecke operator T_%s" % p)

        H = heilbronn.HeilbronnCremona(p)
        # H = heilbronn.HeilbronnMerel(p)
        B = self.manin_basis()
        if rows is not None:
            B = [B[i] for i in rows]

        N = self.level()
        P1 = self.p1list()
        mod2term = self._mod2term
        R = self.manin_gens_to_basis()
        W = R.new_matrix(nrows=len(B), ncols = R.nrows())  # the 0 with given number of rows and cols.
        j = 0
        tm = verbose("Matrix non-reduced", tm)
        for i in B:
            # The following step is where most of the time is spent.
            c,d = P1[i]
            v = H.apply(c,d, N)

            # v is now a list of pairs ((c,d),m), where m is the
            # number of times that (c,d) appears in the image of x
            # under the matrices in H.  Also, the pairs (c,d) are
            # normalized.
            # Let ind(c,d) denote the index of the normalized pair
            # (c,d) in the fixed ordered list of elements of
            # P1(Z/NZ).  Then the list of pairs (ind(c,d), m)
            # obtained from the above list defines a sparse vector
            # s, and the image of x under T_p is the product
            # of s with the matrix R defined above.
            for z, m in v:
                k = P1.index_of_normalized_pair(z[0],z[1])
                if k != -1:
                    f, s = mod2term[k]
                    if s != 0:
                        W[j,f] = W[j,f] + s*m
            j += 1
        tm = verbose("done making non-reduced matrix",tm)
        verbose("start matrix-matrix (%s x %s) times (%s x %s) multiply to get Tp"%(W.nrows(), W.ncols(),
                                                                                         R.nrows(), R.ncols()))
        if hasattr(W, '_matrix_times_matrix_dense'):
            Tp = W._matrix_times_matrix_dense(R)
            verbose("done matrix multiply and computing Hecke operator",tm)
        else:
            Tp = W * R
            tm = verbose("done multiplying",tm)
            Tp = Tp.dense_matrix()
            verbose("done making Hecke operator dense", tm)
        if rows is None:
            self._hecke_matrices[(p,rows)] = Tp
        return Tp

    def boundary_space(self):
        r"""
        Return the space of boundary modular symbols for this space.

        EXAMPLES::

            sage: M = ModularSymbols(100,2)
            sage: M.boundary_space()
            Space of Boundary Modular Symbols for Congruence Subgroup Gamma0(100) of weight 2 over Rational Field
        """
        try:
            return self.__boundary_space
        except AttributeError:
            pass
        self.__boundary_space = boundary.BoundarySpace_wtk_g0(
            self.level(), self.weight(), self.sign(), self.base_ring())
        return self.__boundary_space

    def _hecke_image_of_ith_basis_vector(self, n, i):
        """
        Return `T_n(e_i)`, where `e_i` is the
        `i` th basis vector of this ambient space.

        INPUT:

        - ``n`` -- an integer which should be prime.

        OUTPUT:

        - ``modular symbol`` -- element of this ambient space

        EXAMPLES::

            sage: M = ModularSymbols(43,2,1)
            sage: M._hecke_image_of_ith_basis_vector(2, 0)
            3*(1,0) - 2*(1,33)
            sage: M.hecke_operator(2)(M.0)
            3*(1,0) - 2*(1,33)
            sage: M._hecke_image_of_ith_basis_vector(6, 1)
            4*(1,31) - 3*(1,33) + 3*(1,39)
            sage: M.hecke_operator(6)(M.1)
            4*(1,31) - 3*(1,33) + 3*(1,39)
        """
        c = self.manin_generators()[self.manin_basis()[i]]
        N = self.level()
        I = heilbronn.hecke_images_gamma0_weight2(c.u,c.v,N,[n], self.manin_gens_to_basis())
        return self(I[0])

    def _hecke_images(self, i, v):
        """
        Return images of the `i`-th standard basis vector under the
        Hecke operators `T_p` for all integers in `v`.

        INPUT:


        -  ``i`` - nonnegative integer

        -  ``v`` - a list of positive integer


        OUTPUT:


        -  ``matrix`` - whose rows are the Hecke images


        EXAMPLES::

            sage: M = ModularSymbols(46,2,-1)
            sage: mat = M._hecke_images(1,[3,4,5,6])
            sage: v = M.basis()[1]
            sage: M.T(3)(v).element() == mat[0]
            True
            sage: M.T(4)(v).element() == mat[1]
            True
            sage: M.T(5)(v).element() == mat[2]
            True
            sage: M.T(6)(v).element() == mat[3]
            True
        """
        # Find basis vector for ambient space such that it is not in
        # the kernel of the dual space corresponding to self.
        c = self.manin_generators()[self.manin_basis()[i]]
        N = self.level()
        return heilbronn.hecke_images_gamma0_weight2(c.u,c.v,N, v, self.manin_gens_to_basis())


class ModularSymbolsAmbient_wtk_g1(ModularSymbolsAmbient):
    r"""
    INPUT:


    -  ``level`` - int, the level

    -  ``weight`` - int, the weight = 2

    -  ``sign`` - int, either -1, 0, or 1

    -  ``F`` - field


    EXAMPLES::

        sage: ModularSymbols(Gamma1(17),2)
        Modular Symbols space of dimension 25 for Gamma_1(17) of weight 2 with sign 0 over Rational Field
        sage: [ModularSymbols(Gamma1(7),k).dimension() for k in [2,3,4,5]]
        [5, 8, 12, 16]

    ::

        sage: ModularSymbols(Gamma1(7),3)
        Modular Symbols space of dimension 8 for Gamma_1(7) of weight 3 with sign 0 over Rational Field
        """

    def __init__(self, level, weight, sign, F, custom_init=None, category=None):
        r"""
        Initialize a space of modular symbols for Gamma1(N).

        INPUT:


        -  ``level`` - int, the level

        -  ``weight`` - int, the weight = 2

        -  ``sign`` - int, either -1, 0, or 1

        -  ``F`` - field


        EXAMPLES::

            sage: ModularSymbols(Gamma1(17),2)
            Modular Symbols space of dimension 25 for Gamma_1(17) of weight 2 with sign 0 over Rational Field
            sage: [ModularSymbols(Gamma1(7),k).dimension() for k in [2,3,4,5]]
            [5, 8, 12, 16]

        ::

            sage: M = ModularSymbols(Gamma1(7),3)
        """
        ModularSymbolsAmbient.__init__(self,
                weight=weight,
                group=arithgroup.Gamma1(level),
                sign=sign,
                base_ring=F,
                custom_init=custom_init,
                category=category)

    def _dimension_formula(self):
        r"""
        Return the dimension of this space using the formula.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(7),6)
            sage: M.dimension()
            20
            sage: M._dimension_formula()
            20
        """
        if self.base_ring().characteristic() != 0:
            raise NotImplementedError
        level, weight, sign = self.level(), self.weight(), self.sign()
        if sign != 0:
            return None
        d = 2*self.group().dimension_cusp_forms(weight) + self.group().ncusps()
        if level == 1 and weight % 2:
            return 0
        if weight == 2:
            return d - 1
        if weight % 2 == 0:
            return d

        # TODO: I don't know a formula for dim ModSym_k(Gamma_1(N)) for odd k!!!

        return None

    def _repr_(self):
        r"""
        Return a string representation of this space.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(7),3)
            sage: M # indirect doctest
            Modular Symbols space of dimension 8 for Gamma_1(7) of weight 3 with sign 0 over Rational Field
        """
        return ("Modular Symbols space of dimension %s for Gamma_1(%s) of weight %s with sign %s over %s"
                % (self.dimension(), self.level(), self.weight(), self.sign(), self.base_ring()))

    def _cuspidal_submodule_dimension_formula(self):
        r"""
        Return the dimension of the cuspidal subspace, using the formula.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(11),4)
            sage: M.cuspidal_subspace().dimension()
            20
            sage: M._cuspidal_submodule_dimension_formula()
            20
        """
        if self.sign() == 0:
            m = 2
        else:
            m = 1
        return m * self.group().dimension_cusp_forms(self.weight())

    def _cuspidal_new_submodule_dimension_formula(self):
        r"""
        Return the dimension of the new cuspidal subspace, via the formula.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(22),2)
            sage: M._cuspidal_new_submodule_dimension_formula()
            8
            sage: M.cuspidal_subspace().new_subspace().dimension()
            8
        """
        if self.sign() == 0:
            m = 2
        else:
            m = 1
        return m * self.group().dimension_new_cusp_forms(self.weight())


    def _compute_hecke_matrix_prime_power(self, p, r):
        r"""
        Compute and return the matrix of the Hecke operator `T(p^r)`.

        EXAMPLES::

            sage: m = ModularSymbols(Gamma1(11),2)
            sage: m._compute_hecke_matrix_prime_power(3,4).charpoly('x')
            x^11 - 291*x^10 + 30555*x^9 - 1636145*x^8 + 59637480*x^7 + 1983040928*x^6 - 401988683888*x^5 - 14142158875680*x^4 + 3243232720819520*x^3 - 103658398669404480*x^2 + 197645665452381696*x - 97215957397309696
        """
        return self._compute_hecke_matrix_prime(p**r)

    def _degeneracy_raising_matrix_1(self, M):
        r"""
        Return the matrix of the degeneracy raising map to `M`.

        INPUT:

        - ``M`` -- an ambient space of Gamma1 modular symbols, of level a
          multiple of the level of self

        OUTPUT:

        (matrix) The matrix of the degeneracy raising matrix to the higher level.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(7),3)
            sage: N = ModularSymbols(Gamma1(21), 3)
            sage: M._degeneracy_raising_matrix_1(N)
            8 x 64 dense matrix over Rational Field (use the '.str()' method to see the entries)
            sage: M.dimension()
            8
            sage: N.dimension()
            64
        """
        N = self.level()

        # 1. Find coset representatives H for Gamma_1(M.level()) \ Gamma_1(self.level())
        H = arithgroup.degeneracy_coset_representatives_gamma1(M.level(), N, 1)
        # 2. The map is
        #        [P,pi(g)] |--> sum_{h in H} [P, pi(h*g)]
        #
        MS = matrix_space.MatrixSpace(self.base_ring(), self.dimension(), M.dimension())
        if self.dimension() == 0 or M.dimension() == 0:
            return MS(0)
        rows = []
        B = self.manin_basis()
        syms = self.manin_symbols()
        G = matrix_space.MatrixSpace(ZZ, 2)
        H = [G(h) for h in H]
        for n in B:
            z = M(0)
            s = syms.manin_symbol(n)
            g = G(list(s.lift_to_sl2z(N)))
            i = s.i
            # We apply each matrix in H according to the above formula
            for h in H:
                hg = h*g
                z += M((i, hg[1,0], hg[1,1]))
            rows.append(z.element())

        A = MS(rows)
        return A

    def boundary_space(self):
        r"""
        Return the space of boundary modular symbols for this space.

        EXAMPLES::

            sage: M = ModularSymbols(100,2)
            sage: M.boundary_space()
            Space of Boundary Modular Symbols for Congruence Subgroup Gamma0(100) of weight 2 over Rational Field
        """
        try:
            return self.__boundary_space
        except AttributeError:
            pass
        self.__boundary_space = boundary.BoundarySpace_wtk_g1(
            self.level(), self.weight(), self.sign(), self.base_ring())
        return self.__boundary_space

    def manin_symbols(self):
        r"""
        Return the Manin symbol list of this modular symbol space.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(30),4)
            sage: M.manin_symbols()
            Manin Symbol List of weight 4 for Gamma1(30)
            sage: len(M.manin_symbols())
            1728
        """
        try:
            return self.__manin_symbols
        except AttributeError:
            self.__manin_symbols = ManinSymbolList_gamma1(
                level=self.level(), weight=self.weight())
        return self.__manin_symbols


class ModularSymbolsAmbient_wtk_gamma_h(ModularSymbolsAmbient):
    def __init__(self, group, weight, sign, F, custom_init=None, category=None):
        r"""
        Initialize a space of modular symbols for `\Gamma_H(N)`.

        INPUT:


        -  ``group`` - a congruence subgroup
           `\Gamma_H(N)`.

        -  ``weight`` - int, the weight = 2

        -  ``sign`` - int, either -1, 0, or 1

        -  ``F`` - field


        EXAMPLES::

            sage: ModularSymbols(GammaH(15,[4]),2)
            Modular Symbols space of dimension 9 for Congruence Subgroup Gamma_H(15) with H generated by [4] of weight 2 with sign 0 over Rational Field
        """
        ModularSymbolsAmbient.__init__(self,
                weight=weight, group=group,
                sign=sign, base_ring=F,
                custom_init=custom_init,
                category=category)

    def _dimension_formula(self):
        r"""
        Return None: we have no dimension formulas for `\Gamma_H(N)` spaces.

        EXAMPLES::

            sage: M = ModularSymbols(GammaH(15,[4]),2)
            sage: M.dimension()
            9
            sage: M._dimension_formula()
        """
        return None

    def _repr_(self):
        r"""
        Return a string representation of this space.

        EXAMPLES::

            sage: M = ModularSymbols(GammaH(15,[4]),2)
            sage: M # indirect doctest
            Modular Symbols space of dimension 9 for Congruence Subgroup Gamma_H(15) with H generated by [4] of weight 2 with sign 0 over Rational Field
        """
        return ("Modular Symbols space of dimension %s for %s of weight %s with sign %s over %s"
                % (self.dimension(), self.group(), self.weight(), self.sign(), self.base_ring()))

    def _cuspidal_submodule_dimension_formula(self):
        r"""
        Return None: we have no dimension formulas for `\Gamma_H(N)` spaces.

        EXAMPLES::

            sage: ModularSymbols(GammaH(15,[4]),2)._cuspidal_submodule_dimension_formula() is None
            True

        """
        return None

    def _cuspidal_new_submodule_dimension_formula(self):
        r"""
        Return None: we have no dimension formulas for `\Gamma_H(N)` spaces.

        EXAMPLES::

            sage: ModularSymbols(GammaH(15,[4]),2)._cuspidal_new_submodule_dimension_formula() is None
            True
        """
        return None

    def _compute_hecke_matrix_prime_power(self, p, r):
        r"""
        Return matrix of a prime-power Hecke operator.

        EXAMPLES::

            sage: M = ModularSymbols(GammaH(15,[4]),2)
            sage: M._compute_hecke_matrix_prime_power(2, 3)
            [10  0  5  1  0  0  0  4  0]
            [ 0 10  0  4 10 -5 -5 -4  5]
            [ 5  0 10 -4  0  0  0 -1  0]
            [ 0  0  0 -1  0  0  0 -4  0]
            [ 0  0  0 -7  5  0 10 -3  6]
            [ 0 -5  0 -6  0 10  5 -4 -1]
            [ 0  0  0 -3 10  0  5 -7  6]
            [ 0  0  0 -4  0  0  0 -1  0]
            [ 0  0  0  0  0  0  0  0  3]
            sage: M.hecke_matrix(7)^2  == M.hecke_matrix(49) + 7 * M.diamond_bracket_operator(7).matrix() # indirect doctest
            True
        """
        return self._compute_hecke_matrix_prime(p**r)

    def _degeneracy_raising_matrix_1(self, level):
        r"""
        Return matrix of a degeneracy raising map.

        EXAMPLES::

            sage: ModularSymbols(GammaH(15,[4]),2)._degeneracy_raising_matrix_1(ModularSymbols(GammaH(30, [19]), 2))
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def boundary_space(self):
        r"""
        Return the space of boundary modular symbols for this space.

        EXAMPLES::

            sage: M = ModularSymbols(GammaH(15,[4]),2)
            sage: M.boundary_space()
            Boundary Modular Symbols space for Congruence Subgroup Gamma_H(15) with H generated by [4] of weight 2 over Rational Field
        """
        try:
            return self.__boundary_space
        except AttributeError:
            pass
        self.__boundary_space = boundary.BoundarySpace_wtk_gamma_h(
            self.group(), self.weight(), self.sign(), self.base_ring())
        return self.__boundary_space

    def manin_symbols(self):
        r"""
        Return the Manin symbol list of this modular symbol space.

        EXAMPLES::

            sage: M = ModularSymbols(GammaH(15,[4]),2)
            sage: M.manin_symbols()
            Manin Symbol List of weight 2 for Congruence Subgroup Gamma_H(15) with H generated by [4]
            sage: len(M.manin_symbols())
            96
        """
        try:
            return self.__manin_symbols
        except AttributeError:
            self.__manin_symbols = ManinSymbolList_gamma_h(
                group=self.group(), weight=self.weight())
        return self.__manin_symbols


class ModularSymbolsAmbient_wtk_eps(ModularSymbolsAmbient):
    def __init__(self, eps, weight, sign, base_ring, custom_init=None, category=None):
        """
        Space of modular symbols with given weight, character, base ring and
        sign.

        INPUT:


        -  ``eps`` - dirichlet.DirichletCharacter, the
           "Nebentypus" character.

        -  ``weight`` - int, the weight = 2

        -  ``sign`` - int, either -1, 0, or 1

        - ``base_ring`` - the base ring. It must be possible to change the ring
          of the character to this base ring (not always canonically).


        EXAMPLES::

            sage: eps = DirichletGroup(4).gen(0)
            sage: eps.order()
            2
            sage: ModularSymbols(eps, 2)
            Modular Symbols space of dimension 0 and level 4, weight 2, character [-1], sign 0, over Rational Field
            sage: ModularSymbols(eps, 3)
            Modular Symbols space of dimension 2 and level 4, weight 3, character [-1], sign 0, over Rational Field

        We next create a space with character of order bigger than 2.

        ::

            sage: eps = DirichletGroup(5).gen(0)
            sage: eps     # has order 4
            Dirichlet character modulo 5 of conductor 5 mapping 2 |--> zeta4
            sage: ModularSymbols(eps, 2).dimension()
            0
            sage: ModularSymbols(eps, 3).dimension()
            2

        Here is another example::

            sage: G.<e> = DirichletGroup(5)
            sage: M = ModularSymbols(e,3)
            sage: loads(M.dumps()) == M
            True
        """
        level = eps.modulus()
        ModularSymbolsAmbient.__init__(self,
                weight = weight,
                group = arithgroup.Gamma1(level),
                sign = sign,
                base_ring = base_ring,
                character = eps.change_ring(base_ring),
                custom_init=custom_init,
                category=category)

    def _repr_(self):
        r"""
        Return a string representation of this space.

        EXAMPLES::

            sage: G.<e> = DirichletGroup(5)
            sage: M = ModularSymbols(e,3)
            sage: M # indirect doctest
            Modular Symbols space of dimension 2 and level 5, weight 3, character [zeta4], sign 0, over Cyclotomic Field of order 4 and degree 2
        """
        return ("Modular Symbols space of dimension %s and level %s, weight %s, character %s, sign %s, " + \
                "over %s")%(self.dimension(), self.level(), self.weight(),
                    self.character()._repr_short_(), self.sign(), self.base_ring())


    def _cuspidal_submodule_dimension_formula(self):
        r"""
        Return the dimension for the cuspidal subspace of this space, given by the formula.

        EXAMPLES::

            sage: G.<e> = DirichletGroup(50)
            sage: M = ModularSymbols(e^2,2)
            sage: M.dimension()
            16
            sage: M._cuspidal_submodule_dimension_formula()
            12
        """
        if self.base_ring().characteristic() != 0:
            raise NotImplementedError
        if self.sign() == 0:
            m = 2
        else:
            m = 1
        return m * self.group().dimension_cusp_forms(self.weight(), eps=self.character())

    def _cuspidal_new_submodule_dimension_formula(self):
        r"""
        Return the dimension for the new cuspidal subspace of this space, given by the formula.

        EXAMPLES::

            sage: G.<e> = DirichletGroup(50)
            sage: M = ModularSymbols(e,3)
            sage: M.dimension()
            30
            sage: M._cuspidal_new_submodule_dimension_formula()
            10
        """
        if self.base_ring().characteristic() != 0:
            raise NotImplementedError
        if self.sign() == 0:
            m = 2
        else:
            m = 1
        return m * self.group().dimension_new_cusp_forms(self.weight(), eps=self.character())

    def _matrix_of_operator_on_modular_symbols(self, codomain, R, character_twist=False):
        r"""
        INPUT:


        -  ``self`` - this space of modular symbols

        -  ``codomain`` - space of modular symbols

        -  ``R`` - list of lists [a,b,c,d] of length 4, which
           we view as elements of GL_2(Q).


        OUTPUT: a matrix, which represents the operator

        .. MATH::

                            x \mapsto \sum_{g in R} g.x


        where g.x is the formal linear fractional transformation on modular
        symbols.

        EXAMPLES::

            sage: G.<e> = DirichletGroup(5)
            sage: M = ModularSymbols(e,3)
            sage: M.dimension()
            2
            sage: M._matrix_of_operator_on_modular_symbols(M,HeilbronnCremona(3))
            [ 6  6]
            [ 0 10]
        """
        eps = self.character()
        rows = []
        for b in self.basis():
            v = formal_sum.FormalSum(0, check=False)
            for c, x in b.modular_symbol_rep():
                for g in R:
                    y = x.apply(g)
                    if character_twist:
                        v += y*c*eps(g[0])
                    else:
                        v += y*c
            w = codomain(v).element()
            rows.append(w)
        M = matrix_space.MatrixSpace(self.base_ring(), len(rows), codomain.degree(), sparse=False)
        return M(rows)

    def _degeneracy_raising_matrix_1(self, M):
        r"""
        Return the matrix of the degeneracy raising map to ``M``, which should
        be a space of modular symbols with level a multiple of the level of
        self and with compatible character.

        INPUT:

        - ``M`` -- a space of modular symbols with character, whose level
          should be an integer multiple of the level of self, and whose
          character should be the Dirichlet character at that level obtained by
          extending the character of self.

        The input is *not* sanity-checked in any way -- use with care!

        OUTPUT:

        (matrix) The matrix of the degeneracy raising matrix to the higher level.

        EXAMPLES::

            sage: eps = DirichletGroup(4).gen(0)
            sage: M = ModularSymbols(eps, 3); M
            Modular Symbols space of dimension 2 and level 4, weight 3, character [-1], sign 0, over Rational Field
            sage: M._degeneracy_raising_matrix_1(ModularSymbols(eps.extend(20), 3))
            [ 1  0  0  0 -1 -1  3  1  0  2 -3  0]
            [ 0  5  1 -2 -3  3  0  4 -1  5 -7 -1]
        """
        N = self.level()

        # 1. Find coset representatives H for Gamma_0(M.level()) \ Gamma_0(self.level())
        H = arithgroup.degeneracy_coset_representatives_gamma0(M.level(), N, 1)
        # 2. The map is
        #        [P,pi(g)] |--> sum_{h in H} [P, pi(h*g)]
        #
        MS = matrix_space.MatrixSpace(self.base_ring(), self.dimension(), M.dimension())
        if self.dimension() == 0 or M.dimension() == 0:
            return MS(0)
        rows = []
        B = self.manin_basis()
        syms = self.manin_symbols()
        G = matrix_space.MatrixSpace(ZZ, 2)
        H = [G(h) for h in H]
        eps = self.character()  # note: in my thesis I twisted by eps^(-1), which is definitely a mistake
                                # since twisting by eps gives the right answer and by eps^(-1) does not.
        for n in B:
            z = M(0)
            s = syms.manin_symbol(n)
            g = G(list(s.lift_to_sl2z(N)))
            i = s.i
            # We apply each matrix in H according to the above formula
            for h in H:
                hg = h*g
                z += eps(h[0,0])*M((i, hg[1,0], hg[1,1]))
            rows.append(z.element())
        A = MS(rows)
        return A

    def _dimension_formula(self):
        r"""
        Return None: we have no dimension formula for `\Gamma_H(N)` spaces.

        EXAMPLES::

            sage: eps = DirichletGroup(5).gen(0)
            sage: M = ModularSymbols(eps, 2)
            sage: M.dimension()
            0
            sage: M._dimension_formula()
        """
        return None

    def boundary_space(self):
        r"""
        Return the space of boundary modular symbols for this space.

        EXAMPLES::

            sage: eps = DirichletGroup(5).gen(0)
            sage: M = ModularSymbols(eps, 2)
            sage: M.boundary_space()
            Boundary Modular Symbols space of level 5, weight 2, character [zeta4] and dimension 0 over Cyclotomic Field of order 4 and degree 2
        """
        try:
            return self.__boundary_space
        except AttributeError:
            pass
        self.__boundary_space = boundary.BoundarySpace_wtk_eps(
            self.character(), self.weight(), self.sign())
        return self.__boundary_space

    def manin_symbols(self):
        r"""
        Return the Manin symbol list of this modular symbol space.

        EXAMPLES::

            sage: eps = DirichletGroup(5).gen(0)
            sage: M = ModularSymbols(eps, 2)
            sage: M.manin_symbols()
            Manin Symbol List of weight 2 for Gamma1(5) with character [zeta4]
            sage: len(M.manin_symbols())
            6
        """
        try:
            return self.__manin_symbols
        except AttributeError:
            self.__manin_symbols = ManinSymbolList_character(
                character=self.character(), weight=self.weight())
        return self.__manin_symbols

    def modular_symbols_of_level(self, N):
        r"""
        Return a space of modular symbols with the same parameters as
        this space except with level `N`.

        INPUT:

        - ``N`` (int) -- a positive integer.

        OUTPUT:

        (Modular Symbol space) A space of modular symbols with the
        same defining properties (weight, sign, etc.) as this space
        except with level `N`.

        EXAMPLES::

            sage: eps = DirichletGroup(5).gen(0)
            sage: M = ModularSymbols(eps, 2); M
            Modular Symbols space of dimension 0 and level 5, weight 2, character [zeta4], sign 0, over Cyclotomic Field of order 4 and degree 2
            sage: M.modular_symbols_of_level(15)
            Modular Symbols space of dimension 0 and level 15, weight 2, character [1, zeta4], sign 0, over Cyclotomic Field of order 4 and degree 2
        """
        if self.level() % N == 0:
            eps = self.character().restrict(N)
        elif N % self.level() == 0:
            eps = self.character().extend(N)
        else:
            raise ValueError("The level N (=%s) must be a divisor or multiple of the modulus of the character (=%s)"%(N, self.level()))
        return modsym.ModularSymbols(eps, self.weight(), self.sign(), self.base_ring())

    def modular_symbols_of_sign(self, sign):
        r"""
        Return a space of modular symbols with the same defining
        properties (weight, level, etc.) as this space except with given
        sign.

        INPUT:

        - ``sign`` (int) -- A sign (`+1`, `-1` or `0`).

        OUTPUT:

        (ModularSymbolsAmbient) A space of modular symbols with the
        same defining properties (weight, level, etc.) as this space
        except with given sign.

        EXAMPLES::

            sage: eps = DirichletGroup(5).gen(0)
            sage: M = ModularSymbols(eps, 2); M
            Modular Symbols space of dimension 0 and level 5, weight 2, character [zeta4], sign 0, over Cyclotomic Field of order 4 and degree 2
            sage: M.modular_symbols_of_sign(0) == M
            True
            sage: M.modular_symbols_of_sign(+1)
            Modular Symbols space of dimension 0 and level 5, weight 2, character [zeta4], sign 1, over Cyclotomic Field of order 4 and degree 2
            sage: M.modular_symbols_of_sign(-1)
            Modular Symbols space of dimension 0 and level 5, weight 2, character [zeta4], sign -1, over Cyclotomic Field of order 4 and degree 2

        """
        return modsym.ModularSymbols(self.character(), self.weight(), sign, self.base_ring())

    def modular_symbols_of_weight(self, k):
        r"""
        Return a space of modular symbols with the same defining
        properties (weight, sign, etc.) as this space except with weight
        `k`.

        INPUT:

        - ``k`` (int) -- A positive integer.

        OUTPUT:

        (ModularSymbolsAmbient) A space of modular symbols with the
        same defining properties (level, sign) as this space
        except with given weight.

        EXAMPLES::

            sage: eps = DirichletGroup(5).gen(0)
            sage: M = ModularSymbols(eps, 2); M
            Modular Symbols space of dimension 0 and level 5, weight 2, character [zeta4], sign 0, over Cyclotomic Field of order 4 and degree 2
            sage: M.modular_symbols_of_weight(3)
            Modular Symbols space of dimension 2 and level 5, weight 3, character [zeta4], sign 0, over Cyclotomic Field of order 4 and degree 2
            sage: M.modular_symbols_of_weight(2) == M
            True
        """
        return modsym.ModularSymbols(self.character(), k, self.sign(), self.base_ring())

    def _hecke_images(self, i, v):
        """
        Return images of the `i`-th standard basis vector under the
        Hecke operators `T_p` for all integers in `v`.

        INPUT:

        - ``i`` -- nonnegative integer

        - ``v`` -- a list of positive integer

        OUTPUT:

        - ``matrix`` -- whose rows are the Hecke images

        EXAMPLES::

            sage: G.<e> = DirichletGroup(50,QQ)
            sage: M = ModularSymbols(e^2,2)
            sage: M.dimension()
            15
            sage: M._hecke_images(8,list(range(1,5)))
            [ 0  0  0  0  0  0  0  0  1  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  1  0  0  0  0  1  0  0]
            [ 0  1  0  2  0 -1  1  1  0  0  0  0  0  0  0]
            [ 0  1  1 -1 -1  0 -1  1  1  0  1  2  0 -2  2]
        """
        if self.weight() != 2:
            raise NotImplementedError("hecke images only implemented when the weight is 2")
        chi = self.character()
        # Find basis vector for ambient space such that it is not in
        # the kernel of the dual space corresponding to self.
        c = self.manin_generators()[self.manin_basis()[i]]
        N = self.level()
        if chi.order() > 2:
            return heilbronn.hecke_images_nonquad_character_weight2(c.u,c.v,N,
                                 v, chi, self.manin_gens_to_basis())
        else:
            return heilbronn.hecke_images_quad_character_weight2(c.u,c.v,N,
                                 v, chi, self.manin_gens_to_basis())
        raise NotImplementedError
