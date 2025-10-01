# sage.doctest: needs sage.rings.padics
r"""
Pollack-Stevens' modular symbols spaces

This module contains a class for spaces of modular symbols that use Glenn
Stevens' conventions, as explained in [PS2011]_.

There are two main differences between the modular symbols in this directory
and the ones in :mod:`sage.modular.modsym`:

- There is a shift in the weight: weight `k=0` here corresponds to weight `k=2`
  there.

- There is a duality: these modular symbols are functions from
  `\textrm{Div}^0(P^1(\QQ))` (cohomological objects), the others are formal linear
  combinations of `\textrm{Div}^0(P^1(\QQ))` (homological objects).

EXAMPLES:

First we create the space of modular symbols of weight 0 (`k=2`) and level 11::

    sage: M = PollackStevensModularSymbols(Gamma0(11), 0); M
    Space of modular symbols for Congruence Subgroup Gamma0(11) with sign 0 and values in Sym^0 Q^2

One can also create a space of overconvergent modular symbols, by specifying a prime and a precision::

    sage: M = PollackStevensModularSymbols(Gamma0(11), p = 5, prec_cap = 10, weight = 0); M
    Space of overconvergent modular symbols for Congruence Subgroup Gamma0(11) with sign 0 and values in Space of 5-adic distributions with k=0 action and precision cap 10

Currently not much functionality is available on the whole space, and these
spaces are mainly used as parents for the modular symbols. These can be constructed from the corresponding
classical modular symbols (or even elliptic curves) as follows::

    sage: A = ModularSymbols(13, sign=1, weight=4).decomposition()[0]
    sage: A.is_cuspidal()
    True
    sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
    sage: f = ps_modsym_from_simple_modsym_space(A); f
    Modular symbol of level 13 with values in Sym^2 Q^2
    sage: f.values()
    [(-13, 0, -1),
     (247/2, 13/2, -6),
     (39/2, 117/2, 42),
     (-39/2, 39, 111/2),
     (-247/2, -117, -209/2)]
    sage: f.parent()
    Space of modular symbols for Congruence Subgroup Gamma0(13) with sign 1 and values in Sym^2 Q^2

::

    sage: E = EllipticCurve('37a1')
    sage: phi = E.pollack_stevens_modular_symbol(); phi
    Modular symbol of level 37 with values in Sym^0 Q^2
    sage: phi.values()
    [0, 1, 0, 0, 0, -1, 1, 0, 0]
    sage: phi.parent()
    Space of modular symbols for Congruence Subgroup Gamma0(37) with sign 0 and values in Sym^0 Q^2
"""
# ****************************************************************************
#       Copyright (C) 2012 Robert Pollack <rpollack@math.bu.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.modules.module import Module
from sage.modular.dirichlet import DirichletCharacter
from sage.modular.arithgroup.congroup_gamma0 import Gamma0_constructor as Gamma0
from sage.modular.arithgroup.arithgroup_element import ArithmeticSubgroupElement
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from .fund_domain import ManinRelations
from sage.rings.infinity import infinity as oo
from sage.structure.factory import UniqueFactory

from .distributions import OverconvergentDistributions, Symk
from .modsym import (PSModularSymbolElement, PSModularSymbolElement_symk,
                     PSModularSymbolElement_dist, PSModSymAction)
from .manin_map import ManinMap
from .sigma0 import Sigma0, Sigma0Element


class PollackStevensModularSymbols_factory(UniqueFactory):
    r"""
    Create a space of Pollack-Stevens modular symbols.

    INPUT:

    - ``group`` -- integer or congruence subgroup

    - ``weight`` -- integer `\ge 0`, or ``None``

    - ``sign`` -- integer; -1, 0, 1

    - ``base_ring`` -- ring or ``None``

    - ``p`` -- prime or ``None``

    - ``prec_cap`` -- positive integer or ``None``

    - ``coefficients`` -- the coefficient module (a special type of module,
      typically distributions), or ``None``

    If an explicit coefficient module is given, then the arguments ``weight``,
    ``base_ring``, ``prec_cap``, and ``p`` are redundant and must be ``None``.
    They are only relevant if ``coefficients`` is ``None``, in which case the
    coefficient module is inferred from the other data.

    .. NOTE::

        We emphasize that in the Pollack-Stevens notation, the
        ``weight`` is the usual weight minus 2, so a classical weight
        2 modular form corresponds to a modular symbol of "weight 0".

    EXAMPLES::

        sage: M = PollackStevensModularSymbols(Gamma0(7), weight=0, prec_cap = None); M
        Space of modular symbols for Congruence Subgroup Gamma0(7) with sign 0 and values in Sym^0 Q^2

    An example with an explicit coefficient module::

        sage: D = OverconvergentDistributions(3, 7, prec_cap=10)
        sage: M = PollackStevensModularSymbols(Gamma0(7), coefficients=D); M
        Space of overconvergent modular symbols for Congruence Subgroup Gamma0(7) with sign 0 and values in Space of 7-adic distributions with k=3 action and precision cap 10

    TESTS::

        sage: TestSuite(PollackStevensModularSymbols).run()
    """
    def create_key(self, group, weight=None, sign=0, base_ring=None, p=None, prec_cap=None, coefficients=None):
        r"""
        Sanitize input.

        EXAMPLES::

            sage: D = OverconvergentDistributions(3, 7, prec_cap=10)
            sage: M = PollackStevensModularSymbols(Gamma0(7), coefficients=D) # indirect doctest
        """
        if sign not in (-1, 0, 1):
            raise ValueError("sign must be -1, 0, 1")

        if isinstance(group, (int, Integer)):
            group = Gamma0(group)

        if coefficients is None:
            if isinstance(group, DirichletCharacter):
                character = group.minimize_base_ring()
                group = Gamma0(character.modulus())
                if character.is_trivial():
                    character = None
            else:
                character = None

            if weight is None:
                raise ValueError("you must specify a weight "
                                 "or coefficient module")

            if prec_cap is None:
                coefficients = Symk(weight, base_ring, character)
            else:
                coefficients = OverconvergentDistributions(weight, p, prec_cap, base_ring,
                                             character)
        else:
            if weight is not None or base_ring is not None or p is not None or prec_cap is not None:
                raise ValueError("if coefficients are specified, then weight, "
                                 "base_ring, p, and prec_cap must take their "
                                 "default value None")

        return (group, coefficients, sign)

    def create_object(self, version, key):
        r"""
        Create a space of modular symbols from ``key``.

        INPUT:

        - ``version`` -- the version of the object to create

        - ``key`` -- tuple of parameters, as created by :meth:`create_key`

        EXAMPLES::

            sage: D = OverconvergentDistributions(5, 7, 15)
            sage: M = PollackStevensModularSymbols(Gamma0(7), coefficients=D) # indirect doctest
            sage: M2 = PollackStevensModularSymbols(Gamma0(7), coefficients=D) # indirect doctest
            sage: M is M2
            True
        """
        return PollackStevensModularSymbolspace(*key)


PollackStevensModularSymbols = PollackStevensModularSymbols_factory('PollackStevensModularSymbols')


class PollackStevensModularSymbolspace(Module):
    r"""
    A class for spaces of modular symbols that use Glenn Stevens' conventions.
    This class should not be instantiated directly by the user: this is handled
    by the factory object :class:`PollackStevensModularSymbols_factory`.

    INPUT:

    - ``group`` -- congruence subgroup

    - ``coefficients`` -- a coefficient module

    - ``sign`` -- (default: 0) 0, -1, or 1

    EXAMPLES::

        sage: D = OverconvergentDistributions(2, 11)
        sage: M = PollackStevensModularSymbols(Gamma0(2), coefficients=D); M.sign()
        0
        sage: M = PollackStevensModularSymbols(Gamma0(2), coefficients=D, sign=-1); M.sign()
        -1
        sage: M = PollackStevensModularSymbols(Gamma0(2), coefficients=D, sign=1); M.sign()
        1
    """
    def __init__(self, group, coefficients, sign=0):
        r"""
        INPUT:

        See :class:`PollackStevensModularSymbolspace`

        EXAMPLES::

            sage: D = OverconvergentDistributions(2, 11)
            sage: M = PollackStevensModularSymbols(Gamma0(11), coefficients=D)
            sage: type(M)
            <class 'sage.modular.pollack_stevens.space.PollackStevensModularSymbolspace_with_category'>
            sage: TestSuite(M).run()
        """
        Module.__init__(self, coefficients.base_ring())
        if sign not in [0, -1, 1]:
            # sign must be 0, -1 or 1
            raise ValueError("sign must be 0, -1, or 1")
        self._group = group
        self._coefficients = coefficients
        if coefficients.is_symk():
            self.Element = PSModularSymbolElement_symk
        else:
            self.Element = PSModularSymbolElement_dist
        self._sign = sign
        # should distinguish between Gamma0 and Gamma1...
        self._source = ManinRelations(group.level())

        # Register the action of 2x2 matrices on self.

        if coefficients.is_symk():
            action = PSModSymAction(Sigma0(1), self)
        else:
            action = PSModSymAction(Sigma0(self.prime()), self)

        self._populate_coercion_lists_(action_list=[action])

    def _element_constructor_(self, data):
        r"""
        Construct an element of ``self`` from data.

        EXAMPLES::

            sage: D = OverconvergentDistributions(0, 11)
            sage: M = PollackStevensModularSymbols(Gamma0(11), coefficients=D)
            sage: M(1) # indirect doctest
            Modular symbol of level 11 with values in Space of 11-adic distributions with k=0 action and precision cap 20
        """
        if isinstance(data, PSModularSymbolElement):
            data = data._map
        elif isinstance(data, ManinMap):
            pass
        else:
            # a dict, or a single distribution specifying a constant symbol, etc.
            data = ManinMap(self._coefficients, self._source, data)

        if data._codomain != self._coefficients:
            data = data.extend_codomain(self._coefficients)

        return self.element_class(data, self, construct=True)

    def _coerce_map_from_(self, other):
        r"""
        Used for comparison and coercion.

        EXAMPLES::

            sage: M1 = PollackStevensModularSymbols(Gamma0(11), coefficients=Symk(3))
            sage: M2 = PollackStevensModularSymbols(Gamma0(11), coefficients=Symk(3,Qp(11)))
            sage: M3 = PollackStevensModularSymbols(Gamma0(11), coefficients=Symk(4))
            sage: M4 = PollackStevensModularSymbols(Gamma0(11), coefficients=OverconvergentDistributions(3, 11, 10))
            sage: M1.has_coerce_map_from(M2)
            False
            sage: M2.has_coerce_map_from(M1)
            True
            sage: M1.has_coerce_map_from(M3)
            False
            sage: M1.has_coerce_map_from(M4)
            False
            sage: M2.has_coerce_map_from(M4)
            True
        """
        if isinstance(other, PollackStevensModularSymbolspace):
            return (other.group() == self.group()
                    and self.coefficient_module().has_coerce_map_from(other.coefficient_module()))

        return False

    def _repr_(self):
        r"""
        Return string representation.

        EXAMPLES::

            sage: D = OverconvergentDistributions(2, 11)
            sage: M = PollackStevensModularSymbols(Gamma0(2), coefficients=D)
            sage: M._repr_()
            'Space of overconvergent modular symbols for Congruence Subgroup Gamma0(2) with sign 0 and values in Space of 11-adic distributions with k=2 action and precision cap 20'
        """
        if self.coefficient_module().is_symk():
            s = "Space of modular symbols for "
        else:
            s = "Space of overconvergent modular symbols for "
        s += "%s with sign %s and values in %s" % (self.group(), self.sign(),
                                                   self.coefficient_module())
        return s

    def source(self):
        r"""
        Return the domain of the modular symbols in this space.

        OUTPUT: a :class:`sage.modular.pollack_stevens.fund_domain.PollackStevensModularDomain`

        EXAMPLES::

            sage: D = OverconvergentDistributions(2, 11)
            sage: M = PollackStevensModularSymbols(Gamma0(2), coefficients=D)
            sage: M.source()
            Manin Relations of level 2
        """
        return self._source

    def coefficient_module(self):
        r"""
        Return the coefficient module of this space.

        EXAMPLES::

            sage: D = OverconvergentDistributions(2, 11)
            sage: M = PollackStevensModularSymbols(Gamma0(2), coefficients=D)
            sage: M.coefficient_module()
            Space of 11-adic distributions with k=2 action and precision cap 20
            sage: M.coefficient_module() is D
            True
        """
        return self._coefficients

    def group(self):
        r"""
        Return the congruence subgroup of this space.

        EXAMPLES::

            sage: D = OverconvergentDistributions(2, 5)
            sage: G = Gamma0(23)
            sage: M = PollackStevensModularSymbols(G, coefficients=D)
            sage: M.group()
            Congruence Subgroup Gamma0(23)
            sage: D = Symk(4)
            sage: G = Gamma1(11)
            sage: M = PollackStevensModularSymbols(G, coefficients=D)
            sage: M.group()
            Congruence Subgroup Gamma1(11)
        """
        return self._group

    def sign(self):
        r"""
        Return the sign of this space.

        EXAMPLES::

            sage: D = OverconvergentDistributions(3, 17)
            sage: M = PollackStevensModularSymbols(Gamma(5), coefficients=D)
            sage: M.sign()
            0
            sage: D = Symk(4)
            sage: M = PollackStevensModularSymbols(Gamma1(8), coefficients=D, sign=-1)
            sage: M.sign()
            -1
        """
        return self._sign

    def ngens(self):
        r"""
        Return the number of generators defining this space.

        EXAMPLES::

            sage: D = OverconvergentDistributions(4, 29)
            sage: M = PollackStevensModularSymbols(Gamma1(12), coefficients=D)
            sage: M.ngens()
            5
            sage: D = Symk(2)
            sage: M = PollackStevensModularSymbols(Gamma0(2), coefficients=D)
            sage: M.ngens()
            2
        """
        return len(self._source.indices())

    def ncoset_reps(self):
        r"""
        Return the number of coset representatives defining the domain of the
        modular symbols in this space.

        OUTPUT:

        The number of coset representatives stored in the manin relations.
        (Just the size of `P^1(\ZZ/N\ZZ)`)

        EXAMPLES::

            sage: D = Symk(2)
            sage: M = PollackStevensModularSymbols(Gamma0(2), coefficients=D)
            sage: M.ncoset_reps()
            3
        """
        return len(self._source.reps())

    def level(self):
        r"""
        Return the level `N`, where this space is of level `\Gamma_0(N)`.

        EXAMPLES::

            sage: D = OverconvergentDistributions(7, 11)
            sage: M = PollackStevensModularSymbols(Gamma1(14), coefficients=D)
            sage: M.level()
            14
        """
        return self._source.level()

    def _grab_relations(self):
        r"""
        This is used internally as part of a consistency check.

        EXAMPLES::

            sage: D = OverconvergentDistributions(4, 3)
            sage: M = PollackStevensModularSymbols(Gamma1(13), coefficients=D)
            sage: M._grab_relations()
            [[(1, [1 0]
            [0 1], 0)], [(-1, [-1 -1]
            [ 0 -1], 0)], [(1, [1 0]
            [0 1], 2)], [(1, [1 0]
                [0 1], 3)], [(1, [1 0]
                    [0 1], 4)], [(1, [1 0]
                        [0 1], 5)]]
        """
        S0N = Sigma0(self._source._N)
        v = []
        for r in range(len(self._source.gens())):
            for j in range(len(self._source.reps())):
                R = self._source.relations(j)
                if len(R) == 1 and R[0][2] == self._source.indices(r):
                    if R[0][0] != -1 or R[0][1] != S0N(1):
                        v += [R]
        return v

    def precision_cap(self):
        r"""
        Return the number of moments of each element of this space.

        EXAMPLES::

            sage: D = OverconvergentDistributions(2, 5)
            sage: M = PollackStevensModularSymbols(Gamma1(13), coefficients=D)
            sage: M.precision_cap()
            20
            sage: D = OverconvergentDistributions(3, 7, prec_cap=10)
            sage: M = PollackStevensModularSymbols(Gamma0(7), coefficients=D)
            sage: M.precision_cap()
            10
        """
        # WARNING -- IF YOU ARE WORKING IN SYM^K(Q^2) THIS WILL JUST
        # RETURN K-1.  NOT GOOD
        return self.coefficient_module()._prec_cap

    def weight(self):
        r"""
        Return the weight of this space.

        .. WARNING::

            We emphasize that in the Pollack-Stevens notation, this is
            the usual weight minus 2, so a classical weight 2 modular
            form corresponds to a modular symbol of "weight 0".

        EXAMPLES::

            sage: D = Symk(5)
            sage: M = PollackStevensModularSymbols(Gamma1(7), coefficients=D)
            sage: M.weight()
            5
        """
        return self.coefficient_module()._k

    def prime(self):
        r"""
        Return the prime of this space.

        EXAMPLES::

            sage: D = OverconvergentDistributions(2, 11)
            sage: M = PollackStevensModularSymbols(Gamma(2), coefficients=D)
            sage: M.prime()
            11
        """
        return self.coefficient_module()._p

    def _p_stabilize_parent_space(self, p, new_base_ring):
        r"""
        Return the space of Pollack-Stevens modular symbols of level
        `p N`, with changed base ring.  This is used internally when
        constructing the `p`-stabilization of a modular symbol.

        INPUT:

        - ``p`` -- prime number
        - ``new_base_ring`` -- the base ring of the result

        OUTPUT:

        The space of modular symbols of level `p N`, where `N` is the level
        of this space.

        EXAMPLES::

            sage: D = OverconvergentDistributions(2, 7)
            sage: M = PollackStevensModularSymbols(Gamma(13), coefficients=D)
            sage: M._p_stabilize_parent_space(7, M.base_ring())
            Space of overconvergent modular symbols for Congruence Subgroup
            Gamma(91) with sign 0 and values in Space of 7-adic distributions
            with k=2 action and precision cap 20

            sage: D = OverconvergentDistributions(4, 17)
            sage: M = PollackStevensModularSymbols(Gamma1(3), coefficients=D)
            sage: M._p_stabilize_parent_space(17, Qp(17))
            Space of overconvergent modular symbols for Congruence
            Subgroup Gamma1(51) with sign 0 and values in Space of
            17-adic distributions with k=4 action and precision cap 20
        """
        N = self.level()
        if N % p == 0:
            raise ValueError("the level is not prime to p")
        from sage.modular.arithgroup.congroup_gamma0 import Gamma0_class, Gamma0_constructor as Gamma0
        from sage.modular.arithgroup.congroup_gamma1 import Gamma1_class, Gamma1_constructor as Gamma1
        from sage.modular.arithgroup.congroup_gamma import Gamma_class, Gamma_constructor as Gamma

        G = self.group()
        if isinstance(G, Gamma0_class):
            G = Gamma0(N * p)
        elif isinstance(G, Gamma1_class):
            G = Gamma1(N * p)
        elif isinstance(G, Gamma_class):
            G = Gamma(N * p)
        else:
            raise NotImplementedError
        return PollackStevensModularSymbols(G, coefficients=self.coefficient_module().change_ring(new_base_ring), sign=self.sign())

    def _specialize_parent_space(self, new_base_ring):
        r"""
        Internal function that is used by the specialize method on
        elements.  It returns a space with same parameters as this
        one, but over ``new_base_ring``.

        INPUT:

        - ``new_base_ring`` -- a ring

        OUTPUT: a space of modular symbols to which our space specializes

        EXAMPLES::

            sage: D = OverconvergentDistributions(7, 5)
            sage: M = PollackStevensModularSymbols(Gamma0(2), coefficients=D); M
            Space of overconvergent modular symbols for Congruence Subgroup Gamma0(2) with sign 0 and values in Space of 5-adic distributions with k=7 action and precision cap 20
            sage: M._specialize_parent_space(QQ)
            Space of modular symbols for Congruence Subgroup Gamma0(2) with sign 0 and values in Sym^7 Q^2
            sage: M.base_ring()
            5-adic Ring with capped absolute precision 20
            sage: M._specialize_parent_space(QQ).base_ring()
            Rational Field
        """
        return PollackStevensModularSymbols(self.group(), coefficients=self.coefficient_module().specialize(new_base_ring), sign=self.sign())

    def _lift_parent_space(self, p, M, new_base_ring):
        r"""
        Used internally to lift a space of modular symbols to space of
        overconvergent modular symbols.

        INPUT:

        - ``p`` -- prime
        - ``M`` -- precision cap
        - ``new_base_ring`` -- ring

        OUTPUT: a space of distribution valued modular symbols

        EXAMPLES::

            sage: D = OverconvergentDistributions(4, 17, 2); M = PollackStevensModularSymbols(Gamma1(3), coefficients=D)
            sage: D.is_symk()
            False
            sage: M._lift_parent_space(17, 10, Qp(17))
            Traceback (most recent call last):
            ...
            TypeError: Coefficient module must be a Symk
            sage: PollackStevensModularSymbols(Gamma1(3), weight=1)._lift_parent_space(17,10,Qp(17))
            Space of overconvergent modular symbols for Congruence Subgroup Gamma1(3) with sign 0 and values in Space of 17-adic distributions with k=1 action and precision cap 10
        """
        if self.coefficient_module().is_symk():
            return PollackStevensModularSymbols(self.group(), coefficients=self.coefficient_module().lift(p, M, new_base_ring), sign=self.sign())
        else:
            raise TypeError("Coefficient module must be a Symk")

    def change_ring(self, new_base_ring):
        r"""
        Change the base ring of this space to ``new_base_ring``.

        INPUT:

        - ``new_base_ring`` -- a ring

        OUTPUT: a space of modular symbols over the specified base

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Symk
            sage: D = Symk(4)
            sage: M = PollackStevensModularSymbols(Gamma(6), coefficients=D); M
            Space of modular symbols for Congruence Subgroup Gamma(6) with sign 0 and values in Sym^4 Q^2
            sage: M.change_ring(Qp(5,8))
            Space of modular symbols for Congruence Subgroup Gamma(6) with sign 0 and values in Sym^4 Q_5^2
        """
        return PollackStevensModularSymbols(self.group(), coefficients=self.coefficient_module().change_ring(new_base_ring), sign=self.sign())

    def _an_element_(self):
        r"""
        Return the cusps associated to an element of a congruence subgroup.

        OUTPUT: an element of the modular symbol space

        This returns a "typical" element of this space; in this case
        the constant map sending every element to an element of the
        coefficient module.

        .. WARNING::

            This is not really an element of the space because it does
            not satisfy the Manin relations.

        EXAMPLES::

            sage: D = Symk(4)
            sage: M = PollackStevensModularSymbols(Gamma(6), coefficients=D)
            sage: x = M.an_element(); x       # indirect doctest
            Modular symbol of level 6 with values in Sym^4 Q^2
            sage: x.values()
            [(0, 1, 2, 3, 4), (0, 1, 2, 3, 4), (0, 1, 2, 3, 4)]
            sage: D = Symk(2, Qp(11)); M = PollackStevensModularSymbols(Gamma0(2), coefficients=D)
            sage: x = M.an_element(); x.values()
            [(0, 1 + O(11^20), 2 + O(11^20)), (0, 1 + O(11^20), 2 + O(11^20))]
            sage: x in M
            True
        """
        return self(self.coefficient_module().an_element())

    def random_element(self, M=None):
        r"""
        Return a random overconvergent modular symbol in this space with `M` moments.

        INPUT:

        - ``M`` -- positive integer

        OUTPUT: an element of the modular symbol space with `M` moments

        Returns a random element in this space by randomly choosing
        values of distributions on all but one divisor, and solves the
        difference equation to determine the value on the last
        divisor. ::

            sage: D = OverconvergentDistributions(2, 11)
            sage: M = PollackStevensModularSymbols(Gamma0(11), coefficients=D)
            sage: M.random_element(10)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        # This function still has bugs and is not used in the rest of
        # the package. It is left to be implemented in the future.
        raise NotImplementedError

        if M is None and not self.coefficient_module().is_symk():
            M = self.coefficient_module().precision_cap()

        k = self.coefficient_module()._k
        # p = self.prime()
        manin = self.source()

#        # There must be a problem here with that +1 -- should be
#        # variable depending on a c of some matrix We'll need to
#        # divide by some power of p and so we add extra accuracy
#        # here.
#        if k != 0:
#            MM = M + valuation(k,p) + 1 + M.exact_log(p)
#        else:
#            MM = M + M.exact_log(p) + 1

        # this loop runs thru all of the generators (except
        # (0)-(infty)) and randomly chooses a distribution to assign
        # to this generator (in the 2,3-torsion cases care is taken
        # to satisfy the relevant relation)
        D = {}
        for g in manin.gens():
            D[g] = self.coefficient_module().random_element(M)
            if g in manin.reps_with_two_torsion() and g in manin.reps_with_three_torsion():
                raise ValueError("Level 1 not implemented")
            if g in manin.reps_with_two_torsion():
                gamg = manin.two_torsion_matrix(g)
                D[g] = D[g] - D[g] * gamg
            else:
                if g in manin.reps_with_three_torsion():
                    gamg = manin.three_torsion_matrix(g)
                    D[g] = 2 * D[g] - D[g] * gamg - D[g] * gamg ** 2
                    #            print("post:",D[g])

        # now we compute nu_infty of Prop 5.1 of [PS1]
        t = self.coefficient_module().zero()
        for g in manin.gens()[1:]:
            if (g not in manin.reps_with_two_torsion()) and (g not in manin.reps_with_three_torsion()):
                t += D[g] * manin.gammas[g] - D[g]
            else:
                # this was previously MR.reps_with_two_torsion() but there is no variable MR defined...
                if g in manin.reps_with_two_torsion():
                    t -= D[g]
                else:
                    t -= D[g]

        # If k = 0, then t has total measure zero.  However, this is not true when k != 0
        # (unlike Prop 5.1 of [PS1] this is not a lift of classical symbol).
        # So instead we simply add (const)*mu_1 to some (non-torsion) v[j] to fix this
        # here since (mu_1 |_k ([a,b,c,d]-1))(trivial char) = chi(a) k a^{k-1} c ,
        # we take the constant to be minus the total measure of t divided by (chi(a) k a^{k-1} c)

        if k != 0:
            j = 1
            g = manin.gens()[j]
            while (g in manin.reps_with_two_torsion()) or (g in manin.reps_with_three_torsion()) and (j < len(manin.gens())):
                j = j + 1
                g = manin.gens()[j]
            if j == len(manin.gens()):
                raise ValueError("everything is 2 or 3 torsion!  NOT YET IMPLEMENTED IN THIS CASE")

            gam = manin.gammas[g]
            a = gam.matrix()[0, 0]
            c = gam.matrix()[1, 0]

            if self.coefficient_module()._character is not None:
                chara = self.coefficient_module()._character(a)
            else:
                chara = 1
            err = -t.moment(0) / (chara * k * a ** (k - 1) * c)
            v = [0] * M
            v[1] = 1
            mu_1 = self.base_ring()(err) * self.coefficient_module()(v)
            D[g] += mu_1
            t = t + mu_1 * gam - mu_1

        Id = manin.gens()[0]
        if not self.coefficient_module().is_symk():
            mu = t.solve_difference_equation()
            D[Id] = -mu
        else:
            if self.coefficient_module()._k == 0:
                D[Id] = self.coefficient_module().random_element()
            else:
                raise ValueError("Not implemented for symk with k>0 yet")

        return self(D)


def cusps_from_mat(g):
    r"""
    Return the cusps associated to an element of a congruence subgroup.

    INPUT:

    - ``g`` -- an element of a congruence subgroup or a matrix

    OUTPUT: a tuple of cusps associated to ``g``

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.space import cusps_from_mat
        sage: g = SL2Z.one()
        sage: cusps_from_mat(g)
        (+Infinity, 0)

    You can also just give the matrix of ``g``::

        sage: type(g)
        <class 'sage.modular.arithgroup.arithgroup_element.ArithmeticSubgroupElement'>
        sage: cusps_from_mat(g.matrix())
        (+Infinity, 0)

    Another example::

        sage: from sage.modular.pollack_stevens.space import cusps_from_mat
        sage: g = GammaH(3, [2]).generators()[1].matrix(); g
        [-1  1]
        [-3  2]
        sage: cusps_from_mat(g)
        (1/3, 1/2)
    """
    if isinstance(g, (ArithmeticSubgroupElement, Sigma0Element)):
        g = g.matrix()
    a, b, c, d = g.list()
    if c:
        ac = a / c
    else:
        ac = oo
    if d:
        bd = b / d
    else:
        bd = oo
    return ac, bd


def ps_modsym_from_elliptic_curve(E, sign=0, implementation='eclib'):
    r"""
    Return the overconvergent modular symbol associated to
    an elliptic curve defined over the rationals.

    INPUT:

    - ``E`` -- an elliptic curve defined over the rationals

    - ``sign`` -- the sign (default: 0). If nonzero, returns either
      the plus (if ``sign`` == 1) or the minus (if ``sign`` == -1) modular
      symbol. The default of 0 returns the sum of the plus and minus symbols.

    - ``implementation`` -- either ``'eclib'`` (default) or ``'sage'``. This
      determines which implementation of the underlying classical
      modular symbols is used.

    OUTPUT: the overconvergent modular symbol associated to ``E``

    EXAMPLES::

        sage: E = EllipticCurve('113a1')
        sage: symb = E.pollack_stevens_modular_symbol() # indirect doctest
        sage: symb
        Modular symbol of level 113 with values in Sym^0 Q^2
        sage: symb.values()
        [-1/2, 1, -1, 0, 0, 1, 1, -1, 0, -1, 0, 0, 0, 1, -1, 0, 0, 0, 1, 0, 0]

        sage: E = EllipticCurve([0,1])
        sage: symb = E.pollack_stevens_modular_symbol()
        sage: symb.values()
        [-1/6, 1/3, 1/2, 1/6, -1/6, 1/3, -1/3, -1/2, -1/6, 1/6, 0, -1/6, -1/6]
    """
    if E.base_ring() is not QQ:
        raise ValueError("The elliptic curve must be defined over the "
                         "rationals.")
    sign = Integer(sign)
    if sign not in [0, 1, -1]:
        raise ValueError("The sign must be either 0, 1 or -1")
    N = E.conductor()
    V = PollackStevensModularSymbols(Gamma0(N), 0)
    D = V.coefficient_module()
    manin = V.source()
    # if sage's modular symbols are used we take the
    # normalization given by 'L_ratio' in modular_symbol
    if sign <= 0:
        minus_sym = E.modular_symbol(sign=-1, implementation=implementation)
    if sign >= 0:
        plus_sym = E.modular_symbol(sign=1, implementation=implementation)
    val = {}
    for g in manin.gens():
        ac, bd = cusps_from_mat(g)
        val[g] = D(0)
        if sign >= 0:
            val[g] += D(plus_sym(ac) - plus_sym(bd))
        if sign <= 0:
            val[g] += D(minus_sym(ac) - minus_sym(bd))
    return V(val)


def ps_modsym_from_simple_modsym_space(A, name='alpha'):
    r"""
    Return some choice -- only well defined up a nonzero scalar (!) -- of an
    overconvergent modular symbol that corresponds to ``A``.

    INPUT:

    - ``A`` -- nonzero simple Hecke equivariant new space of modular symbols,
      which need not be cuspidal

    OUTPUT:

    A choice of corresponding overconvergent modular symbols; when dim(A)>1,
    we make an arbitrary choice of defining polynomial for the codomain field.

    EXAMPLES:

    The level 11 example::

        sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
        sage: A = ModularSymbols(11, sign=1, weight=2).decomposition()[0]
        sage: A.is_cuspidal()
        True
        sage: f = ps_modsym_from_simple_modsym_space(A); f
        Modular symbol of level 11 with values in Sym^0 Q^2
        sage: f.values()
        [1, -5/2, -5/2]
        sage: f.weight()         # this is A.weight()-2  !!!!!!
        0

    And the -1 sign for the level 11 example::

        sage: A = ModularSymbols(11, sign=-1, weight=2).decomposition()[0]
        sage: f = ps_modsym_from_simple_modsym_space(A); f.values()
        [0, 1, -1]

    A does not have to be cuspidal; it can be Eisenstein::

        sage: A = ModularSymbols(11, sign=1, weight=2).decomposition()[1]
        sage: A.is_cuspidal()
        False
        sage: f = ps_modsym_from_simple_modsym_space(A); f
        Modular symbol of level 11 with values in Sym^0 Q^2
        sage: f.values()
        [1, 0, 0]

    We create the simplest weight 2 example in which ``A`` has dimension
    bigger than 1::

        sage: A = ModularSymbols(23, sign=1, weight=2).decomposition()[0]
        sage: f = ps_modsym_from_simple_modsym_space(A); f.values()
        [1, 0, 0, 0, 0]
        sage: A = ModularSymbols(23, sign=-1, weight=2).decomposition()[0]
        sage: f = ps_modsym_from_simple_modsym_space(A); f.values()
        [0, 1, -alpha, alpha, -1]
        sage: f.base_ring()
        Number Field in alpha with defining polynomial x^2 + x - 1

    We create the +1 modular symbol attached to the weight 12 modular form ``Delta``::

        sage: A = ModularSymbols(1, sign=+1, weight=12).decomposition()[0]
        sage: f = ps_modsym_from_simple_modsym_space(A); f
        Modular symbol of level 1 with values in Sym^10 Q^2
        sage: f.values()
        [(-1620/691, 0, 1, 0, -9/14, 0, 9/14, 0, -1, 0, 1620/691), (1620/691, 1620/691, 929/691, -453/691, -29145/9674, -42965/9674, -2526/691, -453/691, 1620/691, 1620/691, 0), (0, -1620/691, -1620/691, 453/691, 2526/691, 42965/9674, 29145/9674, 453/691, -929/691, -1620/691, -1620/691)]

    And, the -1 modular symbol attached to ``Delta``::

        sage: A = ModularSymbols(1, sign=-1, weight=12).decomposition()[0]
        sage: f = ps_modsym_from_simple_modsym_space(A); f
        Modular symbol of level 1 with values in Sym^10 Q^2
        sage: f.values()
        [(0, 1, 0, -25/48, 0, 5/12, 0, -25/48, 0, 1, 0), (0, -1, -2, -119/48, -23/12, -5/24, 23/12, 3, 2, 0, 0), (0, 0, 2, 3, 23/12, -5/24, -23/12, -119/48, -2, -1, 0)]

    A consistency check with :meth:`sage.modular.pollack_stevens.space.ps_modsym_from_simple_modsym_space`::

        sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
        sage: E = EllipticCurve('11a')
        sage: f_E = E.pollack_stevens_modular_symbol(); f_E.values()
        [-1/5, 1, 0]
        sage: A = ModularSymbols(11, sign=1, weight=2).decomposition()[0]
        sage: f_plus = ps_modsym_from_simple_modsym_space(A); f_plus.values()
        [1, -5/2, -5/2]
        sage: A = ModularSymbols(11, sign=-1, weight=2).decomposition()[0]
        sage: f_minus = ps_modsym_from_simple_modsym_space(A); f_minus.values()
        [0, 1, -1]

    We find that a linear combination of the plus and minus parts equals the
    Pollack-Stevens symbol attached to ``E``. This illustrates how
    ``ps_modsym_from_simple_modsym_space`` is only well-defined up to a nonzero
    scalar::

        sage: (-1/5)*vector(QQ, f_plus.values()) + (1/2)*vector(QQ, f_minus.values())
        (-1/5, 1, 0)
        sage: vector(QQ, f_E.values())
        (-1/5, 1, 0)

    The next few examples all illustrate the ways in which exceptions are
    raised if A does not satisfy various constraints.

    First, ``A`` must be new::

        sage: A = ModularSymbols(33,sign=1).cuspidal_subspace().old_subspace()
        sage: ps_modsym_from_simple_modsym_space(A)
        Traceback (most recent call last):
        ...
        ValueError: A must be new

    ``A`` must be simple::

        sage: A = ModularSymbols(43,sign=1).cuspidal_subspace()
        sage: ps_modsym_from_simple_modsym_space(A)
        Traceback (most recent call last):
        ...
        ValueError: A must be simple

    ``A`` must have sign -1 or +1 in order to be simple::

        sage: A = ModularSymbols(11).cuspidal_subspace()
        sage: ps_modsym_from_simple_modsym_space(A)
        Traceback (most recent call last):
        ...
        ValueError: A must have sign +1 or -1 (otherwise it is not simple)

    The dimension must be positive::

        sage: A = ModularSymbols(10).cuspidal_subspace(); A
        Modular Symbols subspace of dimension 0 of Modular Symbols space of dimension 3 for Gamma_0(10) of weight 2 with sign 0 over Rational Field
        sage: ps_modsym_from_simple_modsym_space(A)
        Traceback (most recent call last):
        ...
        ValueError: A must have positive dimension

    We check that forms of nontrivial character are getting handled correctly::

        sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
        sage: f = Newforms(Gamma1(13), names='a')[0]
        sage: phi = ps_modsym_from_simple_modsym_space(f.modular_symbols(1))
        sage: phi.hecke(7)
        Modular symbol of level 13 with values in Sym^0 (Number Field in alpha with defining polynomial x^2 + 3*x + 3)^2 twisted by Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -alpha - 1
        sage: phi.hecke(7).values()
        [0, 0, 0, 0, 0]
    """
    if A.dimension() == 0:
        raise ValueError("A must have positive dimension")

    if A.sign() == 0:
        raise ValueError("A must have sign +1 or -1 (otherwise it is"
                         " not simple)")

    if not A.is_new():
        raise ValueError("A must be new")

    if not A.is_simple():
        raise ValueError("A must be simple")

    M = A.ambient_module()
    w = A.dual_eigenvector(name)
    K = w.base_ring()
    chi = A.q_eigenform_character(name)
    V = PollackStevensModularSymbols(chi, A.weight() - 2, base_ring=K, sign=A.sign())
    D = V.coefficient_module()
    # N = V.level()
    k = V.weight()  # = A.weight() - 2
    manin = V.source()
    val = {}
    for g in manin.gens():
        ac, bd = cusps_from_mat(g)
        v = []
        for j in range(k + 1):
            # TODO: The following might be backward: it should be the coefficient of X^j Y^(k-j)
            v.append(w.dot_product(M.modular_symbol([j, ac, bd]).element()) * (-1) ** (k - j))
        val[g] = D(v)
    return V(val)
