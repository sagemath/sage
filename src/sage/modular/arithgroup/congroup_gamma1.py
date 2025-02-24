# sage.doctest: needs sage.libs.pari
r"""
Congruence subgroup `\Gamma_1(N)`
"""

# ****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.misc.cachefunc import cached_method

from sage.misc.misc_c import prod
from .congroup_gammaH import GammaH_class, is_GammaH, GammaH_constructor
from sage.rings.integer_ring import ZZ
from sage.arith.misc import euler_phi as phi, moebius, divisors
from sage.modular.dirichlet import DirichletGroup


def is_Gamma1(x):
    """
    Return ``True`` if x is a congruence subgroup of type Gamma1.

    EXAMPLES::

        sage: from sage.modular.arithgroup.all import is_Gamma1
        sage: is_Gamma1(SL2Z)
        doctest:warning...
        DeprecationWarning: The function is_Gamma1 is deprecated; use 'isinstance(..., Gamma1_class)' instead.
        See https://github.com/sagemath/sage/issues/38035 for details.
        False
        sage: is_Gamma1(Gamma1(13))
        True
        sage: is_Gamma1(Gamma0(6))
        False
        sage: is_Gamma1(GammaH(12, [])) # trick question!
        True
        sage: is_Gamma1(GammaH(12, [5]))
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(38035, "The function is_Gamma1 is deprecated; use 'isinstance(..., Gamma1_class)' instead.")
    # from congroup_sl2z import is_SL2Z
    # return (isinstance(x, Gamma1_class) or is_SL2Z(x))
    return isinstance(x, Gamma1_class)


_gamma1_cache = {}


def Gamma1_constructor(N):
    r"""
    Return the congruence subgroup `\Gamma_1(N)`.

    EXAMPLES::

        sage: Gamma1(5) # indirect doctest
        Congruence Subgroup Gamma1(5)
        sage: G = Gamma1(23)
        sage: G is Gamma1(23)
        True
        sage: G is GammaH(23, [1])
        True
        sage: TestSuite(G).run()
        sage: G is loads(dumps(G))
        True
    """
    if N == 1 or N == 2:
        from .congroup_gamma0 import Gamma0_constructor
        return Gamma0_constructor(N)
    try:
        return _gamma1_cache[N]
    except KeyError:
        _gamma1_cache[N] = Gamma1_class(N)
        return _gamma1_cache[N]


class Gamma1_class(GammaH_class):
    r"""
    The congruence subgroup `\Gamma_1(N)`.

    TESTS::

        sage: [Gamma1(n).genus() for n in prime_range(2,100)]
        [0, 0, 0, 0, 1, 2, 5, 7, 12, 22, 26, 40, 51, 57, 70, 92, 117, 126, 155, 176, 187, 222, 247, 287, 345]
        sage: [Gamma1(n).index() for n in [1..10]]
        [1, 3, 8, 12, 24, 24, 48, 48, 72, 72]

        sage: [Gamma1(n).dimension_cusp_forms() for n in [1..20]]
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 1, 1, 2, 5, 2, 7, 3]
        sage: [Gamma1(n).dimension_cusp_forms(1) for n in [1..20]]
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: [Gamma1(4).dimension_cusp_forms(k) for k in [1..20]]
        [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8]

        sage: Gamma1(23).dimension_cusp_forms(1)
        1
    """

    def __init__(self, level):
        r"""
        The congruence subgroup `\Gamma_1(N)`.

        EXAMPLES::

            sage: G = Gamma1(11); G
            Congruence Subgroup Gamma1(11)
            sage: loads(G.dumps()) == G
            True
        """
        GammaH_class.__init__(self, level, [])

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: Gamma1(133)._repr_()
            'Congruence Subgroup Gamma1(133)'
        """
        return "Congruence Subgroup Gamma1(%s)" % self.level()

    def __reduce__(self):
        """
        Used for pickling ``self``.

        EXAMPLES::

            sage: Gamma1(82).__reduce__()
            (<function Gamma1_constructor at ...>, (82,))
        """
        return Gamma1_constructor, (self.level(),)

    def _latex_(self):
        r"""
        Return the \LaTeX representation of ``self``.

        EXAMPLES::

            sage: Gamma1(3)._latex_()
            '\\Gamma_1(3)'
            sage: latex(Gamma1(3))
            \Gamma_1(3)
        """
        return "\\Gamma_1(%s)" % self.level()

    def is_even(self):
        """
        Return ``True`` precisely if this subgroup contains the matrix -1.

        EXAMPLES::

            sage: Gamma1(1).is_even()
            True
            sage: Gamma1(2).is_even()
            True
            sage: Gamma1(15).is_even()
            False
        """
        return self.level() in [1,2]

    def is_subgroup(self, right):
        """
        Return ``True`` if ``self`` is a subgroup of ``right``.

        EXAMPLES::

            sage: Gamma1(3).is_subgroup(SL2Z)
            True
            sage: Gamma1(3).is_subgroup(Gamma1(5))
            False
            sage: Gamma1(3).is_subgroup(Gamma1(6))
            False
            sage: Gamma1(6).is_subgroup(Gamma1(3))
            True
            sage: Gamma1(6).is_subgroup(Gamma0(2))
            True
            sage: Gamma1(80).is_subgroup(GammaH(40, []))
            True
            sage: Gamma1(80).is_subgroup(GammaH(40, [21]))
            True
        """
        if right.level() == 1:
            return True
        if isinstance(right, GammaH_class):
            return self.level() % right.level() == 0
        else:
            raise NotImplementedError

    @cached_method
    def generators(self, algorithm='farey'):
        r"""
        Return generators for this congruence subgroup. The result is cached.

        INPUT:

        - ``algorithm`` -- string; either ``'farey'`` (default) or
          ``'todd-coxeter'``

        If ``algorithm`` is set to ``'farey'``, then the generators will be
        calculated using Farey symbols, which will always return a *minimal*
        generating set. See :mod:`~sage.modular.arithgroup.farey_symbol` for
        more information.

        If ``algorithm`` is set to ``'todd-coxeter'``, a simpler algorithm
        based on Todd-Coxeter enumeration will be used. This tends to return
        far larger sets of generators.

        EXAMPLES::

            sage: Gamma1(3).generators()
            [
            [1 1]  [ 1 -1]
            [0 1], [ 3 -2]
            ]
            sage: Gamma1(3).generators(algorithm='todd-coxeter')
            [
            [1 1]  [-2  1]  [1 1]  [ 1 -1]  [1 0]  [1 1]  [-5  2]  [ 1  0]
            [0 1], [-3  1], [0 1], [ 0  1], [3 1], [0 1], [12 -5], [-3  1],
            <BLANKLINE>
            [ 1 -1]  [ 1 -1]  [ 4 -1]  [ -5   3]
            [ 3 -2], [ 3 -2], [ 9 -2], [-12   7]
            ]
        """
        if algorithm == "farey":
            return self.farey_symbol().generators()
        elif algorithm == "todd-coxeter":
            from sage.modular.modsym.g1list import G1list
            from .congroup import generators_helper
            level = self.level()
            gen_list = generators_helper(G1list(level), level)
            return [self(g, check=False) for g in gen_list]
        else:
            raise ValueError("Unknown algorithm '%s' (should be either 'farey' or 'todd-coxeter')" % algorithm)

    def _contains_sl2(self, a, b, c, d):
        r"""
        Test whether x is an element of this group.

        EXAMPLES::

            sage: G = Gamma1(5)
            sage: [1, 0, -10, 1] in G
            True
            sage: matrix(ZZ, 2, [6, 1, 5, 1]) in G
            True
            sage: SL2Z.0 in G
            False
            sage: G([1, 1, 6, 7]) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: matrix [1 1]
            [6 7] is not an element of Congruence Subgroup Gamma1(5)
        """
        N = self.level()
        # don't need to check d == 1 mod N as this is automatic from det
        return ((a % N == 1) and (c % N == 0))

    def nu2(self):
        r"""
        Calculate the number of orbits of elliptic points of order 2 for this
        subgroup `\Gamma_1(N)`. This is known to be 0 if N > 2.

        EXAMPLES::

            sage: Gamma1(2).nu2()
            1
            sage: Gamma1(457).nu2()
            0
            sage: [Gamma1(n).nu2() for n in [1..16]]
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        """
        N = self.level()
        if N > 2:
            return 0
        elif N == 2 or N == 1:
            return 1

    def nu3(self):
        r"""
        Calculate the number of orbits of elliptic points of order 3 for this
        subgroup `\Gamma_1(N)`. This is known to be 0 if N > 3.

        EXAMPLES::

            sage: Gamma1(2).nu3()
            0
            sage: Gamma1(3).nu3()
            1
            sage: Gamma1(457).nu3()
            0
            sage: [Gamma1(n).nu3() for n in [1..10]]
            [1, 0, 1, 0, 0, 0, 0, 0, 0, 0]
        """
        N = self.level()
        if N > 3 or N == 2:
            return 0
        else:
            return 1

    def ncusps(self):
        r"""
        Return the number of cusps of this subgroup `\Gamma_1(N)`.

        EXAMPLES::

            sage: [Gamma1(n).ncusps() for n in [1..15]]
            [1, 2, 2, 3, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12, 16]
            sage: [Gamma1(n).ncusps() for n in prime_range(2, 100)]
            [2, 2, 4, 6, 10, 12, 16, 18, 22, 28, 30, 36, 40, 42, 46, 52, 58, 60, 66, 70, 72, 78, 82, 88, 96]
        """
        n = self.level()
        if n <= 4:
            return [None, 1, 2, 2, 3][n]
        return ZZ(sum([phi(d)*phi(n/d)/ZZ(2) for d in n.divisors()]))

    def index(self):
        r"""
        Return the index of ``self`` in the full modular group. This is given by the formula

        .. MATH::

            N^2 \prod_{\substack{p \mid N \\ \text{$p$ prime}}} \left( 1 - \frac{1}{p^2}\right).

        EXAMPLES::

            sage: Gamma1(180).index()
            20736
            sage: [Gamma1(n).projective_index() for n in [1..16]]
            [1, 3, 4, 6, 12, 12, 24, 24, 36, 36, 60, 48, 84, 72, 96, 96]
        """
        return prod([p**(2*e) - p**(2*e-2) for (p,e) in self.level().factor()])

    ##################################################################################
    # Dimension formulas for Gamma1, accepting a Dirichlet character as an argument. #
    ##################################################################################

    def dimension_modular_forms(self, k=2, eps=None, algorithm='CohenOesterle'):
        r"""
        Return the dimension of the space of modular forms for ``self``, or the
        dimension of the subspace corresponding to the given character if one
        is supplied.

        INPUT:

        - ``k`` -- integer (default: 2); the weight

        - ``eps`` -- either ``None`` or a Dirichlet character modulo N, where N is
          the level of this group. If this is ``None``, then the dimension of the
          whole space is returned; otherwise, the dimension of the subspace of
          forms of character eps.

        - ``algorithm`` -- either ``'CohenOesterle'`` (the default) or
          ``'Quer'``. This specifies the method to use in the case of
          nontrivial character: either the Cohen--Oesterle formula as described
          in Stein's book, or by Möbius inversion using the subgroups GammaH (a
          method due to Jordi Quer).

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = CyclotomicField(3)
            sage: eps = DirichletGroup(7*43,K).0^2
            sage: G = Gamma1(7*43)
            sage: G.dimension_modular_forms(2, eps)
            32
            sage: G.dimension_modular_forms(2, eps, algorithm='Quer')
            32

        TESTS:

        Check that :issue:`18436` is fixed::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 + x + 1)
            sage: G = DirichletGroup(13, base_ring=K)
            sage: Gamma1(13).dimension_modular_forms(2, G[1])
            3
            sage: Gamma1(13).dimension_modular_forms(2, G[1], algorithm='Quer')
            3
            sage: Gamma1(39).dimension_modular_forms(2, G[1])
            7
            sage: Gamma1(39).dimension_modular_forms(2, G[1], algorithm='Quer')
            7
        """
        return self.dimension_cusp_forms(k, eps, algorithm) + self.dimension_eis(k, eps, algorithm)

    def dimension_cusp_forms(self, k=2, eps=None, algorithm='CohenOesterle'):
        r"""
        Return the dimension of the space of cusp forms for ``self``, or the
        dimension of the subspace corresponding to the given character if one
        is supplied.

        INPUT:

        - ``k`` -- integer (default: 2); the weight

        - ``eps`` -- either ``None`` or a Dirichlet character modulo N, where N is
          the level of this group. If this is ``None``, then the dimension of the
          whole space is returned; otherwise, the dimension of the subspace of
          forms of character eps.

        - ``algorithm`` -- either ``'CohenOesterle'`` (the default) or
          ``'Quer'``. This specifies the method to use in the case of nontrivial
          character: either the Cohen--Oesterle formula as described in Stein's
          book, or by Möbius inversion using the subgroups GammaH (a method due
          to Jordi Quer). Ignored for weight 1.

        EXAMPLES:

        We compute the same dimension in two different ways ::

            sage: # needs sage.rings.number_field
            sage: K = CyclotomicField(3)
            sage: eps = DirichletGroup(7*43,K).0^2
            sage: G = Gamma1(7*43)

        Via Cohen--Oesterle::

            sage: Gamma1(7*43).dimension_cusp_forms(2, eps)                             # needs sage.rings.number_field
            28

        Via Quer's method::

            sage: Gamma1(7*43).dimension_cusp_forms(2, eps, algorithm='Quer')           # needs sage.rings.number_field
            28

        Some more examples::

            sage: G.<eps> = DirichletGroup(9)
            sage: [Gamma1(9).dimension_cusp_forms(k, eps) for k in [1..10]]
            [0, 0, 1, 0, 3, 0, 5, 0, 7, 0]
            sage: [Gamma1(9).dimension_cusp_forms(k, eps^2) for k in [1..10]]
            [0, 0, 0, 2, 0, 4, 0, 6, 0, 8]

        In weight 1, we can sometimes rule out cusp forms existing via
        Riemann-Roch, but if this does not work, we trigger computation of the
        cusp forms space via Schaeffer's algorithm::

            sage: chi = [u for u in DirichletGroup(40) if u(-1) == -1 and u(21) == 1][0]
            sage: Gamma1(40).dimension_cusp_forms(1, chi)
            0
            sage: G = DirichletGroup(57); chi = (G.0) * (G.1)^6
            sage: Gamma1(57).dimension_cusp_forms(1, chi)
            1
        """
        from .all import Gamma0

        # first deal with special cases

        if eps is None:
            return GammaH_class.dimension_cusp_forms(self, k)

        N = self.level()
        K = eps.base_ring()
        eps = DirichletGroup(N, K)(eps)

        if K.characteristic() != 0:
            raise NotImplementedError('dimension_cusp_forms() is only implemented for rings of characteristic 0')

        if eps.is_trivial():
            return Gamma0(N).dimension_cusp_forms(k)

        if (k <= 0) or ((k % 2) == 1 and eps.is_even()) or ((k % 2) == 0 and eps.is_odd()):
            return ZZ(0)

        if k == 1:
            from sage.modular.modform.weight1 import dimension_wt1_cusp_forms
            return dimension_wt1_cusp_forms(eps)

        # now the main part

        if algorithm == "Quer":
            n = eps.order()
            dim = ZZ(0)
            for d in n.divisors():
                G = GammaH_constructor(N,(eps**d).kernel())
                dim = dim + moebius(d)*G.dimension_cusp_forms(k)
            return dim//phi(n)

        elif algorithm == "CohenOesterle":
            from sage.modular.dims import CohenOesterle
            return ZZ( K(Gamma0(N).index() * (k-1)/ZZ(12)) + CohenOesterle(eps,k) )

        else:  # algorithm not in ["CohenOesterle", "Quer"]:
            raise ValueError("Unrecognised algorithm in dimension_cusp_forms")

    def dimension_eis(self, k=2, eps=None, algorithm='CohenOesterle'):
        r"""
        Return the dimension of the space of Eisenstein series forms for self,
        or the dimension of the subspace corresponding to the given character
        if one is supplied.

        INPUT:

        - ``k`` -- integer (default: 2); the weight

        - ``eps`` -- either ``None`` or a Dirichlet character modulo N, where N is
          the level of this group. If this is ``None``, then the dimension of the
          whole space is returned; otherwise, the dimension of the subspace of
          Eisenstein series of character eps.

        - ``algorithm`` -- either ``'CohenOesterle'`` (the default) or
          ``'Quer'``. This specifies the method to use in the case of nontrivial
          character: either the Cohen--Oesterle formula as described in Stein's
          book, or by Möbius inversion using the subgroups GammaH (a method due
          to Jordi Quer).

        AUTHORS:

        - William Stein - Cohen--Oesterle algorithm

        - Jordi Quer - algorithm based on GammaH subgroups

        - David Loeffler (2009) - code refactoring

        EXAMPLES:

        The following two computations use different algorithms::

            sage: [Gamma1(36).dimension_eis(1,eps) for eps in DirichletGroup(36)]
            [0, 4, 3, 0, 0, 2, 6, 0, 0, 2, 3, 0]
            sage: [Gamma1(36).dimension_eis(1,eps,algorithm='Quer') for eps in DirichletGroup(36)]
            [0, 4, 3, 0, 0, 2, 6, 0, 0, 2, 3, 0]

        So do these::

            sage: [Gamma1(48).dimension_eis(3,eps) for eps in DirichletGroup(48)]
            [0, 12, 0, 4, 0, 8, 0, 4, 12, 0, 4, 0, 8, 0, 4, 0]
            sage: [Gamma1(48).dimension_eis(3,eps,algorithm='Quer') for eps in DirichletGroup(48)]
            [0, 12, 0, 4, 0, 8, 0, 4, 12, 0, 4, 0, 8, 0, 4, 0]
        """
        from .all import Gamma0

        # first deal with special cases

        if eps is None:
            return GammaH_class.dimension_eis(self, k)

        N = self.level()
        K = eps.base_ring()
        eps = DirichletGroup(N, K)(eps)

        if eps.is_trivial():
            return Gamma0(N).dimension_eis(k)

        # Note case of k = 0 and trivial character already dealt with separately, so k <= 0 here is valid:
        if (k <= 0) or ((k % 2) == 1 and eps.is_even()) or ((k % 2) == 0 and eps.is_odd()):
            return ZZ(0)

        if algorithm == "Quer":
            n = eps.order()
            dim = ZZ(0)
            for d in n.divisors():
                G = GammaH_constructor(N,(eps**d).kernel())
                dim = dim + moebius(d)*G.dimension_eis(k)
            return dim//phi(n)

        elif algorithm == "CohenOesterle":
            from sage.modular.dims import CohenOesterle
            j = 2-k
            # We use the Cohen-Oesterle formula in a subtle way to
            # compute dim M_k(N,eps) (see Ch. 6 of William Stein's book on
            # computing with modular forms).
            alpha = -ZZ( K(Gamma0(N).index()*(j-1)/ZZ(12)) + CohenOesterle(eps,j) )
            if k == 1:
                return alpha
            else:
                return alpha - self.dimension_cusp_forms(k, eps)

        else:  # algorithm not in ["CohenOesterle", "Quer"]:
            raise ValueError("Unrecognised algorithm in dimension_eis")

    def dimension_new_cusp_forms(self, k=2, eps=None, p=0, algorithm='CohenOesterle'):
        r"""
        Dimension of the new subspace (or `p`-new subspace) of cusp forms of
        weight `k` and character `\varepsilon`.

        INPUT:

        - ``k`` -- integer (default: 2)

        - ``eps`` -- a Dirichlet character

        - ``p`` -- a prime (default: 0); just the `p`-new subspace if given

        - ``algorithm`` -- either ``'CohenOesterle'`` (the default) or
          ``'Quer'``. This specifies the method to use in the case of nontrivial
          character: either the Cohen--Oesterle formula as described in Stein's
          book, or by Möbius inversion using the subgroups GammaH (a method due
          to Jordi Quer).

        EXAMPLES::

            sage: G = DirichletGroup(9)
            sage: eps = G.0^3
            sage: eps.conductor()
            3
            sage: [Gamma1(9).dimension_new_cusp_forms(k, eps) for k in [2..10]]
            [0, 0, 0, 2, 0, 2, 0, 2, 0]
            sage: [Gamma1(9).dimension_cusp_forms(k, eps) for k in [2..10]]
            [0, 0, 0, 2, 0, 4, 0, 6, 0]
            sage: [Gamma1(9).dimension_new_cusp_forms(k, eps, 3) for k in [2..10]]
            [0, 0, 0, 2, 0, 2, 0, 2, 0]

        Double check using modular symbols (independent calculation)::

            sage: [ModularSymbols(eps,k,sign=1).cuspidal_subspace().new_subspace().dimension()  for k in [2..10]]
            [0, 0, 0, 2, 0, 2, 0, 2, 0]
            sage: [ModularSymbols(eps,k,sign=1).cuspidal_subspace().new_subspace(3).dimension()  for k in [2..10]]
            [0, 0, 0, 2, 0, 2, 0, 2, 0]

        Another example at level 33::

            sage: G = DirichletGroup(33)
            sage: eps = G.1
            sage: eps.conductor()
            11
            sage: [Gamma1(33).dimension_new_cusp_forms(k, G.1) for k in [2..4]]
            [0, 4, 0]
            sage: [Gamma1(33).dimension_new_cusp_forms(k, G.1, algorithm='Quer') for k in [2..4]]
            [0, 4, 0]
            sage: [Gamma1(33).dimension_new_cusp_forms(k, G.1^2) for k in [2..4]]
            [2, 0, 6]
            sage: [Gamma1(33).dimension_new_cusp_forms(k, G.1^2, p=3) for k in [2..4]]
            [2, 0, 6]
        """

        if eps is None:
            return GammaH_class.dimension_new_cusp_forms(self, k, p)

        N = self.level()
        eps = DirichletGroup(N, eps.base_ring())(eps)

        if eps.is_trivial():
            from .all import Gamma0
            return Gamma0(N).dimension_new_cusp_forms(k, p)

        from .congroup_gammaH import mumu

        if p == 0 or N % p != 0 or eps.conductor().valuation(p) == N.valuation(p):
            D = [eps.conductor()*d for d in divisors(N//eps.conductor())]
            return sum([Gamma1_constructor(M).dimension_cusp_forms(k, eps.restrict(M), algorithm)*mumu(N//M) for M in D])
        eps_p = eps.restrict(N//p)
        old = Gamma1_constructor(N//p).dimension_cusp_forms(k, eps_p, algorithm)
        return self.dimension_cusp_forms(k, eps, algorithm) - 2*old
