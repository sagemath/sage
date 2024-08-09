r"""
Congruence subgroup `\Gamma(N)`
"""

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.arith.misc import gcd
from sage.groups.matrix_gps.finitely_generated import MatrixGroup
from sage.matrix.constructor import matrix
from sage.misc.misc_c import prod
from sage.modular.cusps import Cusp
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.integer import GCD_list
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.richcmp import richcmp_method, richcmp

from .congroup_generic import CongruenceSubgroup
from .congroup_sl2z import SL2Z


_gamma_cache = {}
def Gamma_constructor(N):
    r"""
    Return the congruence subgroup `\Gamma(N)`.

    EXAMPLES::

        sage: Gamma(5) # indirect doctest
        Congruence Subgroup Gamma(5)
        sage: G = Gamma(23)
        sage: G is Gamma(23)
        True
        sage: TestSuite(G).run()

    Test global uniqueness::

        sage: G = Gamma(17)
        sage: G is loads(dumps(G))
        True
        sage: G2 = sage.modular.arithgroup.congroup_gamma.Gamma_class(17)
        sage: G == G2
        True
        sage: G is G2
        False
    """
    if N == 1:
        return SL2Z
    try:
        return _gamma_cache[N]
    except KeyError:
        _gamma_cache[N] = Gamma_class(N)
        return _gamma_cache[N]


@richcmp_method
class Gamma_class(CongruenceSubgroup):
    r"""
    The principal congruence subgroup `\Gamma(N)`.
    """
    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: Gamma(133)._repr_()
            'Congruence Subgroup Gamma(133)'
        """
        return "Congruence Subgroup Gamma(%s)" % self.level()

    def _latex_(self):
        r"""
        Return the \LaTeX representation of ``self``.

        EXAMPLES::

            sage: Gamma(20)._latex_()
            '\\Gamma(20)'
            sage: latex(Gamma(20))
            \Gamma(20)
        """
        return "\\Gamma(%s)" % self.level()

    def __reduce__(self):
        """
        Used for pickling ``self``.

        EXAMPLES::

            sage: Gamma(5).__reduce__()
            (<function Gamma_constructor at ...>, (5,))
        """
        return Gamma_constructor, (self.level(),)

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: Gamma(3) == SymmetricGroup(8)                                         # needs sage.groups
            False
            sage: Gamma(3) == Gamma1(3)
            False
            sage: Gamma(5) < Gamma(6)
            True
            sage: Gamma(5) == Gamma(5)
            True
            sage: Gamma(3) == Gamma(3).as_permutation_group()                           # needs sage.groups
            True
        """
        if isinstance(other, Gamma_class):
            return richcmp(self.level(), other.level(), op)
        else:
            return NotImplemented

    def index(self):
        r"""
        Return the index of ``self`` in the full modular group. This is given by

        .. MATH::

          \prod_{\substack{p \mid N \\ \text{$p$ prime}}}\left(p^{3e}-p^{3e-2}\right).

        EXAMPLES::

            sage: [Gamma(n).index() for n in [1..19]]
            [1, 6, 24, 48, 120, 144, 336, 384, 648, 720, 1320, 1152, 2184, 2016, 2880, 3072, 4896, 3888, 6840]
            sage: Gamma(32041).index()
            32893086819240
        """
        return prod([p**(3*e-2)*(p*p-1) for (p,e) in self.level().factor()])

    def _contains_sl2(self, a,b,c,d):
        r"""
        EXAMPLES::

            sage: G = Gamma(5)
            sage: [1, 0, -10, 1] in G
            True
            sage: 1 in G
            True
            sage: SL2Z([26, 5, 5, 1]) in G
            True
            sage: SL2Z([1, 1, 6, 7]) in G
            False
        """
        N = self.level()
        # don't need to check d == 1 as this is automatic from det
        return ((a % N == 1) and (b % N == 0) and (c % N == 0))

    def ncusps(self):
        r"""
        Return the number of cusps of this subgroup `\Gamma(N)`.

        EXAMPLES::

            sage: [Gamma(n).ncusps() for n in [1..19]]
            [1, 3, 4, 6, 12, 12, 24, 24, 36, 36, 60, 48, 84, 72, 96, 96, 144, 108, 180]
            sage: Gamma(30030).ncusps()
            278691840
            sage: Gamma(2^30).ncusps()
            432345564227567616
        """
        n = self.level()
        if n == 1:
            return ZZ(1)
        if n == 2:
            return ZZ(3)
        return prod([p**(2*e) - p**(2*e-2) for (p,e) in n.factor()])//2

    def nirregcusps(self):
        r"""
        Return the number of irregular cusps of ``self``. For principal
        congruence subgroups this is always 0.

        EXAMPLES::

            sage: Gamma(17).nirregcusps()
            0
        """
        return 0

    def _find_cusps(self):
        r"""
        Calculate the reduced representatives of the equivalence classes of
        cusps for this group. Adapted from code by Ron Evans.

        EXAMPLES::

            sage: Gamma(8).cusps() # indirect doctest
            [0, 1/4, 1/3, 3/8, 1/2, 2/3, 3/4, 1, 4/3, 3/2, 5/3, 2, 7/3, 5/2, 8/3, 3, 7/2, 11/3, 4, 14/3, 5, 6, 7, Infinity]
        """
        n = self.level()
        C = [QQ(x) for x in range(n)]

        n0 = n//2
        n1 = (n+1)//2

        for r in range(1, n1):
            if r > 1 and gcd(r,n) == 1:
                C.append(ZZ(r)/ZZ(n))
            if n0 == n/2 and gcd(r,n0) == 1:
                C.append(ZZ(r)/ZZ(n0))

        for s in range(2,n1):
            for r in range(1, 1+n):
                if GCD_list([s,r,n]) == 1:
                    # GCD_list is ~40x faster than gcd, since gcd wastes loads
                    # of time initialising a Sequence type.
                    u,v = _lift_pair(r,s,n)
                    C.append(ZZ(u)/ZZ(v))

        return [Cusp(x) for x in sorted(C)] + [Cusp(1,0)]

    def reduce_cusp(self, c):
        r"""
        Calculate the unique reduced representative of the equivalence of the
        cusp `c` modulo this group. The reduced representative of an
        equivalence class is the unique cusp in the class of the form `u/v`
        with `u, v \ge 0` coprime, `v` minimal, and `u` minimal for that `v`.

        EXAMPLES::

            sage: Gamma(5).reduce_cusp(1/5)
            Infinity
            sage: Gamma(5).reduce_cusp(7/8)
            3/2
            sage: Gamma(6).reduce_cusp(4/3)
            2/3

        TESTS::

            sage: G = Gamma(50)
            sage: all(c == G.reduce_cusp(c) for c in G.cusps())
            True

        We test that :issue:`36163` is fixed::

            sage: Gamma(7).reduce_cusp(Cusp(6,7))
            Infinity
        """
        N = self.level()
        c = Cusp(c)
        u,v = c.numerator() % N, c.denominator() % N
        if (v > N//2) or (2*v == N and u > N//2):
            u,v = -u,-v
        u,v = _lift_pair(u,v,N)
        return Cusp(u,v)

    def are_equivalent(self, x, y, trans=False):
        r"""
        Check if the cusps `x` and `y` are equivalent under the action of this group.

        ALGORITHM: The cusps `u_1 / v_1` and `u_2 / v_2` are equivalent modulo
        `\Gamma(N)` if and only if `(u_1, v_1) = \pm (u_2, v_2) \bmod N`.

        EXAMPLES::

            sage: Gamma(7).are_equivalent(Cusp(2/3), Cusp(5/4))
            True
        """
        if trans:
            return CongruenceSubgroup.are_equivalent(self, x,y,trans=trans)
        N = self.level()
        u1,v1 = (x.numerator() % N, x.denominator() % N)
        u2,v2 = (y.numerator(), y.denominator())

        return ((u1,v1) == (u2 % N, v2 % N)) or ((u1,v1) == (-u2 % N, -v2 % N))

    def nu3(self):
        r"""
        Return the number of elliptic points of order 3 for this arithmetic
        subgroup. Since this subgroup is `\Gamma(N)` for `N \ge 2`, there are
        no such points, so we return 0.

        EXAMPLES::

            sage: Gamma(89).nu3()
            0
        """
        return 0

    # We don't need to override nu2, since the default nu2 implementation knows
    # that nu2 = 0 for odd subgroups.

    def image_mod_n(self):
        r"""
        Return the image of this group modulo `N`, as a subgroup of `SL(2, \ZZ
        / N\ZZ)`. This is just the trivial subgroup.

        EXAMPLES::

            sage: Gamma(3).image_mod_n()
            Matrix group over Ring of integers modulo 3 with 1 generators (
            [1 0]
            [0 1]
            )
        """
        return MatrixGroup([matrix(Zmod(self.level()), 2, 2, 1)])


def is_Gamma(x):
    r"""
    Return ``True`` if x is a congruence subgroup of type Gamma.

    EXAMPLES::

        sage: from sage.modular.arithgroup.all import Gamma_class
        sage: isinstance(Gamma0(13), Gamma_class)
        False
        sage: isinstance(Gamma(4), Gamma_class)
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38035, "The function is_Gamma is deprecated; use 'isinstance(..., Gamma_class)' instead.")
    return isinstance(x, Gamma_class)


def _lift_pair(U,V,N):
    r"""
    Utility function. Given integers ``U, V, N``, with `N \ge 1` and `{\rm
    gcd}(U, V, N) = 1`, return a pair `(u, v)` congruent to `(U, V) \bmod N`,
    such that `{\rm gcd}(u,v) = 1`, `u, v \ge 0`, `v` is as small as possible,
    and `u` is as small as possible for that `v`.

    *Warning*: As this function is for internal use, it does not do a
    preliminary sanity check on its input, for efficiency. It will recover
    reasonably gracefully if ``(U, V, N)`` are not coprime, but only after
    wasting quite a lot of cycles!

    EXAMPLES::

        sage: from sage.modular.arithgroup.congroup_gamma import _lift_pair
        sage: _lift_pair(2,4,7)
        (9, 4)
        sage: _lift_pair(2,4,8) # don't do this
        Traceback (most recent call last):
        ...
        ValueError: (U, V, N) must be coprime
    """
    u = U % N
    v = V % N
    if v == 0:
        if u == 1 or u == N-1:
            return (1,0)
        else:
            v = N
    while gcd(u, v) > 1:
        u = u+N
        if u > N*v:
            raise ValueError("(U, V, N) must be coprime")
    return (u, v)
