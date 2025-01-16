# sage.doctest: needs sage.combinat sage.modules
r"""
Kostka-Foulkes Polynomials

Based on the algorithms in John Stembridge's SF package for Maple
which can be found at http://www.math.lsa.umich.edu/~jrs/maple.html
.
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                     2007 John Stembridge
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

from sage.combinat.partition import _Partitions
from sage.combinat.partitions import ZS1_iterator
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.integer_ring import ZZ


def KostkaFoulkesPolynomial(mu, nu, t=None):
    r"""
    Return the Kostka-Foulkes polynomial `K_{\mu, \nu}(t)`.

    INPUT:

    - ``mu``, ``nu`` -- partitions
    - ``t`` -- an optional parameter (default: ``None``)

    OUTPUT:

    - the Koskta-Foulkes polynomial indexed by partitions ``mu`` and ``nu`` and
      evaluated at the parameter ``t``.  If ``t`` is ``None`` the resulting
      polynomial is in the polynomial ring `\ZZ['t']`.

    EXAMPLES::

        sage: KostkaFoulkesPolynomial([2,2],[2,2])
        1
        sage: KostkaFoulkesPolynomial([2,2],[4])
        0
        sage: KostkaFoulkesPolynomial([2,2],[1,1,1,1])
        t^4 + t^2
        sage: KostkaFoulkesPolynomial([2,2],[2,1,1])
        t
        sage: q = PolynomialRing(QQ,'q').gen()
        sage: KostkaFoulkesPolynomial([2,2],[2,1,1],q)
        q

    TESTS::

        sage: KostkaFoulkesPolynomial([2,4],[2,2])
        Traceback (most recent call last):
        ...
        ValueError: mu must be a partition
        sage: KostkaFoulkesPolynomial([2,2],[2,4])
        Traceback (most recent call last):
        ...
        ValueError: nu must be a partition
        sage: KostkaFoulkesPolynomial([3,2],[2,1])
        Traceback (most recent call last):
        ...
        ValueError: mu and nu must be partitions of the same size
    """
    if mu not in _Partitions:
        raise ValueError("mu must be a partition")
    if nu not in _Partitions:
        raise ValueError("nu must be a partition")
    if sum(mu) != sum(nu):
        raise ValueError("mu and nu must be partitions of the same size")
    return kfpoly(mu, nu, t)


def kfpoly(mu, nu, t=None):
    r"""
    Return the Kostka-Foulkes polynomial `K_{\mu, \nu}(t)`
    by generating all rigging sequences for the shape `\mu`, and then
    selecting those of content `\nu`.

    INPUT:

    - ``mu``, ``nu`` -- partitions
    - ``t`` -- an optional parameter (default: ``None``)

    OUTPUT:

    - the Koskta-Foulkes polynomial indexed by partitions ``mu`` and ``nu`` and
      evaluated at the parameter ``t``.  If ``t`` is ``None`` the resulting polynomial
      is in the polynomial ring `\ZZ['t']`.

    EXAMPLES::

        sage: from sage.combinat.sf.kfpoly import kfpoly
        sage: kfpoly([2,2], [2,1,1])
        t
        sage: kfpoly([4], [2,1,1])
        t^3
        sage: kfpoly([4], [2,2])
        t^2
        sage: kfpoly([1,1,1,1], [2,2])
        0

    TESTS::

        sage: kfpoly([], [])
        1
    """
    if mu == nu:
        return 1

    if t is None:
        t = polygen(ZZ, 't')

    nuc = _Partitions(nu).conjugate()

    f = lambda x: weight(x, t) if x[0] == nuc else 0

    return sum(f(rg) for rg in riggings(mu))


def schur_to_hl(mu, t=None):
    r"""
    Return a dictionary corresponding to `s_\mu` in Hall-Littlewood `P` basis.

    INPUT:

    - ``mu`` -- a partition
    - ``t`` -- an optional parameter (default: the generator from `\ZZ['t']` )

    OUTPUT:

    - a dictionary with the coefficients `K_{\mu\nu}(t)` for `\nu` smaller
      in dominance order than `\mu`

    EXAMPLES::

        sage: from sage.combinat.sf.kfpoly import *
        sage: schur_to_hl([1,1,1])
        {[1, 1, 1]: 1}
        sage: a = schur_to_hl([2,1])
        sage: for mc in sorted(a.items()): print(mc)
        ([1, 1, 1], t^2 + t)
        ([2, 1], 1)
        sage: a = schur_to_hl([3])
        sage: for mc in sorted(a.items()): print(mc)
        ([1, 1, 1], t^3)
        ([2, 1], t)
        ([3], 1)
        sage: a = schur_to_hl([4])
        sage: for mc in sorted(a.items()): print(mc)
        ([1, 1, 1, 1], t^6)
        ([2, 1, 1], t^3)
        ([2, 2], t^2)
        ([3, 1], t)
        ([4], 1)
        sage: a = schur_to_hl([3,1])
        sage: for mc in sorted(a.items()): print(mc)
        ([1, 1, 1, 1], t^5 + t^4 + t^3)
        ([2, 1, 1], t^2 + t)
        ([2, 2], t)
        ([3, 1], 1)
        sage: a = schur_to_hl([2,2])
        sage: for mc in sorted(a.items()): print(mc)
        ([1, 1, 1, 1], t^4 + t^2)
        ([2, 1, 1], t)
        ([2, 2], 1)
        sage: a = schur_to_hl([2,1,1])
        sage: for mc in sorted(a.items()): print(mc)
        ([1, 1, 1, 1], t^3 + t^2 + t)
        ([2, 1, 1], 1)
        sage: a = schur_to_hl([1,1,1,1])
        sage: for mc in sorted(a.items()): print(mc)
        ([1, 1, 1, 1], 1)
        sage: a = schur_to_hl([2,2,2])
        sage: for mc in sorted(a.items()): print(mc)
        ([1, 1, 1, 1, 1, 1], t^9 + t^7 + t^6 + t^5 + t^3)
        ([2, 1, 1, 1, 1], t^4 + t^2)
        ([2, 2, 1, 1], t)
        ([2, 2, 2], 1)
    """
    if mu == []:
        return {mu: 1}
    if t is None:
        t = polygen(ZZ, 't')

    res = {}
    for rg in riggings(mu):
        res[rg[0]] = res.get(rg[0], 0) + weight(rg, t)

    return {key.conjugate(): res[key] for key in res}


def riggings(part):
    r"""
    Generate all possible rigging sequences for a fixed partition ``part``.

    INPUT:

    - ``part`` -- a partition

    OUTPUT: list of riggings associated to the partition ``part``

    EXAMPLES::

        sage: from sage.combinat.sf.kfpoly import *
        sage: riggings([3])
        [[[1, 1, 1]], [[2, 1]], [[3]]]
        sage: riggings([2,1])
        [[[2, 1], [1]], [[3], [1]]]
        sage: riggings([1,1,1])
        [[[3], [2], [1]]]
        sage: riggings([2,2])
        [[[2, 2], [1, 1]], [[3, 1], [1, 1]], [[4], [1, 1]], [[4], [2]]]
        sage: riggings([2,2,2])
        [[[3, 3], [2, 2], [1, 1]],
         [[4, 2], [2, 2], [1, 1]],
         [[5, 1], [2, 2], [1, 1]],
         [[6], [2, 2], [1, 1]],
         [[5, 1], [3, 1], [1, 1]],
         [[6], [3, 1], [1, 1]],
         [[6], [4], [2]]]
    """
    l = len(part)
    res = [ [[],[]] ]
    sa = 0
    for i in sorted(part):
        sa += i
        res = [[new] + nu for nu in res for new in compat(sa, nu[0], nu[1])]

    return [x[:l] for x in res]


def compat(n, mu, nu):
    r"""
    Generate all possible partitions of `n` that can precede `\mu, \nu`
    in a rigging sequence.

    INPUT:

    - ``n`` -- positive integer
    - ``mu``, ``nu`` -- partitions

    OUTPUT: list of partitions

    EXAMPLES::

        sage: from sage.combinat.sf.kfpoly import *
        sage: compat(4, [1], [2,1])
        [[1, 1, 1, 1], [2, 1, 1], [2, 2], [3, 1], [4]]
        sage: compat(3, [1], [2,1])
        [[1, 1, 1], [2, 1], [3]]
        sage: compat(2, [1], [])
        [[2]]
        sage: compat(3, [1], [])
        [[2, 1], [3]]
        sage: compat(3, [2], [1])
        [[3]]
        sage: compat(4, [1,1], [])
        [[2, 2], [3, 1], [4]]
        sage: compat(4, [2], [])
        [[4]]
    """
    l = max(len(mu), len(nu))
    mmu = list(mu) + [0]*(l-len(mu))
    nnu = list(nu) + [0]*(l-len(nu))

    bd = []
    sa = 0
    for i in range(l):
        sa += 2*mmu[i] - nnu[i]
        bd.append(sa)

    for la in ZS1_iterator(n):
        if dom(la, bd):
            return [x.conjugate() for x in _Partitions(la).dominated_partitions()]

    return [] # _Partitions([])


def dom(mup, snu):
    """
    Return ``True`` if ``sum(mu[:i+1]) >= snu[i]`` for all
    ``0 <= i < len(snu)``; otherwise, it returns ``False``.

    INPUT:

    - ``mup`` -- a partition conjugate to ``mu``
    - ``snu`` -- a sequence of positive integers

    OUTPUT: boolean

    EXAMPLES::

        sage: from sage.combinat.sf.kfpoly import *
        sage: dom([3,2,1],[2,4,5])
        True
        sage: dom([3,2,1],[2,4,7])
        False
        sage: dom([3,2,1],[2,6,5])
        False
        sage: dom([3,2,1],[4,4,4])
        False

    TESTS::

        sage: dom([],[])
        True
    """
    if not mup: # mup is empty:
        return not snu # True if and only if snu is empty

    l = len(snu)
    lmup = len(mup)
    # Special case for the largest columns
    if any((k+1)*lmup < snu[k] for k in range(min(mup[-1],l))):
        return False

    pos = mup[-1]
    sa = mup[-1] * lmup
    for i in range(lmup-1, 0, -1):
        for k in range(mup[i-1] - mup[i]):
            if pos >= l: # We've reached the end of snu
                return True
            sa += i
            if sa < snu[pos]:
                return False
            pos += 1
    return all(sa >= snu[j] for j in range(pos, l))


def weight(rg, t=None):
    r"""
    Return the weight of a rigging.

    INPUT:

    - ``rg`` -- a rigging, a list of partitions
    - ``t`` -- an optional parameter, (default: the generator from `\ZZ['t']`)

    OUTPUT: a polynomial in the parameter `t`

    EXAMPLES::

        sage: from sage.combinat.sf.kfpoly import weight
        sage: weight([[2,1], [1]])
        1
        sage: weight([[3], [1]])
        t^2 + t
        sage: weight([[2,1], [3]])
        t^4
        sage: weight([[2, 2], [1, 1]])
        1
        sage: weight([[3, 1], [1, 1]])
        t
        sage: weight([[4], [1, 1]], 2)
        16
        sage: weight([[4], [2]], t=2)
        4
    """
    from sage.combinat.q_analogues import q_binomial
    if t is None:
        t = polygen(ZZ, 't')

    nu = rg + [ [] ]
    l = 1 + max( map(len, nu) )
    nu = [ list(mu) + [0]*l for mu in nu ]
    res = t**int(sum(i * (i-1) // 2 for i in rg[-1]))
    for k in range(1, len(nu)-1):
        sa = 0
        mid = nu[k]
        for i in range( max(len(rg[k]), len(rg[k-1])) ):
            sa += nu[k-1][i] - 2*mid[i] + nu[k+1][i]
            if mid[i] - mid[i+1] + sa >= 0:
                res *= q_binomial(mid[i]-mid[i+1]+sa, sa, t)
            mu = nu[k-1][i] - mid[i]
            res *= t**int(mu * (mu-1) // 2)
    return res
