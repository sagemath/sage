r"""
Sidon sets and their generalizations, Sidon `g`-sets

AUTHORS:

- Martin Raum (07-25-2011)
"""
# ****************************************************************************
#                 Copyright (C) 2011 Martin Raum
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.sets.set import Set
from sage.misc.cachefunc import cached_function
from sage.rings.integer import Integer


def sidon_sets(N, g = 1):
    r"""
    Return the set of all Sidon-`g` sets that have elements less than or equal
    to `N`.

    A Sidon-`g` set is a set of positive integers `A \subset [1, N]` such
    that any integer `M` can be obtain at most `g` times as sums of unordered pairs of
    elements of `A` (the two elements are not necessary distinct):

    .. MATH::

        \#\{ (a_i, a_j) | a_i, a_j \in A, a_i + a_j = M,a_i \leq a_j \} \leq g

    INPUT:

    - `N` -- A positive integer.
    - `g` -- A positive integer (default: `1`).

    OUTPUT:

    - A Sage set with categories whose element are also set of integers.

    EXAMPLES::

        sage: S = sidon_sets(3, 2)
        sage: sorted(S, key=str)
        [{1, 2, 3}, {1, 2}, {1, 3}, {1}, {2, 3}, {2}, {3}, {}]
        sage: S.cardinality()
        8
        sage: S.category()
        Category of finite sets
        sage: sid = S.an_element()
        sage: sid
        {2}
        sage: sid.category()
        Category of finite sets

    TESTS::

        sage: S = sidon_sets(10)
        sage: TestSuite(S).run()
        sage: Set([1,2,4,8,13]) in sidon_sets(13)
        True

    The following piece of code computes the first values of the Sloane sequence
    entitled 'Length of shortest (or optimal) Golomb ruler with n marks' with a
    very dumb algorithm. (sequence identifier A003022)::

        sage: n = 1
        sage: L = []
        sage: for i in range(1,19):
        ....:     nb = max([S.cardinality() for S in sidon_sets(i)])
        ....:     if nb > n:
        ....:         L.append(i-1)
        ....:         n = nb
        sage: L
        [1, 3, 6, 11, 17]

    The following tests check that some generalized Sidon sets satisfy the right
    conditions, using a dumb but exhaustive algorithm::

        sage: from itertools import groupby
        sage: all(all(l <= 3 for l in map(lambda s: len(list(s[1])), groupby(sorted(a + ap for a in sid for ap in sid if a >= ap), lambda s: s))) for sid in sidon_sets(10, 3))
        True
        sage: all(all(l <= 5 for l in map(lambda s: len(list(s[1])), groupby(sorted(a + ap for a in sid for ap in sid if a >= ap), lambda s: s))) for sid in sidon_sets(10, 5))
        True

    Checking of arguments::

        sage: sidon_sets(1,1)
        {{}, {1}}
        sage: sidon_sets(-1,3)
        Traceback (most recent call last):
        ...
        ValueError: N must be a positive integer
        sage: sidon_sets(1, -3)
        Traceback (most recent call last):
        ...
        ValueError: g must be a positive integer
    """
    if not isinstance(N, (int, Integer)) or N < 1:
        raise ValueError( "N must be a positive integer" )
    elif not isinstance(g, (int, Integer)) or g < 1:
        raise ValueError( "g must be a positive integer" )
    return sidon_sets_rec(N, g=g)


# This recursive and cached slave function is mainly here because
# caching the user entry function 'sidon_sets' prevents it from
# appearing in the built documentation.
@cached_function
def sidon_sets_rec(N, g=1):
    r"""
    Return the set of all Sidon-`g` sets that have elements less than or equal
    to `N` without checking the arguments. This internal function should not
    be call directly by user.

    TESTS::

        sage: from sage.combinat.sidon_sets import sidon_sets_rec
        sage: sorted(sidon_sets_rec(3,2), key=str)
        [{1, 2, 3}, {1, 2}, {1, 3}, {1}, {2, 3}, {2}, {3}, {}]
    """
    if N == 1:
        return Set([Set([]), Set([1])])

    pre_sidons = sidon_sets(N - 1, g)
    sidons = set(pre_sidons)
    for psid in pre_sidons:
        psid_shift = Set([n - 1 for n in psid if n != 1] + [N - 1])
        if psid_shift not in pre_sidons:
            continue

        if 1 not in psid:
            add_sid = True
        else:
            add_sid = True
            Np1_count = 0
            for n in psid:
                if N + 1 - n in psid and 2 * n <= N + 1:
                    Np1_count += 1
                    if Np1_count >= g:
                        add_sid = False
                        break

        if add_sid:
            sidons.add(Set(psid.list() + [N]))

    return Set(sidons)
