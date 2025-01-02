r"""
Gray codes

Functions
---------
"""


def product(m):
    r"""
    Iterator over the switch for the iteration of the product
    `[m_0] \times [m_1] \ldots \times [m_k]`.

    The iterator return at each step a pair ``(p,i)`` which corresponds to the
    modification to perform to get the next element. More precisely, one has to
    apply the increment ``i`` at the position ``p``. By construction, the
    increment is either ``+1`` or ``-1``.

    This is algorithm H in [Knu2011]_ Section 7.2.1.1, "Generating All
    `n`-Tuples": loopless reflected mixed-radix Gray generation.

    INPUT:

    - ``m`` -- list or tuple of positive integers that correspond to the size
      of the sets in the product

    EXAMPLES::

        sage: from sage.combinat.gray_codes import product
        sage: l = [0,0,0]
        sage: for p,i in product([3,3,3]):
        ....:     l[p] += i
        ....:     print(l)
        [1, 0, 0]
        [2, 0, 0]
        [2, 1, 0]
        [1, 1, 0]
        [0, 1, 0]
        [0, 2, 0]
        [1, 2, 0]
        [2, 2, 0]
        [2, 2, 1]
        [1, 2, 1]
        [0, 2, 1]
        [0, 1, 1]
        [1, 1, 1]
        [2, 1, 1]
        [2, 0, 1]
        [1, 0, 1]
        [0, 0, 1]
        [0, 0, 2]
        [1, 0, 2]
        [2, 0, 2]
        [2, 1, 2]
        [1, 1, 2]
        [0, 1, 2]
        [0, 2, 2]
        [1, 2, 2]
        [2, 2, 2]
        sage: l = [0,0]
        sage: for i,j in product([2,1]):
        ....:     l[i] += j
        ....:     print(l)
        [1, 0]

    TESTS::

        sage: for t in [[2,2,2],[2,1,2],[3,2,1],[2,1,3]]:
        ....:     assert sum(1 for _ in product(t)) == prod(t)-1
    """
    # n is the length of the element (we ignore sets of size 1)
    n = k = 0

    new_m = []   # will be the set of upper bounds m_i different from 1
    mm = []      # index of each set (we skip sets of cardinality 1)
    for i in m:
        i = int(i)
        if i <= 0:
            raise ValueError("accept only positive integers")
        if i > 1:
            new_m.append(i-1)
            mm.append(k)
            n += 1
        k += 1

    m = new_m
    f = list(range(n + 1))  # focus pointer
    o = [1] * n     # switch +1 or -1
    a = [0] * n     # current element of the product

    j = f[0]
    while j != n:
        f[0] = 0
        oo = o[j]
        a[j] += oo
        if a[j] == 0 or a[j] == m[j]:
            f[j] = f[j+1]
            f[j+1] = j+1
            o[j] = -oo

        yield (mm[j], oo)

        j = f[0]


def combinations(n, t):
    r"""
    Iterator through the switches of the revolving door algorithm.

    The revolving door algorithm is a way to generate all combinations of a set
    (i.e. the subset of given cardinality) in such way that two consecutive
    subsets differ by one element. At each step, the iterator output a pair
    ``(i,j)`` where the item ``i`` has to be removed and ``j`` has to be added.

    The ground set is always `\{0, 1, ..., n-1\}`. Note that ``n`` can be
    infinity in that algorithm.

    See [Knu2011]_ Section 7.2.1.3, "Generating All Combinations".

    INPUT:

    - ``n`` -- integer or ``Infinity``; size of the ground set

    - ``t`` -- integer; size of the subsets

    EXAMPLES::

        sage: from sage.combinat.gray_codes import combinations
        sage: b = [1, 1, 1, 0, 0]
        sage: for i,j in combinations(5,3):
        ....:     b[i] = 0; b[j] = 1
        ....:     print(b)
        [1, 0, 1, 1, 0]
        [0, 1, 1, 1, 0]
        [1, 1, 0, 1, 0]
        [1, 0, 0, 1, 1]
        [0, 1, 0, 1, 1]
        [0, 0, 1, 1, 1]
        [1, 0, 1, 0, 1]
        [0, 1, 1, 0, 1]
        [1, 1, 0, 0, 1]

        sage: s = set([0,1])
        sage: for i,j in combinations(4,2):
        ....:     s.remove(i)
        ....:     s.add(j)
        ....:     print(sorted(s))
        [1, 2]
        [0, 2]
        [2, 3]
        [1, 3]
        [0, 3]

    Note that ``n`` can be infinity::

        sage: c = combinations(Infinity,4)
        sage: s = set([0,1,2,3])
        sage: for _ in range(10):
        ....:     i,j = next(c)
        ....:     s.remove(i); s.add(j)
        ....:     print(sorted(s))
        [0, 1, 3, 4]
        [1, 2, 3, 4]
        [0, 2, 3, 4]
        [0, 1, 2, 4]
        [0, 1, 4, 5]
        [1, 2, 4, 5]
        [0, 2, 4, 5]
        [2, 3, 4, 5]
        [1, 3, 4, 5]
        [0, 3, 4, 5]
        sage: for _ in range(1000):
        ....:     i,j = next(c)
        ....:     s.remove(i); s.add(j)
        sage: sorted(s)
        [0, 4, 13, 14]

    TESTS::

        sage: def check_sets_from_iter(n, k):
        ....:     l = []
        ....:     s = set(range(k))
        ....:     l.append(frozenset(s))
        ....:     for i,j in combinations(n,k):
        ....:         s.remove(i)
        ....:         s.add(j)
        ....:         assert len(s) == k
        ....:         l.append(frozenset(s))
        ....:     assert len(set(l)) == binomial(n,k)
        sage: check_sets_from_iter(9,5)
        sage: check_sets_from_iter(8,5)
        sage: check_sets_from_iter(5,6)
        Traceback (most recent call last):
        ...
        AssertionError: t(=6) must be >=0 and <=n(=5)
    """
    from sage.rings.infinity import Infinity
    t = int(t)
    if n != Infinity:
        n = int(n)
    else:
        n = Infinity
    assert 0 <= t and t <= n, "t(={}) must be >=0 and <=n(={})".format(t,n)
    if t == 0 or t == n:
        return iter([])
    if t % 2:
        return _revolving_door_odd(n,t)
    else:
        return _revolving_door_even(n,t)


def _revolving_door_odd(n, t):
    r"""
    Revolving door switch for odd `t`.

    TESTS::

        sage: from sage.combinat.gray_codes import _revolving_door_odd
        sage: sum(1 for _ in _revolving_door_odd(13,3)) == binomial(13,3) - 1
        True
        sage: sum(1 for _ in _revolving_door_odd(10,5)) == binomial(10,5) - 1
        True
    """
    # note: the numbering of the steps below follows Knuth TAOCP
    c = list(range(t)) + [n]    # the combination (ordered list of numbers of length t+1)

    while True:
        # R3 : easy case
        if c[0] + 1 < c[1]:
            yield c[0], c[0]+1
            c[0] += 1
            continue

        j = 1
        while j < t:
            # R4 : try to decrease c[j]
            # at this point c[j] = c[j-1] + 1
            if c[j] > j:
                yield c[j], j-1
                c[j] = c[j-1]
                c[j-1] = j-1
                break
            j += 1

            # R5 : try to increase c[j]
            # at this point c[j-1] = j-1
            if c[j] + 1 < c[j+1]:
                yield c[j-1], c[j]+1
                c[j-1] = c[j]
                c[j] += 1
                break
            j += 1

        else: # j == t
            break


def _revolving_door_even(n, t):
    r"""
    Revolving door algorithm for even `t`.

    TESTS::

        sage: from sage.combinat.gray_codes import _revolving_door_even
        sage: sum(1 for _ in _revolving_door_even(13,4)) == binomial(13,4) - 1
        True
        sage: sum(1 for _ in _revolving_door_even(12,6)) == binomial(12,6) - 1
        True
    """
    # note: the numbering of the steps below follows Knuth TAOCP

    c = list(range(t)) + [n]    # the combination (ordered list of numbers of length t+1)

    while True:
        # R3 : easy case
        if c[0] > 0:
            yield c[0], c[0]-1
            c[0] -= 1
            continue

        j = 1
        # R5 : try to increase c[j]
        # at this point c[j-1] = j-1
        if c[j] + 1 < c[j+1]:
            yield c[j-1], c[j]+1
            c[j-1] = c[j]
            c[j] += 1
            continue
        j += 1

        while j < t:
            # R4 : try to decrease c[j]
            # at this point c[j] = c[j-1] + 1
            if c[j] > j:
                yield c[j], j-1
                c[j] = c[j-1]
                c[j-1] = j-1
                break
            j += 1

            # R5 : try to increase c[j]
            # at this point c[j-1] = j-1
            if c[j] + 1 < c[j+1]:
                yield c[j-1], c[j] + 1
                c[j-1] = c[j]
                c[j] += 1
                break
            j += 1

        else: # j == t
            break
