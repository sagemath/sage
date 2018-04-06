# the following functions have been moved to sage.combinat.combination
import sage.combinat.combination as combination

def rank(comb, n, check=True):
    """
    Return the rank of ``comb`` in the subsets of ``range(n)`` of size ``k``
    where ``k`` is the length of ``comb``.

    The algorithm used is based on combinadics and James McCaffrey's
    MSDN article. See :wikipedia:`Combinadic`.

    EXAMPLES::

        sage: import sage.combinat.choose_nk as choose_nk
        sage: choose_nk.rank((), 3)
        doctest:...: DeprecationWarning: choose_nk.rank is deprecated and will be removed. Use combination.rank instead
        See http://trac.sagemath.org/18674 for details.
        0
        sage: choose_nk.rank((0,), 3)
        0
        sage: choose_nk.rank((1,), 3)
        1
        sage: choose_nk.rank((2,), 3)
        2
        sage: choose_nk.rank((0,1), 3)
        0
        sage: choose_nk.rank((0,2), 3)
        1
        sage: choose_nk.rank((1,2), 3)
        2
        sage: choose_nk.rank((0,1,2), 3)
        0

        sage: choose_nk.rank((0,1,2,3), 3)
        Traceback (most recent call last):
        ...
        ValueError: len(comb) must be <= n
        sage: choose_nk.rank((0,0), 2)
        Traceback (most recent call last):
        ...
        ValueError: comb must be a subword of (0,1,...,n)

        sage: choose_nk.rank([1,2], 3)
        2
        sage: choose_nk.rank([0,1,2], 3)
        0
    """
    from sage.misc.superseded import deprecation
    deprecation(18674, "choose_nk.rank is deprecated and will be removed. Use combination.rank instead")
    return combination.rank(comb, n, check)

def from_rank(r, n, k):
    r"""
    Returns the combination of rank ``r`` in the subsets of ``range(n)`` of
    size ``k`` when listed in lexicographic order.

    The algorithm used is based on combinadics and James McCaffrey's
    MSDN article. See :wikipedia:`Combinadic`.

    EXAMPLES::

        sage: import sage.combinat.choose_nk as choose_nk
        sage: choose_nk.from_rank(0,3,0)
        doctest:...: DeprecationWarning: choose_nk.from_rank is deprecated and will be removed. Use combination.from_rank instead
        See http://trac.sagemath.org/18674 for details.
        ()
        sage: choose_nk.from_rank(0,3,1)
        (0,)
        sage: choose_nk.from_rank(1,3,1)
        (1,)
        sage: choose_nk.from_rank(2,3,1)
        (2,)
        sage: choose_nk.from_rank(0,3,2)
        (0, 1)
        sage: choose_nk.from_rank(1,3,2)
        (0, 2)
        sage: choose_nk.from_rank(2,3,2)
        (1, 2)
        sage: choose_nk.from_rank(0,3,3)
        (0, 1, 2)
    """
    from sage.misc.superseded import deprecation
    deprecation(18674, "choose_nk.from_rank is deprecated and will be removed. Use combination.from_rank instead")
    return combination.from_rank(r, n, k)

