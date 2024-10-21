# sage_setup: distribution = sagemath-categories
r"""
Subsets satisfying a hereditary property
"""
# ****************************************************************************
#  Copyright (C) 2014 Nathann Cohen <nathann.cohen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************


def subsets_with_hereditary_property(f, X, max_obstruction_size=None, ncpus=1):
    r"""
    Return all subsets `S` of `X` such that `f(S)` is true.

    The boolean function `f` must be decreasing, i.e. `f(S)\Rightarrow f(S')` if
    `S'\subseteq S`.

    This function is implemented to call `f` as few times as possible. More
    precisely, `f` will be called on all sets `S` such that `f(S)` is true, as
    well as on all inclusionwise minimal sets `S` such that `f(S)` is false.

    The problem that this function answers is also known as the learning problem
    on monotone boolean functions, or as computing the set of winning coalitions
    in a simple game.

    INPUT:

    - ``f`` -- boolean function which takes as input a list of elements from
      ``X``

    - ``X`` -- list/iterable

    - ``max_obstruction_size`` -- integer; if you know that there is
      a `k` such that `f(S)` is true if and only if `f(S')` is true
      for all `S'\subseteq S` with `S'\leq k`, set
      ``max_obstruction_size=k``. It may dramatically decrease the
      number of calls to `f`. Set to ``None`` by default, meaning
      `k=|X|`.

    - ``ncpus`` -- number of cpus to use for this computation. Note that
      changing the value from `1` (default) to anything different *enables*
      parallel computations which can have a cost by itself, so it is not
      necessarily a good move. In some cases, however, it is a *great* move. Set
      to ``None`` to automatically detect and use the maximum number of cpus
      available.

      .. NOTE::

          Parallel computations are performed through the
          :func:`~sage.parallel.decorate.parallel` decorator. See its
          documentation for more information, in particular with respect to the
          memory context.

    EXAMPLES:

    Sets whose elements all have the same remainder mod 2::

        sage: from sage.combinat.subsets_hereditary import subsets_with_hereditary_property
        sage: def f(x):
        ....:     return (not x) or all(xx % 2 == x[0] % 2 for xx in x)
        sage: list(subsets_with_hereditary_property(f, range(4)))
        [[], [0], [1], [2], [3], [0, 2], [1, 3]]

    Same, on two threads::

        sage: sorted(subsets_with_hereditary_property(f, range(4), ncpus=2))
        [[], [0], [0, 2], [1], [1, 3], [2], [3]]

    One can use this function to compute the independent sets of a graph. We
    know, however, that in this case the maximum obstructions are the edges, and
    have size 2. We can thus set ``max_obstruction_size=2``, which reduces the
    number of calls to `f` from 91 to 56::

        sage: # needs sage.graphs
        sage: num_calls = 0
        sage: g = graphs.PetersenGraph()
        sage: def is_independent_set(S):
        ....:     global num_calls
        ....:     num_calls += 1
        ....:     return g.subgraph(S).size() == 0
        sage: l1 = list(subsets_with_hereditary_property(is_independent_set,
        ....:                                            g.vertices(sort=False)))
        sage: num_calls
        91
        sage: num_calls = 0
        sage: l2 = list(subsets_with_hereditary_property(is_independent_set,
        ....:                                            g.vertices(sort=False),
        ....:                                            max_obstruction_size=2))
        sage: num_calls
        56
        sage: l1 == l2
        True

    TESTS::

        sage: list(subsets_with_hereditary_property(lambda x: False, range(4)))
        []
        sage: list(subsets_with_hereditary_property(lambda x: len(x)<1, range(4)))
        [[]]
        sage: list(subsets_with_hereditary_property(lambda x: True, range(2)))
        [[], [0], [1], [0, 1]]
    """
    from sage.data_structures.bitset import Bitset
    # About the implementation:
    #
    # 1) We work on X={0,...,n-1} but remember X to return correctly
    #    labelled answers.
    #
    # 2) We maintain a list of sets S such that f(S)=0 (i.e. 'no-sets'), in
    #    order to filter out larger sets for which f is necessarily False.
    #
    # 3) Those sets are stored in an array: bs[i] represents the set of all
    #    no-sets S we found such that i is NOT in S. Why ? Because it makes it
    #    easy to filter out sets: if a set S' whose *complement* is
    #    {i1,i2,...,ik} is such that bs[i1]&bs[i2]&...&bs[ik] is nonempty then
    #    f(S') is necessarily False.
    X_labels = list(X)
    n = len(X_labels)
    X = set(range(n))
    if max_obstruction_size is None:
        max_obstruction_size = n

    bs = [Bitset([], 1) for _ in range(n)]  # collection of no-set
    nforb = 1                              # number of no-sets stored
    current_layer = [[]]                  # all yes-sets of size 'current_size'
    current_size = 0

    def explore_neighbors(s):
        r"""
        Explore the successors of a set s.

        The successors of a set s are all the sets s+[i] where max(s)<i. This
        function returns them all as a partition `(yes_sets,no_sets)`.
        """
        new_yes_sets = []
        new_no_sets = []
        for i in range((s[-1] + 1 if s else 0), n):        # all ways to extend it
            s_plus_i = s + [i]                             # the extended set
            s_plus_i_c = Bitset(s_plus_i, n).complement()  # .. and its complement

            # Filter a no-set using the data collected so far.
            inter = Bitset([], nforb).complement()
            for j in s_plus_i_c:
                inter.intersection_update(bs[j])

            # If we cannot decide yet we must call f(S)
            if not inter:
                if set_size >= max_obstruction_size or f([X_labels[xx] for xx in s_plus_i]):
                    new_yes_sets.append(s_plus_i)
                else:
                    new_no_sets.append(s_plus_i)
        return (new_yes_sets, new_no_sets)

    # The empty set
    if f([]):
        yield []
    else:
        return

    if ncpus != 1:
        from sage.parallel.decorate import parallel
        explore_neighbors_paral = parallel(ncpus=ncpus)(explore_neighbors)

    # All sets of size 0, then size 1, then ...
    set_size = -1
    while current_layer:
        set_size += 1
        new_no_sets = []
        new_yes_sets = []

        if ncpus == 1:
            yes_no_iter = (explore_neighbors(s) for s in current_layer)
        else:
            yes_no_iter = ((yes, no) for (_, (yes, no)) in explore_neighbors_paral(current_layer))

        for yes, no in yes_no_iter:
            new_yes_sets.extend(yes)
            new_no_sets.extend(no)
            for s in yes:
                yield [X_labels[xx] for xx in s]

        current_layer = new_yes_sets

        # Update bs with the new no-sets
        new_nforb = nforb + len(new_no_sets)
        for b in bs:  # resize the bitsets
            b.add(new_nforb)
            b.discard(new_nforb)
        for i, s in enumerate(new_no_sets):  # Fill the new entries
            for j in X.difference(s):
                bs[j].add(i + nforb)
        nforb = new_nforb
        current_size += 1

    # Did we forget to return X itself ?
    #
    # If we did, this was probably the worst choice of algorithm for we computed
    # f(X) for all 2^n sets X, but well...
    if (current_size == len(X) and nforb == 1 and f(X_labels)):
        yield X_labels
