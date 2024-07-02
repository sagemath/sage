# sage.doctest: needs sage.modules

cpdef int _flip_c(W, set positions, list extended_root_conf_indices,
                  int i, side='both') noexcept:
    r"""
    Flip a facet.

    INPUT:

    - ``W`` -- a Coxeter group
    - ``positions`` -- the positions of the elements of the facet
    - ``extended_root_conf_indices`` -- also attached to the facet ?
    - ``i`` -- the position where to flip
    - ``side`` -- string; ``'positive'``, ``'negative'``, or ``'both'`` (default)

    OUTPUT: the new position `j` that has replaced `i`

    EXAMPLES::

        sage: from sage.combinat.subword_complex_c import _flip_c

        sage: # optional - gap3
        sage: W = ReflectionGroup(['A',2])
        sage: w = W.from_reduced_word([1,2,1])
        sage: SC = SubwordComplex([1,2,1,2,1], w)
        sage: F = SC([0, 1])
        sage: _flip_c(W, set([0,1]), F._extended_root_configuration_indices(), 1)
        4
        sage: _flip_c(W, set([0,1]), F._extended_root_configuration_indices(), 0)
        3

        sage: W = CoxeterGroup(['A',2])
        sage: w = W.from_reduced_word([1,2,1])
        sage: SC = SubwordComplex([1,2,1,2,1], w)
        sage: F = SC([0, 1])
        sage: _flip_c(W, set([0,1]), F._extended_root_configuration_indices(), 1)
        4
        sage: _flip_c(W, set([0,1]), F._extended_root_configuration_indices(), 0)
        3
    """
    cdef int r, nr_ref, r_minus, j, k
    cdef list R
    r = extended_root_conf_indices[i]
    nr_ref = W.number_of_reflections()
    r_minus = (r + nr_ref) % (2 * nr_ref)  # get the negative root -r
    j = i
    for k in range(len(extended_root_conf_indices)):
        if extended_root_conf_indices[k] == r_minus and k not in positions and side != "positive":
            j = k
            break
        elif extended_root_conf_indices[k] == r and k not in positions and side != "negative":
            j = k
            break
    positions.remove(i)
    positions.add(j)
    R = list(W.reflections())
    if j != i:
        t = R[min(r, r_minus)]
        for k in range(min(i, j) + 1, max(i, j) + 1):
            extended_root_conf_indices[k] = t.action_on_root_indices(extended_root_conf_indices[k], side='left')
    return j


cpdef list _construct_facets_c(tuple Q, w, int n=-1, int pos=0, int l=-1):
    r"""
    Return the list of facets of the subword complex associated to the
    word `Q` and the element `w` in a Coxeter group `W`.

    EXAMPLES::

        sage: from sage.combinat.subword_complex_c import _construct_facets_c
        sage: W = CoxeterGroup(['A',2])
        sage: w = W.from_reduced_word([1,2])
        sage: _construct_facets_c((2,1), w)
        []
        sage: _construct_facets_c((2,1,2), w)
        [[0]]
        sage: _construct_facets_c((2,1,2,1), w)
        [[0, 3]]
        sage: w = W.from_reduced_word([1,2,1])
        sage: _construct_facets_c((1,2,1,2,1), w)
        [[0, 1], [0, 4], [1, 2], [2, 3], [3, 4]]
    """
    cdef int s
    cdef list X, Y
    if n == -1:
        n = len(Q)
    if l == -1:
        first = True
        l = w.length()
    else:
        first = False

    if l == 0:
        return [list(range(pos, n))]
    elif n < l + pos:
        return []

    s = Q[pos]
    X = _construct_facets_c(Q, w, n=n, pos=pos + 1, l=l)
    for f in X:
        f.append(pos)

    if w.has_left_descent(s):
        Y = _construct_facets_c(Q, w.apply_simple_reflection_left(s),
                                n=n, pos=pos + 1, l=l - 1)
        Y = X + Y
    else:
        Y = X
    if first:
        return sorted([sorted(x) for x in Y])
    else:
        return Y
