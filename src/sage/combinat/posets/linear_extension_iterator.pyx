# cython: binding=True
r"""
Fast linear extension iterator
"""
cimport cython

from copy import copy


def _linear_extension_prepare(D):
    r"""
    The preprocessing routine in Figure 7 of "Generating Linear
    Extensions Fast" by Preusse and Ruskey.

    INPUT:

    - ``D`` -- the Hasse diagram of a poset

    OUTPUT:

    - a triple ``(le, a, b)``, where ``le`` is the first linear
      extension, and ``a`` and ``b`` are lists such that ``a[i]`` and
      ``b[i]`` are minimal elements of ``D`` after removing ``a[:i]``
      and ``b[:i]``.

    TESTS::

        sage: from sage.combinat.posets.linear_extension_iterator import _linear_extension_prepare
        sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })._hasse_diagram
        sage: _linear_extension_prepare(D)
        ([0, 1, 2, 3, 4], [1, 3], [2, 4])
    """
    dag_copy = copy(D)  # this copy is destroyed during preparation
    le = []
    a = []
    b = []

    # the preprocessing routine found in Figure 7 of
    # "Generating Linear Extensions Fast" by
    # Pruesse and Ruskey
    while dag_copy.num_verts() != 0:
        # find all the minimal elements of dag_copy
        minimal_elements = dag_copy.sources()
        if not minimal_elements:
            raise ValueError("the digraph must be acyclic to have linear extensions")
        elif len(minimal_elements) == 1:
            le.append(minimal_elements[0])
            dag_copy.delete_vertex(minimal_elements[0])
        else:
            ap = minimal_elements[0]
            bp = minimal_elements[1]
            a.append(ap)
            b.append(bp)
            le.append(ap)
            le.append(bp)
            dag_copy.delete_vertex(ap)
            dag_copy.delete_vertex(bp)

    return (le, a, b)


@cython.wraparound(False)
@cython.boundscheck(False)
cdef void _linear_extension_switch(list _le, list _a, list _b, list _is_plus, Py_ssize_t i) noexcept:
    """
    This implements the ``Switch`` procedure described on page 7
    of "Generating Linear Extensions Fast" by Pruesse and Ruskey.

    If ``i == -1``, then the sign is changed.  Otherwise, then
    ``_a[i]`` and ``_b[i]`` are transposed.
    """
    cdef Py_ssize_t a_index, b_index
    if i == -1:
        _is_plus[0] = not _is_plus[0]
    else:
        a = _a[i]
        b = _b[i]
        a_index = _le.index(a)
        b_index = _le.index(b)
        _le[a_index] = b
        _le[b_index] = a
        _b[i] = a
        _a[i] = b


@cython.wraparound(False)
@cython.boundscheck(False)
cdef bint _linear_extension_right_a(_D, list _le, list _a, list _b, Py_ssize_t i) noexcept:
    """
    Return ``True`` if and only if ``_a[i]`` is incomparable with the
    element to its right in ``_le`` and the element to the right is
    not ``_b[i]``.

    This is the ``Right`` function described on page 8 of
    "Generating Linear Extensions Fast" by Pruesse and Ruskey.

    ::

        sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })._hasse_diagram             # not tested
        sage: _linear_extension_right_a(D, [0, 1, 2, 4, 3], [1, 4], [2, 3], 0)  # not tested
        False
        sage: _linear_extension_right_a(D, [0, 1, 2, 4, 3], [1, 4], [2, 3], 1)  # not tested
        False
    """
    cdef Py_ssize_t yindex
    x = _a[i]
    yindex = _le.index(x) + 1
    if yindex >= len(_le):
        return False
    y = _le[yindex]
    return y != _b[i] and _D.are_incomparable(x, y)


@cython.wraparound(False)
@cython.boundscheck(False)
cdef bint _linear_extension_right_b(_D, list _le, list _a, list _b, Py_ssize_t i) noexcept:
    """
    Return ``True`` if and only if ``_b[i]`` is incomparable with the
    elements to its right in ``_le``.

    This is the ``Right`` function described on page 8 of
    "Generating Linear Extensions Fast" by Pruesse and Ruskey.

    ::

        sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })._hasse_diagram             # not tested
        sage: _linear_extension_right_b(D, [0, 1, 2, 4, 3], [1, 4], [2, 3], 0)  # not tested
        False
        sage: _linear_extension_right_b(D, [0, 1, 2, 4, 3], [1, 4], [2, 3], 1)  # not tested
        False
    """
    cdef Py_ssize_t yindex
    x = _b[i]
    yindex = _le.index(x) + 1
    if yindex >= len(_le):
        return False
    y = _le[yindex]
    return _D.are_incomparable(x, y)


@cython.wraparound(False)
@cython.boundscheck(False)
def _linear_extension_gen(_D, list _le, list _a, list _b, list _is_plus, Py_ssize_t i):
    """
    This a Python version of the GenLE routine found in Figure 8
    of "Generating Linear Extensions Fast" by Pruesse and Ruskey.

    TESTS::

        sage: from sage.combinat.posets.linear_extension_iterator import _linear_extension_prepare, _linear_extension_gen
        sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })._hasse_diagram
        sage: le, a, b = _linear_extension_prepare(D)
        sage: [e for e in _linear_extension_gen(D, le, a, b, [True], len(a)-1)]         # needs sage.modules
        [[0, 2, 1, 3, 4]]
    """
    cdef int mra, mrb, mla
    cdef Py_ssize_t index, index1
    cdef bint typical
    if i == -1:
        return

    for e in _linear_extension_gen(_D, _le, _a, _b, _is_plus, i-1):
        yield e
    mrb = 0
    typical = False
    while _linear_extension_right_b(_D, _le, _a, _b, i):
        mrb += 1
        # move_right
        index = _le.index(_b[i])
        index1 = index + 1
        _le[index] = _le[index1]
        _le[index1] = _b[i]
        if _is_plus[0]:
            yield _le[:]

        for e in _linear_extension_gen(_D, _le, _a, _b, _is_plus, i-1):
            yield e
        mra = 0
        while _linear_extension_right_a(_D, _le, _a, _b, i):
            typical = True
            mra += 1
            # move_right
            index = _le.index(_a[i])
            index1 = index+1
            _le[index] = _le[index1]
            _le[index1] = _a[i]
            if _is_plus[0]:
                yield _le[:]

            for e in _linear_extension_gen(_D, _le, _a, _b, _is_plus, i-1):
                yield e

        if typical:
            _linear_extension_switch(_le, _a, _b, _is_plus, i-1)
            if _is_plus[0]:
                yield _le[:]

            for e in _linear_extension_gen(_D, _le, _a, _b, _is_plus, i-1):
                yield e
            if mrb % 2 == 1:
                mla = mra - 1
            else:
                mla = mra + 1
            for _ in range(mla):
                # move_left
                index = _le.index(_a[i])
                index1 = index-1
                _le[index] = _le[index1]
                _le[index1] = _a[i]
                if _is_plus[0]:
                    yield _le[:]

                for e in _linear_extension_gen(_D, _le, _a, _b, _is_plus, i-1):
                    yield e

    if typical and (mrb % 2 == 1):
        # move_left
        index = _le.index(_a[i])
        index1 = index-1
        _le[index] = _le[index1]
        _le[index1] = _a[i]
        if _is_plus[0]:
            yield _le[:]
    else:
        _linear_extension_switch(_le, _a, _b, _is_plus, i-1)
        if _is_plus[0]:
            yield _le[:]
    for e in _linear_extension_gen(_D, _le, _a, _b, _is_plus, i-1):
        yield e
    for _ in range(mrb):
        # move_left
        index = _le.index(_b[i])
        index1 = index-1
        _le[index] = _le[index1]
        _le[index1] = _b[i]
        if _is_plus[0]:
            yield _le[:]

        for e in _linear_extension_gen(_D, _le, _a, _b, _is_plus, i-1):
            yield e


def linear_extension_iterator(D):
    """
    Iterate over the linear extensions of the poset.

    The list ``_le`` keeps track of the current linear extensions.  The
    boolean variable ``is_plus`` keeps track of the "sign".

    INPUT:

    - ``D`` -- the Hasse diagram of a poset

    .. WARNING::

        It is assumed that ``D`` is not modified while the linear
        extensions are generated.

    EXAMPLES::

        sage: from sage.combinat.posets.linear_extension_iterator import linear_extension_iterator
        sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })._hasse_diagram
        sage: list(linear_extension_iterator(D))                                        # needs sage.modules
        [[0, 1, 2, 3, 4],
         [0, 2, 1, 3, 4],
         [0, 2, 1, 4, 3],
         [0, 2, 4, 1, 3],
         [0, 1, 2, 4, 3]]

        sage: D = posets.BooleanLattice(3)._hasse_diagram
        sage: len(list(linear_extension_iterator(D)))                                   # needs sage.modules
        48

        sage: D = posets.AntichainPoset(9)._hasse_diagram
        sage: len(list(linear_extension_iterator(D))) == factorial(9)   # long time, needs sage.modules
        True
    """
    _le, _a, _b = _linear_extension_prepare(D)
    _max_pair = len(_a) - 1
    _is_plus = [True]  # this is modified by _linear_extension_switch

    yield _le[:]
    for e in _linear_extension_gen(D, _le, _a, _b, _is_plus, _max_pair):
        yield e
    _linear_extension_switch(_le, _a, _b, _is_plus, _max_pair)
    if _is_plus[0]:
        yield _le[:]
    for e in _linear_extension_gen(D, _le, _a, _b, _is_plus, _max_pair):
        yield e
