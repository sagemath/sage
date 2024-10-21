# sage_setup: distribution = sagemath-gap
# sage.doctest: needs sage.combinat
r"""
Tools for enumeration modulo the action of a permutation group
"""
#*****************************************************************************
#    Copyright (C) 2010-12 Nicolas Borie <nicolas.borie at math dot u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#              The full text of the GPL is available at:
#                    http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.perm_gps.permgroup_element cimport PermutationGroupElement

cpdef list all_children(ClonableIntArray v, int max_part):
    r"""
    Return all the children of an integer vector (:class:`~sage.structure.list_clone.ClonableIntArray`)
    ``v`` in the tree of enumeration by lexicographic order. The children of
    an integer vector ``v`` whose entries have the sum `n` are all integer
    vectors of sum `n+1` which follow ``v`` in the lexicographic order.

    That means this function adds `1` on the last nonzero entries and the
    following ones. For an integer vector `v` such that

    .. MATH::

        v = [ \dots, a , 0, 0 ]  \text{ with } a \neq 0,

    then, the list of children is

    .. MATH::

        [ [ \dots, a+1 , 0, 0 ] , [ \dots, a , 1, 0 ], [ \dots, a , 0, 1 ] ].

    EXAMPLES::

        sage: from sage.combinat.enumeration_mod_permgroup import all_children
        sage: from sage.structure.list_clone_demo import IncreasingIntArrays
        sage: v = IncreasingIntArrays()([1,2,3,4])
        sage: all_children(v, -1)
        [[1, 2, 3, 5]]
    """
    cdef int l, j
    cdef list all_children
    cdef ClonableIntArray child
    l = len(v)
    all_children = []
    j = l - 1
    while v._list[j] == 0 and j >= 1:
        j = j - 1
    for i in range(j,l):
        if (max_part < 0) or ((v[i]+1) <= max_part):
            child = v.clone()
            child._list[i] = v._list[i]+1
            child.set_immutable()
            all_children.append(child)
    return all_children

cpdef int lex_cmp_partial(ClonableIntArray v1, ClonableIntArray v2, int step) noexcept:
    r"""
    Partial comparison of the two lists according the lexicographic
    order. It compares the ``step``-th first entries.

    EXAMPLES::

        sage: from sage.combinat.enumeration_mod_permgroup import lex_cmp_partial
        sage: from sage.structure.list_clone_demo import IncreasingIntArrays
        sage: IA = IncreasingIntArrays()
        sage: lex_cmp_partial(IA([0,1,2,3]),IA([0,1,2,4]),3)
        0
        sage: lex_cmp_partial(IA([0,1,2,3]),IA([0,1,2,4]),4)
        -1
    """
    cdef int i
    if step < 0 or step > v1._len or step > v2._len:
        raise IndexError("list index out of range")

    for i in range(step):
        if v1._list[i] != v2._list[i]:
            break
    if i<step:
        if v1._list[i] > v2._list[i]:
            return 1
        if v1._list[i] < v2._list[i]:
            return -1
    return 0

cpdef int lex_cmp(ClonableIntArray v1, ClonableIntArray v2) noexcept:
    """
    Lexicographic comparison of :class:`~sage.structure.list_clone.ClonableIntArray`.

    INPUT:

    - ``v1``, ``v2`` -- two instances of
      :class:`~sage.structure.list_clone.ClonableIntArray`

    OUTPUT:

    ``-1,0,1``, depending on whether `v_1` is lexicographically smaller, equal, or
    greater than `v_2`.

    EXAMPLES::

        sage: I =  IntegerVectorsModPermutationGroup(SymmetricGroup(5),5)
        sage: I =  IntegerVectorsModPermutationGroup(SymmetricGroup(5))
        sage: J =  IntegerVectorsModPermutationGroup(SymmetricGroup(6))
        sage: v1 = I([2,3,1,2,3], check=False)
        sage: v2 = I([2,3,2,3,2], check=False)
        sage: v3 = J([2,3,1,2,3,1], check=False)
        sage: from sage.combinat.enumeration_mod_permgroup import lex_cmp
        sage: lex_cmp(v1, v1)
        0
        sage: lex_cmp(v1, v2)
        -1
        sage: lex_cmp(v2, v1)
        1
        sage: lex_cmp(v1, v3)
        -1
        sage: lex_cmp(v3, v1)
        1
    """
    cdef int i
    cdef int step = min(v1._len,v2._len)
    for i in range(step):
        if v1._list[i] != v2._list[i]:
            break
    if i<step:
        if v1._list[i] > v2._list[i]:
            return 1
        if v1._list[i] < v2._list[i]:
            return -1
    if v1._len==v2._len:
        return 0
    if v1._len<v2._len:
        return -1
    return 1


cpdef bint is_canonical(list sgs, ClonableIntArray v) except -1:
    r"""
    Return ``True`` if the integer vector `v` is maximal with respect to
    the lexicographic order in its orbit under the action of the
    permutation group whose strong generating system is ``sgs``. Such
    vectors are said to be canonical.

    Let `G` to be the permutation group which admit ``sgs`` as a strong
    generating system.  An integer vector `v` is said to be
    canonical under the action of `G` if and only if:

    .. MATH::

        v = \max_{\text{lex order}} \{g \cdot v | g \in G \}

    EXAMPLES::

        sage: from sage.combinat.enumeration_mod_permgroup import is_canonical
        sage: G = PermutationGroup([[(1,2,3,4)]])
        sage: sgs = G.strong_generating_system()
        sage: from sage.structure.list_clone_demo import IncreasingIntArrays
        sage: IA = IncreasingIntArrays()
        sage: is_canonical(sgs, IA([1,2,3,6]))
        False
    """
    cdef int l, i, comp
    cdef set to_analyse, new_to_analyse
    cdef ClonableIntArray list_test, child
    cdef PermutationGroupElement x
    cdef list transversal
    l = len(v)
    to_analyse = set([v])
    for i in range(l-1):
        new_to_analyse = set()
        transversal = sgs[i]
        for list_test in to_analyse:
            for x in transversal:
                child = x._act_on_array_on_position(list_test)
                comp = lex_cmp_partial(v,child,i+1)
                if comp == -1:
                    return False
                if comp == 0:
                    new_to_analyse.add(child)
        to_analyse = new_to_analyse
    return True


cpdef ClonableIntArray canonical_representative_of_orbit_of(list sgs, ClonableIntArray v):
    r"""
    Return the maximal vector for the lexicographic order living in
    the orbit of `v` under the action of the permutation group whose
    strong generating system is ``sgs``. The maximal vector is also
    called "canonical". Hence, this method returns the canonical
    vector inside the orbit of `v`. If `v` is already canonical,
    the method returns `v`.

    Let `G` to be the permutation group which admits ``sgs`` as a strong
    generating system.  An integer vector `v` is said to be
    canonical under the action of `G` if and only if:

    .. MATH::

        v = \max_{\text{lex order}} \{g \cdot v | g \in G \}

    EXAMPLES::

        sage: from sage.combinat.enumeration_mod_permgroup import canonical_representative_of_orbit_of
        sage: G = PermutationGroup([[(1,2,3,4)]])
        sage: sgs = G.strong_generating_system()
        sage: from sage.structure.list_clone_demo import IncreasingIntArrays
        sage: IA = IncreasingIntArrays()
        sage: canonical_representative_of_orbit_of(sgs, IA([1,2,3,5]))
        [5, 1, 2, 3]
    """
    cdef int l,i,comp
    cdef set to_analyse, new_to_analyse
    cdef ClonableIntArray representative, list_test, child
    cdef PermutationGroupElement x
    representative = v
    l = len(v)
    to_analyse = set([v])
    for i in range(l-1):
        new_to_analyse = set()
        for list_test in to_analyse:
            for x in sgs[i]:
                child = x._act_on_array_on_position(list_test)
                comp = lex_cmp_partial(representative,child,i+1)
                if comp <= 0:
                    new_to_analyse.add(child)
        to_analyse = new_to_analyse
        representative = max(to_analyse)
    return representative

cpdef list canonical_children(list sgs, ClonableIntArray v, int max_part):
    r"""
    Return the canonical children of the integer vector ``v``. This
    function computes all children of the integer vector ``v`` via the
    function :func:`all_children` and returns from this list only
    these which are canonicals identified via the function
    :func:`is_canonical`.

    EXAMPLES::

        sage: from sage.combinat.enumeration_mod_permgroup import canonical_children
        sage: G = PermutationGroup([[(1,2,3,4)]])
        sage: sgs = G.strong_generating_system()
        sage: from sage.structure.list_clone_demo import IncreasingIntArrays
        sage: IA = IncreasingIntArrays()
        sage: canonical_children(sgs, IA([1,2,3,5]), -1)
        []
    """
    cdef ClonableIntArray child
    return [child for child in all_children(v, max_part) if is_canonical(sgs, child)]

cpdef set orbit(list sgs, ClonableIntArray v):
    r"""
    Return the orbit of the integer vector ``v`` under the action of the
    permutation group whose strong generating system is ``sgs``.

    NOTE:

    The returned orbit is a set. In the doctests, we convert it into a
    sorted list.

    EXAMPLES::

        sage: from sage.combinat.enumeration_mod_permgroup import orbit
        sage: G = PermutationGroup([[(1,2,3,4)]])
        sage: sgs = G.strong_generating_system()
        sage: from sage.structure.list_clone_demo import IncreasingIntArrays
        sage: IA = IncreasingIntArrays()
        sage: sorted(orbit(sgs, IA([1,2,3,4])))
        [[1, 2, 3, 4], [2, 3, 4, 1], [3, 4, 1, 2], [4, 1, 2, 3]]
    """
    cdef i,l
    cdef set to_analyse, new_to_analyse
    cdef ClonableIntArray list_test, child
    cdef PermutationGroupElement x
    l = len(v)
    to_analyse = set([v])
    for i in range(l-1):
        new_to_analyse = set()
        for list_test in to_analyse:
            for x in sgs[i]:
                child = x._act_on_array_on_position(list_test)
                new_to_analyse.add(child)
        to_analyse = new_to_analyse
    return to_analyse
