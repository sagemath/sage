# cython: wraparound=False, boundscheck=False
"""
Braid Orbit

Cython function to compute the orbit of the braid moves on a reduced word.
"""
from cysignals.signals cimport sig_check


cpdef set BraidOrbit(list word, list rels):
    r"""
    Return the orbit of ``word`` by all replacements given by ``rels``.

    INPUT:

    - ``word`` -- list of integers

    - ``rels`` -- list of pairs ``(A, B)``, where ``A`` and ``B`` are
      lists of integers the same length

    EXAMPLES::

        sage: from sage.combinat.root_system.braid_orbit import BraidOrbit
        sage: word = [1,2,1,3,2,1]
        sage: rels = [[[2, 1, 2], [1, 2, 1]], [[3, 1], [1, 3]], [[3, 2, 3], [2, 3, 2]]]
        sage: sorted(BraidOrbit(word, rels))
        [(1, 2, 1, 3, 2, 1),
         (1, 2, 3, 1, 2, 1),
         (1, 2, 3, 2, 1, 2),
         (1, 3, 2, 1, 3, 2),
         (1, 3, 2, 3, 1, 2),
         (2, 1, 2, 3, 2, 1),
         (2, 1, 3, 2, 1, 3),
         (2, 1, 3, 2, 3, 1),
         (2, 3, 1, 2, 1, 3),
         (2, 3, 1, 2, 3, 1),
         (2, 3, 2, 1, 2, 3),
         (3, 1, 2, 1, 3, 2),
         (3, 1, 2, 3, 1, 2),
         (3, 2, 1, 2, 3, 2),
         (3, 2, 1, 3, 2, 3),
         (3, 2, 3, 1, 2, 3)]
        sage: len(_)
        16
    """
    cdef int i, l, rel_l, loop_ind, list_len
    cdef tuple left, right, test_word, new_word
    cdef list rel

    l = len(word)
    cdef set words = {tuple(word)}
    cdef list test_words = [tuple(word)]

    rels = rels + [[b, a] for a, b in rels]
    rels = [[tuple(a), tuple(b), len(a)] for a, b in rels]

    loop_ind = 0
    list_len = 1
    while loop_ind < list_len:
        sig_check()
        test_word = <tuple> test_words[loop_ind]
        loop_ind += 1
        for rel in rels:
            left = <tuple> rel[0]
            right = <tuple> rel[1]
            rel_l = <int> rel[2]
            for i in range(l-rel_l+1):
                if pattern_match(test_word, i, left, rel_l):
                    new_word = test_word[:i] + right + test_word[i+rel_l:]
                    if new_word not in words:
                        words.add(new_word)
                        test_words.append(new_word)
                        list_len += 1
    return words


cpdef bint is_fully_commutative(list word, list rels) noexcept:
    r"""
    Check if the braid orbit of ``word`` is using a braid relation.

    INPUT:

    - ``word`` -- list of integers

    - ``rels`` -- list of pairs ``(A, B)``, where ``A`` and ``B`` are
      lists of integers the same length

    EXAMPLES::

        sage: from sage.combinat.root_system.braid_orbit import is_fully_commutative
        sage: rels = [[[2, 1, 2], [1, 2, 1]], [[3, 1], [1, 3]], [[3, 2, 3], [2, 3, 2]]]
        sage: word = [1,2,1,3,2,1]
        sage: is_fully_commutative(word, rels)
        False
        sage: word = [1,2,3]
        sage: is_fully_commutative(word, rels)
        True
    """
    cdef int i, l, rel_l, loop_ind, list_len
    cdef tuple left, right, test_word, new_word
    cdef list rel

    l = len(word)
    cdef set words = {tuple(word)}
    cdef list test_words = [tuple(word)]

    rels = rels + [[b, a] for a, b in rels]
    rels = [[tuple(a), tuple(b), len(a)] for a, b in rels]

    loop_ind = 0
    list_len = 1
    while loop_ind < list_len:
        sig_check()
        test_word = <tuple> test_words[loop_ind]
        loop_ind += 1
        for rel in rels:
            left = <tuple> rel[0]
            right = <tuple> rel[1]
            rel_l = <int> rel[2]
            for i in range(l-rel_l+1):
                if pattern_match(test_word, i, left, rel_l):
                    if rel_l > 2:
                        return False
                    new_word = test_word[:i] + right + test_word[i+rel_l:]
                    if new_word not in words:
                        words.add(new_word)
                        test_words.append(new_word)
                        list_len += 1
    return True


cdef inline bint pattern_match(tuple L, int i, tuple X, int l) noexcept:
    r"""
    Return ``True`` if ``L[i:i+l] == X``.

    Assumes that ``i`` is the length of ``L`` and ``l`` is the
    length of ``X``.
    """
    cdef int ind
    for ind in range(l):
        if <int> L[i+ind] != <int> X[ind]:
            return False
    return True
