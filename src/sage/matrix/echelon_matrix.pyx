r"""
Echelon matrices over finite fields.
"""
# ****************************************************************************
#       Copyright (C) 2014 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.matrix0 cimport Matrix


def reduced_echelon_matrix_iterator(K, k, n, bint sparse=False, bint copy=True, bint set_immutable=False):
    r"""
    An iterator over `(k,n)` reduced echelon matrices over the finite field `K`.

    INPUT:

    - ``K`` -- a finite field

    - ``k`` -- number of rows (or the size of the subspace)

    - ``n`` -- number of columns (or the dimension of the ambient space)

    - ``sparse`` -- boolean (default: ``False``)

    - ``copy`` -- boolean (default: ``True``); if set to ``False`` then
      iterator yields the same matrix over and over (but with different
      entries).  Default is ``True`` which is safer but might be slower.

    - ``set_immutable`` -- boolean; if set to ``True`` then the output matrices
      are immutable. This option automatically turns ``copy`` into ``True``.


    .. NOTE::

        We ensure that the iteration order is so that all matrices with given
        pivot columns are generated consecutively. Furthermore, the order in
        which the pivot columns appear is lexicographic.

        It would be faster to generate the pivots columns following a Gray code.
        There would be only one pivot changing at a time, avoiding the possibly
        expensive ``m0.__copy__()``. However that would modify the generation
        order some functions depend upon.

    EXAMPLES::

        sage: from sage.matrix.echelon_matrix import reduced_echelon_matrix_iterator
        sage: it = reduced_echelon_matrix_iterator(GF(2), 2, 3)
        sage: for m in it:
        ....:     print(m)
        ....:     print(m.pivots())
        ....:     print("*******")
        [1 0 0]
        [0 1 0]
        (0, 1)
        *******
        [1 0 0]
        [0 1 1]
        (0, 1)
        *******
        [1 0 1]
        [0 1 0]
        (0, 1)
        *******
        [1 0 1]
        [0 1 1]
        (0, 1)
        *******
        [1 0 0]
        [0 0 1]
        (0, 2)
        *******
        [1 1 0]
        [0 0 1]
        (0, 2)
        *******
        [0 1 0]
        [0 0 1]
        (1, 2)
        *******

    TESTS:

    Testing cardinalities::

        sage: q = 71
        sage: F = GF(q)
        sage: len(list(reduced_echelon_matrix_iterator(F, 1, 3, copy=False))) == q**2+q+1
        True
        sage: len(list(reduced_echelon_matrix_iterator(F, 2, 3, copy=False))) == q**2+q+1
        True

    Testing options::

        sage: it = reduced_echelon_matrix_iterator(GF(4, 'z'), 2, 4, copy=False)                    # needs sage.rings.finite_rings
        sage: next(it) is next(it)                                                                  # needs sage.rings.finite_rings
        True
        sage: for a in it: pass

        sage: it = reduced_echelon_matrix_iterator(GF(4, 'z'), 2, 4, set_immutable=True)            # needs sage.rings.finite_rings
        sage: all(a.is_immutable() and a.echelon_form() == a for a in it)
        True
    """
    cdef Matrix m0,m,mm
    cdef int i
    n = int(n)
    k = int(k)

    if n < k:
        raise NotImplementedError("echelon matrix with fewer rows than columns "
                "i.e. not full rank) are not implemented")

    from sage.matrix.matrix_space import MatrixSpace
    from itertools import combinations,product

    copy = copy or set_immutable

    m0 = MatrixSpace(K,k,n,sparse=sparse)()
    Klist = list(K)
    K1 = K.one()

    # First, we select which columns will be pivots:
    for pivots in combinations(range(n),k):
        m = m0.__copy__()
        free_positions = []
        for i in range(k):
            m[i,pivots[i]] = K1
            for j in range(pivots[i]+1,n):
                if j not in pivots:
                    free_positions.append((i,j))

        # Next, we fill in those entries that are not
        # determined by the echelon form alone:
        num_free_pos = len(free_positions)
        for v in product(Klist, repeat=num_free_pos):
            for i in range(num_free_pos):
                m[free_positions[i]] = v[i]

            if copy:
                mm = m.__copy__()
                mm.cache('pivots',pivots)
                mm.cache('rank',k)
                mm.cache('in_echelon_form',True)
                if set_immutable:
                    mm.set_immutable()
                yield mm
            else:
                yield m
            del v   # hack: Python itertools reuses the tuple if nobody else uses it
        del pivots  # hack: Python itertools reuses the tuple if nobody else uses it
