r"""
Robinson-Schensted-Knuth correspondence

AUTHORS:

- Travis Scrimshaw (2012-12-07): Initial version

EXAMPLES:

We can perform RSK and the inverse on a variety of objects::

    sage: p = Tableau([[1,2,2],[2]]); q = Tableau([[1,3,3],[2]])
    sage: gp = RSK_inverse(p, q); gp
    [[1, 2, 3, 3], [2, 1, 2, 2]]
    sage: RSK(*gp)
    [[[1, 2, 2], [2]], [[1, 3, 3], [2]]]
    sage: m = RSK_inverse(p, q, 'matrix'); m
    [0 1]
    [1 0]
    [0 2]
    sage: RSK(m)
    [[[1, 2, 2], [2]], [[1, 3, 3], [2]]]

TESTS:

Check that it is a correspondence between all types of input and
the input is preserved::

    sage: t1 = Tableau([[1, 2, 5], [3], [4]])
    sage: t2 = Tableau([[1, 2, 3], [4], [5]])
    sage: gp = RSK_inverse(t1, t2); gp
    [[1, 2, 3, 4, 5], [1, 4, 5, 3, 2]]
    sage: w = RSK_inverse(t1, t2, 'word'); w
    word: 14532
    sage: m = RSK_inverse(t1, t2, 'matrix'); m
    [1 0 0 0 0]
    [0 0 0 1 0]
    [0 0 0 0 1]
    [0 0 1 0 0]
    [0 1 0 0 0]
    sage: p = RSK_inverse(t1, t2, 'permutation'); p
    [1, 4, 5, 3, 2]
    sage: t1
    [[1, 2, 5], [3], [4]]
    sage: t2
    [[1, 2, 3], [4], [5]]
    sage: RSK(*gp) == [t1, t2]
    True
    sage: RSK(w) == [t1, t2]
    True
    sage: RSK(m) == [t1, t2]
    True
    sage: RSK(p) == [t1, t2]
    True
    sage: gp
    [[1, 2, 3, 4, 5], [1, 4, 5, 3, 2]]
    sage: w
    word: 14532
    sage: m
    [1 0 0 0 0]
    [0 0 0 1 0]
    [0 0 0 0 1]
    [0 0 1 0 0]
    [0 1 0 0 0]
    sage: p
    [1, 4, 5, 3, 2]

REFERENCES:

.. [Knu1970] Donald E. Knuth.
   *Permutations, matrices, and generalized Young tableaux*.
   Pacific J. Math. Volume 34, Number 3 (1970), pp. 709-727.
   http://projecteuclid.org/euclid.pjm/1102971948

.. [EG1987] Paul Edelman, Curtis Greene.
   *Balanced Tableaux*.
   Advances in Mathematics 63 (1987), pp. 42-99.
   http://www.sciencedirect.com/science/article/pii/0001870887900636
"""
#*****************************************************************************
#       Copyright (C) 2012 Travis Scrimshaw <tscrim@ucdavis.edu>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.matrix import is_Matrix
from sage.matrix.all import matrix
from itertools import izip
from sage.misc.superseded import deprecated_function_alias

def RSK(obj1=None, obj2=None, insertion='RSK', check_standard=False, **options):
    r"""
    Perform the Robinson-Schensted-Knuth (RSK) correspondence.

    The Robinson-Schensted-Knuth (RSK) correspondence is most naturally stated
    as a bijection between generalized permutations (also known as two-line
    arrays, biwords, ...) and pairs of semi-standard Young tableaux `(P, Q)`
    of identical shape. The tableau `P` is known as the insertion tableau and
    `Q` is known as the recording tableau.

    The basic operation is known as row insertion `P \leftarrow k` (where `P`
    is a given semi-standard Young tableau, and `k` is an integer). Row
    insertion is a recursive algorithm which starts by setting `k_0 = k`,
    and in its `i`-th step inserts the number `k_i` into the `i`-th row of
    `P` by replacing the first integer greater than `k_i` in the row by `k_i`
    and defines `k_{i+1}` as the integer that has been replaced. If no integer
    greater than `k_i` exists in the `i`-th row, then `k_i` is simply appended
    to the row and the algorithm terminates at this point.

    Now the RSK algorithm starts by initializing two semi-standard tableaux
    `P_0` and `Q_0` as empty tableaux. For each nonnegative integer `t`
    starting at `0`, take the pair `(j_t, k_t)` from `p` and set
    `P_{t+1} = P_t \leftarrow k_t`, and define `Q_{t+1}` by adding a new box
    filled with `j_t` to the tableau `Q_t` at the same location the row
    insertion on `P_t` ended (that is to say, adding `j_t` such that
    `P_{t+1}` and `Q_{t+1}` have the same shape). The iterative process stops
    when `t` reaches the size of `p`, and the pair `(P_t, Q_t)` at this point
    is the image of `p` under the Robinson-Schensted-Knuth correspondence.

    This correspondence has been introduced in [Knu1970]_, where it has been
    referred to as "Construction A".

    For more information, see Chapter 7 in [Sta-EC2]_.

    We also note that integer matrices are in bijection with generalized
    permutations. In addition, we can also convert any word `w` (and any
    permutation) to a generalized permutation by considering the top line
    to be `(1, 2, \ldots, n)` where `n` is the length of `w`.

    The optional argument ``insertion`` allows to specify an alternative
    insertion procedure to be used instead of the standard
    Robinson-Schensted-Knuth insertion. If the input is a reduced word of
    a permutation (i.e., a type `A` Coxeter element), one can set
    ``insertion`` to ``'EG'``, which gives Edelman-Greene insertion, an
    algorithm defined in [EG1987]_ Definition 6.20 (where it is referred
    to as Coxeter-Knuth insertion). The Edelman-Greene insertion is similar
    to the standard row insertion except that if `k_i` and `k_i + 1` both
    exist in row `i`, we *only* set `k_{i+1} = k_i + 1` and continue.

    INPUT:

    - ``obj1, obj2`` -- Can be one of the following:

      - A word in an ordered alphabet
      - An integer matrix
      - Two lists of equal length representing a generalized permutation
      - Any object which has a method ``_rsk_iter()`` which returns an iterator
        over the object represented as generalized permutation or a pair of
        lists.

    - ``insertion`` -- (Default: ``'RSK'``) The following types of insertion
      are currently supported:

      - ``'RSK'`` -- Robinson-Schensted-Knuth
      - ``'EG'`` -- Edelman-Greene (only for reduced words of
        permutations/type `A` Coxeter elements)

    - ``check_standard`` -- (Default: ``False``) Check if either of the resulting
      tableaux should be standard tableau.

    EXAMPLES:

    If we only give one line, we treat the top line as being
    `(1, 2, \ldots, n)`::

        sage: RSK([3,3,2,4,1])
        [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]
        sage: RSK(Word([3,3,2,4,1]))
        [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]
        sage: RSK(Word([2,3,3,2,1,3,2,3]))
        [[[1, 2, 2, 3, 3], [2, 3], [3]], [[1, 2, 3, 6, 8], [4, 7], [5]]]

    With a generalized permutation::

        sage: RSK([1, 2, 2, 2], [2, 1, 1, 2])
        [[[1, 1, 2], [2]], [[1, 2, 2], [2]]]
        sage: RSK(Word([1,1,3,4,4]), [1,4,2,1,3])
        [[[1, 1, 3], [2], [4]], [[1, 1, 4], [3], [4]]]
        sage: RSK([1,3,3,4,4], Word([6,2,2,1,7]))
        [[[1, 2, 7], [2], [6]], [[1, 3, 4], [3], [4]]]

    If we give it a matrix::

        sage: RSK(matrix([[0,1],[2,1]]))
        [[[1, 1, 2], [2]], [[1, 2, 2], [2]]]

    We can also give it something looking like a matrix::

        sage: RSK([[0,1],[2,1]])
        [[[1, 1, 2], [2]], [[1, 2, 2], [2]]]

    There are also variations of the insertion algorithm in RSK, and currently
    only Edelman-Greene insertion is supported::

        sage: RSK([2,1,2,3,2], insertion='EG')
        [[[1, 2, 3], [2, 3]], [[1, 3, 4], [2, 5]]]

    We reproduce figure 6.4 in [EG1987]_::

        sage: RSK([2,3,2,1,2,3], insertion='EG')
        [[[1, 2, 3], [2, 3], [3]], [[1, 2, 6], [3, 5], [4]]]

    There is also :func:`~sage.combinat.rsk.RSK_inverse` which performs the
    inverse of the bijection on a pair of semistandard tableaux. We note
    that the inverse function takes 2 separate tableaux inputs, so to compose
    with :func:`~sage.combinat.rsk.RSK`, we need to use the python ``*`` on
    the output::

        sage: RSK_inverse(*RSK([1, 2, 2, 2], [2, 1, 1, 2]))
        [[1, 2, 2, 2], [2, 1, 1, 2]]
        sage: P,Q = RSK([1, 2, 2, 2], [2, 1, 1, 2])
        sage: RSK_inverse(P, Q)
        [[1, 2, 2, 2], [2, 1, 1, 2]]

    TESTS:

    Empty objects::

        sage: RSK(Permutation([]))
        [[], []]
        sage: RSK(Word([]))
        [[], []]
        sage: RSK(matrix([[]]))
        [[], []]
        sage: RSK([], [])
        [[], []]
        sage: RSK([[]])
        [[], []]
        sage: RSK(Word([]), insertion='EG')
        [[], []]
    """
    from sage.combinat.tableau import SemistandardTableau, StandardTableau

    if obj1 is None and obj2 is None:
        if 'matrix' in options:
            obj1 = matrix(options['matrix'])
        else:
            raise ValueError("Invalid input")

    if is_Matrix(obj1):
        obj1 = obj1.rows()
    if len(obj1) == 0:
        return [StandardTableau([]), StandardTableau([])]

    if obj2 is None:
        try:
            itr = obj1._rsk_iter()
        except AttributeError:
            # If this is (something which looks like) a matrix
            #   then build the generalized permutation
            try:
                t = []
                b = []
                for i, row in enumerate(obj1):
                    for j, mult in enumerate(row):
                        if mult > 0:
                            t.extend([i+1]*mult)
                            b.extend([j+1]*mult)
                itr = izip(t, b)
            except TypeError:
                itr = izip(range(1, len(obj1)+1), obj1)
    else:
        if len(obj1) != len(obj2):
            raise ValueError("The two arrays must be the same length")
        # Check it is a generalized permutation
        lt = 0
        lb = 0
        for t,b in izip(obj1, obj2):
            if t < lt or (lt == t and b < lb):
                raise ValueError("Invalid generalized permutation")
            lt = t
            lb = b
        itr = izip(obj1, obj2)

    from bisect import bisect_right
    p = []       #the "insertion" tableau
    q = []       #the "recording" tableau

    use_EG = (insertion == 'EG')

    #For each x in self, insert x into the tableau p.
    lt = 0
    lb = 0
    for i, x in itr:
        for r, qr in izip(p,q):
            if r[-1] > x:
                #Figure out where to insert x into the row r.  The
                #bisect command returns the position of the least
                #element of r greater than x.  We will call it y.
                y_pos = bisect_right(r, x)
                if use_EG and r[y_pos] == x + 1 and y_pos > 0 and x == r[y_pos - 1]:
                    #Special bump: Nothing to do except increment x by 1
                    x += 1
                else:
                    #Switch x and y
                    x, r[y_pos] = r[y_pos], x
            else:
                break
        else:
            #We made through all of the rows of p without breaking
            #so we need to add a new row to p and q.
            r = []; p.append(r)
            qr = []; q.append(qr)

        r.append(x)
        qr.append(i) # Values are always inserted to the right

    if check_standard:
        try:
            P = StandardTableau(p)
        except ValueError:
            P = SemistandardTableau(p)
        try:
            Q = StandardTableau(q)
        except ValueError:
            Q = SemistandardTableau(q)
        return [P, Q]
    return [SemistandardTableau(p), SemistandardTableau(q)]

RobinsonSchenstedKnuth = deprecated_function_alias(15142, RSK)
robinson_schensted_knuth = RSK

def RSK_inverse(p, q, output='array', insertion='RSK'):
    r"""
    Return the generalized permutation corresponding to the pair of tableaux
    `(p,q)` under the inverse of the Robinson-Schensted-Knuth algorithm.

    For more information on the bijeciton, see :func:`RSK`.

    INPUT:

    - ``p``, ``q`` -- Two semi-standard tableaux of the same shape

    - ``output`` -- (Default: ``'array'``) if ``q`` is semi-standard:

      - ``'array'`` -- as a two-line array (i.e. generalized permutation or
        biword)
      -  ``'matrix'`` -- as an integer matrix

      and if ``q`` is standard, we can have the output:

      - ``'word'`` -- as a word

      and additionally if ``p`` is standard, we can also have the output:

      - ``'permutation'`` -- as a permutation

    - ``insertion`` -- (Default: ``RSK``) The insertion algorithm used in the
      bijection. Currently only ``'RSK'`` and ``'EG'`` (Edelman-Greene) are
      supported.

    EXAMPLES:

    If both ``p`` and ``q`` are standard::

        sage: t1 = Tableau([[1, 2, 5], [3], [4]])
        sage: t2 = Tableau([[1, 2, 3], [4], [5]])
        sage: RSK_inverse(t1, t2)
        [[1, 2, 3, 4, 5], [1, 4, 5, 3, 2]]
        sage: RSK_inverse(t1, t2, 'word')
        word: 14532
        sage: RSK_inverse(t1, t2, 'matrix')
        [1 0 0 0 0]
        [0 0 0 1 0]
        [0 0 0 0 1]
        [0 0 1 0 0]
        [0 1 0 0 0]
        sage: RSK_inverse(t1, t2, 'permutation')
        [1, 4, 5, 3, 2]
        sage: RSK_inverse(t1, t1, 'permutation')
        [1, 4, 3, 2, 5]
        sage: RSK_inverse(t2, t2, 'permutation')
        [1, 2, 5, 4, 3]
        sage: RSK_inverse(t2, t1, 'permutation')
        [1, 5, 4, 2, 3]

    If the first tableau is semistandard::

        sage: p = Tableau([[1,2,2],[3]]); q = Tableau([[1,2,4],[3]])
        sage: ret = RSK_inverse(p, q); ret
        [[1, 2, 3, 4], [1, 3, 2, 2]]
        sage: RSK_inverse(p, q, 'word')
        word: 1322

    In general::

        sage: p = Tableau([[1,2,2],[2]]); q = Tableau([[1,3,3],[2]])
        sage: RSK_inverse(p, q)
        [[1, 2, 3, 3], [2, 1, 2, 2]]
        sage: RSK_inverse(p, q, 'matrix')
        [0 1]
        [1 0]
        [0 2]

    Using the Edelman-Greene insertion::

        sage: pq = RSK([2,1,2,3,2], insertion='RSK'); pq
        [[[1, 2, 2], [2, 3]], [[1, 3, 4], [2, 5]]]
        sage: RSK_inverse(*pq, insertion='EG')
        [2, 1, 2, 3, 2]

    .. NOTE::

        Currently the constructor of ``Tableau`` accept as input lists
        that are not even tableaux but only filling of a partition diagram.
        This feature should not be used with ``RSK_inverse``.

    TESTS:

    From empty tableaux::

        sage: RSK_inverse(Tableau([]), Tableau([]))
        [[], []]

    Check that :func:`RSK_inverse` is the inverse of :func:`RSK` on the
    different types of inputs/outputs::

        sage: f = lambda p: RSK_inverse(*RSK(p), output='permutation')
        sage: all(p == f(p) for n in range(7) for p in Permutations(n))
        True
        sage: all(RSK_inverse(*RSK(w), output='word') == w for n in range(4) for w in Words(5, n))
        True
        sage: from sage.combinat.integer_matrices import IntegerMatrices
        sage: M = IntegerMatrices([1,2,2,1], [3,1,1,1])
        sage: all(RSK_inverse(*RSK(m), output='matrix') == m for m in M)
        True

        sage: n = ZZ.random_element(200)
        sage: p = Permutations(n).random_element()
        sage: is_fine = True if p == f(p) else p ; is_fine
        True

    Same for Edelman-Greene (but we are checking only the reduced words that
    can be obtained using the ``reduced_word()`` method from permutations)::

        sage: g = lambda w: RSK_inverse(*RSK(w, insertion='EG'), insertion='EG', output='permutation')
        sage: all(p.reduced_word() == g(p.reduced_word()) for n in range(7) for p in Permutations(n))
        True

        sage: n = ZZ.random_element(200)
        sage: p = Permutations(n).random_element()
        sage: is_fine = True if p == f(p) else p ; is_fine
        True

    Both tableaux must be of the same shape::

        sage: RSK_inverse(Tableau([[1,2,3]]), Tableau([[1,2]]))
        Traceback (most recent call last):
        ...
        ValueError: p(=[[1, 2, 3]]) and q(=[[1, 2]]) must have the same shape
    """
    if p.shape() != q.shape():
        raise ValueError("p(=%s) and q(=%s) must have the same shape"%(p, q))
    from sage.combinat.tableau import SemistandardTableaux
    if p not in SemistandardTableaux():
        raise ValueError("p(=%s) must be a semistandard tableau"%p)

    from bisect import bisect_left
    # Make a copy of p since this is destructive to it
    p_copy = [list(row) for row in p]

    if q.is_standard():
        rev_word = [] # This will be our word in reverse
        d = dict((qij,i) for i, Li in enumerate(q) for qij in Li)
        # d is now a dictionary which assigns to each integer k the
        # number of the row of q containing k.

        use_EG = (insertion == 'EG')

        for i in reversed(d.values()): # Delete last entry from i-th row of p_copy
            x = p_copy[i].pop() # Always the right-most entry
            for row in reversed(p_copy[:i]):
                y_pos = bisect_left(row,x) - 1
                if use_EG and row[y_pos] == x - 1 and y_pos < len(row)-1 and row[y_pos+1] == x:
                    # Nothing to do except decrement x by 1.
                    # (Case 1 on p. 74 of Edelman-Greene [EG1987]_.)
                    x -= 1
                else:
                    # switch x and y
                    x, row[y_pos] = row[y_pos], x
            rev_word.append(x)

        if use_EG:
            return list(reversed(rev_word))
        if output == 'word':
            from sage.combinat.words.word import Word
            return Word(reversed(rev_word))
        if output == 'matrix':
            return to_matrix(list(range(1, len(rev_word)+1)), list(reversed(rev_word)))
        if output == 'array':
            return [list(range(1, len(rev_word)+1)), list(reversed(rev_word))]
        if output == 'permutation':
            if not p.is_standard():
                raise TypeError("p must be standard to have a valid permutation as output")
            from sage.combinat.permutation import Permutation
            return Permutation(reversed(rev_word))
        raise ValueError("Invalid output option")

    # Checks
    if insertion != 'RSK':
        raise NotImplementedError("Only RSK is implemented for non-standard Q")
    if q not in SemistandardTableaux():
        raise ValueError("q(=%s) must be a semistandard tableau"%q)

    upper_row = []
    lower_row = []
    #upper_row and lower_row will be the upper and lower rows of the
    #generalized permutation we get as a result, but both reversed.
    d = {}
    for row, Li in enumerate(q):
        for col, val in enumerate(Li):
            if val in d:
                d[val][col] = row
            else:
                d[val] = {col: row}
    #d is now a double family such that for every integers k and j,
    #the value d[k][j] is the row i such that the (i, j)-th cell of
    #q is filled with k.
    for value, row_dict in reversed(d.items()):
        for i in reversed(row_dict.values()):
            x = p_copy[i].pop() # Always the right-most entry
            for row in reversed(p_copy[:i]):
                y = bisect_left(row,x) - 1
                x, row[y] = row[y], x
            upper_row.append(value)
            lower_row.append(x)

    if output == 'matrix':
        return to_matrix(list(reversed(upper_row)), list(reversed(lower_row)))
    if output == 'array':
        return [list(reversed(upper_row)), list(reversed(lower_row))]
    if output in ['permutation', 'word']:
        raise TypeError("q must be standard to have a %s as valid output"%output)
    raise ValueError("Invalid output option")

RobinsonSchenstedKnuth_inverse = deprecated_function_alias(15142, RSK_inverse)
robinson_schensted_knuth_inverse = RSK_inverse

def to_matrix(t, b):
    r"""
    Return the integer matrix corresponding to a two-line array.

    INPUT:

    - ``t`` -- The top line of the array

    - ``b`` -- The bottom line of the array

    OUTPUT:

    An `m \times n`-matrix (where `m` and `n` are the maximum entries in
    `t` and `b` respectively) whose `(i, j)`-th entry, for any `i` and `j`,
    is the number of all positions `k` satisfying `t_k = i` and `b_k = j`.

    EXAMPLES::

        sage: from sage.combinat.rsk import to_matrix
        sage: to_matrix([1, 1, 3, 3, 4], [2, 3, 1, 1, 3])
        [0 1 1]
        [0 0 0]
        [2 0 0]
        [0 0 1]
    """
    n = len(b)
    if len(t) != n:
        raise ValueError("The two arrays must be the same length")

    # Build the dictionary of entries since the matrix
    #   is typically (very) sparse
    entries = {}
    for i in range(n):
        if (t[i]-1, b[i]-1) in entries:
            entries[(t[i]-1, b[i]-1)] += 1
        else:
            entries[(t[i]-1, b[i]-1)] = 1
    return matrix(entries, sparse=True)

