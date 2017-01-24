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

.. [BKSTY06] \A. Buch, A. Kresch, M. Shimozono, H. Tamvakis, and A. Yong.
   *Stable Grothendieck polynomials and* `K`-*theoretic factor sequences*.
   Math. Ann. **340** Issue 2, (2008), pp. 359--382.
   :arxiv:`math/0601514v1`.
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

from builtins import zip

from sage.matrix.matrix import is_Matrix
from sage.matrix.all import matrix


def RSK(obj1=None, obj2=None, insertion='RSK', check_standard=False, **options):
    r"""
    Perform the Robinson-Schensted-Knuth (RSK) correspondence.

    The Robinson-Schensted-Knuth (RSK) correspondence (also known
    as the RSK algorithm) is most naturally stated as a bijection
    between generalized permutations (also known as two-line arrays,
    biwords, ...) and pairs of semi-standard Young tableaux `(P, Q)`
    of identical shape. The tableau `P` is known as the insertion
    tableau, and `Q` is known as the recording tableau.

    The basic operation is known as row insertion `P \leftarrow k`
    (where `P` is a given semi-standard Young tableau, and `k` is an
    integer). Row insertion is a recursive algorithm which starts by
    setting `k_0 = k`, and in its `i`-th step inserts the number `k_i`
    into the `i`-th row of `P` (we start counting the rows at `0`) by
    replacing the first integer greater than `k_i` in the row by `k_i`
    and defines `k_{i+1}` as the integer that has been replaced. If no
    integer greater than `k_i` exists in the `i`-th row, then `k_i` is
    simply appended to the row and the algorithm terminates at this
    point.

    Now the RSK algorithm, applied to a generalized permutation
    `p = ((j_0, k_0), (j_1, k_1), \ldots, (j_{\ell-1}, k_{\ell-1}))`
    (encoded as a lexicographically sorted list of pairs) starts by
    initializing two semi-standard tableaux `P_0` and `Q_0` as empty
    tableaux. For each nonnegative integer `t` starting at `0`, take
    the pair `(j_t, k_t)` from `p` and set
    `P_{t+1} = P_t \leftarrow k_t`, and define `Q_{t+1}` by adding a
    new box filled with `j_t` to the tableau `Q_t` at the same
    location the row insertion on `P_t` ended (that is to say, adding
    a new box with entry `j_t` such that `P_{t+1}` and `Q_{t+1}` have
    the same shape). The iterative process stops when `t` reaches the
    size of `p`, and the pair `(P_t, Q_t)` at this point is the image
    of `p` under the Robinson-Schensted-Knuth correspondence.

    This correspondence has been introduced in [Knu1970]_, where it has
    been referred to as "Construction A".

    For more information, see Chapter 7 in [Sta-EC2]_.

    We also note that integer matrices are in bijection with generalized
    permutations. Furthermore, we can convert any word `w` (and, in
    particular, any permutation) to a generalized permutation by
    considering the top line to be `(1, 2, \ldots, n)` where `n` is the
    length of `w`.

    The optional argument ``insertion`` allows to specify an alternative
    insertion procedure to be used instead of the standard
    Robinson-Schensted-Knuth insertion. If the input is a reduced word of
    a permutation (i.e., an element of a type-`A` Coxeter group), one can
    set ``insertion`` to ``'EG'``, which gives Edelman-Greene insertion,
    an algorithm defined in [EG1987]_ Definition 6.20 (where it is
    referred to as Coxeter-Knuth insertion). The Edelman-Greene insertion
    is similar to the standard row insertion except that if `k_i` and
    `k_i + 1` both exist in row `i`, we *only* set `k_{i+1} = k_i + 1` and
    continue.

    One can also perform a "Hecke RSK algorithm", defined using the
    Hecke insertion studied in [BKSTY06]_ (but using rows instead of
    columns). The algorithm proceeds similarly to the classical RSK
    algorithm. However, it is not clear in what generality it works;
    thus, following [BKSTY06]_, we shall assume that our biword `p`
    has top line `(1, 2, \ldots, n)` (or, at least, has its top line
    strictly increasing). The Hecke RSK algorithm returns a pair of
    an increasing tableau and a set-valued standard tableau. If
    `p = ((j_0, k_0), (j_1, k_1), \ldots, (j_{\ell-1}, k_{\ell-1}))`,
    then the algorithm recursively constructs pairs
    `(P_0, Q_0), (P_1, Q_1), \ldots, (P_\ell, Q_\ell)` of tableaux.
    The construction of `P_{t+1}` and `Q_{t+1}` from `P_t`, `Q_t`,
    `j_t` and `k_t` proceeds as follows: Set `i = j_t`, `x = k_t`,
    `P = P_t` and `Q = Q_t`. We are going to insert `x` into the
    increasing tableau `P` and update the set-valued "recording
    tableau" `Q` accordingly. As in the classical RSK algorithm, we
    first insert `x` into row `1` of `P`, then into row `2` of the
    resulting tableau, and so on, until the construction terminates.
    The details are different: Suppose we are inserting `x` into
    row `R` of `P`. If (Case 1) there exists an entry `y` in row `R`
    such that `x < y`, then let `y` be the minimal such entry. We
    replace this entry `y` with `x` if the result is still an
    increasing tableau; in either subcase, we then continue
    recursively, inserting `y` into the next row of `P`.
    If, on the other hand, (Case 2) no such `y` exists, then we
    append `x` to the end of `R` if the result is an increasing
    tableau (Subcase 2.1), and otherwise (Subcase 2.2) do nothing.
    Furthermore, in Subcase 2.1, we add the box that we have just
    filled with `x` in `P` to the shape of `Q`, and fill it with
    the one-element set `\{i\}`. In Subcase 2.2, we find the
    bottommost box of the column containing the rightmost box of
    row `R`, and add `i` to the entry of `Q` in this box (this
    entry is a set, since `Q` is a set-valued). In either
    subcase, we terminate the recursion, and set
    `P_{t+1} = P` and `Q_{t+1} = Q`.

    Notice that set-valued tableaux are encoded as tableaux whose
    entries are tuples of positive integers; each such tuple is strictly
    increasing and encodes a set (namely, the set of its entries).

    INPUT:

    - ``obj1, obj2`` -- Can be one of the following:

      - A word in an ordered alphabet
      - An integer matrix
      - Two lists of equal length representing a generalized permutation
      - Any object which has a method ``_rsk_iter()`` which returns an
        iterator over the object represented as generalized permutation or
        a pair of lists.

    - ``insertion`` -- (Default: ``'RSK'``) The following types of insertion
      are currently supported:

      - ``'RSK'`` -- Robinson-Schensted-Knuth
      - ``'EG'`` -- Edelman-Greene (only for reduced words of
        permutations/elements of a type-`A` Coxeter group)
      - ``'hecke'`` -- Hecke insertion (only guaranteed for
        generalized permutations whose top row is strictly increasing)

    - ``check_standard`` -- (Default: ``False``) Check if either of the
      resulting tableaux is a standard tableau, and if so, typecast it
      as such

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

    There are also variations of the insertion algorithm in RSK.
    Here we consider Edelman-Greene insertion::

        sage: RSK([2,1,2,3,2], insertion='EG')
        [[[1, 2, 3], [2, 3]], [[1, 3, 4], [2, 5]]]

    We reproduce figure 6.4 in [EG1987]_::

        sage: RSK([2,3,2,1,2,3], insertion='EG')
        [[[1, 2, 3], [2, 3], [3]], [[1, 2, 6], [3, 5], [4]]]

    Hecke insertion is also supported. We construct Example 2.1
    in :arxiv:`0801.1319v2`::

        sage: w = [5, 4, 1, 3, 4, 2, 5, 1, 2, 1, 4, 2, 4]
        sage: RSK(w, insertion='hecke')
        [[[1, 2, 4, 5], [2, 4, 5], [3, 5], [4], [5]],
         [[(1,), (4,), (5,), (7,)],
          [(2,), (9,), (11, 13)],
          [(3,), (12,)],
          [(6,)],
          [(8, 10)]]]

    There is also :func:`~sage.combinat.rsk.RSK_inverse` which performs
    the inverse of the bijection on a pair of semistandard tableaux. We
    note that the inverse function takes 2 separate tableaux as inputs, so
    to compose with :func:`~sage.combinat.rsk.RSK`, we need to use the
    python ``*`` on the output::

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
        sage: RSK(Word([]), insertion='hecke')
        [[], []]
    """
    from sage.combinat.tableau import SemistandardTableau, StandardTableau

    if insertion == 'hecke':
        return hecke_insertion(obj1, obj2)

    if obj1 is None and obj2 is None:
        if 'matrix' in options:
            obj1 = matrix(options['matrix'])
        else:
            raise ValueError("invalid input")

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
                itr = zip(t, b)
            except TypeError:
                itr = zip(range(1, len(obj1)+1), obj1)
    else:
        if len(obj1) != len(obj2):
            raise ValueError("the two arrays must be the same length")
        # Check it is a generalized permutation
        lt = 0
        lb = 0
        for t,b in zip(obj1, obj2):
            if t < lt or (lt == t and b < lb):
                raise ValueError("invalid generalized permutation")
            lt = t
            lb = b
        itr = zip(obj1, obj2)

    from bisect import bisect_right
    p = []       #the "insertion" tableau
    q = []       #the "recording" tableau

    use_EG = (insertion == 'EG')

    #For each x in self, insert x into the tableau p.
    lt = 0
    lb = 0
    for i, x in itr:
        for r, qr in zip(p,q):
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

robinson_schensted_knuth = RSK

def RSK_inverse(p, q, output='array', insertion='RSK'):
    r"""
    Return the generalized permutation corresponding to the pair of
    tableaux `(p,q)` under the inverse of the Robinson-Schensted-Knuth
    algorithm.

    For more information on the bijection, see :func:`RSK`.

    INPUT:

    - ``p``, ``q`` -- Two semi-standard tableaux of the same shape, or
      (in the case when Hecke insertion is used) an increasing tableau and
      a set-valued tableau of the same shape (see the note below for the
      format of the set-valued tableau)

    - ``output`` -- (Default: ``'array'``) if ``q`` is semi-standard:

      - ``'array'`` -- as a two-line array (i.e. generalized permutation or
        biword)
      -  ``'matrix'`` -- as an integer matrix

      and if ``q`` is standard, we can have the output:

      - ``'word'`` -- as a word

      and additionally if ``p`` is standard, we can also have the output:

      - ``'permutation'`` -- as a permutation

    - ``insertion`` -- (Default: ``RSK``) The insertion algorithm used in the
      bijection. Currently the following are supported:

      - ``'RSK'`` -- Robinson-Schensted-Knuth insertion
      - ``'EG'`` -- Edelman-Greene insertion
      - ``'hecke'`` -- Hecke insertion

    .. NOTE::

        In the case of Hecke insertion, the input variable ``q`` should
        be a set-valued tableau, encoded as a tableau whose entries are
        strictly increasing tuples of positive integers. Each such tuple
        encodes the set of its entries.

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

    Using Edelman-Greene insertion::

        sage: pq = RSK([2,1,2,3,2], insertion='EG'); pq
        [[[1, 2, 3], [2, 3]], [[1, 3, 4], [2, 5]]]
        sage: RSK_inverse(*pq, insertion='EG')
        [2, 1, 2, 3, 2]

    Using Hecke insertion::

        sage: w = [5, 4, 1, 3, 4, 2, 5, 1, 2, 1, 4, 2, 4]
        sage: pq = RSK(w, insertion='hecke')
        sage: RSK_inverse(*pq, insertion='hecke', output='list')
        [5, 4, 1, 3, 4, 2, 5, 1, 2, 1, 4, 2, 4]

    .. NOTE::

        The constructor of ``Tableau`` accepts not only semistandard
        tableaux, but also arbitrary lists that are fillings of a
        partition diagram. (And such lists are used, e.g., for the
        set-valued tableau ``q`` that is passed to
        ``RSK_inverse(p, q, insertion='hecke')``.)
        The user is responsible for ensuring that the tableaux passed to
        ``RSK_inverse`` are of the right types (semistandard, standard,
        increasing, set-valued as needed).

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

    Check that :trac:`20430` is fixed::

        sage: RSK([1,1,1,1,1,1,1,2,2,2,3], [1,1,1,1,1,1,3,2,2,2,1])
        [[[1, 1, 1, 1, 1, 1, 1, 2, 2], [2], [3]],
         [[1, 1, 1, 1, 1, 1, 1, 2, 2], [2], [3]]]
        sage: t = SemistandardTableau([[1, 1, 1, 1, 1, 1, 1, 2, 2], [2], [3]])
        sage: RSK_inverse(t, t, 'array')
        [[1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3],
         [1, 1, 1, 1, 1, 1, 3, 2, 2, 2, 1]]
    """
    if insertion == 'hecke':
        return hecke_insertion_reverse(p, q, output)

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

        for key in sorted(d, reverse=True): # Delete last entry from i-th row of p_copy
            i = d[key]
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
        raise ValueError("invalid output option")

    # Checks
    if insertion != 'RSK':
        raise NotImplementedError("only RSK is implemented for non-standard q")
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
    for value, row_dict in sorted(d.items(), reverse=True, key=lambda x: x[0]):
        for key in sorted(row_dict, reverse=True):
            i = row_dict[key]
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
    raise ValueError("invalid output option")

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

#####################################################################
## Hecke insertion

def hecke_insertion(obj1, obj2=None):
    """
    Return the Hecke insertion of the pair ``[obj1, obj2]``.

    .. SEEALSO::

        :func:`RSK`

    EXAMPLES::

        sage: w = [5, 4, 1, 3, 4, 2, 5, 1, 2, 1, 4, 2, 4]
        sage: RSK(w, insertion='hecke')
        [[[1, 2, 4, 5], [2, 4, 5], [3, 5], [4], [5]],
         [[(1,), (4,), (5,), (7,)],
          [(2,), (9,), (11, 13)],
          [(3,), (12,)],
          [(6,)],
          [(8, 10)]]]
    """
    if obj2 is None:
        obj2 = obj1
        obj1 = list(range(1, len(obj2) + 1))

    from sage.combinat.tableau import SemistandardTableau, Tableau
    from bisect import bisect_right
    p = []       #the "insertion" tableau
    q = []       #the "recording" tableau

    for i, x in zip(obj1, obj2):
        for j,r in enumerate(p):
            if r[-1] > x:
                #Figure out where to insert x into the row r.  The
                #bisect command returns the position of the least
                #element of r greater than x.  We will call it y.
                y_pos = bisect_right(r, x)
                y = r[y_pos]
                # Check to see if we can swap x for y
                if (y_pos == 0 or r[y_pos-1] < x) and (j == 0 or p[j-1][y_pos] < x):
                    r[y_pos] = x
                x = y
            else:
                # We must have len(p[j-1]) > len(r), since x is coming
                # from the previous row.
                if r[-1] < x and (j == 0 or p[j-1][len(r)] < x):
                    # We can add a box to the row
                    r.append(x)
                    q[j].append((i,)) # Values are always inserted to the right
                else:
                    # We must append i to the bottom of this column
                    l = len(r) - 1
                    while j < len(q) and len(q[j]) > l:
                        j += 1
                    q[j-1][-1] = q[j-1][-1] + (i,)
                break
        else:
            #We made through all of the rows of p without breaking
            #so we need to add a new row to p and q.
            p.append([x])
            q.append([(i,)])

    return [SemistandardTableau(p), Tableau(q)]

def hecke_insertion_reverse(p, q, output='array'):
    r"""
    Return the reverse Hecke insertion of ``(p, q)``.

    .. SEEALSO::

        :func:`RSK_inverse`

    EXAMPLES::

        sage: w = [5, 4, 1, 3, 4, 2, 5, 1, 2, 1, 4, 2, 4]
        sage: P,Q = RSK(w, insertion='hecke')
        sage: wp = RSK_inverse(P, Q, insertion='hecke', output='list'); wp
        [5, 4, 1, 3, 4, 2, 5, 1, 2, 1, 4, 2, 4]
        sage: wp == w
        True
    """
    if p.shape() != q.shape():
        raise ValueError("p(=%s) and q(=%s) must have the same shape"%(p, q))
    from sage.combinat.tableau import SemistandardTableaux
    if p not in SemistandardTableaux():
        raise ValueError("p(=%s) must be a semistandard tableau"%p)

    from bisect import bisect_left
    # Make a copy of p and q since this is destructive to it
    p_copy = [list(row) for row in p]
    q_copy = [[list(v) for v in row] for row in q]
    # We shall work on these copies of p and q. Notice that p might get
    # some empty rows in the process; we do not bother pruning them, as
    # they do not matter.

    #upper_row and lower_row will be the upper and lower rows of the
    #generalized permutation we get as a result, but both reversed.
    upper_row = []
    lower_row = []
    d = {}
    for ri, row in enumerate(q):
        for ci, entry in enumerate(row):
            for val in entry:
                if val in d:
                    d[val][ci] = ri
                else:
                    d[val] = {ci: ri}
    #d is now a double family such that for every integers k and j,
    #the value d[k][j] is the row i such that the (i, j)-th cell of
    #q is filled with k.
    for value, row_dict in sorted(d.items(), key=lambda x: -x[0]):
        for i in sorted(row_dict.values(), reverse=True):
            # These are always the right-most entry
            should_be_value = q_copy[i][-1].pop()
            assert value == should_be_value
            if not q_copy[i][-1]:
                # That is, if value was alone in cell q_copy[i][-1].
                q_copy[i].pop()
                x = p_copy[i].pop()
            else:
                x = p_copy[i][-1]
            while i > 0:
                i -= 1
                row = p_copy[i]
                y_pos = bisect_left(row,x) - 1
                y = row[y_pos]
                # Check to see if we can swap x for y
                if ((y_pos == len(row) - 1 or x < row[y_pos+1])
                    and (i == len(p_copy) - 1 or len(p_copy[i+1]) <= y_pos
                         or x < p_copy[i+1][y_pos])):
                    row[y_pos] = x
                x = y
            upper_row.append(value)
            lower_row.append(x)

    if output == 'array':
        return [list(reversed(upper_row)), list(reversed(lower_row))]
    is_standard = (upper_row == list(range(len(upper_row), 0, -1)))
    if output == 'word':
        if not is_standard:
            raise TypeError("q must be standard to have a %s as valid output"%output)
        from sage.combinat.words.word import Word
        return Word(reversed(lower_row))
    if output == 'list':
        if not is_standard:
            raise TypeError("q must be standard to have a %s as valid output"%output)
        return list(reversed(lower_row))
    raise ValueError("invalid output option")

