r"""
Robinson-Schensted-Knuth correspondence

AUTHORS:

- Travis Scrimshaw (2012-12-07): Initial version
- Chaman Agrawal (2019-06-24): Refactoring on the Rule class

Introduction
============

The Robinson-Schensted-Knuth (RSK) correspondence is most naturally
stated as a bijection between generalized permutations (also known
as two-line arrays, biwords, ...) and pairs of semi-standard Young
tableaux `(P, Q)` of identical shape.

The basic operation in the RSK correspondence is a row insertion
`P \leftarrow k` (where `P` is a given semi-standard Young tableau,
and `k` is an integer). Different insertion algorithms have been
implemented for the RSK correspondence and can be specified as
an argument in the function call.

EXAMPLES:

We can perform RSK and its inverse map on a variety of objects::

    sage: p = Tableau([[1,2,2],[2]]); q = Tableau([[1,3,3],[2]])
    sage: gp = RSK_inverse(p, q); gp
    [[1, 2, 3, 3], [2, 1, 2, 2]]
    sage: RSK(*gp) # RSK of a biword
    [[[1, 2, 2], [2]], [[1, 3, 3], [2]]]
    sage: RSK([2,3,2,1,2,3]) # Robinson-Schensted of a word
    [[[1, 2, 2, 3], [2], [3]], [[1, 2, 5, 6], [3], [4]]]
    sage: RSK([2,3,2,1,2,3], insertion=RSK.rules.EG) # Edelman-Greene
    [[[1, 2, 3], [2, 3], [3]], [[1, 2, 6], [3, 5], [4]]]
    sage: m = RSK_inverse(p, q, 'matrix'); m # output as matrix
    [0 1]
    [1 0]
    [0 2]
    sage: RSK(m) # RSK of a matrix
    [[[1, 2, 2], [2]], [[1, 3, 3], [2]]]

Insertions currently available
------------------------------

The following insertion algorithms for RSK correspondence are currently
available:

- RSK insertion (:class:`~sage.combinat.rsk.RuleRSK`)
- Edelman-Greene insertion (:class:`~sage.combinat.rsk.RuleEG`), an algorithm
  defined in [EG1987]_ Definition 6.20 (where it is referred to as
  Coxeter-Knuth insertion).
- Hecke RSK algorithm (:class:`~sage.combinat.rsk.RuleHecke`) , defined
  using the Hecke insertion studied in [BKSTY06]_ (but using rows instead
  of columns).
- Dual RSK insertion (:class:`~sage.combinat.rsk.RuleDualRSK`).
- CoRSK insertion (:class:`~sage.combinat.rsk.RuleCoRSK`), defined in [GR2018v5sol]_.

Implementing your own insertion rule
------------------------------------

The functions :func:`RSK()` and :func:`RSK_inverse()` are written so that it
is easy to implement insertion algorithms you come across in your research.

To implement your own insertion algorithm, you first need to import the
base class for a rule::

    sage: from sage.combinat.rsk import Rule

Using the ``Rule`` class as parent class for your insertion rule,
first implement the insertion and the reverse insertion algorithm
for :func:`RSK()` and :func:`RSK_inverse()` respectively (as methods
``forward_rule`` and ``backward_rule``). If your insertion algorithm
uses the same forward and backward rules as ``RuleRSK``, differing
only in how an entry is inserted into a row, then this is not
necessary, and it suffices to merely implement the
``insertion`` and ``reverse_insertion`` methods.

For more information, see :class:`~sage.combinat.rsk.Rule`.

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
   https://doi.org/10.1016/0001-8708(87)90063-6

.. [BKSTY06] \A. Buch, A. Kresch, M. Shimozono, H. Tamvakis, and A. Yong.
   *Stable Grothendieck polynomials and* `K`-*theoretic factor sequences*.
   Math. Ann. **340** Issue 2, (2008), pp. 359--382.
   :arxiv:`math/0601514v1`.

.. [GR2018v5sol] Darij Grinberg, Victor Reiner.
   *Hopf Algebras In Combinatorics*,
   :arXiv:`1409.8356v5`, available with solutions at
   https://arxiv.org/src/1409.8356v5/anc/HopfComb-v73-with-solutions.pdf
"""

# *****************************************************************************
#       Copyright (C) 2012,2019 Travis Scrimshaw <tcscrims gmail.com>
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
# *****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation

from bisect import bisect_left, bisect_right
from sage.structure.element import is_Matrix
from sage.matrix.all import matrix


class Rule(UniqueRepresentation):
    r"""
    Generic base class for an insertion rule for RSK correspondence.

    An instance of this class should implement a method
    :meth:`insertion` (which can be applied to a letter ``j``
    and a list ``r``, and modifies ``r`` in place by "bumping"
    ``j`` into it appropriately; it then returns the bumped-out
    entry or ``None`` if no such entry exists) and a method
    :meth:`reverse_insertion` (which does the same but for reverse
    bumping).
    It may also implement :meth:`_backward_format_output` and
    :meth:`_forward_format_output` if the RSK correspondence should
    return something other than (semi)standard tableaux (in the
    forward direction) and matrices or biwords (in the backward
    direction).
    The :meth:`to_pairs` method should also be overridden if
    the input for the (forward) RSK correspondence is not the
    usual kind of biwords (i.e., pairs of two `n`-tuples
    `[a_1, a_2, \ldots, a_n]` and `[b_1, b_2, \ldots, b_n]`
    satisfying `(a_1, b_1) \leq (a_2, b_2) \leq \cdots
    \leq (a_n, b_n)` in lexicographic order).
    Finally, it :meth:`forward_rule` and :meth:`backward_rule`
    have to be overridden if the overall structure of the
    RSK correspondence differs from that of classical RSK (see,
    e.g., the case of Hecke insertion, in which a letter bumped
    into a row may change a different row).
    """
    def to_pairs(self, obj1=None, obj2=None, check=True):
        r"""
        Given a valid input for the RSK algorithm, such as
        two `n`-tuples ``obj1`` `= [a_1, a_2, \ldots, a_n]`
        and ``obj2`` `= [b_1, b_2, \ldots, b_n]` forming a biword
        (i.e., satisfying
        `a_1 \leq a_2 \leq \cdots \leq a_n`, and if
        `a_i = a_{i+1}`, then `b_i \leq b_{i+1}`),
        or a matrix ("generalized permutation"), or a single word,
        return the array
        `[(a_1, b_1), (a_2, b_2), \ldots, (a_n, b_n)]`.

        INPUT:

        - ``obj1, obj2`` -- anything representing a biword
          (see the doc of :meth:`forward_rule` for the
          encodings accepted).

        - ``check`` -- (default: ``True``) whether to check
          that ``obj1`` and ``obj2`` actually define a valid
          biword.

        EXAMPLES::

            sage: from sage.combinat.rsk import Rule
            sage: list(Rule().to_pairs([1, 2, 2, 2], [2, 1, 1, 2]))
            [(1, 2), (2, 1), (2, 1), (2, 2)]
            sage: m = Matrix(ZZ, 3, 2, [0,1,1,0,0,2]) ; m
            [0 1]
            [1 0]
            [0 2]
            sage: list(Rule().to_pairs(m))
            [(1, 2), (2, 1), (3, 2), (3, 2)]
        """
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
            if check:
                if len(obj1) != len(obj2):
                    raise ValueError("the two arrays must be the same length")
                # Check it is a generalized permutation (i.e., biword):
                # that is, the pairs of corresponding entries of obj1
                # and obj2 are weakly increasing in lexicographic order.
                lt = 0
                lb = 0
                for t, b in zip(obj1, obj2):
                    if t < lt or (t == lt and b < lb):
                        raise ValueError("invalid generalized permutation")
                    lt = t
                    lb = b
            itr = zip(obj1, obj2)
        return itr

    def forward_rule(self, obj1, obj2, check_standard=False, check=True):
        r"""
        Return a pair of tableaux obtained by applying forward
        insertion to the generalized permutation ``[obj1, obj2]``.

        INPUT:

        - ``obj1, obj2`` -- can be one of the following ways to
          represent a generalized permutation (or, equivalently,
          biword):

          - two lists ``obj1`` and ``obj2`` of equal length,
            to be interpreted as the top row and the bottom row of
            the biword;

          - a matrix ``obj1`` of nonnegative integers, to be
            interpreted as the generalized permutation in matrix
            form (in this case, ``obj2`` is ``None``);

          - a word ``obj1`` in an ordered alphabet, to be
            interpreted as the bottom row of the biword (in this
            case, ``obj2`` is ``None``; the top row of the biword
            is understood to be `(1, 2, \ldots, n)` by default);

          - any object ``obj1`` which has a method ``_rsk_iter()``,
            as long as this method returns an iterator yielding
            pairs of numbers, which then are interperted as top
            entries and bottom entries in the biword (in this case,
            ``obj2`` is ``None``).

        - ``check_standard`` -- (default: ``False``) check if either of the
          resulting tableaux is a standard tableau, and if so, typecast it
          as such

        - ``check`` -- (default: ``True``) whether to check
          that ``obj1`` and ``obj2`` actually define a valid
          biword.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleRSK
            sage: RuleRSK().forward_rule([3,3,2,4,1], None)
            [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]
            sage: RuleRSK().forward_rule([1, 1, 1, 3, 7], None)
            [[[1, 1, 1, 3, 7]], [[1, 2, 3, 4, 5]]]
            sage: RuleRSK().forward_rule([7, 6, 3, 3, 1], None)
            [[[1, 3], [3], [6], [7]], [[1, 4], [2], [3], [5]]]
        """
        itr = self.to_pairs(obj1, obj2, check=check)
        p = []       # the "insertion" tableau
        q = []       # the "recording" tableau
        for i, j in itr:
            for r, qr in zip(p, q):
                j1 = self.insertion(j, r)
                if j1 is None:
                    r.append(j)
                    qr.append(i)  # Values are always inserted to the right
                    break
                else:
                    j = j1
            else:
                # We made through all of the rows of p without breaking
                # so we need to add a new row to p and q.
                p.append([j])
                q.append([i])
        return self._forward_format_output(p, q, check_standard=check_standard)

    def backward_rule(self, p, q, output):
        r"""
        Return the generalized permutation obtained by applying reverse
        insertion to a pair of tableaux ``(p, q)``.

        INPUT:

        - ``p``, ``q`` -- two tableaux of the same shape.

        - ``output`` -- (default: ``'array'``) if ``q`` is semi-standard:

          - ``'array'`` -- as a two-line array (i.e. generalized permutation
            or biword)
          - ``'matrix'`` -- as an integer matrix

          and if ``q`` is standard, we can also have the output:

          - ``'word'`` -- as a word

          and additionally if ``p`` is standard, we can also have the output:

          - ``'permutation'`` -- as a permutation

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleRSK
            sage: t1 = Tableau([[1, 3, 4], [2], [3]])
            sage: t2 = Tableau([[1, 2, 4], [3], [5]])
            sage: RuleRSK().backward_rule(t1, t2, 'array')
            [[1, 2, 3, 4, 5], [3, 3, 2, 4, 1]]
            sage: t1 = Tableau([[1, 1, 1, 3, 7]])
            sage: t2 = Tableau([[1, 2, 3, 4, 5]])
            sage: RuleRSK().backward_rule(t1, t2, 'array')
            [[1, 2, 3, 4, 5], [1, 1, 1, 3, 7]]
            sage: t1 = Tableau([[1, 3], [3], [6], [7]])
            sage: t2 = Tableau([[1, 4], [2], [3], [5]])
            sage: RuleRSK().backward_rule(t1, t2, 'array')
            [[1, 2, 3, 4, 5], [7, 6, 3, 3, 1]]
        """
        from sage.combinat.tableau import SemistandardTableaux
        # Make a copy of p since this is destructive to it
        p_copy = [list(row) for row in p]

        if q.is_standard():
            rev_word = []  # This will be our word in reverse
            d = {qij: i for i, Li in enumerate(q) for qij in Li}
            # d is now a dictionary which assigns to each integer k the
            # number of the row of q containing k.

            for key in sorted(d, reverse=True):
                # Delete last entry from i-th row of p_copy
                i = d[key]
                x = p_copy[i].pop()  # Always the right-most entry
                for row in reversed(p_copy[:i]):
                    x = self.reverse_insertion(x, row)
                rev_word.append(x)
            return self._backward_format_output(rev_word, None, output, p.is_standard(), True)

        if q not in SemistandardTableaux():
            raise ValueError("q(=%s) must be a semistandard tableau" %q)

        # Thus, q is semistandard but not standard.

        upper_row = []
        lower_row = []
        # upper_row and lower_row will be the upper and lower rows of the
        # generalized permutation we get as a result, but both reversed.
        d = {}
        for row, Li in enumerate(q):
            for col, val in enumerate(Li):
                if val in d:
                    d[val][col] = row
                else:
                    d[val] = {col: row}
        # d is now a double family such that for every integers k and j,
        # the value d[k][j] is the row i such that the (i, j)-th cell of
        # q is filled with k.
        for value, row_dict in sorted(d.items(), reverse=True, key=lambda x: x[0]):
            for key in sorted(row_dict, reverse=True):
                i = row_dict[key]
                x = p_copy[i].pop()  # Always the right-most entry
                for row in reversed(p_copy[:i]):
                    x = self.reverse_insertion(x, row)
                lower_row.append(x)
                upper_row.append(value)
        return self._backward_format_output(lower_row, upper_row, output, p.is_standard(), False)

    def _forward_format_output(self, p, q, check_standard):
        r"""
        Return final output of the ``RSK`` correspondence from the
        output of the corresponding ``forward_rule``.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleRSK
            sage: isinstance(RuleRSK()._forward_format_output([[1, 2, 3, 4, 5]]
            ....:            , [[1, 2, 3, 4, 5]], True)[0], StandardTableau)
            True
            sage: isinstance(RuleRSK()._forward_format_output([[1, 2, 3, 4, 5]]
            ....:          , [[1, 2, 3, 4, 5]], False)[0], SemistandardTableau)
            True
            sage: isinstance(RuleRSK()._forward_format_output([[1, 1, 1, 3, 7]]
            ....:          , [[1, 2, 3, 4, 5]], True)[0], SemistandardTableau)
            True
        """
        from sage.combinat.tableau import SemistandardTableau, StandardTableau

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

    def _backward_format_output(self, lower_row, upper_row, output,
                                p_is_standard, q_is_standard):
        r"""
        Return the final output of the ``RSK_inverse`` correspondence
        from the output of the corresponding ``backward_rule``.

        .. NOTE::

            The default implementation of ``backward_rule`` lists
            bumped-out entries in the order in which the reverse
            bumping happens, which is *opposite* to the order of the
            final output.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleRSK
            sage: RuleRSK()._backward_format_output([1, 2, 3, 4], None, 'array', False, True)
            [[1, 2, 3, 4], [4, 3, 2, 1]]
            sage: RuleRSK()._backward_format_output([1, 2, 3, 4], None,
            ....:                                   'matrix', False, True)
            [0 0 0 1]
            [0 0 1 0]
            [0 1 0 0]
            [1 0 0 0]
            sage: RuleRSK()._backward_format_output([1, 2, 3, 4], None, 'word', False, True)
            word: 4321
            sage: RuleRSK()._backward_format_output([3, 2, 1, 1], [2, 1, 1, 1], 'array', False, False)
            [[1, 1, 1, 2], [1, 1, 2, 3]]
            sage: RuleRSK()._backward_format_output([3, 2, 1, 1], [2, 1, 1, 1],
            ....:                                   'matrix', False, False)
            [2 1 0]
            [0 0 1]
            sage: RuleRSK()._backward_format_output([1, 2, 3, 4], None, 'word', False, False)
            Traceback (most recent call last):
            ...
            TypeError: q must be standard to have a word as valid output
            sage: RuleRSK()._backward_format_output([1, 2, 3, 4], None,
            ....:                                     'random_type', False, False)
            Traceback (most recent call last):
            ...
            ValueError: invalid output option
        """
        if q_is_standard:
            if output == 'word':
                from sage.combinat.words.word import Word
                return Word(reversed(lower_row))
            if output == 'matrix':
                return to_matrix(list(range(1, len(lower_row)+1)), list(reversed(lower_row)))
            if output == 'array':
                return [list(range(1, len(lower_row)+1)), list(reversed(lower_row))]
            raise ValueError("invalid output option")

        else:
            if output == 'matrix':
                return to_matrix(list(reversed(upper_row)), list(reversed(lower_row)))
            if output == 'array':
                return [list(reversed(upper_row)), list(reversed(lower_row))]
            if output in ['permutation', 'word']:
                raise TypeError(
                    "q must be standard to have a %s as valid output" %output)
            raise ValueError("invalid output option")

class RuleRSK(Rule):
    r"""
    A rule modeling the classical Robinson-Schensted-Knuth insertion.

    See :func:`RSK` for the definition of this operation.

    EXAMPLES::

        sage: RSK([1, 2, 2, 2], [2, 1, 1, 2], insertion=RSK.rules.RSK)
        [[[1, 1, 2], [2]], [[1, 2, 2], [2]]]
        sage: p = Tableau([[1,2,2],[2]]); q = Tableau([[1,3,3],[2]])
        sage: RSK_inverse(p, q, insertion=RSK.rules.RSK)
        [[1, 2, 3, 3], [2, 1, 2, 2]]

    It is worth noting that in case of ``RSK_inverse`` with
    ``output = 'permutation'`` ``RuleRSK`` returns an object of class
    :class:`~sage.combinat.permutation.Permutation`.
    """
    def insertion(self, j, r):
        r"""
        Insert the letter ``j`` from the second row of the biword
        into the row `r` using classical Schensted insertion,
        if there is bumping to be done.

        The row `r` is modified in place. The bumped-out entry,
        if it exists, is returned.

        .. WARNING::

            This method only changes `r` if bumping occurs.
            Appending `j` to the end of the row should be done
            by the caller.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleRSK
            sage: qr, r = [1,2,3,4,5], [3,3,2,4,8]
            sage: j = RuleRSK().insertion(9, r)
            sage: j is None
            True
            sage: qr, r = [1,2,3,4,5], [3,3,2,4,8]
            sage: j = RuleRSK().insertion(3, r)
            sage: j
            4
        """
        if r[-1] <= j:
            return None  # j needs to be added at the end of the row.
        # Figure out where to insert j into the row r. The
        # bisect command returns the position of the least
        # element of r greater than j.  We will call it y.
        y_pos = bisect_right(r, j)
        # Switch j and y
        j, r[y_pos] = r[y_pos], j
        return j

    def reverse_insertion(self, x, row):
        r"""
        Reverse bump the row ``row`` of the current insertion tableau
        with the number ``x``.

        The row ``row`` is modified in place. The bumped-out entry
        is returned.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleRSK
            sage: r =  [2,3,3,4,8]
            sage: x = RuleRSK().reverse_insertion(4, r); r
            [2, 3, 4, 4, 8]
            sage: x
            3
        """
        y_pos = bisect_left(row, x) - 1
        # switch x and y
        x, row[y_pos] = row[y_pos], x
        return x

    def _backward_format_output(self, lower_row, upper_row, output,
                                p_is_standard, q_is_standard):
        r"""
        Return the final output of the ``RSK_inverse`` correspondence
        from the output of the corresponding ``backward_rule``.

        .. NOTE::

            The default implementation of ``backward_rule`` lists
            bumped-out entries in the order in which the reverse
            bumping happens, which is *opposite* to the order of the
            final output.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleRSK
            sage: RuleRSK()._backward_format_output([1, 2, 3, 4], None,
            ....:                                   'permutation', True, True)
            [4, 3, 2, 1]
            sage: RuleRSK()._backward_format_output([1, 2, 3, 4], None,
            ....:                                   'permutation', False, True)
            Traceback (most recent call last):
            ...
            TypeError: p must be standard to have a valid permutation as output
        """
        if q_is_standard and output == 'permutation':
            if not p_is_standard:
                raise TypeError("p must be standard to have a valid permutation as output")
            from sage.combinat.permutation import Permutation
            return Permutation(reversed(lower_row))
        else:
            return super(RuleRSK, self)._backward_format_output(lower_row, upper_row, output,
                                                                p_is_standard, q_is_standard)


class RuleEG(Rule):
    r"""
    A rule modeling Edelman-Greene insertion.

    For a reduced word of a permutation (i.e., an element of a type `A`
    Coxeter group), one can use Edelman-Greene insertion, an algorithm
    defined in [EG1987]_ Definition 6.20 (where it is referred to as
    Coxeter-Knuth insertion). The Edelman-Greene insertion is similar to the
    standard row insertion except that (using the notations in
    the documentation of :func:`RSK`) if `k_i` and `k_i + 1` both
    exist in row `i`, we *only* set `k_{i+1} = k_i + 1` and continue.

    EXAMPLES:

    Let us reproduce figure 6.4 in [EG1987]_::

        sage: RSK([2,3,2,1,2,3], insertion=RSK.rules.EG)
        [[[1, 2, 3], [2, 3], [3]], [[1, 2, 6], [3, 5], [4]]]

    Some more examples::

        sage: a = [2, 1, 2, 3, 2]
        sage: pq = RSK(a, insertion=RSK.rules.EG); pq
        [[[1, 2, 3], [2, 3]], [[1, 3, 4], [2, 5]]]
        sage: RSK(RSK_inverse(*pq, insertion=RSK.rules.EG, output='matrix'),
        ....:     insertion=RSK.rules.EG)
        [[[1, 2, 3], [2, 3]], [[1, 3, 4], [2, 5]]]
        sage: RSK_inverse(*pq, insertion=RSK.rules.EG)
        [[1, 2, 3, 4, 5], [2, 1, 2, 3, 2]]

    The RSK algorithm (:func:`RSK`) built using the Edelman-Greene
    insertion rule ``RuleEG`` is a bijection from reduced words of
    permutations/elements of a type `A` Coxeter group to pairs
    consisting of an increasing tableau and a standard tableau
    of the same shape (see [EG1987]_ Theorem 6.25).
    The inverse of this bijection is obtained using :func:`RSK_inverse`.
    If the optional parameter ``output = 'permutation'`` is set in
    :func:`RSK_inverse`, then the function returns not the
    reduced word itself but the permutation (of smallest possible
    size) whose reduced word it is (although the order of the
    letters is reverse to the usual Sage convention)::

        sage: w = RSK_inverse(*pq, insertion=RSK.rules.EG, output='permutation'); w
        [4, 3, 1, 2]
        sage: list(reversed(a)) in w.reduced_words()
        True

    TESTS:

    Let us check that :func:`RSK_inverse` is the inverse of :func:`RSK`
    on the different types of inputs/outputs for Edelman-Greene.
    First we can check on the reduced words (specifically, those that can
    be obtained using the ``reduced_word()`` method from permutations)::

        sage: g = lambda w: RSK_inverse(*RSK(w, insertion=RSK.rules.EG),
        ....:                           insertion=RSK.rules.EG, output='word')
        sage: all(p.reduced_word() == list(g(p.reduced_word()))
        ....:     for n in range(7) for p in Permutations(n))
        True

    In case of non-standard tableaux `P, Q`::

        sage: RSK_inverse(*RSK([1, 2, 3, 2, 1], insertion='EG'),
        ....:             insertion='EG')
        [[1, 2, 3, 4, 5], [1, 2, 3, 2, 1]]
        sage: RSK_inverse(*RSK([1, 1, 1, 2], [1, 2, 3, 4],
        ....:             insertion=RSK.rules.EG), insertion=RSK.rules.EG)
        [[1, 1, 1, 2], [1, 2, 3, 4]]
        sage: RSK_inverse(*RSK([1, 2, 3, 3], [2, 1, 2, 2], insertion='EG'),
        ....:             insertion='EG')
        [[1, 2, 3, 3], [2, 1, 2, 2]]

    Since the column reading of the insertion tableau from
    Edelman-Greene insertion gives one of reduced words for the
    original permutation, we can also check for that::

        sage: f = lambda p: [x for row in reversed(p) for x in row]
        sage: g = lambda wd: RSK(wd, insertion=RSK.rules.EG)[0]
        sage: all(p == Permutations(n).from_reduced_word(f(g(wd)))
        ....:     for n in range(6) for p in Permutations(n)
        ....:     for wd in p.reduced_words())
        True
    """
    def insertion(self, j, r):
        r"""
        Insert the letter ``j`` from the second row of the biword
        into the row `r` using Edelman-Greene insertion,
        if there is bumping to be done.

        The row `r` is modified in place. The bumped-out entry,
        if it exists, is returned.

        .. WARNING::

            This method only changes `r` if bumping occurs.
            Appending `j` to the end of the row should be done
            by the caller.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleEG
            sage: qr, r =  [1,2,3,4,5], [3,3,2,4,8]
            sage: j = RuleEG().insertion(9, r)
            sage: j is None
            True
            sage: qr, r = [1,2,3,4,5], [2,3,4,5,8]
            sage: j = RuleEG().insertion(3, r); r
            [2, 3, 4, 5, 8]
            sage: j
            4
            sage: qr, r = [1,2,3,4,5], [2,3,5,5,8]
            sage: j = RuleEG().insertion(3, r); r
            [2, 3, 3, 5, 8]
            sage: j
            5
        """
        if r[-1] <= j:
            return None  # j needs to be added at the end of the row
        # Figure out where to insert j into the row r. The
        # bisect command returns the position of the least
        # element of r greater than j.  We will call it y.
        y_pos = bisect_right(r, j)
        if r[y_pos] == j + 1 and y_pos > 0 and j == r[y_pos - 1]:
            # Special bump: Nothing to do except increment j by 1
            j += 1
        else:
            # Switch j and y
            j, r[y_pos] = r[y_pos], j
        return j

    def reverse_insertion(self, x, row):
        r"""
        Reverse bump the row ``row`` of the current insertion tableau
        with the number ``x``.

        The row ``row`` is modified in place. The bumped-out entry
        is returned.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleEG
            sage: r =  [1,1,1,2,3,3]
            sage: x = RuleEG().reverse_insertion(3, r); r
            [1, 1, 1, 2, 3, 3]
            sage: x
            2
        """
        y_pos = bisect_left(row, x) - 1
        if row[y_pos] == x - 1 and y_pos < len(row)-1 and row[y_pos+1] == x:
            # Nothing to do except decrement x by 1.
            # (Case 1 on p. 74 of Edelman-Greene [EG1987]_.)
            x -= 1
        else:
            # switch x and y
            x, row[y_pos] = row[y_pos], x
        return x

    def _backward_format_output(self, lower_row, upper_row, output,
                                p_is_standard, q_is_standard):
        r"""
        Return the final output of the ``RSK_inverse`` correspondence
        from the output of the corresponding ``backward_rule``.

        .. NOTE::

            The default implementation of ``backward_rule`` lists
            bumped-out entries in the order in which the reverse
            bumping happens, which is *opposite* to the order of the
            final output.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleEG
            sage: RuleEG()._backward_format_output([1, 2, 3, 4], None,
            ....:                                  'permutation', True, False)
            Traceback (most recent call last):
            ...
            TypeError: q must be standard to have a permutation as valid output
        """
        if q_is_standard and output == 'permutation':
            n = 0
            if list(lower_row):
                n = max(list(lower_row)) + 1
            from sage.combinat.permutation import Permutations
            return Permutations(n).from_reduced_word(list(lower_row))
        else:
            return super(RuleEG, self)._backward_format_output(lower_row, upper_row, output,
                                                               p_is_standard, q_is_standard)


class RuleHecke(Rule):
    r"""
    A rule modeling the Hecke insertion algorithm.

    The Hecke RSK algorithm is similar to the classical RSK algorithm,
    but is defined using the Hecke insertion introduced in in
    [BKSTY06]_ (but using rows instead of columns).
    It is not clear in what generality it works; thus, following
    [BKSTY06]_, we shall assume that our biword `p` has top row
    `(1, 2, \ldots, n)` (or, at least, has its top row strictly
    increasing).

    The Hecke RSK algorithm returns a pair of an increasing tableau
    and a set-valued standard tableau. If
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
    entry is a set, since `Q` is set-valued). In either
    subcase, we terminate the recursion, and set
    `P_{t+1} = P` and `Q_{t+1} = Q`.

    Notice that set-valued tableaux are encoded as tableaux whose
    entries are tuples of positive integers; each such tuple is strictly
    increasing and encodes a set (namely, the set of its entries).

    EXAMPLES:

    As an example of Hecke insertion, we reproduce
    Example 2.1 in :arxiv:`0801.1319v2`::

        sage: w = [5, 4, 1, 3, 4, 2, 5, 1, 2, 1, 4, 2, 4]
        sage: P,Q = RSK(w, insertion=RSK.rules.Hecke); [P,Q]
        [[[1, 2, 4, 5], [2, 4, 5], [3, 5], [4], [5]],
         [[(1,), (4,), (5,), (7,)],
          [(2,), (9,), (11, 13)],
          [(3,), (12,)],
          [(6,)],
          [(8, 10)]]]
        sage: wp = RSK_inverse(P, Q, insertion=RSK.rules.Hecke,
        ....:                    output='list'); wp
        [5, 4, 1, 3, 4, 2, 5, 1, 2, 1, 4, 2, 4]
        sage: wp == w
        True
    """
    def forward_rule(self, obj1, obj2, check_standard=False):
        r"""
        Return a pair of tableaux obtained by applying Hecke
        insertion to the generalized permutation ``[obj1, obj2]``.

        INPUT:

        - ``obj1, obj2`` -- can be one of the following ways to
          represent a generalized permutation (or, equivalently,
          biword):

          - two lists ``obj1`` and ``obj2`` of equal length,
            to be interpreted as the top row and the bottom row of
            the biword;

          - a word ``obj1`` in an ordered alphabet, to be
            interpreted as the bottom row of the biword (in this
            case, ``obj2`` is ``None``; the top row of the biword
            is understood to be `(1, 2, \ldots, n)` by default);

        - ``check_standard`` -- (default: ``False``) check if either of the
          resulting tableaux is a standard tableau, and if so, typecast it
          as such

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleHecke
            sage: p, q = RuleHecke().forward_rule([3,3,2,4,1], None);p
            [[1, 4], [2], [3]]
            sage: q
            [[(1, 2), (4,)], [(3,)], [(5,)]]
            sage: isinstance(p, SemistandardTableau)
            True
            sage: isinstance(q, Tableau)
            True
        """
        from sage.combinat.tableau import SemistandardTableau, Tableau

        if obj2 is None:
            obj2 = obj1
            obj1 = list(range(1, len(obj1) + 1))

        p = []       # the "insertion" tableau
        q = []       # the "recording" tableau

        for i, j in zip(obj1, obj2):
            for ir, r in enumerate(p):
                j1 = self.insertion(j, ir, r, p)

                if j1 is None:
                    # We must have len(p[ir-1]) > len(r), since j is coming
                    # from the previous row.
                    if r[-1] < j and (ir == 0 or p[ir-1][len(r)] < j):
                        # We can add a box to the row
                        r.append(j)
                        q[ir].append((i,))  # Values are always inserted to the right
                    else:
                        # We must append i to the bottom of this column
                        l = len(r) - 1
                        while ir < len(q) and len(q[ir]) > l:
                            ir += 1
                        q[ir-1][-1] = q[ir-1][-1] + (i,)
                    break
                else:
                    j = j1
            else:
                # We made through all of the rows of p without breaking
                # so we need to add a new row to p and q.
                p.append([j])
                q.append([(i,)])
        return [SemistandardTableau(p), Tableau(q)]

    def backward_rule(self, p, q, output):
        r"""
        Return the generalized permutation obtained by applying reverse
        Hecke insertion to a pair of tableaux ``(p, q)``.

        INPUT:

        - ``p``, ``q`` -- two tableaux of the same shape.

        -  ``output`` -- (default: ``'array'``) if ``q`` is semi-standard:

          - ``'array'`` -- as a two-line array (i.e. generalized permutation
            or biword)

          and if ``q`` is standard set-valued, we can have the output:

          - ``'word'`` -- as a word
          - ``'list'`` -- as a list

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleHecke
            sage: t1 = Tableau([[1, 4], [2], [3]])
            sage: t2 = Tableau([[(1, 2), (4,)], [(3,)], [(5,)]])
            sage: RuleHecke().backward_rule(t1, t2, 'array')
            [[1, 2, 3, 4, 5], [3, 3, 2, 4, 1]]
            sage: t1 = Tableau([[1, 4], [2, 3]])
            sage: t2 = Tableau([[(1, 2), (4,)], [(3,)], [(5,)]])
            sage: RuleHecke().backward_rule(t1, t2, 'array')
            Traceback (most recent call last):
            ...
            ValueError: p(=[[1, 4], [2, 3]]) and
             q(=[[(1, 2), (4,)], [(3,)], [(5,)]]) must have the same shape
        """
        if p.shape() != q.shape():
            raise ValueError("p(=%s) and q(=%s) must have the same shape" % (p, q))
        from sage.combinat.tableau import SemistandardTableaux
        if p not in SemistandardTableaux():
            raise ValueError("p(=%s) must be a semistandard tableau" % p)

        # Make a copy of p and q since this is destructive to it
        p_copy = [list(row) for row in p]
        q_copy = [[list(v) for v in row] for row in q]
        # We shall work on these copies of p and q. Notice that p might get
        # some empty rows in the process; we do not bother pruning them, as
        # they do not matter.

        # upper_row and lower_row will be the upper and lower rows of the
        # generalized permutation we get as a result, but both reversed.
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
        # d is now a double family such that for every integers k and j,
        # the value d[k][j] is the row i such that the (i, j)-th cell of
        # q is filled with k.
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
                    x = self.reverse_insertion(i, x, row, p_copy)
                lower_row.append(x)
                upper_row.append(value)
        return self._backward_format_output(lower_row, upper_row, output, False, False)

    def insertion(self, j, ir, r, p):
        r"""
        Insert the letter ``j`` from the second row of the biword
        into the row `r` of the increasing tableau `p` using
        Hecke insertion, provided that `r` is the `ir`-th row
        of `p`, and provided that there is bumping to be done.

        The row `r` is modified in place. The bumped-out entry,
        if it exists, is returned.

        .. WARNING::

            This method only changes `r` if bumping occurs.
            Appending `j` to the end of the row should be done
            by the caller.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleHecke
            sage: from bisect import bisect_right
            sage: p, q, r =  [], [], [3,3,8,8,8,9]
            sage: j, ir = 8, 1
            sage: j1 = RuleHecke().insertion(j, ir, r, p)
            sage: j1 == r[bisect_right(r, j)]
            True
        """
        if r[-1] <= j:
            return None  # j needs to be added at the end of the row.
        # Figure out where to insert j into the row r.  The
        # bisect command returns the position of the least
        # element of r greater than j.  We will call it y.
        y_pos = bisect_right(r, j)
        y = r[y_pos]
        # Check to see if we can swap j for y
        if (y_pos == 0 or r[y_pos-1] < j) and (ir == 0 or p[ir-1][y_pos] < j):
            r[y_pos] = j
        j = y
        return j

    def reverse_insertion(self, i, x, row, p):
        r"""
        Reverse bump the row ``row`` of the current insertion tableau
        ``p`` with the number ``x``, provided that ``row`` is the
        `i`-th row of `p`.

        The row ``row`` is modified in place. The bumped-out entry
        is returned.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleHecke
            sage: from bisect import bisect_left
            sage: r =  [2,3,3,4,8,9]
            sage: x, i, p = 9, 1, [1, 2]
            sage: x1 = RuleHecke().reverse_insertion(i, x, r, p)
            sage: x1 == r[bisect_left(r,x) - 1]
            True
        """
        y_pos = bisect_left(row, x) - 1
        y = row[y_pos]
        # Check to see if we can swap x for y
        if ((y_pos == len(row) - 1 or x < row[y_pos+1])
            and (i == len(p) - 1 or len(p[i+1]) <= y_pos
                 or x < p[i+1][y_pos])):
            row[y_pos] = x
        x = y
        return x

    def _backward_format_output(self, lower_row, upper_row, output,
                                p_is_standard, q_is_standard):
        r"""
        Return the final output of the ``RSK_inverse`` correspondence
        from the output of the corresponding ``backward_rule``.

        .. NOTE::

            The default implementation of ``backward_rule`` lists
            bumped-out entries in the order in which the reverse
            bumping happens, which is *opposite* to the order of the
            final output.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleHecke
            sage: RuleHecke()._backward_format_output([1, 1, 3, 9], [1, 2, 3, 4],
            ....:                                     'array', False, False)
            [[4, 3, 2, 1], [9, 3, 1, 1]]
            sage: RuleHecke()._backward_format_output([1, 1, 3, 9], [1, 2, 3, 4],
            ....:                                     'word', False, False)
            Traceback (most recent call last):
            ...
            TypeError: q must be standard to have a word as valid output
            sage: RuleHecke()._backward_format_output([1, 1, 3, 9], [4, 3, 2, 1],
            ....:                                     'word', False, False)
            word: 9311
            sage: RuleHecke()._backward_format_output([1, 1, 3, 9], [1, 2, 3, 4],
            ....:                                     'list', False, False)
            Traceback (most recent call last):
            ...
            TypeError: q must be standard to have a list as valid output
            sage: RuleHecke()._backward_format_output([1, 1, 3, 9], [4, 3, 2, 1],
            ....:                                     'list', False, False)
            [9, 3, 1, 1]
            sage: RuleHecke()._backward_format_output([1, 1, 3, 9], [1, 2, 3, 4],
            ....:                                     'random_type', False, False)
            Traceback (most recent call last):
            ...
            ValueError: invalid output option
        """
        if output == 'array':
            return [list(reversed(upper_row)), list(reversed(lower_row))]
        is_standard = (upper_row == list(range(len(upper_row), 0, -1)))
        if output == 'word':
            if not is_standard:
                raise TypeError(
                    "q must be standard to have a %s as valid output" % output)
            from sage.combinat.words.word import Word
            return Word(reversed(lower_row))
        if output == 'list':
            if not is_standard:
                raise TypeError(
                    "q must be standard to have a %s as valid output" % output)
            return list(reversed(lower_row))
        raise ValueError("invalid output option")

class RuleDualRSK(Rule):
    r"""
    A rule modeling the Dual RSK insertion.

    Dual RSK insertion differs from classical RSK insertion in the
    following ways:

    * The input (in terms of biwords) is no longer an arbitrary biword,
      but rather a strict biword (i.e., a pair of two lists
      `[a_1, a_2, \ldots, a_n]` and `[b_1, b_2, \ldots, b_n]` that
      satisfy the strict inequalities
      `(a_1, b_1) < (a_2, b_2) < \cdots < (a_n, b_n)` in the
      lexicographic order).
      In terms of matrices, this means that the input is not an
      arbitrary matrix with nonnegative integer entries, but rather
      a `\{0, 1\}`-matrix (i.e., a matrix whose entries are `0`'s
      and `1`'s).

    * The output still consists of two tableaux `(P, Q)` of equal
      shapes, but rather than both of them being semistandard, now
      `P` is row-strict (i.e., its transpose is semistandard) while
      `Q` is semistandard.

    * The main difference is in the way bumping works. Namely,
      when a number `k_i` is inserted into the `i`-th row of `P`,
      it bumps out the first integer greater **or equal to** `k_i`
      in this row (rather than greater than `k_i`).

    The RSK and dual RSK algorithms agree for permutation matrices.

    For more information, see Chapter 7, Section 14 in [Sta-EC2]_
    (where dual RSK is called `\mathrm{RSK}^{\ast}`) or the third
    solution to Exercise 2.7.12(a) in [GR2018v5sol]_.

    EXAMPLES::

        sage: RSK([3,3,2,4,1], insertion=RSK.rules.dualRSK)
        [[[1, 4], [2], [3], [3]], [[1, 4], [2], [3], [5]]]
        sage: RSK(Word([3,3,2,4,1]), insertion=RSK.rules.dualRSK)
        [[[1, 4], [2], [3], [3]], [[1, 4], [2], [3], [5]]]
        sage: RSK(Word([2,3,3,2,1,3,2,3]), insertion=RSK.rules.dualRSK)
        [[[1, 2, 3], [2, 3], [2, 3], [3]], [[1, 2, 8], [3, 6], [4, 7], [5]]]

    Using dual RSK insertion with a strict biword::

        sage: RSK([1,1,2,4,4,5],[2,4,1,1,3,2], insertion=RSK.rules.dualRSK)
        [[[1, 2], [1, 3], [2, 4]], [[1, 1], [2, 4], [4, 5]]]
        sage: RSK([1,1,2,3,3,4,5],[1,3,2,1,3,3,2], insertion=RSK.rules.dualRSK)
        [[[1, 2, 3], [1, 2], [3], [3]], [[1, 1, 3], [2, 4], [3], [5]]]
        sage: RSK([1, 2, 2, 2], [2, 1, 2, 4], insertion=RSK.rules.dualRSK)
        [[[1, 2, 4], [2]], [[1, 2, 2], [2]]]
        sage: RSK(Word([1,1,3,4,4]), [1,4,2,1,3], insertion=RSK.rules.dualRSK)
        [[[1, 2, 3], [1], [4]], [[1, 1, 4], [3], [4]]]
        sage: RSK([1,3,3,4,4], Word([6,1,2,1,7]), insertion=RSK.rules.dualRSK)
        [[[1, 2, 7], [1], [6]], [[1, 3, 4], [3], [4]]]

    Using dual RSK insertion with a `\{0, 1\}`-matrix::

        sage: RSK(matrix([[0,1],[1,1]]), insertion=RSK.rules.dualRSK)
        [[[1, 2], [2]], [[1, 2], [2]]]

    We can also give it something looking like a matrix::

        sage: RSK([[0,1],[1,1]], insertion=RSK.rules.dualRSK)
        [[[1, 2], [2]], [[1, 2], [2]]]

    Let us now call the inverse correspondence::

        sage: RSK_inverse(*RSK([1, 2, 2, 2], [2, 1, 2, 3],
        ....:         insertion=RSK.rules.dualRSK),insertion=RSK.rules.dualRSK)
        [[1, 2, 2, 2], [2, 1, 2, 3]]
        sage: P,Q = RSK([1, 2, 2, 2], [2, 1, 2, 3],insertion=RSK.rules.dualRSK)
        sage: RSK_inverse(P, Q, insertion=RSK.rules.dualRSK)
        [[1, 2, 2, 2], [2, 1, 2, 3]]

    When applied to two standard tableaux, reverse dual RSK
    insertion behaves identically to the usual reverse RSK insertion::

        sage: t1 = Tableau([[1, 2, 5], [3], [4]])
        sage: t2 = Tableau([[1, 2, 3], [4], [5]])
        sage: RSK_inverse(t1, t2, insertion=RSK.rules.dualRSK)
        [[1, 2, 3, 4, 5], [1, 4, 5, 3, 2]]
        sage: RSK_inverse(t1, t2, 'word', insertion=RSK.rules.dualRSK)
        word: 14532
        sage: RSK_inverse(t1, t2, 'matrix', insertion=RSK.rules.dualRSK)
        [1 0 0 0 0]
        [0 0 0 1 0]
        [0 0 0 0 1]
        [0 0 1 0 0]
        [0 1 0 0 0]
        sage: RSK_inverse(t1, t2, 'permutation', insertion=RSK.rules.dualRSK)
        [1, 4, 5, 3, 2]
        sage: RSK_inverse(t1, t1, 'permutation', insertion=RSK.rules.dualRSK)
        [1, 4, 3, 2, 5]
        sage: RSK_inverse(t2, t2, 'permutation', insertion=RSK.rules.dualRSK)
        [1, 2, 5, 4, 3]
        sage: RSK_inverse(t2, t1, 'permutation', insertion=RSK.rules.dualRSK)
        [1, 5, 4, 2, 3]

    Let us check that forward and backward dual RSK are mutually
    inverse when the first tableau is merely transpose semistandard::

        sage: p = Tableau([[1,2,2],[1]]); q = Tableau([[1,2,4],[3]])
        sage: ret = RSK_inverse(p, q, insertion=RSK.rules.dualRSK); ret
        [[1, 2, 3, 4], [1, 2, 1, 2]]
        sage: RSK_inverse(p, q, 'word', insertion=RSK.rules.dualRSK)
        word: 1212

    In general for dual RSK::

        sage: p = Tableau([[1,1,2],[1]]); q = Tableau([[1,3,3],[2]])
        sage: RSK_inverse(p, q, insertion=RSK.rules.dualRSK)
        [[1, 2, 3, 3], [1, 1, 1, 2]]
        sage: RSK_inverse(p, q, 'matrix', insertion=RSK.rules.dualRSK)
        [1 0]
        [1 0]
        [1 1]

    TESTS:

    Empty objects::

        sage: RSK(Permutation([]), insertion=RSK.rules.dualRSK)
        [[], []]
        sage: RSK(Word([]), insertion=RSK.rules.dualRSK)
        [[], []]
        sage: RSK(matrix([[]]), insertion=RSK.rules.dualRSK)
        [[], []]
        sage: RSK([], [], insertion=RSK.rules.dualRSK)
        [[], []]
        sage: RSK([[]], insertion=RSK.rules.dualRSK)
        [[], []]

    Check that :func:`RSK_inverse` is the inverse of :func:`RSK` on the
    different types of inputs/outputs::

        sage: RSK_inverse(Tableau([]), Tableau([]),
        ....:                insertion=RSK.rules.dualRSK)
        [[], []]
        sage: f = lambda p: RSK_inverse(*RSK(p, insertion=RSK.rules.dualRSK),
        ....:                output='permutation', insertion=RSK.rules.dualRSK)
        sage: all(p == f(p) for n in range(7) for p in Permutations(n))
        True
        sage: all(RSK_inverse(*RSK(w, insertion=RSK.rules.dualRSK),
        ....:                 output='word', insertion=RSK.rules.dualRSK) == w
        ....:     for n in range(4) for w in Words(5, n))
        True
        sage: from sage.combinat.integer_matrices import IntegerMatrices
        sage: M = IntegerMatrices([1,2,2,1], [3,1,1,1]) #this is probably wrong
        sage: all(RSK_inverse(*RSK(m, insertion=RSK.rules.dualRSK),
        ....:                output='matrix', insertion=RSK.rules.dualRSK) == m
        ....:     for m in M if all(x in [0, 1] for x in m))
        True

        sage: n = ZZ.random_element(200)
        sage: p = Permutations(n).random_element()
        sage: True if p == f(p) else p
        True

    Checking that the tableaux should be of same shape::

        sage: RSK_inverse(Tableau([[1,2,3]]), Tableau([[1,2]]),
        ....:                          insertion=RSK.rules.dualRSK)
        Traceback (most recent call last):
        ...
        ValueError: p(=[[1, 2, 3]]) and q(=[[1, 2]]) must have the same shape
    """
    def to_pairs(self, obj1=None, obj2=None, check=True):
        r"""
        Given a valid input for the dual RSK algorithm, such as
        two `n`-tuples ``obj1`` `= [a_1, a_2, \ldots, a_n]`
        and ``obj2`` `= [b_1, b_2, \ldots, b_n]` forming a strict
        biword (i.e., satisfying
        `a_1 \leq a_2 \leq \cdots \leq a_n`, and if
        `a_i = a_{i+1}`, then `b_i < b_{i+1}`)
        or a `\{0, 1\}`-matrix ("rook placement"), or a
        single word, return the array
        `[(a_1, b_1), (a_2, b_2), \ldots, (a_n, b_n)]`.

        INPUT:

        - ``obj1, obj2`` -- anything representing a strict biword
          (see the doc of :meth:`forward_rule` for the
          encodings accepted).

        - ``check`` -- (default: ``True``) whether to check
          that ``obj1`` and ``obj2`` actually define a valid
          strict biword.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleDualRSK
            sage: list(RuleDualRSK().to_pairs([1, 2, 2, 2], [2, 1, 2, 3]))
            [(1, 2), (2, 1), (2, 2), (2, 3)]
            sage: RuleDualRSK().to_pairs([1, 2, 2, 2], [1, 2, 3, 3])
            Traceback (most recent call last):
            ...
            ValueError: invalid strict biword
            sage: m = Matrix(ZZ, 3, 2, [0,1,1,1,0,1]) ; m
            [0 1]
            [1 1]
            [0 1]
            sage: list(RuleDualRSK().to_pairs(m))
            [(1, 2), (2, 1), (2, 2), (3, 2)]
            sage: m = Matrix(ZZ, 3, 2, [0,1,1,0,0,2]) ; m
            [0 1]
            [1 0]
            [0 2]
            sage: RuleDualRSK().to_pairs(m)
            Traceback (most recent call last):
            ...
            ValueError: dual RSK requires a {0, 1}-matrix
        """
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
                            if mult > 1:
                                raise ValueError("dual RSK requires a {0, 1}-matrix")
                            if mult > 0:
                                t.append(i+1)
                                b.append(j+1)
                    itr = zip(t, b)
                except TypeError:
                    itr = zip(range(1, len(obj1)+1), obj1)
        else:
            if check:
                if len(obj1) != len(obj2):
                    raise ValueError("the two arrays must be the same length")
                # Check it is a generalized permutation (i.e., biword):
                # that is, the pairs of corresponding entries of obj1
                # and obj2 are weakly increasing in lexicographic order.
                lt = 0
                lb = 0
                for t, b in zip(obj1, obj2):
                    if t < lt or (t == lt and b <= lb):
                        raise ValueError("invalid strict biword")
                    lt = t
                    lb = b
            itr = zip(obj1, obj2)
        return itr

    def insertion(self, j, r):
        r"""
        Insert the letter ``j`` from the second row of the biword
        into the row `r` using dual RSK insertion, if there is
        bumping to be done.

        The row `r` is modified in place. The bumped-out entry,
        if it exists, is returned.

        .. WARNING::

            This method only changes `r` if bumping occurs.
            Appending `j` to the end of the row should be done
            by the caller.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleDualRSK
            sage: r = [1, 3, 4, 5]
            sage: j = RuleDualRSK().insertion(4, r); j
            4
            sage: r
            [1, 3, 4, 5]
            sage: r = [1, 2, 3, 6, 7]
            sage: j = RuleDualRSK().insertion(4, r); j
            6
            sage: r
            [1, 2, 3, 4, 7]
            sage: r = [1, 3]
            sage: j = RuleDualRSK().insertion(4, r); j is None
            True
            sage: r
            [1, 3]
        """
        if r[-1] < j:
            return None  # j needs to be added at the end of the row.
        # Figure out where to insert j into the row r.  The
        # bisect command returns the position of the least
        # element of r greater or equal to j.  We will call it y.
        y_pos = bisect_left(r, j)
        # Switch j and y
        j, r[y_pos] = r[y_pos], j
        return j

    def reverse_insertion(self, x, row):
        r"""
        Reverse bump the row ``row`` of the current insertion tableau
        with the number ``x`` using dual RSK insertion.

        The row ``row`` is modified in place. The bumped-out entry
        is returned.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleDualRSK
            sage: r = [1, 2, 4, 6, 7]
            sage: x = RuleDualRSK().reverse_insertion(6, r); r
            [1, 2, 4, 6, 7]
            sage: x
            6
            sage: r = [1, 2, 4, 5, 7]
            sage: x = RuleDualRSK().reverse_insertion(6, r); r
            [1, 2, 4, 6, 7]
            sage: x
            5
        """
        y_pos = bisect_right(row, x) - 1
        # switch x and y
        x, row[y_pos] = row[y_pos], x
        return x

    def _backward_format_output(self, lower_row, upper_row, output, 
                                p_is_standard, q_is_standard):
        r"""
        Return the final output of the ``RSK_inverse`` correspondence
        from the output of the corresponding ``backward_rule``.

        .. NOTE::

            The default implementation of ``backward_rule`` lists
            bumped-out entries in the order in which the reverse
            bumping happens, which is *opposite* to the order of the
            final output.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleDualRSK
            sage: RuleDualRSK()._backward_format_output([1, 2, 3, 4], None,
            ....:                                    'permutation', True, True)
            [4, 3, 2, 1]
            sage: RuleDualRSK()._backward_format_output([1, 2, 3, 4], None,
            ....:                                   'permutation', False, True)
            Traceback (most recent call last):
            ...
            TypeError: p must be standard to have a valid permutation as output
        """
        if q_is_standard and output == 'permutation':
            if not p_is_standard:
                raise TypeError("p must be standard to have a valid permutation as output")
            from sage.combinat.permutation import Permutation
            return Permutation(reversed(lower_row))
        else:
            return super(RuleDualRSK, self)._backward_format_output(lower_row, upper_row, output,
                                                                    p_is_standard, q_is_standard)

    def _forward_format_output(self, p, q, check_standard):
        r"""
        Return final output of the ``RSK`` (here, dual RSK)
        correspondence from the output of the corresponding
        ``forward_rule``.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleDualRSK
            sage: isinstance(RuleDualRSK()._forward_format_output([[1,2,3,4,5]],
            ....:                    [[1,2,3,4,5]], True)[0], StandardTableau)
            True
            sage: isinstance(RuleDualRSK()._forward_format_output([[1,2,3,4,5]],
            ....:                            [[1,2,3,4,5]], False)[0], Tableau)
            True
            sage: isinstance(RuleDualRSK()._forward_format_output([[1,1,1,3,7]],
            ....:                             [[1,2,3,4,5]], True)[0], Tableau)
            True
        """
        from sage.combinat.tableau import Tableau, StandardTableau, SemistandardTableau

        if len(p) == 0:
            return [StandardTableau([]), StandardTableau([])]

        if check_standard:
            try:
                P = StandardTableau(p)
            except ValueError:
                P = Tableau(p)
            try:
                Q = StandardTableau(q)
            except ValueError:
                Q = SemistandardTableau(q)
            return [P, Q]
        return [Tableau(p), SemistandardTableau(q)]


class RuleCoRSK(RuleRSK):
    r"""
    A rule modeling the CoRSK insertion.

    CoRSK insertion differs from classical RSK insertion in the
    following ways:

    * The input (in terms of biwords) is no longer a biword,
      but rather a strict cobiword -- i.e., a pair of two lists
      `[a_1, a_2, \ldots, a_n]` and `[b_1, b_2, \ldots, b_n]` that
      satisfy the strict inequalities
      `(a_1, b_1) \widetilde{<} (a_2, b_2) \widetilde{<} \cdots
      \widetilde{<} (a_n, b_n)`, where
      the binary relation `\widetilde{<}` on pairs of integers
      is defined by having `(u_1, v_1) \widetilde{<} (u_2, v_2)`
      if and only if either `u_1 < u_2` or (`u_1 = u_2` and
      `v_1 > v_2`).
      In terms of matrices, this means that the input is not an
      arbitrary matrix with nonnegative integer entries, but rather
      a `\{0, 1\}`-matrix (i.e., a matrix whose entries are `0`'s
      and `1`'s).

    * The output still consists of two tableaux `(P, Q)` of equal
      shapes, but rather than both of them being semistandard, now
      `Q` is row-strict (i.e., its transpose is semistandard) while
      `P` is semistandard.

    Bumping proceeds in the same way as for RSK insertion.

    The RSK and CoRSK algorithms agree for permutation matrices.

    For more information, see Section A.4 in [Ful1997]_ (specifically,
    construction (1d)) or the second solution to Exercise 2.7.12(a) in
    [GR2018v5sol]_.

    EXAMPLES::

        sage: RSK([1,2,5,3,1], insertion = RSK.rules.coRSK)
        [[[1, 1, 3], [2], [5]], [[1, 2, 3], [4], [5]]]
        sage: RSK(Word([2,3,3,2,1,3,2,3]), insertion = RSK.rules.coRSK)
        [[[1, 2, 2, 3, 3], [2, 3], [3]], [[1, 2, 3, 6, 8], [4, 7], [5]]]
        sage: RSK(Word([3,3,2,4,1]), insertion = RSK.rules.coRSK)
        [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]
        sage: from sage.combinat.rsk import to_matrix
        sage: RSK(to_matrix([1, 1, 3, 3, 4], [3, 2, 2, 1, 3]), insertion = RSK.rules.coRSK)
        [[[1, 2, 3], [2], [3]], [[1, 3, 4], [1], [3]]]

    Using coRSK insertion with a `\{0, 1\}`-matrix::

        sage: RSK(matrix([[0,1],[1,0]]), insertion = RSK.rules.coRSK)
        [[[1], [2]], [[1], [2]]]

    We can also give it something looking like a matrix::

        sage: RSK([[0,1],[1,0]], insertion = RSK.rules.coRSK)
        [[[1], [2]], [[1], [2]]]

    We can also use the inverse correspondence::

        sage: RSK_inverse(*RSK([1, 2, 2, 2], [2, 3, 2, 1],
        ....:         insertion=RSK.rules.coRSK),insertion=RSK.rules.coRSK)
        [[1, 2, 2, 2], [2, 3, 2, 1]]
        sage: P,Q = RSK([1, 2, 2, 2], [2, 3, 2, 1],insertion=RSK.rules.coRSK)
        sage: RSK_inverse(P, Q, insertion=RSK.rules.coRSK)
        [[1, 2, 2, 2], [2, 3, 2, 1]]

    When applied to two standard tableaux, backwards coRSK
    insertion behaves identically to the usual backwards RSK
    insertion::

        sage: t1 = Tableau([[1, 2, 5], [3], [4]])
        sage: t2 = Tableau([[1, 2, 3], [4], [5]])
        sage: RSK_inverse(t1, t2, insertion=RSK.rules.coRSK)
        [[1, 2, 3, 4, 5], [1, 4, 5, 3, 2]]
        sage: RSK_inverse(t1, t2, 'word', insertion=RSK.rules.coRSK)
        word: 14532
        sage: RSK_inverse(t1, t2, 'matrix', insertion=RSK.rules.coRSK)
        [1 0 0 0 0]
        [0 0 0 1 0]
        [0 0 0 0 1]
        [0 0 1 0 0]
        [0 1 0 0 0]
        sage: RSK_inverse(t1, t2, 'permutation', insertion=RSK.rules.coRSK)
        [1, 4, 5, 3, 2]
        sage: RSK_inverse(t1, t1, 'permutation', insertion=RSK.rules.coRSK)
        [1, 4, 3, 2, 5]
        sage: RSK_inverse(t2, t2, 'permutation', insertion=RSK.rules.coRSK)
        [1, 2, 5, 4, 3]
        sage: RSK_inverse(t2, t1, 'permutation', insertion=RSK.rules.coRSK)
        [1, 5, 4, 2, 3]

    For coRSK, the first tableau is semistandard while the second tableau
    is transpose semistandard::

        sage: p = Tableau([[1,2,2],[5]]); q = Tableau([[1,2,4],[3]])
        sage: ret = RSK_inverse(p, q, insertion=RSK.rules.coRSK); ret
        [[1, 2, 3, 4], [1, 5, 2, 2]]
        sage: RSK_inverse(p, q, 'word', insertion=RSK.rules.coRSK)
        word: 1522

    TESTS:

    Empty objects::

        sage: RSK(Permutation([]), insertion=RSK.rules.coRSK)
        [[], []]
        sage: RSK(Word([]), insertion=RSK.rules.coRSK)
        [[], []]
        sage: RSK(matrix([[]]), insertion=RSK.rules.coRSK)
        [[], []]
        sage: RSK([], [], insertion=RSK.rules.coRSK)
        [[], []]
        sage: RSK([[]], insertion=RSK.rules.coRSK)
        [[], []]

    Check that :func:`RSK_inverse` is the inverse of :func:`RSK` on the
    different types of inputs/outputs::

        sage: RSK_inverse(Tableau([]), Tableau([]),
        ....:                insertion=RSK.rules.coRSK)
        [[], []]
        sage: f = lambda p: RSK_inverse(*RSK(p, insertion=RSK.rules.coRSK),
        ....:                output='permutation', insertion=RSK.rules.coRSK)
        sage: all(p == f(p) for n in range(7) for p in Permutations(n))
        True
        sage: all(RSK_inverse(*RSK(w, insertion=RSK.rules.coRSK),
        ....:                 output='word', insertion=RSK.rules.coRSK) == w
        ....:     for n in range(4) for w in Words(5, n))
        True
        sage: from sage.combinat.integer_matrices import IntegerMatrices
        sage: M = IntegerMatrices([1,2,2,1], [3,1,1,1])
        sage: all(RSK_inverse(*RSK(m, insertion=RSK.rules.coRSK),
        ....:                output='matrix', insertion=RSK.rules.coRSK) == m
        ....:     for m in M if all(x in [0, 1] for x in m))
        True

        sage: n = ZZ.random_element(200)
        sage: p = Permutations(n).random_element()
        sage: True if p == f(p) else p
        True

    Checking that the tableaux should be of same shape::

        sage: RSK_inverse(Tableau([[1,2,3]]), Tableau([[1,2]]),
        ....:                          insertion=RSK.rules.dualRSK)
        Traceback (most recent call last):
        ...
        ValueError: p(=[[1, 2, 3]]) and q(=[[1, 2]]) must have the same shape

    Checking that the biword is a strict cobiword::

        sage: RSK([1,2,4,3], [1,2,3,4], insertion=RSK.rules.coRSK)
        Traceback (most recent call last):
        ...
        ValueError: invalid strict cobiword
        sage: RSK([1,2,3,3], [1,2,3,4], insertion=RSK.rules.coRSK)
        Traceback (most recent call last):
        ...
        ValueError: invalid strict cobiword
        sage: RSK([1,2,3,3], [1,2,3,3], insertion=RSK.rules.coRSK)
        Traceback (most recent call last):
        ...
        ValueError: invalid strict cobiword
    """
    def to_pairs(self, obj1=None, obj2=None, check=True):
        r"""
        Given a valid input for the coRSK algorithm, such as
        two `n`-tuples ``obj1`` `= [a_1, a_2, \ldots, a_n]`
        and ``obj2`` `= [b_1, b_2, \ldots, b_n]` forming a
        strict cobiword (i.e., satisfying
        `a_1 \leq a_2 \leq \cdots \leq a_n`, and if
        `a_i = a_{i+1}`, then `b_i > b_{i+1}`),
        or a `\{0, 1\}`-matrix ("rook placement"), or a
        single word, return the array
        `[(a_1, b_1), (a_2, b_2), \ldots, (a_n, b_n)]`.

        INPUT:

        - ``obj1, obj2`` -- anything representing a strict
          cobiword (see the doc of :meth:`forward_rule` for
          the encodings accepted).

        - ``check`` -- (default: ``True``) whether to check
          that ``obj1`` and ``obj2`` actually define a valid
          strict cobiword.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleCoRSK
            sage: list(RuleCoRSK().to_pairs([1, 2, 2, 2], [2, 3, 2, 1]))
            [(1, 2), (2, 3), (2, 2), (2, 1)]
            sage: RuleCoRSK().to_pairs([1, 2, 2, 2], [1, 2, 3, 3])
            Traceback (most recent call last):
            ...
            ValueError: invalid strict cobiword
            sage: m = Matrix(ZZ, 3, 2, [0,1,1,1,0,1]) ; m
            [0 1]
            [1 1]
            [0 1]
            sage: list(RuleCoRSK().to_pairs(m))
            [(1, 2), (2, 2), (2, 1), (3, 2)]
            sage: m = Matrix(ZZ, 3, 2, [0,1,1,0,0,2]) ; m
            [0 1]
            [1 0]
            [0 2]
            sage: RuleCoRSK().to_pairs(m)
            Traceback (most recent call last):
            ...
            ValueError: coRSK requires a {0, 1}-matrix
        """
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
                        for j, mult in reversed(list(enumerate(row))):
                            if mult > 1:
                                raise ValueError("coRSK requires a {0, 1}-matrix")
                            if mult > 0:
                                t.append(i+1)
                                b.append(j+1)
                    itr = zip(t, b)
                except TypeError:
                    itr = zip(range(1, len(obj1)+1), obj1)
        else:
            if check:
                if len(obj1) != len(obj2):
                    raise ValueError("the two arrays must be the same length")
                # Check it is a generalized permutation (i.e., biword):
                # that is, the pairs of corresponding entries of obj1
                # and obj2 are weakly increasing in lexicographic order.
                lt = 0
                lb = 0
                for t, b in zip(obj1, obj2):
                    if t < lt or (t == lt and b >= lb):
                        raise ValueError("invalid strict cobiword")
                    lt = t
                    lb = b
            itr = zip(obj1, obj2)
        return itr

    def _forward_format_output(self, p, q, check_standard):
        r"""
        Return final output of the ``RSK`` correspondence from the
        output of the corresponding ``forward_rule``.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleCoRSK
            sage: isinstance(RuleCoRSK()._forward_format_output([[1,2,3,4,5]],
            ....:                    [[1,2,3,4,5]], True)[0], StandardTableau)
            True
            sage: isinstance(RuleCoRSK()._forward_format_output([[1,2,3,4,5]],
            ....:                            [[1,2,3,4,5]], False)[0], SemistandardTableau)
            True
            sage: isinstance(RuleCoRSK()._forward_format_output([[1,1,1,3,7]],
            ....:                             [[1,2,3,4,5]], True)[1], Tableau)
            True
            sage: isinstance(RuleCoRSK()._forward_format_output([[1,1,1,3,7]],
            ....:                             [[1,2,3,4,5]], False)[1], Tableau)
            True
        """
        from sage.combinat.tableau import SemistandardTableau, StandardTableau, Tableau

        if check_standard:
            try:
                P = StandardTableau(p)
            except ValueError:
                P = SemistandardTableau(p)
            try:
                Q = StandardTableau(q)
            except ValueError:
                Q = Tableau(q)
            return [P, Q]
        return [SemistandardTableau(p), Tableau(q)]

    def backward_rule(self, p, q, output):
        r"""
        Return the strict cobiword obtained by applying reverse
        coRSK insertion to a pair of tableaux ``(p, q)``.

        INPUT:

        - ``p``, ``q`` -- two tableaux of the same shape.

        - ``output`` -- (default: ``'array'``) if ``q`` is row-strict:

          - ``'array'`` -- as a two-line array (i.e. strict cobiword)
          - ``'matrix'`` -- as a `\{0, 1\}`-matrix

          and if ``q`` is standard, we can have the output:

          - ``'word'`` -- as a word

          and additionally if ``p`` is standard, we can also have the output:

          - ``'permutation'`` -- as a permutation

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleCoRSK
            sage: t1 = Tableau([[1, 1, 2], [2, 3], [4]])
            sage: t2 = Tableau([[1, 4, 5], [1, 4], [2]])
            sage: RuleCoRSK().backward_rule(t1, t2, 'array')
            [[1, 1, 2, 4, 4, 5], [4, 2, 1, 3, 1, 2]]
        """
        # Make a copy of p since this is destructive to it
        p_copy = [list(row) for row in p]

        if q.is_standard():
            rev_word = []  # This will be our word in reverse
            d = {qij: i for i, Li in enumerate(q) for qij in Li}
            # d is now a dictionary which assigns to each integer k the
            # number of the row of q containing k.

            for key in sorted(d, reverse=True):
                # Delete last entry from i-th row of p_copy
                i = d[key]
                x = p_copy[i].pop()  # Always the right-most entry
                for row in reversed(p_copy[:i]):
                    x = self.reverse_insertion(x, row)
                rev_word.append(x)
            return self._backward_format_output(rev_word, None, output, p.is_standard(), True)

        upper_row = []
        lower_row = []
        # upper_row and lower_row will be the upper and lower rows of the
        # strict cobiword we get as a result, but both reversed.
        d = {}
        for row, Li in enumerate(q):
            for val in Li:
                if val in d:
                    d[val].append(row)
                else:
                    d[val] = [row]
        # d is now a dictionary which assigns to each integer k the
        # list of the rows of q containing k.
        for value, row_list in sorted(d.items(), reverse=True, key=lambda x: x[0]):
            for i in sorted(row_list, reverse=True):
                x = p_copy[i].pop()  # Always the right-most entry
                for row in reversed(p_copy[:i]):
                    x = self.reverse_insertion(x, row)
                lower_row.append(x)
                upper_row.append(value)
        return self._backward_format_output(lower_row, upper_row, output, p.is_standard(), False)


class InsertionRules(object):
    r"""
    Catalog of rules for RSK-like insertion algorithms.
    """
    RSK = RuleRSK
    EG = RuleEG
    Hecke = RuleHecke
    dualRSK = RuleDualRSK
    coRSK = RuleCoRSK

#####################################################################

def RSK(obj1=None, obj2=None, insertion=InsertionRules.RSK, check_standard=False, **options):
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
    simply appended to the row and the algorithm terminates at this point.

    A *generalized permutation* (or *biword*) is a list
    `((j_0, k_0), (j_1, k_1), \ldots, (j_{\ell-1}, k_{\ell-1}))`
    of pairs such that the letters `j_0, j_1, \ldots, j_{\ell-1}`
    are weakly increasing (that is,
    `j_0 \leq j_1 \leq \cdots \leq j_{\ell-1}`), whereas the letters
    `k_i` satisfy `k_i \leq k_{i+1}` whenever `j_i = j_{i+1}`.
    The `\ell`-tuple `(j_0, j_1, \ldots, j_{\ell-1})` is called the
    *top line* of this generalized permutation,
    whereas the `\ell`-tuple `(k_0, k_1, \ldots, k_{\ell-1})` is
    called its *bottom line*.

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
    considering the top row to be `(1, 2, \ldots, n)` where `n` is the
    length of `w`.

    The optional argument ``insertion`` allows to specify an alternative
    insertion procedure to be used instead of the standard
    Robinson-Schensted-Knuth insertion.

    INPUT:

    - ``obj1, obj2`` -- can be one of the following:

      - a word in an ordered alphabet (in this case, ``obj1`` is said
        word, and ``obj2`` is ``None``)
      - an integer matrix
      - two lists of equal length representing a generalized permutation
        (namely, the lists `(j_0, j_1, \ldots, j_{\ell-1})` and
        `(k_0, k_1, \ldots, k_{\ell-1})` represent the generalized
        permutation
        `((j_0, k_0), (j_1, k_1), \ldots, (j_{\ell-1}, k_{\ell-1}))`)
      - any object which has a method ``_rsk_iter()`` which returns an
        iterator over the object represented as generalized permutation or
        a pair of lists (in this case, ``obj1`` is said object,
        and ``obj2`` is ``None``).

    - ``insertion`` -- (default: ``RSK.rules.RSK``) the following types
      of insertion are currently supported:

      - ``RSK.rules.RSK`` (or ``'RSK'``) -- Robinson-Schensted-Knuth
        insertion (:class:`~sage.combinat.rsk.RuleRSK`)
      - ``RSK.rules.EG`` (or ``'EG'``) -- Edelman-Greene insertion
        (only for reduced words of permutations/elements of a type `A`
        Coxeter group) (:class:`~sage.combinat.rsk.RuleEG`)
      - ``RSK.rules.Hecke`` (or ``'hecke'``) -- Hecke insertion (only
        guaranteed for generalized permutations whose top row is strictly
        increasing) (:class:`~sage.combinat.rsk.RuleHecke`)
      - ``RSK.rules.dualRSK`` (or ``'dualRSK'``) -- Dual RSK insertion 
        (only for strict biwords) (:class:`~sage.combinat.rsk.RuleDualRSK`)
      - ``RSK.rules.coRSK`` (or ``'coRSK'``) -- CoRSK insertion (only 
        for strict cobiwords) (:class:`~sage.combinat.rsk.RuleCoRSK`)

    - ``check_standard`` -- (default: ``False``) check if either of the
      resulting tableaux is a standard tableau, and if so, typecast it
      as such

    For precise information about constraints on the input and output,
    as well as the definition of the algorithm (if it is not standard
    RSK), see the particular :class:`~sage.combinat.rsk.Rule` class.

    EXAMPLES:

    If we only input one row, it is understood that the top row
    should be `(1, 2, \ldots, n)`::

        sage: RSK([3,3,2,4,1])
        [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]
        sage: RSK(Word([3,3,2,4,1]))
        [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]
        sage: RSK(Word([2,3,3,2,1,3,2,3]))
        [[[1, 2, 2, 3, 3], [2, 3], [3]], [[1, 2, 3, 6, 8], [4, 7], [5]]]

    We can provide a generalized permutation::

        sage: RSK([1, 2, 2, 2], [2, 1, 1, 2])
        [[[1, 1, 2], [2]], [[1, 2, 2], [2]]]
        sage: RSK(Word([1,1,3,4,4]), [1,4,2,1,3])
        [[[1, 1, 3], [2], [4]], [[1, 1, 4], [3], [4]]]
        sage: RSK([1,3,3,4,4], Word([6,2,2,1,7]))
        [[[1, 2, 7], [2], [6]], [[1, 3, 4], [3], [4]]]

    We can provide a matrix::

        sage: RSK(matrix([[0,1],[2,1]]))
        [[[1, 1, 2], [2]], [[1, 2, 2], [2]]]

    We can also provide something looking like a matrix::

        sage: RSK([[0,1],[2,1]])
        [[[1, 1, 2], [2]], [[1, 2, 2], [2]]]

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
        sage: RSK(Word([]), insertion=RSK.rules.EG)
        [[], []]
        sage: RSK(Word([]), insertion=RSK.rules.Hecke)
        [[], []]

    """
    if isinstance(insertion, str):
        if insertion == 'RSK':
            insertion = RSK.rules.RSK
        elif insertion == 'EG':
            insertion = RSK.rules.EG
        elif insertion == 'hecke':
            insertion = RSK.rules.Hecke
        elif insertion == 'dualRSK':
            insertion = RSK.rules.dualRSK
        elif insertion == 'coRSK':
            insertion = RSK.rules.coRSK
        else:
            raise ValueError("invalid input")

    rule = insertion()
    if not isinstance(rule, Rule):
        raise TypeError("the insertion must be an instance of Rule")

    if obj1 is None and obj2 is None:
        if 'matrix' in options:
            obj1 = matrix(options['matrix'])
        else:
            raise ValueError("invalid input")

    if is_Matrix(obj1):
        obj1 = obj1.rows()

    output = rule.forward_rule(obj1, obj2, check_standard)
    return output


robinson_schensted_knuth = RSK
RSK.rules = InsertionRules

def RSK_inverse(p, q, output='array', insertion=InsertionRules.RSK):
    r"""
    Return the generalized permutation corresponding to the pair of
    tableaux `(p, q)` under the inverse of the Robinson-Schensted-Knuth
    correspondence.

    For more information on the bijection, see :func:`RSK`.

    INPUT:

    - ``p``, ``q`` -- two semi-standard tableaux of the same shape, or
      (in the case when Hecke insertion is used) an increasing tableau and
      a set-valued tableau of the same shape (see the note below for the
      format of the set-valued tableau)

    - ``output`` -- (default: ``'array'``) if ``q`` is semi-standard:

      - ``'array'`` -- as a two-line array (i.e. generalized permutation or
        biword)
      - ``'matrix'`` -- as an integer matrix

      and if ``q`` is standard, we can also have the output:

      - ``'word'`` -- as a word

      and additionally if ``p`` is standard, we can also have the output:

      - ``'permutation'`` -- as a permutation

    - ``insertion`` -- (default: ``RSK.rules.RSK``) the insertion algorithm
      used in the bijection. Currently the following are supported:

      - ``RSK.rules.RSK`` (or ``'RSK'``) -- Robinson-Schensted-Knuth
        insertion (:class:`~sage.combinat.rsk.RuleRSK`)
      - ``RSK.rules.EG`` (or ``'EG'``) -- Edelman-Greene insertion
        (only for reduced words of permutations/elements of a type `A`
        Coxeter group) (:class:`~sage.combinat.rsk.RuleEG`)
      - ``RSK.rules.Hecke`` (or ``'hecke'``) -- Hecke insertion (only
        guaranteed for generalized permutations whose top row is strictly
        increasing) (:class:`~sage.combinat.rsk.RuleHecke`)
      - ``RSK.rules.dualRSK`` (or ``'dualRSK'``) -- Dual RSK insertion 
        (only for strict biwords) (:class:`~sage.combinat.rsk.RuleDualRSK`)
      - ``RSK.rules.coRSK`` (or ``'coRSK'``) -- CoRSK insertion (only 
        for strict cobiwords) (:class:`~sage.combinat.rsk.RuleCoRSK`)

    For precise information about constraints on the input and
    output, see the particular :class:`~sage.combinat.rsk.Rule` class.

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

    Using Hecke insertion::

        sage: w = [5, 4, 3, 1, 4, 2, 5, 5]
        sage: pq = RSK(w, insertion=RSK.rules.Hecke)
        sage: RSK_inverse(*pq, insertion=RSK.rules.Hecke, output='list')
        [5, 4, 3, 1, 4, 2, 5, 5]

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
        sage: all(RSK_inverse(*RSK(w), output='word') == w for n in range(4)
        ....:                                            for w in Words(5, n))
        True
        sage: from sage.combinat.integer_matrices import IntegerMatrices
        sage: M = IntegerMatrices([1,2,2,1], [3,1,1,1])
        sage: all(RSK_inverse(*RSK(m), output='matrix') == m for m in M)
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
    if isinstance(insertion, str):
        if insertion == 'RSK':
            insertion = RSK.rules.RSK
        elif insertion == 'EG':
            insertion = RSK.rules.EG
        elif insertion == 'hecke':
            insertion = RSK.rules.Hecke
        elif insertion == 'dualRSK':
            insertion = RSK.rules.dualRSK
        elif insertion == 'coRSK':
            insertion = RSK.rules.coRSK
        else:
            raise ValueError("invalid input")

    rule = insertion()
    if not isinstance(rule, Rule):
        raise TypeError("the insertion must be an instance of Rule")

    if p.shape() != q.shape():
        raise ValueError("p(=%s) and q(=%s) must have the same shape" %(p, q))

    answer = rule.backward_rule(p, q, output)
    return answer

robinson_schensted_knuth_inverse = RSK_inverse

def to_matrix(t, b):
    r"""
    Return the integer matrix corresponding to a two-line array.

    INPUT:

    - ``t`` -- the top row of the array
    - ``b`` -- the bottom row of the array

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
        raise ValueError("the two arrays must be the same length")

    # Build the dictionary of entries since the matrix
    #   is typically (very) sparse
    entries = {}
    for i in range(n):
        pos = (t[i]-1, b[i]-1)
        if pos in entries:
            entries[pos] += 1
        else:
            entries[pos] = 1
    return matrix(entries, sparse=True)

