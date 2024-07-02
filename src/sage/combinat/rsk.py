# sage.doctest: needs sage.combinat sage.modules
r"""
Robinson-Schensted-Knuth correspondence

AUTHORS:

- Travis Scrimshaw (2012-12-07): Initial version
- Chaman Agrawal (2019-06-24): Refactoring on the Rule class
- Matthew Lancellotti (2018): initial version of super RSK
- Jianping Pan, Wencin Poh, Anne Schilling (2020-08-31): initial version of RuleStar

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

- RSK insertion (:class:`~sage.combinat.rsk.RuleRSK`).
- Edelman-Greene insertion (:class:`~sage.combinat.rsk.RuleEG`), an algorithm
  defined in [EG1987]_ Definition 6.20 (where it is referred to as
  Coxeter-Knuth insertion).
- Hecke RSK algorithm (:class:`~sage.combinat.rsk.RuleHecke`) , defined
  using the Hecke insertion studied in [BKSTY06]_ (but using rows instead
  of columns).
- Dual RSK insertion (:class:`~sage.combinat.rsk.RuleDualRSK`).
- CoRSK insertion (:class:`~sage.combinat.rsk.RuleCoRSK`), defined in
  [GR2018v5sol]_.
- Super RSK insertion (:class:`~sage.combinat.rsk.RuleSuperRSK`), a
  combination of row and column insertions defined in [Muth2019]_.
- Star insertion (:class:`~sage.combinat.rsk.RuleStar`), defined in [MPPS2020]_.

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
   :doi:`10.1016/0001-8708(87)90063-6`

.. [BKSTY06] \A. Buch, A. Kresch, M. Shimozono, H. Tamvakis, and A. Yong.
   *Stable Grothendieck polynomials and* `K`-*theoretic factor sequences*.
   Math. Ann. **340** Issue 2, (2008), pp. 359--382.
   :arxiv:`math/0601514v1`.

.. [GR2018v5sol] Darij Grinberg, Victor Reiner.
   *Hopf Algebras In Combinatorics*,
   :arxiv:`1409.8356v5`, available with solutions at
   https://arxiv.org/src/1409.8356v5/anc/HopfComb-v73-with-solutions.pdf
"""

# *****************************************************************************
#       Copyright (C) 2012,2019 Travis Scrimshaw <tcscrims at gmail.com>
#                          2019 Chaman Agrawal <chaman.ag at gmail.com>
#                          2019 Matthew Lancellotti <mareoraft at gmail.com>
#                          2020 Jianping Pan <jppan at math dot ucdavis dot edu>
#                               Wencin Poh <wpoh at ucdavis dot edu>
#                               Anne Schilling <anne at math dot ucdavis dot edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from bisect import bisect_left, bisect_right

from sage.misc.lazy_import import lazy_import
from sage.rings.integer_ring import ZZ
from sage.structure.element import Matrix
from sage.structure.unique_representation import UniqueRepresentation

lazy_import('sage.matrix.constructor', 'matrix')


class Rule(UniqueRepresentation):
    r"""
    Generic base class for an insertion rule for an RSK-type correspondence.

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

        - ``obj1``, ``obj2`` -- anything representing a biword
          (see the doc of :meth:`forward_rule` for the
          encodings accepted)

        - ``check`` -- boolean (default: ``True``); whether to check that
          ``obj1`` and ``obj2`` actually define a valid biword

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

        - ``obj1``, ``obj2`` -- can be one of the following ways to
          represent a generalized permutation (or, equivalently,
          biword):

          - two lists ``obj1`` and ``obj2`` of equal length,
            to be interpreted as the top row and the bottom row of
            the biword

          - a matrix ``obj1`` of nonnegative integers, to be
            interpreted as the generalized permutation in matrix
            form (in this case, ``obj2`` is ``None``)

          - a word ``obj1`` in an ordered alphabet, to be
            interpreted as the bottom row of the biword (in this
            case, ``obj2`` is ``None``; the top row of the biword
            is understood to be `(1, 2, \ldots, n)` by default)

          - any object ``obj1`` which has a method ``_rsk_iter()``,
            as long as this method returns an iterator yielding
            pairs of numbers, which then are interperted as top
            entries and bottom entries in the biword (in this case,
            ``obj2`` is ``None``)

        - ``check_standard`` -- boolean (default: ``False``); check if either
          of the resulting tableaux is a standard tableau, and if so, typecast
          it as such

        - ``check`` -- boolean (default: ``True``); whether to check
          that ``obj1`` and ``obj2`` actually define a valid
          biword

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

        - ``p``, ``q`` -- two tableaux of the same shape

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
            raise ValueError("q(=%s) must be a semistandard tableau" % q)

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
            sage: isinstance(RuleRSK()._forward_format_output([[1, 2, 3, 4, 5]],
            ....:            [[1, 2, 3, 4, 5]], True)[0], StandardTableau)
            True
            sage: isinstance(RuleRSK()._forward_format_output([[1, 2, 3, 4, 5]],
            ....:            [[1, 2, 3, 4, 5]], False)[0], SemistandardTableau)
            True
            sage: isinstance(RuleRSK()._forward_format_output([[1, 1, 1, 3, 7]],
            ....:            [[1, 2, 3, 4, 5]], True)[0], SemistandardTableau)
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
            sage: RuleRSK()._backward_format_output([1, 2, 3, 4], None, 'matrix', False, True)
            [0 0 0 1]
            [0 0 1 0]
            [0 1 0 0]
            [1 0 0 0]
            sage: RuleRSK()._backward_format_output([1, 2, 3, 4], None, 'word', False, True)
            word: 4321
            sage: RuleRSK()._backward_format_output([3, 2, 1, 1], [2, 1, 1, 1], 'array', False, False)
            [[1, 1, 1, 2], [1, 1, 2, 3]]
            sage: RuleRSK()._backward_format_output([3, 2, 1, 1], [2, 1, 1, 1], 'matrix', False, False)
            [2 1 0]
            [0 0 1]
            sage: RuleRSK()._backward_format_output([1, 2, 3, 4], None, 'word', False, False)
            Traceback (most recent call last):
            ...
            TypeError: q must be standard to have a word as valid output
            sage: RuleRSK()._backward_format_output([1, 2, 3, 4], None, 'random_type', False, False)
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
                    "q must be standard to have a %s as valid output" % output)
            raise ValueError("invalid output option")


class RuleRSK(Rule):
    r"""
    Rule for the classical Robinson-Schensted-Knuth insertion.

    See :func:`RSK` for the definition of this operation.

    EXAMPLES::

        sage: RSK([1, 2, 2, 2], [2, 1, 1, 2], insertion=RSK.rules.RSK)
        [[[1, 1, 2], [2]], [[1, 2, 2], [2]]]
        sage: p = Tableau([[1,2,2],[2]]); q = Tableau([[1,3,3],[2]])
        sage: RSK_inverse(p, q, insertion=RSK.rules.RSK)
        [[1, 2, 3, 3], [2, 1, 2, 2]]
    """

    def insertion(self, j, r):
        r"""
        Insert the letter ``j`` from the second row of the biword
        into the row `r` using classical Schensted insertion,
        if there is bumping to be done.

        The row `r` is modified in place if bumping occurs. The bumped-out
        entry, if it exists, is returned.

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
        return super()._backward_format_output(lower_row, upper_row, output,
                                               p_is_standard, q_is_standard)


class RuleEG(Rule):
    r"""
    Rule for Edelman-Greene insertion.

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

        The row `r` is modified in place if bumping occurs. The bumped-out
        entry, if it exists, is returned.

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
            return super()._backward_format_output(lower_row, upper_row, output,
                                                   p_is_standard, q_is_standard)


class RuleHecke(Rule):
    r"""
    Rule for Hecke insertion.

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

        - ``obj1``, ``obj2`` -- can be one of the following ways to
          represent a generalized permutation (or, equivalently,
          biword):

          - two lists ``obj1`` and ``obj2`` of equal length,
            to be interpreted as the top row and the bottom row of
            the biword

          - a word ``obj1`` in an ordered alphabet, to be
            interpreted as the bottom row of the biword (in this
            case, ``obj2`` is ``None``; the top row of the biword
            is understood to be `(1, 2, \ldots, n)` by default)

        - ``check_standard`` -- boolean (default: ``False``); check if either
          of the resulting tableaux is a standard tableau, and if so, typecast
          it as such

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

        - ``p``, ``q`` -- two tableaux of the same shape

        - ``output`` -- (default: ``'array'``) if ``q`` is semi-standard:

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

        The row `r` is modified in place if bumping occurs. The bumped-out
        entry, if it exists, is returned.

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
    Rule for dual RSK insertion.

    Dual RSK insertion differs from classical RSK insertion in the
    following ways:

    * The input (in terms of biwords) is no longer an arbitrary biword,
      but rather a strict biword (i.e., a pair of two lists
      `[a_1, a_2, \ldots, a_n]` and `[b_1, b_2, \ldots, b_n]` that
      satisfy the strict inequalities
      `(a_1, b_1) < (a_2, b_2) < \cdots < (a_n, b_n)` in
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
        biword (i.e., satisfying `a_1 \leq a_2 \leq \cdots \leq a_n`,
        and if `a_i = a_{i+1}`, then `b_i < b_{i+1}`) or a
        `\{0, 1\}`-matrix ("rook placement"), or a single word, return
        the array `[(a_1, b_1), (a_2, b_2), \ldots, (a_n, b_n)]`.

        INPUT:

        - ``obj1``, ``obj2`` -- anything representing a strict biword
          (see the doc of :meth:`forward_rule` for the
          encodings accepted)

        - ``check`` -- boolean (default: ``True``); whether to check
          that ``obj1`` and ``obj2`` actually define a valid
          strict biword

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

        The row `r` is modified in place if bumping occurs. The bumped-out
        entry, if it exists, is returned.

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
            return super()._backward_format_output(lower_row, upper_row, output,
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
    Rule for coRSK insertion.

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

    The RSK and coRSK algorithms agree for permutation matrices.

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

        - ``obj1``, ``obj2`` -- anything representing a strict
          cobiword (see the doc of :meth:`forward_rule` for
          the encodings accepted)

        - ``check`` -- boolean (default: ``True``); whether to check
          that ``obj1`` and ``obj2`` actually define a valid
          strict cobiword

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

        - ``p``, ``q`` -- two tableaux of the same shape

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


class RuleSuperRSK(RuleRSK):
    r"""
    Rule for super RSK insertion.

    Super RSK is based on `\epsilon`-insertion, a combination of
    row and column classical RSK insertion.

    Super RSK insertion differs from the classical RSK insertion in the
    following ways:

    * The input (in terms of biwords) is no longer an arbitrary biword,
      but rather a restricted super biword (i.e., a pair of two lists
      `[a_1, a_2, \ldots, a_n]` and `[b_1, b_2, \ldots, b_n]` that
      contains entries with even and odd parity and pairs with mixed
      parity entries do not repeat).

    * The output still consists of two tableaux `(P, Q)` of equal
      shapes, but rather than both of them being semistandard, now
      they are semistandard super tableaux.

    * The main difference is in the way bumping works. Instead of having
      only row bumping super RSK uses `\epsilon`-insertion, a combination
      of classical RSK bumping along the rows and a dual RSK like bumping
      (i.e. when a number `k_i` is inserted into the `i`-th row of `P`, it
      bumps out the first integer greater **or equal to** `k_i` in the column)
      along the column.

    EXAMPLES::

        sage: RSK([1], [1], insertion='superRSK')
        [[[1]], [[1]]]
        sage: RSK([1, 2], [1, 3], insertion='superRSK')
        [[[1, 3]], [[1, 2]]]
        sage: RSK([1, 2, 3], [1, 3, "3p"], insertion='superRSK')
        [[[1, 3], [3']], [[1, 2], [3]]]
        sage: RSK([1, 3, "3p", "2p"], insertion='superRSK')
        [[[1, 3', 3], [2']], [[1', 1, 2'], [2]]]
        sage: RSK(["1p", "2p", 2, 2, "3p", "3p", 3, 3],
        ....:     ["1p", 1, "2p", 2, "3p", "3p", "3p", 3], insertion='superRSK')
        [[[1', 2, 3', 3], [1, 3'], [2'], [3']], [[1', 2, 3', 3], [2', 3'], [2], [3]]]
        sage: P = SemistandardSuperTableau([[1, '3p', 3], ['2p']])
        sage: Q = SemistandardSuperTableau([['1p', 1, '2p'], [2]])
        sage: RSK_inverse(P, Q, insertion=RSK.rules.superRSK)
        [[1', 1, 2', 2], [1, 3, 3', 2']]

    We apply super RSK on Example 5.1 in [Muth2019]_::

        sage: P,Q = RSK(["1p", "2p", 2, 2, "3p", "3p", 3, 3],
        ....:           ["3p", 1, 2, 3, "3p", "3p", "2p", "1p"], insertion='superRSK')
        sage: (P, Q)
        ([[1', 2', 3', 3], [1, 2, 3'], [3']], [[1', 2, 2, 3'], [2', 3, 3], [3']])
        sage: ascii_art((P, Q))
        (  1' 2' 3'  3   1'  2  2 3' )
        (   1  2 3'      2'  3  3    )
        (  3'         ,  3'          )
        sage: RSK_inverse(P, Q, insertion=RSK.rules.superRSK)
        [[1', 2', 2, 2, 3', 3', 3, 3], [3', 1, 2, 3, 3', 3', 2', 1']]

    Example 6.1 in [Muth2019]_::

        sage: P,Q = RSK(["1p", "2p", 2, 2, "3p", "3p", 3, 3],
        ....:           ["3p", 1, 2, 3, "3p", "3p", "2p", "1p"], insertion='superRSK')
        sage: ascii_art((P, Q))
        (  1' 2' 3'  3   1'  2  2 3' )
        (   1  2 3'      2'  3  3    )
        (  3'         ,  3'          )
        sage: RSK_inverse(P, Q, insertion=RSK.rules.superRSK)
        [[1', 2', 2, 2, 3', 3', 3, 3], [3', 1, 2, 3, 3', 3', 2', 1']]

        sage: P,Q = RSK(["1p", 1, "2p", 2, "3p", "3p", "3p", 3],
        ....:           [3, "2p", 3, 2, "3p", "3p", "1p", 2], insertion='superRSK')
        sage: ascii_art((P, Q))
        (  1'  2  2 3'   1' 2' 3'  3 )
        (  2'  3  3       1  2 3'    )
        (  3'         ,  3'          )
        sage: RSK_inverse(P, Q, insertion=RSK.rules.superRSK)
        [[1', 1, 2', 2, 3', 3', 3', 3], [3, 2', 3, 2, 3', 3', 1', 2]]

    Let us now call the inverse correspondence::

        sage: P, Q = RSK([1, 2, 2, 2], [2, 1, 2, 3],
        ....:            insertion=RSK.rules.superRSK)
        sage: RSK_inverse(P, Q, insertion=RSK.rules.superRSK)
        [[1, 2, 2, 2], [2, 1, 2, 3]]

    When applied to two tableaux with only even parity elements, reverse super
    RSK insertion behaves identically to the usual reversel RSK insertion::

        sage: t1 = Tableau([[1, 2, 5], [3], [4]])
        sage: t2 = Tableau([[1, 2, 3], [4], [5]])
        sage: RSK_inverse(t1, t2, insertion=RSK.rules.RSK)
        [[1, 2, 3, 4, 5], [1, 4, 5, 3, 2]]
        sage: t1 = SemistandardSuperTableau([[1, 2, 5], [3], [4]])
        sage: t2 = SemistandardSuperTableau([[1, 2, 3], [4], [5]])
        sage: RSK_inverse(t1, t2, insertion=RSK.rules.superRSK)
        [[1, 2, 3, 4, 5], [1, 4, 5, 3, 2]]

    TESTS:

    Empty objects::

        sage: RSK(Word([]), insertion=RSK.rules.superRSK)
        [[], []]
        sage: RSK([], [], insertion=RSK.rules.superRSK)
        [[], []]

    Check that :func:`RSK_inverse` is the inverse of :func:`RSK` on the
    different types of inputs/outputs::

        sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
        sage: RSK_inverse(SemistandardSuperTableau([]),
        ....:             SemistandardSuperTableau([]), insertion=RSK.rules.superRSK)
        [[], []]
        sage: f = lambda p: RSK_inverse(*RSK(p, insertion=RSK.rules.superRSK),
        ....:                           insertion=RSK.rules.superRSK)
        sage: all(p == f(p)[1] for n in range(5) for p in Permutations(n))
        True

        sage: SST = StandardSuperTableaux([3,2,1])
        sage: f = lambda P, Q: RSK(*RSK_inverse(P, Q, insertion=RSK.rules.superRSK),
        ....:                      insertion=RSK.rules.superRSK)
        sage: all([P, Q] == f(P, Q) for n in range(7) for la in Partitions(n)
        ....:     for P in StandardSuperTableaux(la) for Q in StandardSuperTableaux(la))
        True

    Checking that tableaux should be of same shape::

        sage: RSK_inverse(SemistandardSuperTableau([[1, 2, 3]]),
        ....:             SemistandardSuperTableau([[1, 2]]),
        ....:             insertion=RSK.rules.superRSK)
        Traceback (most recent call last):
        ...
        ValueError: p(=[[1, 2, 3]]) and q(=[[1, 2]]) must have the same shape
    """

    def to_pairs(self, obj1=None, obj2=None, check=True):
        r"""
        Given a valid input for the super RSK algorithm, such as
        two `n`-tuples ``obj1`` `= [a_1, a_2, \ldots, a_n]`
        and ``obj2`` `= [b_1, b_2, \ldots, b_n]` forming a restricted
        super biword (i.e., entries with even and odd parity and no
        repetition of corresponding pairs with mixed parity entries)
        return the array `[(a_1, b_1), (a_2, b_2), \ldots, (a_n, b_n)]`.

        INPUT:

        - ``obj1``, ``obj2`` -- anything representing a restricted super biword
          (see the doc of :meth:`forward_rule` for the
          encodings accepted)

        - ``check`` -- boolean (default: ``True``); whether to check
          that ``obj1`` and ``obj2`` actually define a valid
          restricted super biword

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleSuperRSK
            sage: list(RuleSuperRSK().to_pairs([2, '1p', 1],[1, 1, '1p']))
            [(2, 1), (1', 1), (1, 1')]
            sage: list(RuleSuperRSK().to_pairs([1, '1p', '2p']))
            [(1', 1), (1, 1'), (2', 2')]
            sage: list(RuleSuperRSK().to_pairs([1, 1], ['1p', '1p']))
            Traceback (most recent call last):
            ...
            ValueError: invalid restricted superbiword
        """
        from sage.combinat.shifted_primed_tableau import PrimedEntry
        # Initializing itr for itr = None case
        itr = None
        if obj2 is None:
            try:
                itr = obj1._rsk_iter()
            except AttributeError:
                # set recording list (obj1) to default value [1', 1, 2', 2, ...]
                obj2, obj1 = obj1, []
                a = ZZ.one() / ZZ(2)
                for i in range(len(obj2)):
                    obj1.append(a)
                    a = a + ZZ.one() / ZZ(2)
        else:
            if check:
                if len(obj1) != len(obj2):
                    raise ValueError("the two arrays must be the same length")
                mixed_parity = []
                # Check it is a restricted superbiword: that is,
                # the entries can have even or odd parity, but repetition of
                # the pairs of corresponding entries of obj1
                # and obj2 with mixed-parity is not allowed
                for t, b in zip(obj1, obj2):
                    if PrimedEntry(t).is_primed() != PrimedEntry(b).is_primed():
                        if (t, b) in mixed_parity:
                            raise ValueError("invalid restricted superbiword")
                        else:
                            mixed_parity.append((t, b))
        # Since the _rsk_iter() gives unprimed entries
        # We will create obj1 and obj2 from it.
        if itr:
            obj1, obj2 = [], []
            for i, j in itr:
                obj1.append(i)
                obj2.append(j)
        # Converting entries of obj1 and obj2 to PrimedEntry
        for i in range(len(obj1)):
            obj1[i] = PrimedEntry(obj1[i])
            obj2[i] = PrimedEntry(obj2[i])
        return zip(obj1, obj2)

    def _get_col(self, t, col_index):
        r"""
        Return the column as a list of a given tableau ``t`` (list of lists)
        at index ``col_index`` (indexing  starting from zero).

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleSuperRSK
            sage: t = [[1,2,3,4], [5,6,7,8], [9,10]];
            sage: RuleSuperRSK()._get_col(t, 0)
            [1, 5, 9]
            sage: RuleSuperRSK()._get_col(t, 2)
            [3, 7]
        """
        # t is the tableau (list of lists)
        # compute how many rows will contribute to the col
        num_rows_long_enough = 0
        for row in t:
            if len(row) > col_index:
                num_rows_long_enough += 1
            else:
                break
        # create the col
        col = [t[row_index][col_index] for row_index in range(num_rows_long_enough)]
        return col

    def _set_col(self, t, col_index, col):
        r"""
        Set the column of a given tableau ``t`` (list of lists) at
        index ``col_index`` (indexing  starting from zero) as  ``col``.

        .. NOTE::

            If ``length(col)`` is greater than the corresponding column in
            tableau ``t`` then only those rows of ``t`` will be set which
            have ``length(row) <= col_index``. Similarly if ``length(col)``
            is less than the corresponding column in tableau ``t`` then only
            those entries of the corresponding column in ``t`` which have row
            index less than ``length(col)`` will be set, rest will remain
            unchanged.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleSuperRSK
            sage: t = [[1,2,3,4], [5,6,7,8], [9,10]]
            sage: col = [1, 2, 3, 4]
            sage: RuleSuperRSK()._set_col(t, 0, col); t
            [[1, 2, 3, 4], [2, 6, 7, 8], [3, 10], [4]]
            sage: col = [1]
            sage: RuleSuperRSK()._set_col(t, 2, col); t
            [[1, 2, 1, 4], [2, 6, 7, 8], [3, 10], [4]]
        """
        # overwrite a column in tableau t (list of lists) with col
        for row_index, val in enumerate(col):
            # add a box/node if necessary
            if row_index == len(t):
                t.append([])
            if col_index == len(t[row_index]):
                t[row_index].append(None)
            # set value
            t[row_index][col_index] = val

    def forward_rule(self, obj1, obj2, check_standard=False, check=True):
        r"""
        Return a pair of tableaux obtained by applying forward
        insertion to the restricted super biword ``[obj1, obj2]``.

        INPUT:

        - ``obj1``, ``obj2`` -- can be one of the following ways to
          represent a generalized permutation (or, equivalently,
          biword):

          - two lists ``obj1`` and ``obj2`` of equal length,
            to be interpreted as the top row and the bottom row of
            the biword

          - a word ``obj1`` in an ordered alphabet, to be
            interpreted as the bottom row of the biword (in this
            case, ``obj2`` is ``None``; the top row of the biword
            is understood to be `(1, 2, \ldots, n)` by default)

          - any object ``obj1`` which has a method ``_rsk_iter()``,
            as long as this method returns an iterator yielding
            pairs of numbers, which then are interperted as top
            entries and bottom entries in the biword (in this case,
            ``obj2`` is ``None``)

        - ``check_standard`` -- boolean (default: ``False``); check if either
          of the resulting tableaux is a standard super tableau, and if so,
          typecast it as such

        - ``check`` -- boolean (default: ``True``); whether to check
          that ``obj1`` and ``obj2`` actually define a valid
          restricted super biword

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleSuperRSK
            sage: p, q = RuleSuperRSK().forward_rule([1, 2], [1, 3]); p
            [[1, 3]]
            sage: q
            [[1, 2]]
            sage: isinstance(p, SemistandardSuperTableau)
            True
            sage: isinstance(q, SemistandardSuperTableau)
            True
        """
        itr = self.to_pairs(obj1, obj2, check=check)
        p = []       # the "insertion" tableau
        q = []       # the "recording" tableau
        for i, j in itr:
            # loop
            row_index = -1
            col_index = -1
            epsilon = 1 if i.is_primed() else 0
            while True:
                if i.is_primed() == j.is_primed():
                    # row insertion
                    row_index += 1
                    if row_index == len(p):
                        p.append([j])
                        q.append([i])
                        break
                    else:
                        j1, col_index = self.insertion(j, p[row_index], epsilon=epsilon)
                        if j1 is None:
                            p[row_index].append(j)
                            q[row_index].append(i)
                            break
                        else:
                            j = j1
                else:
                    # column insertion
                    col_index += 1
                    if not p or col_index == len(p[0]):
                        self._set_col(p, col_index, [j])
                        self._set_col(q, col_index, [i])
                        break
                    else:
                        # retrieve column
                        c = self._get_col(p, col_index)
                        j1, row_index = self.insertion(j, c, epsilon=epsilon)
                        if j1 is None:
                            c.append(j)
                            self._set_col(p, col_index, c)
                            if col_index == 0:
                                q.append([])
                            q[row_index].append(i)
                            break
                        else:
                            j = j1
                        self._set_col(p, col_index, c)
        return self._forward_format_output(p, q, check_standard=check_standard)

    def insertion(self, j, r, epsilon=0):
        r"""
        Insert the letter ``j`` from the second row of the biword
        into the row ``r`` using dual RSK insertion or classical
        Schensted insertion depending on the value of ``epsilon``,
        if there is bumping to be done.

        The row `r` is modified in place if bumping occurs. The bumped-out
        entry, if it exists, is returned.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleSuperRSK
            sage: from bisect import bisect_left, bisect_right
            sage: r = [1, 3, 3, 3, 4]
            sage: j = 3
            sage: j, y_pos = RuleSuperRSK().insertion(j, r, epsilon=0); r
            [1, 3, 3, 3, 3]
            sage: j
            4
            sage: y_pos
            4
            sage: r = [1, 3, 3, 3, 4]
            sage: j = 3
            sage: j, y_pos = RuleSuperRSK().insertion(j, r, epsilon=1); r
            [1, 3, 3, 3, 4]
            sage: j
            3
            sage: y_pos
            1
        """
        bisect = bisect_right if epsilon == 0 else bisect_left

        if (r[-1] < j) or (r[-1] == j and epsilon == 0):
            return None, len(r)  # j needs to be added at the end of the list r.
        # Figure out where to insert j into the list r. The
        # bisect command returns the position of the least
        # element of r greater than j.  We will call it y.
        y_pos = bisect(r, j)
        # Switch j and y
        j, r[y_pos] = r[y_pos], j
        return j, y_pos

    def _forward_format_output(self, p, q, check_standard):
        r"""
        Return final output of the ``RSK`` (here, super RSK)
        correspondence from the output of the corresponding
        ``forward_rule``.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleSuperRSK
            sage: isinstance(RuleSuperRSK()._forward_format_output(
            ....:           [['1p', 1, '2p']], [['1p', '1', '2p']], True)[0],
            ....:           StandardSuperTableau)
            True
            sage: isinstance(RuleSuperRSK()._forward_format_output(
            ....:                   [[1, '2p', 3]], [[1, 2, 3]], False)[0],
            ....:                   SemistandardSuperTableau)
            True
            sage: isinstance(RuleSuperRSK()._forward_format_output(
            ....:                       [[1, 1, 3]], [[1, 2, 3]], True)[0],
            ....:                       SemistandardSuperTableau)
            True
        """
        from sage.combinat.tableau import StandardTableau
        from sage.combinat.super_tableau import SemistandardSuperTableau, StandardSuperTableau

        if not p:
            return [StandardTableau([]), StandardTableau([])]
        if check_standard:
            try:
                P = StandardSuperTableau(p)
            except ValueError:
                P = SemistandardSuperTableau(p)
            try:
                Q = StandardSuperTableau(q)
            except ValueError:
                Q = SemistandardSuperTableau(q)
            return [P, Q]
        return [SemistandardSuperTableau(p), SemistandardSuperTableau(q)]

    def backward_rule(self, p, q, output='array'):
        r"""
        Return the restricted super biword obtained by applying reverse
        super RSK insertion to a pair of tableaux ``(p, q)``.

        INPUT:

        - ``p``, ``q`` -- two tableaux of the same shape

        - ``output`` -- (default: ``'array'``) if ``q`` is row-strict:

          - ``'array'`` -- as a two-line array (i.e. restricted super biword)

          and if ``q`` is standard, we can have the output:

          - ``'word'`` -- as a word

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleSuperRSK
            sage: t1 = SemistandardSuperTableau([['1p', '3p', '4p'], [2], [3]])
            sage: t2 = SemistandardSuperTableau([[1, 2, 4], [3], [5]])
            sage: RuleSuperRSK().backward_rule(t1, t2, 'array')
            [[1, 2, 3, 4, 5], [4', 3, 3', 2, 1']]
            sage: t1 = SemistandardSuperTableau([[1, 3], ['3p']])
            sage: t2 = SemistandardSuperTableau([[1, 2], [3]])
            sage: RuleSuperRSK().backward_rule(t1, t2, 'array')
            [[1, 2, 3], [1, 3, 3']]
        """
        p_copy = [list(row) for row in p]
        upper_row = []
        lower_row = []
        # upper_row and lower_row will be the upper and lower rows of the
        # generalized permutation we get as a result, but both reversed
        d = {}
        for row, Li in enumerate(q):
            for col, val in enumerate(Li):
                if val in d and col in d[val]:
                    d[val][col].append(row)
                elif val not in d:
                    d[val] = {col: [row]}
                else:
                    d[val][col] = [row]
        # d is now a double family such that for every integers k and j,
        # the value d[k][j] is the list of rows i such that the (i, j)-th
        # cell of q is filled with k.
        for value, iter_dict in sorted(d.items(), reverse=True, key=lambda x: x[0]):
            epsilon = 1 if value.is_primed() else 0
            if epsilon == 1:
                iter_copy = dict(iter_dict)
                iter_dict = {}
                for k, v in iter_copy.items():
                    for vi in v:
                        if vi in iter_dict:
                            iter_dict[vi].append(k)
                        else:
                            iter_dict[vi] = [k]
            for key in sorted(iter_dict, reverse=True):
                for rows in iter_dict[key]:
                    row_index, col_index = (rows, key) if epsilon == 0 else (key, rows)
                    x = p_copy[row_index].pop()  # Always the right-most entry
                    while True:
                        if value.is_primed() == x.is_primed():
                            # row bumping
                            row_index -= 1
                            if row_index < 0:
                                break
                            x, col_index = self.reverse_insertion(x, p_copy[row_index], epsilon=epsilon)
                        else:
                            # column bumping
                            col_index -= 1
                            if col_index < 0:
                                break
                            c = self._get_col(p_copy, col_index)
                            x, row_index = self.reverse_insertion(x, c, epsilon=epsilon)
                            self._set_col(p_copy, col_index, c)
                    upper_row.append(value)
                    lower_row.append(x)
        return self._backward_format_output(lower_row, upper_row, output, q.is_standard())

    def reverse_insertion(self, x, row, epsilon=0):
        r"""
        Reverse bump the row ``row`` of the current insertion tableau
        with the number ``x`` using dual RSK insertion or classical
        Schensted insertion depending on the value of `epsilon`.

        The row ``row`` is modified in place. The bumped-out entry
        is returned along with the bumped position.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleSuperRSK
            sage: from bisect import bisect_left, bisect_right
            sage: r = [1, 3, 3, 3, 4]
            sage: j = 2
            sage: j, y = RuleSuperRSK().reverse_insertion(j, r, epsilon=0); r
            [2, 3, 3, 3, 4]
            sage: j
            1
            sage: y
            0
            sage: r = [1, 3, 3, 3, 4]
            sage: j = 3
            sage: j, y = RuleSuperRSK().reverse_insertion(j, r, epsilon=0); r
            [3, 3, 3, 3, 4]
            sage: j
            1
            sage: y
            0
            sage: r = [1, 3, 3, 3, 4]
            sage: j = (3)
            sage: j, y = RuleSuperRSK().reverse_insertion(j, r, epsilon=1); r
            [1, 3, 3, 3, 4]
            sage: j
            3
            sage: y
            3
        """
        bisect = bisect_left if epsilon == 0 else bisect_right
        y_pos = bisect(row, x) - 1
        # switch x and y
        x, row[y_pos] = row[y_pos], x
        return x, y_pos

    def _backward_format_output(self, lower_row, upper_row, output,
                                q_is_standard):
        r"""
        Return the final output of the ``RSK_inverse`` correspondence
        from the output of the corresponding ``backward_rule``.

        .. NOTE::

            The default implementation of ``backward_rule`` lists
            bumped-out entries in the order in which the reverse
            bumping happens, which is *opposite* to the order of the
            final output.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleSuperRSK
            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: RuleSuperRSK()._backward_format_output([PrimedEntry('1p'),
            ....:       PrimedEntry(1), PrimedEntry('3p'), PrimedEntry(9)],
            ....:       [PrimedEntry(1), PrimedEntry('2p'), PrimedEntry('3p'),
            ....:       PrimedEntry(4)], 'array', False)
            [[4, 3', 2', 1], [9, 3', 1, 1']]
            sage: RuleSuperRSK()._backward_format_output([PrimedEntry(1),
            ....:       PrimedEntry('2p'), PrimedEntry('3p'), PrimedEntry(4)],
            ....:       [PrimedEntry('1p'), PrimedEntry(1), PrimedEntry('2p'),
            ....:       PrimedEntry(2)], 'word', True)
            word: 4,3',2',1
            sage: RuleSuperRSK()._backward_format_output([PrimedEntry(1),
            ....:       PrimedEntry(2), PrimedEntry(3), PrimedEntry(4)],
            ....:       [PrimedEntry('1p'), PrimedEntry(1), PrimedEntry('2p'),
            ....:       PrimedEntry(2)], 'word', True)
            word: 4321
        """
        if output == 'array':
            return [list(reversed(upper_row)), list(reversed(lower_row))]
        if output == 'word':
            if q_is_standard:
                from sage.combinat.words.word import Word
                return Word(reversed(lower_row))
            else:
                raise TypeError("q must be standard to have a %s as "
                                "valid output" % output)
        raise ValueError("invalid output option")


class RuleStar(Rule):
    r"""
    Rule for `\star`-insertion.

    The `\star`-insertion is similar to the classical RSK algorithm
    and is defined in [MPPS2020]_. The bottom row of the increasing
    Hecke biword is a word in the 0-Hecke monoid that is fully
    commutative. When inserting a letter `x` into a row `R`, there
    are three cases:

    - Case 1: If `R` is empty or `x > \max(R)`, append `x` to row `R`
      and terminate.

    - Case 2: Otherwise if `x` is not in `R`, locate the smallest `y` in `R`
      with `y > x`. Bump `y` with `x` and insert `y` into the next row.

    - Case 3: Otherwise, if `x` is in `R`, locate the smallest `y` in `R` with
      `y \leq x` and interval `[y,x]` contained in `R`. Row `R` remains
      unchanged and `y` is to be inserted into the next row.

    The `\star`-insertion returns a pair consisting a conjugate of a
    semistandard tableau and a semistandard tableau. It is a bijection from the
    collection of all increasing Hecke biwords whose bottom row is a fully
    commutative word to pairs (P, Q) of tableaux of the same shape such that
    P is conjugate semistandard, Q is semistandard and the row reading word of
    P is fully commutative [MPPS2020]_.

    EXAMPLES:

    As an example of `\star`-insertion, we reproduce Example 28 in [MPPS2020]_::

        sage: from sage.combinat.rsk import RuleStar
        sage: p,q = RuleStar().forward_rule([1,1,2,2,4,4], [1,3,2,4,2,4])
        sage: ascii_art(p, q)
          1  2  4  1  1  2
          1  4     2  4
          3        4
        sage: line1,line2 = RuleStar().backward_rule(p, q)
        sage: line1,line2
        ([1, 1, 2, 2, 4, 4], [1, 3, 2, 4, 2, 4])
        sage: RSK_inverse(p, q, output='DecreasingHeckeFactorization', insertion='Star')
        (4, 2)()(4, 2)(3, 1)

        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
        sage: h = DecreasingHeckeFactorization([[4, 2], [], [4, 2], [3, 1]])
        sage: RSK_inverse(*RSK(h,insertion='Star'),insertion='Star',
        ....:             output='DecreasingHeckeFactorization')
        (4, 2)()(4, 2)(3, 1)
        sage: p,q = RSK(h, insertion='Star')
        sage: ascii_art(p, q)
        1  2  4  1  1  2
        1  4     2  4
        3        4
        sage: RSK_inverse(p, q, insertion='Star')
        [[1, 1, 2, 2, 4, 4], [1, 3, 2, 4, 2, 4]]
        sage: f = RSK_inverse(p, q, output='DecreasingHeckeFactorization', insertion='Star')
        sage: f == h
        True

    .. WARNING::

        When ``output`` is set to ``'DecreasingHeckeFactorization'``, the
        inverse of `\star`-insertion of `(P,Q)` returns a decreasing
        factorization whose number of factors is the maximum entry of `Q`::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            sage: h1 = DecreasingHeckeFactorization([[],[3,1],[1]]); h1
            ()(3, 1)(1)
            sage: P,Q = RSK(h1, insertion='Star')
            sage: ascii_art(P, Q)
              1  3  1  2
              1     2
            sage: h2 = RSK_inverse(P, Q, insertion='Star',
            ....: output='DecreasingHeckeFactorization'); h2
            (3, 1)(1)

    TESTS:

    Check that :func:`RSK` is the inverse of :func:`RSK_inverse` for various
    outputs/inputs::

        sage: from sage.combinat.partition import Partitions_n
        sage: shapes = [shape for n in range(7) for shape in Partitions_n(n)]
        sage: row_reading = lambda T: [x for row in reversed(T) for x in row]
        sage: from sage.monoids.hecke_monoid import HeckeMonoid
        sage: H = HeckeMonoid(SymmetricGroup(4+1))
        sage: from sage.combinat import permutation
        sage: reduce = lambda w: permutation.from_reduced_word(H.from_reduced_word(w).reduced_word())
        sage: fc = lambda w: not reduce(w).has_pattern([3,2,1])
        sage: FC_tabs = [T for shape in shapes
        ....:                  for T in SemistandardTableaux(shape, max_entry=4)
        ....:                      if fc(row_reading(T.conjugate()))]
        sage: Checks = []
        sage: for T in FC_tabs:  # long time
        ....:    shape = T.shape().conjugate()
        ....:    P = T.conjugate()
        ....:    Checks += [all((P,Q) == tuple(RSK(*RSK_inverse(P, Q,
        ....:               insertion='Star', output='array'),
        ....:               insertion='Star'))
        ....:               for Q in SemistandardTableaux(shape, max_entry=5))]
        sage: all(Checks)
        True
        sage: Checks = []
        sage: for T in FC_tabs:  # long time
        ....:    shape = T.shape().conjugate()
        ....:    P = T.conjugate()
        ....:    Checks += [all((P,Q) == tuple(RSK(RSK_inverse(P, Q,
        ....:               insertion='Star', output='DecreasingHeckeFactorization'),
        ....:               insertion='Star'))
        ....:               for Q in SemistandardTableaux(shape, max_entry=5))]
        sage: all(Checks)
        True
        sage: Checks = []
        sage: for T in FC_tabs:
        ....:    shape = T.shape().conjugate()
        ....:    P = T.conjugate()
        ....:    for Q in StandardTableaux(shape, max_entry=5):
        ....:        Checks += [(P,Q) == tuple(RSK(RSK_inverse(P, Q,
        ....:                   insertion='Star', output='word'),
        ....:                   insertion='Star'))]
        sage: all(Checks)
        True

    Check that :func:`RSK_inverse` is the inverse of :func:`RSK` on arrays
    and words::

        sage: S = SymmetricGroup(3+1)
        sage: from sage.combinat import permutation
        sage: FC = [x
        ....:       for x in S
        ....:           if (not permutation.from_reduced_word(
        ....:           x.reduced_word()).has_pattern([3,2,1]) and
        ....:           x.reduced_word())]
        sage: Triples = [(w, factors, ex)
        ....:            for w in FC
        ....:                for factors in range(2, 5+1)
        ....:                    for ex in range(4)]
        sage: Checks = []
        sage: for t in Triples:
        ....:     B = crystals.FullyCommutativeStableGrothendieck(*t)
        ....:     Checks += [all(b.to_increasing_hecke_biword() ==
        ....:                RSK_inverse(*RSK(
        ....:                *b.to_increasing_hecke_biword(),
        ....:                insertion='Star'), insertion='Star')
        ....:                for b in B)]
        sage: all(Checks)
        True

        sage: from sage.monoids.hecke_monoid import HeckeMonoid
        sage: Checks = []
        sage: H = HeckeMonoid(SymmetricGroup(3+1))
        sage: reduce = lambda w: permutation.from_reduced_word(H.from_reduced_word(w).reduced_word())
        sage: fc = lambda w: not reduce(w).has_pattern([3,2,1])
        sage: words = [w for n in range(10) for w in Words(3, n) if fc(w)]
        sage: all([all(w == RSK_inverse(*RSK(w, insertion='Star'),
        ....:          insertion='Star', output='word') for w in words)])
        True
    """

    def forward_rule(self, obj1, obj2=None, check_braid=True):
        r"""
        Return a pair of tableaux obtained by applying forward insertion
        to the increasing Hecke biword ``[obj1, obj2]``.

        INPUT:

        - ``obj1``, ``obj2`` -- can be one of the following ways to represent a
          biword (or, equivalently, an increasing 0-Hecke factorization) that
          is fully commutative:

          - two lists ``obj1`` and ``obj2`` of equal length, to be
            interpreted as the top row and the bottom row of the biword.

          - a word ``obj1`` in an ordered alphabet, to be interpreted as
            the bottom row of the biword (in this case, ``obj2`` is ``None``;
            the top row of the biword is understood to be `(1,2,\ldots,n)`
            by default).

          - a DecreasingHeckeFactorization ``obj1``, the whose increasing
            Hecke biword will be interpreted as the bottom row; the top row is
            understood to be the indices of the factors for each letter in
            this biword.

        - ``check_braid`` -- boolean (default: ``True``); indicator to validate
          that input is associated to a fully commutative word in the 0-Hecke
          monoid, validation is performed if set to ``True``. Οtherwise, this
          validation is ignored.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleStar
            sage: p,q = RuleStar().forward_rule([1,1,2,3,3], [2,3,3,1,3]); p,q
            ([[1, 3], [2, 3], [2]], [[1, 1], [2, 3], [3]])
            sage: p,q = RuleStar().forward_rule([2,3,3,1,3]); p,q
            ([[1, 3], [2, 3], [2]], [[1, 2], [3, 5], [4]])
            sage: p,q = RSK([1,1,2,3,3], [2,3,3,1,3], insertion=RSK.rules.Star); p,q
            ([[1, 3], [2, 3], [2]], [[1, 1], [2, 3], [3]])

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            sage: h = DecreasingHeckeFactorization([[3, 1], [3], [3, 2]])
            sage: p,q = RSK(h, insertion=RSK.rules.Star); p,q
            ([[1, 3], [2, 3], [2]], [[1, 1], [2, 3], [3]])

        TESTS:

        Empty objects::

            sage: from sage.combinat.rsk import RuleStar
            sage: p,q = RuleStar().forward_rule([]); p,q
            ([], [])

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            sage: h = DecreasingHeckeFactorization([[],[]])
            sage: p,q = RuleStar().forward_rule(h); p,q
            ([], [])

        Invalid inputs::

            sage: p,q = RuleStar().forward_rule([1,1,2,3,3], [2,2,3,1,3])
            Traceback (most recent call last):
            ...
            ValueError: [1, 1, 2, 3, 3], [2, 2, 3, 1, 3] is not an increasing factorization
            sage: p,q = RuleStar().forward_rule([1,1,2,2,4,4], [1,3,2,4,1,3])
            Traceback (most recent call last):
            ...
            ValueError: the Star insertion is not defined for non-fully commutative words
        """
        if obj2 is None and obj1 is not None:
            from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            if not isinstance(obj1, DecreasingHeckeFactorization):
                obj2 = obj1
                obj1 = list(range(1, len(obj1)+1))
            else:
                h = obj1
                obj1 = sum([[h.factors-i]*len(h.value[i]) for i in reversed(range(h.factors))], [])
                obj2 = [i for f in h.value[::-1] for i in reversed(f)]
        if len(obj1) != len(obj2):
            raise ValueError(f"{obj1} and {obj2} have different number of elements")
        for i in range(len(obj1)-1):
            if obj1[i] > obj1[i+1] or (obj1[i] == obj1[i+1] and obj2[i] >= obj2[i+1]):
                raise ValueError(f"{obj1}, {obj2} is not an increasing factorization")
        if check_braid:
            N = max(obj2)+1 if obj2 else 1
            from sage.monoids.hecke_monoid import HeckeMonoid
            from sage.groups.perm_gps.permgroup_named import SymmetricGroup
            H = HeckeMonoid(SymmetricGroup(N))
            h = H.from_reduced_word(obj2)
            from sage.combinat import permutation
            p = permutation.from_reduced_word(h.reduced_word())
            if p.has_pattern([3, 2, 1]):
                raise ValueError("the Star insertion is not defined for non-fully commutative words")

        p = []  # the "insertion" tableau
        q = []  # the "recording" tableau
        for i, j in zip(obj1, obj2):
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
        from sage.combinat.tableau import Tableau, SemistandardTableau
        p = Tableau(p)
        q = SemistandardTableau(q)
        return [p, q]

    def backward_rule(self, p, q, output='array'):
        r"""
        Return the increasing Hecke biword obtained by applying reverse
        `\star`-insertion to a pair of tableaux ``(p, q)``.

        INPUT:

        - ``p``, ``q`` -- two tableaux of the same shape, where ``p`` is the
          conjugate of a semistandard tableau, whose reading word is fully
          commutative and ``q`` is a semistandard tableau.

        - ``output`` -- (default: ``'array'``) if ``q`` is semi-standard:

          - ``'array'`` -- as a two-line array (i.e. generalized permutation
            or biword) that is an increasing Hecke biword
          - ``'DecreasingHeckeFactorization'`` -- as a decreasing
            factorization in the 0-Hecke monoid

          and if ``q`` is standard:

          - ``'word'`` -- as a (possibly non-reduced) word in the 0-Hecke
            monoid

        .. WARNING::

            When output is 'DecreasingHeckeFactorization', the number of factors
            in the output is the largest number in ``obj1``.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleStar
            sage: p,q = RuleStar().forward_rule([1,1,2,2,4,4], [1,3,2,4,2,4])
            sage: ascii_art(p, q)
              1  2  4  1  1  2
              1  4     2  4
              3        4
            sage: line1,line2 = RuleStar().backward_rule(p, q); line1,line2
            ([1, 1, 2, 2, 4, 4], [1, 3, 2, 4, 2, 4])
            sage: RuleStar().backward_rule(p, q, output = 'DecreasingHeckeFactorization')
            (4, 2)()(4, 2)(3, 1)

        TESTS:

        Empty objects::

            sage: RuleStar().backward_rule(Tableau([]), Tableau([]))
            [[], []]
            sage: RuleStar().backward_rule(Tableau([]), Tableau([]), output='word')
            word:
            sage: RuleStar().backward_rule(Tableau([]), Tableau([]),output='DecreasingHeckeFactorization')
            ()
        """
        from sage.combinat.tableau import SemistandardTableaux
        if p.shape() != q.shape():
            raise ValueError("p(=%s) and q(=%s) must have the same shape" % (p, q))
        if q not in SemistandardTableaux():
            raise ValueError("q(=%s) must be a semistandard tableau" % q)
        if p.conjugate() not in SemistandardTableaux():
            raise ValueError("the conjugate of p(=%s) must be a semistandard tableau" % p.conjugate())

        row_reading = [ele for row in reversed(p) for ele in row]
        N = max(row_reading + [0]) + 1
        from sage.monoids.hecke_monoid import HeckeMonoid
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        H = HeckeMonoid(SymmetricGroup(N))
        h = H.from_reduced_word(row_reading)
        from sage.combinat import permutation
        w = permutation.from_reduced_word(h.reduced_word())
        if w.has_pattern([3, 2, 1]):
            raise ValueError(f"the row reading word of the insertion tableau {p} is not fully-commutative")

        p_copy = p.to_list()
        line1 = []
        line2 = []
        d = {}
        for i, row in enumerate(q):
            for j, val in enumerate(row):
                if val in d:
                    d[val][j] = i
                else:
                    d[val] = {j: i}
        # d is now a double family such that for all integers k and j,
        # d[k][j]=i means that (i, j)-th cell of q is filled with k.
        for value, row_dict in sorted(d.items(), key=lambda x: -x[0]):
            for j in sorted(row_dict, reverse=True):
                i = row_dict[j]
                x = p_copy[i].pop()  # Always the right-most entry
                while i > 0:
                    i -= 1
                    row = p_copy[i]
                    x = self.reverse_insertion(x, row)
                line2.append(x)
                line1.append(value)
        return self._backward_format_output(line1[::-1], line2[::-1], output)

    def insertion(self, b, r):
        r"""
        Insert the letter ``b`` from the second row of the biword into the row
        ``r`` using `\star`-insertion defined in [MPPS2020]_.

        The row `r` is modified in place if bumping occurs and `b` is not in
        row `r`. The bumped-out entry, if it exists, is returned.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleStar
            sage: RuleStar().insertion(3, [1,2,4,5])
            4
            sage: RuleStar().insertion(3, [1,2,3,5])
            1
            sage: RuleStar().insertion(6, [1,2,3,5]) is None
            True
        """
        if r[-1] < b:
            return None  # append b to the end of row r
        if b in r:
            k = b
            while k in r:
                k -= 1
            k += 1
        else:
            y_pos = bisect_right(r, b)
            k = r[y_pos]
            r[y_pos] = b
        return k

    def reverse_insertion(self, x, r):
        """
        Reverse bump the row ``r`` of the current insertion tableau ``p``
        with number ``x``, provided that ``r`` is the ``i``-th row of ``p``.

        The row ``r`` is modified in place. The bumped-out entry is returned.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleStar
            sage: RuleStar().reverse_insertion(4, [1,2,3,5])
            3
            sage: RuleStar().reverse_insertion(1, [1,2,3,5])
            3
            sage: RuleStar().reverse_insertion(5, [1,2,3,5])
            5
        """
        if x in r:
            y = x
            while y in r:
                y += 1
            y -= 1
        else:
            y_pos = bisect_left(r, x) - 1
            y = r[y_pos]
            r[y_pos] = x
        x = y
        return x

    def _backward_format_output(self, obj1, obj2, output):
        r"""
        Return the final output of the ``RSK_inverse`` correspondence from
        the output of the corresponding ``backward_rule``.

        EXAMPLES::

            sage: from sage.combinat.rsk import RuleStar
            sage: RuleStar()._backward_format_output([1, 1, 2, 2, 4, 4], [1, 3, 2, 4, 2, 4], 'array')
            [[1, 1, 2, 2, 4, 4], [1, 3, 2, 4, 2, 4]]
            sage: RuleStar()._backward_format_output([1, 1, 2, 2, 4, 4], [1, 3, 2, 4, 2, 4], 'DecreasingHeckeFactorization')
            (4, 2)()(4, 2)(3, 1)
            sage: RuleStar()._backward_format_output([1, 2, 3, 4, 5, 6], [4, 2, 4, 2, 3, 1], 'word')
            word: 424231
        """
        if len(obj1) != len(obj2):
            raise ValueError(f"{obj1} and {obj2} are of different lengths")
        if output == 'array':
            return [obj1, obj2]
        elif output == 'word':
            if obj1 == list(range(1, len(obj1)+1)):
                from sage.combinat.words.word import Word
                return Word(obj2)
            else:
                raise TypeError("upper row must be standard")
        elif output == 'DecreasingHeckeFactorization':
            from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            obj1.reverse()
            obj2.reverse()
            df = []
            for j in range(len(obj1)):
                if j == 0:
                    df.append([])
                if j > 0 and obj1[j] < obj1[j-1]:
                    for _ in range(obj1[j-1]-obj1[j]):
                        df.append([])
                df[-1].append(obj2[j])
            if obj1:
                for a in range(obj1[-1]-1):
                    df.append([])
            # If biword is empty, return a decreasing factorization with 1 factor
            else:
                df.append([])
            return DecreasingHeckeFactorization(df)


class InsertionRules:
    r"""
    Catalog of rules for RSK-like insertion algorithms.
    """
    RSK = RuleRSK
    EG = RuleEG
    Hecke = RuleHecke
    dualRSK = RuleDualRSK
    coRSK = RuleCoRSK
    superRSK = RuleSuperRSK
    Star = RuleStar

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

    - ``obj1``, ``obj2`` -- can be one of the following:

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
      - ``RSK.rules.Hecke`` -- (or ``'hecke'``) Hecke insertion (only
        guaranteed for generalized permutations whose top row is strictly
        increasing) (:class:`~sage.combinat.rsk.RuleHecke`)
      - ``RSK.rules.dualRSK`` (or ``'dualRSK'``) -- Dual RSK insertion
        (only for strict biwords) (:class:`~sage.combinat.rsk.RuleDualRSK`)
      - ``RSK.rules.coRSK`` (or ``'coRSK'``) -- CoRSK insertion (only
        for strict cobiwords) (:class:`~sage.combinat.rsk.RuleCoRSK`)
      - ``RSK.rules.superRSK`` (or ``'super'``) -- Super RSK insertion (only for
        restricted super biwords) (:class:`~sage.combinat.rsk.RuleSuperRSK`)
      - ``RSK.rules.Star`` (or ``'Star'``) -- `\star`-insertion (only for
        fully commutative words in the 0-Hecke monoid)
        (:class:`~sage.combinat.rsk.RuleStar`)

    - ``check_standard`` -- boolean (default: ``False``); check if either of
      the resulting tableaux is a standard tableau, and if so, typecast it
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
        elif insertion == 'superRSK':
            insertion = RSK.rules.superRSK
        elif insertion == 'Star':
            insertion = RSK.rules.Star
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

    if isinstance(obj1, Matrix):
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
      - ``RSK.rules.superRSK`` (or ``'super'``) -- Super RSK insertion (only for
        restricted super biwords) (:class:`~sage.combinat.rsk.RuleSuperRSK`)
      - ``RSK.rules.Star`` (or ``'Star'``) -- `\star`-insertion (only for
        fully commutative words in the 0-Hecke monoid)
        (:class:`~sage.combinat.rsk.RuleStar`)

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

    Check that :issue:`20430` is fixed::

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
        elif insertion == 'superRSK':
            insertion = RSK.rules.superRSK
        elif insertion == 'Star':
            insertion = RSK.rules.Star
        else:
            raise ValueError("invalid input")

    rule = insertion()
    if not isinstance(rule, Rule):
        raise TypeError("the insertion must be an instance of Rule")

    if p.shape() != q.shape():
        raise ValueError(f"p(={p}) and q(={q}) must have the same shape")

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
