r"""
`k`-regular Sequences

An introduction and formal definition of `k`-regular sequences can be
found, for example, on the :wikipedia:`k-regular_sequence` or in
[AS2003]_.


.. WARNING::

    As this code is experimental, warnings are thrown when a
    `k`-regular sequence space is created for the first time in a
    session (see :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: Seq2 = kRegularSequenceSpace(2, ZZ)
        doctest:...: FutureWarning: This class/method/function is
        marked as experimental. It, its functionality or its interface
        might change without a formal deprecation.
        See http://trac.sagemath.org/21202 for details.


Examples
========

Binary sum of digits
--------------------

The binary sum of digits `S(n)` of a nonnegative integer `n` satisfies
`S(2n) = S(n)` and `S(2n+1) = S(n) + 1`. We model this by the following::

    sage: Seq2 = kRegularSequenceSpace(2, ZZ)
    sage: S = Seq2((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [1, 1]])),
    ....:          left=vector([0, 1]), right=vector([1, 0]))
    sage: S
    2-regular sequence 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, ...
    sage: all(S[n] == sum(n.digits(2)) for n in srange(10))
    True

Number of odd entries in Pascal's triangle
------------------------------------------

Let us consider the number of odd entries in the first `n` rows
of Pascals's triangle::

    sage: @cached_function
    ....: def u(n):
    ....:     if n <= 1:
    ....:         return n
    ....:     return 2*u(floor(n/2)) + u(ceil(n/2))
    sage: tuple(u(n) for n in srange(10))
    (0, 1, 3, 5, 9, 11, 15, 19, 27, 29)

There is a `2`-recursive sequence describing the numbers above as well::

    sage: U = Seq2((Matrix([[3, 2], [0, 1]]), Matrix([[2, 0], [1, 3]])),
    ....:          left=vector([0, 1]), right=vector([1, 0])).transposed()
    sage: all(U[n] == u(n) for n in srange(30))
    True


Various
=======

.. SEEALSO::

    :mod:`recognizable series <sage.combinat.recognizable_series>`,
    :mod:`sage.rings.cfinite_sequence`,
    :mod:`sage.combinat.binary_recurrence_sequences`.

AUTHORS:

- Daniel Krenn (2016, 2021)
- Gabriel F. Lipnik (2021)

ACKNOWLEDGEMENT:

- Daniel Krenn is supported by the
  Austrian Science Fund (FWF): P 24644-N26.
- Gabriel F. Lipnik is supported by the
  Austrian Science Fund (FWF): W 1230.


Classes and Methods
===================
"""
#*****************************************************************************
#       Copyright (C) 2016 Daniel Krenn <dev@danielkrenn.at>
#                     2021 Gabriel F. Lipnik <dev@gabriellipnik.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .recognizable_series import RecognizableSeries
from .recognizable_series import RecognizableSeriesSpace
from sage.misc.cachefunc import cached_function, cached_method


class kRegularSequence(RecognizableSeries):
    def __init__(self, parent, mu, left=None, right=None):
        r"""
        A `k`-regular sequence.

        INPUT:

        - ``parent`` -- an instance of :class:`kRegularSequenceSpace`

        - ``mu`` -- a family of square matrices, all of which have the
          same dimension. The indices of this family are `0,...,k-1`.
          ``mu`` may be a list or tuple of cardinality `k`
          as well. See also
          :meth:`~sage.combinat.recognizable_series.RecognizableSeries.mu`.

        - ``left`` -- (default: ``None``) a vector.
          When evaluating the sequence, this vector is multiplied
          from the left to the matrix product. If ``None``, then this
          multiplication is skipped.

        - ``right`` -- (default: ``None``) a vector.
          When evaluating the sequence, this vector is multiplied
          from the right to the matrix product. If ``None``, then this
          multiplication is skipped.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: S = Seq2((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:          vector([0, 1]), vector([1, 0])).transposed(); S
            2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...

        We can access the coefficients of a sequence by
        ::

            sage: S[5]
            11

        or iterating over the first, say `10`, by
        ::

            sage: from itertools import islice
            sage: list(islice(S, 10))
            [0, 1, 3, 5, 9, 11, 15, 19, 27, 29]

        .. SEEALSO::

            :doc:`k-regular sequence <k_regular_sequence>`,
            :class:`kRegularSequenceSpace`.
        """
        super(kRegularSequence, self).__init__(
            parent=parent, mu=mu, left=left, right=right)

    def _repr_(self):
        r"""
        Return a representation string of this `k`-regular sequence.

        OUTPUT:

        A string

        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: s = Seq2((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:           vector([0, 1]), vector([1, 0])).transposed()
            sage: repr(s)  # indirect doctest
            '2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...'
        """
        from sage.misc.lazy_list import lazy_list_formatter
        return lazy_list_formatter(
            self,
            name='{}-regular sequence'.format(self.parent().k),
            opening_delimiter='', closing_delimiter='',
            preview=10)

    @cached_method
    def __getitem__(self, n, **kwds):
        r"""
        Return the `n`-th entry of this sequence.

        INPUT:

        - ``n`` -- a nonnegative integer

        OUTPUT:

        An element of the universe of the sequence

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: S = Seq2((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
            ....:          left=vector([0, 1]), right=vector([1, 0]))
            sage: S[7]
            3

        TESTS::

            sage: S[-1]
            Traceback (most recent call last):
            ...
            ValueError: value -1 of index is negative

        ::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: W = Seq2.indices()
            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: S = Seq2((M0, M1), [0, 1], [1, 1])
            sage: S._mu_of_word_(W(0.digits(2))) == M0
            True
            sage: S._mu_of_word_(W(1.digits(2))) == M1
            True
            sage: S._mu_of_word_(W(3.digits(2))) == M1^2
            True
        """
        return self.coefficient_of_word(self.parent()._n_to_index_(n), **kwds)

    def __iter__(self):
        r"""
        Return an iterator over the coefficients of this sequence.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: S = Seq2((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
            ....:          left=vector([0, 1]), right=vector([1, 0]))
            sage: from itertools import islice
            sage: tuple(islice(S, 10))
             (0, 1, 1, 2, 1, 2, 2, 3, 1, 2)

        TESTS::

            sage: it = iter(S)
            sage: iter(it) is it
            True
            sage: iter(S) is not it
            True
        """
        from itertools import count
        return iter(self[n] for n in count())


class kRegularSequenceSpace(RecognizableSeriesSpace):
    r"""
    The space of `k`-regular Sequences over the given ``coefficients``.

    INPUT:

    - ``k`` -- an integer at least `2` specifying the base

    - ``coefficient_ring`` -- a (semi-)ring.

    - ``category`` -- (default: ``None``) the category of this
      space

    EXAMPLES::

        sage: kRegularSequenceSpace(2, ZZ)
        Space of 2-regular sequences over Integer Ring
        sage: kRegularSequenceSpace(3, ZZ)
        Space of 3-regular sequences over Integer Ring

    .. SEEALSO::

        :doc:`k-regular sequence <k_regular_sequence>`,
        :class:`kRegularSequence`.
    """
    Element = kRegularSequence

    @classmethod
    def __normalize__(cls, k, coefficient_ring, **kwds):
        r"""
        Normalizes the input in order to ensure a unique
        representation.

        For more information see :class:`kRegularSequenceSpace`.

        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2.category()
            Category of modules over Integer Ring
            sage: Seq2.alphabet()
            {0, 1}
        """
        from sage.arith.srange import srange
        nargs = super(kRegularSequenceSpace, cls).__normalize__(
            coefficient_ring, alphabet=srange(k), **kwds)
        return (k,) + nargs

    def __init__(self, k, *args, **kwds):
        r"""
        See :class:`kRegularSequenceSpace` for details.

        INPUT:

        - ``k`` -- an integer at least `2` specifying the base

        Other input arguments are passed on to
        :meth:`~sage.combinat.recognizable_series.RecognizableSeriesSpace.__init__`.

        TESTS::

            sage: kRegularSequenceSpace(2, ZZ)
            Space of 2-regular sequences over Integer Ring
            sage: kRegularSequenceSpace(3, ZZ)
            Space of 3-regular sequences over Integer Ring

        .. SEEALSO::

            :doc:`k-regular sequence <k_regular_sequence>`,
            :class:`kRegularSequence`.
        """
        self.k = k
        super(kRegularSequenceSpace, self).__init__(*args, **kwds)

    def _repr_(self):
        r"""
        Return a representation string of this `k`-regular sequence space.

        OUTPUT:

        A string

        TESTS::

            sage: repr(kRegularSequenceSpace(2, ZZ))  # indirect doctest
            'Space of 2-regular sequences over Integer Ring'
        """
        return 'Space of {}-regular sequences over {}'.format(self.k, self.base())

    def _n_to_index_(self, n):
        r"""
        Convert `n` to an index usable by the underlying
        recognizable series.

        INPUT:

        - ``n`` -- a nonnegative integer

        OUTPUT:

        A word

        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2._n_to_index_(6)
            word: 011
            sage: Seq2._n_to_index_(-1)
            Traceback (most recent call last):
            ...
            ValueError: value -1 of index is negative
        """
        from sage.rings.integer_ring import ZZ
        n = ZZ(n)
        W = self.indices()
        try:
            return W(n.digits(self.k))
        except OverflowError:
            raise ValueError('value {} of index is negative'.format(n)) from None

    def from_recurrence(self, *args, **kwds):
        r"""
        Construct a `k`-regular sequence that fulfills the recurrence relations
        given in ``equations``.

        INPUT:

        Positional arguments:

        - ``equations`` -- A list of equations where the elements have
          either the form

          - `f(k^M n + r) = c_{r,l} f(k^m n + l) + c_{r,l + 1} f(k^m n
            + l + 1) + ... + c_{r,u} f(k^m n + u)` for some integers
            `0 \leq r < k^M`, `M > m \geq 0` and `l \leq u`, and some
            coefficients `c_{r,j}` from the (semi)ring ``coefficents``
            of the corresponding :class:`kRegularSequenceSpace`, valid
            for all integers `n \geq \text{offset}` for some integer
            `\text{offset} \geq \max(-l/k^m, 0)` (default: ``0``), and
            there is an equation of this form (with the same
            parameters `M` and `m`) for all `r`

          or the form

          - ``f(k) == t`` for some integer ``k`` and some ``t`` from the (semi)ring
            ``coefficients``.

          The recurrence relations above uniquely determine a `k`-regular sequence;
          see [HKL2021]_ for further information.

        - ``function`` -- symbolic function ``f`` occurring in the equations

        - ``var`` -- symbolic variable (``n`` in the above description of ``equations``)

        Keyword-only argument:

        - ``offset`` -- an integer (default: ``0``). See explanation for ``equations`` above.

        OUTPUT:

        A :class:`kRegularSequence`.

        EXAMPLES:

        Stern--Brocot Sequence::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: var('n')
            n
            sage: function('f')
            f
            sage: Seq2.from_recurrence([
            ....:     f(2*n) == f(n), f(2*n + 1) == f(n) + f(n + 1),
            ....:     f(0) == 0, f(1) == 1], f, n)
            2-regular sequence 0, 1, 1, 2, 1, 3, 2, 3, 1, 4, ...

        Number of Odd Entries in Pascal's Triangle::

            sage: Seq2.from_recurrence([
            ....:     f(2*n) == 3*f(n), f(2*n + 1) == 2*f(n) + f(n + 1),
            ....:     f(0) == 0, f(1) == 1], f, n)
            2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...

        Number of Unbordered Factors in the Thue--Morse Sequence::

            sage: Seq2.from_recurrence([
            ....:     f(8*n) == 2*f(4*n),
            ....:     f(8*n + 1) == f(4*n + 1),
            ....:     f(8*n + 2) == f(4*n + 1) + f(4*n + 3),
            ....:     f(8*n + 3) == -f(4*n + 1) + f(4*n + 2),
            ....:     f(8*n + 4) == 2*f(4*n + 2),
            ....:     f(8*n + 5) == f(4*n + 3),
            ....:     f(8*n + 6) == -f(4*n + 1) + f(4*n + 2) + f(4*n + 3),
            ....:     f(8*n + 7) == 2*f(4*n + 1) + f(4*n + 3),
            ....:     f(0) == 1, f(1) == 2, f(2) == 2, f(3) == 4, f(4) == 2,
            ....:     f(5) == 4, f(6) == 6, f(7) == 0, f(8) == 4, f(9) == 4,
            ....:     f(10) == 4, f(11) == 4, f(12) == 12, f(13) == 0, f(14) == 4,
            ....:     f(15) == 4, f(16) == 8, f(17) == 4, f(18) == 8, f(19) == 0,
            ....:     f(20) == 8, f(21) == 4, f(22) == 4, f(23) == 8], f, n, offset=3)
            2-regular sequence 1, 2, 2, 4, 2, 4, 6, 0, 4, 4, ...

        Number of Non-Zero Elements in the Generalized Pascal's Triangle (see [LRS2017]_)::

            sage: Seq2 = kRegularSequenceSpace(2, QQ)
            sage: Seq2.from_recurrence([
            ....:     f(4*n) == 5/3*f(2*n) - 1/3*f(2*n + 1),
            ....:     f(4*n + 1) == 4/3*f(2*n) + 1/3*f(2*n + 1),
            ....:     f(4*n + 2) == 1/3*f(2*n) + 4/3*f(2*n + 1),
            ....:     f(4*n + 3) == -1/3*f(2*n) + 5/3*f(2*n + 1),
            ....:     f(0) == 1, f(1) == 2], f, n)
            2-regular sequence 1, 2, 3, 3, 4, 5, 5, 4, 5, 7, ...

        TESTS::

            sage: Seq2.from_recurrence([
            ....:     f(4*n) == f(2*n),
            ....:     f(4*n + 1) == f(2*n),
            ....:     f(4*n + 2) == f(2*n),
            ....:     f(4*n + 3) == f(2*n + 1024),
            ....:     f(0) == 1, f(1) == 1], f, n, offset=2)
            Traceback (most recent call last):
            ...
            ValueError: Initial values for arguments in [2, ..., 2044] are missing.

        ::

            sage: S = Seq2.from_recurrence([
            ....:     f(4*n) == f(2*n),
            ....:     f(4*n + 1) == f(2*n),
            ....:     f(4*n + 2) == f(2*n),
            ....:     f(4*n + 3) == f(2*n + 16),
            ....:     f(0) == 1, f(1) == 1, f(2) == 2, f(3) == 3, f(4) == 4,
            ....:     f(5) == 5, f(6) == 6, f(7) == 7, f(16) == 4, f(18) == 4,
            ....:     f(20) == 4, f(22) == 4, f(24) == 6, f(26) == 6, f(28) == 6],
            ....:     f, n, offset=2)
            sage: all([S[4*i] == S[2*i] and
            ....:      S[4*i + 1] == S[2*i] and
            ....:      S[4*i + 2] == S[2*i] and
            ....:      S[4*i + 3] == S[2*i + 16] for i in srange(2, 100)])
            True

        ::

            sage: S = Seq2.from_recurrence([
            ....:     f(4*n) == f(2*n),
            ....:     f(4*n + 1) == f(2*n),
            ....:     f(4*n + 2) == f(2*n),
            ....:     f(4*n + 3) == f(2*n - 16),
            ....:     f(0) == 1, f(1) == 1, f(2) == 2, f(3) == 3, f(4) == 4,
            ....:     f(5) == 5, f(6) == 6, f(7) == 7, f(8) == 8, f(9) == 9,
            ....:     f(10) == 10, f(11) == 11, f(12) == 12, f(13) == 13,
            ....:     f(14) == 14, f(15) == 15, f(16) == 16, f(17) == 17,
            ....:     f(18) == 18, f(19) == 19, f(20) == 20, f(21) == 21,
            ....:     f(22) == 22, f(23) == 23, f(24) == 24, f(25) == 25,
            ....:     f(26) == 26, f(27) == 27, f(28) == 28, f(29) == 29,
            ....:     f(30) == 30, f(31) == 31], f, n, offset=8)
            sage: all([S[4*i] == S[2*i] and
            ....:      S[4*i + 1] == S[2*i] and
            ....:      S[4*i + 2] == S[2*i] and
            ....:      S[4*i + 3] == S[2*i - 16] for i in srange(8, 100)])
            True

        Same test with different variable and function names::

            sage: var('m')
            m
            sage: function('g')
            g
            sage: T = Seq2.from_recurrence([
            ....:     g(4*m) == g(2*m),
            ....:     g(4*m + 1) == g(2*m),
            ....:     g(4*m + 2) == g(2*m),
            ....:     g(4*m + 3) == g(2*m - 16),
            ....:     g(0) == 1, g(1) == 1, g(2) == 2, g(3) == 3, g(4) == 4,
            ....:     g(5) == 5, g(6) == 6, g(7) == 7, g(8) == 8, g(9) == 9,
            ....:     g(10) == 10, g(11) == 11, g(12) == 12, g(13) == 13,
            ....:     g(14) == 14, g(15) == 15, g(16) == 16, g(17) == 17,
            ....:     g(18) == 18, g(19) == 19, g(20) == 20, g(21) == 21,
            ....:     g(22) == 22, g(23) == 23, g(24) == 24, g(25) == 25,
            ....:     g(26) == 26, g(27) == 27, g(28) == 28, g(29) == 29,
            ....:     g(30) == 30, g(31) == 31], g, m, offset=8)
            sage: all([S[i] == T[i] for i in srange(1000)])
            True

        Zero-sequence with non-zero initial values::

            sage: Seq2.from_recurrence([
            ....:     f(2*n) == 0, f(2*n + 1) == 0,
            ....:     f(0) == 1, f(1) == 1, f(2) == 2, f(3) == 3], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Initial value for argument 0 does not match with the given recurrence relations.

        ::

            sage: Seq2.from_recurrence([
            ....:     f(2*n) == 0, f(2*n + 1) == 0,
            ....:     f(0) == 1, f(1) == 1, f(2) == 2, f(3) == 3], f, n, offset=2)
            2-regular sequence 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, ...
        """
        RP = RecurrenceParser(self.k, self.coefficient_ring())
        mu, left, right = RP(*args, **kwds)
        return self(mu, left, right)


class RecurrenceParser(object):
    r"""
    A parser for symbolic recurrence relations that allow
    the construction of a `k`-linear representation
    for the sequence satisfying these recurrence relations.

    This is used by :meth:`kRegularSequenceSpace.from_recurrence`
    to construct a :class:`kRegularSequence`.
    """
    def __init__(self, k, coefficient_ring):
        r"""
        See :class:`RecurrenceParser`.

        INPUT:

        - ``k`` -- an integer at least `2` specifying the base

        - ``coefficient_ring`` -- a ring.

        These are the same parameters used when creating
        a :class:`kRegularSequenceSpace`.

        TESTS::

            sage: from sage.combinat.k_regular_sequence import RecurrenceParser
            sage: RecurrenceParser(2, ZZ)
            <sage.combinat.k_regular_sequence.RecurrenceParser object at 0x...>
        """
        self.k = k
        self.coefficient_ring = coefficient_ring

    def parse_recurrence(self, equations, function, var):
        r"""
        Parse recurrence relations as admissible in :meth:`kRegularSequenceSpace.from_recurrence`.

        INPUT:

        - ``equations`` -- see :meth:`kRegularSequenceSpace.from_recurrence`

        - ``function`` -- see :meth:`kRegularSequenceSpace.from_recurrence`

        - ``var`` -- see :meth:`kRegularSequenceSpace.from_recurrence`

        OUTPUT:

        A tuple consisting of

        - ``M``, ``m`` -- parameters of the recursive sequences,
          see [HKL2021]_, Definition 3.1

        - ``coeffs`` -- a dictionary mapping ``(r, j)`` to the coefficients
          `c_{r, j}` as given in [HKL2021]_, Equation (3.1)

        - ``initial_values`` -- a dictionary mapping integers ``n`` to the
          ``n``-th value of the sequence

        EXAMPLES::

            sage: from sage.combinat.k_regular_sequence import RecurrenceParser
            sage: RP = RecurrenceParser(2, ZZ)
            sage: var('n')
            n
            sage: function('f')
            f
            sage: RP.parse_recurrence([
            ....:     f(4*n) == f(2*n) + 2*f(2*n + 1) + 3*f(2*n - 2),
            ....:     f(4*n + 1) == 4*f(2*n) + 5*f(2*n + 1) + 6*f(2*n - 2),
            ....:     f(4*n + 2) == 7*f(2*n) + 8*f(2*n + 1) + 9*f(2*n - 2),
            ....:     f(4*n + 3) == 10*f(2*n) + 11*f(2*n + 1) + 12*f(2*n - 2),
            ....:     f(0) == 1, f(1) == 2, f(2) == 1], f, n)
            (2, 1, {(0, -2): 3, (0, 0): 1, (0, 1): 2, (1, -2): 6, (1, 0): 4,
            (1, 1): 5, (2, -2): 9, (2, 0): 7, (2, 1): 8, (3, -2): 12, (3, 0): 10,
            (3, 1): 11}, {0: 1, 1: 2, 2: 1})

        Stern--Brocot Sequence::

            sage: RP.parse_recurrence([
            ....:    f(2*n) == f(n), f(2*n + 1) == f(n) + f(n + 1),
            ....:    f(0) == 0, f(1) == 1], f, n)
            (1, 0, {(0, 0): 1, (1, 0): 1, (1, 1): 1}, {0: 0, 1: 1})

        .. SEEALSO::

            :meth:`kRegularSequenceSpace.from_recurrence`

        TESTS:

        The following tests check that the equations are well-formed::

            sage: RP.parse_recurrence([], f, n)
            Traceback (most recent call last):
            ...
            ValueError: List of recurrence equations is empty.

        ::

            sage: RP.parse_recurrence([f(4*n + 1)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: f(4*n + 1) is not an equation with ==.

        ::

            sage: RP.parse_recurrence([42], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 42 is not a symbolic expression.

        ::

            sage: RP.parse_recurrence([f(2*n) + 1 == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(2*n) + 1 in the equation f(2*n) + 1 == f(n) is
            not an evaluation of f.

        ::

            sage: RP.parse_recurrence([f(2*n, 5) == 3], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(2*n, 5) in the equation f(2*n, 5) == 3 does not
            have one argument.

        ::

            sage: RP.parse_recurrence([f() == 3], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f() in the equation f() == 3 does not have one
            argument.

        ::

            sage: RP.parse_recurrence([f(1/n + 1) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(1/n + 1) in the equation f(1/n + 1) == f(n):
            1/n + 1 is not a polynomial in n with integer coefficients.

        ::

            sage: RP.parse_recurrence([f(2*n + 1/2) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(2*n + 1/2) in the equation f(2*n + 1/2) == f(n):
            2*n + 1/2 is not a polynomial in n with integer coefficients.

        ::

            sage: RP.parse_recurrence([f(4*n^2) == f(2*n^2)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(4*n^2) in the equation f(4*n^2) == f(2*n^2):
            4*n^2 is not a polynomial in n of degree smaller than 2.

        ::

            sage: RP.parse_recurrence([f(42) == 1/2], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Initial value 1/2 given by the equation f(42) == (1/2)
            is not in Integer Ring.

        ::

            sage: RP.parse_recurrence([f(42) == 0, f(42) == 1], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Initial value f(42) is given twice.

        ::

            sage: RP.parse_recurrence([f(42) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Initial value f(n) given by the equation f(42) == f(n)
            is not in Integer Ring.

        ::

            sage: RP.parse_recurrence([f(4*n) == f(n), f(2*n) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(2*n) in the equation f(2*n) == f(n): 2 does not
            equal 4. Expected subsequence modulo 4 as in another equation, got
            subsequence modulo 2.

        ::

            sage: RP.parse_recurrence([f(3*n + 1) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(3*n + 1) in the equation f(3*n + 1) == f(n):
            3 is not a power of 2.

        ::

            sage: RP.parse_recurrence([f(n + 1) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(n + 1) in the equation f(n + 1) == f(n):
            1 is less than 2. Modulus must be at least 2.

        ::

            sage: RP.parse_recurrence([f(2*n) == f(n), f(2*n) == 0], f, n)
            Traceback (most recent call last):
            ...
            ValueError: There are more than one recurrence relation for f(2*n).

        ::

            sage: RP.parse_recurrence([f(2*n + 2) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(2*n + 2) in the equation f(2*n + 2) == f(n):
            remainder 2 is not smaller than modulus 2.

        ::

            sage: RP.parse_recurrence([f(2*n - 1) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(2*n - 1) in the equation f(2*n - 1) == f(n):
            remainder -1 is smaller than 0.

        ::

            sage: RP.parse_recurrence([f(2*n) == 2*n], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term 2*n in the equation f(2*n) == 2*n does not
            contain f.

        ::

            sage: RP.parse_recurrence([f(2*n) == 1/2*f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term 1/2*f(n) in the equation f(2*n) == 1/2*f(n):
            1/2 is not a valid coefficient since it is not in Integer Ring.

        ::

            sage: RP.parse_recurrence([f(2*n) == 1/f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 1/f(n) is not a valid right hand side.

        ::

            sage: RP.parse_recurrence([f(2*n) == 2*n*f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 2*n*f(n) is not a valid right hand side.

        ::

            sage: RP.parse_recurrence([f(2*n) == 2*f(n, 5)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(n, 5) in the equation f(2*n) == 2*f(n, 5)
            has more than one argument.

        ::

            sage: RP.parse_recurrence([f(2*n) == 2*f()], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f() in the equation f(2*n) == 2*f() has no argument.

        ::

            sage: RP.parse_recurrence([f(2*n) == 1/f(n) + 2*f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term 1/f(n) in the equation f(2*n) == 1/f(n) + 2*f(n)
            is not a valid summand.

        ::

            sage: RP.parse_recurrence([f(2*n) == 2*f(1/n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(1/n) in the equation f(2*n) == 2*f(1/n):
            1/n is not a polynomial in n with integer coefficients.

        ::

            sage: RP.parse_recurrence([f(2*n) == f(n + 1/2)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(n + 1/2) in the equation f(2*n) == f(n + 1/2):
            n + 1/2 is not a polynomial in n with integer coefficients.

        ::

            sage: RP.parse_recurrence([f(2*n) == f(1/2*n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(1/2*n) in the equation f(2*n) == f(1/2*n):
            1/2*n is not a polynomial in n with integer coefficients.

        ::

            sage: RP.parse_recurrence([f(2*n) == f(n^2 + 1)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(n^2 + 1) in the equation f(2*n) == f(n^2 + 1):
            polynomial n^2 + 1 does not have degree 1.

        ::

            sage: RP.parse_recurrence([f(2*n) == f(1)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(1) in the equation f(2*n) == f(1):
            polynomial 1 does not have degree 1.

        ::

            sage: RP.parse_recurrence([f(4*n) == f(2*n) + f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(n) in the equation f(4*n) == f(2*n) + f(n):
            1 does not equal 2. Expected subsequence modulo 2 as in another
            summand or equation, got subsequence modulo 1.

        ::

            sage: RP.parse_recurrence([f(4*n) == f(2*n), f(4*n + 1) == f(n)],
            ....:     f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(n) in the equation f(4*n + 1) == f(n): 1 does not
            equal 2. Expected subsequence modulo 2 as in another summand or
            equation, got subsequence modulo 1.

        ::

            sage: RP.parse_recurrence([f(4*n) == f(3*n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(3*n) in the equation f(4*n) == f(3*n): 3 is not
            a power of 2.

        ::

            sage: RP.parse_recurrence([f(2*n) == f(4*n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(4*n) in the equation f(2*n) == f(4*n):
            4 is not smaller than 2.

        ::

            sage: RP.parse_recurrence([f(2*n) == f(2*n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(2*n) in the equation f(2*n) == f(2*n):
            2 is not smaller than 2.

        ::

            sage: RP.parse_recurrence([f(2*n) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Recurrence relations for [f(2*n + 1)] are missing.

        ::

            sage: RP.parse_recurrence([f(4*n) == f(n), f(4*n + 3) == 0], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Recurrence relations for [f(4*n + 1), f(4*n + 2)]
            are missing.

        ::

            sage: RP.parse_recurrence([f(42) == 0], f, n)
            Traceback (most recent call last):
            ...
            ValueError: No recurrence relations are given.

        ::

            sage: RP.parse_recurrence(
            ....:     [f(4*n + r) == f(n) for r in srange(4)], f, n)
            (2, 0, {(0, 0): 1, (1, 0): 1, (2, 0): 1, (3, 0): 1}, {})

        ::

            sage: RP.parse_recurrence(
            ....:     [f(8*n) == f(n)] +
            ....:     [f(8*n + r) == f(2*n) for r in srange(1,8)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Term f(2*n) in the equation f(8*n + 1) == f(2*n):
            2 does not equal 1. Expected subsequence modulo 1 as in another
            summand or equation, got subsequence modulo 2.

        Finally, also for the zero-sequence the output is as expected::

            sage: RP.parse_recurrence([f(2*n) == 0, f(2*n + 1) == 0], f, n)
            (1, 0, {}, {})
        """
        from sage.arith.srange import srange
        from sage.functions.log import log
        from sage.rings.integer_ring import ZZ
        from sage.symbolic.operators import add_vararg, mul_vararg, operator

        k = self.k
        coefficient_ring = self.coefficient_ring
        M = None
        m = None
        coeffs = {}
        initial_values = {}
        remainders = set()

        def parse_multiplication(op, eq):
            operands = op.operands()
            assert op.operator() == mul_vararg and len(operands) == 2
            if operands[1].operator() == function:
                return [operands[0], operands[1]]
            elif operands[0].operator() == function:
                return [operands[1], operands[0]]
            else:
                raise ValueError('Term %s in the equation %s '
                                 'does not contain %s.'
                                 % (op, eq, function)) from None

        def parse_one_summand(summand, eq):
            if summand.operator() == mul_vararg:
                coeff, op = parse_multiplication(summand, eq)
            elif summand.operator() == function:
                coeff, op = 1, summand
            else:
                raise ValueError('Term %s in the equation %s is not a valid summand.'
                                 % (summand, eq)) from None
            if len(op.operands()) > 1:
                raise ValueError('Term %s in the equation %s has more than one argument.'
                                 % (op, eq)) from None
            elif len(op.operands()) == 0:
                raise ValueError('Term %s in the equation %s has no argument.'
                                 % (op, eq)) from None
            try:
                poly = ZZ[var](op.operands()[0])
            except TypeError:
                raise ValueError('Term %s in the equation %s: '
                                 '%s is not a polynomial in %s with integer coefficients.'
                                 % (op, eq, op.operands()[0], var)) from None
            if poly.degree() != 1:
                raise ValueError("Term %s in the equation %s: "
                                 "polynomial %s does not have degree 1."
                                 % (op, eq, poly)) from None
            d, base_power_m = list(poly)
            m = log(base_power_m, base=k)
            return [coeff, m, d]

        if not equations:
            raise ValueError("List of recurrence equations is empty.") from None

        for eq in equations:
            try:
                if eq.operator() != operator.eq:
                    raise ValueError("%s is not an equation with ==."
                                     % eq) from None
            except AttributeError:
                raise ValueError("%s is not a symbolic expression."
                                 % eq) from None
            left_side, right_side = eq.operands()
            if left_side.operator() != function:
                raise ValueError("Term %s in the equation %s is not an evaluation of %s."
                                 % (left_side, eq, function)) from None
            if  len(left_side.operands()) != 1:
                raise ValueError("Term %s in the equation %s does not have "
                                 "one argument."
                                 % (left_side, eq)) from None
            try:
                polynomial_left = ZZ[var](left_side.operands()[0])
            except TypeError:
                raise ValueError("Term %s in the equation %s: "
                                 "%s is not a polynomial in %s with "
                                 "integer coefficients."
                                 % (left_side, eq,
                                    left_side.operands()[0], var)) from None
            if polynomial_left.degree()  > 1:
                raise ValueError("Term %s in the equation %s: "
                                 "%s is not a polynomial in %s of degree smaller than 2."
                                 % (left_side, eq, polynomial_left, var)) from None
            if polynomial_left in ZZ:
                if right_side in coefficient_ring:
                    if (polynomial_left in initial_values.keys() and
                        initial_values[polynomial_left] != right_side):
                        raise ValueError("Initial value %s is given twice."
                                         % (function(polynomial_left))) from None
                    initial_values.update({polynomial_left: right_side})
                else:
                    raise ValueError("Initial value %s given by the equation %s "
                                     "is not in %s."
                                     % (right_side, eq, coefficient_ring)) from None
            else:
                [r, base_power_M] = list(polynomial_left)
                M_new = log(base_power_M, base=k)
                if M and M != M_new:
                    raise ValueError(("Term {0} in the equation {1}: "
                                      "{2} does not equal {3}. Expected "
                                      "subsequence modulo {3} as in another "
                                      "equation, got subsequence modulo {2}.").format(
                                          left_side, eq,
                                          base_power_M, k**M)) from None
                elif not M:
                    M = M_new
                    if M not in ZZ:
                        raise ValueError("Term %s in the equation %s: "
                                         "%s is not a power of %s."
                                         % (left_side, eq,
                                            base_power_M, k)) from None
                    if M < 1:
                        raise ValueError(("Term {0} in the equation {1}: "
                                          "{2} is less than {3}. Modulus must "
                                          "be at least {3}.").format(
                                              left_side, eq,
                                              base_power_M, k)) from None
                if r in remainders:
                    raise ValueError("There are more than one recurrence relation for %s."
                                     % (left_side,)) from None
                if r >= k**M:
                    raise ValueError("Term %s in the equation %s: "
                                     "remainder %s is not smaller than modulus %s."
                                     % (left_side, eq, r, k**M)) from None
                elif r < 0:
                    raise ValueError("Term %s in the equation %s: "
                                     "remainder %s is smaller than 0."
                                     % (left_side, eq, r)) from None
                else:
                    remainders.add(r)
                if right_side != 0:
                    if (len(right_side.operands()) == 1 and right_side.operator() == function
                        or right_side.operator() == mul_vararg and len(right_side.operands()) == 2):
                        summands = [right_side]
                    elif right_side.operator() == add_vararg:
                        summands = right_side.operands()
                    else:
                        raise ValueError("%s is not a valid right hand side."
                                         % (right_side,)) from None
                    for summand in summands:
                        coeff, new_m, d = parse_one_summand(summand, eq)
                        if coeff not in coefficient_ring:
                            raise ValueError("Term %s in the equation %s: "
                                             "%s is not a valid coefficient "
                                             "since it is not in %s."
                                             % (summand, eq, coeff, coefficient_ring)) from None
                        if m is not None and m != new_m:
                            raise ValueError(("Term {0} in the equation {1}: "
                                              "{2} does not equal {3}. Expected "
                                              "subsequence modulo {3} as in another "
                                              "summand or equation, got subsequence "
                                              "modulo {2}.").format(
                                                  summand, eq,
                                                  k**new_m, k**m)) from None
                        elif m is None:
                            m = new_m
                            if m not in ZZ:
                                raise ValueError("Term %s in the equation %s: "
                                                 "%s is not a power of %s."
                                                 % (summand, eq,
                                                    k**m, k)) from None
                            if M <= m:
                                raise ValueError("Term %s in the equation %s: "
                                                 "%s is not smaller than %s."
                                                 % (summand, eq,
                                                    k**m, k**M)) from None
                        coeffs.update({(r, d): coeff})

        if not M:
            raise ValueError("No recurrence relations are given.") from None
        elif M and m is None: # for the zero sequence
            m = M - 1

        missing_remainders = [rem for rem in srange(k**M)
                              if rem not in remainders]
        if missing_remainders:
            raise ValueError("Recurrence relations for %s are missing."
                             % ([function(k**M*var + rem)
                                 for rem in missing_remainders],))

        return (M, m, coeffs, initial_values)

    def parameters(self, M, m, coeffs, initial_values, offset):
        r"""
        Determine parameters from recurrence relations as admissible in
        :meth:`kRegularSequenceSpace.from_recurrence`.

        INPUT:

        - ``M``, ``m``, ``offset`` -- parameters of the recursive sequences,
          see [HKL2021]_, Definition 3.1, as well as :meth:`kRegularSequenceSpace.from_recurrence`

        - ``coeffs`` -- a dictionary where ``coeffs[(r, j)]`` is the
          coefficient `c_{r,j}` as given in :meth:`kRegularSequenceSpace.from_recurrence`.
          If ``coeffs[(r, j)]`` is not given for some ``r`` and ``j``,
          then it is assumed to be zero.

        - ``initial_values`` -- a dictionary mapping integers ``n`` to the
          ``n``-th value of the sequence

        OUTPUT:

        A namedtuple ``recurrence_rules`` consisting of

        - ``M``, ``m``, ``l``, ``u``, ``offset`` -- parameters of the recursive
          sequences, see [HKL2021]_, Definition 3.1

        - ``ll``, ``uu``, ``n1``, ``dim`` -- parameters and dimension of the
          resulting linear representation, see [HKL2021]_, Theorem A

        - ``coeffs`` -- a dictionary mapping ``(r, j)`` to the coefficients
          `c_{r, j}` as given in [HKL2021]_, Equation (3.1).
          If ``coeffs[(r, j)]`` is not given for some ``r`` and ``j``,
          then it is assumed to be zero.

        - ``initial_values`` -- a dictionary mapping integers ``n`` to the
          ``n``-th value of the sequence

        EXAMPLES::

            sage: from sage.combinat.k_regular_sequence import RecurrenceParser
            sage: RP = RecurrenceParser(2, ZZ)
            sage: RP.parameters(2, 1,
            ....: {(0, -2): 3, (0, 0): 1, (0, 1): 2, (1, -2): 6, (1, 0): 4,
            ....: (1, 1): 5, (2, -2): 9, (2, 0): 7, (2, 1): 8, (3, -2): 12,
            ....: (3, 0): 10, (3, 1): 11}, {0: 1, 1: 2, 2: 1, 3: 4}, 0)
            recurrence_rules(M=2, m=1, l=-2, u=1, ll=-6, uu=3, dim=14,
            coeffs={(0, -2): 3, (0, 0): 1, (0, 1): 2, (1, -2): 6, (1, 0): 4,
            (1, 1): 5, (2, -2): 9, (2, 0): 7, (2, 1): 8, (3, -2): 12,
            (3, 0): 10, (3, 1): 11}, initial_values={0: 1, 1: 2, 2: 1, 3: 4,
            4: 12, 5: 30, 6: 48, 7: 66, 8: 75, 9: 204, 10: 333, 11: 462,
            12: 216, 13: 594, -6: 0, -5: 0, -4: 0, -3: 0, -2: 0, -1: 0},
            offset=1, n1=3)

        .. SEEALSO::

            :meth:`kRegularSequenceSpace.from_recurrence`

        TESTS::

            sage: RP.parameters(1, 0, {(0, 0): 1}, {}, 0)
            Traceback (most recent call last):
            ...
            ValueError: No initial values are given.

        ::

            sage: RP.parameters(1, 0,
            ....: {(0, 0): 1, (1, 0): 1, (1, 1): 1}, {0: 1/2, 1: 2*i}, 0)
            Traceback (most recent call last):
            ...
            ValueError: Initial values for arguments in [0, 1] are not in Integer Ring.

        ::

            sage: RP.parameters(1, 0, {(0, 0): 1},
            ....: {0: 1, 1: 0}, 0)
            recurrence_rules(M=1, m=0, l=0, u=0, ll=0, uu=0, dim=1,
            coeffs={(0, 0): 1}, initial_values={0: 1, 1: 0}, offset=0, n1=0)

        Finally, also for the zero-sequence the output is as expected::

            sage: RP.parameters(1, 0, {}, {0: 0}, 0)
            recurrence_rules(M=1, m=0, l=0, u=0, ll=0, uu=0, dim=1,
            coeffs={}, initial_values={0: 0}, offset=0, n1=0)

        ::

            sage: RP.parameters(1, 0,
            ....: {(0, 0): 0, (1, 1): 0}, {0: 0}, 0)
            recurrence_rules(M=1, m=0, l=0, u=0, ll=0, uu=0, dim=1,
            coeffs={(0, 0): 0, (1, 1): 0}, initial_values={0: 0},
            offset=0, n1=0)
        """
        from collections import namedtuple

        from sage.functions.other import ceil, floor

        coefficient_ring = self.coefficient_ring
        k = self.k
        keys_coeffs = coeffs.keys()
        indices_right = [key[1] for key in keys_coeffs if coeffs[key]]

        if not indices_right: # the sequence is the zero sequence
            l = 0
            u = 0
        else:
            l = min(indices_right)
            u = max(indices_right)

        if offset < max(0, -l/k**m):
            offset = max(0, ceil(-l/k**m))

        ll = (floor((l*k**(M-m) - k**M + 1)/(k**(M-m) - 1)) + 1)*(l < 0)
        uu = max([ceil((u*k**(M-m) + k**M - k**m)/(k**(M-m) - 1)) - 1, k**m - 1])
        n1 = offset - floor(ll/k**M)
        dim = (k**M - 1)/(k - 1) + (M - m)*(uu - ll - k**m + 1) + n1

        if not initial_values:
            raise ValueError("No initial values are given.") from None
        keys_initial = initial_values.keys()
        values_not_in_ring = [n for n in keys_initial
                              if initial_values[n] not in coefficient_ring]
        if values_not_in_ring:
            raise ValueError("Initial values for arguments in %s are not in %s."
                             % (values_not_in_ring, coefficient_ring)) from None

        last_value_needed = max(
            k**(M-1) - k**m + uu + (n1 > 0)*k**(M-1)*(k*(n1 - 1) + k - 1), # for matrix W
            k**m*offset + u,
            max(keys_initial))
        initial_values = self.values(
            M=M, m=m, l=l, u=u, ll=ll, coeffs=coeffs,
            initial_values=initial_values, last_value_needed=last_value_needed,
            offset=offset)

        recurrence_rules = namedtuple('recurrence_rules',
                                      ['M', 'm', 'l', 'u', 'll', 'uu', 'dim',
                                       'coeffs', 'initial_values', 'offset', 'n1'])

        return recurrence_rules(M=M, m=m, l=l, u=u, ll=ll, uu=uu, dim=dim,
                                coeffs=coeffs, initial_values=initial_values,
                                offset=offset, n1=n1)

    def values(self, *, M, m, l, u, ll, coeffs,
               initial_values, last_value_needed, offset):
        r"""
        Determine enough values of the corresponding recursive sequence by
        applying the recurrence relations given in :meth:`kRegularSequenceSpace.from_recurrence`
        to the values given in ``initial_values``.

        INPUT:

        - ``M``, ``m``, ``l``, ``u``, ``offset`` -- parameters of the
          recursive sequences, see [HKL2021]_, Definition 3.1

        - ``ll`` -- parameter of the resulting linear representation,
          see [HKL2021]_, Theorem A

        - ``coeffs`` -- a dictionary where ``coeffs[(r, j)]`` is the
          coefficient `c_{r,j}` as given in :meth:`kRegularSequenceSpace.from_recurrence`.
          If ``coeffs[(r, j)]`` is not given for some ``r`` and ``j``,
          then it is assumed to be zero.

        - ``initial_values`` -- a dictionary mapping integers ``n`` to the
          ``n``-th value of the sequence

        - ``last_value_needed`` -- last initial value which is needed to
          determine the linear representation

        OUTPUT:

        A dictionary mapping integers ``n`` to the ``n``-th value of the
        sequence for all ``n`` up to ``last_value_needed``.

        EXAMPLES:

        Stern--Brocot Sequence::

            sage: from sage.combinat.k_regular_sequence import RecurrenceParser
            sage: RP = RecurrenceParser(2, ZZ)
            sage: RP.values(M=1, m=0, l=0, u=1, ll=0,
            ....:     coeffs={(0, 0): 1, (1, 0): 1, (1, 1): 1},
            ....:     initial_values={0: 0, 1: 1, 2: 1}, last_value_needed=20,
            ....:     offset=0)
            {0: 0, 1: 1, 2: 1, 3: 2, 4: 1, 5: 3, 6: 2, 7: 3, 8: 1, 9: 4, 10: 3,
            11: 5, 12: 2, 13: 5, 14: 3, 15: 4, 16: 1, 17: 5, 18: 4, 19: 7, 20: 3}

        .. SEEALSO::

            :meth:`kRegularSequenceSpace.from_recurrence`

        TESTS:

        For the equations `f(2n) = f(n)` and `f(2n + 1) = f(n) + f(n + 1)`::

            sage: RP.values(M=1, m=0, l=0, u=1, ll=0,
            ....:     coeffs={(0, 0): 1, (1, 0): 1, (1, 1): 1},
            ....:     initial_values={0: 0, 1: 2}, last_value_needed=20,
            ....:     offset=0)
            {0: 0, 1: 2, 2: 2, 3: 4, 4: 2, 5: 6, 6: 4, 7: 6, 8: 2, 9: 8, 10: 6,
            11: 10, 12: 4, 13: 10, 14: 6, 15: 8, 16: 2, 17: 10, 18: 8, 19: 14,
            20: 6}

        ::

            sage: RP.values(M=1, m=0, l=0, u=1, ll=0,
            ....:     coeffs={(0, 0): 1, (1, 0): 1, (1, 1): 1},
            ....:     initial_values={}, last_value_needed=20, offset=0)
            Traceback (most recent call last):
            ...
            ValueError: Initial values for arguments in [0, 1] are missing.

        ::

            sage: RP.values(M=1, m=0, l=0, u=1, ll=0,
            ....:     coeffs={(0, 0): 1, (1, 0): 1, (1, 1): 1},
            ....:     initial_values={0: 0}, last_value_needed=20, offset=0)
            Traceback (most recent call last):
            ...
            ValueError: Initial values for arguments in [1] are missing.

        ::

            sage: RP.values(M=1, m=0, l=0, u=1, ll=0,
            ....:     coeffs={(0, 0): 1, (1, 0): 1, (1, 1): 1},
            ....:     initial_values={0: 0, 2: 1}, last_value_needed=20,
            ....:     offset=0)
            Traceback (most recent call last):
            ...
            ValueError: Initial values for arguments in [1] are missing.

        ::

            sage: RP.values(M=1, m=0, l=0, u=1, ll=0,
            ....:     coeffs={(0, 0): 1, (1, 0): 1, (1, 1): 1},
            ....:     initial_values={0: 0, 1: 2, 2:0}, last_value_needed=20,
            ....:     offset=0)
            Traceback (most recent call last):
            ...
            ValueError: Initial value for argument 2 does not match with the given
            recurrence relations.

        ::

            sage: RP.values(M=1, m=0, l=-2, u=2, ll=-2,
            ....:     coeffs={(0, -2): 1, (0, 2): 1, (1, -2): 1, (1, 2): 1},
            ....:     initial_values={0: 0, 1: 2, 2: 4, 3: 3, 4: 2},
            ....:     last_value_needed=20, offset=2)
            {-2: 0, -1: 0, 0: 0, 1: 2, 2: 4, 3: 3, 4: 2, 5: 2, 6: 4, 7: 4,
            8: 8, 9: 8, 10: 7, 11: 7, 12: 10, 13: 10, 14: 10, 15: 10, 16: 11,
            17: 11, 18: 11, 19: 11, 20: 18}

        Finally, also for the zero-sequence the output is as expected::

            sage: RP.values(M=1, m=0, l=0, u=0, ll=0,
            ....:     coeffs={}, initial_values={}, last_value_needed=10,
            ....:     offset=0)
            {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0}

        ::

            sage: RP.values(M=1, m=0, l=0, u=0, ll=0,
            ....:     coeffs={(0, 0): 0, (1, 1): 0}, initial_values={},
            ....:     last_value_needed=10, offset=0)
            {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0}
        """
        from sage.arith.srange import srange
        from sage.rings.integer_ring import ZZ

        k = self.k
        keys_initial = initial_values.keys()

        values = {n: None if n not in keys_initial else initial_values[n]
                  for n in srange(last_value_needed + 1)}
        missing_values = []

        @cached_function
        def coeff(r, k):
            try:
                return coeffs[(r, k)]
            except KeyError:
                return 0

        def f(n):
            f_n = values[n]
            if f_n is not None and f_n != "pending":
                return f_n
            elif f_n == "pending":
                missing_values.append(n)
                return 0
            else:
                values.update({n: "pending"})
                q, r = ZZ(n).quo_rem(k**M)
                if q < offset:
                    missing_values.append(n)
                return sum([coeff(r, j)*f(k**m*q + j)
                            for j in srange(l, u + 1)
                            if coeff(r, j)])

        for n in srange(last_value_needed + 1):
            values.update({n: f(n)})

        if missing_values:
            raise ValueError("Initial values for arguments in %s are missing."
                             % (list(set(missing_values)),)) from None

        for n in keys_initial:
            q, r = ZZ(n).quo_rem(k**M)
            if (q >= offset and
                values[n] != sum([coeff(r, j)*values[k**m*q + j]
                                  for j in srange(l, u + 1)])):
                raise ValueError("Initial value for argument %s does not match with "
                                 "the given recurrence relations."
                                 % (n,)) from None

        values.update({n: 0 for n in srange(ll, 0)})

        return values

    @cached_method
    def ind(self, M, m, ll, uu):
        r"""
        Determine the index operator corresponding to the recursive
        sequence as defined in [HKL2021]_.

        INPUT:

        - ``M``, ``m`` -- parameters of the recursive sequences,
          see [HKL2021]_, Definition 3.1

        - ``ll``, ``uu`` -- parameters of the resulting linear representation,
          see [HKL2021]_, Theorem A

        OUTPUT:

        A dictionary which maps both row numbers to subsequence parameters and
        vice versa, i.e.,

        - ``ind[i]`` -- a pair ``(j, d)`` representing the sequence `x(k^j n + d)`
          in the `i`-th component (0-based) of the resulting linear representation,

        - ``ind[(j, d)]`` -- the (0-based) row number of the sequence
          `x(k^j n + d)` in the linear representation.

        EXAMPLES::

            sage: from sage.combinat.k_regular_sequence import RecurrenceParser
            sage: RP = RecurrenceParser(2, ZZ)
            sage: RP.ind(3, 1, -3, 3)
            {(0, 0): 0, (1, -1): 3, (1, -2): 2, (1, -3): 1,
            (1, 0): 4, (1, 1): 5, (1, 2): 6, (1, 3): 7, (2, -1): 10,
            (2, -2): 9, (2, -3): 8, (2, 0): 11, (2, 1): 12, (2, 2): 13,
            (2, 3): 14, (2, 4): 15, (2, 5): 16, 0: (0, 0), 1: (1, -3),
            10: (2, -1), 11: (2, 0), 12: (2, 1), 13: (2, 2), 14: (2, 3),
            15: (2, 4), 16: (2, 5), 2: (1, -2), 3: (1, -1), 4: (1, 0),
            5: (1, 1), 6: (1, 2), 7: (1, 3), 8: (2, -3), 9: (2, -2)}

        .. SEEALSO::

            :meth:`kRegularSequenceSpace.from_recurrence`
        """
        from sage.arith.srange import srange

        k = self.k
        ind = {}

        pos = 0
        for j in srange(m):
            for d in srange(k**j):
                ind.update({(j, d): pos, pos: (j, d)})
                pos += 1
        for j in srange(m, M):
            for d in srange(ll, k**j - k**m + uu + 1):
                ind.update({(j, d): pos, pos: (j, d)})
                pos += 1

        return ind

    def v_eval_n(self, recurrence_rules, n):
        r"""
        Return the vector `v(n)` as given in [HKL2021]_, Theorem A.

        INPUT:

        - ``recurrence_rules`` -- a namedtuple generated by
          :meth:`parameters`

        - ``n`` -- an integer

        OUTPUT:

        A vector.

        EXAMPLES:

        Stern--Brocot Sequence::

            sage: from sage.combinat.k_regular_sequence import RecurrenceParser
            sage: RP = RecurrenceParser(2, ZZ)
            sage: SB_rules = RP.parameters(
            ....:     1, 0, {(0, 0): 1, (1, 0): 1, (1, 1): 1},
            ....:     {0: 0, 1: 1, 2: 1}, 0)
            sage: RP.v_eval_n(SB_rules, 0)
            (0, 1, 1)

        .. SEEALSO::

            :meth:`kRegularSequenceSpace.from_recurrence`
        """
        from sage.arith.srange import srange
        from sage.modules.free_module_element import vector

        k = self.k
        M = recurrence_rules.M
        m = recurrence_rules.m
        ll = recurrence_rules.ll
        uu = recurrence_rules.uu
        dim = recurrence_rules.dim - recurrence_rules.n1
        initial_values = recurrence_rules.initial_values
        ind = self.ind(M, m, ll, uu)

        return vector(
            [initial_values[k**ind[i][0]*n + ind[i][1]] for i in srange(dim)])

    def matrix(self, recurrence_rules, rem, correct_offset=True):
        r"""
        Construct the matrix for remainder ``rem`` of the linear
        representation of the sequence represented by ``recurrence_rules``.

        INPUT:

        - ``recurrence_rules`` -- a namedtuple generated by
          :meth:`parameters`

        - ``rem`` -- an integer between ``0`` and ``k - 1``

        - ``correct_offset`` -- (default: ``True``) a boolean. If
          ``True``, then the resulting linear representation has no
          offset.  See [HKL2021]_ for more information.

        OUTPUT:

        A matrix.

        EXAMPLES:

        The following example illustrates how the coefficients in the
        right-hand sides of the recurrence relations correspond to the entries of
        the matrices. ::

            sage: from sage.combinat.k_regular_sequence import RecurrenceParser
            sage: RP = RecurrenceParser(2, ZZ)
            sage: var('n')
            n
            sage: function('f')
            f
            sage: M, m, coeffs, initial_values = RP.parse_recurrence([
            ....:     f(8*n) == -1*f(2*n - 1) + 1*f(2*n + 1),
            ....:     f(8*n + 1) == -11*f(2*n - 1) + 10*f(2*n) + 11*f(2*n + 1),
            ....:     f(8*n + 2) == -21*f(2*n - 1) + 20*f(2*n) + 21*f(2*n + 1),
            ....:     f(8*n + 3) == -31*f(2*n - 1) + 30*f(2*n) + 31*f(2*n + 1),
            ....:     f(8*n + 4) == -41*f(2*n - 1) + 40*f(2*n) + 41*f(2*n + 1),
            ....:     f(8*n + 5) == -51*f(2*n - 1) + 50*f(2*n) + 51*f(2*n + 1),
            ....:     f(8*n + 6) == -61*f(2*n - 1) + 60*f(2*n) + 61*f(2*n + 1),
            ....:     f(8*n + 7) == -71*f(2*n - 1) + 70*f(2*n) + 71*f(2*n + 1),
            ....:     f(0) == 0, f(1) == 1, f(2) == 2, f(3) == 3, f(4) == 4,
            ....:     f(5) == 5, f(6) == 6, f(7) == 7], f, n)
            sage: rules = RP.parameters(
            ....:     M, m, coeffs, initial_values, 0)
            sage: RP.matrix(rules, 0, False)
            [  0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0]
            [  0 -51  50  51   0   0   0   0   0   0   0   0   0   0   0   0   0]
            [  0 -61  60  61   0   0   0   0   0   0   0   0   0   0   0   0   0]
            [  0 -71  70  71   0   0   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0  -1   0   1   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -11  10  11   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -21  20  21   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -31  30  31   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -41  40  41   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -51  50  51   0   0   0   0   0   0   0   0   0   0   0]
            sage: RP.matrix(rules, 1, False)
            [  0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1]
            [  0   0   0 -11  10  11   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -21  20  21   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -31  30  31   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -41  40  41   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -51  50  51   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -61  60  61   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -71  70  71   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0  -1   0   1   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0 -11  10  11   0   0   0   0   0   0   0   0   0]

        Stern--Brocot Sequence::

            sage: SB_rules = RP.parameters(
            ....:     1, 0, {(0, 0): 1, (1, 0): 1, (1, 1): 1},
            ....:     {0: 0, 1: 1, 2: 1}, 0)
            sage: RP.matrix(SB_rules, 0)
            [1 0 0]
            [1 1 0]
            [0 1 0]
            sage: RP.matrix(SB_rules, 1)
            [1 1 0]
            [0 1 0]
            [0 1 1]

        Number of Unbordered Factors in the Thue--Morse Sequence::

            sage: M, m, coeffs, initial_values = RP.parse_recurrence([
            ....:     f(8*n) == 2*f(4*n),
            ....:     f(8*n + 1) == f(4*n + 1),
            ....:     f(8*n + 2) == f(4*n + 1) + f(4*n + 3),
            ....:     f(8*n + 3) == -f(4*n + 1) + f(4*n + 2),
            ....:     f(8*n + 4) == 2*f(4*n + 2),
            ....:     f(8*n + 5) == f(4*n + 3),
            ....:     f(8*n + 6) == -f(4*n + 1) + f(4*n + 2) + f(4*n + 3),
            ....:     f(8*n + 7) == 2*f(4*n + 1) + f(4*n + 3),
            ....:     f(0) == 1, f(1) == 2, f(2) == 2, f(3) == 4, f(4) == 2,
            ....:     f(5) == 4, f(6) == 6, f(7) == 0, f(8) == 4, f(9) == 4,
            ....:     f(10) == 4, f(11) == 4, f(12) == 12, f(13) == 0, f(14) == 4,
            ....:     f(15) == 4, f(16) == 8, f(17) == 4, f(18) == 8, f(19) == 0,
            ....:     f(20) == 8, f(21) == 4, f(22) == 4, f(23) == 8], f, n)
            sage: UB_rules = RP.parameters(
            ....:     M, m, coeffs, initial_values, 3)
            sage: RP.matrix(UB_rules, 0)
            [ 0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  2  0  0  0  0  0  0  0  0  0 -1  0  0]
            [ 0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  1  0  1  0  0  0  0  0  0 -4  0  0]
            [ 0  0  0  0 -1  1  0  0  0  0  0  0  0  4  2  0]
            [ 0  0  0  0  0  2  0  0  0  0  0  0  0 -2  0  0]
            [ 0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0 -1  1  1  0  0  0  0  0  0  2  2  0]
            [ 0  0  0  0  2  0  1  0  0  0  0  0  0 -8 -4 -4]
            [ 0  0  0  0  0  0  0  2  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0]
            sage: RP.matrix(UB_rules, 1)
            [ 0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  2  0  0  0  0  0  0  0 -2  0  0]
            [ 0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0 -1  1  1  0  0  0  0  0  0  2  2  0]
            [ 0  0  0  0  2  0  1  0  0  0  0  0  0 -8 -4 -4]
            [ 0  0  0  0  0  0  0  2  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0 -1  1  0  0  0  2  0  0]
            [ 0  0  0  0  0  0  0  0  0  2  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]

        .. SEEALSO::

            :meth:`kRegularSequenceSpace.from_recurrence`
        """
        from sage.arith.srange import srange
        from sage.matrix.constructor import Matrix
        from sage.matrix.special import block_matrix, zero_matrix
        from sage.modules.free_module_element import vector

        coefficient_ring = self.coefficient_ring
        k = self.k
        M = recurrence_rules.M
        m = recurrence_rules.m
        l = recurrence_rules.l
        ll = recurrence_rules.ll
        uu = recurrence_rules.uu
        dim = recurrence_rules.dim
        n1 = recurrence_rules.n1
        dim_without_corr = dim - n1
        coeffs = recurrence_rules.coeffs
        ind = self.ind(M, m, ll, uu)

        @cached_function
        def coeff(r, k):
            try:
                return coeffs[(r, k)]
            except KeyError:
                return 0

        def entry(i, kk):
            j, d = ind[i]
            if j < M - 1:
                return int(kk == ind[(j + 1, k**j*rem + d)])
            else:
                rem_d = k**(M-1)*rem + (d%k**M)
                dd = d // k**M
                if rem_d < k**M:
                    lambd = l - ind[(m, (k**m)*dd + l)]
                    return coeff(rem_d, kk + lambd)
                else:
                    lambd = l - ind[(m, k**m*dd + k**m + l)]
                    return coeff(rem_d - k**M, kk + lambd)

        mat = Matrix(coefficient_ring, dim_without_corr, dim_without_corr, entry)

        if n1 == 0 or not correct_offset:
            return mat
        else:
            W = Matrix(coefficient_ring, dim_without_corr, 0)
            for i in srange(n1):
                W = W.augment(
                    self.v_eval_n(recurrence_rules, k*i + rem) -
                    mat*self.v_eval_n(recurrence_rules, i))

            J = Matrix(coefficient_ring, 0, n1)
            for i in srange(n1):
                J = J.stack(vector([int(j*k == i - rem) for j in srange(n1)]))

            Z = zero_matrix(coefficient_ring, n1, dim_without_corr)
            return block_matrix([[mat, W], [Z, J]], subdivide=False)

    def left(self, recurrence_rules):
        r"""
        Construct the vector ``left`` of the linear representation of
        recursive sequences.

        INPUT:

        - ``recurrence_rules`` -- a namedtuple generated by
          :meth:`parameters`; it only needs to contain a field
          ``dim`` (a positive integer)

        OUTPUT:

        A vector.

        EXAMPLES::

            sage: from collections import namedtuple
            sage: from sage.combinat.k_regular_sequence import RecurrenceParser
            sage: RP = RecurrenceParser(2, ZZ)
            sage: RRD = namedtuple('recurrence_rules_dim', ['dim'])
            sage: recurrence_rules = RRD(dim=5)
            sage: RP.left(recurrence_rules)
            (1, 0, 0, 0, 0)

        .. SEEALSO::

            :meth:`kRegularSequenceSpace.from_recurrence`
        """
        from sage.modules.free_module_element import vector
        dim = recurrence_rules.dim

        return vector([1] + (dim - 1)*[0])

    def right(self, recurrence_rules):
        r"""
        Construct the vector ``right`` of the linear
        representation of the sequence induced by ``recurrence_rules``.

        INPUT:

        - ``recurrence_rules`` -- a namedtuple generated by
          :meth:`parameters`

        OUTPUT:

        A vector.

        .. SEEALSO::

            :meth:`kRegularSequenceSpace.from_recurrence`

        TESTS:

        Stern--Brocot Sequence::

            sage: from sage.combinat.k_regular_sequence import RecurrenceParser
            sage: RP = RecurrenceParser(2, ZZ)
            sage: var('n')
            n
            sage: function('f')
            f
            sage: SB_rules = RP.parameters(
            ....:     1, 0, {(0, 0): 1, (1, 0): 1, (1, 1): 1},
            ....:     {0: 0, 1: 1, 2: 1}, 0)
            sage: RP.right(SB_rules)
            (0, 1, 1)

        Number of Unbordered Factors in the Thue--Morse Sequence::

            sage: M, m, coeffs, initial_values = RP.parse_recurrence([
            ....:     f(8*n) == 2*f(4*n),
            ....:     f(8*n + 1) == f(4*n + 1),
            ....:     f(8*n + 2) == f(4*n + 1) + f(4*n + 3),
            ....:     f(8*n + 3) == -f(4*n + 1) + f(4*n + 2),
            ....:     f(8*n + 4) == 2*f(4*n + 2),
            ....:     f(8*n + 5) == f(4*n + 3),
            ....:     f(8*n + 6) == -f(4*n + 1) + f(4*n + 2) + f(4*n + 3),
            ....:     f(8*n + 7) == 2*f(4*n + 1) + f(4*n + 3),
            ....:     f(0) == 1, f(1) == 2, f(2) == 2, f(3) == 4, f(4) == 2,
            ....:     f(5) == 4, f(6) == 6, f(7) == 0, f(8) == 4, f(9) == 4,
            ....:     f(10) == 4, f(11) == 4, f(12) == 12, f(13) == 0, f(14) == 4,
            ....:     f(15) == 4, f(16) == 8, f(17) == 4, f(18) == 8, f(19) == 0,
            ....:     f(20) == 8, f(21) == 4, f(22) == 4, f(23) == 8], f, n)
            sage: UB_rules = RP.parameters(
            ....:     M, m, coeffs, initial_values, 3)
            sage: RP.right(UB_rules)
            (1, 1, 2, 1, 2, 2, 4, 2, 4, 6, 0, 4, 4, 1, 0, 0)
        """
        from sage.modules.free_module_element import vector

        n1 = recurrence_rules.n1
        right = self.v_eval_n(recurrence_rules, 0)

        if n1 >= 1:
            right = vector(list(right) + [1] + (n1 - 1)*[0])

        return right

    def __call__(self, equations, function, var, *, offset=0):
        r"""
        Construct a `k`-linear representation that fulfills the recurrence relations
        given in ``equations``.

        This is the main method of :class:`RecurrenceParser` and
        is called by :meth:`kRegularSequenceSpace.from_recurrence`
        to construct a :class:`kRegularSequence`.

        INPUT:

        All parameters are explained in the high-level method
        :meth:`kRegularSequenceSpace.from_recurrence`

        OUTPUT:

        A linear representation ``(left, mu, right)``.

        Many examples can be found in
        :meth:`kRegularSequenceSpace.from_recurrence`.

        TESTS::

            sage: from sage.combinat.k_regular_sequence import RecurrenceParser
            sage: RP = RecurrenceParser(2, ZZ)
            sage: var('n')
            n
            sage: function('f')
            f
            sage: RP([f(2*n) == f(n), f(2*n + 1) == f(n) + f(n + 1),
            ....:     f(0) == 0, f(1) == 1], f, n)
            ([
              [1 0 0]  [1 1 0]
              [1 1 0]  [0 1 0]
              [0 1 0], [0 1 1]
             ],
             (1, 0, 0),
             (0, 1, 1))
        """
        from sage.arith.srange import srange

        k = self.k
        M, m, coeffs, initial_values = self.parse_recurrence(equations, function, var)
        recurrence_rules = self.parameters(M, m, coeffs, initial_values, offset)

        mu = [self.matrix(recurrence_rules, rem)
              for rem in srange(k)]

        left = self.left(recurrence_rules)
        right = self.right(recurrence_rules)

        return (mu, left, right)
