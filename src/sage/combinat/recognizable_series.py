r"""
Recognizable Series

Let `A` be an alphabet and `K` a semiring. Then a formal series `S`
with coefficients in `K` and indices in the words `A^*` is called
recognizable if it has a linear representation, i.e., there exists

- a nonnegative integer `n`

and there exist

- two vectors `\mathit{left}` and `\mathit{right}` of dimension `n` and

- a morphism of monoids `\mu` from `A^*` to `n\times n` matrices over `K`

such that the coefficient corresponding to a word `w\in A^*` equals

.. MATH::

    \mathit{left} \, \mu(w) \, \mathit{right}.

.. NOTE::

    Whenever a minimization (:meth:`~RecognizableSeries.minimized`) of
    a series needs to be computed, it is required that `K` is a field.
    In particular, minimization is called before checking if a series is
    nonzero.

.. WARNING::

    As this code is experimental, warnings are thrown when a
    recognizable series space is created for the first time in a
    session (see :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
        doctest:...: FutureWarning: This class/method/function is
        marked as experimental. It, its functionality or its interface
        might change without a formal deprecation.
        See http://trac.sagemath.org/21202 for details.


Various
=======

.. SEEALSO::

    :mod:`k-regular sequence <sage.combinat.k_regular_sequence>`,
    :mod:`sage.rings.cfinite_sequence`,
    :mod:`sage.combinat.binary_recurrence_sequences`.


AUTHORS:

- Daniel Krenn (2016, 2021)


ACKNOWLEDGEMENT:

- Daniel Krenn is supported by the
  Austrian Science Fund (FWF): P 24644-N26.


Classes and Methods
===================
"""
#*****************************************************************************
#       Copyright (C) 2016 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from functools import wraps

from sage.misc.cachefunc import cached_method
from sage.misc.superseded import experimental
from sage.structure.element import ModuleElement
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class PrefixClosedSet(object):
    def __init__(self, words):
        r"""
        A prefix-closed set.

        Creation of this prefix-closed set is interactive
        iteratively.

        INPUT:

        - ``words`` -- a class of words
          (instance of :class:`~sage.combinat.words.words.Words`)

        EXAMPLES::

            sage: from sage.combinat.recognizable_series import PrefixClosedSet
            sage: P = PrefixClosedSet(Words([0, 1], infinite=False)); P
            [word: ]

            sage: P = PrefixClosedSet.create_by_alphabet([0, 1]); P
            [word: ]

        See :meth:`iterate_possible_additions` for further examples.
        """
        self.words = words
        self.elements = [self.words([])]

    @classmethod
    def create_by_alphabet(cls, alphabet):
        r"""
        A prefix-closed set

        This is a convenience method for the
        creation of prefix-closed sets by specifying an alphabet.

        INPUT:

        - ``alphabet`` -- finite words over this ``alphabet``
          will used

        EXAMPLES::

            sage: from sage.combinat.recognizable_series import PrefixClosedSet
            sage: P = PrefixClosedSet.create_by_alphabet([0, 1]); P
            [word: ]
        """
        from sage.combinat.words.words import Words
        return cls(Words(alphabet, infinite=False))

    def __repr__(self):
        r"""
        A representation string of this prefix-closed set

        OUTPUT:

        A string

        EXAMPLES::

            sage: from sage.combinat.recognizable_series import PrefixClosedSet
            sage: P = PrefixClosedSet.create_by_alphabet([0, 1])
            sage: repr(P)  # indirect doctest
            '[word: ]'
        """
        return repr(self.elements)

    def add(self, w, check=True):
        r"""
        Add a word to this prefix-closed set.

        INPUT:

        - ``w`` -- a word

        - ``check`` -- boolean (default: ``True``). If set, then it is verified
          whether all proper prefixes of ``w`` are already in this
          prefix-closed set.

        OUTPUT:

        Nothing, but a
        :python:`RuntimeError<library/exceptions.html#exceptions.ValueError>`
        is raised if the check fails.

        EXAMPLES::

            sage: from sage.combinat.recognizable_series import PrefixClosedSet
            sage: P = PrefixClosedSet.create_by_alphabet([0, 1])
            sage: W = P.words
            sage: P.add(W([0])); P
            [word: , word: 0]
            sage: P.add(W([0, 1])); P
            [word: , word: 0, word: 01]
            sage: P.add(W([1, 1]))
            Traceback (most recent call last):
            ...
            ValueError: Cannot add as not all prefixes of 11 are included yet.
        """
        if check and any(p not in self.elements
                         for p in w.prefixes_iterator()
                         if p != w):
            raise ValueError('Cannot add as not all prefixes of '
                             '{} are included yet.'.format(w))
        self.elements.append(w)

    def iterate_possible_additions(self):
        r"""
        Return an iterator over all elements including possible new elements.

        OUTPUT:

        An iterator

        EXAMPLES::

            sage: from sage.combinat.recognizable_series import PrefixClosedSet
            sage: P = PrefixClosedSet.create_by_alphabet([0, 1]); P
            [word: ]
            sage: for n, p in enumerate(P.iterate_possible_additions()):
            ....:     print('{}?'.format(p))
            ....:     if n in (0, 2, 3, 5):
            ....:         P.add(p)
            ....:         print('...added')
            0?
            ...added
            1?
            00?
            ...added
            01?
            ...added
            000?
            001?
            ...added
            010?
            011?
            0010?
            0011?
            sage: P.elements
            [word: , word: 0, word: 00, word: 01, word: 001]

        Calling the iterator once more, returns all elements::

            sage: list(P.iterate_possible_additions())
            [word: 0,
             word: 1,
             word: 00,
             word: 01,
             word: 000,
             word: 001,
             word: 010,
             word: 011,
             word: 0010,
             word: 0011]

        The method :meth:`iterate_possible_additions` is roughly equivalent to
        ::

            sage: list(p + a
            ....:      for p in P.elements
            ....:      for a in P.words.iterate_by_length(1))
            [word: 0,
             word: 1,
             word: 00,
             word: 01,
             word: 000,
             word: 001,
             word: 010,
             word: 011,
             word: 0010,
             word: 0011]

        However, the above does not allow to add elements during iteration,
        whereas :meth:`iterate_possible_additions` does.
        """
        n = 0
        it = self.words.iterate_by_length(1)
        while n < len(self.elements):
            try:
                nn = next(it)
                yield self.elements[n] + nn #next(it)
            except StopIteration:
                n += 1
                it = self.words.iterate_by_length(1)

    def prefix_set(self):
        r"""
        Return the set of minimal (with respect to prefix ordering) elements
        of the complement of this prefix closed set.

        See also Proposition 2.3.1 of [BR2010a]_.

        OUTPUT:

        A list

        EXAMPLES::

            sage: from sage.combinat.recognizable_series import PrefixClosedSet
            sage: P = PrefixClosedSet.create_by_alphabet([0, 1]); P
            [word: ]
            sage: for n, p in enumerate(P.iterate_possible_additions()):
            ....:     if n in (0, 1, 2, 4, 6):
            ....:         P.add(p)
            sage: P
            [word: , word: 0, word: 1, word: 00, word: 10, word: 000]
            sage: P.prefix_set()
            [word: 01, word: 11, word: 001, word: 100,
             word: 101, word: 0000, word: 0001]
        """
        return [p + a
                for p in self.elements
                for a in self.words.iterate_by_length(1)
                if p + a not in self.elements]


def minimize_result(operation):
    r"""
    A decorator for operations that enables control of
    automatic minimization on the result.

    INPUT:

    - ``operation`` -- a method

    OUTPUT:

    A method with the following additional argument:

    - ``minimize`` -- (default: ``None``) a boolean or ``None``.
      If ``True``, then :meth:`minimized` is called after the operation,
      if ``False``, then not. If this argument is ``None``, then
      the default specified by the parent's ``minimize_results`` is used.

    TESTS::

        sage: from sage.combinat.recognizable_series import minimize_result
        sage: class P(object):
        ....:     pass
        sage: p = P()
        sage: class S(object):
        ....:     def __init__(self, s):
        ....:         self.s = s
        ....:     def __repr__(self):
        ....:         return self.s
        ....:     def parent(self):
        ....:         return p
        ....:     def minimized(self):
        ....:         return S(self.s + ' minimized')
        ....:     @minimize_result
        ....:     def operation(self):
        ....:         return S(self.s + ' result')

        sage: p.minimize_results = True
        sage: S('some').operation()
        some result minimized
        sage: S('some').operation(minimize=True)
        some result minimized
        sage: S('some').operation(minimize=False)
        some result

        sage: p.minimize_results = False
        sage: S('some').operation()
        some result
        sage: S('some').operation(minimize=True)
        some result minimized
        sage: S('some').operation(minimize=False)
        some result
    """
    @wraps(operation)
    def minimized(self, *args, **kwds):
        minimize = kwds.pop('minimize', None)
        if minimize is None:
            minimize = self.parent().minimize_results

        result = operation(self, *args, **kwds)

        if minimize:
            result = result.minimized()

        return result

    return minimized


class RecognizableSeries(ModuleElement):
    def __init__(self, parent, mu, left, right):
        r"""
        A recognizable series.

        - ``parent`` -- an instance of :class:`RecognizableSeriesSpace`

        - ``mu`` -- a family of square matrices, all of which have the
          same dimension.
          The indices of this family are the elements of the alphabet.
          ``mu`` may be a list or tuple of the same cardinality as the
          alphabet as well. See also :meth:`mu <mu>`.

        - ``left`` -- a vector. When evaluating a
          coefficient, this vector is multiplied from the left to the
          matrix obtained from :meth:`mu <mu>` applying on a word.
          See also :meth:`left <left>`.

        - ``right`` -- a vector. When evaluating a
          coefficient, this vector is multiplied from the right to the
          matrix obtained from :meth:`mu <mu>` applying on a word.
          See also :meth:`right <right>`.

        When created via the parent :class:`RecognizableSeriesSpace`, then
        the following option is available.

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: S = Rec((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:         vector([0, 1]), vector([1, 0])).transposed(); S
            [1] + 3*[01] + [10] + 5*[11] + 9*[001] + 3*[010] + ...

        We can access coefficients by
        ::

            sage: W = Rec.indices()
            sage: S[W([0, 0, 1])]
            9

        .. SEEALSO::

            :doc:`recognizable series <recognizable_series>`,
            :class:`RecognizableSeriesSpace`.

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, (0,1))
            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: Rec((M0, M1), (0, 1), (1, 1))
            [] + [0] + 3*[1] + [00] + 3*[01] + 3*[10] + 5*[11] + [000] + 3*[001] + 3*[010] + ...
        """
        super(RecognizableSeries, self).__init__(parent=parent)

        from sage.modules.free_module_element import vector
        from sage.sets.family import Family

        A = self.parent().alphabet()
        if isinstance(mu, (list, tuple)):
            mu = dict(zip(A, mu))
        mu = Family(mu)
        if not mu.is_finite():
            raise NotImplementedError('mu is not a finite family of matrices.')

        self._left_ = vector(left)
        self._mu_ = mu
        self._right_ = vector(right)

    @property
    def mu(self):
        r"""
        When evaluating a coefficient, this is applied on each letter
        of a word; the result is a matrix.
        This extends :meth:`mu <mu>` to words over the parent's
        :meth:`~RecognizableSeriesSpace.alphabet`.

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: S = Rec((M0, M1), [0, 1], [1, 1])
            sage: S.mu[0] == M0 and S.mu[1] == M1
            True
        """
        return self._mu_

    @property
    def left(self):
        r"""
        When evaluating a coefficient, this vector is multiplied from
        the left to the matrix obtained from :meth:`mu <mu>` applied on a
        word.

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: Rec((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:     vector([0, 1]), vector([1, 0])).transposed().left
            (1, 0)
        """
        return self._left_

    @property
    def right(self):
        r"""
        When evaluating a coefficient, this vector is multiplied from
        the right to the matrix obtained from :meth:`mu <mu>` applied on a
        word.

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: Rec((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:     vector([0, 1]), vector([1, 0])).transposed().right
            (0, 1)
        """
        return self._right_

    def linear_representation(self):
        r"""
        Return the linear representation of this series.

        OUTPUT:

        A triple ``(left, mu, right)`` containing
        the vectors :meth:`left <left>` and :meth:`right <right>`,
        and the family of matrices :meth:`mu <mu>`.

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: Rec((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:     vector([0, 1]), vector([1, 0])
            ....:    ).transposed().linear_representation()
            ((1, 0),
             Finite family {0: [3 0]
                               [6 1],
                            1: [ 0  1]
                               [-6  5]},
             (0, 1))
        """
        return (self.left, self.mu, self.right)

    def _repr_(self, latex=False):
        r"""
        A representation string for this recognizable series.

        INPUT:

        - ``latex`` -- (default: ``False``) a boolean. If set, then
          LaTeX-output is returned.

        OUTPUT:

        A string

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: S = Rec((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:         vector([0, 1]), vector([1, 0])).transposed()
            sage: repr(S)  # indirect doctest
            '[1] + 3*[01] + [10] + 5*[11] + 9*[001] + 3*[010] + ...'

        TESTS::

            sage: S = Rec((Matrix([[0]]), Matrix([[0]])),
            ....:         vector([1]), vector([1]))
            sage: repr(S)  # indirect doctest
            '[] + ...'

            sage: S = Rec((Matrix([[0, 1], [0, 0]]), Matrix([[0, 0], [0, 0]])),
            ....:         vector([0, 1]), vector([1, 0]))
            sage: repr(S)  # indirect doctest
            '0 + ...'
        """
        if self.is_trivial_zero():
            return '0'

        from itertools import islice

        if latex:
            from sage.misc.latex import latex as latex_repr
            fr = latex_repr
            fs = latex_repr
            times = ' '
        else:
            fr = repr
            fs = str
            times = '*'

        def summand(w, c):
            if c == 1:
                return '[{w}]'.format(w=fs(w))
            return '{c}{times}[{w}]'.format(c=fr(c), times=times, w=fs(w))

        def all_coefficients():
            number_of_zeros = 0
            for w in self.parent().indices():
                c = self[w]
                if c != 0:
                    number_of_zeros = 0
                    yield (w, self[w])
                else:
                    number_of_zeros += 1
                if number_of_zeros >= 100:
                    return

        coefficients = islice(all_coefficients(), 10)

        s = ' + '.join(summand(w, c)
                       for w, c in coefficients)
        s = s.replace('+ -', '- ')
        if not s:
            s = '0'
        return s + ' + ...'

    def _latex_(self):
        r"""
        A LaTeX-representation string for this recognizable series.

        OUTPUT:

        A string

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: S = Rec((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:         vector([0, 1]), vector([1, 0])).transposed()
            sage: latex(S)  # indirect doctest
            [1] + 3 [01] + [10] + 5 [11] + 9 [001] + 3 [010]
                + 15 [011] + [100] + 11 [101] + 5 [110] + ...
        """
        return self._repr_(latex=True)

    @cached_method
    def coefficient_of_word(self, w, multiply_left=True, multiply_right=True):
        r"""
        Return the coefficient to word `w` of this series.

        INPUT:

        - ``w`` -- a word over the parent's
          :meth:`~RecognizableSeriesSpace.alphabet`

        - ``multiply_left`` -- (default: ``True``) a boolean. If ``False``,
          then multiplication by :meth:`left <left>` is skipped.

        - ``multiply_right`` -- (default: ``True``) a boolean. If ``False``,
          then multiplication by :meth:`right <right>` is skipped.

        OUTPUT:

        An element in the parent's
        :meth:`~RecognizableSeriesSpace.coefficient_ring`

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: W = Rec.indices()
            sage: S = Rec((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
            ....:          left=vector([0, 1]), right=vector([1, 0]))
            sage: S[W(7.digits(2))]  # indirect doctest
            3

        TESTS::

            sage: w = W(6.digits(2))
            sage: S.coefficient_of_word(w)
            2
            sage: S.coefficient_of_word(w, multiply_left=False)
            (-1, 2)
            sage: S.coefficient_of_word(w, multiply_right=False)
            (2, 3)
            sage: S.coefficient_of_word(w, multiply_left=False, multiply_right=False)
            [-1 -2]
            [ 2  3]
        """
        result = self._mu_of_word_(w)
        if multiply_left:
            result = self.left * result
        if multiply_right:
            result = result * self.right
        return result

    __getitem__ = coefficient_of_word

    @cached_method
    def _mu_of_empty_word_(self):
        r"""
        Return :meth:`mu <mu>` applied on the empty word.

        OUTPUT:

        A matrix

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: W = Rec.indices()
            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: S = Rec({W([0]): M0, W([1]): M1}, [0, 1], [1, 1])
            sage: S._mu_of_empty_word_()
            [1 0]
            [0 1]
            sage: I = Matrix([[1, 0], [0, 1]])
            sage: T = Rec({W([]): I, W([0]): M0, W([1]): M1}, [0, 1], [1, 1])
            sage: T._mu_of_empty_word_()
            [1 0]
            [0 1]
            sage: _ is I
            True
        """
        eps = self.parent().indices()()
        try:
            return self.mu[eps]
        except KeyError:
            return next(iter(self.mu)).parent().one()

    @cached_method
    def _mu_of_word_(self, w):
        r"""
        Return :meth:`mu <mu>` applied on the word `w`.

        INPUT:

        - ``w`` -- a word over the parent's
          :meth:`~RecognizableSeriesSpace.alphabet`

        OUTPUT:

        A matrix

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: W = Rec.indices()
            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: S = Rec((M0, M1), [0, 1], [1, 1])
            sage: S._mu_of_word_(W([0])) == M0
            True
            sage: S._mu_of_word_(W([1])) == M1
            True
            sage: S._mu_of_word_(W(3.digits(2))) == M1^2
            True

        ::

            sage: S._mu_of_word_(-1)
            Traceback (most recent call last):
            ...
            ValueError: Index -1 is not in Finite words over {0, 1}.
        """
        W = self.parent().indices()
        if w not in W:
            raise ValueError('Index {} is not in {}.'.format(w, W))
        from sage.misc.misc_c import prod
        return prod((self.mu[a] for a in w), z=self._mu_of_empty_word_())

    def __iter__(self):
        r"""
        Return an iterator over pairs ``(index, coefficient)``.

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: S = Rec((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
            ....:         left=vector([0, 1]), right=vector([1, 0]))
            sage: from itertools import islice
            sage: list(islice(S, 10))
            [(word: , 0),
             (word: 0, 0),
             (word: 1, 1),
             (word: 00, 0),
             (word: 01, 1),
             (word: 10, 1),
             (word: 11, 2),
             (word: 000, 0),
             (word: 001, 1),
             (word: 010, 1)]
            sage: list(islice((s for s in S if s[1] != 0), 10))
            [(word: 1, 1),
             (word: 01, 1),
             (word: 10, 1),
             (word: 11, 2),
             (word: 001, 1),
             (word: 010, 1),
             (word: 011, 2),
             (word: 100, 1),
             (word: 101, 2),
             (word: 110, 2)]

            sage: S = Rec((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
            ....:         left=vector([1, 0]), right=vector([1, 0]))
            sage: list(islice((s for s in S if s[1] != 0), 10))
            [(word: , 1),
             (word: 0, 1),
             (word: 00, 1),
             (word: 11, -1),
             (word: 000, 1),
             (word: 011, -1),
             (word: 101, -1),
             (word: 110, -1),
             (word: 111, -2),
             (word: 0000, 1)]

        TESTS::

            sage: it = iter(S)
            sage: iter(it) is it
            True
            sage: iter(S) is not it
            True
        """
        return iter((w, self[w]) for w in self.parent().indices())

    def is_trivial_zero(self):
        r"""
        Return whether this recognizable series is trivially equal to
        zero (without any :meth:`minimization <minimized>`).

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:     left=vector([0, 1]), right=vector([1, 0])).is_trivial_zero()
            False
            sage: Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:     left=vector([0, 0]), right=vector([1, 0])).is_trivial_zero()
            True
            sage: Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:     left=vector([0, 1]), right=vector([0, 0])).is_trivial_zero()
            True

        The following two differ in the coefficient of the empty word::

            sage: Rec((Matrix([[0, 0], [0, 0]]), Matrix([[0, 0], [0, 0]])),
            ....:     left=vector([0, 1]), right=vector([1, 0])).is_trivial_zero()
            True
            sage: Rec((Matrix([[0, 0], [0, 0]]), Matrix([[0, 0], [0, 0]])),
            ....:     left=vector([1, 1]), right=vector([1, 1])).is_trivial_zero()
            False

        TESTS::

            sage: Rec.zero().is_trivial_zero()
            True

        The following is zero, but not trivially zero::

            sage: S = Rec((Matrix([[1, 0], [0, 0]]), Matrix([[1, 0], [0, 0]])),
            ....:         left=vector([0, 1]), right=vector([1, 0]))
            sage: S.is_trivial_zero()
            False
            sage: S.is_zero()
            True
        """
        return not self.left or not self.right or \
            (all(not self.mu[a] for a in self.parent().alphabet()) and
             not self[self.parent().indices()()])

    def __bool__(self):
        r"""
        Return whether this recognizable series is nonzero.

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: bool(Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:          left=vector([0, 1]), right=vector([1, 0])))
            False
            sage: bool(Rec((Matrix([[0, 0], [0, 0]]), Matrix([[0, 0], [0, 0]])),
            ....:          left=vector([0, 1]), right=vector([1, 0])))
            False
            sage: bool(Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:          left=vector([0, 0]), right=vector([1, 0])))
            False
            sage: bool(Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:          left=vector([0, 1]), right=vector([0, 0])))
            False

        ::

            sage: S = Rec((Matrix([[1, 0], [0, 0]]), Matrix([[1, 0], [0, 0]])),
            ....:         left=vector([0, 1]), right=vector([1, 0]))
            sage: bool(S)
            False
        """
        if self.is_trivial_zero():
            return False
        try:
            M = self.minimized()
        except ValueError:
            pass
        else:
            if M.is_trivial_zero():
                return False
        return True

    __nonzero__ = __bool__

    def transposed(self):
        r"""
        Return the transposed series.

        OUTPUT:

        A :class:`RecognizableSeries`

        Each of the matrices in :meth:`mu <mu>` is transposed. Additionally
        the vectors :meth:`left <left>` and :meth:`right <right>` are switched.

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: S = Rec((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:         vector([0, 1]), vector([1, 0])).transposed()
            sage: S
            [1] + 3*[01] + [10] + 5*[11] + 9*[001] + 3*[010]
                + 15*[011] + [100] + 11*[101] + 5*[110] + ...
            sage: S.mu[0], S.mu[1], S.left, S.right
            (
            [3 0]  [ 0  1]
            [6 1], [-6  5], (1, 0), (0, 1)
            )
            sage: T = S.transposed()
            sage: T
            [1] + [01] + 3*[10] + 5*[11] + [001] + 3*[010]
                + 5*[011] + 9*[100] + 11*[101] + 15*[110] + ...
            sage: T.mu[0], T.mu[1], T.left, T.right
            (
            [3 6]  [ 0 -6]
            [0 1], [ 1  5], (0, 1), (1, 0)
            )
        """
        return self.parent()(self.mu.map(lambda M: M.transpose()),
                             left=self.right,
                             right=self.left)

    @cached_method
    def minimized(self):
        r"""
        Return a recognizable series equivalent to this series, but
        with a minimized linear representation.

        The coefficients of the involved matrices need be in a field.
        If this is not the case, then the coefficients are
        automatically coerced to their fraction field.

        OUTPUT:

        A :class:`RecognizableSeries`

        ALGORITHM:

        This method implements the minimization algorithm presented in
        Chapter 2 of [BR2010a]_.

        EXAMPLES::

            sage: from itertools import islice
            sage: from six.moves import zip
            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])

            sage: S = Rec((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:         vector([0, 1]), vector([1, 0])).transposed()
            sage: S
            [1] + 3*[01] + [10] + 5*[11] + 9*[001] + 3*[010]
                + 15*[011] + [100] + 11*[101] + 5*[110] + ...
            sage: M = S.minimized()
            sage: M.mu[0], M.mu[1], M.left, M.right
            (
            [3 0]  [ 0  1]
            [6 1], [-6  5], (1, 0), (0, 1)
            )
            sage: all(c == d and v == w
            ....:     for (c, v), (d, w) in islice(zip(iter(S), iter(M)), 20))
            True

            sage: S = Rec((Matrix([[2, 0], [1, 1]]), Matrix([[2, 0], [2, 1]])),
            ....:         vector([1, 0]), vector([1, 1]))
            sage: S
            [] + 2*[0] + 2*[1] + 4*[00] + 4*[01] + 4*[10] + 4*[11]
               + 8*[000] + 8*[001] + 8*[010] + ...
            sage: M = S.minimized()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([2], [2], (1), (1))
            sage: all(c == d and v == w
            ....:     for (c, v), (d, w) in islice(zip(iter(S), iter(M)), 20))
            True
        """
        return self._minimized_right_()._minimized_left_()

    def _minimized_right_(self):
        r"""
        Return a recognizable series equivalent to this series, but
        with a right minimized linear representation.

        OUTPUT:

        A :class:`RecognizableSeries`

        See :meth:`minimized` for details.

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: S = Rec((Matrix([[0, 0], [0, 0]]), Matrix([[0, 0], [0, 0]])),
            ....:         vector([1, 1]), vector([1, 1]))
            sage: M = S._minimized_right_()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([0], [0], (2), (1))
        """
        return self.transposed()._minimized_left_().transposed()

    def _minimized_left_(self):
        r"""
        Return a recognizable series equivalent to this series, but
        with a left minimized linear representation.

        OUTPUT:

        A :class:`RecognizableSeries`

        See :meth:`minimized` for details.

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: S = Rec((Matrix([[0, 0], [0, 0]]), Matrix([[0, 0], [0, 0]])),
            ....:         vector([1, 1]), vector([1, 1]))
            sage: M = S._minimized_left_()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([0], [0], (1), (2))
            sage: M = S.minimized()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([0], [0], (1), (2))

        ::

            sage: S = Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:         vector([1, -1]), vector([1, 1]))._minimized_left_()
            sage: S.mu[0], S.mu[1], S.left, S.right
            ([1], [1], (1), (0))
            sage: M = S.minimized()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([], [], (), ())

            sage: S = Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:         vector([1, 1]), vector([1, -1]))
            sage: M = S._minimized_left_()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([1], [1], (1), (0))
            sage: M = S.minimized()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([], [], (), ())

            sage: S = Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:         left=vector([0, 1]), right=vector([1, 0]))
            sage: M = S._minimized_left_()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([1], [1], (1), (0))
            sage: M = S.minimized()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([], [], (), ())
        """
        from sage.matrix.constructor import Matrix
        from sage.modules.free_module_element import vector
        from sage.rings.integer_ring import ZZ

        pcs = PrefixClosedSet(self.parent().indices())
        left = self.coefficient_of_word(pcs.elements[0], multiply_right=False)
        if left.is_zero():
            return self.parent().zero()
        Left = [left]
        for p in pcs.iterate_possible_additions():
            left = self.coefficient_of_word(p,
                                            multiply_left=True,
                                            multiply_right=False)
            try:
                Matrix(Left).solve_left(left)
            except ValueError:
                # no solution found
                pcs.add(p)
                Left.append(left)
        P = pcs.elements
        C = pcs.prefix_set()

        ML = Matrix(Left)

        def alpha(c):
            return ML.solve_left(self.coefficient_of_word(c, multiply_right=False))

        mu_prime = []
        for a in self.parent().alphabet():
            a = self.parent().indices()([a])
            M = Matrix([alpha(c) if c in C else tuple(ZZ(c==q) for q in P)
                        for c in (p + a for p in P)])
            mu_prime.append(M)

        left_prime = vector([ZZ(1)] + (len(P)-1)*[ZZ(0)])
        right_prime = vector(self.coefficient_of_word(p) for p in P)

        P = self.parent()
        return P.element_class(P, mu_prime, left_prime, right_prime)

    def dimension(self):
        r"""
        Return the dimension of this recognizable series.

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:     left=vector([0, 1]), right=vector([1, 0])).dimension()
            2
        """
        return self.mu.first().nrows()

    @minimize_result
    def _add_(self, other):
        r"""
        Return the sum of this recognizable series and the ``other``
        recognizable series.

        INPUT:

        - ``other`` -- a :class:`RecognizableSeries` with the same parent
          as this recognizable series

        - ``minimize`` -- (default: ``None``) a boolean or ``None``.
          If ``True``, then :meth:`minimized` is called after the operation,
          if ``False``, then not. If this argument is ``None``, then
          the default specified by the parent's ``minimize_results`` is used.

        OUTPUT:

        A :class:`RecognizableSeries`.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: E = Seq2((Matrix([[0, 1], [0, 1]]), Matrix([[0, 0], [0, 1]])),
            ....:          vector([1, 0]), vector([1, 1]))
            sage: E
            2-regular sequence 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, ...
            sage: O = Seq2((Matrix([[0, 0], [0, 1]]), Matrix([[0, 1], [0, 1]])),
            ....:          vector([1, 0]), vector([0, 1]))
            sage: O
            2-regular sequence 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, ...
            sage: I = E + O  # indirect doctest
            sage: I
            2-regular sequence 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
            sage: I.linear_representation()
            ((1),
             Finite family {0: [1],
                            1: [1]},
             (1))
        """
        from sage.modules.free_module_element import vector
        P = self.parent()

        result = P.element_class(
            P,
            dict((a, self.mu[a].block_sum(other.mu[a])) for a in P.alphabet()),
            vector(tuple(self.left) + tuple(other.left)),
            vector(tuple(self.right) + tuple(other.right)))

        return result

    def _neg_(self):
        r"""
        Return the additive inverse of this recognizable series.

        OUTPUT:

        A :class:`RecognizableSeries`

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: E = Seq2((Matrix([[0, 1], [0, 1]]), Matrix([[0, 0], [0, 1]])),
            ....:          vector([1, 0]), vector([1, 1]))
            sage: -E
            2-regular sequence -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, ...
            sage: Z = E - E
            sage: Z.is_trivial_zero()
            True
        """
        P = self.parent()
        return P.element_class(P, self.mu, -self.left, self.right)

    def _rmul_(self, other):
        r"""
        Multiply this recognizable series from the right
        by an element ``other`` of its coefficient (semi-)ring.

        INPUT:

        - ``other`` -- an element of the coefficient (semi-)ring

        OUTPUT:

        A :class:`RecognizableSeries`

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: E = Seq2((Matrix([[0, 1], [0, 1]]), Matrix([[0, 0], [0, 1]])),
            ....:          vector([1, 0]), vector([1, 1]))
            sage: M = 2 * E  # indirect doctest
            sage: M
            2-regular sequence 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, ...
            sage: M.linear_representation()
            ((2, 0),
             Finite family {0: [0 1]
                               [0 1],
                            1: [0 0]
                               [0 1]},
             (1, 1))

        TESTS::

            sage: 1 * E is E
            True

        We test that ``_rmul_`` and ``_lmul_`` are actually called::

            sage: def print_name(f):
            ....:     def f_with_printed_name(*args, **kwds):
            ....:         print(f.__name__)
            ....:         return f(*args, **kwds)
            ....:     return f_with_printed_name

            sage: E._rmul_ = print_name(E._rmul_)
            sage: E._lmul_ = print_name(E._lmul_)
            sage: 2 * E
            _rmul_
            2-regular sequence 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, ...
            sage: E * 2
            _lmul_
            _lmul_
            2-regular sequence 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, ...
        """
        if other.is_one():
            return self
        P = self.parent()
        return P.element_class(P, self.mu, other*self.left, self.right)

    def _lmul_(self, other):
        r"""
        Multiply this recognizable series from the left
        by an element ``other`` of its coefficient (semi-)ring.

        INPUT:

        - ``other`` -- an element of the coefficient (semi-)ring

        OUTPUT:

        A :class:`RecognizableSeries`

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: E = Seq2((Matrix([[0, 1], [0, 1]]), Matrix([[0, 0], [0, 1]])),
            ....:          vector([1, 0]), vector([1, 1]))
            sage: M = E * 2 # indirect doctest
            sage: M
            2-regular sequence 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, ...
            sage: M.linear_representation()
            ((1, 0),
             Finite family {0: [0 1]
                               [0 1],
                            1: [0 0]
                               [0 1]},
             (2, 2))

        TESTS::

            sage: E * 1 is E
            True

        The following is not tested, as `MS^i` for integers `i` does
        not work, thus ``vector([m])`` fails. (See :trac:`21317` for
        details.)

        ::

            sage: MS = MatrixSpace(ZZ,2,2)
            sage: Rec = RecognizableSeriesSpace(MS, [0, 1])
            sage: m = MS.an_element()
            sage: S = Rec((Matrix([[m]]), Matrix([[m]])),  # not tested
            ....:         vector([m]), vector([m]))
            sage: S  # not tested
            sage: M = m * S  # not tested indirect doctest
            sage: M  # not tested
            sage: M.linear_representation()  # not tested
        """
        if other.is_one():
            return self
        P = self.parent()
        return P.element_class(P, self.mu, self.left, self.right*other)

    @minimize_result
    def hadamard_product(self, other):
        r"""
        Return the Hadamard product of this recognizable series
        and the ``other`` recognizable series, i.e., multiply the two
        series coefficient-wise.

        INPUT:

        - ``other`` -- a :class:`RecognizableSeries` with the same parent
          as this recognizable series

        - ``minimize`` -- (default: ``None``) a boolean or ``None``.
          If ``True``, then :meth:`minimized` is called after the operation,
          if ``False``, then not. If this argument is ``None``, then
          the default specified by the parent's ``minimize_results`` is used.

        OUTPUT:

        A :class:`RecognizableSeries`

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)

            sage: E = Seq2((Matrix([[0, 1], [0, 1]]), Matrix([[0, 0], [0, 1]])),
            ....:          vector([1, 0]), vector([1, 1]))
            sage: E
            2-regular sequence 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, ...

            sage: O = Seq2((Matrix([[0, 0], [0, 1]]), Matrix([[0, 1], [0, 1]])),
            ....:          vector([1, 0]), vector([0, 1]))
            sage: O
            2-regular sequence 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, ...

            sage: C = Seq2((Matrix([[2, 0], [2, 1]]), Matrix([[0, 1], [-2, 3]])),
            ....:          vector([1, 0]), vector([0, 1]))
            sage: C
            2-regular sequence 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, ...

        ::

            sage: CE = C.hadamard_product(E)
            sage: CE
            2-regular sequence 0, 0, 2, 0, 4, 0, 6, 0, 8, 0, ...
            sage: CE.linear_representation()
            ((1, 0, 0),
             Finite family {0: [0 1 0]
                               [0 2 0]
                               [0 2 1],
                            1: [ 0  0  0]
                               [ 0  0  1]
                               [ 0 -2  3]},
             (0, 0, 2))

            sage: Z = E.hadamard_product(O)
            sage: Z
            2-regular sequence 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            sage: Z.linear_representation()
            ((),
             Finite family {0: [],
                            1: []},
             ())

        TESTS::

            sage: EC = E.hadamard_product(C, minimize=False)
            sage: EC
            2-regular sequence 0, 0, 2, 0, 4, 0, 6, 0, 8, 0, ...
            sage: EC.linear_representation()
            ((1, 0, 0, 0),
             Finite family {0: [0 0 2 0]
                               [0 0 2 1]
                               [0 0 2 0]
                               [0 0 2 1],
                            1: [ 0  0  0  0]
                               [ 0  0  0  0]
                               [ 0  0  0  1]
                               [ 0  0 -2  3]},
             (0, 1, 0, 1))
            sage: MEC = EC.minimized()
            sage: MEC
            2-regular sequence 0, 0, 2, 0, 4, 0, 6, 0, 8, 0, ...
            sage: MEC.linear_representation()
            ((1, 0, 0),
             Finite family {0: [0 1 0]
                               [0 2 0]
                               [0 2 1],
                            1: [ 0  0  0]
                               [ 0  0  1]
                               [ 0 -2  3]},
             (0, 0, 2))
        """
        from sage.matrix.constructor import Matrix
        from sage.modules.free_module_element import vector
        P = self.parent()

        def tensor_product(left, right):
            T = left.tensor_product(right)
            T.subdivide()
            return T
        result = P.element_class(
            P,
            dict((a, tensor_product(self.mu[a], other.mu[a]))
                 for a in P.alphabet()),
            vector(tensor_product(Matrix(self.left), Matrix(other.left))),
            vector(tensor_product(Matrix(self.right), Matrix(other.right))))

        return result


class RecognizableSeriesSpace(UniqueRepresentation, Parent):
    r"""
    The space of recognizable series on the given alphabet and
    with the given coefficients.

    INPUT:

    - ``coefficient_ring`` -- a (semi-)ring

    - ``alphabet`` -- a tuple, list or
      :class:`~sage.sets.totally_ordered_finite_set.TotallyOrderedFiniteSet`.
      If specified, then the ``indices`` are the
      finite words over this ``alphabet``.
      ``alphabet`` and ``indices`` cannot be specified
      at the same time.

    - ``indices`` -- a SageMath-parent of finite words over an alphabet.
      ``alphabet`` and ``indices`` cannot be specified
      at the same time.

    - ``category`` -- (default: ``None``) the category of this
      space

    EXAMPLES:

    We create a recognizable series that counts the number of ones in each word::

        sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
        sage: Rec
        Space of recognizable series on {0, 1} with coefficients in Integer Ring
        sage: Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 1], [0, 1]])),
        ....:     vector([1, 0]), vector([0, 1]))
        [1] + [01] + [10] + 2*[11] + [001] + [010] + 2*[011] + [100] + 2*[101] + 2*[110] + ...

    All of the following examples create the same space::

        sage: Rec1 = RecognizableSeriesSpace(ZZ, [0, 1])
        sage: Rec1
        Space of recognizable series on {0, 1} with coefficients in Integer Ring
        sage: Rec2 = RecognizableSeriesSpace(coefficient_ring=ZZ, alphabet=[0, 1])
        sage: Rec2
        Space of recognizable series on {0, 1} with coefficients in Integer Ring
        sage: Rec3 = RecognizableSeriesSpace(ZZ, indices=Words([0, 1], infinite=False))
        sage: Rec3
        Space of recognizable series on {0, 1} with coefficients in Integer Ring

    .. SEEALSO::

        :doc:`recognizable series <recognizable_series>`,
        :class:`RecognizableSeries`.
    """
    Element = RecognizableSeries

    @staticmethod
    def __classcall__(cls, *args, **kwds):
        r"""
        Prepare normalizing the input in order to ensure a
        unique representation.

        For more information see :class:`RecognizableSeriesSpace`
        and :meth:`__normalize__`.

        TESTS::

            sage: Rec1 = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: Rec1
            Space of recognizable series on {0, 1} with coefficients in Integer Ring
            sage: Rec2 = RecognizableSeriesSpace(coefficient_ring=ZZ, alphabet=[0, 1])
            sage: Rec2
            Space of recognizable series on {0, 1} with coefficients in Integer Ring
            sage: Rec3 = RecognizableSeriesSpace(ZZ, indices=Words([0, 1], infinite=False))
            sage: Rec3
            Space of recognizable series on {0, 1} with coefficients in Integer Ring
            sage: Rec1 is Rec2 is Rec3
            True
        """
        return super(RecognizableSeriesSpace, cls).__classcall__(
            cls, *cls.__normalize__(*args, **kwds))

    @classmethod
    def __normalize__(cls,
                      coefficient_ring=None,
                      alphabet=None, indices=None,
                      category=None,
                      minimize_results=True):
        r"""
        Normalizes the input in order to ensure a unique
        representation.

        For more information see :class:`RecognizableSeriesSpace`.

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])  # indirect doctest
            sage: Rec.category()
            Category of modules over Integer Ring
            sage: RecognizableSeriesSpace([0, 1], [0, 1])
            Traceback (most recent call last):
            ...
            ValueError: Coefficient ring [0, 1] is not a semiring.

        ::

            sage: W = Words([0, 1], infinite=False)
            sage: RecognizableSeriesSpace(ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Specify either 'alphabet' or 'indices'.
            sage: RecognizableSeriesSpace(ZZ, alphabet=[0, 1], indices=W)
            Traceback (most recent call last):
            ...
            ValueError: Specify either 'alphabet' or 'indices'.
            sage: RecognizableSeriesSpace(alphabet=[0, 1])
            Traceback (most recent call last):
            ...
            ValueError: No coefficient ring specified.
            sage: RecognizableSeriesSpace(ZZ, indices=Words(ZZ))
            Traceback (most recent call last):
            ...
            NotImplementedError: Alphabet is not finite.
        """
        if (alphabet is None) == (indices is None):
            raise ValueError("Specify either 'alphabet' or 'indices'.")

        if indices is None:
            from sage.combinat.words.words import Words
            indices = Words(alphabet, infinite=False)
        if not indices.alphabet().is_finite():
            raise NotImplementedError('Alphabet is not finite.')

        if coefficient_ring is None:
            raise ValueError('No coefficient ring specified.')
        from sage.categories.semirings import Semirings
        if coefficient_ring not in Semirings():
            raise ValueError(
                'Coefficient ring {} is not a semiring.'.format(coefficient_ring))

        from sage.categories.modules import Modules
        category = category or Modules(coefficient_ring)

        return (coefficient_ring, indices, category, minimize_results)

    @experimental(trac_number=21202)
    def __init__(self, coefficient_ring, indices, category, minimize_results):
        r"""
        See :class:`RecognizableSeriesSpace` for details.

        INPUT:

        - ``coefficients`` -- a (semi-)ring

        - ``indices`` -- a SageMath-parent of finite words over an alphabet

        - ``category`` -- (default: ``None``) the category of this
          space

        - ``minimize_results`` -- (default: ``True``) a boolean. If set, then
          :meth:`RecognizableSeries.minimized` is automatically called
          after performing operations.

        TESTS::

            sage: RecognizableSeriesSpace(ZZ, [0, 1])
            Space of recognizable series on {0, 1} with coefficients in Integer Ring
        """
        self._indices_ = indices
        self._minimize_results_ = minimize_results
        super(RecognizableSeriesSpace, self).__init__(
            category=category, base=coefficient_ring)

    def alphabet(self):
        r"""
        Return the alphabet of this recognizable series space.

        OUTPUT:

        A totally ordered set

        EXAMPLES::

            sage: RecognizableSeriesSpace(ZZ, [0, 1]).alphabet()
            {0, 1}

        TESTS::

            sage: type(RecognizableSeriesSpace(ZZ, [0, 1]).alphabet())
            <class 'sage.sets.totally_ordered_finite_set.TotallyOrderedFiniteSet_with_category'>
        """
        return self.indices().alphabet()

    def indices(self):
        r"""
        Return the indices of the recognizable series.

        OUTPUT:

        The set of finite words over the alphabet

        EXAMPLES::

            sage: RecognizableSeriesSpace(ZZ, [0, 1]).indices()
            Finite words over {0, 1}
        """
        return self._indices_

    def coefficient_ring(self):
        r"""
        Return the coefficients of this recognizable series space.

        OUTPUT:

        A (semi-)ring

        EXAMPLES::

            sage: RecognizableSeriesSpace(ZZ, [0, 1]).coefficient_ring()
            Integer Ring
        """
        return self.base()

    @property
    def minimize_results(self):
        r"""
        A boolean indicating whether
        :meth:`RecognizableSeries.minimized` is automatically called
        after performing operations.

        TESTS::

            sage: RecognizableSeriesSpace(ZZ, [0, 1]).minimize_results
            True
            sage: RecognizableSeriesSpace(ZZ, [0, 1], minimize_results=True).minimize_results
            True
            sage: RecognizableSeriesSpace(ZZ, [0, 1], minimize_results=False).minimize_results
            False
        """
        return self._minimize_results_

    def _repr_(self):
        r"""
        Return a representation string of this recognizable sequence
        space.

        OUTPUT:

        A string

        TESTS::

            sage: repr(RecognizableSeriesSpace(ZZ, [0, 1]))  # indirect doctest
            'Space of recognizable series on {0, 1} with coefficients in Integer Ring'
        """
        return 'Space of recognizable series on {} ' \
               'with coefficients in {}'.format(self.alphabet(),
                                                self.coefficient_ring())

    @cached_method
    def one_hadamard(self):
        r"""
        Return the identity with respect to the
        :meth:`~RecognizableSeries.hadamard_product`, i.e. the
        coefficient-wise multiplication.

        OUTPUT:

        A :class:`RecognizableSeries`

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: Rec.one_hadamard()
            [] + [0] + [1] + [00] + [01] + [10]
               + [11] + [000] + [001] + [010] + ...

        TESTS::

            sage: Rec.one_hadamard() is Rec.one_hadamard()
            True
        """
        from sage.matrix.constructor import Matrix
        from sage.modules.free_module_element import vector

        one = self.coefficient_ring()(1)
        return self(dict((a, Matrix([[one]])) for a in self.alphabet()),
                    vector([one]), vector([one]))

    def _element_constructor_(self, data,
                              left=None, right=None):
        r"""
        Return a recognizable series.

        See :class:`RecognizableSeriesSpace` for details.

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])

            sage: Rec.zero()
            0
            sage: type(_)
            <class 'sage.combinat.recognizable_series.RecognizableSeriesSpace_with_category.element_class'>

        ::

            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: S = Rec((M0, M1), [0, 1], [1, 1])
            sage: Rec(S) is S
            True

            sage: Rec((M0, M1))
            Traceback (most recent call last):
            ...
            ValueError: Left or right vector is None.
            sage: Rec((M0, M1), [0, 1])
            Traceback (most recent call last):
            ...
            ValueError: Left or right vector is None.
            sage: Rec((M0, M1), left=[0, 1])
            Traceback (most recent call last):
            ...
            ValueError: Left or right vector is None.
            sage: Rec((M0, M1), right=[0, 1])
            Traceback (most recent call last):
            ...
            ValueError: Left or right vector is None.
        """
        if isinstance(data, int) and data == 0:
            from sage.matrix.constructor import Matrix
            from sage.modules.free_module_element import vector
            from sage.sets.family import Family

            return self.element_class(
                self, Family(self.alphabet(), lambda a: Matrix()),
                vector([]), vector([]))

        if type(data) == self.element_class and data.parent() == self:
            element = data

        elif isinstance(data, RecognizableSeries):
            element = self.element_class(self, data.mu, data.left, data.right)

        else:
            mu = data
            if left is None or right is None:
                raise ValueError('Left or right vector is None.')

            element = self.element_class(self, mu, left, right)

        return element
