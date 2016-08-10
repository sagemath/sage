r"""
`k`-regular Sequences

An introduction and formal definition of `k`-regular sequences can be
found, for example, on the :wikipedia:`k-regular_sequence` or in
[AS2003]_.


Examples
========

Binary sum of digits
--------------------

::

    sage: Seq2 = kRegularSequenceSpace(2, ZZ)
    sage: S = Seq2((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
    ....:          initial=vector([0, 1]), selection=vector([1, 0]))
    sage: S
    2-regular sequence 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, ...
    sage: all(S[n] == sum(n.digits(2)) for n in srange(10))
    True

Number of odd entries in Pascal's triangle
------------------------------------------

::

    sage: @cached_function
    ....: def u(n):
    ....:     if n <= 1:
    ....:         return n
    ....:     return 2*u(floor(n/2)) + u(ceil(n/2))
    sage: tuple(u(n) for n in srange(10))
    (0, 1, 3, 5, 9, 11, 15, 19, 27, 29)

    sage: U = Seq2((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
    ....:          initial=vector([0, 1]), selection=vector([1, 0]), transpose=True)
    sage: all(U[n] == u(n) for n in srange(30))
    True


Various
=======

.. SEEALSO::

    :doc:`sage/rings/cfinite_sequence`,
    :doc:`sage/combinat/binary_recurrence_sequences`.

REFERENCES:

.. [AS2003] Jean-Paul Allouche, Jeffrey Shallit,
   *Automatic Sequences: Theory, Applications, Generalizations*,
   Cambridge University Press, 2003.

AUTHORS:

- Daniel Krenn (2016)

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

import itertools
from sage.misc.cachefunc import cached_function, cached_method


def pad_right(T, length, zero=0):
    r"""
    Pad ``T`` to the right by ``zero``s to have
    at least the given ``length``.

    INPUT:

    - ``T`` -- A tuple, list or other iterable.

    - ``length`` -- a nonnegative integer.

    - ``zero`` -- (default: ``0``) the elements to pad with.

    OUTPUT:

    An object of the same type as ``T``.

    EXAMPLES::

        sage: from sage.combinat.k_regular_sequence import pad_right
        sage: pad_right((1,2,3), 10)
        (1, 2, 3, 0, 0, 0, 0, 0, 0, 0)
        sage: pad_right((1,2,3), 2)
        (1, 2, 3)

    TESTS::

        sage: pad_right([1,2,3], 10)
        [1, 2, 3, 0, 0, 0, 0, 0, 0, 0]
    """
    return T + type(T)(zero for _ in xrange(length - len(T)))


def value(D, k):
    r"""
    Return the value of the expansion with digits `D` in base `k`, i.e.

    .. MATH::

        \sum_{0\leq j < \operator{len}D} D[j] k^j.

    INPUT:

    - ``D`` -- a tuple or other iterable.

    - ``k`` -- the base.

    OUTPUT:

    An element in the common parent of the base `k` and of the entries
    of `D`.

    EXAMPLES::

        sage: from sage.combinat.k_regular_sequence import value
        sage: value(42.digits(7), 7)
        42
    """
    return sum(d * k**j for j, d in enumerate(D))


def split_interlace(n, k, p):
    r"""
    Split each digit in the `k`-ary expansion of `n` into `p` parts and
    return the value of the expansion obtained by each of these parts.

    INPUT:

    - ``n`` -- an integer.

    - ``k`` -- an integer specifying the base.

    - ``p`` -- a positive integer specifying in how many parts
      the input ``n`` is split. This has to be a divisor of ``k``.

    OUTPUT:

    A tuple of integers.

    EXAMPLES::

        sage: from sage.combinat.k_regular_sequence import split_interlace
        sage: [(n, split_interlace(n, 4, 2)) for n in srange(20)]
        [(0, (0, 0)), (1, (1, 0)), (2, (0, 1)), (3, (1, 1)),
         (4, (2, 0)), (5, (3, 0)), (6, (2, 1)), (7, (3, 1)),
         (8, (0, 2)), (9, (1, 2)), (10, (0, 3)), (11, (1, 3)),
         (12, (2, 2)), (13, (3, 2)), (14, (2, 3)), (15, (3, 3)),
         (16, (4, 0)), (17, (5, 0)), (18, (4, 1)), (19, (5, 1))]
        sage: [(n, split_interlace(n, 6, 3)) for n in srange(9)]
        [(0, (0, 0, 0)), (1, (1, 0, 0)), (2, (0, 1, 0)),
         (3, (1, 1, 0)), (4, (0, 0, 1)), (5, (1, 0, 1)),
         (6, (2, 0, 0)), (7, (3, 0, 0)), (8, (2, 1, 0))]

    TESTS::

        sage: split_interlace(42, 4, 3)
        Traceback (most recent call last):
        ...
        ValueError: p=3 is not a divisor of k=4.
    """
    if k % p != 0:
        raise ValueError('p={} is not a divisor of k={}.'.format(p, k))
    ki = k // p
    return tuple(value(D, ki)
                 for D in zip(*(d.digits(ki, padto=p)
                                for d in n.digits(k, padto=1))))


from sage.structure.element import Element


class kRegularSequence(Element):

    def __init__(self, parent, matrices, initial=None, selection=None,
                 output_function=None, transpose=False):
        r"""
        A `k`-regular sequence.

        INPUT:

        - ``parent`` -- an instance of :class:`kRegularSequenceSpace`.

        - ``matrices`` -- a tuple or other iterable of square matrices,
          all of which have the same dimension.

        - ``initial`` -- (default: ``None``) a vector.
          When evaluating the sequence, this vector is multiplied
          from the left to the matrix product. If ``None``, then this
          multiplication is skipped.

        - ``selection`` -- (default: ``None``) a vector.
          When evaluating the sequence, this vector is multiplied
          from the left to the matrix product. If ``None``, then this
          multiplication is skipped.

        - ``output_function`` -- (default: ``None``) a function, which is
          applied after evaluating the sequence. This may be used to
          extract the value of a 1x1 matrix.

        - ``transpose`` -- (default: ``False``) a boolean. If set, then
          each of the ``matrices``. Additionally the vectors ``initial``
          and ``selection`` are switched and (if possible)
          transposed as well.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:      vector([0, 1]), vector([1, 0]),
            ....:      transpose=True)
            2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...

        Using an output function::

            sage: Seq2((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:      Matrix([[0, 1]]), Matrix([[1], [0]]),
            ....:      lambda o: o[0, 0], transpose=True)
            2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...
        """
        super(kRegularSequence, self).__init__(parent=parent)

        def tr(M):
            try:
                return M.transpose() if transpose else M
            except AttributeError:
                return M

        self.matrices = tuple(tr(M) for M in matrices)
        self.k = len(self.matrices)
        self.d = self.matrices[0].nrows()
        if not all(M.dimensions() == (self.d, self.d) for M in self.matrices):
            raise ValueError  # TODO

        if not transpose:
            self.initial = initial
            self.selection = selection
        else:
            self.initial = tr(selection)
            self.selection = tr(initial)

        if output_function is None:
            self.output_function = lambda o: o
        else:
            self.output_function = output_function


    def _repr_(self):
        r"""
        Return a representation string of this `k`-regular sequence

        OUTPUT:

        A string

        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: s = Seq2((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:           vector([0, 1]), vector([1, 0]), transpose=True)
            sage: repr(s)  # indirect doctest
            '2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...'
        """
        from sage.misc.lazy_list import lazy_list_formatter
        return lazy_list_formatter(
            [self[n] for n in xrange(11)],  # once slicing works, use self
            name='{}-regular sequence'.format(self.parent().k),
            opening_delimiter='', closing_delimiter='',
            preview=10)


    def info(self):
        r"""
        Displays the matrices of the `k`-linear representation, the initial
        vector and the selection vector.

        OUTPUT:

        Nothing; printing to standard output.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2.guess(lambda n: sum(n.digits(2))).info()
            matrices:
            (
            [1 0]  [ 0 -1]
            [0 1], [ 1  2]
            )
            initial:
            (0, 1)
            selection:
            (1, 0)
        """
        from sys import displayhook
        print('matrices:')
        displayhook(self.matrices)
        print('initial:')
        displayhook(self.initial)
        print('selection:')
        displayhook(self.selection)


    @cached_method
    def __getitem__(self, n):
        r"""
        Return the `n`th entry of this sequence.

        INPUT:

        - ``n`` -- a nonnegative integer.

        OUTPUT:

        An element of the universe of the sequence.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: S = Seq2((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
            ....:          initial=vector([0, 1]), selection=vector([1, 0]))
            sage: S[7]
            3
        """
        result = self._product_of_matrices_(n)
        if self.initial is not None:
            result = self.initial * result
        if self.selection is not None:
            result = result * self.selection
        return self.output_function(result)


    @cached_method
    def _product_of_matrices_(self, m):
        r"""
        Return the product of matrices according to the `k`-ary
        digit expansion of `m`.

        INPUT:

        - ``m`` -- a nonnegative integer.

        OUTPUT:

        A matrix.

        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: S = Seq2((M0, M1))
            sage: S._product_of_matrices_(0) == M0
            True
            sage: S._product_of_matrices_(1) == M1
            True
            sage: S._product_of_matrices_(3) == M1^2
            True

        ::

            sage: S._product_of_matrices_(-1)
            Traceback (most recent call last):
            ...
            ValueError: m=-1 is not a nonnegative integer.
         """
        k = self.parent().k
        if m < 0:
            raise ValueError('m={} is not a nonnegative integer.'.format(m))
        if 0 <= m < k:
            return self.matrices[m]
        n = m // k
        r = m - n*k
        return self.matrices[r] * self._product_of_matrices_(n)


from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

class kRegularSequenceSpace(UniqueRepresentation, Parent):

    Element = kRegularSequence

    def __init__(self, k, base, category=None):
        r"""
        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2
            Space of 2-regular sequences over Integer Ring
            sage: Seq2.category()
            Category of sets
        """
        from sage.categories.sets_cat import Sets
        self.k = k
        super(kRegularSequenceSpace, self).__init__(
            category=category or Sets(), base=base)


    def _repr_(self):
        return 'Space of {}-regular sequences over {}'.format(self.k, self.base())


    def guess(self, f, n_max=None, d_max=None, domain=None, sequence=None):
        r"""

        EXAMPLES:

        Binary sum of digits::

            sage: @cached_function
            ....: def s(n):
            ....:     if n == 0:
            ....:         return 0
            ....:     return s(n//2) + ZZ(is_odd(n))
            sage: all(s(n) == sum(n.digits(2)) for n in srange(10))
            True
            sage: [s(n) for n in srange(10)]
            [0, 1, 1, 2, 1, 2, 2, 3, 1, 2]

        ::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: import logging
            sage: logging.basicConfig(level=logging.INFO)
            sage: S1 = Seq2.guess(s)
            INFO:...:including f_{1*m+0}
            INFO:...:M_0: f_{2*m+0} = (1) * X_m
            INFO:...:including f_{2*m+1}
            INFO:...:M_1: f_{2*m+1} = (0, 1) * X_m
            INFO:...:M_0: f_{4*m+1} = (0, 1) * X_m
            INFO:...:M_1: f_{4*m+3} = (-1, 2) * X_m
            sage: S1.info()
            matrices:
            (
            [1 0]  [ 0 -1]
            [0 1], [ 1  2]
            )
            initial:
            (0, 1)
            selection:
            (1, 0)
            sage: logging.shutdown(); _ = reload(logging)

        ::

            sage: C = Seq2((Matrix([[1]]), Matrix([[1]])), vector([1]), vector([1])); C
            2-regular sequence 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
            sage: S2 = Seq2.guess(s, sequence=C)
            sage: S2.info()
            matrices:
            (
            [1 0]  [1 1]
            [0 1], [0 1]
            )
            initial:
            (1, 0)
            selection:
            (0, 1)

        TESTS::

            sage: Seq2.guess(lambda n: 2, sequence=C).info()
            matrices:
            ([1], [1])
            initial:
            (1)
            selection:
            (2)
        """
        import logging
        logger = logging.getLogger(__name__)

        from sage.arith.srange import srange, xsrange
        from sage.matrix.constructor import Matrix
        from sage.misc.mrange import cantor_product
        from sage.modules.free_module_element import vector

        k = self.k
        if n_max is None:
            n_max = 100
        if d_max is None:
            d_max = 10
        if domain is None:
            domain = self.base()  # TODO
        if sequence is None:
            matrices = [[] for _ in srange(k)]
            class ES(object):
                def __getitem__(self, m):
                    return tuple()
            sequence = ES()
        else:
            matrices = [M.rows() for M in sequence.matrices]
            sequence = sequence.parent()(sequence.matrices,
                                         initial=sequence.initial)

        zero = domain(0)
        one = domain(1)

        def values(m, lines):
            return tuple(sequence[m]) + tuple(f(k**t_R * m + r_R) for t_R, r_R, s_R in lines)

        @cached_function(key=lambda lines: len(lines))  # we assume that existing lines are not changed (we allow appending of new lines)
        def some_inverse_U_matrix(lines):
            d = len(sequence[0]) + len(lines)

            for m_indices in cantor_product(xsrange(n_max), repeat=d, min_slope=1):
                U = Matrix(domain, d, d, [values(m, lines) for m in m_indices]).transpose()
                try:
                    return U.inverse(), m_indices
                except ZeroDivisionError:
                    pass
            else:
                raise RuntimeError

        def guess_linear_dependence(t_L, r_L, lines):
            iU, m_indices = some_inverse_U_matrix(lines)
            X_L = vector(f(k**t_L * m + r_L) for m in m_indices)
            return X_L * iU

        def verify_linear_dependence(t_L, r_L, linear_dependence, lines):
            return all(f(k**t_L * m + r_L) ==
                       linear_dependence * vector(values(m, lines))
                       for m in xsrange(0, (n_max - r_L) // k**t_L + 1))

        def find_linear_dependence(t_L, r_L, lines):
            linear_dependence = guess_linear_dependence(t_L, r_L, lines)
            if not verify_linear_dependence(t_L, r_L, linear_dependence, lines):
                raise ValueError
            return linear_dependence

        selection = None
        if sequence[0]:
            try:
                solution = find_linear_dependence(0, 0, [])
            except ValueError:
                pass
            else:
                selection = vector(solution)

        to_branch = []
        lines = []
        def include(line):
            to_branch.append(line)
            lines.append(line)
            t, r, s = line
            logger.info('including f_{%s*m+%s}', k**t, r)

        if selection is None:
            line_L = (0, 0, 0)  # entries (t, r, s) --> k**t * m + r, belong to M_s
            include(line_L)
            selection = vector((len(sequence[0]) + len(lines)-1)*(zero,) + (one,))

        while to_branch:
            line_R = to_branch.pop(0)
            t_R, r_R, s_R = line_R
            if t_R >= d_max:
                raise RuntimeError

            t_L = t_R + 1
            for s_L in srange(k):
                r_L = k**t_R * s_L + r_R
                line_L = t_L, r_L, s_L

                try:
                    solution = find_linear_dependence(t_L, r_L, lines)
                except ValueError:
                    include(line_L)
                    solution = (len(lines)-1)*(zero,) + (one,)
                logger.info('M_%s: f_{%s*m+%s} = %s * X_m',
                            s_L, k**t_L, r_L, solution)
                matrices[s_L].append(solution)

        d = len(sequence[0]) + len(lines)
        matrices = tuple(Matrix(domain, [pad_right(tuple(row), d, zero=zero) for row in M]).transpose()
                         for M in matrices)
        initial = vector(values(0, lines))
        selection = vector(pad_right(tuple(selection), d, zero=zero))
        return self(matrices, initial, selection)
