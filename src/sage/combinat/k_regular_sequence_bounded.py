r"""
Boundedness of `k`-Regular Sequences

This module contains a collection of algorithms to check for boundedness.
This is done
- based on eigenvalues and
- by the criterion presented in [MS1977].

Various
=======

.. SEEALSO::

    :mod:`sage.combinat.k_regular_sequence`,
    :mod:`recognizable series <sage.combinat.recognizable_series>`,
    :mod:`sage.rings.cfinite_sequence`,
    :mod:`sage.combinat.binary_recurrence_sequences`.

REFERENCES:

.. [MS1977] Arnaldo Mandel, Imre Simon,
   *On Finite Semigroups of Matrices*,
   Theoretical Computer Science 5 (101-111),
   North-Holland Publishing Company

AUTHORS:

- Gabriel Lipnik (2017)

ACKNOWLEDGEMENT:

- Gabriel Lipnik is supported by the
  Austrian Science Fund (FWF): P 24644-N26.
"""
#*****************************************************************************
#       Copyright (C) 2017 Gabriel Lipnik <devel@gabriellipnik.at>
#
# This program is free software: You can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


def multiply_reduce(A, B):
    r"""
    Return the matrix `A\cdot B` with entries `\min{(A\cdot B)_{ij},2}`.

    This function is used in
    :func:`<sage.combinat.k_regular_sequence_bounded.mandel_simon_algorithm>`.

    INPUT:

    - ``A`` -- an `m \times n` matrix
    - ``B`` -- an `n \times p` matrix

    OUTPUT:

    An `m \times p` matrix with entries `\min{(A\cdot B)_{ij},2}`.

    EXAMPLES::

        sage: from sage.combinat.k_regular_sequence_bounded import multiply_reduce
        sage: A = Matrix([[2, 0], [0, 2]])
        sage: B = Matrix([[-2, 0], [0, 2]])
        sage: A*B
        [-4  0]
        [ 0  4]
        sage: multiply_reduce(A, B)
        [-4  0]
        [ 0  2]

    ::

        sage: A = Matrix([[1, 2, 3], [-1, -2, -3], [1, 2, 3]])
        sage: B = Matrix([[1, 2, 3], [2, 3, 4], [1, 2, 3]])
        sage: A*B
        [  8  14  20]
        [ -8 -14 -20]
        [  8  14  20]
        sage: multiply_reduce(A, B)
        [  2   2   2]
        [ -8 -14 -20]
        [  2   2   2]
    """
    return (A*B).apply_map(lambda m: min(m, 2))


def construct_phi(matrices):
    r"""
    Return the set `\phi(S)` as defined in [MS1977].

    INPUT:

    - ``matrices`` -- a list of non-negative square matrices
      in the same dimension

    OUTPUT:

    A list of matrices.

    EXAMPLES::

        sage: from sage.combinat.k_regular_sequence_bounded import construct_phi
        sage: L = [Matrix([[2, 2], [1, 3]]), Matrix([[0, 2], [1, 1]])]
        sage: construct_phi(L)
        [
        [2 2]  [0 2]  [2 2]
        [1 2], [1 1], [2 2]
        ]

    ::

        sage: L = [Matrix([[42, 1, 0], [5, 0, 1], [0, 0, 1]]), Matrix([[0, 1, 1],
        ....: [4, 1, 1], [1, 2, 2]]), Matrix([[5, 1, 1], [1, 7, 1], [0, 1, 32]])]
        sage: construct_phi(L)
        [
        [2 1 0]  [0 1 1]  [2 1 1]  [2 2 1]  [2 2 2]  [2 2 2]  [2 2 2]  [2 2 2]
        [2 0 1]  [2 1 1]  [1 2 1]  [2 2 1]  [1 2 2]  [2 2 2]  [2 2 2]  [2 2 2]
        [0 0 1], [1 2 2], [0 1 2], [0 0 1], [1 2 2], [0 1 2], [0 0 1], [1 2 2],

        [2 0 2]  [2 2 2]  [1 2 2]  [2 2 2]  [2 2 2]  [2 2 2]  [2 2 2]  [2 1 2]
        [2 2 2]  [2 2 2]  [2 2 2]  [2 1 2]  [2 2 2]  [2 1 2]  [2 2 2]  [2 2 2]
        [2 1 2], [2 2 2], [2 2 2], [2 0 2], [2 1 2], [2 1 2], [2 0 2], [2 2 2]
        ]

    Tests::

        sage: L = [Matrix([[20, 1, 0], [2, 0, 0], [117, 0, 8]]),
        ....: Matrix([[0, 2, 1], [1, 0, 0], [1,1,2]]), Matrix([[8, 1, 0],
        ....: [0, 0, 3], [0, 1, 0]])]
        sage: construct_phi(L)
        [
        [2 1 0]  [0 2 1]  [2 1 0]  [2 2 0]  [1 2 2]  [2 2 2]  [2 2 2]  [2 0 2]
        [2 0 0]  [1 0 0]  [0 0 2]  [2 2 0]  [0 2 2]  [2 2 0]  [2 2 2]  [2 1 0]
        [2 0 2], [1 1 2], [0 1 0], [2 2 2], [2 2 2], [2 2 0], [2 2 2], [2 1 2],

        [2 1 2]  [0 1 2]  [2 2 2]  [2 2 2]  [2 2 0]  [2 2 2]  [2 2 2]  [2 2 2]
        [0 2 1]  [2 1 0]  [2 2 0]  [1 2 2]  [2 2 2]  [2 0 2]  [2 1 2]  [0 1 2]
        [2 2 2], [2 2 2], [2 2 2], [2 2 2], [2 2 2], [2 2 2], [2 2 2], [2 2 2],

        [2 2 0]  [1 2 2]  [2 2 2]  [2 2 0]  [2 2 2]  [2 1 2]  [2 2 2]  [2 2 2]
        [2 0 2]  [2 2 2]  [0 2 0]  [2 2 2]  [2 2 2]  [2 2 2]  [2 2 2]  [2 2 2]
        [2 0 0], [1 0 0], [0 0 2], [2 2 0], [0 2 2], [2 1 0], [0 2 1], [2 1 0],

        [2 2 2]  [2 2 2]  [2 2 2]  [2 2 2]  [2 2 2]  [2 2 2]  [2 2 2]  [2 2 2]
        [2 2 2]  [2 2 2]  [2 2 2]  [2 2 2]  [2 2 2]  [2 0 0]  [2 0 0]  [0 0 2]
        [2 2 0], [1 2 2], [2 0 2], [2 1 2], [0 1 2], [2 0 2], [2 2 2], [0 2 0],

        [2 2 2]  [2 2 2]  [2 2 2]  [2 1 2]  [2 0 2]  [0 2 2]  [2 0 2]
        [0 2 2]  [2 0 2]  [2 2 2]  [2 0 2]  [2 2 0]  [2 2 2]  [2 2 2]
        [2 2 2], [2 0 0], [2 0 0], [2 2 2], [2 2 2], [2 2 2], [2 2 2]
        ]
    """
    from sage.arith.srange import srange
    length = len(matrices)
    phi = list(M.apply_map(lambda m: min(m, 2)) for M in matrices)
    for counter in range(1000000):
        for A in matrices:
            for B in phi:
                prod = multiply_reduce(A, B)
                if prod not in phi:
                    phi.append(prod)
        if len(phi) == length:
            return phi
        else:
            length = len(phi)
    raise RuntimeError('Phi too large.')


def is_integer_valued(matrices):
    r"""
    Return whether every matrix in ``matrices`` is integer-valued.

    INPUT:

    - ``matrices`` -- a list of square matrices in the same dimension

    OUTPUT:

    A boolean.

    EXAMPLES::

        sage: from sage.combinat.k_regular_sequence_bounded import is_integer_valued
        sage: matrices = [Matrix([[1, 2], [-1, 0]]), Matrix([[42, -42], [0, 0]])]
        sage: is_integer_valued(matrices)
        True

    ::

        sage: matrices = [Matrix([[1, pi], [-1, 0]])]
        sage: is_integer_valued(matrices)
        False

    ::

        sage: matrices = [Matrix([[1, 1/2], [2/4, 0]])]
        sage: is_integer_valued(matrices)
        False

    ::

        sage: matrices = [Matrix([[1, 4/2], [-1, 0]])]
        sage: is_integer_valued(matrices)
        True
    """
    from sage.matrix.matrix_space import MatrixSpace
    from sage.rings.integer_ring import ZZ
    M = MatrixSpace(ZZ, matrices[0].nrows(), matrices[0].ncols())
    return all(mat in M for mat in matrices)


def is_non_negative(matrices):
    r"""
    Return whether every matrix in ``matrices`` is non-negative.

    INPUT:

    - ``matrices`` -- a list of square matrices in the same dimension

    OUTPUT:

    A boolean.

    EXAMPLES::

        sage: from sage.combinat.k_regular_sequence_bounded import is_non_negative
        sage: matrices = [Matrix([[1, 2], [1, 0]]), Matrix([[42, -42], [0, 0]])]
        sage: is_non_negative(matrices)
        False

    ::

        sage: matrices = [Matrix([[0]])]
        sage: is_non_negative(matrices)
        True

    ::

        sage: matrices = [Matrix([[1, 1/2], [2/4, 0]])]
        sage: is_non_negative(matrices)
        True
    """
    return all(min(mat.list()) >= 0 for mat in matrices)


def is_bounded_via_mandel_simon_algorithm(matrices):
    r"""
    Return whether the semigroup generated whether the semigroup of all
    possible products of ``matrices`` is finite/bounded,
    using a criterion in [MS1977].

    INPUT:

    - ``matrices`` -- a list of non-negative square matrices
    in the same dimension

    OUTPUT:

    A boolean.

    EXAMPLES::

        sage: from sage.combinat.k_regular_sequence_bounded import is_bounded_via_mandel_simon_algorithm
        sage: J = [Matrix([[1, 0, 1], [0, 1, 1], [0, 0, 0]])]
        sage: is_bounded_via_mandel_simon_algorithm(J)
        True

    ::

        sage: from sage.combinat.k_regular_sequence_bounded import is_bounded_via_mandel_simon_algorithm
        sage: K = [Matrix([[1, 1], [1, 1]])]
        sage: is_bounded_via_mandel_simon_algorithm(K)
        False

    ::

        sage: L = [Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 0]])]
        sage: is_bounded_via_mandel_simon_algorithm(L)
        True

    ::

        sage: M = [Matrix([[1, 0], [0, 2]]), Matrix([[1, 0], [0, 0]])]
        sage: is_bounded_via_mandel_simon_algorithm(M)
        False

    Non-integer-valued Input::

        sage: N = [Matrix([[0.5, 0], [1, 0]])]
        sage: is_bounded_via_mandel_simon_algorithm(N)
        Traceback (most recent call last):
        ...
        ValueError: Not all matrices are integer-valued.
    """
    if not is_integer_valued(matrices):
        raise ValueError('Not all matrices are integer-valued.')

    phi = construct_phi(matrices)
    return not any(multiply_reduce(M, M) == M and not M**2 == M**3
                   for M in phi)


def has_bounded_matrix_powers(matrices):
    r"""
    Return whether `M^n` is bounded for `n \to \infty`
    for all `M` in ``matrices``.

    INPUT:

    - ``matrices`` -- a list of square matrices

    ALGORITHM:

    Eigenvalues are used for the check.

    OUTPUT:

    A boolean.

    TEST::

        max(abs(eigenvalues)) = 1, multiplicities okay:
        sage: from sage.combinat.k_regular_sequence_bounded import has_bounded_matrix_powers
        sage: matrices = [Matrix([[-1, 1, 1], [-1, 1, 1], [1, -1, 1]]),
        ....:             Matrix([[-1, 1, 1], [-1, 0, 0], [1, 1, 1]])]
        sage: has_bounded_matrix_powers(matrices)
        True

    ::

        max(abs(eigenvalues)) > 1:
        sage: matrices = [Matrix([[1, 1], [1/2, -1]])]
        sage: has_bounded_matrix_powers(matrices)
        False

    ::

        max(abs(eigenvalues)) = 1, multiplicities different:
        sage: matrices = [Matrix([[0,1],[1,1]])]
        sage: has_bounded_matrix_powers(matrices)
        False

    ::

        max(abs(eigenvalues)) < 1:
        sage: matrices = [Matrix([[1, -1], [1/2, -1]])]
        sage: has_bounded_matrix_powers(matrices)
        True
    """
    from sage.matrix.constructor import Matrix
    from sage.arith.srange import srange

    return all(abs(eVn[0]) < 1 or
                (abs(eVn[0]) == 1 and len(eVn[1]) == eVn[2])
                for mat in matrices
                for eVn in mat.eigenvectors_right())


def make_positive(matrices):
    r"""
    Return a list of non-negative matrices

    INPUT:

    - ``matrices`` -- a list of matrices where every matrix is either
      non-negative or non-positive

    OUTPUT:

    A list of matrices containing every non-negative matrix of ``matrices``,
    and `-M` if `M` is a non-positive matrix of ``matrices``.

    EXAMPLES::

        sage: from sage.combinat.k_regular_sequence_bounded import make_positive
        sage: matrices = [Matrix([[1, 2], [1, 0]]), Matrix([[42, 42], [0, 0]])]
        sage: make_positive(matrices)
        [
        [1 2]  [42 42]
        [1 0], [ 0  0]
        ]

    ::

        sage: matrices = [Matrix([[1, 2], [1, 0]]), Matrix([[-42, -42], [0, 0]])]
        sage: make_positive(matrices)
        [
        [1 2]  [42 42]
        [1 0], [ 0  0]
        ]

    ::

        sage: matrices = [Matrix([[1, 2], [1, 0]]), Matrix([[42, -42], [0, 0]])]
        sage: make_positive(matrices)
        Traceback (most recent call last):
        ...
        ValueError: There is a matrix which is neither non-negative nor non-positive.
    """
    from sage.arith.srange import srange
    def do(mat):
        if is_non_negative(mat):
            return mat
        elif is_non_negative(-mat):
            return -mat
        else:
            raise ValueError('There is a matrix which is neither non-negative nor non-positive.')

    return list(do(mat) for mat in matrices)


def k_regular_sequence_is_bounded(S):
    r"""
    Return whether this `k`-regular sequence is bounded.

    INPUT:

    - ``S`` -- a `k`-regular sequence

    OUTPUT:

    A boolean.

    EXAMPLES:

    Thue-Morse Sequence::

        sage: from sage.combinat.k_regular_sequence_bounded import k_regular_sequence_is_bounded
        sage: Seq2 = kRegularSequenceSpace(2, ZZ)
        doctest:...: FutureWarning: This class/method/function is
        marked as experimental. It, its functionality or its interface
        might change without a formal deprecation.
        See http://trac.sagemath.org/21202 for details.
        sage: TM = Seq2([Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [1, 0]])],
        ....: left=vector([1, 0]), right=vector([0, 1]))
        sage: k_regular_sequence_is_bounded(TM)
        True

    Binary Sum of Digits::

        sage: SD = Seq2([Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])],
        ....: left=vector([0, 1]), right=vector([1, 0]))
        sage: k_regular_sequence_is_bounded(SD)
        False

    Sequence of All Natural Numbers::

        sage: N = Seq2([Matrix([[2, 0], [2, 1]]), Matrix([[0, 1], [-2, 3]])],
        ....: left=vector([1, 0]), right=vector([0, 1]))
        sage: k_regular_sequence_is_bounded(N)
        False

    Indicator Function of Even Integers::

        sage: E = Seq2([Matrix([[0, 1], [0, 1]]), Matrix([[0, 0], [0, 1]])],
        ....: left=vector([1, 0]), right=vector([1, 1]))
        sage: k_regular_sequence_is_bounded(E)
        True

    Indicator Function of Odd Integers::

        sage: O = Seq2([Matrix([[0, 0], [0, 1]]), Matrix([[0, 1], [0, 1]])],
        ....: left=vector([1, 0]), right=vector([0, 1]))
        sage: k_regular_sequence_is_bounded(O)
        True

    Number of Odd Entries in Pascal's Triangle::

        sage: U = Seq2([Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])],
        ....: left=vector([0, 1]), right=vector([1, 0]), transpose=True)
        sage: k_regular_sequence_is_bounded(U)
        False

    Counting '10' in the Binary Representation::

        sage: C = Seq2([Matrix([[0, 1, 0, 0], [0, 0, 0, 1], [-1, 0, 1, 1],
        ....: [0, 0, 0, 1]]), Matrix([[0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0],
        ....: [-1, 0, 1, 1]])], left=vector([1, 0, 0, 0]),
        ....: right=vector([0, 0, 1, 0]))
        sage: k_regular_sequence_is_bounded(C)
        False

    Numbers Starting with '10'::

        sage: D = Seq2([Matrix([[0, 1, 0, 0], [0, 0, 1, 0], [0, -2, 3, 0],
        ....: [0, -2, 2, 1]]), Matrix([[2, 0, 0, 0], [0, 0, 0, 1],
        ....: [0, 2, 0, 1], [0, -2, 0, 3]])], left=vector([1, 0, 0, 0]),
        ....: right=vector([2, 2, 2, 5]))
        sage: k_regular_sequence_is_bounded(D)
        False

    Signum Function::

        sage: S = Seq2([Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 1]])],
        ....: left=vector([1, 0]), right=vector([0, 1]))
        sage: k_regular_sequence_is_bounded(S)
        True

    Number of Digits from the Right to the First '1'::

        sage: S = Seq2([Matrix([[0, 1, 0], [-1, 2, 0], [0, 0, 1]]),
        ....: Matrix([[0, 0, 1], [0, 0, 2], [0, 0, 1]])],
        ....: left=vector([1, 0, 0]), right=vector([0, 0, 1]))
        sage: k_regular_sequence_is_bounded(S)
        Traceback (most recent call last):
        ...
        RuntimeError: It is not decidable with this implementation whether the
        sequence is bounded or not.

    Find another Representation or a Sequence which is equivalent in Boundedness::

        sage: S1 = S.subsequence(a = 1, b = 1)
        sage: k_regular_sequence_is_bounded(S1)
        False

    TESTS::

        sage: S = Seq2((Matrix([[0, 1, 0], [0, 0, 1], [-1, 2, 0]]),
        ....: Matrix([[-1, 0, 0], [-3/4, -1/4, 3/4], [-1/4, 1/4, -3/4]])),
        ....: left=vector([1, 0, 0]), right=vector([-4, -4, -4]))
        sage: k_regular_sequence_is_bounded(S)
        False

    ::

        sage: S = Seq2((Matrix([[1,0], [1,0]]), Matrix([[0,1],[1,0]])),
        ....:          left = vector([1,1]), right = vector([1,0]))
        No handlers could be found for logger "sage.combinat.k_regular_sequence"
        sage: k_regular_sequence_is_bounded(S)
        True
    """
    from sage.arith.srange import srange

    matrices = list(S.mu)
    length = len(matrices)
    try:
        if is_bounded_via_mandel_simon_algorithm(make_positive(matrices)):
            return True
    except ValueError:
        pass

    matrices = list(S.minimized().mu)
    matricesWithout = matrices[1:]
    if not has_bounded_matrix_powers(matricesWithout):
        return False

    try:
        if is_bounded_via_mandel_simon_algorithm(make_positive(matrices)):
            return True
    except ValueError:
        pass

    matricesProd = list(ell*em for ell in matrices for em in matrices[1:]
                        if ell != em)
    if not has_bounded_matrix_powers(matricesProd):
        return False

    raise RuntimeError('It is not decidable with this implementation ' +
                       'whether the sequence is bounded or not.')
