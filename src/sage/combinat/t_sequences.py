# sage_setup: distribution = sagemath-combinat
r"""
T-sequences

T-sequences are tuples of four (-1, 0, 1) sequences of length `t` where
for every `i` exactly one sequence has a nonzero entry at index `i`
and for which the nonperiodic autocorrelation function is equal to zero
(i.e. they are complementary). See Definition 7.5 of [Seb2017]_.

These can be constructed from Turyn sequences. In particular,
if Turyn sequences of length `l` exists, there will be T-sequences
of length `4l-1` and `2l-1`.

Turyn sequences are tuples of four (-1, +1) sequences `X, U, Y, V` of length
`l`, `l`, `l-1`, `l-1` with nonperiodic autocorrelation equal to zero and
the additional constraints that:

* the first element of `X` is 1
* the last element of `X` is -1
* the last element of `U` is 1

The nonperiodic autocorrelation of a familiy of sequences
`X=\{A_1, A_2, ..., A_n\}` is defined as (see Definition 7.2 of [Seb2017]_):

.. MATH::

    N_X(j) = \sum_{i=1}^{n-j}(a_{1,i}a_{1,i+j} + a_{2,i}a_{2,i+j} + ... + a_{n,i}a_{n,i+j})

AUTHORS:

- Matteo Cati (2022-11-16): initial version
"""

# ***************************************************************************
#       Copyright (C) 2022 Matteo Cati matteo.cati@keble.ox.ac.uk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sequence import Sequence


def _nonperiodic_autocorrelation(sequences, j):
    r"""
    Compute the nonperiodic autocorrelation of a familiy of sequences.

    Namely, given a family of sequences `X` it computes:

    .. MATH::

        N_X(j) = \sum_{i=1}^{n-j}(a_{1,i}a_{1,i+j} + a_{2,i}a_{2,i+j} + ... + a_{n,i}a_{n,i+j})

    INPUT:

    - ``sequences`` -- either a single sequence or a list of sequences for
      which we want to compute the nonperiodic autocorrelation

    - ``j`` -- integer; the parameter `j` used when calculating the nonperiodic
      autocorrelation
    """
    if not isinstance(sequences[0], list):
        sequences = [sequences]

    t = len(sequences[0])
    result = 0
    for i in range(t - j):
        for seq in sequences:
            result += seq[i] * seq[i + j]
    return result


def is_skew(seq, verbose=False):
    r"""
    Check if the given sequence is skew.

    A sequence `X=\{x_1, x_2, ...,x_n\}` is defined skew (according to Definition
    7.4 of [Seb2017]_) if `n` is even and `x_i = -x_{n-i+1}`.

    INPUT:

    - ``seq`` -- the sequence that should be checked

    - ``verbose`` -- boolean (default: ``False``); if ``True`` the function
      will be verbose when the sequences do not satisfy the contraints

    EXAMPLES::

        sage: from sage.combinat.t_sequences import is_skew
        sage: is_skew([1, -1, 1, -1, 1, -1])
        True
        sage: is_skew([1, -1, -1, -1], verbose=True)
        Constraint not satisfied at index 1
        False

    TESTS::

        sage: is_skew([1, -1, -1])
        False
        sage: is_skew([1, -1, -1, 1, -1], verbose=True)
        Sequence should be of even length
        False
    """
    n = len(seq)

    if n % 2:
        if verbose:
            print('Sequence should be of even length')
        return False

    for i in range(n):
        if seq[i] != -seq[n-i-1]:
            if verbose:
                print(f'Constraint not satisfied at index {i}')
            return False
    return True


def is_symmetric(seq, verbose=False) -> bool:
    r"""
    Check if the given sequence is symmetric.

    A sequence `X=\{x_1, x_2, ...,x_n\}` is defined symmetric (according to Definition
    7.4 of [Seb2017]_) if `n` is odd and `x_i = x_{n-i+1}`.

    INPUT:

    - ``seq`` -- the sequence that should be checked

    - ``verbose`` -- boolean (default: ``False``); if ``True`` the function will be
      verbose when the sequences do not satisfy the contraints

    EXAMPLES::

        sage: from sage.combinat.t_sequences import is_symmetric
        sage: is_symmetric([1, -1, 1, -1, 1])
        True
        sage: is_symmetric([1, -1, 1, 1, 1], verbose=True)
        Constraint not satisfied at index 1
        False

    TESTS::

        sage: is_symmetric([1, -1, -1, 1])
        False
        sage: is_symmetric([1, -1, -1, 1], verbose=True)
        Sequence should be of odd length
        False
    """
    n = len(seq)

    if n % 2 == 0:
        if verbose:
            print('Sequence should be of odd length')
        return False

    for i in range(n):
        if seq[i] != seq[n-i-1]:
            if verbose:
                print(f'Constraint not satisfied at index {i}')
            return False
    return True


def is_T_sequences_set(sequences, verbose=False):
    r"""
    Check if a family of sequences is composed of T-sequences.

    Given 4 (-1, 0, +1) sequences, they will be T-sequences if
    (Definition 7.4 of [Seb2017]_):

    * they have all the same length `t`
    * for each index `i`, exactly one sequence is nonzero at `i`
    * the nonperiodic autocorrelation is equal to `0`

    INPUT:

    - ``sequences`` -- list of four sequences

    - ``verbose`` -- boolean (default: ``False``); if ``True`` the function will be
      verbose when the sequences do not satisfy the contraints

    EXAMPLES::

        sage: from sage.combinat.t_sequences import is_T_sequences_set
        sage: seqs = [[1, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, -1], [0, 0, 0, 0, 0]]
        sage: is_T_sequences_set(seqs)
        True
        sage: seqs = [[1, 1, 0, 1, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, -1], [0, 0, 0, 0, 0]]
        sage: is_T_sequences_set(seqs, verbose=True)
        There should be exactly a nonzero element at every index, found 2 such
         elements at index 3
        False


    TESTS::

        sage: seqs = [[1, 1, 0, 1, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, -1], [0, 0, 0, 0, 0]]
        sage: is_T_sequences_set(seqs)
        False
        sage: seqs = [[1, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, -1, -1], [0, 0, 0, 0, 0]]
        sage: is_T_sequences_set(seqs, verbose=True)
        Nonperiodic autocorrelation should always be zero, found 2 for parameter 1
        False
        sage: is_T_sequences_set([[1, 0, ], [0, -1, 0], [0, 0, 1]])
        False
        sage: seqs = [[1, 2, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, -1, -1], [0, 0, 0, 0, 0]]
        sage: is_T_sequences_set(seqs, verbose=True)
        Elements should be in (-1, 0, +1), but 2 was found at index 1
        False
    """
    if len(sequences) != 4:
        if verbose:
            print(f"T-Sequence should contain 4 sequences, found {len(sequences)} instead")
        return False

    t = len(sequences[0])

    for i in range(t):
        tot = 0
        for seq in sequences:
            if seq[i] not in [-1, 0, 1]:
                if verbose:
                    print(f"Elements should be in (-1, 0, +1), but {seq[i]} was found at index {i}")
                return False
            tot += abs(seq[i])
        if tot != 1:
            if verbose:
                print(f"There should be exactly a nonzero element at every index, found {tot} such elements at index {i}")
            return False

    for j in range(1, t):
        autocorr = _nonperiodic_autocorrelation(sequences, j)
        if autocorr != 0:
            if verbose:
                print(f"Nonperiodic autocorrelation should always be zero, found {autocorr} for parameter {j}")
            return False

    return True


def turyn_sequences_smallcases(l, existence=False):
    r"""
    Construction of Turyn sequences for small values of `l`.

    The data is taken from [Seb2017]_ and [CRSKKY1989]_.

    INPUT:

    - ``l`` -- integer; the length of the Turyn sequences

    - ``existence`` -- boolean (default: ``False``); if ``True``, only return
      whether the Turyn sequences are available for the given length

    EXAMPLES:

    By default, this method returns the four Turyn sequences ::

        sage: from sage.combinat.t_sequences import turyn_sequences_smallcases
        sage: turyn_sequences_smallcases(4)
        [[1, 1, -1, -1], [1, 1, -1, 1], [1, 1, 1], [1, -1, 1]]

    If we pass the ``existence`` flag, the method will return a boolean ::

        sage: turyn_sequences_smallcases(4, existence=True)
        True

    TESTS::

        sage: turyn_sequences_smallcases(17)
        Traceback (most recent call last):
        ...
        ValueError: Turyn sequence of length 17 is not implemented yet.
        sage: turyn_sequences_smallcases(17, existence=True)
        False
    """
    db = {
        2: [[1, -1], [1, 1], [1], [1]],
        3: [[1, 1, 1], [1, 1, -1], [1, -1], [1, -1]],
        4: [[1, 1, -1, -1], [1, 1, -1, 1], [1, 1, 1], [1, -1, 1]],
        5: [[1, 1, -1, 1, 1], [1, 1, 1, 1, -1], [1, 1, -1, -1], [1, -1, 1, -1]],
        6: [[1, 1, 1, -1, -1, -1], [1, 1, -1, 1, -1, 1], [1, 1, -1, 1, 1], [1, 1, -1, 1, 1]],
        7: [[1, 1, 1, -1, 1, 1, 1], [1, 1, -1, -1, -1, 1, -1], [1, 1, -1, 1, -1, -1], [1, 1, -1, 1, -1, -1]],
        8: [[1, 1, -1, 1, -1, 1, -1, -1], [1, 1, 1, 1, -1, -1, -1, 1], [1, 1, 1, -1, 1, 1, 1], [1, -1, -1, 1, -1, -1, 1]],
        13: [[1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 1, 1], [1, 1, 1, -1, -1, 1, -1, 1, -1, -1, 1, 1, -1],
             [1, 1, 1, -1, 1, 1, -1, -1, 1, -1, -1, -1], [1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1]],
        15: [[1, 1, -1, 1, 1, 1, -1, 1, -1, 1, 1, 1, -1, 1, 1], [1, 1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1, 1, -1],
             [1, 1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1, -1], [1, -1, -1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1, -1]],
    }

    if existence:
        return l in db

    if l not in db:
        raise ValueError(f"Turyn sequence of length {l} is not implemented yet.")

    return list(map(Sequence, db[l]))


def T_sequences_construction_from_base_sequences(base_sequences, check=True):
    r"""
    Construct T-sequences of length `2n+p` from base sequences of length `n+p, n+p, n, n`.

    Given base sequences `A, B, C, D`, the T-sequences are constructed as described in
    [KTR2005]_:

    .. MATH::

        \begin{aligned}
        T_1 &= \frac{1}{2}(A+B); 0_{n} \\
        T_2 &= \frac{1}{2}(A-B); 0_{n} \\
        T_3 &= 0_{n+p} + \frac{1}{2}(C+D) \\
        T_4 &= 0_{n+p} + \frac{1}{2}(C-D)
        \end{aligned}

    INPUT:

    - ``base_sequences`` -- the base sequences that should be used to construct
      the T-sequences

    - ``check`` -- boolean (default: ``True``); check that the sequences
      created are T-sequences before returning them

    EXAMPLES::

        sage: from sage.combinat.t_sequences import turyn_sequences_smallcases, T_sequences_construction_from_base_sequences
        sage: seqs = turyn_sequences_smallcases(4)
        sage: T_sequences_construction_from_base_sequences(seqs)
        [[1, 1, -1, 0, 0, 0, 0],
        [0, 0, 0, -1, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 1],
        [0, 0, 0, 0, 0, 1, 0]]

    TESTS::

        sage: from sage.combinat.t_sequences import base_sequences_construction, is_T_sequences_set
        sage: seqs = turyn_sequences_smallcases(4)
        sage: is_T_sequences_set(T_sequences_construction_from_base_sequences(seqs))
        True
        sage: T_sequences_construction_from_base_sequences([[1, -1], [-1, 1], [1]])
        Traceback (most recent call last):
        ...
        AssertionError
        sage: X = [1,1,-1,1,-1,1,-1,1]
        sage: Y = [1,-1,-1,-1,-1,-1,-1,1]
        sage: Z = [1,-1,-1,1,1,1,1,-1]
        sage: W = [1,1,1,-1,1,1,-1]
        sage: base_seqs = base_sequences_construction([X, Y, Z, W])
        sage: is_T_sequences_set(T_sequences_construction_from_base_sequences(base_seqs))
        True
    """

    assert len(base_sequences) == 4

    A, B, C, D = base_sequences
    n = len(C)
    p = len(A)-n

    assert len(A) == len(B) == len(C)+p == len(D)+p

    def seq_sum(seq1, seq2):
        return [(a+b)//2 for (a, b) in zip(seq1, seq2)]

    def seq_subtract(seq1, seq2):
        return [(a-b)//2 for (a, b) in zip(seq1, seq2)]

    def zero_seq(n):
        return [0 for _ in range(n)]

    X1 = Sequence(seq_sum(A, B) + zero_seq(n))
    X2 = Sequence(seq_subtract(A, B) + zero_seq(n))
    X3 = Sequence(zero_seq(n+p) + seq_sum(C, D))
    X4 = Sequence(zero_seq(n+p) + seq_subtract(C, D))

    res = [X1, X2, X3, X4]
    if check:
        assert is_T_sequences_set(res)
    return res


def T_sequences_construction_from_turyn_sequences(turyn_sequences, check=True):
    r"""
    Construct T-sequences of length `4l-1` from Turyn sequences of length `l`.

    Given Turyn sequences `X, U, Y, V`, the T-sequences are constructed as described in
    theorem 7.7 of [Seb2017]_:

    .. MATH::

        \begin{aligned}
        T_1 &= 1; 0_{4l-2} \\
        T_2 &= 0; X/Y; 0_{2l-1} \\
        T_3 &= 0_{2l}; U/0_{l-2} \\
        T_4 &= 0_{2l} + 0_{l}/V
        \end{aligned}

    INPUT:

    - ``turyn_sequences`` -- the Turyn sequences that should be used to
      construct the T-sequences

    - ``check`` -- boolean (default: ``True``); check that the sequences
      created are T-sequences before returning them

    EXAMPLES::

        sage: from sage.combinat.t_sequences import turyn_sequences_smallcases, T_sequences_construction_from_turyn_sequences, is_T_sequences_set
        sage: seqs = turyn_sequences_smallcases(4)
        sage: T_sequences_construction_from_turyn_sequences(seqs)
        [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, -1, 1, -1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -1, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 1, 0]]

    TESTS::

        sage: seqs = turyn_sequences_smallcases(4)
        sage: is_T_sequences_set(T_sequences_construction_from_turyn_sequences(seqs))
        True
        sage: T_sequences_construction_from_turyn_sequences([[1, -1], [-1, 1], [1]])
        Traceback (most recent call last):
        ...
        AssertionError
    """

    assert len(turyn_sequences) == 4

    X, U, Y, V = turyn_sequences
    l = len(X)

    assert len(X) == len(U) == len(Y)+1 == len(V)+1

    def zero_seq(n):
        return [0 for _ in range(n)]

    def interleave(seq1, seq2):
        res = []
        for i in range(len(seq1) + len(seq2)):
            if i % 2 == 0:
                res.append(seq1[i//2])
            else:
                res.append(seq2[i//2])
        return res

    X1 = Sequence([1] + zero_seq(4*l-2))
    X2 = Sequence([0] + interleave(X, Y) + zero_seq(2*l-1))
    X3 = Sequence(zero_seq(2*l) + interleave(U, zero_seq(l-1)))
    X4 = Sequence(zero_seq(2*l) + interleave(zero_seq(l), V))

    res = [X1, X2, X3, X4]
    if check:
        assert is_T_sequences_set(res)
    return res


def T_sequences_smallcases(t, existence=False, check=True):
    r"""
    Construct T-sequences for some small values of `t`.

    This method will try to use the constructions defined in
    :func:`T_sequences_construction_from_base_sequences` and
    :func:`T_sequences_construction_from_turyn_sequences`
    together with the Turyn sequences stored in :func:`turyn_sequences_smallcases`,
    or base sequences created by :func:`base_sequences_smallcases`.

    This function contains also some T-sequences taken directly from [CRSKKY1989]_.

    INPUT:

    - ``t`` -- integer; the length of the T-sequences to construct

    - ``existence`` -- boolean (default: ``False``); if ``True``, this method
      only returns whether a T-sequences of the given size can be constructed

    - ``check`` -- boolean (default: ``True``); check that the sequences are
      T-sequences before returning them

    EXAMPLES:

    By default, this method returns the four T-sequences ::

        sage: from sage.combinat.t_sequences import T_sequences_smallcases, is_T_sequences_set
        sage: T_sequences_smallcases(9)
        [[1, 1, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, -1, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, -1],
        [0, 0, 0, 0, 0, 0, 1, -1, 0]]

    If the existence flag is passed, the method returns a boolean ::

        sage: T_sequences_smallcases(9, existence=True)
        True

    TESTS::

        sage: T_sequences_smallcases(66)
        Traceback (most recent call last):
        ...
        ValueError: T Sequences of length 66 not yet implemented.
        sage: is_T_sequences_set(T_sequences_smallcases(47))
        True
        sage: is_T_sequences_set(T_sequences_smallcases(11))
        True
        sage: T_sequences_smallcases(69, existence=True)
        False
    """
    db = {
        47: [
            [1,-1,-1,0,0,-1,1,-1]+[0]*8+[1,-1,-1,0,0,-1,-1]+[0]*24,
            [0,0,0,-1,1,0,0,0,-1,-1,-1,1,1,1,1,1,0,0,0,1,-1,0,0,1]+[0]*23,
            [0]*26+[-1,0,1,0,0,0,0,1,-1,1,1,1,0,0,0,0,1,0,-1,0,0],
            [0]*24 + [1,1,0,-1,0,-1,1,1,-1,0,0,0,0,0,-1,1,-1,-1,0,-1,0,-1,1]
        ],
        65: [
            [0]*33+[1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,-1,-1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1,1],
            [0]*32+[1]+[0]*32,
            [1]*5+[-1,-1,1,1,-1,1,-1,1,1]+[-1]*7+[1,1,-1,1,-1,1,-1,1,1,-1,-1]+[0]*33,
            [0]*65
        ],
        93: [
            [0,-1,0,0,-1,1,0,-1,1,0,1,1,0,0,1,1,1,0,0,-1,0,-1,1,1,1,-1,0,1,0,0,1]+[0]*33+[1,1,0,0,1,0,0,-1,0,0,-1,1,0,1]+[0]*15,
            [-1,0,-1,1,0,0,1,0,0,-1,0,0,-1,-1,0,0,0,-1,1,0,1]+[0]*5+[-1,0,1,1]+[0]*32+[1,1,0,0,1,1,0,1,-1,0,1,-1,0,0,-1]+[0]*16,
            [0]*32+[1,0,0,1,-1,0,1,-1,0,-1,-1,0,0,-1,-1,1,0,0,-1,0,-1,1,1,1,-1,0,1,0,0,1]+[0]*17+[1,1,0,-1]+[0]*5+[1,0,1,-1,0],
            [0]*31+[1,0,1,-1,0,0,-1,0,0,1,0,0,1,1,0,0,0,-1,1,0,1]+[0]*5+[-1,0,1,1]+[0]*17+[-1,0,0,-1,0,1,-1,-1,-1,1,0,1,0,0,-1]
        ]
    }

    if t in db:
        if existence:
            return True
        sequences = list(map(Sequence, db[t]))
        if check:
            assert is_T_sequences_set(sequences)
        return sequences
    if (t+1) % 2 == 0 and turyn_sequences_smallcases((t+1)//2, existence=True):
        if existence:
            return True
        turyn_seqs = turyn_sequences_smallcases((t+1)//2)
        return T_sequences_construction_from_base_sequences(turyn_seqs, check=check)

    if (t+1) % 4 == 0 and turyn_sequences_smallcases((t+1)//4, existence=True):
        if existence:
            return True
        turyn_seqs = turyn_sequences_smallcases((t+1)//4)
        return T_sequences_construction_from_turyn_sequences(turyn_seqs, check=check)

    for p in range(1, t):
        n = (t-p)//2
        if (t-p) % 2 == 0 and base_sequences_smallcases(n, p, existence=True):
            if existence:
                return True
            base_seqs = base_sequences_smallcases(n, p, check=False)
            return T_sequences_construction_from_base_sequences(base_seqs, check=check)

    if existence:
        return False
    raise ValueError(f'T Sequences of length {t} not yet implemented.')


def base_sequences_construction(turyn_type_seqs, check=True):
    r"""
    Construct base sequences of length `2n-1, 2n-1, n, n` from Turyn type
    sequences of length `n,n,n,n-1`.

    Given Turyn type sequences `X, Y, Z, W` of length `n,n,n,n-1`, Theorem 1 of
    [KTR2005]_  shows that the following are base sequences of length
    `2n-1, 2n-1, n, n`:

    .. MATH::

        \begin{aligned}
        A &= Z;W \\
        B &= Z; -W \\
        C &= X \\
        D &= Y
        \end{aligned}

    INPUT:

    - ``turyn_type_seqs`` -- the list of 4 Turyn type sequences that should be
      used to construct the base sequences

    - ``check`` -- boolean (default: ``True``); check that the resulting
      sequences are base sequences before returning them

    OUTPUT: list containing the four base sequences

    EXAMPLES::

        sage: from sage.combinat.t_sequences import base_sequences_construction
        sage: X = [1,1,-1,1,-1,1,-1,1]
        sage: Y = [1,-1,-1,-1,-1,-1,-1,1]
        sage: Z = [1,-1,-1,1,1,1,1,-1]
        sage: W = [1,1,1,-1,1,1,-1]
        sage: base_sequences_construction([X, Y, Z, W])
        [[1, -1, -1, 1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, -1],
        [1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, -1, -1, 1],
        [1, 1, -1, 1, -1, 1, -1, 1],
        [1, -1, -1, -1, -1, -1, -1, 1]]

    TESTS::

        sage: base_sequences_construction([[1, -1], [1], [1], [-1]])
        Traceback (most recent call last):
        ...
        AssertionError

    .. SEEALSO::

        :func:`is_base_sequences_tuple`
    """
    assert len(turyn_type_seqs) == 4
    X, Y, Z, W = turyn_type_seqs

    assert len(X) == len(Y) == len(Z) == len(W)+1

    A = Sequence(Z + W)
    B = Sequence(Z + [-el for el in W])
    C = X
    D = Y

    if check:
        assert is_base_sequences_tuple([A, B, C, D])
    return [A, B, C, D]


def is_base_sequences_tuple(base_sequences, verbose=False):
    r"""Check if the given sequences are base sequences.

    Four (-1, +1) sequences `A, B, C, D` of length `n+p, n+p, n, n` are called
    base sequences if for all `j \ge 1`:

    .. MATH::

        N_A(j)+N_B(j)+N_C(j)+N_D(j) = 0

    where `N_X(j)` is the nonperiodic autocorrelation (See definition in [KTR2005]_).

    INPUT:

    - ``base_sequences`` -- the list of 4 sequences that should be checked

    - ``verbose`` -- boolean (default: ``False``); if ``True`` the function
      will be verbose when the sequences do not satisfy the contraints

    EXAMPLES::

        sage: from sage.combinat.t_sequences import is_base_sequences_tuple
        sage: seqs = [[1, -1, -1, 1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, -1],[1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, -1, -1, 1],[1, 1, -1, 1, -1, 1, -1, 1],[1, -1, -1, -1, -1, -1, -1, 1]]
        sage: is_base_sequences_tuple(seqs)
        True

    If verbose is true, the function will be verbose ::

        sage: seqs = [[1, -1], [1, 1], [-1], [2]]
        sage: is_base_sequences_tuple(seqs, verbose=True)
        Base sequences should only contiain -1, +1, found 2
        False

    TESTS:

        sage: seqs = [[1, -1], [1], [-1]]
        sage: is_base_sequences_tuple(seqs)
        False
        sage: seqs = [[1, -1], [1, -1], [-1], [1]]
        sage: is_base_sequences_tuple(seqs)
        False
        sage: seqs = [[1, -1], [1, 1], [-1], [2]]
        sage: is_base_sequences_tuple(seqs)
        False
        sage: seqs = [[1, -1], [1], [-1], [1]]
        sage: is_base_sequences_tuple(seqs)
        False

    .. SEEALSO::

        :func:`base_sequences_construction`
    """
    if len(base_sequences) != 4:
        if verbose:
            print(f'Base sequences should be 4, found {len(base_sequences)}')
        return False
    A, B, C, D = base_sequences
    n = len(C)
    p = len(A) - len(C)
    if not (len(A) == len(B) == len(C)+p == len(D)+p):
        if verbose:
            print(f'Base sequences should have length n+p, n+p, n, n, found {len(A)}, {len(B)}, {len(C)}, {len(D)}')
        return False

    for seq in base_sequences:
        for el in seq:
            if abs(el) != 1:
                if verbose:
                    print(f'Base sequences should only contiain -1, +1, found {el}')
                return False

    for j in range(1, n+p):
        autocorr = _nonperiodic_autocorrelation(A, j) + _nonperiodic_autocorrelation(B, j) + _nonperiodic_autocorrelation(C, j) + _nonperiodic_autocorrelation(D, j)
        if autocorr != 0:
            if verbose:
                print(f"Nonperiodic autocorrelation should always be zero, found {autocorr} for parameter {j}")
            return False

    return True


def turyn_type_sequences_smallcases(n, existence=False):
    r"""
    Construction of Turyn type sequences for small values of `n`.

    The data is taken from [KTR2005]_ for `n= 36`, and from [BDKR2013]_ for `n\le 32`.

    INPUT:

    - ``n`` -- integer; the length of the Turyn type sequences

    - ``existence`` -- boolean (default: ``False``); if ``True``, only return
      whether the Turyn type sequences are available for the given length

    EXAMPLES:

    By default, this method returns the four Turyn type sequences ::

        sage: from sage.combinat.t_sequences import turyn_type_sequences_smallcases
        sage: turyn_type_sequences_smallcases(4)
        [[1, 1, 1, 1], [1, 1, -1, 1], [1, 1, -1, -1], [1, -1, 1]]

    If we pass the ``existence`` flag, the method will return a boolean ::

        sage: turyn_type_sequences_smallcases(4, existence=True)
        True

    TESTS::

        sage: turyn_type_sequences_smallcases(17)
        Traceback (most recent call last):
        ...
        ValueError: Turyn type sequences of length 17 are not implemented yet.
        sage: turyn_type_sequences_smallcases(17, existence=True)
        False

    ALGORITHM:

    The Turyn type sequences are stored in hexadecimal format.
    Given `n` hexadecimal digits `h_1, h_2,...,h_n`, it is possible to get the
    Turyn type sequences by converting each `h_i` (`1 \le i \le n-1`) into a
    four digits binary number. Then, the `j`-th binary digit is `0` if the `i`-th
    number in the `j`-th sequence is `1`, and it is `1` if the number in the
    sequence is -1.

    For the `n`-th digit, it should be converted to a 3 digits binary number, and
    then the same mapping as before can be used (see also [BDKR2013]_).
    """
    def convertLists(hexstring):
        seqs = [Sequence([]), Sequence([]), Sequence([]), Sequence([])]
        for c in hexstring[:-1]:
            binary = bin(int(c, 16))[2:].zfill(4)
            for i in range(4):
                if binary[i] == '0':
                    seqs[i].append(1)
                else:
                    seqs[i].append(-1)
        last = bin(int(hexstring[-1], 16))[2:].zfill(3)
        for i in range(3):
            if last[i] == '0':
                seqs[i].append(1)
            else:
                seqs[i].append(-1)
        return seqs

    db = {
        2: '01',
        4: '0161',
        6: '006d61',
        8: '06e5c4d1',
        10: '0001f4a961',
        12: '0004f90bc961',
        14: '00036ac71c7651',
        16: '0000778e52de5561',
        18: '00006758b30d1e9a51',
        20: '000038e2739c7a0b6951',
        22: '00000f702c71a9ad565961',
        24: '00000b7c2cb2bc4b6cd9a961',
        26: '000000ff0f846f1ca5a5aa9551',
        28: '0000067cde3e50639ab46135aa51',
        30: '000000f70b106f9d427a25e9a96951',
        32: '00000138f64f1c1e77844f26d95a5961',
        36: '060989975b685d8fc80750b21c0212eceb26',
    }

    if existence:
        return n in db

    if n not in db:
        raise ValueError(f"Turyn type sequences of length {n} are not implemented yet.")

    return convertLists(db[n])


def base_sequences_smallcases(n, p, existence=False, check=True):
    r"""Construct base sequences of length `n+p, n+p, n, n` from available data.

    The function uses the construction :func:`base_sequences_construction`,
    together with Turyn type sequences from :func:`turyn_type_sequences_smallcases`
    to construct base sequences with `p = n-1`.

    Furthermore, this function uses also Turyn sequences (i.e. base sequences
    with `p=1`) from :func:`turyn_sequences_smallcases`.

    INPUT:

    - ``n`` -- integer; the length of the last two base sequences

    - ``p`` -- integer; `n+p` will be the length of the first two base
      sequences

    - ``existence`` -- boolean (default: ``False``); if ``True``, the function
      will only check whether the base sequences can be constructed

    - ``check`` -- boolean (default: ``True``); check that the resulting
      sequences are base sequences before returning them

    OUTPUT:

    If ``existence`` is ``False``, the function returns a list containing the
    four base sequences, or raises an error if the base sequences cannot be
    constructed. If ``existence`` is ``True``, the function returns a boolean,
    which is ``True`` if the base sequences can be constructed and ``False``
    otherwise.

    EXAMPLES::

        sage: from sage.combinat.t_sequences import base_sequences_smallcases
        sage: base_sequences_smallcases(8, 7)
        [[1, -1, -1, 1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, -1],
        [1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, -1, -1, 1],
        [1, 1, -1, 1, -1, 1, -1, 1],
        [1, -1, -1, -1, -1, -1, -1, 1]]

    If ``existence`` is ``True``, the function returns a boolean ::

        sage: base_sequences_smallcases(8, 7, existence=True)
        True
        sage: base_sequences_smallcases(7, 5, existence=True)
        False

    TESTS::

        sage: base_sequences_smallcases(7, 5)
        Traceback (most recent call last):
        ...
        ValueError: Base sequences of order 12, 12, 7, 7 not yet implemented.
        sage: seqs = base_sequences_smallcases(16, 15)
        sage: len(seqs[0]) == len(seqs[1]) == 16+15
        True
        sage: len(seqs[2]) == len(seqs[3]) == 16
        True
    """

    if existence:
        return p == n-1 and turyn_type_sequences_smallcases(n, existence=True)

    if p == n-1 and turyn_type_sequences_smallcases(n, existence=True):
        if existence:
            return True
        turyn_type_seqs = turyn_type_sequences_smallcases(n)
        return base_sequences_construction(turyn_type_seqs, check=check)
    if p == 1 and turyn_sequences_smallcases(n + p, existence=True):
        if existence:
            return True
        return turyn_sequences_smallcases(n + p)

    raise ValueError(f'Base sequences of order {n+p}, {n+p}, {n}, {n} not yet implemented.')
