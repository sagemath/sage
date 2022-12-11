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

The nonperiodic autocorrelation of a familiy of sequences `X=\{A_1, A_2, ..., A_n\}` is defined as
(see Definition 7.2 of [Seb2017]_):

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


def  _nonperiodic_autocorrelation(sequences, j):
    r"""
    Compute the nonperiodic autocorrelation of a familiy of sequences.

    Namely, given a family of sequences `X` it computes:

    .. MATH::

        N_X(j) = \sum_{i=1}^{n-j}(a_{1,i}a_{1,i+j} + a_{2,i}a_{2,i+j} + ... + a_{n,i}a_{n,i+j})

    INPUT:
    
    - ``sequences`` -- either a single sequence or a list of sequences for which we want
      to compute the nonperiodic autocorrelation.
    
    - ``j`` -- integer, the parameter `j` used when calculating the nonperiodic autocorrelation.
    """
    if not isinstance(sequences[0], list):
        sequences = [sequences]
    
    t = len(sequences[0])
    result = 0
    for i in range(t-j):
        for seq in sequences:
            result += seq[i]*seq[i+j]
    return result

def is_skew(seq, verbose=False):
    r"""
    Check if the given sequence is skew.

    A sequence `X=\{x_1, x_2, ...,x_n\}` is defined skew (according to Definition 
    7.4 of [Seb2017]_) if `n` is even and `x_i = -x_{n-i+1}`.

    INPUT:

    - ``seq`` -- the sequence that should be checked.

    - ``verbose`` -- a boolean (default false). If true the function will be verbose
      when the sequences do not satisfy the contraints.

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

    if n%2 == 1:
        if verbose:
            print('Sequence should be of even length')
        return False

    for i in range(n):
        if seq[i] != -seq[n-i-1]:
            if verbose:
                print(f'Constraint not satisfied at index {i}')
            return False
    return True

def is_symmetric(seq, verbose=False):
    r"""
    Check if the given sequence is symmetric.

    A sequence `X=\{x_1, x_2, ...,x_n\}` is defined symmetric (according to Definition 
    7.4 of [Seb2017]_) if `n` is odd and `x_i = x_{n-i+1}`.

    INPUT:
    
    - ``seq`` -- the sequence that should be checked.

    - ``verbose`` -- a boolean (default false). If true the function will be verbose
      when the sequences do not satisfy the contraints.

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

    if n%2 == 0:
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
    
    - ``sequences`` -- a list of four sequences.

    - ``verbose`` -- a boolean (default false). If true the function will be verbose
      when the sequences do not satisfy the contraints.

    EXAMPLES::

        sage: from sage.combinat.t_sequences import is_T_sequences_set
        sage: seqs = [[1, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, -1], [0, 0, 0, 0, 0]]
        sage: is_T_sequences_set(seqs)
        True
        sage: seqs = [[1, 1, 0, 1, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, -1], [0, 0, 0, 0, 0]]
        sage: is_T_sequences_set(seqs, verbose=True)
        There should  be exactly a nonzero element at every index, found 2 such elemnents at index 3
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
                print(f"There should  be exactly a nonzero element at every index, found {tot} such elemnents at index {i}")
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

    - ``l`` -- integer, the length of the Turyn sequences.

    - ``existence`` -- boolean (default False). If true, only return whether the 
      Turyn sequences are available for the given length.

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

def T_sequences_construction_I(turyn_sequences, check=True):
    r"""
    Construct T-sequences of length `2l-1` from Turyn sequences of length `l`.

    Given Turyn sequences `X, U, Y, V`, the T-sequences are constructed as described in
    Theorem 7.7 of [Seb2017]_:

    .. MATH::

        \begin{aligned}
        T_1 &= \frac{1}{2}(X+U); 0_{l-1} \\
        T_2 &= \frac{1}{2}(X-U); 0_{l-1} \\
        T_3 &= 0_{l} + \frac{1}{2}(Y+V) \\
        T_4 &= 0_{l} + \frac{1}{2}(Y-V)
        \end{aligned}

    INPUT:

    - ``turyn_sequences`` -- the Turyn sequences that should be used to construct the T-sequences.

    - ``check`` -- boolean, if true (default) checks that the sequences created are T-sequences before returning them.

    EXAMPLES::

        sage: from sage.combinat.t_sequences import turyn_sequences_smallcases, T_sequences_construction_I, is_T_sequences_set
        sage: seqs = turyn_sequences_smallcases(4)
        sage: T_sequences_construction_I(seqs)
        [[1, 1, -1, 0, 0, 0, 0],
        [0, 0, 0, -1, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 1],
        [0, 0, 0, 0, 0, 1, 0]]

    TESTS::

        sage: seqs = turyn_sequences_smallcases(4)
        sage: is_T_sequences_set(T_sequences_construction_I(seqs))
        True
        sage: T_sequences_construction_I([[1, -1], [-1, 1], [1]])
        Traceback (most recent call last):
        ...
        AssertionError

    .. NOTE::

        The construction detailed in [Seb2017]_ contains a typo. The first two sequences are
        defined as `T_1 = \frac{1}{2}(X+Y); 0_{l-1}` and `T_2 = \frac{1}{2}(X-Y); 0_{l-1}`,
        but the correct formulas are `T_1 = \frac{1}{2}(X+U); 0_{l-1}` and
        `T_2 = \frac{1}{2}(X-U); 0_{l-1}`. A very similar typo can be seen in the proof of
        Theorem 9 in [CRSKKY1989]_ as well.
    """

    assert len(turyn_sequences) == 4

    X, U, Y, V = turyn_sequences
    l = len(X)

    assert len(X) == len(U) == len(Y)+1 == len(V)+1

    def seq_sum(seq1, seq2):
        return [(a+b)//2 for (a, b) in zip(seq1, seq2)]
    
    def seq_subtract(seq1, seq2):
        return [(a-b)//2 for (a, b) in zip(seq1, seq2)]

    def zero_seq(n):
        return [0 for _ in range(n)]
    
    X1 = Sequence(seq_sum(X, U) + zero_seq(l-1))
    X2 =  Sequence(seq_subtract(X, U) + zero_seq(l-1))
    X3 =  Sequence(zero_seq(l) + seq_sum(Y, V))
    X4 =  Sequence(zero_seq(l) + seq_subtract(Y, V))

    res = [X1, X2, X3, X4]
    if check:
        assert is_T_sequences_set(res)
    return res

def T_sequences_construction_II(turyn_sequences, check=True):
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

    - ``turyn_sequences`` -- the Turyn sequences that should be used to construct the T-sequences .

    - ``check`` -- boolean, if true (default) checks that the sequences created are T-sequences before returning them.

    EXAMPLES::

        sage: from sage.combinat.t_sequences import turyn_sequences_smallcases, T_sequences_construction_II, is_T_sequences_set
        sage: seqs = turyn_sequences_smallcases(4)
        sage: T_sequences_construction_II(seqs)
        [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, -1, 1, -1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -1, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 1, 0]]

    TESTS::
        sage: seqs = turyn_sequences_smallcases(4)
        sage: is_T_sequences_set(T_sequences_construction_II(seqs))
        True
        sage: T_sequences_construction_II([[1, -1], [-1, 1], [1]])
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
            if i%2 == 0:
                res.append(seq1[i//2])
            else:
                res.append(seq2[i//2])
        return res

    X1 = Sequence([1]+ zero_seq(4*l-2))
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
    :func:`sage.combinat.t_sequences.T_sequences_construction_I` and 
    :func:`sage.combinat.t_sequences.T_sequences_construction_II` 
    together with the Turyn sequences stored in :func:`sage.combinat.t_sequences.turyn_sequences_smallcases`, 
    or some T-sequences taken directly from [CRSKKY1989]_.

    INPUT:

    - ``t`` -- integer, the length of the T-sequences to construct.

    - ``existence`` -- boolean (default false). If true, thsi method only returns whether a T-sequences of 
      the given size can be constructed.

    - ``check`` -- boolean, if true (default) check that the sequences are T-sequences before returning them.

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
            [0]*24+ [1,1,0,-1,0,-1,1,1,-1,0,0,0,0,0,-1,1,-1,-1,0,-1,0,-1,1]
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
        sequences =  list(map(Sequence, db[t]))
        if check:
            assert is_T_sequences_set(sequences)
        return sequences
    if (t+1) %2 == 0 and turyn_sequences_smallcases((t+1)//2, existence=True):
        if existence:
            return True
        turyn_seqs = turyn_sequences_smallcases((t+1)//2)
        return T_sequences_construction_I(turyn_seqs, check=check)
    
    if (t+1)%4 == 0 and turyn_sequences_smallcases((t+1)//4, existence=True):
        if existence:
            return True
        turyn_seqs = turyn_sequences_smallcases((t+1)//4)
        return T_sequences_construction_II(turyn_seqs, check=check)
    
    if existence:
        return False
    raise ValueError(f'T Sequences of length {t} not yet implemented.')
