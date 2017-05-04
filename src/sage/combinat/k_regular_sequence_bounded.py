r"""
Boundedness of `k`-regular sequences

AUTHORS:

- Gabriel Lipnik (2017)
"""

#*****************************************************************************
#       Copyright (C) 2017 Gabriel Lipnik <galipnik@edu.aau.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_function

def red_mult(A, B):
    r"""
    Return the matrix `A\cdot B` with entries `\min{(A\cdot B)_{ij},2}`.

    INPUT: 

    - ``A`` -- an `m \times n` matrix
    - ``B`` -- an `n \times p` matrix

    OUTPUT: 

    An `m \times p` matrix with entries `\min{(A\cdot B)_{ij},2}`

    EXAMPLES::

        sage: from sage.combinat.k_regular_sequence_bounded import red_mult
        sage: A = Matrix([[2, 0], [0, 2]])
        sage: B = Matrix([[-2, 0], [0, 2]])
        sage: A*B
        [-4  0]
        [ 0  4]        
        sage: red_mult(A, B)
        [-4  0]
        [ 0  2]

    ::

        sage: A = Matrix([[1, 2, 3], [-1, -2, -3], [1, 2, 3]])
        sage: B = Matrix([[1, 2, 3], [2, 3, 4], [1, 2, 3]])
        sage: A*B
        [  8  14  20]
        [ -8 -14 -20]
        [  8  14  20]
        sage: red_mult(A, B)
        [  2   2   2]
        [ -8 -14 -20]
        [  2   2   2]
    """
    return (A*B).apply_map(lambda m: min(m, 2))


def construct_phi(matrices):
    r"""
    Return the set `\phi(\text{span}\ L)` with `\phi` as defined in [###paper###]
    
    INPUT:
    
    - ``L`` -- a list of non-negative square matrices with the same dimension

    OUTPUT:
    
    The set `\phi(L)`

    EXAMPLES::
        sage: from sage.combinat.k_regular_sequence_bounded import construct_phi
        sage: L = [Matrix([[2, 2], [1, 3]]), Matrix([[0, 2], [1, 1]])]
        sage: construct_phi(L)
        [
        [0 2]  [2 2]  [2 2]
        [1 1], [2 2], [1 2]
        ]
    """
    redMatrices = list(ell.apply_map(lambda m: min(m,2)) for ell in matrices[1:])
    redLength = len(redMatrices)
    counter = 1
    while(True):
        for A in redMatrices:
            for B in matrices:
                prod = red_mult(B, A)
                if prod not in redMatrices:
                    redMatrices.append(prod)
        if len(redMatrices) == redLength:
            return redMatrices
        else:
            redLength = len(redMatrices)
        counter = counter + 1
        if counter > 1000000:
            raise RuntimeError('while loop too long...')


def mandel_simon_algorithm(matrices):
    from sage.arith.srange import srange
    from sage.matrix.constructor import Matrix
    from sage.matrix.matrix_space import MatrixSpace
    from sage.rings.integer_ring import ZZ

    length = len(matrices)
    if not all(matrices[i] in MatrixSpace(ZZ, matrices[i].nrows(), matrices[i].ncols()) for i in srange(1, length)):
        raise ValueError('Matrix is not integer-valued.')

    def is_non_negative(M):
        return min(M.list()) >= 0

    posMatrices = list()
    for i in srange(length):
        M = matrices[i]
        if is_non_negative(M):
            posMatrices.append(M)
        elif is_non_negative(-M):
            posMatrices.append(-M)
        else:
            raise ValueError('M[{}] is neither non-negative nor non-positive.'.format(i))

    constrMatrices = construct_phi(posMatrices)
    return not any(red_mult(M, M) == red_mult(red_mult(M, M), M) and not M**2 == M**3
                   for M in constrMatrices)


def check_eigenvalues(matrices):
    from sage.matrix.constructor import Matrix
    from sage.arith.srange import srange

    for mat in matrices:
        eigen = mat.eigenvectors_right()
        l = len(eigen)
        isBounded = all(abs(eigen[i][0]) < 1 or (abs(eigen[i][0]) == 1 and
                        len(eigen[i][1]) == eigen[i][2]) for i in srange(0, l))
        if isBounded == False:
            return False
    return True

@cached_function
def k_regular_sequence_is_bounded(seq):
    r"""
    Return whether this `k`-regular sequence is bounded or not (if decidable with this implementation)

    INPUT:
    
    - ``S`` -- a `k`-regular sequence

    EXAMPLES::
        sage: from sage.combinat.k_regular_sequence_bounded import k_regular_sequence_is_bounded
        sage: Seq2 = kRegularSequenceSpace(2, ZZ)
        doctest:...: FutureWarning: This class/method/function is
        marked as experimental. It, its functionality or its interface
        might change without a formal deprecation.
        See http://trac.sagemath.org/21202 for details.
        sage: S = Seq2((Matrix([[0, 1, 0], [0, 0, 1], [-1, 2, 0]]), Matrix([[-1, 0, 0], 
        ....: [-3/4, -1/4, 3/4], [-1/4, 1/4, -3/4]])), left=vector([1, 0, 0]), right=vector([-4, -4, -4]))
        sage: k_regular_sequence_is_bounded(S)
        ev2
        False
    
    """
    from sage.arith.srange import srange
    matrices = seq.mu
    length = len(matrices)
    matricesWithout = list(matrices[i] for i in srange(1, length))
    try:
        if mandel_simon_algorithm(matricesWithout):
            print 'ms1'
            return True
    except ValueError:
        print 'ms1err'
        pass
    
    matrices = seq.minimized().mu
    matricesWithout = list(matrices[i] for i in srange(1, length))
    if not check_eigenvalues(matricesWithout):
        print 'ev1'
        return False

    try:
        print 'ms2'
        return mandel_simon_algorithm(matricesWithout)
    except ValueError:
        print 'm2err'
        pass

    matricesProd2 = list(matrices[i]*matrices[j] for i in srange(length) for j in srange(1, length) if i != j)
    if not check_eigenvalues(matricesProd):
        print 'ev2'
        return False

    raise RuntimeError('It is not decidable with this implementation whether the sequence is bounded or not.')
