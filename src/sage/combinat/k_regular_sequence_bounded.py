from __future__ import absolute_import

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
        sage: A = Matrix([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
        sage: red_mult(A, A)
        [2 0 0]
        [0 2 0]
        [0 0 2]
    """
    return (A*B).apply_map(lambda m: min(m, 2))


def construct_phi(matrices):
    r"""
    Return the set `\phi(L)` with `\phi` as defined in [###paper###]
    
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

                                  
def is_bounded(S):
    r"""
    Return ``True`` if S is bounded and ``False`` if S is not bounded

    INPUT:
    
    - ``S`` -- a `k`-regular sequence
    """

    L = S.minimized().mu
    le = len(L)
    for i in srange(1, le):
        M = L[i]
        es = M.eigenvectors_left()
        l = len(es)
        for i in srange(0, l):
            if abs(es[i][0]) > 1 or (abs(es[i][0]) == 1 and len(es[i][1]) < es[i][2]):
                return False

    mandelSimon = True
    posL = []
    l = len(L)
    for i in srange(0, l):
        M = L[i]
        if max(M._list()) <= 0:
            posL.append(Matrix(numpy.absolute(numpy.array(M))))
        elif min(M._list()) >= 0:
            posL.append(M)
        elif max(M._list()) > 0 and min(M._list()) < 0:
            mandelSimon = False
            break

    if mandelSimon == True:
        print 'maldelSimon = True'
        return mandel_simon_algorithm(posL)
    else:
        raise Exception( "Not decidable with this implementation.")
    from sage.arith.srange import srange
    
