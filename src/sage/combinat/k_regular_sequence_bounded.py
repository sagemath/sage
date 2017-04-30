
from sage.all_cmdline import *   # import sage library
import numpy # needed!

def reduce_matrix_max_2(A):
    r"""
    Return the matrix `B` with entries `B_{ij} = \max{A_{ij},2}`.

    INPUT: 

    - ``A`` -- a matrix 

    OUTPUT: 

    A matrix with entries `\max{A_{ij},2}`

    EXAMPLES::
        sage: from sage.combinat.k_automatic_sequence import reduce_matrix_max_2
        sage: A = Matrix([[0, 1, 2], [3, 4, 5], [-6, -7, -8]])
        sage: reduce_matrix_max_2(A)
        [ 0  1  2]
        [ 2  2  2]
        [-6 -7 -8]
    """

    n = A.nrows()
    m = A.ncols()
    B = zero_matrix(n, m)

    for i in srange(n):
        for j in srange(m):
            # TO THINK: negative values?
            if A[i, j] > 2:
                B[i, j] = 2 
            else:
                B[i, j] = A[i, j]
    return B

def red_mult(A, B):
    return reduce_matrix_max_2(A*B)

def construct_phi(L):
    r"""
    Return the set `\phi(L)` with `\phi` as defined in []
    
    INPUT:
    
    - ``L`` -- a list of square matrices with the same dimension

    OUTPUT:
    
    The set `\phi(L)`

    EXAMPLES::
        sage: from sage.combinat.k_automatic_sequence import construct_phi
        sage: L = [Matrix([[2, 2], [1, 3]]), Matrix([[0, 2], [1, 1]])]
        sage: construct_phi(L)
        [
        [0 2]  [2 2]  [2 2]
        [1 1], [2 2], [1 2]
        ]
    """
    
    l = len(L)
    redL = []
    for i in srange(1, l):
        redL.append(reduce_matrix_max_2(L[i]))
    n = len(redL)
    while(True): 
        for A in redL:
            for B in L:
                C = red_mult(A, B)
                if C not in redL:
                    redL.append(C)
        if len(redL) == n:
            return redL
        else:
            n = len(redL)

    
def mandel_simon_algorithm(L):
    A = construct_phi(L)
    for M in A:
        if red_mult(M, M) == red_mult(red_mult(M, M), M) and not M**2 == M**3:
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
    
