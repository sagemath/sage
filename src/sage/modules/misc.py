# sage_setup: distribution = sagemath-modules
"""
Miscellaneous module-related functions

AUTHORS:

- William Stein (2007-11-18)
"""

###############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
###############################################################################

from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ

# Function below could be replicated into
# sage.matrix.matrix_integer_dense.Matrix_integer_dense.is_LLL_reduced
# which is its only current use (2011-02-26).  Then this could
# be deprecated and this file removed.


def gram_schmidt(B):
    r"""
    Return the Gram-Schmidt orthogonalization of the entries in the list
    B of vectors, along with the matrix mu of Gram-Schmidt coefficients.

    Note that the output vectors need not have unit length. We do this
    to avoid having to extract square roots.

    .. NOTE::

        Use of this function is discouraged.  It fails on linearly
        dependent input and its output format is not as natural as it
        could be.  Instead, see :meth:`sage.matrix.matrix2.Matrix2.gram_schmidt`
        which is safer and more general-purpose.

    EXAMPLES::

        sage: B = [vector([1,2,1/5]), vector([1,2,3]), vector([-1,0,0])]
        sage: from sage.modules.misc import gram_schmidt
        sage: G, mu = gram_schmidt(B)
        sage: G
        [(1, 2, 1/5), (-1/9, -2/9, 25/9), (-4/5, 2/5, 0)]
        sage: G[0] * G[1]
        0
        sage: G[0] * G[2]
        0
        sage: G[1] * G[2]
        0
        sage: mu
        [      0       0       0]
        [   10/9       0       0]
        [-25/126    1/70       0]
        sage: a = matrix([])
        sage: a.gram_schmidt()
        ([], [])
        sage: a = matrix([[],[],[],[]])
        sage: a.gram_schmidt()
         ([], [])

    Linearly dependent input leads to a zero dot product in a denominator.
    This shows that :issue:`10791` is fixed. ::

        sage: from sage.modules.misc import gram_schmidt
        sage: V = [vector(ZZ,[1,1]), vector(ZZ,[2,2]), vector(ZZ,[1,2])]
        sage: gram_schmidt(V)
        Traceback (most recent call last):
        ...
        ValueError: linearly dependent input for module version of Gram-Schmidt

    TESTS::

        sage: from sage.modules.misc import gram_schmidt
        sage: V = []
        sage: gram_schmidt(V)
        ([], [])
        sage: V = [vector(ZZ,[0])]
        sage: gram_schmidt(V)
        Traceback (most recent call last):
        ...
        ValueError: linearly dependent input for module version of Gram-Schmidt
    """
    from sage.modules.free_module_element import vector
    if len(B) == 0 or len(B[0]) == 0:
        return B, matrix(ZZ, 0, 0, [])
    n = len(B)
    Bstar = [B[0]]
    K = B[0].base_ring().fraction_field()
    zero = vector(K, len(B[0]))
    if Bstar[0] == zero:
        raise ValueError("linearly dependent input for module version of Gram-Schmidt")
    mu = matrix(K, n, n)
    for i in range(1, n):
        for j in range(i):
            mu[i, j] = B[i].dot_product(Bstar[j]) / (Bstar[j].dot_product(Bstar[j]))
        Bstar.append(B[i] - sum(mu[i, j] * Bstar[j] for j in range(i)))
        if Bstar[i] == zero:
            raise ValueError("linearly dependent input for module version of Gram-Schmidt")
    return Bstar, mu
