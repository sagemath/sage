# sage.doctest: needs sage.modules
"""
Tests for chain complexes, simplicial complexes, etc.
"""

from sage.misc.random_testing import random_testing
from sage.misc.prandom import randint
from sage.matrix.constructor import random_matrix
from sage.homology.chain_complex import ChainComplex
from sage.rings.integer_ring import ZZ
from sage.topology.simplicial_complex_examples import RandomComplex


def random_chain_complex(level=1):
    """
    Return a random chain complex, defined by specifying a single
    random matrix in a random degree, with differential of degree
    either 1 or -1.  The matrix is randomly sparse or dense.

    - ``level`` -- positive integer (default: 1); measure of complexity. The larger
      this is, the larger the matrix can be, and the larger its degree can be
      in the chain complex.

    EXAMPLES::

        sage: from sage.homology.tests import random_chain_complex
        sage: C = random_chain_complex()
        sage: C  # random
        Chain complex with at most ... nonzero terms over Integer Ring
        sage: len(C.nonzero_degrees()) in [0, 1, 2]
        True
        sage: C.degree_of_differential() in [-1, 1]
        True
    """
    bound = 50 * level
    nrows = randint(0, bound)
    ncols = randint(0, bound)
    sparseness = bool(randint(0, 1))
    mat = random_matrix(ZZ, nrows, ncols, sparse=sparseness)
    dim = randint(-bound, bound)
    deg = 2 * randint(0, 1) - 1  # -1 or 1
    return ChainComplex({dim: mat}, degree=deg)


def random_simplicial_complex(level=1, p=0.5):
    """
    Return a random simplicial complex.

    - ``level`` -- positive integer (default: 1); measure of complexity. The
      larger this is, the more vertices and therefore the larger the possible
      dimension of the complex.
    - ``p`` -- float between 0 and 1 (default: 0.5); probability, passed on to
      ``simplicial_complexes.RandomComplex``

    EXAMPLES::

        sage: from sage.homology.tests import random_simplicial_complex
        sage: X = random_simplicial_complex()
        sage: X  # random
        Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6, 7) and 31 facets
        sage: X.dimension() < 11
        True
    """
    n = randint(2, 4 * level)
    dim = randint(1, n)
    return RandomComplex(n, dim, p)
