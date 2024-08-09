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


<<<<<<< HEAD
@random_testing
def test_random_chain_complex(level=1, trials=1, verbose=False):
    """
    Compute the homology of a random chain complex with and without
    CHomP, and compare the results.  If they are not the same, raise
    an error.

    This function is deprecated: see :issue:`33777`.

    :param level: measure of complexity of the chain complex -- see
      :func:`random_chain_complex`
    :type level: positive integer; optional, default 1
    :param trials: number of trials to conduct
    :type trials: positive integer; optional, default 1
    :param verbose: if ``True``, print verbose messages
    :type verbose: boolean; optional, default ``False``

    EXAMPLES::

        sage: from sage.homology.tests import test_random_chain_complex
        sage: test_random_chain_complex(trials=2)  # optional - CHomP
        doctest:...: DeprecationWarning: the CHomP interface is deprecated; hence so is this function
        See https://github.com/sagemath/sage/issues/33777 for details.
        done
    """
    deprecation(33777, 'the CHomP interface is deprecated; hence so is this function')
    for i in range(trials):
        C = random_chain_complex(level=level)
        for d in C.differential():
            chomp = C.homology(d, verbose=verbose)
            no_chomp = C.homology(d, algorithm='no_chomp', verbose=verbose)
            if chomp != no_chomp:
                print("Homology in dimension %s according to CHomP: %s" % (d, chomp))
                print("Homology in dimension %s according to Sage: %s" % (d, no_chomp))
                print("Chain complex: %s" % C.differential())
                raise ValueError
    print("done")


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


@random_testing
def test_random_simplicial_complex(level=1, trials=1, verbose=False):
    """
    Compute the homology of a random simplicial complex with and
    without CHomP, and compare the results.  If they are not the same,
    raise an error.

    :param level: measure of complexity of the simplicial complex --
      see :func:`random_simplicial_complex`
    :type level: positive integer; optional, default 1
    :param trials: number of trials to conduct
    :type trials: positive integer; optional, default 1
    :param verbose: if ``True``, print verbose messages
    :type verbose: boolean; optional, default ``False``

    This gets pretty slow if ``level`` is more than 3.

    EXAMPLES::

        sage: from sage.homology.tests import test_random_simplicial_complex
        sage: test_random_simplicial_complex(trials=2)  # optional - CHomP
        doctest:...: DeprecationWarning: the CHomP interface is deprecated; hence so is this function
        See https://github.com/sagemath/sage/issues/33777 for details.
        done
    """
    deprecation(33777, 'the CHomP interface is deprecated; hence so is this function')
    for i in range(trials):
        X = random_simplicial_complex(level=level)
        chomp = X.homology(verbose=verbose)
        no_chomp = X.homology(algorithm='no_chomp', verbose=verbose)
        if chomp != no_chomp:
            print("Homology according to CHomP: %s" % chomp)
            print("Homology according to Sage: %s" % no_chomp)
            print("Simplicial complex: %s" % X)
            print("Its chain complex: %s" % X.chain_complex())
            raise ValueError
    print("done")
=======
>>>>>>> origin/develop
