"""
Tests for chain complexes, simplicial complexes, etc.

These test whether CHomP gives the same answers as Sage's built-in
homology calculator.

TESTS::

    sage: from sage.homology.tests import test_random_chain_complex
    sage: test_random_chain_complex(trials=20)  # optional - CHomP
    sage: test_random_chain_complex(level=2, trials=20)  # optional - CHomP
    sage: test_random_chain_complex(level=3, trials=20)  # long time # optional - CHomP

    sage: from sage.homology.tests import test_random_simplicial_complex
    sage: test_random_simplicial_complex(level=1, trials=20)  # optional - CHomP
    sage: test_random_simplicial_complex(level=2, trials=20)  # optional - CHomP
    sage: test_random_simplicial_complex(level=5/2, trials=10)  # long time # optional - CHomP
"""
from sage.misc.random_testing import random_testing

def random_chain_complex(level=1):
    """
    Return a random chain complex, defined by specifying a single
    random matrix in a random degree, with differential of degree
    either 1 or -1.  The matrix is randomly sparse or dense.

    :param level: measure of complexity: the larger this is, the
      larger the matrix can be, and the larger its degree can be in
      the chain complex.
    :type level: positive integer; optional, default 1

    EXAMPLES::

        sage: from sage.homology.tests import random_chain_complex
        sage: C = random_chain_complex()
        sage: C
        Chain complex with at most 2 nonzero terms over Integer Ring
        sage: C.degree_of_differential() # random: either 1 or -1
        1
    """
    from sage.misc.prandom import randint
    from sage.matrix.constructor import random_matrix
    from sage.homology.chain_complex import ChainComplex
    from sage.rings.integer_ring import ZZ
    bound = 50*level
    nrows = randint(0, bound)
    ncols = randint(0, bound)
    sparseness = bool(randint(0, 1))
    mat = random_matrix(ZZ, nrows, ncols, sparse=sparseness)
    dim = randint(-bound, bound)
    deg = 2 * randint(0, 1) - 1  # -1 or 1
    return ChainComplex({dim: mat}, degree = deg)

@random_testing
def test_random_chain_complex(level=1, trials=1, verbose=False):
    """
    Compute the homology of a random chain complex with and without
    CHomP, and compare the results.  If they are not the same, raise
    an error.

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
    """
    for i in range(trials):
        C = random_chain_complex(level=level)
        for d in C.differential():
            chomp = C.homology(d, verbose=verbose)
            no_chomp = C.homology(d, algorithm='no_chomp', verbose=verbose)
            if chomp != no_chomp:
                print "Homology in dimension %s according to CHomP: %s" % (d, chomp)
                print "Homology in dimension %s according to Sage: %s" % (d, no_chomp)
                print "Chain complex: %s" % C.differential()
                raise ValueError

def random_simplicial_complex(level=1, p=0.5):
    """
    Return a random simplicial complex.

    :param level: measure of complexity: the larger this is, the more
      vertices and therefore the larger the possible dimension of the
      complex.
    :type level: positive integer; optional, default 1
    :param p: probability, passed on to ``simplicial_complexes.RandomComplex``
    :type p: float between 0 and 1; optional; default 0.5

    EXAMPLES::

        sage: from sage.homology.tests import random_simplicial_complex
        sage: X = random_simplicial_complex()
        sage: X # random
        Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6, 7) and 31 facets
        sage: X.dimension() < 11
        True
    """
    from sage.misc.prandom import randint
    from sage.homology.examples import simplicial_complexes
    n = randint(2, 4*level)
    dim = randint(1, n)
    return simplicial_complexes.RandomComplex(n, dim, p)

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
    """
    for i in range(trials):
        X = random_simplicial_complex(level=level)
        chomp = X.homology(verbose=verbose)
        no_chomp = X.homology(algorithm='no_chomp', verbose=verbose)
        if chomp != no_chomp:
            print "Homology according to CHomP: %s" % chomp
            print "Homology according to Sage: %s" % no_chomp
            print "Simplicial complex: %s" % X
            print "Its chain complex: %s" % X.chain_complex()
            raise ValueError
