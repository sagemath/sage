# sage_setup: distribution = sagemath-objects
r"""
Random Numbers with Python API

AUTHORS:
    -- Carl Witty (2008-03): new file

This module has the same functions as the Python standard module
\module{random}, but uses the current \sage random number state from
\module{sage.misc.randstate} (so that it can be controlled by the same
global random number seeds).

The functions here are less efficient than the functions in \module{random},
because they look up the current random number state on each call.

If you are going to be creating many random numbers in a row, it is
better to use the functions in \module{sage.misc.randstate} directly.

Here is an example:

(The imports on the next two lines are not necessary, since
\function{randrange} and \function{current_randstate} are both available
by default at the \code{sage:} prompt; but you would need them
to run these examples inside a module.) ::

    sage: from sage.misc.prandom import randrange
    sage: from sage.misc.randstate import current_randstate
    sage: def test1():
    ....:    return sum([randrange(100) for i in range(100)])
    sage: def test2():
    ....:    randrange = current_randstate().python_random().randrange
    ....:    return sum([randrange(100) for i in range(100)])

Test2 will be slightly faster than test1, but they give the same answer::

    sage: with seed(0): test1()
    5169
    sage: with seed(0): test2()
    5169
    sage: with seed(1): test1()
    5097
    sage: with seed(1): test2()
    5097
    sage: timeit('test1()') # random
    625 loops, best of 3: 590 us per loop
    sage: timeit('test2()') # random
    625 loops, best of 3: 460 us per loop

The docstrings for the functions in this file are mostly copied from
Python's \file{random.py}, so those docstrings are "Copyright (c)
2001, 2002, 2003, 2004, 2005, 2006, 2007 Python Software Foundation;
All Rights Reserved" and are available under the terms of the
Python Software Foundation License Version 2.
"""

# We deliberately omit "seed" and several other seed-related functions...
# setting seeds should only be done through sage.misc.randstate .

from sage.misc.randstate import current_randstate


def _pyrand():
    r"""
    A tiny private helper function to return an instance of
    random.Random from the current \sage random number state.
    Only for use in prandom.py; other modules should use
    current_randstate().python_random().

    EXAMPLES::

        sage: set_random_seed(0)
        sage: from sage.misc.prandom import _pyrand
        sage: _pyrand()
        <...random.Random object at 0x...>
        sage: _pyrand().getrandbits(10)
        114
    """
    return current_randstate().python_random()


def getrandbits(k):
    r"""
    getrandbits(k) -> x.  Generates a long int with k random bits.

    EXAMPLES::

        sage: getrandbits(10) in range(2^10)
        True
        sage: getrandbits(200) in range(2^200)
        True
        sage: getrandbits(4) in range(2^4)
        True
    """
    return _pyrand().getrandbits(k)


def randrange(start, stop=None, step=1):
    r"""
    Choose a random item from range(start, stop[, step]).

    This fixes the problem with randint() which includes the
    endpoint; in Python this is usually not what you want.

    EXAMPLES::

        sage: s = randrange(0, 100, 11)
        sage: 0 <= s < 100
        True
        sage: s % 11
        0

        sage: 5000 <= randrange(5000, 5100) < 5100
        True
        sage: s = [randrange(0, 2) for i in range(15)]
        sage: all(t in [0, 1] for t in s)
        True

        sage: s = randrange(0, 1000000, 1000)
        sage: 0 <= s < 1000000
        True
        sage: s % 1000
        0
        sage: -100 <= randrange(-100, 10) < 10
        True
    """
    return _pyrand().randrange(start, stop, step)


def randint(a, b):
    r"""
    Return random integer in range [a, b], including both end points.

    EXAMPLES::

        sage: s = [randint(0, 2) for i in range(15)]; s  # random
        [0, 1, 0, 0, 1, 0, 2, 0, 2, 1, 2, 2, 0, 2, 2]
        sage: all(t in [0, 1, 2] for t in s)
        True
        sage: -100 <= randint(-100, 10) <= 10
        True
    """
    return _pyrand().randint(a, b)


def choice(seq):
    r"""
    Choose a random element from a non-empty sequence.

    EXAMPLES::

        sage: s = [choice(list(primes(10, 100))) for i in range(5)]; s  # random        # needs sage.libs.pari
        [17, 47, 11, 31, 47]
        sage: all(t in primes(10, 100) for t in s)                                      # needs sage.libs.pari
        True
    """
    return _pyrand().choice(seq)


def shuffle(x):
    r"""
    x, random=random.random -> shuffle list x in place; return None.

    Optional arg random is a 0-argument function returning a random
    float in [0.0, 1.0); by default, the sage.misc.random.random.

    EXAMPLES::

        sage: shuffle([1 .. 10])
    """
    return _pyrand().shuffle(x)


def sample(population, k):
    r"""
    Choose k unique random elements from a population sequence.

    Return a new list containing elements from the population while
    leaving the original population unchanged.  The resulting list is
    in selection order so that all sub-slices will also be valid random
    samples.  This allows raffle winners (the sample) to be partitioned
    into grand prize and second place winners (the subslices).

    Members of the population need not be hashable or unique.  If the
    population contains repeats, then each occurrence is a possible
    selection in the sample.

    To choose a sample in a range of integers, use xrange as an
    argument (in Python 2) or range (in Python 3).  This is especially
    fast and space efficient for sampling from a large population:
    sample(range(10000000), 60)

    EXAMPLES::

        sage: from sage.misc.misc import is_sublist
        sage: l = ["Here", "I", "come", "to", "save", "the", "day"]
        sage: s = sample(l, 3); s  # random
        ['Here', 'to', 'day']
        sage: is_sublist(sorted(s), sorted(l))
        True
        sage: len(s)
        3

        sage: s = sample(range(2^30), 7); s  # random
        [357009070, 558990255, 196187132, 752551188, 85926697, 954621491, 624802848]
        sage: len(s)
        7
        sage: all(t in range(2^30) for t in s)
        True
    """
    return _pyrand().sample(population, k)


def random():
    r"""
    Get the next random number in the range [0.0, 1.0).

    EXAMPLES::

        sage: sample = [random() for i in [1 .. 4]]; sample  # random
        [0.111439293741037, 0.5143475134191677, 0.04468968524815642, 0.332490606442413]
        sage: all(0.0 <= s <= 1.0 for s in sample)
        True
    """
    return _pyrand().random()


def uniform(a, b):
    r"""
    Get a random number in the range [a, b).

    Equivalent to \code{a + (b-a) * random()}.

    EXAMPLES::

        sage: s = uniform(0, 1); s  # random
        0.111439293741037
        sage: 0.0 <= s <= 1.0
        True

        sage: s = uniform(e, pi); s  # random                                           # needs sage.symbolic
        0.5143475134191677*pi + 0.48565248658083227*e
        sage: bool(e <= s <= pi)                                                        # needs sage.symbolic
        True
    """
    return _pyrand().uniform(a, b)


def betavariate(alpha, beta):
    r"""
    Beta distribution.

    Conditions on the parameters are alpha > 0 and beta > 0.
    Returned values range between 0 and 1.

    EXAMPLES::

        sage: s = betavariate(0.1, 0.9); s  # random
        9.75087916621299e-9
        sage: 0.0 <= s <= 1.0
        True

        sage: s = betavariate(0.9, 0.1); s  # random
        0.941890400939253
        sage: 0.0 <= s <= 1.0
        True
    """
    return _pyrand().betavariate(alpha, beta)


def expovariate(lambd):
    r"""
    Exponential distribution.

    lambd is 1.0 divided by the desired mean.  (The parameter would be
    called "lambda", but that is a reserved word in Python.)  Returned
    values range from 0 to positive infinity.

    EXAMPLES::

        sage: sample = [expovariate(0.001) for i in range(3)]; sample  # random
        [118.152309288166, 722.261959038118, 45.7190543690470]
        sage: all(s >= 0.0 for s in sample)
        True

        sage: sample = [expovariate(1.0) for i in range(3)]; sample  # random
        [0.404201816061304, 0.735220464997051, 0.201765578600627]
        sage: all(s >= 0.0 for s in sample)
        True

        sage: sample = [expovariate(1000) for i in range(3)]; sample  # random
        [0.0012068700332283973, 8.340929747302108e-05, 0.00219877067980605]
        sage: all(s >= 0.0 for s in sample)
        True
    """
    return _pyrand().expovariate(lambd)


def gammavariate(alpha, beta):
    r"""
    Gamma distribution.  (Not the gamma function.)

    Conditions on the parameters are alpha > 0 and beta > 0.

    EXAMPLES::

        sage: sample = gammavariate(1.0, 3.0); sample  # random
        6.58282586130638
        sage: sample > 0
        True
        sage: sample = gammavariate(3.0, 1.0); sample  # random
        3.07801512341612
        sage: sample > 0
        True
    """
    return _pyrand().gammavariate(alpha, beta)


def gauss(mu, sigma):
    r"""
    Gaussian distribution.

    mu is the mean, and sigma is the standard deviation.  This is
    slightly faster than the normalvariate() function, but is not
    thread-safe.

    EXAMPLES::

       sage: [gauss(0, 1) for i in range(3)]  # random
       [0.9191011757657915, 0.7744526756246484, 0.8638996866800877]
       sage: [gauss(0, 100) for i in range(3)]  # random
       [24.916051749154448, -62.99272061579273, -8.1993122536718...]
       sage: [gauss(1000, 10) for i in range(3)]  # random
       [998.7590700045661, 996.1087338511692, 1010.1256817458031]
    """
    return _pyrand().gauss(mu, sigma)


def lognormvariate(mu, sigma):
    r"""
    Log normal distribution.

    If you take the natural logarithm of this distribution, you'll get a
    normal distribution with mean mu and standard deviation sigma.
    mu can have any value, and sigma must be greater than zero.

    EXAMPLES::

        sage: [lognormvariate(100, 10) for i in range(3)]  # random
        [2.9410355688290246e+37, 2.2257548162070125e+38, 4.142299451717446e+43]
    """
    return _pyrand().lognormvariate(mu, sigma)


def normalvariate(mu, sigma):
    r"""
    Normal distribution.

    mu is the mean, and sigma is the standard deviation.

    EXAMPLES::

       sage: [normalvariate(0, 1) for i in range(3)]  # random
       [-1.372558980559407, -1.1701670364898928, 0.04324100555110143]
       sage: [normalvariate(0, 100) for i in range(3)]  # random
       [37.45695875041769, 159.6347743233298, 124.1029321124009]
       sage: [normalvariate(1000, 10) for i in range(3)]  # random
       [1008.5303090383741, 989.8624892644895, 985.7728921150242]
    """
    return _pyrand().normalvariate(mu, sigma)


def vonmisesvariate(mu, kappa):
    r"""
    Circular data distribution.

    mu is the mean angle, expressed in radians between 0 and 2*pi, and
    kappa is the concentration parameter, which must be greater than or
    equal to zero.  If kappa is equal to zero, this distribution reduces
    to a uniform random angle over the range 0 to 2*pi.

    EXAMPLES::

        sage: sample = [vonmisesvariate(1.0r, 3.0r) for i in range(1, 5)]; sample  # random
        [0.898328639355427, 0.6718030007041281, 2.0308777524813393, 1.714325253725145]
        sage: all(s >= 0.0 for s in sample)
        True
    """
    return _pyrand().vonmisesvariate(mu, kappa)


def paretovariate(alpha):
    r"""
    Pareto distribution.  alpha is the shape parameter.

    EXAMPLES::

        sage: sample = [paretovariate(3) for i in range(1, 5)]; sample  # random
        [1.0401699394233033, 1.2722080162636495, 1.0153564009379579, 1.1442323078983077]
        sage: all(s >= 1.0 for s in sample)
        True
    """
    return _pyrand().paretovariate(alpha)


def weibullvariate(alpha, beta):
    r"""
    Weibull distribution.

    alpha is the scale parameter and beta is the shape parameter.

    EXAMPLES::

        sage: sample = [weibullvariate(1, 3) for i in range(1, 5)]; sample  # random
        [0.49069775546342537, 0.8972185564611213, 0.357573846531942, 0.739377255516847]
        sage: all(s >= 0.0 for s in sample)
        True
    """
    return _pyrand().weibullvariate(alpha, beta)
