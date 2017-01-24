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

        sage: from sage.misc.prandom import _pyrand
        sage: _pyrand()
        <random.Random object at 0x...>
        sage: _pyrand().getrandbits(10)
        114L
    """
    return current_randstate().python_random()

def getrandbits(k):
    r"""
    getrandbits(k) -> x.  Generates a long int with k random bits.

    EXAMPLES::

        sage: getrandbits(10)
        114L
        sage: getrandbits(200)
        1251230322675596703523231194384285105081402591058406420468435L
        sage: getrandbits(10)
        533L
    """
    return _pyrand().getrandbits(k)

def randrange(start, stop=None, step=1):
    r"""
    Choose a random item from range(start, stop[, step]).

    This fixes the problem with randint() which includes the
    endpoint; in Python this is usually not what you want.

    EXAMPLES::

        sage: randrange(0, 100, 11)
        11
        sage: randrange(5000, 5100)
        5051
        sage: [randrange(0, 2) for i in range(15)]
        [0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1]
        sage: randrange(0, 1000000, 1000)
        486000
        sage: randrange(-100, 10)
        -56
    """
    return _pyrand().randrange(start, stop, step)

def randint(a, b):
    r"""
    Return random integer in range [a, b], including both end points.

    EXAMPLES::

        sage: [randint(0, 2) for i in range(15)]
        [0, 1, 0, 0, 1, 0, 2, 0, 2, 1, 2, 2, 0, 2, 2]
        sage: randint(-100, 10)
        -46
    """
    return _pyrand().randint(a, b)

def choice(seq):
    r"""
    Choose a random element from a non-empty sequence.

    EXAMPLES::

        sage: [choice(list(primes(10, 100))) for i in range(5)]
        [17, 47, 11, 31, 47]
    """
    return _pyrand().choice(seq)

def shuffle(x, random=None):
    r"""
    x, random=random.random -> shuffle list x in place; return None.

    Optional arg random is a 0-argument function returning a random
    float in [0.0, 1.0); by default, the sage.misc.random.random.

    EXAMPLES::

        sage: shuffle([1 .. 10])
    """
    if random is None:
        random = _pyrand().random
    return _pyrand().shuffle(x, random)

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

        sage: sample(["Here", "I", "come", "to", "save", "the", "day"], 3)
        ['Here', 'to', 'day']
        sage: from six.moves import range
        sage: sample(range(2^30), 7)
        [357009070, 558990255, 196187132, 752551188, 85926697, 954621491, 624802848]
    """
    return _pyrand().sample(population, k)

def random():
    r"""
    Get the next random number in the range [0.0, 1.0).

    EXAMPLES::

        sage: [random() for i in [1 .. 4]]
        [0.111439293741037, 0.5143475134191677, 0.04468968524815642, 0.332490606442413]
    """
    return _pyrand().random()

def uniform(a, b):
    r"""
    Get a random number in the range [a, b).

    Equivalent to \code{a + (b-a) * random()}.

    EXAMPLES::

        sage: uniform(0, 1)
        0.111439293741037
        sage: uniform(e, pi)
        0.5143475134191677*pi + 0.48565248658083227*e
        sage: RR(_)
        2.93601069876846
    """
    return _pyrand().uniform(a, b)

def betavariate(alpha, beta):
    r"""
    Beta distribution.

    Conditions on the parameters are alpha > 0 and beta > 0.
    Returned values range between 0 and 1.

    EXAMPLES::

        sage: betavariate(0.1, 0.9)
        9.75087916621299e-9
        sage: betavariate(0.9, 0.1)
        0.941890400939253
    """
    return _pyrand().betavariate(alpha, beta)

def expovariate(lambd):
    r"""
    Exponential distribution.

    lambd is 1.0 divided by the desired mean.  (The parameter would be
    called "lambda", but that is a reserved word in Python.)  Returned
    values range from 0 to positive infinity.

    EXAMPLES::

        sage: [expovariate(0.001) for i in range(3)]
        [118.152309288166, 722.261959038118, 45.7190543690470]
        sage: [expovariate(1.0) for i in range(3)]
        [0.404201816061304, 0.735220464997051, 0.201765578600627]
        sage: [expovariate(1000) for i in range(3)]
        [0.0012068700332283973, 8.340929747302108e-05, 0.00219877067980605]
    """
    return _pyrand().expovariate(lambd)

def gammavariate(alpha, beta):
    r"""
    Gamma distribution.  Not the gamma function!

    Conditions on the parameters are alpha > 0 and beta > 0.

    EXAMPLES::

        sage: gammavariate(1.0, 3.0)
        6.58282586130638
        sage: gammavariate(3.0, 1.0)
        3.07801512341612
    """
    return _pyrand().gammavariate(alpha, beta)

def gauss(mu, sigma):
    r"""
    Gaussian distribution.

    mu is the mean, and sigma is the standard deviation.  This is
    slightly faster than the normalvariate() function, but is not
    thread-safe.

    EXAMPLES::

       sage: [gauss(0, 1) for i in range(3)]
       [0.9191011757657915, 0.7744526756246484, 0.8638996866800877]
       sage: [gauss(0, 100) for i in range(3)]
       [24.916051749154448, -62.99272061579273, -8.1993122536718...]
       sage: [gauss(1000, 10) for i in range(3)]
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

        sage: [lognormvariate(100, 10) for i in range(3)]
        [2.9410355688290246e+37, 2.2257548162070125e+38, 4.142299451717446e+43]
    """
    return _pyrand().lognormvariate(mu, sigma)

def normalvariate(mu, sigma):
    r"""
    Normal distribution.

    mu is the mean, and sigma is the standard deviation.

    EXAMPLES::

       sage: [normalvariate(0, 1) for i in range(3)]
       [-1.372558980559407, -1.1701670364898928, 0.04324100555110143]
       sage: [normalvariate(0, 100) for i in range(3)]
       [37.45695875041769, 159.6347743233298, 124.1029321124009]
       sage: [normalvariate(1000, 10) for i in range(3)]
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

        sage: [vonmisesvariate(1.0r, 3.0r) for i in range(1, 5)]  # abs tol 1e-12
        [0.898328639355427, 0.6718030007041281, 2.0308777524813393, 1.714325253725145]
    """
    return _pyrand().vonmisesvariate(mu, kappa)

def paretovariate(alpha):
    r"""
    Pareto distribution.  alpha is the shape parameter.

    EXAMPLES::

        sage: [paretovariate(3) for i in range(1, 5)]
        [1.0401699394233033, 1.2722080162636495, 1.0153564009379579, 1.1442323078983077]
    """
    return _pyrand().paretovariate(alpha)

def weibullvariate(alpha, beta):
    r"""
    Weibull distribution.

    alpha is the scale parameter and beta is the shape parameter.

    EXAMPLES::

        sage: [weibullvariate(1, 3) for i in range(1, 5)]
        [0.49069775546342537, 0.8972185564611213, 0.357573846531942, 0.739377255516847]
    """
    return _pyrand().weibullvariate(alpha, beta)
