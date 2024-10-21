# sage_setup: distribution = sagemath-modules
r"""
Probability Distributions

This module provides three types of probability distributions:

- :class:`RealDistribution`: various real-valued probability distributions.

- :class:`SphericalDistribution`: uniformly distributed points on the
  surface of an `n-1` sphere in `n` dimensional euclidean space.

- :class:`GeneralDiscreteDistribution`: user-defined discrete distributions.

AUTHORS:

- Josh Kantor (2007-02): first version

- William Stein (2007-02): rewrite of docs, conventions, etc.

- Carlo Hamalainen (2008-08): full doctest coverage, more documentation,
  GeneralDiscreteDistribution, misc fixes.

- Kwankyu Lee (2010-05-29): F-distribution support.

REFERENCES:

    GNU gsl library, General discrete distributions
    http://www.gnu.org/software/gsl/manual/html_node/General-Discrete-Distributions.html

    GNU gsl library, Random number distributions
    http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Distributions.html
"""
# ****************************************************************************
#       Copyright (C) 2004, 2005, 2006 Joshua Kantor <kantor.jm@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from cysignals.memory cimport sig_malloc, sig_free

from sage.libs.gsl.all cimport *
import sage.misc.prandom as random
import sage.rings.real_double
from sage.modules.free_module_element import vector

# TODO: Add more distributions available in gsl
# available but not currently wrapped are laplace, cauchy, landau, logistic.

cdef enum:
    uniform
    gaussian
    rayleigh
    lognormal
    pareto
    t
    F
    chisquared
    exppow
    weibull
    beta
    exponential
    gamma

cdef class ProbabilityDistribution:
    r"""
    Concrete probability distributions should be derived from this
    abstract class.
    """

    def __init__(self):
        r"""
        To be implemented by a derived class::

            sage: P = sage.probability.probability_distribution.ProbabilityDistribution()
        """

        pass

    def get_random_element(self):
        r"""
        To be implemented by a derived class::

            sage: P = sage.probability.probability_distribution.ProbabilityDistribution()
            sage: P.get_random_element()
            Traceback (most recent call last):
            ...
            NotImplementedError: implement in derived class
        """

        raise NotImplementedError("implement in derived class")

    def generate_histogram_data(self, num_samples=1000, bins=50):
        r"""
        Compute a histogram of the probability distribution.

        INPUT:

        - ``num_samples`` -- (optional) number of times to sample from
          the probability distribution

        - ``bins`` -- (optional) number of bins to divide the samples
          into

        OUTPUT:

        - a tuple. The first element of the tuple is a list of length
          ``bins``, consisting of the normalised histogram of the random
          samples. The second list is the bins.

        EXAMPLES::

            sage: set_random_seed(0)
            sage: from sage.probability.probability_distribution import GeneralDiscreteDistribution
            sage: P = [0.3, 0.4, 0.3]
            sage: X = GeneralDiscreteDistribution(P)
            sage: h, b = X.generate_histogram_data(bins=10)                             # needs sage.plot
            sage: h  # rel tol 1e-08                                                    # needs sage.plot
            [1.6299999999999999,
             0.0,
             0.0,
             0.0,
             0.0,
             1.9049999999999985,
             0.0,
             0.0,
             0.0,
             1.4650000000000003]
            sage: b                                                                     # needs sage.plot
            [0.0,
             0.2,
             0.4,
             0.6000000000000001,
             0.8,
             1.0,
             1.2000000000000002,
             1.4000000000000001,
             1.6,
             1.8,
             2.0]
        """
        import numpy as np
        if int(np.version.short_version[0]) > 1:
            np.set_printoptions(legacy="1.25")
        import pylab
        ell = [float(self.get_random_element()) for _ in range(num_samples)]
        S = pylab.hist(ell, bins, density=True)
        return [list(S[0]), list(S[1])]

    def generate_histogram_plot(self, name, num_samples=1000, bins=50):
        r"""
        Save the histogram from :func:`generate_histogram_data() <sage.libs.gsl.ProbabilityDistribution.generate_histogram_data>`
        to a file.

        INPUT:

        - ``name`` -- file to save the histogram plot (as a PNG)

        - ``num_samples`` -- (optional) number of times to sample from
          the probability distribution

        - ``bins`` -- (optional) number of bins to divide the samples
          into

        EXAMPLES:

        This saves the histogram plot to a temporary file::

            sage: from sage.probability.probability_distribution import GeneralDiscreteDistribution
            sage: import tempfile
            sage: P = [0.3, 0.4, 0.3]
            sage: X = GeneralDiscreteDistribution(P)
            sage: with tempfile.NamedTemporaryFile() as f:                              # needs sage.plot
            ....:     X.generate_histogram_plot(f.name)
        """
        import pylab
        ell = [float(self.get_random_element()) for _ in range(num_samples)]
        pylab.hist(ell, bins, density=True)
        pylab.savefig(name)


cdef class SphericalDistribution(ProbabilityDistribution):
    r"""
    This class is capable of producing random points uniformly distributed
    on the surface of an `(n-1)`-sphere in `n`-dimensional euclidean space. The
    dimension `n` is selected via the keyword ``dimension``. The random
    number generator which drives it can be selected using the keyword
    ``rng``. Valid choices are ``'default'`` which uses the Mersenne-Twister,
    ``'luxury'`` which uses RANDLXS, and ``'taus'`` which uses the tausworth
    generator. The default dimension is ``3``.

    EXAMPLES::

        sage: T = SphericalDistribution()
        sage: s = T.get_random_element()
        sage: s.norm()  # rel tol 1e-14
        1.0
        sage: len(s)
        3
        sage: T = SphericalDistribution(dimension=4, rng='luxury')
        sage: s = T.get_random_element()
        sage: s.norm()  # rel tol 1e-14
        1.0
        sage: len(s)
        4

    TESTS:

    Make sure that repeated initializations are randomly seeded
    (:issue:`9770`)::

        sage: Xs = [tuple(SphericalDistribution(2).get_random_element()) for _ in range(1000)]
        sage: len(set(Xs)) > 2^^32
        True
    """

    cdef gsl_rng *r
    cdef gsl_rng_type *T
    cdef long int seed
    cdef Py_ssize_t dimension
    cdef double* vec

    def __init__(self, dimension=3, rng='default', seed=None):
        r"""
        EXAMPLES::

            sage: T = SphericalDistribution()
            sage: T.get_random_element().norm()  # rel tol 1e-14
            1.0

        TESTS:

        Until :issue:`15089` a value of the ``seed`` keyword
        besides ``None`` was ignored. We check here that setting
        a seed is effective. ::

            sage: T = SphericalDistribution(seed=876)
            sage: one = [T.get_random_element() for _ in range(10)]
            sage: T = SphericalDistribution(seed=876)
            sage: two = [T.get_random_element() for _ in range(10)]
            sage: T = SphericalDistribution(seed=123)
            sage: three = [T.get_random_element() for _ in range(10)]
            sage: one == two
            True
            sage: one == three
            False
        """
        gsl_rng_env_setup()
        self.set_random_number_generator(rng)
        self.r = gsl_rng_alloc(self.T)
        if seed is None:
            seed = random.randint(1, 2**31)
        self.set_seed(seed)
        self.dimension = dimension
        self.vec = <double *>sig_malloc(self.dimension*(sizeof(double)))

    def set_seed(self, seed):
        r"""
        Set the seed for the underlying random number generator.

        EXAMPLES::

            sage: T = SphericalDistribution(seed=0)
            sage: T.set_seed(100)
        """
        gsl_rng_set(self.r, seed)
        self.seed = seed

    def set_random_number_generator(self, rng='default'):
        r"""
        Set the gsl random number generator to be one of ``default``,
        ``luxury``, or ``taus``.

        EXAMPLES::

            sage: T = SphericalDistribution()
            sage: T.set_random_number_generator('default')
            sage: T.set_seed(0)
            sage: T.get_random_element()  # rel tol 4e-16
            (0.07961564104639995, -0.05237671627581255, 0.9954486572862178)
            sage: T.set_random_number_generator('luxury')
            sage: T.set_seed(0)
            sage: T.get_random_element()  # rel tol 4e-16
            (0.07961564104639995, -0.05237671627581255, 0.9954486572862178)
        """
        if rng == 'default':
            self.T = gsl_rng_default
        elif rng == 'luxury':
            self.T = gsl_rng_ranlxd2
        elif rng == 'taus':
            self.T = gsl_rng_taus2
        else:
            raise TypeError("Not a valid random number generator")

    def __dealloc__(self):
        if self.r != NULL:
            gsl_rng_free(self.r)
        sig_free(self.vec)

    def get_random_element(self):
        r"""
        Get a random sample from the probability distribution.

        EXAMPLES::

            sage: T = SphericalDistribution(seed=0)
            sage: T.get_random_element()  # rel tol 4e-16
            (0.07961564104639995, -0.05237671627581255, 0.9954486572862178)
        """
        cdef int i
        v = [0]*self.dimension
        gsl_ran_dir_nd(self.r, self.dimension, self.vec)
        for i in range(self.dimension):
            v[i] = self.vec[i]
        return vector(sage.rings.real_double.RDF, v)  # This could be made more efficient by directly constructing the vector, TODO.

    def reset_distribution(self):
        r"""
        This method resets the distribution.

        EXAMPLES::

            sage: T = SphericalDistribution(seed=0)
            sage: [T.get_random_element() for _ in range(4)]  # rel tol 4e-16
            [(0.07961564104639995, -0.05237671627581255, 0.9954486572862178),
             (0.4123599490593727, 0.5606817859360097, -0.7180495855658982),
             (-0.9619860891623148, -0.2726473494040498, -0.015690351211529927),
             (0.5674297579435619, -0.011206783800420301, -0.8233455397322326)]
            sage: T.reset_distribution()
            sage: [T.get_random_element() for _ in range(4)]  # rel tol 4e-16
            [(0.07961564104639995, -0.05237671627581255, 0.9954486572862178),
             (0.4123599490593727, 0.5606817859360097, -0.7180495855658982),
             (-0.9619860891623148, -0.2726473494040498, -0.015690351211529927),
             (0.5674297579435619, -0.011206783800420301, -0.8233455397322326)]
        """
        if self.r != NULL:
            gsl_rng_free(self.r)
        self.r = gsl_rng_alloc(self.T)
        self.set_seed(self.seed)
#        gsl_rng_env_setup()

cdef class RealDistribution(ProbabilityDistribution):
    r"""
    The :class:`RealDistribution` class provides a number of routines for sampling
    from and analyzing and visualizing probability distributions.
    For precise definitions of the distributions and their parameters
    see the gsl reference manuals chapter on random number generators
    and probability distributions.

    EXAMPLES:

    Uniform distribution on the interval ``[a, b]``::

        sage: a = 0
        sage: b = 2
        sage: T = RealDistribution('uniform', [a, b])
        sage: a <= T.get_random_element() <= b
        True
        sage: T.distribution_function(0)
        0.5
        sage: T.cum_distribution_function(1)
        0.5
        sage: T.cum_distribution_function_inv(.5)
        1.0

    The gaussian distribution takes 1 parameter ``sigma``. The standard
    gaussian distribution has ``sigma = 1``::

        sage: sigma = 1
        sage: T = RealDistribution('gaussian', sigma)
        sage: s = T.get_random_element()
        sage: s.parent()
        Real Double Field
        sage: T.distribution_function(0)
        0.3989422804014327
        sage: T.cum_distribution_function(1)
        0.8413447460685429
        sage: T.cum_distribution_function_inv(.5)
        0.0

    The rayleigh distribution has 1 parameter ``sigma``::

        sage: sigma = 3
        sage: T = RealDistribution('rayleigh', sigma)
        sage: s = T.get_random_element()
        sage: s >= 0
        True
        sage: s.parent()
        Real Double Field
        sage: T.distribution_function(0)
        0.0
        sage: T.cum_distribution_function(1)
        0.054040531093234534
        sage: T.cum_distribution_function_inv(.5)
        3.532230067546424...

    The lognormal distribution has two parameters ``sigma``
    and ``zeta``::

        sage: zeta = 0
        sage: sigma = 1
        sage: T = RealDistribution('lognormal', [zeta, sigma])
        sage: s = T.get_random_element()
        sage: s >= 0
        True
        sage: s.parent()
        Real Double Field
        sage: T.distribution_function(0)
        0.0
        sage: T.cum_distribution_function(1)
        0.5
        sage: T.cum_distribution_function_inv(.5)
        1.0

    The pareto distribution has two parameters ``a``, and ``b``::

        sage: a = 1
        sage: b = 1
        sage: T = RealDistribution('pareto', [a, b])
        sage: s = T.get_random_element()
        sage: s >= b
        True
        sage: s.parent()
        Real Double Field
        sage: T.distribution_function(0)
        0.0
        sage: T.cum_distribution_function(1)
        0.0
        sage: T.cum_distribution_function_inv(.5)
        2.0

    The t-distribution has one parameter ``nu``::

        sage: nu = 1
        sage: T = RealDistribution('t', nu)
        sage: s = T.get_random_element()
        sage: s.parent()
        Real Double Field
        sage: T.distribution_function(0)      # rel tol 1e-15
        0.3183098861837906
        sage: T.cum_distribution_function(1)  # rel tol 1e-15
        0.75
        sage: T.cum_distribution_function_inv(.5)
        0.0

    The F-distribution has two parameters ``nu1`` and ``nu2``::

        sage: nu1 = 9; nu2 = 17
        sage: F = RealDistribution('F', [nu1,nu2])
        sage: s = F.get_random_element()
        sage: s >= 0
        True
        sage: s.parent()
        Real Double Field
        sage: F.distribution_function(1)  # rel tol 1e-14
        0.6695025505192798
        sage: F.cum_distribution_function(3.68)  # rel tol 1e-14
        0.9899717772300652
        sage: F.cum_distribution_function_inv(0.99)  # rel tol 1e-14
        3.682241524045864

    The chi-squared distribution has one parameter ``nu``::

        sage: nu = 1
        sage: T = RealDistribution('chisquared', nu)
        sage: s = T.get_random_element()
        sage: s >= 0
        True
        sage: s.parent()
        Real Double Field
        sage: T.distribution_function(0)
        +infinity
        sage: T.cum_distribution_function(1)  # rel tol 1e-14
        0.6826894921370856
        sage: T.cum_distribution_function_inv(.5)  # rel tol 1e-14
        0.45493642311957305

    The exponential power distribution has two parameters ``a`` and
    ``b``::

        sage: a = 1
        sage: b = 2.5
        sage: T = RealDistribution('exppow', [a, b])
        sage: s = T.get_random_element()
        sage: s.parent()
        Real Double Field
        sage: T.distribution_function(0)  # rel tol 1e-14
        0.5635302489930136
        sage: T.cum_distribution_function(1)  # rel tol 1e-14
        0.940263052542855

    The beta distribution has two parameters ``a`` and ``b``::

        sage: a = 2
        sage: b = 2
        sage: T = RealDistribution('beta', [a, b])
        sage: s = T.get_random_element()
        sage: 0 <= s <= 1
        True
        sage: s.parent()
        Real Double Field
        sage: T.distribution_function(0)
        0.0
        sage: T.cum_distribution_function(1)
        1.0

    The exponential distribution has one parameter ``mu``::

        sage: mu = 2
        sage: T = RealDistribution('exponential', mu)
        sage: s = T.get_random_element()
        sage: 0 <= s
        True
        sage: s.parent()
        Real Double Field
        sage: T.distribution_function(0)
        0.5

    The gamma distribution has two parameters ``a`` and ``b``::

        sage: a = 2
        sage: b = 2
        sage: T = RealDistribution('gamma', [a, b])
        sage: s = T.get_random_element()
        sage: 0 <= s
        True
        sage: s.parent()
        Real Double Field
        sage: T.distribution_function(0)
        0.0

    The weibull distribution has two parameters ``a`` and ``b``::

        sage: a = 1
        sage: b = 1
        sage: T = RealDistribution('weibull', [a, b])
        sage: s = T.get_random_element()
        sage: s >= 0
        True
        sage: s.parent()
        Real Double Field
        sage: T.distribution_function(0)
        1.0
        sage: T.cum_distribution_function(1)
        0.6321205588285577
        sage: T.cum_distribution_function_inv(.5)
        0.6931471805599453

    It is possible to select which random number generator drives the
    sampling as well as the seed.  The default is the Mersenne
    twister. Also available are the RANDLXS algorithm and the
    Tausworthe generator (see the gsl reference manual for more
    details). These are all supposed to be simulation quality
    generators. For RANDLXS use ``rng='luxury'`` and for
    tausworth use ``rng='taus'``::

         sage: T = RealDistribution('gaussian', 1, rng='luxury', seed=10)

    To change the seed at a later time use ``set_seed``::

         sage: T.set_seed(100)

    TESTS:

    Make sure that repeated initializations are randomly seeded
    (:issue:`9770`)::

        sage: Xs = [RealDistribution('gaussian', 1).get_random_element() for _ in range(1000)]
        sage: len(set(Xs)) > 2^^32
        True
    """
    cdef gsl_rng_type *T
    cdef gsl_rng *r
    cdef int distribution_type
    cdef double* parameters
    cdef long int seed
    cdef object name
    # cdef double (*generator_1)(gsl_rng*)
    # cdef double (*generator_2)(gsl_rng*, double)
    # cdef _get_random_element_c(self)

    def __init__(self, type='uniform', parameters=None,
                 rng='default', seed=None):
        r"""
        EXAMPLES::

            sage: T = RealDistribution('gaussian', 1, seed=0)
            sage: T.get_random_element()  # rel tol 4e-16
            0.13391860811867587

        TESTS:

        Until :issue:`15089` a value of the ``seed`` keyword
        besides ``None`` was ignored. We check here that setting
        a seed is effective. ::

            sage: T = RealDistribution("beta",[1.6,4.3], seed=876)
            sage: one = [T.get_random_element() for _ in range(10)]
            sage: T = RealDistribution("beta",[1.6,4.3], seed=876)
            sage: two = [T.get_random_element() for _ in range(10)]
            sage: T = RealDistribution("beta",[1.6,4.3], seed=123)
            sage: three = [T.get_random_element() for _ in range(10)]
            sage: one == two
            True
            sage: one == three
            False
        """
        gsl_rng_env_setup()
        self.parameters = NULL
        self.set_random_number_generator(rng)
        self.r = gsl_rng_alloc(self.T)
        if seed is None:
            seed = random.randint(1, 2**31)
        self.set_seed(seed)
        self.name = " "
        self.set_distribution(type, parameters)

    def set_seed(self, seed):
        r"""
        Set the seed for the underlying random number generator.

        EXAMPLES::

            sage: T = RealDistribution('gaussian', 1, rng='luxury', seed=10)
            sage: T.set_seed(100)
        """

        gsl_rng_set(self.r, seed)
        self.seed = seed

    def set_random_number_generator(self, rng='default'):
        r"""
        Set the gsl random number generator to be one of ``'default'``,
        ``'luxury'``, or ``'taus'``.

        EXAMPLES::

            sage: T = SphericalDistribution()
            sage: T.set_random_number_generator('default')
            sage: T.set_seed(0)
            sage: T.get_random_element()  # rel tol 4e-16
            (0.07961564104639995, -0.05237671627581255, 0.9954486572862178)
            sage: T.set_random_number_generator('luxury')
            sage: T.set_seed(0)
            sage: T.get_random_element()  # rel tol 4e-16
            (0.07961564104639995, -0.05237671627581255, 0.9954486572862178)
        """
        if rng == 'default':
            self.T = gsl_rng_default
        elif rng == 'luxury':
            self.T = gsl_rng_ranlxd2
        elif rng == 'taus':
            self.T = gsl_rng_taus2
        else:
            raise TypeError("Not a valid random number generator")

    def __dealloc__(self):
        if self.r != NULL:
            gsl_rng_free(self.r)
        sig_free(self.parameters)

    def __str__(self):
        r"""
        Return the name of the current distribution.

        EXAMPLES::

            sage: T = RealDistribution('gaussian', 1)
            sage: str(T)
            'gaussian'
            sage: T = RealDistribution('beta', [2, 2])
            sage: str(T)
            'beta'
        """
        return self.name

    def get_random_element(self):
        r"""
        Get a random sample from the probability distribution.

        EXAMPLES::

            sage: T = RealDistribution('gaussian', 1, seed=0)
            sage: T.get_random_element()  # rel tol 4e-16
            0.13391860811867587
        """
        cdef double result
        if self.distribution_type == uniform:
            result = gsl_ran_flat(self.r, self.parameters[0], self.parameters[1])
#            result = gsl_rng_uniform(self.r)
        elif self.distribution_type == gaussian:
            result = gsl_ran_gaussian(self.r, self.parameters[0])
        elif self.distribution_type == rayleigh:
            result = gsl_ran_rayleigh(self.r, self.parameters[0])
        elif self.distribution_type == lognormal:
            result = gsl_ran_lognormal(self.r, self.parameters[0], self.parameters[1])
        elif self.distribution_type == pareto:
            result = gsl_ran_pareto(self.r, self.parameters[0], self.parameters[1])
        elif self.distribution_type == t:
            result = gsl_ran_tdist(self.r, self.parameters[0])
        elif self.distribution_type == F:
            result = gsl_ran_fdist(self.r, self.parameters[0], self.parameters[1])
        elif self.distribution_type == chisquared:
            result = gsl_ran_chisq(self.r, self.parameters[0])
        elif self.distribution_type == exppow:
            result = gsl_ran_exppow(self.r, self.parameters[0], self.parameters[1])
        elif self.distribution_type == weibull:
            result = gsl_ran_weibull(self.r, self.parameters[0], self.parameters[1])
        elif self.distribution_type == beta:
            result = gsl_ran_beta(self.r, self.parameters[0], self.parameters[1])
        elif self.distribution_type == exponential:
            result = gsl_ran_exponential(self.r, self.parameters[0])
        elif self.distribution_type == gamma:
            result = gsl_ran_gamma(self.r, self.parameters[0], self.parameters[1])
        else:
            raise TypeError("Not a supported probability distribution")

        return sage.rings.real_double.RDF(result)

    def set_distribution(self, name='uniform', parameters=None):
        r"""
        This method can be called to change the current probability distribution.

        EXAMPLES::

            sage: T = RealDistribution('gaussian', 1)
            sage: T.set_distribution('gaussian', 1)
            sage: T.set_distribution('pareto', [0, 1])
        """
        sig_free(self.parameters)

        if parameters is None:
            parameters = []

        if name == 'uniform':
            self.distribution_type = uniform
            for x in parameters:
                try:
                    float(x)
                except Exception:
                    raise TypeError("Uniform distribution requires parameters coercible to float")
            self.parameters = <double*>sig_malloc(sizeof(double)*2)
            self.parameters[0] = parameters[0]
            self.parameters[1] = parameters[1]
        elif name == 'gaussian':
            try:
                float(parameters)
            except Exception:
                raise TypeError("gaussian distribution requires parameter sigma coercible to float")
            self.parameters = <double*>sig_malloc(sizeof(double))
            self.parameters[0] = float(parameters)
            self.distribution_type = gaussian
        elif name == 'pareto':
            if len(parameters) != 2:
                raise TypeError("pareto distribution has two parameters")
            try:
                map(float, parameters)
            except Exception:
                raise TypeError("parameters must be coercible to float")
            self.parameters = <double*>sig_malloc(sizeof(double)*2)
            self.parameters[0] = float(parameters[0])
            self.parameters[1] = float(parameters[1])
            self.distribution_type = pareto
        elif name == 'rayleigh':
            try:
                float(parameters)
            except Exception:
                raise TypeError("rayleigh distribution requires parameter sigma coercible to float")
            self.parameters = <double*>sig_malloc(sizeof(double))
            self.parameters[0] = float(parameters)
            self.distribution_type = rayleigh
        elif name == 'lognormal':
            if len(parameters) != 2:
                raise TypeError("Lognormal distribution requires two parameters")
            for x in parameters:
                try:
                    float(x)
                except Exception:
                    raise TypeError("Lognormal distribution requires real parameters")
            self.parameters = <double*>sig_malloc(sizeof(double)*2)
            self.parameters[0] = float(parameters[0])
            self.parameters[1] = float(parameters[1])
            self.distribution_type = lognormal
        elif name == 't':
            try:
                float(parameters)
            except Exception:
                raise TypeError("parameter to t distribution must be coercible to float")
            self.parameters = <double*>sig_malloc(sizeof(double))
            self.parameters[0] = float(parameters)
            self.distribution_type = t
        elif name == 'F':
            if len(parameters) != 2:
                raise TypeError("F-distribution requires two real parameters")
            try:
                map(float, parameters)
            except Exception:
                raise TypeError("F-distribution requires real parameters")
            self.parameters = <double *>sig_malloc(sizeof(double)*2)
            self.parameters[0] = float(parameters[0])
            self.parameters[1] = float(parameters[1])
            self.distribution_type = F
        elif name == 'chisquared':
            try:
                float(parameters)
            except Exception:
                raise TypeError("parameters to t distribution must be coercible to float")
            self.parameters = <double *>sig_malloc(sizeof(double))
            self.parameters[0] = float(parameters)
            self.distribution_type = chisquared
        elif name == 'exppow':
            if len(parameters) != 2:
                raise TypeError("exponential power distribution requires two parameters")
            for x in parameters:
                try:
                    float(x)
                except Exception:
                    raise TypeError("exponential power distribution requires real parameters")
            self.parameters = <double*>sig_malloc(sizeof(double)*2)
            self.parameters[0] = float(parameters[0])
            self.parameters[1] = float(parameters[1])
            self.distribution_type = exppow
        elif name == 'weibull':
            if len(parameters) != 2:
                raise TypeError("weibull distribution requires two real parameters")
            try:
                map(float, parameters)
            except Exception:
                raise TypeError("weibull distribution requires real parameters")
            self.parameters = <double *>sig_malloc(sizeof(double)*2)
            self.parameters[0] = float(parameters[0])
            self.parameters[1] = float(parameters[1])
            self.distribution_type = weibull
        elif name == 'beta':
            if len(parameters) != 2:
                raise TypeError("beta distribution requires two real parameters")
            try:
                map(float, parameters)
            except Exception:
                raise TypeError("beta distribution requires real parameters")
            self.parameters = <double *>sig_malloc(sizeof(double)*2)
            self.parameters[0] = float(parameters[0])
            self.parameters[1] = float(parameters[1])
            self.distribution_type = beta
        elif name == 'exponential':
            try:
                float(parameters)
            except Exception:
                raise TypeError("exponential distribution requires parameter mu coercible to float")
            self.parameters = <double*>sig_malloc(sizeof(double))
            self.parameters[0] = float(parameters)
            self.distribution_type = exponential
        elif name == 'gamma':
            if len(parameters) != 2:
                raise TypeError("gamma distribution requires two real parameters")
            try:
                map(float, parameters)
            except Exception:
                raise TypeError("gamma distribution requires real parameters")
            self.parameters = <double *>sig_malloc(sizeof(double)*2)
            self.parameters[0] = float(parameters[0])
            self.parameters[1] = float(parameters[1])
            self.distribution_type = gamma
        else:
            raise TypeError("Not a supported probability distribution")

        self.name = name

    # def _get_random_element_c():

    def reset_distribution(self):
        r"""
        Reset the distribution.

        EXAMPLES::

            sage: T = RealDistribution('gaussian', 1, seed=10)
            sage: [T.get_random_element() for _ in range(10)]  # rel tol 4e-16
            [-0.7460999595745819, -0.004644606626413462, -0.8720538317207641, 0.6916259921666037, 2.67668674666043, 0.6325002813661014, -0.7974263521959355, -0.5284976893366636, 1.1353119849528792, 0.9912505673230749]
            sage: T.reset_distribution()
            sage: [T.get_random_element() for _ in range(10)]  # rel tol 4e-16
            [-0.7460999595745819, -0.004644606626413462, -0.8720538317207641, 0.6916259921666037, 2.67668674666043, 0.6325002813661014, -0.7974263521959355, -0.5284976893366636, 1.1353119849528792, 0.9912505673230749]
        """
        if self.r != NULL:
            gsl_rng_free(self.r)
        self.r = gsl_rng_alloc(self.T)
        self.set_seed(self.seed)
#        gsl_rng_env_setup()

    def distribution_function(self, x):
        r"""
        Evaluate the distribution function of the
        probability distribution at ``x``.

        EXAMPLES::

            sage: T = RealDistribution('uniform', [0, 2])
            sage: T.distribution_function(0)
            0.5
            sage: T.distribution_function(1)
            0.5
            sage: T.distribution_function(1.5)
            0.5
            sage: T.distribution_function(2)
            0.0
        """
        if self.distribution_type == uniform:
            return sage.rings.real_double.RDF(gsl_ran_flat_pdf(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == gaussian:
            return sage.rings.real_double.RDF(gsl_ran_gaussian_pdf(x, self.parameters[0]))
        elif self.distribution_type == rayleigh:
            return sage.rings.real_double.RDF(gsl_ran_rayleigh_pdf(x, self.parameters[0]))
        elif self.distribution_type == lognormal:
            return sage.rings.real_double.RDF(gsl_ran_lognormal_pdf(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == pareto:
            return sage.rings.real_double.RDF(gsl_ran_pareto_pdf(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == t:
            return sage.rings.real_double.RDF(gsl_ran_tdist_pdf(x, self.parameters[0]))
        elif self.distribution_type == F:
            return sage.rings.real_double.RDF(gsl_ran_fdist_pdf(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == chisquared:
            return sage.rings.real_double.RDF(gsl_ran_chisq_pdf(x, self.parameters[0]))
        elif self.distribution_type == exppow:
            return sage.rings.real_double.RDF(gsl_ran_exppow_pdf(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == weibull:
            return sage.rings.real_double.RDF(gsl_ran_weibull_pdf(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == beta:
            return sage.rings.real_double.RDF(gsl_ran_beta_pdf(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == exponential:
            return sage.rings.real_double.RDF(gsl_ran_exponential_pdf(x, self.parameters[0]))
        elif self.distribution_type == gamma:
            return sage.rings.real_double.RDF(gsl_ran_gamma_pdf(x, self.parameters[0], self.parameters[1]))
        else:
            raise TypeError("Not a supported probability distribution")

    def cum_distribution_function(self, x):
        r"""
        Evaluate the cumulative distribution function of
        the probability distribution at ``x``.

        EXAMPLES::

            sage: T = RealDistribution('uniform', [0, 2])
            sage: T.cum_distribution_function(1)
            0.5
        """
        if self.distribution_type == uniform:
            return sage.rings.real_double.RDF(gsl_cdf_flat_P(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == gaussian:
            return sage.rings.real_double.RDF(gsl_cdf_gaussian_P(x, self.parameters[0]))
        elif self.distribution_type == rayleigh:
            return sage.rings.real_double.RDF(gsl_cdf_rayleigh_P(x, self.parameters[0]))
        elif self.distribution_type == lognormal:
            return sage.rings.real_double.RDF(gsl_cdf_lognormal_P(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == pareto:
            return sage.rings.real_double.RDF(gsl_cdf_pareto_P(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == t:
            return sage.rings.real_double.RDF(gsl_cdf_tdist_P(x, self.parameters[0]))
        elif self.distribution_type == F:
            return sage.rings.real_double.RDF(gsl_cdf_fdist_P(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == chisquared:
            return sage.rings.real_double.RDF(gsl_cdf_chisq_P(x, self.parameters[0]))
        elif self.distribution_type == exppow:
            return sage.rings.real_double.RDF(gsl_cdf_exppow_P(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == weibull:
            return sage.rings.real_double.RDF(gsl_cdf_weibull_P(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == beta:
            return sage.rings.real_double.RDF(gsl_cdf_beta_P(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == exponential:
            return sage.rings.real_double.RDF(gsl_cdf_exponential_P(x, self.parameters[0]))
        elif self.distribution_type == gamma:
            return sage.rings.real_double.RDF(gsl_cdf_gamma_P(x, self.parameters[0], self.parameters[1]))
        else:
            raise TypeError("Not a supported probability distribution")

    def cum_distribution_function_inv(self, x):
        r"""
        Evaluate the inverse of the cumulative distribution
        distribution function of the probability distribution at ``x``.

        EXAMPLES::

            sage: T = RealDistribution('uniform', [0, 2])
            sage: T.cum_distribution_function_inv(.5)
            1.0
        """
        if self.distribution_type == uniform:
            return sage.rings.real_double.RDF(gsl_cdf_flat_Pinv(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == gaussian:
            return sage.rings.real_double.RDF(gsl_cdf_gaussian_Pinv(x, self.parameters[0]))
        elif self.distribution_type == rayleigh:
            return sage.rings.real_double.RDF(gsl_cdf_rayleigh_Pinv(x, self.parameters[0]))
        elif self.distribution_type == lognormal:
            return sage.rings.real_double.RDF(gsl_cdf_lognormal_Pinv(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == pareto:
            return sage.rings.real_double.RDF(gsl_cdf_pareto_Pinv(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == t:
            return sage.rings.real_double.RDF(gsl_cdf_tdist_Pinv(x, self.parameters[0]))
        elif self.distribution_type == F:
            return sage.rings.real_double.RDF(gsl_cdf_fdist_Pinv(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == chisquared:
            return sage.rings.real_double.RDF(gsl_cdf_chisq_Pinv(x, self.parameters[0]))
        elif self.distribution_type == exppow:
            raise NotImplementedError("gsl does not provide inverse for exponential power")
#            return sage.rings.real_double.RDF(gsl_cdf_exppow_Pinv(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == weibull:
            return sage.rings.real_double.RDF(gsl_cdf_weibull_Pinv(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == beta:
            return sage.rings.real_double.RDF(gsl_cdf_beta_Pinv(x, self.parameters[0], self.parameters[1]))
        elif self.distribution_type == exponential:
            return sage.rings.real_double.RDF(gsl_cdf_exponential_Pinv(x, self.parameters[0]))
        elif self.distribution_type == gamma:
            return sage.rings.real_double.RDF(gsl_cdf_gamma_Pinv(x, self.parameters[0], self.parameters[1]))
        else:
            raise TypeError("Not a supported probability distribution")

    def plot(self, *args, **kwds):
        r"""
        Plot the distribution function for the probability
        distribution. Parameters to :func:`sage.plot.plot.plot` can be
        passed through ``*args`` and ``**kwds``.

        EXAMPLES::

            sage: T = RealDistribution('uniform', [0, 2])
            sage: P = T.plot()                                                          # needs sage.plot
        """
        from sage.plot.plot import plot
        return plot(self.distribution_function, *args, **kwds)


cdef class GeneralDiscreteDistribution(ProbabilityDistribution):
    r"""
    Create a discrete probability distribution.

    INPUT:

    - ``P`` -- list of probabilities; the list will automatically be
      normalised if ``sum(P)`` is not equal to 1

    - ``rng`` -- (optional) random number generator to use; may be
      one of ``'default'``, ``'luxury'``, or ``'taus'``

    - ``seed`` -- (optional) seed to use with the random number
      generator

    OUTPUT: a probability distribution where the probability of selecting
    ``x`` is ``P[x]``.

    EXAMPLES:

    Construct a ``GeneralDiscreteDistribution`` with the probability
    distribution `P` where `P(0) = 0.3`, `P(1) = 0.4`, `P(2) = 0.3`::

        sage: P = [0.3, 0.4, 0.3]
        sage: X = GeneralDiscreteDistribution(P)
        sage: X.get_random_element() in (0, 1, 2)
        True

    Checking the distribution of samples::

        sage: P = [0.3, 0.4, 0.3]
        sage: counts = [0] * len(P)
        sage: X = GeneralDiscreteDistribution(P)
        sage: nr_samples = 10000
        sage: for _ in range(nr_samples):
        ....:     counts[X.get_random_element()] += 1
        sage: [1.0*x/nr_samples for x in counts]  # abs tol 3e-2
        [0.3, 0.4, 0.3]

    The distribution probabilities will automatically be normalised::

        sage: P = [0.1, 0.3]
        sage: X = GeneralDiscreteDistribution(P, seed=0)
        sage: counts = [0, 0]
        sage: for _ in range(10000):
        ....:     counts[X.get_random_element()] += 1
        sage: float(counts[1]/counts[0])
        3.042037186742118

    TESTS:

    Make sure that repeated initializations are randomly seeded
    (:issue:`9770`)::

        sage: P = [0.001] * 1000
        sage: Xs = [GeneralDiscreteDistribution(P).get_random_element() for _ in range(1000)]
        sage: len(set(Xs)) > 2^^32
        True

    The distribution probabilities must be nonnegative::

        sage: GeneralDiscreteDistribution([0.1, -0.1])
        Traceback (most recent call last):
        ...
        ValueError: The distribution probabilities must be nonnegative
    """
    cdef gsl_rng_type * T
    cdef gsl_rng * r
    cdef gsl_ran_discrete_t *dist
    cdef long seed

    def __init__(self, P, rng='default', seed=None):
        r"""
        Given a list of probabilities P construct an instance of a gsl
        discrete random variable generator.

        EXAMPLES::

            sage: P = [0.3, 0.4, 0.3]
            sage: X = GeneralDiscreteDistribution(P)
            sage: assert X.get_random_element() in range(len(P))

        TESTS:

        Until :issue:`15089` a value of the ``seed`` keyword
        besides ``None`` was ignored. We check here that setting
        a seed is effective. ::

            sage: P = [0.2, 0.3, 0.1, 0.4]
            sage: T = GeneralDiscreteDistribution(P, seed=876)
            sage: one = [T.get_random_element() for _ in range(50)]
            sage: T = GeneralDiscreteDistribution(P, seed=876)
            sage: two = [T.get_random_element() for _ in range(50)]
            sage: T = GeneralDiscreteDistribution(P, seed=123)
            sage: three = [T.get_random_element() for _ in range(50)]
            sage: one == two
            True
            sage: one == three
            False

        Testing that :issue:`24416` is fixed for when entries are larger
        than `2^{1024}`::

            sage: from collections import Counter
            sage: X = GeneralDiscreteDistribution([1,2,2^1024])
            sage: Counter(X.get_random_element() for _ in range(100))
            Counter({2: 100})
        """
        gsl_rng_env_setup()
        self.set_random_number_generator(rng)
        self.r = gsl_rng_alloc(self.T)
        if seed is None:
            seed = random.randint(1, 2**31)
        self.set_seed(seed)

        cdef int n
        n = len(P)

        s = sum(P)
        if s != 1:
            P = [p/s for p in P]

        cdef double *P_vec
        P_vec = <double *> sig_malloc(n*(sizeof(double)))

        cdef int i
        for i in range(n):
            if P[i] < 0:
                raise ValueError("The distribution probabilities must "
                                 "be nonnegative")
            P_vec[i] = P[i]

        self.dist = gsl_ran_discrete_preproc(n, P_vec)

        sig_free(P_vec)

    def set_seed(self, seed):
        r"""
        Set the seed to be used by the random number generator.

        EXAMPLES::

            sage: X = GeneralDiscreteDistribution([0.3, 0.4, 0.3])
            sage: X.set_seed(1)
            sage: X.get_random_element()
            1
        """
        gsl_rng_set(self.r, seed)
        self.seed = seed

    def set_random_number_generator(self, rng='default'):
        r"""
        Set the random number generator to be used by gsl.

        EXAMPLES::

            sage: X = GeneralDiscreteDistribution([0.3, 0.4, 0.3])
            sage: X.set_random_number_generator('taus')
        """
        if rng == 'default':
            self.T = gsl_rng_default
        elif rng == 'luxury':
            self.T = gsl_rng_ranlxd2
        elif rng == 'taus':
            self.T = gsl_rng_taus2
        else:
            raise TypeError("Not a valid random number generator")

    def __dealloc__(self):
        if self.r != NULL:
            gsl_rng_free(self.r)

        if self.dist != NULL:
            gsl_ran_discrete_free(self.dist)

    def get_random_element(self):
        r"""
        Get a random sample from the probability distribution.

        EXAMPLES::

            sage: P = [0.3, 0.4, 0.3]
            sage: X = GeneralDiscreteDistribution(P)
            sage: all(X.get_random_element() in (0,1,2) for _ in range(10))
            True
            sage: isinstance(X.get_random_element(), sage.rings.integer.Integer)
            True
        """
        return sage.rings.integer.Integer(gsl_ran_discrete(self.r, self.dist))

    def reset_distribution(self):
        r"""
        This method resets the distribution.

        EXAMPLES::

            sage: T = GeneralDiscreteDistribution([0.1, 0.3, 0.6])
            sage: T.set_seed(0)
            sage: [T.get_random_element() for _ in range(10)]
            [2, 2, 2, 2, 2, 1, 2, 2, 1, 2]
            sage: T.reset_distribution()
            sage: [T.get_random_element() for _ in range(10)]
            [2, 2, 2, 2, 2, 1, 2, 2, 1, 2]
        """
        if self.r != NULL:
            gsl_rng_free(self.r)
        self.r = gsl_rng_alloc(self.T)
        self.set_seed(self.seed)
