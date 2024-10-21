# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.symbolic
#
# distutils: sources = sage/stats/distributions/dgs_gauss_mp.c sage/stats/distributions/dgs_gauss_dp.c sage/stats/distributions/dgs_bern.c
# distutils: depends = sage/stats/distributions/dgs_gauss.h sage/stats/distributions/dgs_bern.h sage/stats/distributions/dgs_misc.h
# distutils: extra_compile_args = -D_XOPEN_SOURCE=600

r"""
Discrete Gaussian Samplers over the Integers

This class realizes oracles which returns integers proportionally to
`\exp(-(x-c)^2/(2σ^2))`. All oracles are implemented using rejection sampling.
See :func:`DiscreteGaussianDistributionIntegerSampler.__init__` for which algorithms are
available.

AUTHORS:

- Martin Albrecht (2014-06-28): initial version

EXAMPLES:

We construct a sampler for the distribution `D_{3,c}` with width `σ=3` and center `c=0`::

    sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
    sage: sigma = 3.0
    sage: D = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)

We ask for 100000 samples::

    sage: from collections import defaultdict
    sage: counter = defaultdict(Integer)
    sage: n = 0
    sage: def add_samples(i):
    ....:     global counter, n
    ....:     for _ in range(i):
    ....:         counter[D()] += 1
    ....:         n += 1

    sage: add_samples(100000)

These are sampled with a probability proportional to `\exp(-x^2/18)`. More
precisely we have to normalise by dividing by the overall probability over all
integers. We use the fact that hitting anything more than 6 standard deviations
away is very unlikely and compute::

    sage: bound = (6*sigma).floor()
    sage: norm_factor = sum([exp(-x^2/(2*sigma^2)) for x in range(-bound,bound+1)])
    sage: norm_factor
    7.519...

With this normalisation factor, we can now test if our samples follow the
expected distribution::

    sage: expected = lambda x : ZZ(round(n*exp(-x^2/(2*sigma^2))/norm_factor))
    sage: observed = lambda x : counter[x]

    sage: add_samples(10000)
    sage: while abs(observed(0)*1.0/expected(0)   - 1.0) > 5e-2: add_samples(10000)
    sage: while abs(observed(4)*1.0/expected(4)   - 1.0) > 5e-2: add_samples(10000)
    sage: while abs(observed(-10)*1.0/expected(-10) - 1.0) > 5e-2: add_samples(10000)  # long time

We construct an instance with a larger width::

    sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
    sage: sigma = 127
    sage: D = DiscreteGaussianDistributionIntegerSampler(sigma=sigma, algorithm='uniform+online')

ask for 100000 samples::

    sage: from collections import defaultdict
    sage: counter = defaultdict(Integer)
    sage: n = 0
    sage: def add_samples(i):
    ....:     global counter, n
    ....:     for _ in range(i):
    ....:         counter[D()] += 1
    ....:         n += 1

    sage: add_samples(100000)

and check if the proportions fit::

    sage: expected = lambda x, y: (
    ....:     exp(-x^2/(2*sigma^2))/exp(-y^2/(2*sigma^2)).n())
    sage: observed = lambda x, y: float(counter[x])/counter[y]

    sage: while not all(v in counter for v in (0, 1, -100)): add_samples(10000)

    sage: while abs(expected(0, 1) - observed(0, 1)) > 2e-1: add_samples(10000)
    sage: while abs(expected(0, -100) - observed(0, -100)) > 2e-1: add_samples(10000)

We construct a sampler with `c\%1 != 0`::

    sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
    sage: sigma = 3
    sage: D = DiscreteGaussianDistributionIntegerSampler(sigma=sigma, c=1/2)
    sage: s = 0
    sage: n = 0
    sage: def add_samples(i):
    ....:     global s, n
    ....:     for _ in range(i):
    ....:         s += D()
    ....:         n += 1
    ....:
    sage: add_samples(100000)
    sage: while abs(float(s)/n - 0.5) > 5e-2: add_samples(10000)

REFERENCES:

- [DDLL2013]_
"""
#******************************************************************************
#
#                        DGS - Discrete Gaussian Samplers
#
# Copyright (c) 2014, Martin Albrecht  <martinralbrecht+dgs@googlemail.com>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# The views and conclusions contained in the software and documentation are
# those of the authors and should not be interpreted as representing official
# policies, either expressed or implied, of the FreeBSD Project.
#*****************************************************************************/

from cysignals.signals cimport sig_on, sig_off

from sage.rings.real_mpfr cimport RealNumber, RealField
from sage.libs.mpfr cimport mpfr_set, MPFR_RNDN
from sage.rings.integer cimport Integer
from sage.misc.randstate cimport randstate, current_randstate

from sage.stats.distributions.dgs cimport dgs_disc_gauss_mp_init, dgs_disc_gauss_mp_clear, dgs_disc_gauss_mp_flush_cache
from sage.stats.distributions.dgs cimport dgs_disc_gauss_dp_init, dgs_disc_gauss_dp_clear, dgs_disc_gauss_dp_flush_cache
from sage.stats.distributions.dgs cimport DGS_DISC_GAUSS_UNIFORM_TABLE, DGS_DISC_GAUSS_UNIFORM_ONLINE, DGS_DISC_GAUSS_UNIFORM_LOGTABLE, DGS_DISC_GAUSS_SIGMA2_LOGTABLE

cdef class DiscreteGaussianDistributionIntegerSampler(SageObject):
    r"""
    A Discrete Gaussian Sampler using rejection sampling.

    .. automethod:: __init__
    .. automethod:: __call__
    """

    # We use tables for σt ≤ table_cutoff
    table_cutoff = 10**6

    def __init__(self, sigma, c=0, tau=6, algorithm=None, precision='mp'):
        r"""
        Construct a new sampler for a discrete Gaussian distribution.

        INPUT:

        - ``sigma`` -- samples `x` are accepted with probability proportional to
          `\exp(-(x-c)²/(2σ²))`

        - ``c`` -- the mean of the distribution. The value of ``c`` does not have
          to be an integer. However, some algorithms only support integer-valued
          ``c`` (default: ``0``)

        - ``tau`` -- samples outside the range `(⌊c⌉-⌈στ⌉,...,⌊c⌉+⌈στ⌉)` are
          considered to have probability zero. This bound applies to algorithms which
          sample from the uniform distribution (default: ``6``)

        - ``algorithm`` -- see list below (default: ``'uniform+table'`` for
           `σt` bounded by ``DiscreteGaussianDistributionIntegerSampler.table_cutoff`` and
           ``'uniform+online'`` for bigger `στ`)

        - ``precision`` -- either ``'mp'`` for multi-precision where the actual
          precision used is taken from sigma or ``'dp'`` for double precision. In
          the latter case results are not reproducible. (default: ``'mp'``)

        ALGORITHMS:

        - ``'uniform+table'`` -- classical rejection sampling, sampling from the
          uniform distribution and accepted with probability proportional to
          `\exp(-(x-c)²/(2σ²))` where `\exp(-(x-c)²/(2σ²))` is precomputed and
          stored in a table. Any real-valued `c` is supported.

        - ``'uniform+logtable'`` -- samples are drawn from a uniform distribution and
          accepted with probability proportional to `\exp(-(x-c)²/(2σ²))` where
          `\exp(-(x-c)²/(2σ²))` is computed using logarithmically many calls to
          Bernoulli distributions. See [DDLL2013]_ for details.  Only
          integer-valued `c` are supported.

        - ``'uniform+online'`` -- samples are drawn from a uniform distribution and
          accepted with probability proportional to `\exp(-(x-c)²/(2σ²))` where
          `\exp(-(x-c)²/(2σ²))` is computed in each invocation. Typically this
          is very slow.  See [DDLL2013]_ for details.  Any real-valued `c` is
          accepted.

        - ``'sigma2+logtable'`` -- samples are drawn from an easily samplable
          distribution with `σ = k·σ_2` with `σ_2 = \sqrt{1/(2\log 2)}` and accepted
          with probability proportional to `\exp(-(x-c)²/(2σ²))` where
          `\exp(-(x-c)²/(2σ²))` is computed using  logarithmically many calls to Bernoulli
          distributions (but no calls to `\exp`). See [DDLL2013]_ for details. Note that this
          sampler adjusts `σ` to match `k·σ_2` for some integer `k`.
          Only integer-valued `c` are supported.

        EXAMPLES::

            sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
            sage: DiscreteGaussianDistributionIntegerSampler(3.0, algorithm='uniform+online')
            Discrete Gaussian sampler over the Integers with sigma = 3.000000 and c = 0.000000
            sage: DiscreteGaussianDistributionIntegerSampler(3.0, algorithm='uniform+table')
            Discrete Gaussian sampler over the Integers with sigma = 3.000000 and c = 0.000000
            sage: DiscreteGaussianDistributionIntegerSampler(3.0, algorithm='uniform+logtable')
            Discrete Gaussian sampler over the Integers with sigma = 3.000000 and c = 0.000000

        Note that ``'sigma2+logtable'`` adjusts `σ`::

            sage: DiscreteGaussianDistributionIntegerSampler(3.0, algorithm='sigma2+logtable')
            Discrete Gaussian sampler over the Integers with sigma = 3.397287 and c = 0.000000

        TESTS:

        We are testing invalid inputs::

            sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
            sage: DiscreteGaussianDistributionIntegerSampler(-3.0)
            Traceback (most recent call last):
            ...
            ValueError: sigma must be > 0.0 but got -3.000000

            sage: DiscreteGaussianDistributionIntegerSampler(3.0, tau=-1)
            Traceback (most recent call last):
            ...
            ValueError: tau must be >= 1 but got -1

            sage: DiscreteGaussianDistributionIntegerSampler(3.0, tau=2, algorithm='superfastalgorithmyouneverheardof')
            Traceback (most recent call last):
            ...
            ValueError: Algorithm 'superfastalgorithmyouneverheardof' not supported by class 'DiscreteGaussianDistributionIntegerSampler'

            sage: DiscreteGaussianDistributionIntegerSampler(3.0, c=1.5, algorithm='sigma2+logtable')
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'uniform+logtable' requires c%1 == 0

        We are testing correctness for multi-precision::

            sage: def add_samples(i):
            ....:     global mini, maxi, s, n
            ....:     for _ in range(i):
            ....:         x = D()
            ....:         s += x
            ....:         maxi = max(maxi, x)
            ....:         mini = min(mini, x)
            ....:         n += 1

            sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
            sage: D = DiscreteGaussianDistributionIntegerSampler(1.0, c=0, tau=2)
            sage: mini = 1000; maxi = -1000; s = 0; n = 0
            sage: add_samples(2^16)
            sage: while mini != 0 - 2*1.0 or maxi != 0 + 2*1.0 or abs(float(s)/n) >= 0.01:
            ....:     add_samples(2^16)

            sage: D = DiscreteGaussianDistributionIntegerSampler(1.0, c=2.5, tau=2)
            sage: mini = 1000; maxi = -1000; s = 0; n = 0
            sage: add_samples(2^16)
            sage: while mini != 2 - 2*1.0 or maxi != 2 + 2*1.0 or abs(float(s)/n - 2.45) >= 0.01:
            ....:     add_samples(2^16)

            sage: D = DiscreteGaussianDistributionIntegerSampler(1.0, c=2.5, tau=6)
            sage: mini = 1000; maxi = -1000; s = 0; n = 0
            sage: add_samples(2^18)
            sage: while mini > 2 - 4*1.0 or maxi < 2 + 5*1.0 or abs(float(s)/n - 2.5) >= 0.01:  # long time
            ....:     add_samples(2^18)

        We are testing correctness for double precision::

            sage: def add_samples(i):
            ....:     global mini, maxi, s, n
            ....:     for _ in range(i):
            ....:         x = D()
            ....:         s += x
            ....:         maxi = max(maxi, x)
            ....:         mini = min(mini, x)
            ....:         n += 1

            sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
            sage: D = DiscreteGaussianDistributionIntegerSampler(1.0, c=0, tau=2, precision='dp')
            sage: mini = 1000; maxi = -1000; s = 0; n = 0
            sage: add_samples(2^16)
            sage: while mini != 0 - 2*1.0 or maxi != 0 + 2*1.0 or abs(float(s)/n) >= 0.05:
            ....:     add_samples(2^16)

            sage: D = DiscreteGaussianDistributionIntegerSampler(1.0, c=2.5, tau=2, precision='dp')
            sage: mini = 1000; maxi = -1000; s = 0; n = 0
            sage: add_samples(2^16)
            sage: while mini != 2 - 2*1.0 or maxi != 2 + 2*1.0 or abs(float(s)/n - 2.45) >= 0.01:
            ....:     add_samples(2^16)

            sage: D = DiscreteGaussianDistributionIntegerSampler(1.0, c=2.5, tau=6, precision='dp')
            sage: mini = 1000; maxi = -1000; s = 0; n = 0
            sage: add_samples(2^16)
            sage: while mini > -1 or maxi < 6 or abs(float(s)/n - 2.5) >= 0.1:
            ....:     add_samples(2^16)

        We plot a histogram::

            sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
            sage: D = DiscreteGaussianDistributionIntegerSampler(17.0)
            sage: S = [D() for _ in range(2^16)]
            sage: list_plot([(v,S.count(v)) for v in set(S)])  # long time
            Graphics object consisting of 1 graphics primitive

        These generators cache random bits for performance reasons. Hence, resetting
        the seed of the PRNG might not have the expected outcome. You can flush this cache with
        ``_flush_cache()``::

            sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
            sage: D = DiscreteGaussianDistributionIntegerSampler(3.0)
            sage: sage.misc.randstate.set_random_seed(0); D()
            3
            sage: sage.misc.randstate.set_random_seed(0); D()
            3
            sage: sage.misc.randstate.set_random_seed(0); D._flush_cache(); D()
            3

            sage: D = DiscreteGaussianDistributionIntegerSampler(3.0)
            sage: sage.misc.randstate.set_random_seed(0); D()
            3
            sage: sage.misc.randstate.set_random_seed(0); D()
            3
            sage: sage.misc.randstate.set_random_seed(0); D()
            -3
        """
        if sigma <= 0.0:
            raise ValueError("sigma must be > 0.0 but got %f"%sigma)

        if tau < 1:
            raise ValueError("tau must be >= 1 but got %d"%tau)

        if algorithm is None:
            if sigma*tau <= DiscreteGaussianDistributionIntegerSampler.table_cutoff:
                algorithm = "uniform+table"
            else:
                algorithm = "uniform+online"

        algorithm_str = algorithm

        if algorithm == "uniform+table":
            algorithm = DGS_DISC_GAUSS_UNIFORM_TABLE
        elif algorithm == "uniform+online":
            algorithm = DGS_DISC_GAUSS_UNIFORM_ONLINE
        elif algorithm == "uniform+logtable":
            if (c%1):
                raise ValueError("algorithm 'uniform+logtable' requires c%1 == 0")
            algorithm = DGS_DISC_GAUSS_UNIFORM_LOGTABLE
        elif algorithm == "sigma2+logtable":
            if (c%1):
                raise ValueError("algorithm 'uniform+logtable' requires c%1 == 0")
            algorithm = DGS_DISC_GAUSS_SIGMA2_LOGTABLE
        else:
            raise ValueError("Algorithm '%s' not supported by class 'DiscreteGaussianDistributionIntegerSampler'"%(algorithm))

        if precision == "mp":
            if not isinstance(sigma, RealNumber):
                RR = RealField()
                sigma = RR(sigma)

            if not isinstance(c, RealNumber):
                c = sigma.parent()(c)
            sig_on()
            self._gen_mp = dgs_disc_gauss_mp_init((<RealNumber>sigma).value, (<RealNumber>c).value, tau, algorithm)
            sig_off()
            self._gen_dp = NULL
            self.sigma = sigma.parent()(0)
            mpfr_set(self.sigma.value, self._gen_mp.sigma, MPFR_RNDN)
            self.c = c
        elif precision == "dp":
            RR = RealField()
            if not isinstance(sigma, RealNumber):
                sigma = RR(sigma)
            sig_on()
            self._gen_dp = dgs_disc_gauss_dp_init(sigma, c, tau, algorithm)
            sig_off()
            self._gen_mp = NULL
            self.sigma = RR(sigma)
            self.c = RR(c)
        else:
            raise ValueError("Parameter precision '%s' not supported."%precision)

        self.tau = Integer(tau)
        self.algorithm = algorithm_str

    def _flush_cache(self):
        r"""
        Flush the internal cache of random bits.

        EXAMPLES::

            sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler

            sage: f = lambda: sage.misc.randstate.set_random_seed(0)

            sage: f()
            sage: D = DiscreteGaussianDistributionIntegerSampler(30.0)
            sage: [D() for k in range(16)]
            [21, 23, 37, 6, -64, 29, 8, -22, -3, -10, 7, -43, 1, -29, 25, 38]

            sage: f()
            sage: D = DiscreteGaussianDistributionIntegerSampler(30.0)
            sage: l = []
            sage: for i in range(16):
            ....:     f(); l.append(D())
            sage: l
            [21, 21, 21, 21, -21, 21, 21, -21, -21, -21, 21, -21, 21, -21, 21, 21]

            sage: f()
            sage: D = DiscreteGaussianDistributionIntegerSampler(30.0)
            sage: l = []
            sage: for i in range(16):
            ....:     f(); D._flush_cache(); l.append(D())
            sage: l
            [21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21]
        """
        if self._gen_mp:
            dgs_disc_gauss_mp_flush_cache(self._gen_mp)
        if self._gen_dp:
            dgs_disc_gauss_dp_flush_cache(self._gen_dp)

    def __dealloc__(self):
        r"""
        TESTS::

            sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
            sage: D = DiscreteGaussianDistributionIntegerSampler(3.0, algorithm='uniform+online')
            sage: del D
        """
        if self._gen_mp:
            dgs_disc_gauss_mp_clear(self._gen_mp)
        if self._gen_dp:
            dgs_disc_gauss_dp_clear(self._gen_dp)

    def __call__(self):
        r"""
        Return a new sample.

        EXAMPLES::

            sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
            sage: DiscreteGaussianDistributionIntegerSampler(3.0, algorithm='uniform+online')()  # random
            -3
            sage: DiscreteGaussianDistributionIntegerSampler(3.0, algorithm='uniform+table')()  # random
            3

        TESTS::

            sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
            sage: DiscreteGaussianDistributionIntegerSampler(3.0, algorithm='uniform+logtable', precision='dp')() # random output
            13
        """
        cdef randstate rstate
        cdef Integer rop
        if self._gen_mp:
            rstate = current_randstate()
            rop = Integer()
            sig_on()
            self._gen_mp.call(rop.value, self._gen_mp, rstate.gmp_state)
            sig_off()
            return rop
        else:
            sig_on()
            r = self._gen_dp.call(self._gen_dp)
            sig_off()
            return Integer(r)

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
            sage: repr(DiscreteGaussianDistributionIntegerSampler(3.0, 2))
            'Discrete Gaussian sampler over the Integers with sigma = 3.000000 and c = 2.000000'
        """
        return f"Discrete Gaussian sampler over the Integers with sigma = {self.sigma:.6f} and c = {self.c:.6f}"
