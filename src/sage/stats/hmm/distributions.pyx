# sage_setup: distribution = sagemath-modules
# sage.doctest: optional - numpy
r"""
Distributions used in implementing Hidden Markov Models

These distribution classes are designed specifically for HMM's and not
for general use in statistics. For example, they have fixed or
non-fixed status, which only make sense relative to being used in a
hidden Markov model.

AUTHOR:

- William Stein, 2010-03
"""

#############################################################################
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

from cpython.object cimport PyObject_RichCompare

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double sqrt(double)

import math
cdef double sqrt2pi = sqrt(2*math.pi)

from sage.misc.randstate cimport current_randstate, randstate
from sage.stats.time_series cimport TimeSeries


cdef double random_normal(double mean, double std, randstate rstate) noexcept:
    r"""
    Return a floating point number chosen from the normal distribution
    with given mean and standard deviation, using the given randstate.
    The computation uses the box muller algorithm.

    INPUT:

    - ``mean`` -- float; the mean
    - ``std`` -- float; the standard deviation
    - ``rstate`` -- randstate; the random number generator state

    OUTPUT:

    - double
    """
    # Ported from http://users.tkk.fi/~nbeijar/soft/terrain/source_o2/boxmuller.c
    # This the box muller algorithm.
    # Client code can get the current random state from:
    #         cdef randstate rstate = current_randstate()
    cdef double x1, x2, w, y1, y2
    while True:
        x1 = 2*rstate.c_rand_double() - 1
        x2 = 2*rstate.c_rand_double() - 1
        w = x1*x1 + x2*x2
        if w < 1: break
    w = sqrt( (-2*log(w))/w )
    y1 = x1 * w
    return mean + y1*std

# Abstract base class for distributions used for hidden Markov models.

cdef class Distribution:
    r"""
    A distribution.
    """
    def sample(self, n=None):
        r"""
        Return either a single sample (the default) or `n` samples from
        this probability distribution.

        INPUT:

        - ``n`` -- ``None`` or a positive integer

        OUTPUT: a single sample if `n` is 1; otherwise many samples

        EXAMPLES:

        This method must be defined in a derived class::

            sage: import sage.stats.hmm.distributions
            sage: sage.stats.hmm.distributions.Distribution().sample()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def prob(self, x):
        r"""
        The probability density function evaluated at `x`.

        INPUT:

        - ``x`` -- object

        OUTPUT: float

        EXAMPLES:

        This method must be defined in a derived class::

            sage: import sage.stats.hmm.distributions
            sage: sage.stats.hmm.distributions.Distribution().prob(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def plot(self, *args, **kwds):
        r"""
        Return a plot of the probability density function.

        INPUT:

        - ``args`` and ``kwds``, passed to the Sage :func:`plot` function

        OUTPUT: a :class:`Graphics` object

        EXAMPLES::

            sage: P = hmm.GaussianMixtureDistribution([(.2,-10,.5),(.6,1,1),(.2,20,.5)])
            sage: P.plot(-10,30)                                                        # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        from sage.plot.all import plot
        return plot(self.prob, *args, **kwds)

cdef class GaussianMixtureDistribution(Distribution):
    r"""
    A probability distribution defined by taking a weighted linear
    combination of Gaussian distributions.

    EXAMPLES::

        sage: P = hmm.GaussianMixtureDistribution([(.3,1,2),(.7,-1,1)]); P
        0.3*N(1.0,2.0) + 0.7*N(-1.0,1.0)
        sage: P[0]
        (0.3, 1.0, 2.0)
        sage: P.is_fixed()
        False
        sage: P.fix(1)
        sage: P.is_fixed(0)
        False
        sage: P.is_fixed(1)
        True
        sage: P.unfix(1)
        sage: P.is_fixed(1)
        False
    """
    def __init__(self, B, eps=1e-8, bint normalize=True):
        r"""
        INPUT:

        - ``B`` -- list of triples ``(c_i, mean_i, std_i)``, where
          the ``c_i`` and ``std_i`` are positive and the sum of the
          ``c_i`` is `1`

        - ``eps`` -- positive real number; any standard deviation in B
          less than eps is replaced by eps

        - ``normalize`` -- if ``True``, ensure that the ``c_i`` are nonnegative

        EXAMPLES::

            sage: hmm.GaussianMixtureDistribution([(.3,1,2),(.7,-1,1)])
            0.3*N(1.0,2.0) + 0.7*N(-1.0,1.0)
            sage: hmm.GaussianMixtureDistribution([(1,-1,0)], eps=1e-3)
            1.0*N(-1.0,0.001)
        """
        B = [[(c if c >= 0 else 0), mu, (std if std > 0 else eps)]
             for c, mu, std in B]
        if not B:
            raise ValueError("must specify at least one component of the mixture model")
        cdef double s
        if normalize:
            s = sum([a[0] for a in B])
            if s != 1:
                if s == 0:
                    s = 1.0/len(B)
                    for a in B:
                        a[0] = s
                else:
                    for a in B:
                        a[0] /= s
        self.c0 = TimeSeries([c/(sqrt2pi*std) for c, _, std in B])
        self.c1 = TimeSeries([-1.0/(2*std*std) for _, _, std in B])
        self.param = TimeSeries(sum([list(x) for x in B], []))
        self.fixed = IntList(self.c0._length)

    def __getitem__(self, Py_ssize_t i):
        r"""
        Return triple (coefficient, mu, std).

        INPUT:

        - ``i`` -- integer

        OUTPUT: triple of floats

        EXAMPLES::

            sage: P = hmm.GaussianMixtureDistribution([(.2,-10,.5),(.6,1,1),(.2,20,.5)])
            sage: P[0]
            (0.2, -10.0, 0.5)
            sage: P[2]
            (0.2, 20.0, 0.5)
            sage: [-1]
            [-1]
            sage: P[-1]
            (0.2, 20.0, 0.5)
            sage: P[3]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: P[-4]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        if i < 0: i += self.param._length//3
        if i < 0 or i >= self.param._length//3:
            raise IndexError("index out of range")
        return self.param._values[3*i], self.param._values[3*i+1], self.param._values[3*i+2]

    def __reduce__(self):
        r"""
        Used in pickling.

        EXAMPLES::

            sage: G = hmm.GaussianMixtureDistribution([(.1,1,2), (.9,0,1)])
            sage: loads(dumps(G)) == G
            True
        """
        return unpickle_gaussian_mixture_distribution_v1, (
            self.c0, self.c1, self.param, self.fixed)

    def __richcmp__(self, other, op):
        r"""
        EXAMPLES::

            sage: G = hmm.GaussianMixtureDistribution([(.1,1,2), (.9,0,1)])
            sage: H = hmm.GaussianMixtureDistribution([(.3,1,2), (.7,1,5)])
            sage: G < H
            True
            sage: H > G
            True
            sage: G == H
            False
            sage: G == G
            True
        """
        if not isinstance(other, GaussianMixtureDistribution):
            return NotImplemented
        return PyObject_RichCompare(self.__reduce__()[1],
                                    other.__reduce__()[1], op)

    def __len__(self):
        r"""
        Return the number of components of this GaussianMixtureDistribution.

        EXAMPLES::

            sage: len(hmm.GaussianMixtureDistribution([(.2,-10,.5),(.6,1,1),(.2,20,.5)]))
            3
        """
        return self.c0._length

    cpdef is_fixed(self, i=None):
        r"""
        Return whether or not this :class:`GaussianMixtureDistribution` is
        fixed when using Baum-Welch to update the corresponding HMM.

        INPUT:

        - ``i`` -- ``None`` (default) or integer; if given, only return
          whether the `i`-th component is fixed

        EXAMPLES::

            sage: P = hmm.GaussianMixtureDistribution([(.2,-10,.5),(.6,1,1),(.2,20,.5)])
            sage: P.is_fixed()
            False
            sage: P.is_fixed(0)
            False
            sage: P.fix(0); P.is_fixed()
            False
            sage: P.is_fixed(0)
            True
            sage: P.fix(); P.is_fixed()
            True
        """
        if i is None:
            return bool(self.fixed.prod())
        else:
            return bool(self.fixed[i])

    def fix(self, i=None):
        r"""
        Set that this :class:`GaussianMixtureDistribution` (or its `i`-th
        component) is fixed when using Baum-Welch to update
        the corresponding HMM.

        INPUT:

        - ``i`` -- ``None`` (default) or integer; if given, only fix the
          `i`-th component

        EXAMPLES::

            sage: P = hmm.GaussianMixtureDistribution([(.2,-10,.5),(.6,1,1),(.2,20,.5)])
            sage: P.fix(1); P.is_fixed()
            False
            sage: P.is_fixed(1)
            True
            sage: P.fix(); P.is_fixed()
            True
        """
        cdef int j
        if i is None:
            for j in range(self.c0._length):
                self.fixed[j] = 1
        else:
            self.fixed[i] = 1

    def unfix(self, i=None):
        r"""
        Set that this :class:`GaussianMixtureDistribution` (or its `i`-th
        component) is not fixed when using Baum-Welch to update the
        corresponding HMM.

        INPUT:

        - ``i`` -- ``None`` (default) or integer; if given, only fix the
          `i`-th component

        EXAMPLES::

            sage: P = hmm.GaussianMixtureDistribution([(.2,-10,.5),(.6,1,1),(.2,20,.5)])
            sage: P.fix(1); P.is_fixed(1)
            True
            sage: P.unfix(1); P.is_fixed(1)
            False
            sage: P.fix(); P.is_fixed()
            True
            sage: P.unfix(); P.is_fixed()
            False
        """
        cdef int j
        if i is None:
            for j in range(self.c0._length):
                self.fixed[j] = 0
        else:
            self.fixed[i] = 0

    def __repr__(self):
        r"""
        Return string representation of this mixed Gaussian distribution.

        EXAMPLES::

            sage: hmm.GaussianMixtureDistribution([(.2,-10,.5),(.6,1,1),(.2,20,.5)]).__repr__()
            '0.2*N(-10.0,0.5) + 0.6*N(1.0,1.0) + 0.2*N(20.0,0.5)'
        """
        return ' + '.join("%s*N(%s,%s)" % x for x in self)

    def sample(self, n=None):
        r"""
        Return a single sample from this distribution (by default), or
        if `n>1`, return a :class:`TimeSeries` of samples.

        INPUT:

        - ``n`` -- integer or ``None`` (default: ``None``)

        OUTPUT:

        - float if ``n`` is ``None`` (default); otherwise a :class:`TimeSeries`

        EXAMPLES::

            sage: P = hmm.GaussianMixtureDistribution([(.2,-10,.5),(.6,1,1),(.2,20,.5)])
            sage: type(P.sample())
            <class 'float'>
            sage: l = P.sample(1)
            sage: len(l)
            1
            sage: type(l)
            <class 'sage.stats.time_series.TimeSeries'>
            sage: l = P.sample(5)
            sage: len(l)
            5
            sage: type(l)
            <class 'sage.stats.time_series.TimeSeries'>
            sage: l = P.sample(0)
            sage: len(l)
            0
            sage: type(l)
            <class 'sage.stats.time_series.TimeSeries'>
            sage: P.sample(-3)
            Traceback (most recent call last):
            ...
            ValueError: n must be nonnegative
        """
        cdef randstate rstate = current_randstate()
        cdef Py_ssize_t i
        cdef TimeSeries T
        if n is None:
            return self._sample(rstate)
        else:
            _n = n
            if _n < 0:
                raise ValueError("n must be nonnegative")
            T = TimeSeries(_n)
            for i in range(_n):
                T._values[i] = self._sample(rstate)
            return T

    cdef double _sample(self, randstate rstate) noexcept:
        r"""
        Used internally to compute a sample from this distribution quickly.

        INPUT:

        - ``rstate`` -- a randstate object

        OUTPUT:

        - double
        """
        cdef double accum, r
        cdef int n
        accum = 0
        r = rstate.c_rand_double()

        # See the remark in hmm.pyx about using GSL to remove this
        # silly way of sampling from a discrete distribution.
        for n in range(self.c0._length):
            accum += self.param._values[3*n]
            if r <= accum:
                return random_normal(self.param._values[3*n+1], self.param._values[3*n+2], rstate)
        raise RuntimeError("invalid probability distribution")

    cpdef double prob(self, double x) noexcept:
        r"""
        Return the probability of `x`.

        Since this is a continuous distribution, this is defined to be
        the limit of the p's such that the probability of [x,x+h] is p*h.

        INPUT:

        - ``x`` -- float

        OUTPUT: float

        EXAMPLES::

            sage: P = hmm.GaussianMixtureDistribution([(.2,-10,.5),(.6,1,1),(.2,20,.5)])
            sage: P.prob(.5)
            0.21123919605857971
            sage: P.prob(-100)
            0.0
            sage: P.prob(20)
            0.1595769121605731
        """
        # The tricky-looking code below is a fast version of this:
        #       return sum([c/(sqrt(2*math.pi)*std) * \
        #                  exp(-(x-mean)*(x-mean)/(2*std*std)) for
        #                  c, mean, std in self.B])
        cdef double s=0, mu
        cdef int n
        for n in range(self.c0._length):
            mu = self.param._values[3*n+1]
            s += self.c0._values[n]*exp((x-mu)*(x-mu)*self.c1._values[n])
        return s

    cpdef double prob_m(self, double x, int m) noexcept:
        r"""
        Return the probability of `x` using just the `m`-th summand.

        INPUT:

        - ``x`` -- float
        - ``m`` -- integer

        OUTPUT: float

        EXAMPLES::

            sage: P = hmm.GaussianMixtureDistribution([(.2,-10,.5),(.6,1,1),(.2,20,.5)])
            sage: P.prob_m(.5, 0)
            2.7608117680508...e-97
            sage: P.prob_m(.5, 1)
            0.21123919605857971
            sage: P.prob_m(.5, 2)
            0.0
        """
        cdef double s, mu
        if m < 0 or m >= self.param._length//3:
            raise IndexError("index out of range")
        mu = self.param._values[3*m+1]
        return self.c0._values[m]*exp((x-mu)*(x-mu)*self.c1._values[m])


def unpickle_gaussian_mixture_distribution_v1(TimeSeries c0, TimeSeries c1,
                                              TimeSeries param, IntList fixed):
    r"""
    Used in unpickling :class:`GaussianMixtureDistribution` objects.

    EXAMPLES::

        sage: P = hmm.GaussianMixtureDistribution([(.2,-10,.5),(.6,1,1),(.2,20,.5)])
        sage: loads(dumps(P)) == P          # indirect doctest
        True
    """
    cdef GaussianMixtureDistribution G = GaussianMixtureDistribution.__new__(GaussianMixtureDistribution)
    G.c0 = c0
    G.c1 = c1
    G.param = param
    G.fixed = fixed
    return G
