# sage.doctest: needs numpy
r"""
Continuous Emission Hidden Markov Models

AUTHOR:

- William Stein, 2010-03
"""

#*****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.object cimport PyObject_RichCompare
from libc.math cimport log, sqrt, exp, isnormal, isfinite, M_PI
cdef double sqrt2pi = sqrt(2*M_PI)
from cysignals.signals cimport sig_on, sig_off

from sage.misc.flatten  import flatten
from sage.structure.element import Matrix

from sage.stats.time_series cimport TimeSeries
from sage.stats.intlist cimport IntList

from sage.stats.hmm.hmm cimport HiddenMarkovModel
from sage.stats.hmm.util cimport HMM_Util
from sage.stats.hmm.distributions cimport GaussianMixtureDistribution

cdef HMM_Util util = HMM_Util()

from sage.misc.randstate cimport current_randstate, randstate


# TODO: DELETE THIS FUNCTION WHEN MOVE Gaussian stuff to distributions.pyx!!! (next version)
cdef double random_normal(double mean, double std, randstate rstate) noexcept:
    r"""
    Return a number chosen randomly with given mean and standard deviation.

    INPUT:

    - ``mean`` -- double
    - ``std`` -- double, standard deviation
    - ``rstate`` -- a randstate object

    OUTPUT:

    a double
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

cdef class GaussianHiddenMarkovModel(HiddenMarkovModel):
    r"""
    Gaussian emissions Hidden Markov Model.

    INPUT:

    - ``A`` -- matrix; the `N \times N` transition matrix
    - ``B`` -- list of pairs ``(mu, sigma)`` that define the distributions
    - ``pi`` -- initial state probabilities
    - ``normalize`` -- boolean (default: ``True``)

    EXAMPLES:

    We illustrate the primary functions with an example 2-state Gaussian HMM::

        sage: m = hmm.GaussianHiddenMarkovModel([[.1,.9],[.5,.5]],
        ....:                                   [(1,1), (-1,1)],
        ....:                                   [.5,.5]); m
        Gaussian Hidden Markov Model with 2 States
        Transition matrix:
        [0.1 0.9]
        [0.5 0.5]
        Emission parameters:
        [(1.0, 1.0), (-1.0, 1.0)]
        Initial probabilities: [0.5000, 0.5000]

    We query the defining transition matrix, emission parameters, and
    initial state probabilities::

        sage: m.transition_matrix()
        [0.1 0.9]
        [0.5 0.5]
        sage: m.emission_parameters()
        [(1.0, 1.0), (-1.0, 1.0)]
        sage: m.initial_probabilities()
        [0.5000, 0.5000]

    We obtain a sample sequence with 10 entries in it, and compute the
    logarithm of the probability of obtaining this sequence, given the
    model::

        sage: obs = m.sample(5); obs  # random
        [-1.6835, 0.0635, -2.1688, 0.3043, -0.3188]
        sage: log_likelihood = m.log_likelihood(obs)
        sage: counter = 0
        sage: n = 0
        sage: def add_samples(i):
        ....:     global counter, n
        ....:     for _ in range(i):
        ....:         n += 1
        ....:         obs2 = m.sample(5)
        ....:         if all(abs(obs2[i] - obs[i]) < 0.25 for i in range(5)):
        ....:             counter += 1

        sage: add_samples(10000)
        sage: while abs(log_likelihood - log(counter*1.0/n/0.5^5)) < 0.1:
        ....:     add_samples(10000)

    We compute the Viterbi path, and probability that the given path
    of states produced obs::

        sage: m.viterbi(obs)  # random
        ([1, 0, 1, 0, 1], -8.714092684611794)

    We use the Baum-Welch iterative algorithm to find another model
    for which our observation sequence is more likely::

        sage: try:
        ....:     p, s = m.baum_welch(obs)
        ....:     assert p > log_likelihood
        ....:     assert (1 <= s <= 500)
        ....: except RuntimeError:
        ....:     pass

    Notice that running Baum-Welch changed our model::

        sage: m  # random
        Gaussian Hidden Markov Model with 2 States
        Transition matrix:
        [   0.4154981366185841     0.584501863381416]
        [   0.9999993174253741 6.825746258991804e-07]
        Emission parameters:
        [(0.4178882427119503, 0.5173109664360919),
         (-1.5025208631331122, 0.5085512836055119)]
        Initial probabilities: [0.0000, 1.0000]
    """
    cdef TimeSeries B, prob
    cdef int n_out

    def __init__(self, A, B, pi, bint normalize=True):
        r"""
        Create a Gaussian emissions HMM with transition probability
        matrix `A`, normal emissions given by `B`, and initial state
        probability distribution ``pi``.

        INPUT:

        - ``A`` -- list of lists or a square `N \times N` matrix, whose
          `(i,j)` entry gives the probability of transitioning from
          state `i` to state `j`.

        - ``B`` -- list of `N` pairs ``(mu, std)``, where if ``B[i]=(mu,std)``,
          then the probability distribution associated with state `i`
          normal with mean ``mu`` and standard deviation ``std``.

        - ``pi`` -- the probabilities of starting in each initial
          state, i.e., ``pi[i]`` is the probability of starting in
          state `i`.

        - ``normalize`` -- boolean (default: ``True``); if given, input is
          normalized to define valid probability distributions,
          e.g., the entries of `A` are made nonnegative and the rows
          sum to 1.

        EXAMPLES::

            sage: hmm.GaussianHiddenMarkovModel([[.1,.9],[.5,.5]], [(1,1), (-1,1)], [.5,.5])
            Gaussian Hidden Markov Model with 2 States
            Transition matrix:
            [0.1 0.9]
            [0.5 0.5]
            Emission parameters:
            [(1.0, 1.0), (-1.0, 1.0)]
            Initial probabilities: [0.5000, 0.5000]

        We input a model in which both `A` and ``pi`` have to be
        renormalized to define valid probability distributions::

            sage: hmm.GaussianHiddenMarkovModel([[-1,.7],[.3,.4]], [(1,1), (-1,1)], [-1,.3])  # rel tol 3e-14
            Gaussian Hidden Markov Model with 2 States
            Transition matrix:
            [                0.0                 1.0]
            [0.42857142857142855  0.5714285714285714]
            Emission parameters:
            [(1.0, 1.0), (-1.0, 1.0)]
            Initial probabilities: [0.0000, 1.0000]

        Bad things can happen::

            sage: hmm.GaussianHiddenMarkovModel([[-1,.7],[.3,.4]], [(1,1), (-1,1)], [-1,.3],
            ....:                               normalize=False)
            Gaussian Hidden Markov Model with 2 States
            Transition matrix:
            [-1.0  0.7]
            [ 0.3  0.4]
            ...
        """
        self.pi = util.initial_probs_to_TimeSeries(pi, normalize)
        self.N = len(self.pi)
        self.A = util.state_matrix_to_TimeSeries(A, self.N, normalize)

        # B should be a matrix of N rows, with column 0 the mean and 1
        # the standard deviation.
        if isinstance(B, Matrix):
            B = B.list()
        else:
            B = flatten(B)
        self.B = TimeSeries(B)
        self.probability_init()

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` and ``other``, which must both be GaussianHiddenMarkovModel's.

        EXAMPLES::

            sage: m = hmm.GaussianHiddenMarkovModel([[1]], [(0,1)], [1])
            sage: n = hmm.GaussianHiddenMarkovModel([[1]], [(1,1)], [1])
            sage: m < n
            True
            sage: m == m
            True
            sage: n > m
            True
            sage: n < m
            False
        """
        if not isinstance(other, GaussianHiddenMarkovModel):
            return NotImplemented
        return PyObject_RichCompare(self.__reduce__()[1],
                                    other.__reduce__()[1], op)

    def __getitem__(self, Py_ssize_t i):
        r"""
        Return the mean and standard distribution for the `i`-th state.

        INPUT:

        - ``i`` -- integer

        OUTPUT: 2 floats

        EXAMPLES::

            sage: m = hmm.GaussianHiddenMarkovModel([[.1,.9],[.5,.5]], [(1,.5), (-2,.3)], [.5,.5])
            sage: m[0]
            (1.0, 0.5)
            sage: m[1]
            (-2.0, 0.3)
            sage: m[-1]
            (-2.0, 0.3)
            sage: m[3]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: m[-3]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        if i < 0:
            i += self.N
        if i < 0 or i >= self.N:
            raise IndexError('index out of range')

        # TODO: change to be a normal distribution class (next version)
        return self.B[2*i], self.B[2*i+1]

    def __reduce__(self):
        r"""
        Used in pickling.

        EXAMPLES::

            sage: G = hmm.GaussianHiddenMarkovModel([[1]], [(0,1)], [1])
            sage: loads(dumps(G)) == G
            True
        """
        return unpickle_gaussian_hmm_v1, \
               (self.A, self.B, self.pi, self.prob, self.n_out)

    def emission_parameters(self):
        r"""
        Return the parameters that define the normal distributions
        associated to all of the states.

        OUTPUT:

        a list ``B`` of pairs ``B[i] = (mu, std)``, such that the
        distribution associated to state `i` is normal with mean
        ``mu`` and standard deviation ``std``.

        EXAMPLES::

            sage: M = hmm.GaussianHiddenMarkovModel([[.1,.9],[.5,.5]],
            ....:                                   [(1,.5), (-1,3)],
            ....:                                   [.1,.9])
            sage: M.emission_parameters()
            [(1.0, 0.5), (-1.0, 3.0)]
        """
        cdef Py_ssize_t i
        from sage.rings.real_double import RDF
        return [(RDF(self.B[2*i]),RDF(self.B[2*i+1])) for i in range(self.N)]

    def __repr__(self):
        r"""
        Return string representation.

        EXAMPLES::

            sage: hmm.GaussianHiddenMarkovModel([[.1,.9],[.5,.5]], [(1,.5), (-1,3)], [.1,.9]).__repr__()
            'Gaussian Hidden Markov Model with 2 States\nTransition matrix:\n[0.1 0.9]\n[0.5 0.5]\nEmission parameters:\n[(1.0, 0.5), (-1.0, 3.0)]\nInitial probabilities: [0.1000, 0.9000]'
        """
        s = "Gaussian Hidden Markov Model with %s States"%self.N
        s += '\nTransition matrix:\n%s'%self.transition_matrix()
        s += '\nEmission parameters:\n%s'%self.emission_parameters()
        s += '\nInitial probabilities: %s'%self.initial_probabilities()
        return s

    def generate_sequence(self, Py_ssize_t length, starting_state=None):
        r"""
        Return a sample of the given length from this HMM.

        INPUT:

        - ``length`` -- positive integer
        - ``starting_state`` -- integer (or ``None``); if specified then
          generate a sequence using this model starting with the given state
          instead of the initial probabilities to determine the
          starting state.

        OUTPUT:

        - an :class:`IntList` or list of emission symbols
        - :class:`TimeSeries` of emissions

        EXAMPLES::

            sage: m = hmm.GaussianHiddenMarkovModel([[.1,.9],[.5,.5]],
            ....:                                   [(1,.5), (-1,3)],
            ....:                                   [.1,.9])
            sage: m.generate_sequence(5)  # random
            ([-3.0505, 0.5317, -4.5065, 0.6521, 1.0435], [1, 0, 1, 0, 1])
            sage: m.generate_sequence(0)
            ([], [])
            sage: m.generate_sequence(-1)
            Traceback (most recent call last):
            ...
            ValueError: length must be nonnegative

        Verify numerically that the starting state is 0 with probability about 0.1::

            sage: counter = 0
            sage: n = 0
            sage: def add_samples(i):
            ....:     global counter, n
            ....:     for i in range(i):
            ....:         n += 1
            ....:         if m.generate_sequence(1)[1][0] == 0:
            ....:             counter += 1

            sage: add_samples(10^5)
            sage: while abs(counter*1.0 / n - 0.1) > 0.01: add_samples(10^5)

        Example in which the starting state is 0 (see :issue:`11452`)::

            sage: set_random_seed(23);  m.generate_sequence(2)
            ([0.6501, -2.0151], [0, 1])

        Force a starting state of 1 even though as we saw above it would be 0::

            sage: set_random_seed(23);  m.generate_sequence(2, starting_state=1)
            ([-3.1491, -1.0244], [1, 1])
        """
        if length < 0:
            raise ValueError("length must be nonnegative")

        # Create Integer lists for states and observations
        cdef IntList states = IntList(length)
        cdef TimeSeries obs = TimeSeries(length)
        if length == 0:
            return states, obs

        # Setup variables, including random state.
        cdef Py_ssize_t i, j
        cdef randstate rstate = current_randstate()
        cdef int q = 0
        cdef double r, accum

        # Choose the starting state.
        # See the remark in hmm.pyx about how this should get
        # replaced by some general fast discrete distribution code.
        if starting_state is None:
            r = rstate.c_rand_double()
            accum = 0
            for i in range(self.N):
                if r < self.pi._values[i] + accum:
                    q = i
                    break
                else:
                    accum += self.pi._values[i]
        else:
            q = starting_state
            if q < 0 or q>= self.N:
                raise ValueError("starting state must be between 0 and %s"%(self.N-1))

        states._values[0] = q
        obs._values[0] = self.random_sample(q, rstate)

        cdef double* row
        cdef int O
        sig_on()
        for i in range(1, length):
            accum = 0
            row = self.A._values + q*self.N
            r = rstate.c_rand_double()
            for j in range(self.N):
                if r < row[j] + accum:
                    q = j
                    break
                else:
                    accum += row[j]
            states._values[i] = q
            obs._values[i] = self.random_sample(q, rstate)
        sig_off()

        return obs, states

    cdef probability_init(self):
        r"""
        Used internally to compute caching information that makes
        certain computations in the Baum-Welch algorithm faster.  This
        function has no input or output.
        """
        self.prob = TimeSeries(2*self.N)
        cdef int i
        for i in range(self.N):
            self.prob[2*i] = 1.0/(sqrt2pi*self.B[2*i+1])
            self.prob[2*i+1] = -1.0/(2*self.B[2*i+1]*self.B[2*i+1])

    cdef double random_sample(self, int state, randstate rstate) noexcept:
        r"""
        Return a random sample from the normal distribution associated
        to the given state.

        This is only used internally, and no bounds or other error
        checking is done, so calling this improperly can lead to seg
        faults.

        INPUT:

        - ``state`` -- integer
        - ``rstate`` -- randstate instance

        OUTPUT:

        double
        """
        return random_normal(self.B._values[state*2], self.B._values[state*2+1], rstate)

    cdef double probability_of(self, int state, double observation) noexcept:
        r"""
        Return a useful continuous analogue of "the probability b_j(o)"
        of seeing the given observation given that we're in the given
        state j (=state).

        The distribution is a normal distribution, and we're asking
        about the probability of a particular point being observed;
        the probability of a particular point is 0, which is not
        useful.  Thus we instead consider the limit p = prob([o,o+d])/d
        as d goes to 0.  There is a simple closed form formula for p,
        derived in the source code.  Note that p can be bigger than 1;
        for example, if we set observation=mean in the closed formula
        we get p=1/(sqrt(2*pi)*std), so p>1 when std<1/sqrt(2*pi).

        INPUT:

        - ``state`` -- integer
        - ``observation`` -- double

        OUTPUT:

        double
        """
        # The code below is an optimized version of the following code:
        #     cdef double mean = self.B._values[2*state], \
        #                 std = self.B._values[2*state+1]
        #     return 1/(sqrt2pi*std) * \
        #             exp(-(observation-mean)*(observation-mean)/(2*std*std))
        #
        # Here is how to use Sage to verify that the above formula computes
        # the limit claimed above:
        #
        # var('x,d,obs,mean,std')
        # n = 1/sqrt(2*pi*std^2) * exp(-(x-mean)^2/(2*std^2))
        # assume(std>0); assume(d>0)
        # m = n.integrate(x,obs,obs+d)/d
        # p = SR(m.limit(d=0).simplify_full())
        # q = 1/(sqrt(2*pi)*std) * exp(-(obs-mean)*(obs-mean)/(2*std*std))
        # bool(p==q)  # outputs True

        cdef double x = observation - self.B._values[2*state] # observation - mean
        return self.prob._values[2*state] * exp(x*x*self.prob._values[2*state+1])

    def log_likelihood(self, obs):
        r"""
        Return the logarithm of a continuous analogue of the
        probability that this model produced the given observation
        sequence.

        Note that the "continuous analogue of the probability" above can
        be bigger than 1, hence the logarithm can be positive.

        INPUT:

        - ``obs`` -- sequence of observations

        OUTPUT: float

        EXAMPLES::

            sage: m = hmm.GaussianHiddenMarkovModel([[.1,.9],[.5,.5]],
            ....:                                   [(1,.5), (-1,3)],
            ....:                                   [.1,.9])
            sage: m.log_likelihood([1,1,1])
            -4.297880766072486
            sage: s = m.sample(20)
            sage: -80 < m.log_likelihood(s) < -20
            True
        """
        if len(obs) == 0:
            return 1.0
        if not isinstance(obs, TimeSeries):
            obs = TimeSeries(obs)
        return self._forward_scale(obs)

    def _forward_scale(self, TimeSeries obs):
        r"""
        Memory-efficient implementation of the forward algorithm (with scaling).

        INPUT:

        - ``obs`` -- integer list of observation states

        OUTPUT:

        float -- the log of the probability that the model
        produced this sequence

        EXAMPLES::

            sage: m = hmm.GaussianHiddenMarkovModel([[.1,.9],[.5,.5]], [(1,.5), (-1,3)], [.1,.9])
            sage: m._forward_scale(stats.TimeSeries([1,-1,-1,1]))
            -7.641988207069133
        """
        cdef Py_ssize_t i, j, t, T = len(obs)

        # The running sum of the log probabilities
        cdef double log_probability = 0, sum, a

        cdef TimeSeries alpha  = TimeSeries(self.N), \
                        alpha2 = TimeSeries(self.N)

        # Initialization
        sum = 0
        for i in range(self.N):
            a = self.pi[i] * self.probability_of(i, obs._values[0])
            alpha[i] = a
            sum += a

        log_probability = log(sum)
        alpha.rescale(1/sum)

        # Induction
        cdef double s
        for t in range(1, T):
            sum = 0
            for j in range(self.N):
                s = 0
                for i in range(self.N):
                    s += alpha._values[i] * self.A._values[i*self.N + j]
                a = s * self.probability_of(j, obs._values[t])
                alpha2._values[j] = a
                sum += a

            log_probability += log(sum)
            for j in range(self.N):
                alpha._values[j] = alpha2._values[j] / sum

        # Termination
        return log_probability

    def viterbi(self, obs):
        r"""
        Determine "the" hidden sequence of states that is most likely
        to produce the given sequence ``obs`` of observations, along with
        the probability that this hidden sequence actually produced
        the observation.

        INPUT:

        - ``obs`` -- sequence of emitted ints or symbols

        OUTPUT:

        - ``list`` -- "the" most probable sequence of hidden states, i.e.,
          the Viterbi path

        - ``float`` -- log of probability that the observed sequence
          was produced by the Viterbi sequence of states

        EXAMPLES:

        We find the optimal state sequence for a given model::

            sage: m = hmm.GaussianHiddenMarkovModel([[0.5,0.5],[0.5,0.5]],
            ....:                                   [(0,1),(10,1)],
            ....:                                   [0.5,0.5])
            sage: m.viterbi([0,1,10,10,1])
            ([0, 0, 1, 1, 0], -9.0604285688230...)

        Another example in which the most likely states change based
        on the last observation::

            sage: m = hmm.GaussianHiddenMarkovModel([[.1,.9],[.5,.5]],
            ....:                                   [(1,.5), (-1,3)],
            ....:                                   [.1,.9])
            sage: m.viterbi([-2,-1,.1,0.1])
            ([1, 1, 0, 1], -9.61823698847639...)
            sage: m.viterbi([-2,-1,.1,0.3])
            ([1, 1, 1, 0], -9.566023653378513)
        """
        cdef TimeSeries _obs
        if not isinstance(obs, TimeSeries):
            _obs = TimeSeries(obs)
        else:
            _obs = obs

        # The algorithm is the same as _viterbi above, except
        # we take the logarithms of everything first, and add
        # instead of multiply.
        cdef Py_ssize_t t, T = _obs._length
        cdef IntList state_sequence = IntList(T)
        if T == 0:
            return state_sequence, 0.0

        cdef int i, j, N = self.N

        # delta[i] is the maximum of the probabilities over all
        # paths ending in state i.
        cdef TimeSeries delta = TimeSeries(N), delta_prev = TimeSeries(N)

        # We view psi as an N x T matrix of ints. The quantity
        #          psi[N*t + j]
        # is a most probable hidden state at time t-1, given
        # that j is a most probable state at time j.
        cdef IntList psi = IntList(N * T)  # initialized to 0 by default

        # Log Preprocessing
        cdef TimeSeries A = self.A.log()
        cdef TimeSeries pi = self.pi.log()

        # Initialization:
        for i in range(N):
            delta._values[i] = pi._values[i] + log(self.probability_of(i, _obs._values[0]))

        # Recursion:
        cdef double mx, tmp, minus_inf = float('-inf')
        cdef int index

        for t in range(1, T):
            delta_prev, delta = delta, delta_prev
            for j in range(N):
                # Compute delta_t[j] = max_i(delta_{t-1}(i) a_{i,j}) * b_j(_obs[t])
                mx = minus_inf
                index = -1
                for i in range(N):
                    tmp = delta_prev._values[i] + A._values[i*N+j]
                    if tmp > mx:
                        mx = tmp
                        index = i
                delta._values[j] = mx + log(self.probability_of(j, _obs._values[t]))
                psi._values[N*t + j] = index

        # Termination:
        mx, index = delta.max(index=True)

        # Path (state sequence) backtracking:
        state_sequence._values[T-1] = index
        t = T-2
        while t >= 0:
            state_sequence._values[t] = psi._values[N*(t+1) + state_sequence._values[t+1]]
            t -= 1

        return state_sequence, mx

    cdef TimeSeries _backward_scale_all(self, TimeSeries obs, TimeSeries scale):
        r"""
        This function returns the matrix beta_t(i), and is used
        internally as part of the Baum-Welch algorithm.

        The quantity beta_t(i) is the probability of observing the
        sequence obs[t+1:] assuming that the model is in state i at
        time t.

        INPUT:

        - ``obs`` -- TimeSeries
        - ``scale`` -- TimeSeries

        OUTPUT:

        - TimeSeries beta such that beta_t(i) = beta[t*N + i]
        - scale is also changed by this function
        """
        cdef Py_ssize_t t, T = obs._length
        cdef int N = self.N, i, j
        cdef double s
        cdef TimeSeries beta = TimeSeries(N*T, initialize=False)

        # 1. Initialization
        for i in range(N):
            beta._values[(T-1)*N + i] = 1 / scale._values[T-1]

        # 2. Induction
        t = T-2
        while t >= 0:
            for i in range(N):
                s = 0
                for j in range(N):
                    s += self.A._values[i*N+j] * \
                         self.probability_of(j, obs._values[t+1]) * beta._values[(t+1)*N+j]
                beta._values[t*N + i] = s/scale._values[t]
            t -= 1
        return beta

    cdef _forward_scale_all(self, TimeSeries obs):
        r"""
        Return scaled values alpha_t(i), the sequence of scalings, and
        the log probability.

        The quantity alpha_t(i) is the probability of observing the
        sequence obs[:t+1] assuming that the model is in state i at
        time t.

        INPUT:

        - ``obs`` -- TimeSeries

        OUTPUT:

        - TimeSeries alpha with alpha_t(i) = alpha[t*N + i]
        - TimeSeries scale with scale[t] the scaling at step t
        - ``float`` -- log_probability of the observation sequence
          being produced by the model
        """
        cdef Py_ssize_t i, j, t, T = len(obs)
        cdef int N = self.N

        # The running some of the log probabilities
        cdef double log_probability = 0, sum, a

        cdef TimeSeries alpha = TimeSeries(N*T, initialize=False)
        cdef TimeSeries scale = TimeSeries(T, initialize=False)

        # Initialization
        sum = 0
        for i in range(self.N):
            a = self.pi._values[i] * self.probability_of(i, obs._values[0])
            alpha._values[0*N + i] = a
            sum += a

        scale._values[0] = sum
        log_probability = log(sum)
        for i in range(self.N):
            alpha._values[0*N + i] /= sum

        # Induction
        # The code below is just an optimized version of:
        #     alpha = (alpha * A).pairwise_product(B[O[t+1]])
        #     alpha = alpha.scale(1/alpha.sum())
        # along with keeping track of the log of the scaling factor.
        cdef double s
        for t in range(1,T):
            sum = 0
            for j in range(N):
                s = 0
                for i in range(N):
                    s += alpha._values[(t-1)*N + i] * self.A._values[i*N + j]
                a = s * self.probability_of(j, obs._values[t])
                alpha._values[t*N + j] = a
                sum += a

            log_probability += log(sum)
            scale._values[t] = sum
            for j in range(self.N):
                alpha._values[t*N + j] /= sum

        # Termination
        return alpha, scale, log_probability

    cdef TimeSeries _baum_welch_xi(self, TimeSeries alpha, TimeSeries beta, TimeSeries obs):
        r"""
        Used internally to compute the scaled quantity xi_t(i,j)
        appearing in the Baum-Welch reestimation algorithm.

        INPUT:

        - ``alpha`` -- TimeSeries as output by the scaled forward algorithm
        - ``beta`` -- TimeSeries as output by the scaled backward algorithm
        - ``obs`` -- TimeSeries of observations

        OUTPUT:

        TimeSeries xi such that xi[t*N*N + i*N + j] = xi_t(i,j).
        """
        cdef int i, j, N = self.N
        cdef double sum
        cdef Py_ssize_t t, T = alpha._length//N
        cdef TimeSeries xi = TimeSeries(T*N*N, initialize=False)
        for t in range(T-1):
            sum = 0.0
            for i in range(N):
                for j in range(N):
                    xi._values[t*N*N+i*N+j] = alpha._values[t*N+i]*beta._values[(t+1)*N+j]*\
                       self.A._values[i*N+j] * self.probability_of(j, obs._values[t+1])
                    sum += xi._values[t*N*N+i*N+j]
            for i in range(N):
                for j in range(N):
                    xi._values[t*N*N+i*N+j] /= sum
        return xi

    def baum_welch(self, obs, int max_iter=500, double log_likelihood_cutoff=1e-4,
                   double min_sd=0.01, bint fix_emissions=False, bint v=False):
        r"""
        Given an observation sequence ``obs``, improve this HMM using the
        Baum-Welch algorithm to increase the probability of observing ``obs``.

        INPUT:

        - ``obs`` -- a time series of emissions

        - ``max_iter`` -- integer (default: 500); maximum number
          of Baum-Welch steps to take

        - ``log_likelihood_cutoff`` -- positive float (default: 1e-4);
          the minimal improvement in likelihood with respect to
          the last iteration required to continue. Relative value
          to log likelihood.

        - ``min_sd`` -- positive float (default: 0.01); when
          reestimating, the standard deviation of emissions is not
          allowed to be less than ``min_sd``.

        - ``fix_emissions`` -- boolean (default: ``False``); if ``True``, do not
          change emissions when updating

        OUTPUT:

        changes the model in place, and returns the log
        likelihood and number of iterations.

        EXAMPLES::

            sage: m = hmm.GaussianHiddenMarkovModel([[.1,.9],[.5,.5]],
            ....:                                   [(1,.5), (-1,3)],
            ....:                                   [.1,.9])
            sage: m.log_likelihood([-2,-1,.1,0.1])
            -8.858282215986275
            sage: m.baum_welch([-2,-1,.1,0.1])
            (4.534646052182..., 7)
            sage: m.log_likelihood([-2,-1,.1,0.1])
            4.534646052182...
            sage: m  # rel tol 3e-14
            Gaussian Hidden Markov Model with 2 States
            Transition matrix:
            [   0.9999999992430161 7.569839394440382e-10]
            [  0.49998462791192644    0.5000153720880736]
            Emission parameters:
            [(0.09999999999999999, 0.01), (-1.4999508147591902, 0.5000710504895474)]
            Initial probabilities: [0.0000, 1.0000]

        We illustrate bounding the standard deviation below.  Note that above we had
        different emission parameters when the ``min_sd`` was the default of 0.01::

            sage: m = hmm.GaussianHiddenMarkovModel([[.1,.9],[.5,.5]],
            ....:                                   [(1,.5), (-1,3)],
            ....:                                   [.1,.9])
            sage: m.baum_welch([-2,-1,.1,0.1], min_sd=1)
            (-4.07939572755..., 32)
            sage: m.emission_parameters()
            [(-0.2663018798..., 1.0), (-1.99850979..., 1.0)]

        We watch the log likelihoods of the model converge, step by step::

            sage: m = hmm.GaussianHiddenMarkovModel([[.1,.9],[.5,.5]],
            ....:                                   [(1,.5), (-1,3)],
            ....:                                   [.1,.9])
            sage: v = m.sample(10)
            sage: l = stats.TimeSeries([m.baum_welch(v, max_iter=1)[0]
            ....:                       for _ in range(len(v))])
            sage: all(l[i] <= l[i+1] + 0.0001 for i in range(9))
            True
            sage: l  # random
            [-20.1167, -17.7611, -16.9814, -16.9364, -16.9314,
             -16.9309, -16.9309, -16.9309, -16.9309, -16.9309]

        We illustrate fixing emissions::

            sage: m = hmm.GaussianHiddenMarkovModel([[.1,.9],[.9,.1]],
            ....:                                   [(1,2),(-1,.5)],
            ....:                                   [.3,.7])
            sage: set_random_seed(0); v = m.sample(100)
            sage: m.baum_welch(v,fix_emissions=True)
            (-164.72944548204..., 23)
            sage: m.emission_parameters()
            [(1.0, 2.0), (-1.0, 0.5)]
            sage: m = hmm.GaussianHiddenMarkovModel([[.1,.9],[.9,.1]],
            ....:                                   [(1,2),(-1,.5)],
            ....:                                   [.3,.7])
            sage: m.baum_welch(v)
            (-162.854370397998..., 49)
            sage: m.emission_parameters()  # rel tol 3e-14
            [(1.2722419172602375, 2.371368751761901),
             (-0.9486174675179113, 0.5762360385123765)]
        """
        if not isinstance(obs, TimeSeries):
            obs = TimeSeries(obs)
        cdef TimeSeries _obs = obs
        cdef TimeSeries alpha, beta, scale, gamma, xi
        cdef double log_probability, log_probability0, log_probability_prev, delta
        cdef int i, j, k, N, n_iterations
        cdef Py_ssize_t t, T
        cdef double denominator_A, numerator_A, denominator_B, numerator_mean, numerator_std

        # Initialization
        alpha, scale, log_probability0 = self._forward_scale_all(_obs)
        if not isfinite(log_probability0):
            return (0.0, 0)
        log_probability = log_probability0
        beta = self._backward_scale_all(_obs, scale)
        gamma = self._baum_welch_gamma(alpha, beta)
        xi = self._baum_welch_xi(alpha, beta, _obs)
        log_probability_prev = log_probability
        N = self.N
        n_iterations = 0
        T = len(_obs)

        # Re-estimation
        while True:
            # Reestimate
            for i in range(N):
                if not isfinite(gamma._values[0*N+i]):
                    # Before raising an error, leave self in a valid state.
                    util.normalize_probability_TimeSeries(self.pi, 0, self.pi._length)
                    raise RuntimeError("impossible to compute gamma during reestimation")
                self.pi._values[i] = gamma._values[0*N+i]

            # Update the probabilities pi to define a valid discrete distribution
            util.normalize_probability_TimeSeries(self.pi, 0, self.pi._length)

            # Reestimate transition matrix and emission probabilities in
            # each state.
            for i in range(N):
                # Compute the updated transition matrix
                denominator_A = 0.0
                for t in range(T-1):
                    denominator_A += gamma._values[t*N+i]
                if not isnormal(denominator_A):
                    raise RuntimeError("unable to re-estimate transition matrix")
                for j in range(N):
                    numerator_A = 0.0
                    for t in range(T-1):
                        numerator_A += xi._values[t*N*N+i*N+j]
                    self.A._values[i*N+j] = numerator_A / denominator_A

                # Rescale the i-th row of the transition matrix to be
                # a valid stochastic matrix:
                util.normalize_probability_TimeSeries(self.A, i*N, (i+1)*N)

                if not fix_emissions:
                    denominator_B = denominator_A + gamma._values[(T-1)*N + i]
                    if not isnormal(denominator_B):
                        raise RuntimeError("unable to re-estimate emission probabilities")

                    numerator_mean = 0.0
                    numerator_std = 0.0
                    for t in range(T):
                        numerator_mean += gamma._values[t*N + i] * _obs._values[t]
                        numerator_std  += gamma._values[t*N + i] * \
                               (_obs._values[t] - self.B._values[2*i])*(_obs._values[t] - self.B._values[2*i])
                    # restimated mean
                    self.B._values[2*i] = numerator_mean / denominator_B
                    # restimated standard deviation
                    self.B._values[2*i+1] = sqrt(numerator_std / denominator_B)
                    if self.B._values[2*i+1] < min_sd:
                        self.B._values[2*i+1] = min_sd
                    self.probability_init()

            n_iterations += 1
            if n_iterations >= max_iter: break

            # Initialization for next iteration
            alpha, scale, log_probability0 = self._forward_scale_all(_obs)

            if not isfinite(log_probability0): break
            log_probability = log_probability0
            beta = self._backward_scale_all(_obs, scale)
            gamma = self._baum_welch_gamma(alpha, beta)
            xi = self._baum_welch_xi(alpha, beta, _obs)

            # Compute the difference between the log probability of
            # two iterations.
            delta = log_probability - log_probability_prev
            log_probability_prev = log_probability

            # If the log probability does not change by more than delta,
            # then terminate
            if delta >= 0 and delta <= log_likelihood_cutoff:
                break

        return log_probability, n_iterations


cdef class GaussianMixtureHiddenMarkovModel(GaussianHiddenMarkovModel):
    r"""
    Gaussian mixture Hidden Markov Model.

    INPUT:

    - ``A`` -- matrix; the `N \times N` transition matrix

    - ``B`` -- list of mixture definitions for each state.  Each
      state may have a varying number of gaussians with selection
      probabilities that sum to 1 and encoded as ``(p, (mu,sigma))``

    - ``pi`` -- initial state probabilities

    - ``normalize`` -- boolean (default: ``True``); if given, input is
      normalized to define valid probability distributions,
      e.g., the entries of `A` are made nonnegative and the rows
      sum to 1, and the probabilities in ``pi`` are normalized.

    EXAMPLES::

        sage: A = [[0.5,0.5],[0.5,0.5]]
        sage: B = [[(0.9,(0.0,1.0)), (0.1,(1,10000))],[(1,(1,1)), (0,(0,0.1))]]
        sage: hmm.GaussianMixtureHiddenMarkovModel(A, B, [1,0])
        Gaussian Mixture Hidden Markov Model with 2 States
        Transition matrix:
        [0.5 0.5]
        [0.5 0.5]
        Emission parameters:
        [0.9*N(0.0,1.0) + 0.1*N(1.0,10000.0), 1.0*N(1.0,1.0) + 0.0*N(0.0,0.1)]
        Initial probabilities: [1.0000, 0.0000]

    TESTS:

    If a standard deviation is 0, it is normalized to be slightly bigger than 0.::

        sage: hmm.GaussianMixtureHiddenMarkovModel([[1]], [[(1,(0,0))]], [1])
        Gaussian Mixture Hidden Markov Model with 1 States
        Transition matrix:
        [1.0]
        Emission parameters:
        [1.0*N(0.0,1e-08)]
        Initial probabilities: [1.0000]

    We test that number of emission distributions must be the same as the number of states::

        sage: hmm.GaussianMixtureHiddenMarkovModel([[1]], [], [1])
        Traceback (most recent call last):
        ...
        ValueError: number of GaussianMixtures must be the same as number of entries of pi

        sage: hmm.GaussianMixtureHiddenMarkovModel([[1]], [[]], [1])
        Traceback (most recent call last):
        ...
        ValueError: must specify at least one component of the mixture model

    We test that the number of initial probabilities must equal the number of states::

        sage: hmm.GaussianMixtureHiddenMarkovModel([[1]], [[]], [1,2])
        Traceback (most recent call last):
        ...
        ValueError: number of entries of transition matrix A must be the square of the number of entries of pi
    """

    cdef object mixture # mixture

    def __init__(self, A, B, pi=None, bint normalize=True):
        r"""
        Initialize a Gaussian mixture hidden Markov model.

        EXAMPLES::

            sage: hmm.GaussianMixtureHiddenMarkovModel([[.9,.1],[.4,.6]], [[(.4,(0,1)), (.6,(1,0.1))],[(1,(0,1))]], [.7,.3])
            Gaussian Mixture Hidden Markov Model with 2 States
            Transition matrix:
            [0.9 0.1]
            [0.4 0.6]
            Emission parameters:
            [0.4*N(0.0,1.0) + 0.6*N(1.0,0.1), 1.0*N(0.0,1.0)]
            Initial probabilities: [0.7000, 0.3000]
        """
        self.pi = util.initial_probs_to_TimeSeries(pi, normalize)
        self.N = len(self.pi)
        self.A = util.state_matrix_to_TimeSeries(A, self.N, normalize)
        if self.N*self.N != len(self.A):
            raise ValueError("number of entries of transition matrix A must be the square of the number of entries of pi")

        self.mixture = [b if isinstance(b, GaussianMixtureDistribution) else
                        GaussianMixtureDistribution([flatten(x) for x in b])
                        for b in B]
        if len(self.mixture) != self.N:
            raise ValueError("number of GaussianMixtures must be the same as number of entries of pi")

    def __repr__(self):
        r"""
        Return string representation.

        EXAMPLES::

            sage: hmm.GaussianMixtureHiddenMarkovModel([[.9,.1],[.4,.6]], [[(.4,(0,1)), (.6,(1,0.1))],[(1,(0,1))]], [.7,.3]).__repr__()
            'Gaussian Mixture Hidden Markov Model with 2 States\nTransition matrix:\n[0.9 0.1]\n[0.4 0.6]\nEmission parameters:\n[0.4*N(0.0,1.0) + 0.6*N(1.0,0.1), 1.0*N(0.0,1.0)]\nInitial probabilities: [0.7000, 0.3000]'
        """
        s = "Gaussian Mixture Hidden Markov Model with %s States"%self.N
        s += '\nTransition matrix:\n%s'%self.transition_matrix()
        s += '\nEmission parameters:\n%s'%self.emission_parameters()
        s += '\nInitial probabilities: %s'%self.initial_probabilities()
        return s

    def __reduce__(self):
        r"""
        Used in pickling.

        EXAMPLES::

            sage: m = hmm.GaussianMixtureHiddenMarkovModel([[1]], [[(.4,(0,1)), (.6,(1,0.1))]], [1])
            sage: loads(dumps(m)) == m
            True
        """
        return unpickle_gaussian_mixture_hmm_v1, \
               (self.A, self.B, self.pi, self.mixture)

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` and ``other``, which must both be
        ``GaussianMixtureHiddenMarkovModel``s.

        EXAMPLES::

            sage: m = hmm.GaussianMixtureHiddenMarkovModel([[1]], [[(.4,(0,1)), (.6,(1,0.1))]], [1])
            sage: n = hmm.GaussianMixtureHiddenMarkovModel([[1]], [[(.5,(0,1)), (.5,(1,0.1))]], [1])
            sage: m < n
            True
            sage: m == m
            True
            sage: n > m
            True
            sage: n < m
            False
        """
        if not isinstance(other, GaussianMixtureHiddenMarkovModel):
            return NotImplemented
        return PyObject_RichCompare(self.__reduce__()[1],
                                    other.__reduce__()[1], op)

    def __getitem__(self, Py_ssize_t i):
        r"""
        Return the Gaussian mixture distribution associated to the
        i-th state.

        INPUT:

        - ``i`` -- integer

        OUTPUT: a Gaussian mixture distribution object

        EXAMPLES::

            sage: m = hmm.GaussianMixtureHiddenMarkovModel([[.9,.1],[.4,.6]], [[(.4,(0,1)), (.6,(1,0.1))],[(1,(0,1))]], [.7,.3])
            sage: m[0]
            0.4*N(0.0,1.0) + 0.6*N(1.0,0.1)
            sage: m[1]
            1.0*N(0.0,1.0)

        Negative indexing works::

            sage: m[-1]
            1.0*N(0.0,1.0)

        Bounds are checked::

            sage: m[2]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: m[-3]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        if i < 0:
            i += self.N
        if i < 0 or i >= self.N:
            raise IndexError('index out of range')
        return self.mixture[i]

    def emission_parameters(self):
        r"""
        Return a list of all the emission distributions.

        OUTPUT: list of Gaussian mixtures

        EXAMPLES::

            sage: m = hmm.GaussianMixtureHiddenMarkovModel([[.9,.1],[.4,.6]],
            ....:                                          [[(.4,(0,1)), (.6,(1,0.1))], [(1,(0,1))]],
            ....:                                          [.7,.3])
            sage: m.emission_parameters()
            [0.4*N(0.0,1.0) + 0.6*N(1.0,0.1), 1.0*N(0.0,1.0)]
        """
        return list(self.mixture)

    cdef double random_sample(self, int state, randstate rstate) noexcept:
        r"""
        Return a random sample from the normal distribution associated
        to the given state.

        This is only used internally, and no bounds or other error
        checking is done, so calling this improperly can lead to seg
        faults.

        INPUT:

        - ``state`` -- integer
        - ``rstate`` -- randstate instance

        OUTPUT:

        double
        """
        cdef GaussianMixtureDistribution G = self.mixture[state]
        return G._sample(rstate)

    cdef double probability_of(self, int state, double observation) noexcept:
        r"""
        Return the probability b_j(o) of see the given observation o
        (=observation) given that we're in the given state j (=state).

        This is a continuous probability, so this really returns a
        number p such that the probability of a value in the interval
        [o,o+d] is p*d.

        INPUT:

        - ``state`` -- integer
        - ``observation`` -- double

        OUTPUT:

        double
        """
        cdef GaussianMixtureDistribution G = self.mixture[state]
        return G.prob(observation)

    cdef TimeSeries _baum_welch_mixed_gamma(self, TimeSeries alpha, TimeSeries beta,
                                            TimeSeries obs, int j):
        r"""
        Let gamma_t(j,m) be the m-component (in the mixture) of the
        probability of being in state j at time t, given the
        observation sequence.  This function outputs a TimeSeries v
        such that v[m*T + t] gives gamma_t(j, m) where T is the number
        of time steps.

        INPUT:

        - ``alpha`` -- TimeSeries
        - ``beta`` -- TimeSeries
        - ``obs`` -- TimeSeries
        - ``j`` -- integer

        OUTPUT:

        TimeSeries
        """
        cdef int i, k, m, N = self.N
        cdef Py_ssize_t t, T = alpha._length//N

        cdef double numer, alpha_minus, P, s, prob
        cdef GaussianMixtureDistribution G = self.mixture[j]
        cdef int M = len(G)
        cdef TimeSeries mixed_gamma = TimeSeries(T*M)

        for t in range(T):
            prob = self.probability_of(j, obs._values[t])
            if prob == 0:
                # If the probability of observing obs[t] in state j is 0, then
                # all of the m-mixture components *have* to automatically be 0,
                # since prob is the sum of those and they are all nonnegative.
                for m in range(M):
                    mixed_gamma._values[m*T + t] = 0
            else:
                # Compute the denominator we used when scaling gamma.
                # The key thing is that this is consistent between
                # gamma and mixed_gamma.
                P = 0
                for k in range(N):
                    P += alpha._values[t*N+k]*beta._values[t*N+k]

                # Divide out the total probability, so we can multiply back in
                # the m-components of the probability.
                alpha_minus = alpha._values[t*N + j] / prob
                for m in range(M):
                    numer = alpha_minus * G.prob_m(obs._values[t], m) * beta._values[t*N + j]
                    mixed_gamma._values[m*T + t] = numer / P

        return mixed_gamma

    def baum_welch(self, obs, int max_iter=1000, double log_likelihood_cutoff=1e-12,
                   double min_sd=0.01, bint fix_emissions=False):
        r"""
        Given an observation sequence ``obs``, improve this HMM using the
        Baum-Welch algorithm to increase the probability of observing ``obs``.

        INPUT:

        - ``obs`` -- a time series of emissions
        - ``max_iter`` -- integer (default: 1000); maximum number
          of Baum-Welch steps to take
        - ``log_likelihood_cutoff`` -- positive float (default: 1e-12);
          the minimal improvement in likelihood with respect to
          the last iteration required to continue. Relative value
          to log likelihood.
        - ``min_sd`` -- positive float (default: 0.01); when
          reestimating, the standard deviation of emissions is not
          allowed to be less than ``min_sd``.
        - ``fix_emissions`` -- boolean (default: ``False``); if ``True``, do not
          change emissions when updating

        OUTPUT:

        changes the model in place, and returns the log
        likelihood and number of iterations.

        EXAMPLES::

            sage: m = hmm.GaussianMixtureHiddenMarkovModel(
            ....:         [[.9,.1],[.4,.6]],
            ....:         [[(.4,(0,1)), (.6,(1,0.1))], [(1,(0,1))]],
            ....:         [.7,.3])
            sage: set_random_seed(0); v = m.sample(10); v
            [0.3576, -0.9365, 0.9449, -0.6957, 1.0217,
             0.9644, 0.9987, -0.5950, -1.0219, 0.6477]
            sage: m.log_likelihood(v)
            -8.31408655939536...
            sage: m.baum_welch(v)
            (2.18905068682..., 15)
            sage: m.log_likelihood(v)
            2.18905068682...
            sage: m  # rel tol 6e-12
            Gaussian Mixture Hidden Markov Model with 2 States
            Transition matrix:
            [   0.8746363339773399   0.12536366602266016]
            [                  1.0 1.451685202290174e-40]
            Emission parameters:
            [0.500161629343*N(-0.812298726239,0.173329026744)
              + 0.499838370657*N(0.982433690378,0.029719932009),
             1.0*N(0.503260056832,0.145881515324)]
            Initial probabilities: [0.0000, 1.0000]

        We illustrate bounding the standard deviation below.  Note that above we had
        different emission parameters when the min_sd was the default of 0.01::

            sage: m = hmm.GaussianMixtureHiddenMarkovModel(
            ....:         [[.9,.1],[.4,.6]],
            ....:         [[(.4,(0,1)), (.6,(1,0.1))], [(1,(0,1))]],
            ....:         [.7,.3])
            sage: m.baum_welch(v, min_sd=1)
            (-12.617885761692..., 1000)
            sage: m.emission_parameters()  # rel tol 6e-12
            [0.503545634447*N(0.200166509595,1.0) + 0.496454365553*N(0.200166509595,1.0),
             1.0*N(0.0543433426535,1.0)]

        We illustrate fixing all emissions::

            sage: m = hmm.GaussianMixtureHiddenMarkovModel(
            ....:         [[.9,.1],[.4,.6]],
            ....:         [[(.4,(0,1)), (.6,(1,0.1))], [(1,(0,1))]],
            ....:         [.7,.3])
            sage: set_random_seed(0); v = m.sample(10)
            sage: m.baum_welch(v, fix_emissions=True)
            (-7.58656858997..., 36)
            sage: m.emission_parameters()
            [0.4*N(0.0,1.0) + 0.6*N(1.0,0.1),
             1.0*N(0.0,1.0)]
        """
        if not isinstance(obs, TimeSeries):
            obs = TimeSeries(obs)
        cdef TimeSeries _obs = obs
        cdef TimeSeries alpha, beta, scale, gamma, mixed_gamma, mixed_gamma_m, xi
        cdef double log_probability, log_probability0, log_probability_prev, delta
        cdef int i, j, k, m, N, n_iterations
        cdef Py_ssize_t t, T
        cdef double denominator_A, numerator_A, denominator_B, numerator_mean, numerator_std, \
             numerator_c, c, mu, std, numer, denom, new_mu, new_std, new_c, s
        cdef GaussianMixtureDistribution G

        # Initialization
        alpha, scale, log_probability0 = self._forward_scale_all(_obs)
        if not isfinite(log_probability0):
            return (0.0, 0)
        log_probability = log_probability0
        beta = self._backward_scale_all(_obs, scale)
        gamma = self._baum_welch_gamma(alpha, beta)
        xi = self._baum_welch_xi(alpha, beta, _obs)
        log_probability_prev = log_probability
        N = self.N
        n_iterations = 0
        T = len(_obs)

        # Re-estimation
        while True:

            # Reestimate frequency of state i in time t=0
            for i in range(N):
                if not isfinite(gamma._values[0*N+i]):
                    # Before raising an error, leave self in a valid state.
                    util.normalize_probability_TimeSeries(self.pi, 0, self.pi._length)
                    raise RuntimeError("impossible to compute gamma during reestimation")
                self.pi._values[i] = gamma._values[0*N+i]

            # Update the probabilities pi to define a valid discrete distribution
            util.normalize_probability_TimeSeries(self.pi, 0, self.pi._length)

            # Reestimate transition matrix and emission probabilities in
            # each state.
            for i in range(N):
                # Reestimate the state transition matrix
                denominator_A = 0.0
                for t in range(T-1):
                    denominator_A += gamma._values[t*N+i]
                if not isnormal(denominator_A):
                    raise RuntimeError("unable to re-estimate pi (1)")
                for j in range(N):
                    numerator_A = 0.0
                    for t in range(T-1):
                        numerator_A += xi._values[t*N*N+i*N+j]
                    self.A._values[i*N+j] = numerator_A / denominator_A

                # Rescale the i-th row of the transition matrix to be
                # a valid stochastic matrix:
                util.normalize_probability_TimeSeries(self.A, i*N, (i+1)*N)

                ########################################################################
                # Re-estimate the emission probabilities
                ########################################################################
                G = self.mixture[i]
                if not fix_emissions and not G.is_fixed():
                    mixed_gamma = self._baum_welch_mixed_gamma(alpha, beta, _obs, i)
                    new_G = []
                    for m in range(len(G)):
                        if G.fixed._values[m]:
                            new_G.append(G[m])
                            continue

                        # Compute re-estimated mu_{j,m}
                        numer = 0
                        denom = 0
                        for t in range(T):
                            numer += mixed_gamma._values[m*T + t] * _obs._values[t]
                            denom += mixed_gamma._values[m*T + t]
                        new_mu = numer / denom

                        # Compute re-estimated standard deviation
                        numer = 0
                        mu = G[m][1]
                        for t in range(T):
                            numer += mixed_gamma._values[m*T + t] * \
                                     (_obs._values[t] - mu)*(_obs._values[t] - mu)

                        new_std = sqrt(numer / denom)
                        if new_std < min_sd:
                            new_std = min_sd

                        # Compute re-estimated weighting coefficient
                        new_c = denom
                        s = 0
                        for t in range(T):
                            s += gamma._values[t*N + i]
                        new_c /= s

                        new_G.append((new_c,new_mu,new_std))

                    self.mixture[i] = GaussianMixtureDistribution(new_G)

            n_iterations += 1
            if n_iterations >= max_iter: break

            ########################################################################
            # Initialization for next iteration
            ########################################################################
            alpha, scale, log_probability0 = self._forward_scale_all(_obs)
            if not isfinite(log_probability0): break
            log_probability = log_probability0
            beta = self._backward_scale_all(_obs, scale)
            gamma = self._baum_welch_gamma(alpha, beta)
            xi = self._baum_welch_xi(alpha, beta, _obs)

            # Compute the difference between the log probability of
            # two iterations.
            delta = log_probability - log_probability_prev
            log_probability_prev = log_probability

            # If the log probability does not improve by more than
            # delta, then terminate
            if delta >= 0 and delta <= log_likelihood_cutoff:
                break

        return log_probability, n_iterations


##################################################
# For Unpickling
##################################################

# We keep the _v0 function for backwards compatible.
def unpickle_gaussian_hmm_v0(A, B, pi, name):
    r"""
    EXAMPLES::

        sage: m = hmm.GaussianHiddenMarkovModel([[1]], [(0,1)], [1])
        sage: sage.stats.hmm.chmm.unpickle_gaussian_hmm_v0(m.transition_matrix(), m.emission_parameters(), m.initial_probabilities(), 'test')
        Gaussian Hidden Markov Model with 1 States
        Transition matrix:
        [1.0]
        Emission parameters:
        [(0.0, 1.0)]
        Initial probabilities: [1.0000]
    """
    return GaussianHiddenMarkovModel(A,B,pi)


def unpickle_gaussian_hmm_v1(A, B, pi, prob, n_out):
    r"""
    EXAMPLES::

        sage: m = hmm.GaussianHiddenMarkovModel([[1]], [(0,1)], [1])
        sage: loads(dumps(m)) == m   # indirect test
        True
    """
    cdef GaussianHiddenMarkovModel m = GaussianHiddenMarkovModel.__new__(GaussianHiddenMarkovModel)
    m.A = A
    m.B = B
    m.pi = pi
    m.prob = prob
    m.n_out = n_out
    return m


def unpickle_gaussian_mixture_hmm_v1(A, B, pi, mixture):
    r"""
    EXAMPLES::

        sage: m = hmm.GaussianMixtureHiddenMarkovModel([[1]], [[(.4,(0,1)), (.6,(1,0.1))]], [1])
        sage: loads(dumps(m)) == m   # indirect test
        True
    """
    cdef GaussianMixtureHiddenMarkovModel m = GaussianMixtureHiddenMarkovModel.__new__(GaussianMixtureHiddenMarkovModel)
    m.A = A
    m.B = B
    m.pi = pi
    m.mixture = mixture
    return m
