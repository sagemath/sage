# sage_setup: distribution = sagemath-modules
# sage.doctest: needs numpy sage.modules
r"""
Hidden Markov Models

This is a complete pure-Cython optimized implementation of Hidden
Markov Models.  It fully supports Discrete, Gaussian, and Mixed
Gaussian emissions.

The best references for the basic HMM algorithms implemented here are:

-  Tapas Kanungo's "Hidden Markov Models"

-  Jackson's HMM tutorial:
    http://personal.ee.surrey.ac.uk/Personal/P.Jackson/tutorial/

LICENSE: Some of the code in this file is based on reading Kanungo's
GPLv2+ implementation of discrete HMM's, hence the present code must
be licensed with a GPLv2+ compatible license.

AUTHOR:

   - William Stein, 2010-03
"""

# ****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libc.math cimport log
from cysignals.signals cimport sig_on, sig_off

from sage.stats.time_series cimport TimeSeries
from sage.structure.element import Matrix
from sage.matrix.constructor import matrix
from sage.misc.randstate cimport current_randstate, randstate
from cpython.object cimport PyObject_RichCompare

from sage.stats.hmm.util cimport HMM_Util

cdef HMM_Util util = HMM_Util()

###########################################

cdef class HiddenMarkovModel:
    r"""
    Abstract base class for all Hidden Markov Models.
    """
    def initial_probabilities(self):
        r"""
        Return the initial probabilities as a :class:`TimeSeries` of
        length `N`, where `N` is the number of states of the Markov model.

        EXAMPLES::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]],
            ....:                                   [[0.1,0.9],[0.5,0.5]],
            ....:                                   [.2,.8])
            sage: pi = m.initial_probabilities(); pi
            [0.2000, 0.8000]
            sage: type(pi)
            <... 'sage.stats.time_series.TimeSeries'>

        The returned time series is a copy, so changing it does not
        change the model::

            sage: pi[0] = .1; pi[1] = .9
            sage: m.initial_probabilities()
            [0.2000, 0.8000]

        Some other models::

            sage: m = hmm.GaussianHiddenMarkovModel([[.1,.9],[.5,.5]],
            ....:                                   [(1,1), (-1,1)],
            ....:                                   [.1,.9])
            sage: m.initial_probabilities()
            [0.1000, 0.9000]
            sage: m = hmm.GaussianMixtureHiddenMarkovModel(
            ....:         [[.9,.1],[.4,.6]],
            ....:         [[(.4,(0,1)), (.6,(1,0.1))], [(1,(0,1))]],
            ....:         [.7,.3])
            sage: m.initial_probabilities()
            [0.7000, 0.3000]
        """
        return TimeSeries(self.pi)

    def transition_matrix(self):
        r"""
        Return the state transition matrix.

        OUTPUT: a Sage matrix with real double precision (RDF) entries.

        EXAMPLES::

            sage: M = hmm.DiscreteHiddenMarkovModel([[0.7,0.3],[0.9,0.1]],
            ....:                                   [[0.5,.5],[.1,.9]],
            ....:                                   [0.3,0.7])
            sage: T = M.transition_matrix(); T
            [0.7 0.3]
            [0.9 0.1]

        The returned matrix is mutable, but changing it does not
        change the transition matrix for the model::

            sage: T[0,0] = .1; T[0,1] = .9
            sage: M.transition_matrix()
            [0.7 0.3]
            [0.9 0.1]

        Transition matrices for other types of models::

            sage: M = hmm.GaussianHiddenMarkovModel([[.1,.9],[.5,.5]],
            ....:                                   [(1,1), (-1,1)],
            ....:                                   [.5,.5])
            sage: M.transition_matrix()
            [0.1 0.9]
            [0.5 0.5]
            sage: M = hmm.GaussianMixtureHiddenMarkovModel(
            ....:         [[.9,.1],[.4,.6]],
            ....:         [[(.4,(0,1)), (.6,(1,0.1))],[(1,(0,1))]],
            ....:         [.7,.3])
            sage: M.transition_matrix()
            [0.9 0.1]
            [0.4 0.6]
        """
        from sage.matrix.constructor import matrix
        from sage.rings.real_double import RDF
        return matrix(RDF, self.N, self.A.list())

    def graph(self, eps=1e-3):
        r"""
        Create a weighted directed graph from the transition matrix,
        not including any edge with a probability less than ``eps``.

        INPUT:

        - ``eps`` -- nonnegative real number

        OUTPUT: a :class:`DiGraph`

        EXAMPLES::

            sage: m = hmm.DiscreteHiddenMarkovModel([[.3,0,.7],[0,0,1],[.5,.5,0]],
            ....:                                   [[.5,.5,.2]]*3,
            ....:                                   [1/3]*3)
            sage: G = m.graph(); G                                                      # needs sage.graphs
            Looped digraph on 3 vertices
            sage: G.edges(sort=True)                                                    # needs sage.graphs
            [(0, 0, 0.3), (0, 2, 0.7), (1, 2, 1.0), (2, 0, 0.5), (2, 1, 0.5)]
            sage: G.plot()                                                              # needs sage.graphs sage.plot
            Graphics object consisting of 11 graphics primitives
        """
        cdef int i, j
        m = self.transition_matrix()
        for i in range(self.N):
            for j in range(self.N):
                if m[i,j] < eps:
                    m[i,j] = 0
        from sage.graphs.digraph import DiGraph
        return DiGraph(m, weighted=True)

    def sample(self, Py_ssize_t length, number=None, starting_state=None):
        r"""
        Return number samples from this HMM of given length.

        INPUT:

        - ``length`` -- positive integer
        - ``number`` -- (default: ``None``) if given, compute list of this many
          sample sequences
        - ``starting_state`` -- integer (or ``None``); if specified, generate a
          sequence using this model starting with the given state instead of
          the initial probabilities to determine the starting state

        OUTPUT:

        - if ``number`` is not given, return a single :class:`TimeSeries`.
        - if ``number`` is given, return a list of :class:`TimeSeries`.

        EXAMPLES::

            sage: set_random_seed(0)
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.1,0.9]],
            ....:                                   [[1,0],[0,1]],
            ....:                                   [0,1])
            sage: print(a.sample(10, 3))
            [[1, 0, 1, 1, 1, 1, 0, 1, 1, 1],
             [1, 1, 0, 0, 1, 1, 1, 1, 1, 1],
             [1, 1, 1, 1, 0, 1, 0, 1, 1, 1]]
            sage: a.sample(15)
            [1, 1, 1, 1, 0 ... 1, 1, 1, 1, 1]
            sage: a.sample(3, 1)
            [[1, 1, 1]]
            sage: list(a.sample(1000)).count(0)
            88

        If the emission symbols are set::

            sage: set_random_seed(0)
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.5,0.5],[0.1,0.9]],
            ....:                                   [[1,0],[0,1]], [0,1],
            ....:                                   ['up', 'down'])
            sage: a.sample(10)
            ['down', 'up', 'down', 'down', 'down', 'down', 'up', 'up', 'up', 'up']

        Force a starting state::

            sage: set_random_seed(0); a.sample(10, starting_state=0)
            ['up', 'up', 'down', 'down', 'down', 'down', 'up', 'up', 'up', 'up']
        """
        if number is None:
            return self.generate_sequence(length, starting_state=starting_state)[0]

        cdef Py_ssize_t i
        return [self.generate_sequence(length, starting_state=starting_state)[0] for i in range(number)]

    #########################################################
    # Some internal functions used for various general
    # HMM algorithms.
    #########################################################
    cdef TimeSeries _baum_welch_gamma(self, TimeSeries alpha, TimeSeries beta):
        r"""
        Used internally to compute the scaled quantity gamma_t(j)
        appearing in the Baum-Welch reestimation algorithm.

        The quantity gamma_t(j) is the (scaled!) probability of being
        in state j at time t, given the observation sequence.

        INPUT:

        - ``alpha`` -- :class:`TimeSeries` as output by the scaled forward algorithm
        - ``beta`` -- :class:`TimeSeries` as output by the scaled backward algorithm

        OUTPUT: :class:`TimeSeries` gamma such that gamma[t*N+j] is gamma_t(j).
        """
        cdef int j, N = self.N
        cdef Py_ssize_t t, T = alpha._length//N
        cdef TimeSeries gamma = TimeSeries(alpha._length, initialize=False)
        cdef double denominator
        sig_on()
        for t in range(T):
            denominator = 0
            for j in range(N):
                gamma._values[t*N + j] = alpha._values[t*N + j] * beta._values[t*N + j]
                denominator += gamma._values[t*N + j]
            for j in range(N):
                gamma._values[t*N + j] /= denominator
        sig_off()
        return gamma


cdef class DiscreteHiddenMarkovModel(HiddenMarkovModel):
    r"""
    A discrete Hidden Markov model implemented using double precision
    floating point arithmetic.

    INPUT:

    - ``A`` -- list of lists or a square `N \times N` matrix, whose
      `(i,j)` entry gives the probability of transitioning from
      state `i` to state `j`.

    - ``B`` -- list of `N` lists or a matrix with `N` rows, such that
      `B[i,k]` gives the probability of emitting symbol `k` while
      in state `i`.

    - ``pi`` -- the probabilities of starting in each initial
      state, i.e., ``pi[i]`` is the probability of starting in
      state `i`.

    - ``emission_symbols`` -- ``None`` or list (default: ``None``); if
      None, the emission_symbols are the ints ``[0..N-1]``, where `N`
      is the number of states.  Otherwise, they are the entries
      of the list ``emissions_symbols``, which must all be hashable.

    - ``normalize`` -- boolean (default: ``True``); if given, input is
      normalized to define valid probability distributions,
      e.g., the entries of `A` are made nonnegative and the rows
      sum to 1, and the probabilities in ``pi`` are normalized.

    EXAMPLES::

        sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]],
        ....:                                   [[0.1,0.9],[0.5,0.5]],
        ....:                                   [.5,.5]); m
        Discrete Hidden Markov Model with 2 States and 2 Emissions
        Transition matrix:
        [0.4 0.6]
        [0.1 0.9]
        Emission matrix:
        [0.1 0.9]
        [0.5 0.5]
        Initial probabilities: [0.5000, 0.5000]
        sage: m.log_likelihood([0,1,0,1,0,1])
        -4.66693474691329...
        sage: m.viterbi([0,1,0,1,0,1])
        ([1, 1, 1, 1, 1, 1], -5.378832842208748)
        sage: m.baum_welch([0,1,0,1,0,1])
        (0.0, 22)
        sage: m  # rel tol 1e-10
        Discrete Hidden Markov Model with 2 States and 2 Emissions
        Transition matrix:
        [1.0134345614745788e-70                    1.0]
        [                   1.0 3.9974352713558623e-19]
        Emission matrix:
        [ 7.380221566254936e-54                    1.0]
        [                   1.0 3.9974352626002193e-19]
        Initial probabilities: [0.0000, 1.0000]
        sage: m.sample(10)
        [0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
        sage: m.graph().plot()                                                          # needs sage.plot
        Graphics object consisting of 6 graphics primitives

    A 3-state model that happens to always outputs 'b'::

        sage: m = hmm.DiscreteHiddenMarkovModel([[1/3]*3]*3, [[0,1,0]]*3, [1/3]*3, ['a','b','c'])
        sage: m.sample(10)
        ['b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b']
    """
    cdef TimeSeries B
    cdef int n_out
    cdef object _emission_symbols, _emission_symbols_dict

    def __init__(self, A, B, pi, emission_symbols=None, bint normalize=True):
        r"""
        Create a discrete emissions HMM with transition probability
        matrix A, emission probabilities given by B, initial state
        probabilities pi, and given emission symbols (which default
        to the first few nonnegative integers).

        EXAMPLES::

            sage: hmm.DiscreteHiddenMarkovModel([.5,0,-1,.5], [[1],[1]],[.5,.5]).transition_matrix()
            [1.0 0.0]
            [0.0 1.0]
            sage: hmm.DiscreteHiddenMarkovModel([0,0,.1,.9], [[1],[1]],[.5,.5]).transition_matrix()
            [0.5 0.5]
            [0.1 0.9]
            sage: hmm.DiscreteHiddenMarkovModel([-1,-2,.1,.9], [[1],[1]],[.5,.5]).transition_matrix()
            [0.5 0.5]
            [0.1 0.9]
            sage: hmm.DiscreteHiddenMarkovModel([1,2,.1,1.2], [[1],[1]],[.5,.5]).transition_matrix()
            [ 0.3333333333333333  0.6666666666666666]
            [0.07692307692307693   0.923076923076923]
        """
        self.pi = util.initial_probs_to_TimeSeries(pi, normalize)
        self.N = len(self.pi)
        self.A = util.state_matrix_to_TimeSeries(A, self.N, normalize)
        self._emission_symbols = emission_symbols
        if self._emission_symbols is not None:
            self._emission_symbols_dict = dict([(y,x) for x,y in enumerate(emission_symbols)])

        if not isinstance(B, Matrix):
            B = matrix(B)
        if B.nrows() != self.N:
            raise ValueError("number of rows of B must equal number of states")
        self.B = TimeSeries(B.list())
        self.n_out = B.ncols()
        if emission_symbols is not None and len(emission_symbols) != self.n_out:
            raise ValueError("number of emission symbols must equal number of output states")
        cdef Py_ssize_t i
        if normalize:
            for i in range(self.N):
                util.normalize_probability_TimeSeries(self.B, i*self.n_out, (i+1)*self.n_out)

    def __reduce__(self):
        r"""
        Used in pickling.

        EXAMPLES::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.0,1.0],[1,1]], [0,1], ['a','b'])
            sage: loads(dumps(m)) == m
            True
        """
        return unpickle_discrete_hmm_v1, \
               (self.A, self.B, self.pi, self.n_out, self._emission_symbols, self._emission_symbols_dict)

    def __richcmp__(self, other, op):
        r"""
        EXAMPLES::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.0,1.0],[0.5,0.5]], [.5,.5])
            sage: n = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.0,1.0],[0.5,0.5]], [1,0])
            sage: m == n
            False
            sage: m == m
            True
            sage: m < n
            True
            sage: n < m
            False
        """
        if not isinstance(other, DiscreteHiddenMarkovModel):
            return NotImplemented
        return PyObject_RichCompare(self.__reduce__()[1],
                                    other.__reduce__()[1], op)

    def emission_matrix(self):
        r"""
        Return the matrix whose `i`-th row specifies the emission
        probability distribution for the `i`-th state.

        More precisely,
        the `i,j` entry of the matrix is the probability of the Markov
        model outputting the `j`-th symbol when it is in the `i`-th state.

        OUTPUT: a Sage matrix with real double precision (RDF) entries.

        EXAMPLES::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]],
            ....:                                   [[0.1,0.9],[0.5,0.5]],
            ....:                                   [.5,.5])
            sage: E = m.emission_matrix(); E
            [0.1 0.9]
            [0.5 0.5]

        The returned matrix is mutable, but changing it does not
        change the transition matrix for the model::

            sage: E[0,0] = 0; E[0,1] = 1
            sage: m.emission_matrix()
            [0.1 0.9]
            [0.5 0.5]
        """
        from sage.matrix.constructor import matrix
        from sage.rings.real_double import RDF
        return matrix(RDF, self.N, self.n_out, self.B.list())

    def __repr__(self):
        r"""
        Return string representation of this discrete hidden Markov model.

        EXAMPLES::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.1,0.9],[0.5,0.5]], [.2,.8])
            sage: m.__repr__()
            'Discrete Hidden Markov Model with 2 States and 2 Emissions\nTransition matrix:\n[0.4 0.6]\n[0.1 0.9]\nEmission matrix:\n[0.1 0.9]\n[0.5 0.5]\nInitial probabilities: [0.2000, 0.8000]'
        """
        s = "Discrete Hidden Markov Model with %s States and %s Emissions"%(
            self.N, self.n_out)
        s += '\nTransition matrix:\n%s'%self.transition_matrix()
        s += '\nEmission matrix:\n%s'%self.emission_matrix()
        s += '\nInitial probabilities: %s'%self.initial_probabilities()
        if self._emission_symbols is not None:
            s += '\nEmission symbols: %s'%self._emission_symbols
        return s

    def _emission_symbols_to_IntList(self, obs):
        r"""
        Internal function used to convert a list of emission symbols to an :class:`IntList`.

        INPUT:

        - ``obs`` -- list of objects

        OUTPUT: an IntList

        EXAMPLES::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.1,0.9],[0.5,0.5]], [.2,.8], ['smile', 'frown'])
            sage: m._emission_symbols_to_IntList(['frown','smile'])
            [1, 0]
        """
        d = self._emission_symbols_dict
        return IntList([d[x] for x in obs])

    def _IntList_to_emission_symbols(self, obs):
        r"""
        Internal function used to convert a list of emission symbols to an IntList.

        INPUT:

        - ``obs`` -- list of objects

        OUTPUT: an IntList

        EXAMPLES::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.1,0.9],[0.5,0.5]], [.2,.8], ['smile', 'frown'])
            sage: m._IntList_to_emission_symbols(stats.IntList([0,0,1,0]))
            ['smile', 'smile', 'frown', 'smile']
        """
        d = self._emission_symbols
        return [d[x] for x in obs]

    def log_likelihood(self, obs, bint scale=True):
        r"""
        Return the logarithm of the probability that this model produced the given
        observation sequence.  Thus the output is a nonpositive number.

        INPUT:

        - ``obs`` -- sequence of observations

        - ``scale`` -- boolean (default: ``True``); if ``True``, use rescaling
          to overoid loss of precision due to the very limited
          dynamic range of floats.  You should leave this as ``True``
          unless the ``obs`` sequence is very small.

        EXAMPLES::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]],
            ....:                                   [[0.1,0.9],[0.5,0.5]],
            ....:                                   [.2,.8])
            sage: m.log_likelihood([0, 1, 0, 1, 1, 0, 1, 0, 0, 0])
            -7.3301308009370825
            sage: m.log_likelihood([0, 1, 0, 1, 1, 0, 1, 0, 0, 0], scale=False)
            -7.330130800937082
            sage: m.log_likelihood([])
            0.0

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]],
            ....:                                   [[0.1,0.9],[0.5,0.5]],
            ....:                                   [.2,.8], ['happy','sad'])
            sage: m.log_likelihood(['happy','happy'])
            -1.6565295199679506
            sage: m.log_likelihood(['happy','sad'])
            -1.4731602941415523

        Overflow from not using the scale option::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]],
            ....:                                   [[0.1,0.9],[0.5,0.5]],
            ....:                                   [.2,.8])
            sage: m.log_likelihood([0,1]*1000, scale=True)
            -1433.820666652728
            sage: m.log_likelihood([0,1]*1000, scale=False)
            -inf
        """
        if len(obs) == 0:
            return 0.0
        if self._emission_symbols is not None:
            obs = self._emission_symbols_to_IntList(obs)
        if not isinstance(obs, IntList):
            obs = IntList(obs)
        if scale:
            return self._forward_scale(obs)
        else:
            return self._forward(obs)

    def _forward(self, IntList obs):
        r"""
        Memory-efficient implementation of the forward algorithm, without
        scaling.

        INPUT:

        - ``obs`` -- integer list of observation states

        OUTPUT: ``float`` -- the log of the probability that the model produced
        this sequence

        EXAMPLES::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.1,0.9],[0.5,0.5]], [.2,.8])
            sage: m._forward(stats.IntList([0,1]*10))
            -14.378663512219456

        The forward algorithm computes the log likelihood::

            sage: m.log_likelihood(stats.IntList([0,1]*10), scale=False)
            -14.378663512219456

        But numerical overflow will happen (without scaling) for long sequences::

            sage: m._forward(stats.IntList([0,1]*1000))
            -inf
        """
        if obs.max() > self.N or obs.min() < 0:
            raise ValueError("invalid observation sequence, since it includes unknown states")

        cdef Py_ssize_t i, j, t, T = len(obs)

        cdef TimeSeries alpha  = TimeSeries(self.N), \
                        alpha2 = TimeSeries(self.N)

        # Initialization
        for i in range(self.N):
            alpha[i] = self.pi[i] * self.B[i*self.n_out + obs._values[0]]
        alpha[i] = self.pi[i] * self.B[i*self.n_out + obs._values[0]]

        # Induction
        cdef double s
        for t in range(1, T):
            for j in range(self.N):
                s = 0
                for i in range(self.N):
                    s += alpha._values[i] * self.A._values[i*self.N + j]
                alpha2._values[j] = s * self.B._values[j*self.n_out+obs._values[t]]
            for j in range(self.N):
                alpha._values[j] = alpha2._values[j]

        # Termination
        return log(alpha.sum())

    def _forward_scale(self, IntList obs):
        r"""
        Memory-efficient implementation of the forward algorithm, with scaling.

        INPUT:

        - ``obs`` -- integer list of observation states

        OUTPUT: ``float`` -- the log of the probability that the model produced
        this sequence

        EXAMPLES::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.1,0.9],[0.5,0.5]], [.2,.8])
            sage: m._forward_scale(stats.IntList([0,1]*10))
            -14.378663512219454

        The forward algorithm computes the log likelihood::

            sage: m.log_likelihood(stats.IntList([0,1]*10))
            -14.378663512219454

        Note that the scale algorithm ensures that numerical overflow
        won't happen for long sequences like it does with the forward
        non-scaled algorithm::

            sage: m._forward_scale(stats.IntList([0,1]*1000))
            -1433.820666652728
            sage: m._forward(stats.IntList([0,1]*1000))
            -inf

        A random sequence produced by the model is more likely::

            sage: set_random_seed(0); v = m.sample(1000)
            sage: m._forward_scale(v)
            -686.8753189365056
        """
        # This is just like self._forward(obs) above, except at every step of the
        # algorithm, we rescale the vector alpha so that the sum of
        # the entries in alpha is 1.  Then, at the end of the
        # algorithm, the sum of probabilities (alpha.sum()) is of
        # course just 1.  However, the true probability that we want
        # is the product of the numbers that we divided by when
        # rescaling, since at each step of the iteration the next term
        # depends linearly on the previous one.  Instead of returning
        # the product, we return the sum of the logs, which avoid
        # numerical error.
        cdef Py_ssize_t i, j, t, T = len(obs)

        # The running some of the log probabilities
        cdef double log_probability = 0, sum, a

        cdef TimeSeries alpha  = TimeSeries(self.N), \
                        alpha2 = TimeSeries(self.N)

        # Initialization
        sum = 0
        for i in range(self.N):
            a = self.pi[i] * self.B[i*self.n_out + obs._values[0]]
            alpha[i] = a
            sum += a

        log_probability = log(sum)
        alpha.rescale(1/sum)

        # Induction
        # The code below is just an optimized version of:
        #     alpha = (alpha * A).pairwise_product(B[O[t+1]])
        #     alpha = alpha.scale(1/alpha.sum())
        # along with keeping track of the log of the scaling factor.
        cdef double s
        for t in range(1, T):
            sum = 0
            for j in range(self.N):
                s = 0
                for i in range(self.N):
                    s += alpha._values[i] * self.A._values[i*self.N + j]
                a = s * self.B._values[j*self.n_out+obs._values[t]]
                alpha2._values[j] = a
                sum += a

            log_probability += log(sum)
            for j in range(self.N):
                alpha._values[j] = alpha2._values[j] / sum

        # Termination
        return log_probability

    def generate_sequence(self, Py_ssize_t length, starting_state=None):
        r"""
        Return a sample of the given length from this HMM.

        INPUT:

        - ``length`` -- positive integer
        - ``starting_state`` -- integer (or ``None``); if specified, generate a
          sequence using this model starting with the given state instead of
          the initial probabilities to determine the starting state

        OUTPUT:

        - an :class:`IntList` or list of emission symbols
        - :class:`IntList` of the actual states the model was in when
          emitting the corresponding symbols

        EXAMPLES:

        In this example, the emission symbols are not set::

            sage: set_random_seed(0)
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.1,0.9]],
            ....:                                   [[1,0],[0,1]],
            ....:                                   [0,1])
            sage: a.generate_sequence(5)
            ([1, 0, 1, 1, 1], [1, 0, 1, 1, 1])
            sage: list(a.generate_sequence(1000)[0]).count(0)
            90

        Here the emission symbols are set::

            sage: set_random_seed(0)
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.5,0.5],[0.1,0.9]],
            ....:                                   [[1,0],[0,1]],
            ....:                                   [0,1], ['up', 'down'])
            sage: a.generate_sequence(5)
            (['down', 'up', 'down', 'down', 'down'], [1, 0, 1, 1, 1])

        Specify the starting state::

            sage: set_random_seed(0); a.generate_sequence(5, starting_state=0)
            (['up', 'up', 'down', 'down', 'down'], [0, 0, 1, 1, 1])
        """
        if length < 0:
            raise ValueError("length must be nonnegative")

        # Create Integer lists for states and observations
        cdef IntList states = IntList(length)
        cdef IntList obs = IntList(length)
        if length == 0:
            # A special case
            if self._emission_symbols is None:
                return states, obs
            else:
                return states, []

        # Setup variables, including random state.
        cdef Py_ssize_t i, j
        cdef randstate rstate = current_randstate()
        cdef int q = 0
        cdef double r, accum
        r = rstate.c_rand_double()

        # Now choose random variables from our discrete distribution.

        # This standard naive algorithm has complexity that is linear
        # in the number of states.  It might be possible to replace it
        # by something that is more efficient.  However, make sure to
        # refactor this code into distribution.pyx before doing so.
        # Note that state switching involves the same algorithm
        # (below).  Use GSL as described here to make this O(1):
        #    http://www.gnu.org/software/gsl/manual/html_node/General-Discrete-Distributions.html

        # Choose initial state:
        if starting_state is None:
            accum = 0
            for i in range(self.N):
                if r < self.pi._values[i] + accum:
                    q = i
                else:
                    accum += self.pi._values[i]
        else:
            q = starting_state
            if q < 0 or q >= self.N:
                raise ValueError("starting state must be between 0 and %s"%(self.N-1))

        states._values[0] = q
        # Generate a symbol from state q
        obs._values[0] = self._gen_symbol(q, rstate.c_rand_double())

        cdef double* row
        cdef int O
        sig_on()
        for i in range(1, length):
            # Choose next state
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
            # Generate symbol from this new state q
            obs._values[i] = self._gen_symbol(q, rstate.c_rand_double())
        sig_off()

        if self._emission_symbols is None:
            # No emission symbol mapping
            return obs, states
        else:
            # Emission symbol mapping, so change our intlist into a list of symbols
            return self._IntList_to_emission_symbols(obs), states

    cdef int _gen_symbol(self, int q, double r) noexcept:
        r"""
        Generate a symbol in state q using the randomly chosen
        floating point number r, which should be between 0 and 1.

        INPUT:

        - ``q`` -- nonnegative integer, which specifies a state
        - ``r`` -- a real number between 0 and 1

        OUTPUT: a nonnegative int

        EXAMPLES::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.1,0.9],[0.9,0.1]], [.2,.8])
            sage: set_random_seed(0)
            sage: m.generate_sequence(10)  # indirect test
            ([0, 1, 0, 0, 0, 0, 1, 0, 1, 1], [1, 0, 1, 1, 1, 1, 0, 0, 0, 0])
        """
        cdef Py_ssize_t j
        cdef double a, accum = 0
        # See the comments above about switching to GSL for this; also note
        # that this should get factored out into a call to something
        # defined in distributions.pyx.
        for j in range(self.n_out):
            a = self.B._values[q*self.n_out + j]
            if r < a + accum:
                return j
            else:
                accum += a
        # None of the values was selected: shouldn't this only happen if the
        # distribution is invalid?   Anyway, this will get factored out.
        return self.n_out - 1

    def viterbi(self, obs, log_scale=True):
        r"""
        Determine "the" hidden sequence of states that is most likely
        to produce the given sequence seq of observations, along with
        the probability that this hidden sequence actually produced
        the observation.

        INPUT:

        - ``seq`` -- sequence of emitted ints or symbols

        - ``log_scale`` -- boolean (default: ``True``); whether to scale the
          sequence in order to avoid numerical overflow

        OUTPUT:

        - ``list`` -- "the" most probable sequence of hidden states, i.e.,
          the Viterbi path

        - ``float`` -- log of probability that the observed sequence
          was produced by the Viterbi sequence of states

        EXAMPLES::

            sage: a = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.1,0.9]],
            ....:                                   [[0.9,0.1],[0.1,0.9]],
            ....:                                   [0.5,0.5])
            sage: a.viterbi([1,0,0,1,0,0,1,1])
            ([1, 0, 0, 1, ..., 0, 1, 1], -11.06245322477221...)

        We predict the state sequence when the emissions are 3/4 and 'abc'.::

            sage: a = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.1,0.9]],
            ....:                                   [[0.9,0.1],[0.1,0.9]],
            ....:                                   [0.5,0.5], [3/4, 'abc'])

        Note that state 0 is common below, despite the model trying hard to
        switch to state 1::

            sage: a.viterbi([3/4, 'abc', 'abc'] + [3/4]*10)
            ([0, 1, 1, 0, 0 ... 0, 0, 0, 0, 0], -25.299405845367794)
        """
        if self._emission_symbols is not None:
            obs = self._emission_symbols_to_IntList(obs)
        elif not isinstance(obs, IntList):
            obs = IntList(obs)
        if log_scale:
            return self._viterbi_scale(obs)
        else:
            return self._viterbi(obs)

    cpdef _viterbi(self, IntList obs):
        r"""
        Used internally to compute the viterbi path, without
        rescaling.  This can be useful for short sequences.

        INPUT:

        - ``obs`` -- IntList

        OUTPUT:

        - IntList (most likely state sequence)

        - log of probability that the observed sequence was
          produced by the Viterbi sequence of states.

        EXAMPLES::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.9,0.1]], [[1,0],[0,1]], [.2,.8])
            sage: m._viterbi(stats.IntList([1]*5))
            ([1, 1, 1, 1, 1], -9.433483923290392)
            sage: m._viterbi(stats.IntList([0]*5))
            ([0, 0, 0, 0, 0], -10.819778284410283)

        Long sequences will overflow::

            sage: m._viterbi(stats.IntList([0]*1000))
            ([0, 0, 0, 0, 0 ... 0, 0, 0, 0, 0], -inf)
        """
        cdef Py_ssize_t t, T = obs._length
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

        # Initialization:
        for i in range(N):
            delta._values[i] = self.pi._values[i] * self.B._values[self.n_out*i + obs._values[0]]

        # Recursion:
        cdef double mx, tmp
        cdef int index
        for t in range(1, T):
            delta_prev, delta = delta, delta_prev
            for j in range(N):
                # delta_t[j] = max_i(delta_{t-1}(i) a_{i,j}) * b_j(obs[t])
                mx = -1
                index = -1
                for i in range(N):
                    tmp = delta_prev._values[i]*self.A._values[i*N+j]
                    if tmp > mx:
                        mx = tmp
                        index = i
                delta._values[j] = mx * self.B._values[self.n_out*j + obs._values[t]]
                psi._values[N*t + j] = index

        # Termination:
        mx, index = delta.max(index=True)

        # Path (state sequence) backtracking:
        state_sequence._values[T-1] = index
        t = T-2
        while t >= 0:
            state_sequence._values[t] = psi._values[N*(t+1) + state_sequence._values[t+1]]
            t -= 1

        return state_sequence, log(mx)

    cpdef _viterbi_scale(self, IntList obs):
        r"""
        Used internally to compute the viterbi path with rescaling.

        INPUT:

        - ``obs`` -- IntList

        OUTPUT:

        - IntList (most likely state sequence)

        - log of probability that the observed sequence was
          produced by the Viterbi sequence of states.

        EXAMPLES::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.9,0.1]], [[.5,.5],[0,1]], [.2,.8])
            sage: m._viterbi_scale(stats.IntList([1]*10))
            ([1, 0, 1, 0, 1, 0, 1, 0, 1, 0], -4.637124095034373)

        Long sequences should not overflow::

            sage: m._viterbi_scale(stats.IntList([1]*1000))
            ([1, 0, 1, 0, 1 ... 0, 1, 0, 1, 0], -452.05188897345...)
        """
        # The algorithm is the same as _viterbi above, except
        # we take the logarithms of everything first, and add
        # instead of multiply.
        cdef Py_ssize_t t, T = obs._length
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
        cdef TimeSeries B = self.B.log()
        cdef TimeSeries pi = self.pi.log()

        # Initialization:
        for i in range(N):
            delta._values[i] = pi._values[i] + B._values[self.n_out*i + obs._values[0]]

        # Recursion:
        cdef double mx, tmp, minus_inf = float('-inf')
        cdef int index

        for t in range(1, T):
            delta_prev, delta = delta, delta_prev
            for j in range(N):
                # delta_t[j] = max_i(delta_{t-1}(i) a_{i,j}) * b_j(obs[t])
                mx = minus_inf
                index = -1
                for i in range(N):
                    tmp = delta_prev._values[i] + A._values[i*N+j]
                    if tmp > mx:
                        mx = tmp
                        index = i
                delta._values[j] = mx + B._values[self.n_out*j + obs._values[t]]
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

    cdef TimeSeries _backward_scale_all(self, IntList obs, TimeSeries scale):
        r"""
        Return the scaled matrix of values `\beta_t(i)` that appear in
        the backtracking algorithm.  This function is used internally
        by the Baum-Welch algorithm.

        The matrix is stored as a TimeSeries T, such that
        `\beta_t(i) = T[t*N + i]` where N is the number of states of
        the Hidden Markov Model.

        The quantity beta_t(i) is the probability of observing the
        sequence obs[t+1:] assuming that the model is in state i at
        time t.

        INPUT:

        - ``obs`` -- IntList
        - ``scale`` -- series that is *changed* in place, so that
          after calling this function, scale[t] is value that is
          used to scale each of the `\beta_t(i)`.

        OUTPUT:

        - a TimeSeries of values beta_t(i).
        - the input object scale is modified
        """
        cdef Py_ssize_t t, T = obs._length
        cdef int N = self.N, i, j
        cdef double s
        cdef TimeSeries beta = TimeSeries(N*T, initialize=False)

        # 1. Initialization
        for i in range(N):
            beta._values[(T-1)*N + i] = 1/scale._values[T-1]

        # 2. Induction
        t = T-2
        while t >= 0:
            for i in range(N):
                s = 0
                for j in range(N):
                    s += self.A._values[i*N+j] * \
                         self.B._values[j*self.n_out+obs._values[t+1]] * beta._values[(t+1)*N+j]
                beta._values[t*N + i] = s/scale._values[t]
            t -= 1
        return beta

    cdef _forward_scale_all(self, IntList obs):
        r"""
        Return scaled values alpha_t(i), the sequence of scalings, and
        the log probability.

        INPUT:

        - ``obs`` -- IntList

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
            a = self.pi._values[i] * self.B._values[i*self.n_out + obs._values[0]]
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
                a = s * self.B._values[j*self.n_out + obs._values[t]]
                alpha._values[t*N + j] = a
                sum += a

            log_probability += log(sum)
            scale._values[t] = sum
            for j in range(self.N):
                alpha._values[t*N + j] /= sum

        # Termination
        return alpha, scale, log_probability

    cdef TimeSeries _baum_welch_xi(self, TimeSeries alpha, TimeSeries beta, IntList obs):
        r"""
        Used internally to compute the scaled quantity xi_t(i,j)
        appearing in the Baum-Welch reestimation algorithm.

        INPUT:

        - ``alpha`` -- TimeSeries as output by the scaled forward algorithm
        - ``beta`` -- TimeSeries as output by the scaled backward algorithm
        - ``obs `` -- IntList of observations

        OUTPUT:

        - TimeSeries xi such that xi[t*N*N + i*N + j] = xi_t(i,j).
        """
        cdef int i, j, N = self.N
        cdef double sum
        cdef Py_ssize_t t, T = alpha._length//N
        cdef TimeSeries xi = TimeSeries(T*N*N, initialize=False)
        sig_on()
        for t in range(T-1):
            sum = 0.0
            for i in range(N):
                for j in range(N):
                    xi._values[t*N*N + i*N + j] = alpha._values[t*N + i]*beta._values[(t+1)*N + j]*\
                       self.A._values[i*N + j] * self.B._values[j*self.n_out + obs._values[t+1]]
                    sum += xi._values[t*N*N + i*N + j]
            for i in range(N):
                for j in range(N):
                    xi._values[t*N*N + i*N + j] /= sum
        sig_off()
        return xi

    def baum_welch(self, obs, int max_iter=100, double log_likelihood_cutoff=1e-4, bint fix_emissions=False):
        r"""
        Given an observation sequence obs, improve this HMM using the
        Baum-Welch algorithm to increase the probability of observing obs.

        INPUT:

        - ``obs`` -- list of emissions

        - ``max_iter`` -- integer (default: 100); maximum number
          of Baum-Welch steps to take

        - ``log_likelihood_cutoff`` -- positive float (default: 1e-4);
          the minimal improvement in likelihood with respect to
          the last iteration required to continue. Relative value
          to log likelihood.

        - ``fix_emissions`` -- boolean (default: ``False``); if ``True``, do
          not change emissions when updating

        OUTPUT:

        changes the model in place, and returns the log
        likelihood and number of iterations.

        EXAMPLES::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.9,0.1]],
            ....:                                   [[.5,.5],[0,1]],
            ....:                                   [.2,.8])
            sage: m.baum_welch([1,0]*20, log_likelihood_cutoff=0)
            (0.0, 4)
            sage: m  # rel tol 1e-14
            Discrete Hidden Markov Model with 2 States and 2 Emissions
            Transition matrix:
            [1.3515269707707603e-51                    1.0]
            [                   1.0                    0.0]
            Emission matrix:
            [                  1.0 6.462537138850569e-52]
            [                  0.0                   1.0]
            Initial probabilities: [0.0000, 1.0000]

        The following illustrates how Baum-Welch is only a local
        optimizer, i.e., the above model is far more likely to produce
        the sequence [1,0]*20 than the one we get below::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.5,0.5],[0.5,0.5]],
            ....:                                   [[.5,.5],[.5,.5]],
            ....:                                   [.5,.5])
            sage: m.baum_welch([1,0]*20, log_likelihood_cutoff=0)
            (-27.725887222397784, 1)
            sage: m
            Discrete Hidden Markov Model with 2 States and 2 Emissions
            Transition matrix:
            [0.5 0.5]
            [0.5 0.5]
            Emission matrix:
            [0.5 0.5]
            [0.5 0.5]
            Initial probabilities: [0.5000, 0.5000]

        We illustrate fixing emissions::

            sage: m = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.9,0.1]],
            ....:                                   [[.5,.5],[.2,.8]],
            ....:                                   [.2,.8])
            sage: set_random_seed(0); v = m.sample(100)
            sage: m.baum_welch(v,fix_emissions=True)
            (-66.98630856918774, 100)
            sage: m.emission_matrix()
            [0.5 0.5]
            [0.2 0.8]
            sage: m = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.9,0.1]],
            ....:                                   [[.5,.5],[.2,.8]],
            ....:                                   [.2,.8])
            sage: m.baum_welch(v)
            (-66.782360659293..., 100)
            sage: m.emission_matrix()  # rel tol 1e-14
            [ 0.5303085748626447 0.46969142513735535]
            [ 0.2909775550173978  0.7090224449826023]
        """
        if self._emission_symbols is not None:
            obs = self._emission_symbols_to_IntList(obs)
        elif not isinstance(obs, IntList):
            obs = IntList(obs)
        cdef IntList _obs = obs
        cdef TimeSeries alpha, beta, scale, gamma, xi
        cdef double log_probability, log_probability_prev, delta
        cdef int i, j, k, N, n_iterations
        cdef Py_ssize_t t, T
        cdef double denominator_A, numerator_A, denominator_B, numerator_B

        # Initialization
        alpha, scale, log_probability = self._forward_scale_all(_obs)
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
                self.pi._values[i] = gamma._values[0*N+i]

            # Reestimate transition matrix and emissions probability in
            # each state.
            for i in range(N):
                denominator_A = 0.0
                for t in range(T-1):
                    denominator_A += gamma._values[t*N+i]
                for j in range(N):
                    numerator_A = 0.0
                    for t in range(T-1):
                        numerator_A += xi._values[t*N*N+i*N+j]
                    self.A._values[i*N+j] = numerator_A / denominator_A

                if not fix_emissions:
                    denominator_B = denominator_A + gamma._values[(T-1)*N + i]
                    for k in range(self.n_out):
                        numerator_B = 0.0
                        for t in range(T):
                            if _obs._values[t] == k:
                                numerator_B += gamma._values[t*N + i]
                        self.B._values[i*self.n_out + k] = numerator_B / denominator_B

            # Initialization
            alpha, scale, log_probability = self._forward_scale_all(_obs)
            beta = self._backward_scale_all(_obs, scale)
            gamma = self._baum_welch_gamma(alpha, beta)
            xi = self._baum_welch_xi(alpha, beta, _obs)

            # Compute the difference between the log probability of
            # two iterations.
            delta = log_probability - log_probability_prev
            log_probability_prev = log_probability
            n_iterations += 1

            # If the log probability does not change by more than
            # delta, then terminate
            if delta <= log_likelihood_cutoff or n_iterations >= max_iter:
                break

        return log_probability, n_iterations


# Keep this -- it's for backwards compatibility with the GHMM based implementation
def unpickle_discrete_hmm_v0(A, B, pi, emission_symbols, name):
    r"""
    TESTS::

        sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.0,1.0],[0.5,0.5]], [1,0])
        sage: sage.stats.hmm.hmm.unpickle_discrete_hmm_v0(m.transition_matrix(), m.emission_matrix(), m.initial_probabilities(), ['a','b'], 'test model')
        Discrete Hidden Markov Model with 2 States and 2 Emissions...
    """
    return DiscreteHiddenMarkovModel(A, B, pi, emission_symbols, normalize=False)


def unpickle_discrete_hmm_v1(A, B, pi, n_out, emission_symbols, emission_symbols_dict):
    r"""
    Return a :class:`DiscreteHiddenMarkovModel`, restored from the arguments.

    This function is used internally for unpickling.

    TESTS::

        sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.0,1.0],[0.5,0.5]], [1,0],['a','b'])
        sage: m2 = loads(dumps(m)) # indirect doctest
        sage: m2 == m
        True

    Test that :issue:`15711` has been resolved::

        sage: str(m2) == str(m)
        True
        sage: m2.log_likelihood('baa'*2) == m.log_likelihood('baa'*2)
        True
    """
    cdef DiscreteHiddenMarkovModel m = DiscreteHiddenMarkovModel.__new__(DiscreteHiddenMarkovModel)
    m.A = A
    m.B = B
    m.pi = pi
    m.N = len(pi)
    m.n_out = n_out
    m._emission_symbols = emission_symbols
    m._emission_symbols_dict = emission_symbols_dict
    return m
