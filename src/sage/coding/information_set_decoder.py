# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Information-set decoding for linear codes

Information-set decoding is a probabilistic decoding strategy that
essentially tries to guess `k` correct positions in the received word,
where `k` is the dimension of the code. A codeword agreeing with the
received word on the guessed position can easily be computed, and their
difference is one possible error vector. A "correct" guess is assumed when
this error vector has low Hamming weight.

This simple algorithm is not very efficient in itself, but there are numerous
refinements to the strategy that make it very capable over rather large codes.
Still, the decoding algorithm is exponential in dimension of the code and the
log of the field size.

The ISD strategy requires choosing how many errors is deemed acceptable. One
choice could be `d/2`, where `d` is the minimum distance of the code, but
sometimes `d` is not known, or sometimes more errors are expected. If one
chooses anything above `d/2`, the algorithm does not guarantee to return a
nearest codeword.

AUTHORS:

- David Lucas, Johan Rosenkilde, Yann Laigle-Chapuy (2016-02, 2017-06): initial
  version
"""

#******************************************************************************
#       Copyright (C) 2017 David Lucas <david.lucas@inria.fr>
#                          Johan Rosenkilde <jsrn@jsrn.dk>
#                          Yann Laigle-Chapuy
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.arith.misc import binomial
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.modules.free_module_element import free_module_element as vector
from sage.structure.sage_object import SageObject
from .decoder import Decoder


def _format_decoding_interval(decoding_interval):
    r"""
    Format the decoding interval of an ISD decoder when calling ``_repr_`` or
    ``_latex_``.

    EXAMPLES::

        sage: from sage.coding.information_set_decoder import _format_decoding_interval
        sage: _format_decoding_interval((0,3))
        'up to 3'
        sage: _format_decoding_interval((2,3))
        'between 2 and 3'
        sage: _format_decoding_interval((3,3))
        'exactly 3'
    """
    if decoding_interval[0] == 0:
        return "up to {0}".format(decoding_interval[1])
    if decoding_interval[0] == decoding_interval[1]:
        return "exactly {0}".format(decoding_interval[0])
    return "between {0} and {1}".format(decoding_interval[0], decoding_interval[1])


class InformationSetAlgorithm(SageObject):
    r"""
    Abstract class for algorithms for
    :class:`sage.coding.information_set_decoder.LinearCodeInformationSetDecoder`.

    To sub-class this class, override ``decode`` and ``calibrate``, and call the
    super constructor from ``__init__``.

    INPUT:

    - ``code`` -- a linear code for which to decode

    - ``number_errors`` -- integer; the maximal number of errors to accept as
      correct decoding. An interval can also be specified by giving a pair of
      integers, where both end values are taken to be in the interval.

    - ``algorithm_name`` -- a name for the specific ISD algorithm used (used for
      printing)

    - ``parameters`` -- (optional) A dictionary for setting the parameters of
      this ISD algorithm. Note that sanity checking this dictionary for the
      individual sub-classes should be done in the sub-class constructor.

    EXAMPLES::

        sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
        sage: LeeBrickellISDAlgorithm(codes.GolayCode(GF(2)), (0,4))
        ISD Algorithm (Lee-Brickell) for [24, 12, 8] Extended Golay code over GF(2)
         decoding up to 4 errors

    A minimal working example of how to sub-class::

        sage: from sage.coding.information_set_decoder import InformationSetAlgorithm
        sage: from sage.coding.decoder import DecodingError
        sage: class MinimalISD(InformationSetAlgorithm):
        ....:     def __init__(self, code, decoding_interval):
        ....:         super().__init__(code, decoding_interval, "MinimalISD")
        ....:     def calibrate(self):
        ....:         self._parameters = { } # calibrate parameters here
        ....:         self._time_estimate = 10.0  # calibrated time estimate
        ....:     def decode(self, r):
        ....:         # decoding algorithm here
        ....:         raise DecodingError("I failed")
        sage: MinimalISD(codes.GolayCode(GF(2)), (0,4))
        ISD Algorithm (MinimalISD) for [24, 12, 8] Extended Golay code over GF(2)
         decoding up to 4 errors
    """

    def __init__(self, code, decoding_interval, algorithm_name, parameters=None):
        r"""
        TESTS::

            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: LeeBrickellISDAlgorithm(codes.GolayCode(GF(2)), (0,4))
            ISD Algorithm (Lee-Brickell) for [24, 12, 8] Extended Golay code over GF(2) decoding up to 4 errors
        """
        self._code = code
        self._decoding_interval = decoding_interval
        self._algorithm_name = algorithm_name
        if parameters:
            self._parameters = parameters
            self._parameters_specified = True
        else:
            self._parameters_specified = False

    def name(self):
        r"""
        Return the name of this ISD algorithm.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: A = LeeBrickellISDAlgorithm(C, (0,2))
            sage: A.name()
            'Lee-Brickell'
        """
        return self._algorithm_name

    def decode(self, r):
        r"""
        Decode a received word using this ISD decoding algorithm.

        Must be overridden by sub-classes.

        EXAMPLES::

            sage: M = matrix(GF(2), [[1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0],\
                                     [0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1],\
                                     [0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0],\
                                     [0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1],\
                                     [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1]])
            sage: C = codes.LinearCode(M)
            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: A = LeeBrickellISDAlgorithm(C, (2,2))
            sage: r = vector(GF(2), [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            sage: A.decode(r)
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        """
        raise NotImplementedError

    def time_estimate(self):
        """
        Estimate for how long this ISD algorithm takes to perform a single decoding.

        The estimate is for a received word whose number of errors is within the
        decoding interval of this ISD algorithm.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: A = LeeBrickellISDAlgorithm(C, (0,2))
            sage: A.time_estimate() #random
            0.0008162108571427874
        """
        if not hasattr(self, "_time_estimate"):
            self.calibrate()
        return self._time_estimate

    def calibrate(self):
        """
        Uses test computations to estimate optimal values for any parameters
        this ISD algorithm may take.

        Must be overridden by sub-classes.

        If ``self._parameters_specified`` is ``False``, this method shall set
        ``self._parameters`` to the best parameters estimated. It shall always
        set ``self._time_estimate`` to the time estimate of using
        ``self._parameters``.

        EXAMPLES::

            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: C = codes.GolayCode(GF(2))
            sage: A = LeeBrickellISDAlgorithm(C, (0,3))
            sage: A.calibrate()
            sage: A.parameters() #random
            {'search_size': 1}
        """
        raise NotImplementedError

    def code(self):
        r"""
        Return the code associated to this ISD algorithm.

        EXAMPLES::

            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: C = codes.GolayCode(GF(2))
            sage: A = LeeBrickellISDAlgorithm(C, (0,3))
            sage: A.code()
            [24, 12, 8] Extended Golay code over GF(2)
        """
        return self._code

    def decoding_interval(self):
        r"""
         A pair of integers specifying the interval of number of errors this
         ISD algorithm will attempt to correct.

         The interval includes both end values.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: A = LeeBrickellISDAlgorithm(C, (0,2))
            sage: A.decoding_interval()
            (0, 2)
        """
        return self._decoding_interval

    def parameters(self):
        """
        Return any parameters this ISD algorithm uses.

        If the parameters have not already been set, efficient values will first
        be calibrated and returned.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: A = LeeBrickellISDAlgorithm(C, (0,4), search_size=3)
            sage: A.parameters()
            {'search_size': 3}

        If not set, calibration will determine a sensible value::

            sage: A = LeeBrickellISDAlgorithm(C, (0,4))
            sage: A.parameters() #random
            {'search_size': 1}
        """
        if not hasattr(self, "_parameters"):
            self.calibrate()
        return self._parameters

    def __eq__(self, other):
        r"""
        Test equality between ISD algorithm objects.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: A = LeeBrickellISDAlgorithm(C, (0,4))
            sage: A == LeeBrickellISDAlgorithm(C, (0,4))
            True
            sage: A == LeeBrickellISDAlgorithm(C, (0,5))
            False
            sage: other_search = 1 if A.parameters()['search_size'] != 1 else 2
            sage: A == LeeBrickellISDAlgorithm(C, (0,4), search_size=other_search)
            False

        ISD Algorithm objects can be equal only if they have both calibrated
        the parameters, or if they both had it set and to the same value::

            sage: A2 = LeeBrickellISDAlgorithm(C, (0,4), search_size=A.parameters()['search_size'])
            sage: A == A2
            False
            sage: A2 == LeeBrickellISDAlgorithm(C, (0,4), search_size=A.parameters()['search_size'])
            True
        """
        return isinstance(other, self.__class__)\
                and self.code() == other.code()\
                and self.decoding_interval() == other.decoding_interval()\
                and self._parameters_specified == other._parameters_specified\
                and (not self._parameters_specified or self.parameters() == other.parameters())

    def __hash__(self):
        r"""
        Return the hash value of ``self``.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: A = LeeBrickellISDAlgorithm(C, (0,4))
            sage: hash(A) #random
            5884357732955478461
            sage: C2 = codes.GolayCode(GF(3))
            sage: A2 = LeeBrickellISDAlgorithm(C2, (0,4))
            sage: hash(A) != hash(A2)
            True
        """
        return hash(str(self))

    def _repr_(self):
        r"""
        Return a string representation of this ISD algorithm.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: A = LeeBrickellISDAlgorithm(C, (0,4))
            sage: A
            ISD Algorithm (Lee-Brickell) for [24, 12, 8] Extended Golay code over GF(2) decoding up to 4 errors
        """
        return "ISD Algorithm ({}) for {} decoding {} errors ".format(self._algorithm_name, self.code(), _format_decoding_interval(self.decoding_interval()))

    def _latex_(self):
        r"""
        Return a latex representation of this ISD algorithm.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: A = LeeBrickellISDAlgorithm(C, (0,4))
            sage: latex(A)
            \textnormal{ISD Algorithm (Lee-Brickell) for }[24, 12, 8] \textnormal{ Extended Golay Code over } \Bold{F}_{2} \textnormal{decoding up to 4 errors}
        """
        return "\\textnormal{{ISD Algorithm ({}) for }}{} \\textnormal{{decoding {} errors}}".format(self._algorithm_name, self.code()._latex_(), _format_decoding_interval(self.decoding_interval()))


class LeeBrickellISDAlgorithm(InformationSetAlgorithm):
    r"""
    The Lee-Brickell algorithm for information-set decoding.

    For a description of the information-set decoding paradigm (ISD), see
    :class:`sage.coding.information_set_decoder.LinearCodeInformationSetDecoder`.

    This implements the Lee-Brickell variant of ISD, see [LB1988]_ for the
    original binary case, and [Pet2010]_ for the `q`-ary extension.

    Let `C` be a `[n, k]`-linear code over `\GF{q}`, and let `r \in \GF{q}^{n}` be
    a received word in a transmission. We seek the codeword whose Hamming
    distance from `r` is minimal. Let `p` and `w` be integers, such that `0\leq
    p\leq w`, Let `G` be a generator matrix of `C`, and for any set of indices
    `I`, we write `G_{I}` for the matrix formed by the columns of `G` indexed by
    `I`. The Lee-Brickell ISD loops the following until it is successful:

        1. Choose an information set `I` of `C`.
        2. Compute `r' = r - r_{I} G_I^{-1}  G`
        3. Consider every size-`p` subset of `I`, `\{a_1, \dots, a_p\}`.
           For each `m = (m_1, \dots, m_p) \in \GF{q}^{p}`, compute
           the error vector `e = r' - \sum_{i=1}^{p} m_i g_{a_i}`,
        4. If `e` has a Hamming weight at most `w`, return `r-e`.

    INPUT:

    - ``code`` -- a linear code for which to decode

    - ``decoding_interval`` -- a pair of integers specifying an interval of
      number of errors to correct; includes both end values

    - ``search_size`` -- (optional) the size of subsets to use on step 3 of the
      algorithm as described above. Usually a small number. It has to be at most
      the largest allowed number of errors. A good choice will be approximated
      if this option is not set; see
      :meth:`sage.coding.LeeBrickellISDAlgorithm.calibrate`
      for details.

    EXAMPLES::

        sage: C = codes.GolayCode(GF(2))
        sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
        sage: A = LeeBrickellISDAlgorithm(C, (0,4)); A
        ISD Algorithm (Lee-Brickell) for [24, 12, 8] Extended Golay code over GF(2)
         decoding up to 4 errors

        sage: C = codes.GolayCode(GF(2))
        sage: A = LeeBrickellISDAlgorithm(C, (2,3)); A
        ISD Algorithm (Lee-Brickell) for [24, 12, 8] Extended Golay code over GF(2)
         decoding between 2 and 3 errors
    """
    def __init__(self, code, decoding_interval, search_size=None):
        r"""
        TESTS:

        If ``search_size`` is not a positive integer, or is bigger than the
        decoding radius, an error will be raised::

            sage: C = codes.GolayCode(GF(2))
            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: LeeBrickellISDAlgorithm(C, (1, 3), search_size=-1)
            Traceback (most recent call last):
            ...
            ValueError: The search size parameter has to be a positive integer

            sage: LeeBrickellISDAlgorithm(C, (1, 3), search_size=4)
            Traceback (most recent call last):
            ...
            ValueError: The search size parameter has to be at most the maximal number of allowed errors
        """
        if search_size is not None:
            if not isinstance(search_size, (Integer, int)) or search_size < 0:
                raise ValueError("The search size parameter has to be a positive integer")
            if search_size > decoding_interval[1]:
                raise ValueError("The search size parameter has to be at most"
                                 " the maximal number of allowed errors")
            super().__init__(code, decoding_interval, "Lee-Brickell",
                             parameters={'search_size': search_size})
            self._parameters_specified = True
        else:
            self._parameters_specified = False
            super().__init__(code, decoding_interval, "Lee-Brickell")

    def decode(self, r):
        r"""
        The Lee-Brickell algorithm as described in the class doc.

        Note that either parameters must be given at construction time or
        :meth:`sage.coding.information_set_decoder.InformationSetAlgorithm.calibrate()`
        should be called before calling this method.

        INPUT:

        - ``r`` -- a received word, i.e. a vector in the ambient space of
          :meth:`decoder.Decoder.code`

        OUTPUT: a codeword whose distance to `r` satisfies ``self.decoding_interval()``.

        EXAMPLES::

            sage: M = matrix(GF(2), [[1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0],\
                                     [0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1],\
                                     [0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0],\
                                     [0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1],\
                                     [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1]])
            sage: C = codes.LinearCode(M)
            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: A = LeeBrickellISDAlgorithm(C, (2,2))
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), 2)
            sage: r = Chan(c)
            sage: c_out = A.decode(r)
            sage: (r - c).hamming_weight() == 2
            True
        """
        import itertools
        from sage.misc.prandom import sample
        C = self.code()
        n, k = C.length(), C.dimension()
        tau = self.decoding_interval()
        p = self.parameters()['search_size']
        F = C.base_ring()
        G = C.generator_matrix()
        Fstar = F.list()[1:]
        while True:
            # step 1.
            I = sample(range(n), k)
            Gi = G.matrix_from_columns(I)
            try:
                Gi_inv = Gi.inverse()
            except ZeroDivisionError:
                # I was not an information set
                continue
            Gt = Gi_inv * G
            # step 2.
            y = r - vector([r[i] for i in I]) * Gt
            g = Gt.rows()
            # step 3.
            for pi in range(p + 1):
                for A in itertools.combinations(range(k), pi):
                    for m in itertools.product(Fstar, repeat=pi):
                        e = y - sum(m[i] * g[A[i]] for i in range(pi))
                        errs = e.hamming_weight()
                        if tau[0] <= errs <= tau[1]:
                            return r - e

    def calibrate(self):
        r"""
        Run some test computations to estimate the optimal search size.

        Let `p` be the search size. We should simply choose `p` such that the
        average expected time is minimal. The algorithm succeeds when it chooses
        an information set with at least `k - p` correct positions, where `k` is
        the dimension of the code and `p` the search size. The expected number
        of trials we need before this occurs is:

        .. MATH::

            \binom{n}{k}/(\rho \sum_{i=0}^p \binom{n-\tau}{k-i} \binom{\tau}{i})

        Here `\rho` is the fraction of `k` subsets of indices which are
        information sets. If `T` is the average time for steps 1 and 2
        (including selecting `I` until an information set is found), while `P(i)`
        is the time for the body of the ``for``-loop in step 3 for `m` of weight
        `i`, then each information set trial takes roughly time `T +
        \sum_{i=0}^{p} P(i) \binom{k}{i} (q-1)^i`, where `\GF{q}` is the base
        field.

        The values `T` and `P` are here estimated by running a few test
        computations similar to those done by the decoding algorithm.
        We don't explicitly estimate `\rho`.

        OUTPUT: does not output anything but sets private fields used by
        :meth:`sage.coding.information_set_decoder.InformationSetAlgorithm.parameters()`
        and
        :meth:`sage.coding.information_set_decoder.InformationSetAlgorithm.time_estimate()`.

        EXAMPLES::

            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: C = codes.GolayCode(GF(2))
            sage: A = LeeBrickellISDAlgorithm(C, (0,3)); A
            ISD Algorithm (Lee-Brickell) for [24, 12, 8] Extended Golay code over GF(2)
             decoding up to 3 errors
            sage: A.calibrate()
            sage: A.parameters() #random
            {'search_size': 1}
            sage: A.time_estimate() #random
            0.0008162108571427874

        If we specify the parameter at construction time, calibrate does not override this choice::

            sage: A = LeeBrickellISDAlgorithm(C, (0,3), search_size=2); A
            ISD Algorithm (Lee-Brickell) for [24, 12, 8] Extended Golay code over GF(2)
             decoding up to 3 errors
            sage: A.parameters()
            {'search_size': 2}
            sage: A.calibrate()
            sage: A.parameters()
            {'search_size': 2}
            sage: A.time_estimate() #random
            0.0008162108571427874
        """
        from sage.matrix.special import random_matrix
        from sage.misc.prandom import sample, randint
        from sage.modules.free_module_element import random_vector
        from time import process_time

        C = self.code()
        G = C.generator_matrix()
        n, k = C.length(), C.dimension()
        tau = self.decoding_interval()[1]
        F = C.base_ring()
        q = F.cardinality()
        Fstar = F.list()[1:]

        def time_information_set_steps():
            before = process_time()
            while True:
                I = sample(range(n), k)
                Gi = G.matrix_from_columns(I)
                try:
                    Gi_inv = Gi.inverse()
                except ZeroDivisionError:
                    continue
                return process_time() - before

        def time_search_loop(p):
            y = random_vector(F, n)
            g = random_matrix(F, p, n).rows()
            scalars = [  [ Fstar[randint(0,q-2)] for i in range(p) ]
                             for s in range(100) ]
            before = process_time()
            for m in scalars:
                e = y - sum(m[i]*g[i] for i in range(p))
            return (process_time() - before) / 100.
        T = sum([ time_information_set_steps() for s in range(5) ]) / 5.
        P = [ time_search_loop(p) for p in range(tau+1) ]

        def compute_estimate(p):
            iters = 1. * binomial(n, k) / \
                sum( binomial(n-tau, k-i)*binomial(tau,i) for i in range(p+1) )
            estimate = iters*(T +
                sum(P[pi] * (q-1)**pi * binomial(k, pi) for pi in range(p+1) ))
            return estimate

        if self._parameters_specified:
            self._time_estimate = compute_estimate(self._parameters['search_size'])
        else:
            self._calibrate_select([ compute_estimate(p) for p in range(tau+1) ])

    def _calibrate_select(self, estimates):
        r"""
        Internal method used by ``self.calibrate()``.

        Given the timing estimates, select the best parameter and set the
        appropriate private fields.

        INPUT:

        - ``estimates`` -- list of time estimates, for the search size set to the
          index of the list entry

        OUTPUT: none, but sets the private fields `self._parameters` and
        `self._time_estimate`.

        TESTS::

            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: C = codes.GolayCode(GF(2))
            sage: A = LeeBrickellISDAlgorithm(C, (0,3)); A
            ISD Algorithm (Lee-Brickell) for [24, 12, 8] Extended Golay code over GF(2) decoding up to 3 errors
            sage: A._calibrate_select([ 1.0, 2.0, 3.0, 0.5, 0.6, 1.0 ])
            sage: A._time_estimate
            0.500000000000000
            sage: A._parameters
            {'search_size': 3}
        """
        search_size = 0
        for p in range(1, len(estimates)):
            if estimates[p] < estimates[search_size]:
                search_size = p
        self._parameters = { 'search_size': search_size }
        self._time_estimate = estimates[search_size]


class LinearCodeInformationSetDecoder(Decoder):
    r"""
    Information-set decoder for any linear code.

    Information-set decoding is a probabilistic decoding strategy that
    essentially tries to guess `k` correct positions in the received word,
    where `k` is the dimension of the code. A codeword agreeing with the
    received word on the guessed position can easily be computed, and their
    difference is one possible error vector. A "correct" guess is assumed when
    this error vector has low Hamming weight.

    The ISD strategy requires choosing how many errors is deemed acceptable. One
    choice could be `d/2`, where `d` is the minimum distance of the code, but
    sometimes `d` is not known, or sometimes more errors are expected. If one
    chooses anything above `d/2`, the algorithm does not guarantee to return a
    nearest codeword.

    This simple algorithm is not very efficient in itself, but there are numerous
    refinements to the strategy. Specifying which strategy to use among those
    that Sage knows is done using the ``algorithm`` keyword. If this is not set,
    an efficient choice will be made for you.

    The various ISD algorithms all need to select a number of parameters. If you
    choose a specific algorithm to use, you can pass these parameters as named
    parameters directly to this class' constructor. If you don't, efficient
    choices will be calibrated for you.

    .. WARNING::

        If there is no codeword within the specified decoding distance, then the
        decoder may never terminate, or it may raise a
        :exc:`sage.coding.decoder.DecodingError` exception, depending on the ISD
        algorithm used.

    INPUT:

    - ``code`` -- a linear code for which to decode

    - ``number_errors`` -- integer; the maximal number of errors to accept as
      correct decoding. An interval can also be specified by giving a pair of
      integers, where both end values are taken to be in the interval.

    - ``algorithm`` -- (optional) the string name of the ISD algorithm to
      employ. If this is not set, an appropriate one will be chosen.
      A constructed
      :class:`sage.coding.information_set_decoder.InformationSetAlgorithm`
      object may also be given. In this case ``number_errors`` must match that
      of the passed algorithm.

    - ``**kwargs`` -- (optional) any number of named arguments passed on to the
      ISD algorithm. Such are usually not required, and they can only be set if
      ``algorithm`` is set to a specific algorithm. See the documentation for
      each individual ISD algorithm class for information on any named arguments
      they may accept. The easiest way to access this documentation is to first
      construct the decoder without passing any named arguments, then accessing
      the ISD algorithm using
      :meth:`sage.coding.information_set_decoder.LinearCodeInformationSetDecoder.algorithm`,
      and then reading the `?` help on the constructed object.

    EXAMPLES:

    The principal way to access this class is through the
    :meth:`sage.code.linear_code.AbstractLinearCode.decoder` method::

        sage: C = codes.GolayCode(GF(3))
        sage: D = C.decoder('InformationSet', 2); D
        Information-set decoder (Lee-Brickell) for [12, 6, 6] Extended Golay code over GF(3)
         decoding up to 2 errors

    You can specify which algorithm you wish to use, and you should do so in
    order to pass special parameters to it::

        sage: C = codes.GolayCode(GF(3))
        sage: D2 = C.decoder('InformationSet', 2, algorithm='Lee-Brickell', search_size=2); D2
        Information-set decoder (Lee-Brickell) for [12, 6, 6] Extended Golay code over GF(3)
         decoding up to 2 errors
        sage: D2.algorithm()
        ISD Algorithm (Lee-Brickell) for [12, 6, 6] Extended Golay code over GF(3)
         decoding up to 2 errors
        sage: D2.algorithm().parameters()
        {'search_size': 2}

    If you specify an algorithm which is not known, you get a friendly error message::

        sage: C.decoder('InformationSet', 2, algorithm="NoSuchThing")
        Traceback (most recent call last):
        ...
        ValueError: Unknown ISD algorithm 'NoSuchThing'.
        The known algorithms are ['Lee-Brickell'].

    You can also construct an ISD algorithm separately and pass that. This is
    mostly useful if you write your own ISD algorithms::

        sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
        sage: A = LeeBrickellISDAlgorithm(C, (0, 2))
        sage: D = C.decoder('InformationSet', 2, algorithm=A); D
        Information-set decoder (Lee-Brickell) for [12, 6, 6] Extended Golay code over GF(3)
         decoding up to 2 errors

    When passing an already constructed ISD algorithm, you can't also pass
    parameters to the ISD algorithm when constructing the decoder::

        sage: C.decoder('InformationSet', 2, algorithm=A, search_size=2)
        Traceback (most recent call last):
        ...
        ValueError: ISD algorithm arguments are not allowed
        when supplying a constructed ISD algorithm

    We can also information-set decode non-binary codes::

        sage: C = codes.GolayCode(GF(3))
        sage: D = C.decoder('InformationSet', 2); D
        Information-set decoder (Lee-Brickell) for [12, 6, 6] Extended Golay code over GF(3)
         decoding up to 2 errors

    There are two other ways to access this class::

        sage: D = codes.decoders.LinearCodeInformationSetDecoder(C, 2); D
        Information-set decoder (Lee-Brickell) for [12, 6, 6] Extended Golay code over GF(3)
         decoding up to 2 errors

        sage: from sage.coding.information_set_decoder import LinearCodeInformationSetDecoder
        sage: D = LinearCodeInformationSetDecoder(C, 2); D
        Information-set decoder (Lee-Brickell) for [12, 6, 6] Extended Golay code over GF(3)
         decoding up to 2 errors
    """
    def __init__(self, code, number_errors, algorithm=None, **kwargs):
        r"""
        TESTS:

        ``number_errors`` has to be either a list of Integers/ints, a tuple of Integers/ints,
        or an Integer/int::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder('InformationSet', "aa")
            Traceback (most recent call last):
            ...
            ValueError: number_errors should be an integer or a pair of integers

        If ``number_errors`` is passed as a list/tuple, it has to contain only
        two values, the first one being at most the second one::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder('InformationSet', (4, 2))
            Traceback (most recent call last):
            ...
            ValueError: number_errors should be a positive integer or a valid interval within the positive integers

        You cannot ask the decoder to correct more errors than the code length::

            sage: D = C.decoder('InformationSet', 25)
            Traceback (most recent call last):
            ...
            ValueError: The provided number of errors should be at most the code's length

        If ``algorithm`` is not set, additional parameters cannot be passed to
        the ISD algorithm::

            sage: D = C.decoder('InformationSet', 2, search_size=2)
            Traceback (most recent call last):
            ...
            ValueError: Additional arguments to an information-set decoder algorithm are only allowed if a specific algorithm is selected by setting the algorithm keyword

        If ``algorithm`` is set to a constructed ISD algorithm, additional
        parameters cannot be passed to the ISD algorithm::

            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: A = LeeBrickellISDAlgorithm(C, (0, 2))
            sage: D = C.decoder('InformationSet', 2, A, search_size=3)
            Traceback (most recent call last):
            ...
            ValueError: ISD algorithm arguments are not allowed when supplying a constructed ISD algorithm

        If ``algorithm`` is set to a constructed
        :class:`sage.coding.information_set_decoder.InformationSetAlgorithm`,
        then ``number_errors`` must match that of the algorithm::

            sage: C = codes.GolayCode(GF(2))
            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: A = LeeBrickellISDAlgorithm(C, (0, 2))
            sage: D = C.decoder('InformationSet', 2, A); D
            Information-set decoder (Lee-Brickell) for [24, 12, 8] Extended Golay code over GF(2) decoding up to 2 errors
            sage: D = C.decoder('InformationSet', (0,2), A); D
            Information-set decoder (Lee-Brickell) for [24, 12, 8] Extended Golay code over GF(2) decoding up to 2 errors
            sage: D = C.decoder('InformationSet', 3, A); D
            Traceback (most recent call last):
            ...
            ValueError: number_errors must match that of the passed ISD algorithm
        """
        if isinstance(number_errors, (Integer, int)):
            number_errors = (0, number_errors)
        if isinstance(number_errors, (tuple, list)) and len(number_errors) == 2 \
            and number_errors[0] in ZZ and number_errors[1] in ZZ:
            if 0 > number_errors[0] or number_errors[0] > number_errors[1]:
                raise ValueError(
                        "number_errors should be a positive integer or"
                        " a valid interval within the positive integers")
            if number_errors[1] > code.length():
                raise ValueError("The provided number of errors should be at"
                                 " most the code's length")
        else:
            raise ValueError("number_errors should be an integer or a pair of integers")

        self._number_errors = number_errors

        super().__init__(code, code.ambient_space(), code._default_encoder_name)

        if algorithm is None:
            if kwargs:
                raise ValueError("Additional arguments to an information-set decoder"
                                " algorithm are only allowed if a specific"
                                " algorithm is selected by setting the algorithm"
                                " keyword")
            algorithm = "Lee-Brickell"
        algorithm_names = LinearCodeInformationSetDecoder.known_algorithms(dictionary=True)

        if isinstance(algorithm, InformationSetAlgorithm):
            if kwargs:
                raise ValueError("ISD algorithm arguments are not allowed when"
                                " supplying a constructed ISD algorithm")
            if number_errors != algorithm.decoding_interval():
                raise ValueError("number_errors must match that of the passed"
                                " ISD algorithm")
            self._algorithm = algorithm
        elif algorithm in algorithm_names:
            self._algorithm = algorithm_names[algorithm](code, number_errors, **kwargs)
        else:
            raise ValueError("Unknown ISD algorithm '{}'."
                            " The known algorithms are {}."
                            .format(algorithm, sorted(algorithm_names)))

    _known_algorithms = {
        "Lee-Brickell": LeeBrickellISDAlgorithm
        }

    @staticmethod
    def known_algorithms(dictionary=False):
        r"""
        Return the list of ISD algorithms that Sage knows.

        Passing any of these to the constructor of
        :class:`sage.coding.information_set_decoder.LinearCodeInformationSetDecoder`
        will make the ISD decoder use that algorithm.

        INPUT:

        - ``dictionary`` -- boolean (default: ``False``); if set to ``True``,
          return a ``dict`` mapping decoding algorithm name to its class

        OUTPUT: list of strings or a ``dict`` from string to ISD algorithm class

        EXAMPLES::

            sage: from sage.coding.information_set_decoder import LinearCodeInformationSetDecoder
            sage: sorted(LinearCodeInformationSetDecoder.known_algorithms())
            ['Lee-Brickell']
        """
        if dictionary:
            return LinearCodeInformationSetDecoder._known_algorithms
        else:
            return LinearCodeInformationSetDecoder._known_algorithms.keys()

    def algorithm(self):
        r"""
        Return the ISD algorithm used by this ISD decoder.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder('InformationSet', (2,4), "Lee-Brickell")
            sage: D.algorithm()
            ISD Algorithm (Lee-Brickell) for [24, 12, 8] Extended Golay code over GF(2)
             decoding between 2 and 4 errors
        """
        return self._algorithm

    def decode_to_code(self, r):
        r"""
        Decodes a received word with respect to the associated code of this decoder.

        .. WARNING::

            If there is no codeword within the decoding radius of this decoder, this
            method may never terminate, or it may raise a
            :exc:`sage.coding.decoder.DecodingError` exception, depending on the ISD
            algorithm used.

        INPUT:

        - ``r`` -- a vector in the ambient space of :meth:`decoder.Decoder.code`

        OUTPUT: a codeword of :meth:`decoder.Decoder.code`

        EXAMPLES::

            sage: M = matrix(GF(2), [[1,0,0,0,0,0,1,0,1,0,1,1,0,0,1],\
                                     [0,1,0,0,0,1,1,1,1,0,0,0,0,1,1],\
                                     [0,0,1,0,0,0,0,1,0,1,1,1,1,1,0],\
                                     [0,0,0,1,0,0,1,0,1,0,0,0,1,1,0],\
                                     [0,0,0,0,1,0,0,0,1,0,1,1,0,1,0]])
            sage: C = LinearCode(M)
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), 2)
            sage: r = Chan(c)
            sage: D = C.decoder('InformationSet', 2)
            sage: c == D.decode_to_code(r)
            True

        Information-set decoding a non-binary code::

            sage: C = codes.GolayCode(GF(3)); C
            [12, 6, 6] Extended Golay code over GF(3)
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), 2)
            sage: r = Chan(c)
            sage: D = C.decoder('InformationSet', 2)
            sage: c == D.decode_to_code(r)
            True

        Let's take a bigger example, for which syndrome decoding or
        nearest-neighbor decoding would be infeasible: the `[59, 30]` Quadratic
        Residue code over `\GF{3}` has true minimum distance 17, so we can
        correct 8 errors::

            sage: C = codes.QuadraticResidueCode(59, GF(3))
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), 2)
            sage: r = Chan(c)
            sage: D = C.decoder('InformationSet', 8)
            sage: c == D.decode_to_code(r) # long time
            True
        """
        C = self.code()
        if r in C:
            return r
        return self.algorithm().decode(r)

    def decoding_radius(self):
        r"""
        Return the maximal number of errors this decoder can decode.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder('InformationSet', 2)
            sage: D.decoding_radius()
            2
        """
        return self._number_errors[1]

    def decoding_interval(self):
        r"""
        A pair of integers specifying the interval of number of errors this
        decoder will attempt to correct.

        The interval includes both end values.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder('InformationSet', 2)
            sage: D.decoding_interval()
            (0, 2)
        """
        return self._number_errors

    def _repr_(self):
        r"""
        Return a string representation of this decoding algorithm.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder('InformationSet', 2)
            sage: D
            Information-set decoder (Lee-Brickell) for [24, 12, 8] Extended Golay code over GF(2) decoding up to 2 errors
        """
        return "Information-set decoder ({}) for {} decoding {} errors ".format(self.algorithm().name(), self.code(), _format_decoding_interval(self.decoding_interval()))

    def _latex_(self):
        r"""
        Return a latex representation of this decoding algorithm.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: from sage.coding.information_set_decoder import LeeBrickellISDAlgorithm
            sage: D = C.decoder('InformationSet', 2)
            sage: latex(D)
            \textnormal{Information-set decoder (Lee-Brickell) for }[24, 12, 8] \textnormal{ Extended Golay Code over } \Bold{F}_{2} \textnormal{decoding up to 2 errors}
        """
        return "\\textnormal{{Information-set decoder ({}) for }}{} \\textnormal{{decoding {} errors}}".format(self.algorithm().name(), self.code()._latex_(), _format_decoding_interval(self.decoding_interval()))


LinearCodeInformationSetDecoder._decoder_type = {"hard-decision",
    "probabilistic", "not-always-closest", "bounded-distance", "might-fail"}
