# -*- coding: utf-8 -*-
r"""
Discrete Gaussian Samplers over Lattices

This file implements oracles which return samples from a lattice following a
discrete Gaussian distribution. That is, if `σ` is big enough relative to the
provided basis, then vectors are returned with a probability proportional to
`\exp(-|x-c|_2^2/(2σ^2))`. More precisely lattice vectors in `x ∈ Λ` are
returned with probability:

    `\exp(-|x-c|_2^2/(2σ²))/(∑_{x ∈ Λ} \exp(-|x|_2^2/(2σ²)))`

AUTHORS:

- Martin Albrecht (2014-06-28): initial version

EXAMPLES::

    sage: from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
    sage: D = DiscreteGaussianDistributionLatticeSampler(ZZ^10, 3.0)
    sage: D(), D(), D()  # random
    ((3, 0, -5, 0, -1, -3, 3, 3, -7, 2), (4, 0, 1, -2, -4, -4, 4, 0, 1, -4), (-3, 0, 4, 5, 0, 1, 3, 2, 0, -1))
    sage: a = D()
    sage: a.parent()
    Ambient free module of rank 10 over the principal ideal domain Integer Ring
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

from sage.functions.log import exp
from sage.rings.real_mpfr import RealField
from sage.rings.real_mpfr import RR
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from .discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
from sage.structure.sage_object import SageObject
from sage.misc.functional import sqrt
from sage.symbolic.constants import pi
from sage.matrix.constructor import matrix
from sage.modules.free_module import FreeModule
from sage.modules.free_module_element import vector


def _iter_vectors(n, lower, upper, step=None):
    r"""
    Iterate over all integer vectors of length ``n`` between ``lower`` and
    ``upper`` bound.

    INPUT:

    - ``n`` - length, integer ``>0``,
    - ``lower`` - lower bound (inclusive), integer ``< upper``.
    - ``upper`` - upper bound (exclusive), integer ``> lower``.
    - ``step`` - used for recursion, ignore.

    EXAMPLES::

      sage: from sage.stats.distributions.discrete_gaussian_lattice import _iter_vectors
      sage: list(_iter_vectors(2, -1, 2))
      [(-1, -1), (0, -1), (1, -1), (-1, 0), (0, 0), (1, 0), (-1, 1), (0, 1), (1, 1)]

    """
    if step is None:
        if ZZ(lower) >= ZZ(upper):
            raise ValueError("Expected lower < upper, but got %d >= %d" % (lower, upper))
        if ZZ(n) <= 0:
            raise ValueError("Expected n>0 but got %d <= 0" % n)
        step = n

    assert step > 0
    if step == 1:
        for x in range(lower, upper):
            v = vector(ZZ, n)
            v[0] = x
            yield v
        return
    else:
        for x in range(lower, upper):
            for v in _iter_vectors(n, lower, upper, step - 1):
                v[step - 1] = x
                yield v


class DiscreteGaussianDistributionLatticeSampler(SageObject):
    r"""
    GPV sampler for Discrete Gaussians over Lattices.

    EXAMPLES::

        sage: from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
        sage: D = DiscreteGaussianDistributionLatticeSampler(ZZ^10, 3.0); D
        Discrete Gaussian sampler with σ = 3.000000, c=(0, 0, 0, 0, 0, 0, 0, 0, 0, 0) over lattice with basis
        <BLANKLINE>
        [1 0 0 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0 0 0]
        [0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 1 0 0 0 0 0]
        [0 0 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 1 0 0 0]
        [0 0 0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 0 0 1 0]
        [0 0 0 0 0 0 0 0 0 1]


    We plot a histogram::

        sage: from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
        sage: import warnings
        sage: warnings.simplefilter('ignore', UserWarning)
        sage: D = DiscreteGaussianDistributionLatticeSampler(identity_matrix(2), 3.0)
        sage: S = [D() for _ in range(2^12)]
        sage: l = [vector(v.list() + [S.count(v)]) for v in set(S)]
        sage: list_plot3d(l, point_list=True, interpolation='nn')
        Graphics3d Object

    REFERENCES:

    - [GPV2008]_

    .. automethod:: __init__
    .. automethod:: __call__
    """
    @staticmethod
    def compute_precision(precision, sigma):
        r"""
        Compute precision to use.

        INPUT:

        - ``precision`` - an integer `>= 53` nor ``None``.
        - ``sigma`` - if ``precision`` is ``None`` then the precision of
          ``sigma`` is used.

        EXAMPLES::

            sage: from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
            sage: DiscreteGaussianDistributionLatticeSampler.compute_precision(100, RR(3))
            100
            sage: DiscreteGaussianDistributionLatticeSampler.compute_precision(100, RealField(200)(3))
            100
            sage: DiscreteGaussianDistributionLatticeSampler.compute_precision(100, 3)
            100
            sage: DiscreteGaussianDistributionLatticeSampler.compute_precision(None, RR(3))
            53
            sage: DiscreteGaussianDistributionLatticeSampler.compute_precision(None, RealField(200)(3))
            200
            sage: DiscreteGaussianDistributionLatticeSampler.compute_precision(None, 3)
            53

        """
        if precision is None:
            try:
                precision = ZZ(sigma.precision())
            except AttributeError:
                return 53
        precision = max(53, precision)
        return precision

    def _normalisation_factor_zz(self, tau=None, prec=None):
        r"""
        This function returns an approximation of `∑_{x ∈ B}
        \exp(-|x|_2^2/(2σ²))`, i.e. the normalisation factor such that the sum
        over all probabilities is 1 for `B`, via Poisson summation.

        INPUT:

        - ``tau`` -- (default: ``None``) all vectors `v` with `|v|_2^2 ≤ τ·σ`
          are enumerated; if none is provided, enumerate vectors with
          increasing norm until the sum converges to given precision. For high
          dimension lattice, this is recommended.

        - ``prec`` -- (default: ``None``) Passed to :meth:`compute_precision`

        EXAMPLES::

            sage: from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
            sage: n = 3; sigma = 1.0
            sage: D = DiscreteGaussianDistributionLatticeSampler(ZZ^n, sigma)
            sage: f = D.f
            sage: c = D._normalisation_factor_zz(); c
            15.7496...

            sage: from collections import defaultdict
            sage: counter = defaultdict(Integer)
            sage: m = 0
            sage: def add_samples(i):
            ....:     global counter, m
            ....:     for _ in range(i):
            ....:         counter[D()] += 1
            ....:         m += 1

            sage: v = vector(ZZ, n, (0, 0, 0))
            sage: v.set_immutable()
            sage: while v not in counter: add_samples(1000)

            sage: while abs(m*f(v)*1.0/c/counter[v] - 1.0) >= 0.1: add_samples(1000)

            sage: v = vector(ZZ, n, (-1, 2, 3))
            sage: v.set_immutable()
            sage: while v not in counter: add_samples(1000)

            sage: while abs(m*f(v)*1.0/c/counter[v] - 1.0) >= 0.2: add_samples(1000)  # long time

            sage: D = DiscreteGaussianDistributionLatticeSampler(ZZ^8, 0.5)
            sage: D._normalisation_factor_zz(tau=3)
            3.1653...
            sage: D._normalisation_factor_zz()
            6.8249...
            sage: D = DiscreteGaussianDistributionLatticeSampler(ZZ^8, 1000)
            sage: round(D._normalisation_factor_zz(prec=100))
            1558545456544038969634991553

            sage: M = Matrix(ZZ, [[1, 3, 0], [-2, 5, 1], [3, -4, 2]])
            sage: D = DiscreteGaussianDistributionLatticeSampler(M, 1.7)
            sage: D._normalisation_factor_zz() # long time
            7247.1975...

            sage: M = Matrix(ZZ, [[1, 3, 0], [-2, 5, 1]])
            sage: D = DiscreteGaussianDistributionLatticeSampler(M, 3)
            sage: D._normalisation_factor_zz()
            Traceback (most recent call last):
            ...
            NotImplementedError: Basis must be a square matrix for now.

            sage: c = vector([3, 7, 1])
            sage: D = DiscreteGaussianDistributionLatticeSampler(M, 3, c = c)
        """

        # If σ > 1:
        # We use the Fourier transform g(t) of f(x) = exp(-k^2 / 2σ^2), but
        # taking the norm of vector t^2 as input, and with norm_factor factored.
        # If σ ≤ 1:
        # The formula in docstring converges quickly since it has -1 / σ^2 in
        # the exponent
        def f(x):
            # Fun fact: If you remove this R() and delay the call to return,
            # It might give an error due to precision error. For example,
            # RR(1 + 100 * exp(-5.0 * pi^2)) == 0
            if sigma > 1:
                return R(exp(-pi**2 * (2 * sigma**2) * x))
            return R(exp(-x / (2 * sigma**2)))

        if self.B.nrows() != self.B.ncols():
            raise NotImplementedError("Basis must be a square matrix for now.")

        if not self._c_in_lattice:
            raise NotImplementedError("Lattice must contain 0 for now.")

        n = self.B.nrows()

        sigma = self._sigma
        prec = DiscreteGaussianDistributionLatticeSampler.compute_precision(
            prec, sigma
        )
        R = RealField(prec=prec)
        if sigma > 1:
            det = self.B.det()
            norm_factor = (sigma * sqrt(2 * pi))**n / det
        else:
            det = 1
            norm_factor = 1

        # qfrep computes theta series of a quadratic form, which is *half* the
        # generating function of number of vectors with given norm (and no 0)
        Q = self.B * self.B.T
        if tau is not None:
            freq = Q.__pari__().qfrep(tau * sigma, 0)
            res = R(1)
            for x, fq in enumerate(freq):
                res += 2 * ZZ(fq) * f((x + 1) / det**n)
            return R(norm_factor * res)

        res = R(1)
        bound = 0
        # There might still be precision issue but whatever
        while True:
            bound += 1
            cnt = ZZ(Q.__pari__().qfrep(bound, 0)[bound - 1])
            inc = 2 * cnt * f(bound / det**n)
            if cnt > 0 and res == res + inc:
                return R(norm_factor * res)
            res += inc

    def __init__(self, B, sigma=1, c=None, precision=None):
        r"""
        Construct a discrete Gaussian sampler over the lattice `Λ(B)`
        with parameter ``sigma`` and center `c`.

        INPUT:

        - ``B`` -- a basis for the lattice, one of the following:

          - an integer matrix,
          - an object with a ``matrix()`` method, e.g. ``ZZ^n``, or
          - an object where ``matrix(B)`` succeeds, e.g. a list of vectors.

        - ``sigma`` -- Gaussian parameter, one of the following:

          - a real number `σ > 0`,
          - a positive definite matrix `Σ`, or
          - any matrix-like ``B``, equivalent to ``Σ = BBᵀ``

        - ``c`` -- center `c`, any vector in `\ZZ^n` is supported, but `c ∈ Λ(B)` is faster.
        - ``precision`` -- bit precision `≥ 53`.

        EXAMPLES::

            sage: from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
            sage: n = 2; sigma = 3.0
            sage: D = DiscreteGaussianDistributionLatticeSampler(ZZ^n, sigma)
            sage: f = D.f
            sage: c = D._normalisation_factor_zz(); c
            56.5486677646...

            sage: from collections import defaultdict
            sage: counter = defaultdict(Integer)
            sage: m = 0
            sage: def add_samples(i):
            ....:     global counter, m
            ....:     for _ in range(i):
            ....:         counter[D()] += 1
            ....:         m += 1

            sage: v = vector(ZZ, n, (-3, -3))
            sage: v.set_immutable()
            sage: while v not in counter: add_samples(1000)
            sage: while abs(m*f(v)*1.0/c/counter[v] - 1.0) >= 0.1: add_samples(1000)

            sage: v = vector(ZZ, n, (0, 0))
            sage: v.set_immutable()
            sage: while v not in counter: add_samples(1000)
            sage: while abs(m*f(v)*1.0/c/counter[v] - 1.0) >= 0.1: add_samples(1000)

            sage: from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
            sage: qf = QuadraticForm(matrix(3, [2, 1, 1,  1, 2, 1,  1, 1, 2]))
            sage: D = DiscreteGaussianDistributionLatticeSampler(qf, 3.0); D
            Discrete Gaussian sampler with σ = 3.000000, c=(0, 0, 0) over lattice with basis
            <BLANKLINE>
            [2 1 1]
            [1 2 1]
            [1 1 2]
            sage: D().parent() is D.c.parent()
            True
        """
        precision = DiscreteGaussianDistributionLatticeSampler.compute_precision(precision, sigma)

        self._RR = RealField(precision)
        # Check if sigma is a (real) number or a scaled identity matrix
        self.is_spherical = True
        try:
            self._sigma = self._RR(sigma)
        except TypeError as e:
            print("error:", e)
            self._sigma = matrix(self._RR, sigma)
            if self._sigma == self._sigma[0, 0]:
                self._sigma = self._RR(self._sigma[0, 0])
            else:
                self.is_spherical = False

        try:
            B = matrix(B)
        except (TypeError, ValueError):
            pass

        try:
            B = B.matrix()
        except AttributeError:
            pass

        self.B = B
        self._G = B.gram_schmidt()[0]

        try:
            c = vector(ZZ, B.ncols(), c)
        except TypeError:
            try:
                c = vector(QQ, B.ncols(), c)
            except TypeError:
                c = vector(RR, B.ncols(), c)

        self._c = c

        self.f = lambda x: exp(-(vector(ZZ, B.ncols(), x) - c).norm() ** 2 / (2 * self._sigma ** 2))

        if self.is_spherical:
            # deal with trivial case first, it is common
            if self._G == 1 and self._c == 0:
                self._c_in_lattice = True
                D = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)
                self.D = tuple([D for _ in range(self.B.nrows())])
                self.VS = FreeModule(ZZ, B.nrows())
                return

            else:
                self._c_in_lattice = False
                try:
                    w = B.solve_left(c)
                    if w in ZZ ** B.nrows():
                        self._c_in_lattice = True
                        D = []
                        for i in range(self.B.nrows()):
                            sigma_ = self._sigma / self._G[i].norm()
                            D.append(DiscreteGaussianDistributionIntegerSampler(sigma=sigma_))
                        self.D = tuple(D)
                        self.VS = FreeModule(ZZ, B.nrows())
                except ValueError:
                    pass
        else:
            # TODO: Precompute basis of sqrt(sigma), change _str_
            raise NotImplementedError

    def __call__(self):
        r"""
        Return a new sample.

        EXAMPLES::

            sage: from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
            sage: D = DiscreteGaussianDistributionLatticeSampler(ZZ^3, 3.0, c=(1,0,0))
            sage: L = [D() for _ in range(2^12)]
            sage: mean_L = sum(L) / len(L)
            sage: norm(mean_L.n() - D.c) < 0.25
            True

            sage: D = DiscreteGaussianDistributionLatticeSampler(ZZ^3, 3.0, c=(1/2,0,0))
            sage: L = [D() for _ in range(2^12)]  # long time
            sage: mean_L = sum(L) / len(L)        # long time
            sage: norm(mean_L.n() - D.c) < 0.25   # long time
            True
        """
        if self._c_in_lattice:
            v = self._call_in_lattice()
        else:
            v = self._call()
        v.set_immutable()
        return v

    @property
    def sigma(self):
        r"""
        Gaussian parameter `σ`.

        Samples from this sampler will have expected norm `\sqrt{n}σ` where `n`
        is the dimension of the lattice.

        EXAMPLES::

            sage: from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
            sage: D = DiscreteGaussianDistributionLatticeSampler(ZZ^3, 3.0, c=(1,0,0))
            sage: D.sigma
            3.00000000000000
        """
        return self._sigma

    @property
    def c(self):
        r"""Center `c`.

        Samples from this sampler will be centered at `c`.

        EXAMPLES::

            sage: from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
            sage: D = DiscreteGaussianDistributionLatticeSampler(ZZ^3, 3.0, c=(1,0,0)); D
            Discrete Gaussian sampler with σ = 3.000000, c=(1, 0, 0) over lattice with basis
            <BLANKLINE>
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: D.c
            (1, 0, 0)
        """
        return self._c

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
            sage: D = DiscreteGaussianDistributionLatticeSampler(ZZ^3, 3.0, c=(1,0,0)); D
            Discrete Gaussian sampler with σ = 3.000000, c=(1, 0, 0) over lattice with basis
            <BLANKLINE>
            [1 0 0]
            [0 1 0]
            [0 0 1]

        """
        # beware of unicode character in ascii string !
        return "Discrete Gaussian sampler with σ = %f, c=%s over lattice with basis\n\n%s" % (self._sigma, self._c, self.B)

    def _call_in_lattice(self):
        r"""
        Return a new sample assuming `c ∈ Λ(B)`.

        EXAMPLES::

            sage: from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
            sage: D = DiscreteGaussianDistributionLatticeSampler(ZZ^3, 3.0, c=(1,0,0))
            sage: L = [D._call_in_lattice() for _ in range(2^12)]
            sage: mean_L = sum(L) / len(L)
            sage: norm(mean_L.n() - D.c) < 0.25
            True

        .. note::

           Do not call this method directly, call :func:`DiscreteGaussianDistributionLatticeSampler.__call__` instead.
        """
        w = self.VS([d() for d in self.D], check=False)
        return w * self.B + self._c

    def _call(self):
        """
        Return a new sample.

        EXAMPLES::

            sage: from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
            sage: D = DiscreteGaussianDistributionLatticeSampler(ZZ^3, 3.0, c=(1/2,0,0))
            sage: L = [D._call() for _ in range(2^12)]  # long time
            sage: mean_L = sum(L) / len(L)              # long time
            sage: norm(mean_L.n() - D.c) < 0.25         # long time
            True

        .. note::

           Do not call this method directly, call :func:`DiscreteGaussianDistributionLatticeSampler.__call__` instead.
        """
        v = 0
        c, sigma, B = self._c, self._sigma, self.B

        m = self.B.nrows()

        for i in range(m - 1, -1, -1):
            b_ = self._G[i]
            c_ = c.dot_product(b_) / b_.dot_product(b_)
            sigma_ = sigma / b_.norm()
            assert sigma_ > 0
            z = DiscreteGaussianDistributionIntegerSampler(sigma=sigma_, c=c_, algorithm="uniform+online")()
            c = c - z * B[i]
            v = v + z * B[i]
        return v

    def _call_non_spherical(self):
        """
        Non-spherical sampler

        .. note::

           Do not call this method directly, call :func:`DiscreteGaussianDistributionLatticeSampler.__call__` instead.
        """
        # TODO: Implement
        raise NotImplementedError
