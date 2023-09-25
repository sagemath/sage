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

- Gareth Ma (2023-09-22): implement non-spherical sampling

EXAMPLES::

    sage: from sage.stats.all import DGL
    sage: D = DGL(ZZ^10, 3.0)
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
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.stats.distributions.discrete_gaussian_integer import DGI
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.misc.functional import sqrt
from sage.misc.prandom import normalvariate
from sage.misc.verbose import verbose
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

        sage: from sage.stats.all import DGL
        sage: D = DGL(ZZ^10, 3.0); D
        Discrete Gaussian sampler with Gaussian parameter σ = 3.00000000000000, c=(0, 0, 0, 0, 0, 0, 0, 0, 0, 0) over lattice with basis
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

        sage: from sage.stats.all import DGL
        sage: import warnings
        sage: warnings.simplefilter('ignore', UserWarning)
        sage: D = DGL(identity_matrix(2), 3.0)
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

            sage: from sage.stats.all import DGL
            sage: DGL.compute_precision(100, RR(3))
            100
            sage: DGL.compute_precision(100, RealField(200)(3))
            100
            sage: DGL.compute_precision(100, 3)
            100
            sage: DGL.compute_precision(None, RR(3))
            53
            sage: DGL.compute_precision(None, RealField(200)(3))
            200
            sage: DGL.compute_precision(None, 3)
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

            sage: from sage.stats.all import DGL
            sage: n = 3; sigma = 1.0
            sage: D = DGL(ZZ^n, sigma)
            sage: f = D.f
            sage: nf = D._normalisation_factor_zz(); nf
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

            sage: while abs(m*f(v)*1.0/nf/counter[v] - 1.0) >= 0.1: add_samples(1000)

            sage: v = vector(ZZ, n, (-1, 2, 3))
            sage: v.set_immutable()
            sage: while v not in counter: add_samples(1000)

            sage: while abs(m*f(v)*1.0/nf/counter[v] - 1.0) >= 0.2: add_samples(1000)  # long time

            sage: D = DGL(ZZ^8, 0.5)
            sage: D._normalisation_factor_zz(tau=3)
            3.1653...
            sage: D._normalisation_factor_zz()
            6.8249...
            sage: D = DGL(ZZ^8, 1000)
            sage: round(D._normalisation_factor_zz(prec=100))
            1558545456544038969634991553

            sage: M = Matrix(ZZ, [[1, 3, 0], [-2, 5, 1], [3, -4, 2]])
            sage: D = DGL(M, 1.7)
            sage: D._normalisation_factor_zz() # long time
            7247.1975...

            sage: M = Matrix(ZZ, [[1, 3, 0], [-2, 5, 1]])
            sage: D = DGL(M, 3)
            sage: D._normalisation_factor_zz()
            Traceback (most recent call last):
            ...
            NotImplementedError: Basis must be a square matrix for now.
        """

        # If σ > 1:
        # We use the Fourier transform g(t) of f(x) = exp(-k^2 / 2σ^2), but
        # taking the norm of vector t^2 as input, and with norm_factor factored.
        # If σ ≤ 1:
        # The formula in docstring converges quickly since it has -1 / σ^2 in
        # the exponent
        def f_or_hat(x):
            # Fun fact: If you remove this R() and delay the call to return,
            # It might give an error due to precision error. For example,
            # RR(1 + 100 * exp(-5.0 * pi^2)) == 0

            if sigma > 1:
                return R(exp(-pi**2 * (2 * sigma**2) * x))

            return R(exp(-x / (2 * sigma**2)))

        if not self.is_spherical:
            # NOTE: This is only a poor approximation placeholder.
            # It should be easy to implement, since the Fourier transform
            # is essentially the same, but I can't figure out how to
            # tweak the `.qfrep` call below correctly. TODO.
            from warnings import warn
            warn("Note: `_normalisation_factor_zz` has not been properly "\
                 "implemented for non-spherical distributions.")
            import itertools
            from sage.functions.log import log
            basis = self.B.LLL()
            base = vector(ZZ, [v.round() for v in basis.solve_left(self.c)])
            BOUND = max(1, (self._RR(log(10**4, self.n)).ceil() - 1) // 2)
            if BOUND > 10:
                BOUND = 10
            coords = itertools.product(range(-BOUND, BOUND + 1), repeat=self.n)
            return sum(self.f((vector(u) + base) * self.B) for u in coords)

        if self.B.nrows() != self.B.ncols():
            raise NotImplementedError("Basis must be a square matrix for now.")

        if self.is_spherical and not self._c_in_lattice:
            raise NotImplementedError("Lattice must contain 0 for now.")

        if self.B.base_ring() not in ZZ:
            raise NotImplementedError("Lattice must be integral for now.")

        sigma = self._sigma
        prec = DGL.compute_precision(
            prec, sigma
        )
        R = RealField(prec=prec)
        if sigma > 1:
            det = self.B.det()
            norm_factor = (sigma * sqrt(2 * pi))**self.n / det
        else:
            det = 1
            norm_factor = 1

        # qfrep computes theta series of a quadratic form, which is *half* the
        # generating function of number of vectors with given norm (and no 0)
        Q = self.Q
        if tau is not None:
            freq = Q.__pari__().qfrep(tau * sigma, 0)
            res = R(1)
            for x, fq in enumerate(freq):
                res += 2 * ZZ(fq) * f_or_hat((x + 1) / det**self.n)
            return R(norm_factor * res)

        res = R(1)
        bound = 0
        # There might still be precision issue but whatever
        while True:
            bound += 1
            cnt = ZZ(Q.__pari__().qfrep(bound, 0)[bound - 1])
            inc = 2 * cnt * f_or_hat(bound / det**self.n)
            if cnt > 0 and res == res + inc:
                return R(norm_factor * res)
            res += inc

    @cached_method
    def _maximal_r(self):
        r"""
        This function computes the largest value `r > 0` such that `Σ - r²BBᵀ`
        is positive definite.

        This is equivalent to finding `λ₁(Σ / Q) = 1 / λₙ(Q / Σ)`, which is done
        via the Power iteration method.

        EXAMPLES::

            sage: from sage.stats.all import DGL
            sage: n = 3
            sage: Sigma = Matrix(ZZ, [[5, -2, 4], [-2, 10, -5], [4, -5, 5]])
            sage: c = vector(ZZ, [7, 2, 5])
            sage: D = DGL(ZZ^n, Sigma, c)
            sage: r = D._maximal_r(); r
            0.58402...
            sage: e_vals = (D.sigma - r^2 * D.Q).eigenvalues()
            sage: assert all(e_val >= -1e-12 for e_val in e_vals)
        """
        # TODO: Write doctest
        if self.is_spherical:
            raise RuntimeError("You have encountered a bug. File it! :)")

        Q = self.Q.change_ring(self._RR) / self._sigma.change_ring(self._RR)
        v = Q[0].change_ring(self._RR)
        cnt = 0
        while cnt < 10000:
            nv = (Q * v).normalized()
            if (nv - v).norm() < 1e-12:
                break
            v = nv
            cnt += 1
        res = (v[0] / (Q * v)[0]).sqrt()
        return res

    def _randomise(self, v):
        r"""
        Randomly round to the latice coset `\ZZ + v` with Gaussian parameter `r`

        REFERENCES:

        - [Pei2010]_, Section 4.1
        """
        return vector(ZZ, [DGI(self.r, c=vi)() for vi in v])

    def __init__(self, B, sigma=1, c=None, r=None, precision=None):
        r"""
        Construct a discrete Gaussian sampler over the lattice `Λ(B)`
        with parameter ``sigma`` and center `c`.

        INPUT:

        - ``B`` -- a basis for the lattice, one of the following:

          - an integer matrix,
          - an object with a ``matrix()`` method, e.g. ``ZZ^n``, or
          - an object where ``matrix(B)`` succeeds, e.g. a list of vectors.

        - ``sigma`` -- Gaussian parameter, one of the following:

          - a real number `σ > 0` (spherical),
          - a positive definite matrix `Σ` (non-spherical), or
          - any matrix-like ``S``, equivalent to ``Σ = SSᵀ``

        - ``c`` -- (default: None) center `c`, any vector in `\ZZ^n` is
          supported, but `c ∈ Λ(B)` is faster.

        - ``r`` -- (default: None) rounding parameter `r` as defined in
          [Pei2010]_. Ignored for spherical Gaussian parameter. If not provided,
          set to be the maximal possible such that Σ - r²BBᵀ is positive
          definite.

        - ``precision`` -- bit precision `≥ 53`.

        EXAMPLES::

            sage: from sage.stats.all import DGL
            sage: n = 2; sigma = 3.0
            sage: D = DGL(ZZ^n, sigma)
            sage: f = D.f
            sage: nf = D._normalisation_factor_zz(); nf
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
            sage: while abs(m*f(v)*1.0/nf/counter[v] - 1.0) >= 0.1: add_samples(1000)

            sage: counter = defaultdict(Integer)
            sage: v = vector(ZZ, n, (0, 0))
            sage: v.set_immutable()
            sage: while v not in counter: add_samples(1000)
            sage: while abs(m*f(v)*1.0/nf/counter[v] - 1.0) >= 0.1: add_samples(1000)

        The sampler supports non-spherical covariance in the form of a Gram
        matrix.

            sage: n = 3
            sage: Sigma = Matrix(ZZ, [[5, -2, 4], [-2, 10, -5], [4, -5, 5]])
            sage: c = vector(ZZ, [7, 2, 5])
            sage: D = DGL(ZZ^n, Sigma, c)
            sage: nf = D._normalisation_factor_zz(); nf # This has not been properly implemented
            63.76927...
            sage: while v not in counter: add_samples(1000)
            sage: while abs(m*f(v)*1.0/nf/counter[v] - 1.0) >= 0.1: add_samples(1000)

        The non-spherical sampler supports offline computation to speed up
        sampling. This will be useful when changing the center `c` is supported.
        The difference is more significant for larger matrices. For 128x128 the
        author of this sentence see a 4x speedup (86s -> 20s).

            sage: D.offline_samples = []
            sage: T = 2**12
            sage: L = [D() for _ in range(T)] # 560ms
            sage: D.add_offline_samples(T)    # 150ms
            sage: L = [D() for _ in range(T)] # 370ms

        We can also initialise with matrix-like objects:

            sage: from sage.stats.all import DGL
            sage: qf = QuadraticForm(matrix(3, [2, 1, 1,  1, 2, 1,  1, 1, 2]))
            sage: D = DGL(qf, 3.0); D
            Discrete Gaussian sampler with Gaussian parameter σ = 3.00000000000000, c=(0, 0, 0) over lattice with basis
            <BLANKLINE>
            [2 1 1]
            [1 2 1]
            [1 1 2]
            sage: D().parent() is D.c.parent()
            True
        """
        precision = DGL.compute_precision(precision, sigma)

        self._RR = RealField(precision)
        # Check if sigma is a (real) number or a scaled identity matrix
        self.is_spherical = True
        try:
            self._sigma = self._RR(sigma)
        except TypeError:
            self._sigma = matrix(self._RR, sigma)
            # Will it be "annoying" if a matrix Sigma has different behaviour
            # sometimes? There should be a parameter in the consrtuctor
            if self._sigma == self._sigma[0, 0]:
                self._sigma = self._RR(self._sigma[0, 0])
            else:
                if not self._sigma.is_positive_definite():
                    raise RuntimeError(f"Sigma(={self._sigma}) is not positive definite")
                self.is_spherical = False

        # TODO: Support taking a basis for the covariance
        try:
            B = matrix(B)
        except (TypeError, ValueError):
            pass

        try:
            B = B.matrix()
        except AttributeError:
            pass

        self.n = B.ncols()
        self.B = B
        self.Q = B * B.T
        self._G = B.gram_schmidt()[0]
        self._c_in_lattice = False

        try:
            c = vector(ZZ, self.n, c)
        except TypeError:
            try:
                c = vector(QQ, self.n, c)
            except TypeError:
                c = vector(self._RR, self.n, c)

        self._c = c

        if self.is_spherical:
            # deal with trivial case first, it is common
            if self._G == 1 and self.c == 0:
                self._c_in_lattice = True
                D = DGI(sigma=sigma)
                self.D = tuple([D for _ in range(self.B.nrows())])
                self.VS = FreeModule(ZZ, B.nrows())
                return

            else:
                try:
                    w = B.solve_left(c)
                    if w in ZZ ** B.nrows():
                        self._c_in_lattice = True
                        D = []
                        for i in range(self.B.nrows()):
                            sigma_ = self._sigma / self._G[i].norm()
                            D.append(DGI(sigma=sigma_))
                        self.D = tuple(D)
                        self.VS = FreeModule(ZZ, B.nrows())
                except ValueError:
                    pass
        else:
            # Variables Sigma2 and r are from [Pei2010]_
            # TODO: B is implicitly assumed to be full-rank for the
            # non-spherical case. Remove this assumption :)

            # Offline samples of B⁻¹D₁
            self.offline_samples = []
            self.B_inv = B.inverse()
            self.sigma_inv = self.sigma.inverse()

            if r is None:
                # Compute the maximal r such that (Sigma - r^2 * Q) > 0
                r = self._maximal_r() * 0.9999
            r = self._RR(r)

            Sigma2 = self._sigma - r**2 * self.Q
            try:
                self.r = r
                verbose(f"Computing Cholesky decomposition of a {Sigma2.dimensions()} matrix")
                self.B2 = Sigma2.cholesky().T
                self.B2_B_inv = self.B2 * self.B_inv
            except ValueError:
                raise ValueError("Σ₂ is not positive definite. Is your "\
                                 f"r(={r}) too large? It should be at most "\
                                 f"{self._maximal_r()}")


    def __call__(self):
        r"""
        Return a new sample.

        EXAMPLES::

            sage: from sage.stats.all import DGL
            sage: D = DGL(ZZ^3, 3.0, c=(1,0,0))
            sage: L = [D() for _ in range(2^12)]
            sage: mean_L = sum(L) / len(L)
            sage: norm(mean_L.n() - D.c) < 0.25
            True

            sage: D = DGL(ZZ^3, 3.0, c=(1/2,0,0))
            sage: L = [D() for _ in range(2^12)]  # long time
            sage: mean_L = sum(L) / len(L)        # long time
            sage: norm(mean_L.n() - D.c) < 0.25   # long time
            True
        """
        if not self.is_spherical:
            v = self._call_non_spherical()
        elif self._c_in_lattice:
            v = self._call_in_lattice()
        else:
            v = self._call()
        v.set_immutable()
        return v

    @property
    def f(self):
        r"""
        Returns closure that computes the Gaussian `\rho_{\Lambda, c, \Sigma}`.
        """
        def g(x):
            try:
                x = vector(ZZ, self.n, x)
            except TypeError:
                try:
                    x = vector(QQ, self.n, x)
                except TypeError:
                    x = vector(self._RR, self.n, x)
            x -= self.c
            if self.is_spherical:
                return exp(-x.norm() ** 2 / (2 * self.sigma**2))
            return exp(-x * self.sigma_inv * x / 2)

        return g

    @property
    def sigma(self):
        r"""
        Gaussian parameter `σ`.

        If σ is a real number, samples from this sampler will have expected norm
        `\sqrt{n}σ` where `n` is the dimension of the lattice.

        EXAMPLES::

            sage: from sage.stats.all import DGL
            sage: D = DGL(ZZ^3, 3.0, c=(1,0,0))
            sage: D.sigma
            3.00000000000000
        """
        return self._sigma

    @property
    def c(self):
        r"""
        Center `c`.

        Samples from this sampler will be centered at `c`.

        EXAMPLES::

            sage: from sage.stats.all import DGL
            sage: D = DGL(ZZ^3, 3.0, c=(1,0,0)); D
            Discrete Gaussian sampler with Gaussian parameter σ = 3.00000000000000, c=(1, 0, 0) over lattice with basis
            <BLANKLINE>
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: D.c
            (1, 0, 0)
        """
        return self._c

    @c.setter
    def c(self, _):
        r"""
        Modifies center `c`

        EXAMPLES::

            sage: from sage.stats.all import DGL
            sage: D = DGL(ZZ^3, 3.0, c=(1,0,0))
            sage: D.c = 5
            Traceback (most recent call last):
            ...
            NotImplementedError: Modifying c is not yet supported!
        """
        # TODO: Isolate code to set `c` here, so that the offline part of
        # non-spherical sampling can be effectively utilised
        raise NotImplementedError("Modifying c is not yet supported!")

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from sage.stats.all import DGL
            sage: D = DGL(ZZ^3, 3.0, c=(1,0,0)); D
            Discrete Gaussian sampler with Gaussian parameter σ = 3.00000000000000, c=(1, 0, 0) over lattice with basis
            <BLANKLINE>
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: Sigma = Matrix(ZZ, [[10, -6, 1], [-6, 5, -1], [1, -1, 2]])
            sage: D = DGL(ZZ^3, Sigma); D
            Discrete Gaussian sampler with Gaussian parameter Σ =
            [ 10.0000000000000 -6.00000000000000  1.00000000000000]
            [-6.00000000000000  5.00000000000000 -1.00000000000000]
            [ 1.00000000000000 -1.00000000000000  2.00000000000000], c=(0, 0, 0) over lattice with basis
            <BLANKLINE>
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        # beware of unicode character in ascii string !
        if self.is_spherical:
            sigma_str = f"σ = {self._sigma}"
        else:
            sigma_str = f"Σ =\n{self._sigma}"
        return f"Discrete Gaussian sampler with Gaussian parameter {sigma_str}, c={self.c} over lattice with basis\n\n{self.B}"


    def _call_in_lattice(self):
        r"""
        Return a new sample assuming `c ∈ Λ(B)`.

        EXAMPLES::

            sage: from sage.stats.all import DGL
            sage: D = DGL(ZZ^3, 3.0, c=(1,0,0))
            sage: L = [D._call_in_lattice() for _ in range(2^12)]
            sage: mean_L = sum(L) / len(L)
            sage: norm(mean_L.n() - D.c) < 0.25
            True

        .. note::

           Do not call this method directly, call :func:`DGL.__call__` instead.
        """
        w = self.VS([d() for d in self.D], check=False)
        return w * self.B + self.c

    def _call(self):
        """
        Return a new sample.

        EXAMPLES::

            sage: from sage.stats.all import DGL
            sage: D = DGL(ZZ^3, 3.0, c=(1/2,0,0))
            sage: L = [D._call() for _ in range(2^12)]
            sage: mean_L = sum(L) / len(L)
            sage: norm(mean_L.n() - D.c) < 0.25
            True

        .. note::

           Do not call this method directly, call :func:`DGL.__call__` instead.
        """
        v = 0
        c, sigma, B = self._c, self._sigma, self.B

        m = self.B.nrows()

        for i in range(m - 1, -1, -1):
            b_ = self._G[i]
            c_ = c.dot_product(b_) / b_.dot_product(b_)
            sigma_ = sigma / b_.norm()
            assert sigma_ > 0
            z = DGI(sigma=sigma_, c=c_, algorithm="uniform+online")()
            c = c - z * B[i]
            v = v + z * B[i]
        return v


    def add_offline_samples(self, cnt=1):
        """
        Precompute samples from B⁻¹D₁ to be used in :meth:`_call_non_spherical`

        EXAMPLES::

            sage: from sage.stats.all import DGL
            sage: Sigma = Matrix([[5, -2, 4], [-2, 10, -5], [4, -5, 5]])
            sage: D = DGL(ZZ^3, Sigma)
            sage: assert not D.is_spherical
            sage: D.add_offline_samples(2^12)
            sage: L = [D() for _ in range(2^12)] # Takes less time
        """
        # Just to document the difference with [Pei2010]_, in the paper (Algo 1)
        # he samples from Λ + c, but we instead sample from Λ with distribution
        # sampled at c (D_{Λ, c}), but that's the same as c + D_{Λ - c}
        # Also, we use row notation instead of column notation. Sorry.
        for _ in range(cnt):
            coord = [normalvariate(mu=0, sigma=1) for _ in range(self.n)]
            self.offline_samples.append(vector(self._RR, coord) * self.B2_B_inv)

    def _call_non_spherical(self):
        """
        Return a new sample.

        EXAMPLES::

            sage: from sage.stats.all import DGL
            sage: Sigma = Matrix([[5, -2, 4], [-2, 10, -5], [4, -5, 5]])
            sage: D = DGL(ZZ^3, Sigma, c=(1/2,0,0))
            sage: L = [D._call_non_spherical() for _ in range(2^12)]
            sage: mean_L = sum(L) / len(L)
            sage: norm(mean_L.n() - D.c) < 0.25
            True

        .. note::

           Do not call this method directly, call :func:`DGL.__call__` instead.
        """
        if len(self.offline_samples) == 0:
            self.add_offline_samples()
        vec = self.c * self.B_inv - self.offline_samples.pop()
        return self._randomise(vec) * self.B


DGL = DiscreteGaussianDistributionLatticeSampler
