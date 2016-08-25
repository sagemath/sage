# -*- coding: utf-8 -*-
r"""
Discrete Gaussian Samplers for `\ZZ[x]`

This class realizes oracles which returns polynomials in `\ZZ[x]`
where each coefficient is sampled independently with a probability
proportional to `\exp(-(x-c)²/(2σ²))`.

AUTHORS:

- Martin Albrecht, Robert Fitzpatrick, Daniel Cabracas, Florian Göpfert,
  Michael Schneider: initial version

EXAMPLE::

    sage: from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianDistributionPolynomialSampler
    sage: sigma = 3.0; n=1000
    sage: l = [DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], 64, sigma)() for _ in xrange(n)]
    sage: l = map(lambda f: vector(f).norm().n(), l)
    sage: mean(l), sqrt(64)*sigma
    (23.83..., 24.0...)

"""
from __future__ import absolute_import
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

from sage.rings.all import RealField, RR, ZZ
from .discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
from sage.structure.sage_object import SageObject

class DiscreteGaussianDistributionPolynomialSampler(SageObject):
    r"""
    Discrete Gaussian sampler for polynomials.

    EXAMPLE::

        sage: from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianDistributionPolynomialSampler
        sage: DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], 8, 3.0)()
        3*x^7 + 3*x^6 - 3*x^5 - x^4 - 5*x^2 + 3
        sage: gs = DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], 8, 3.0)
        sage: [gs() for _ in xrange(3)]
        [4*x^7 + 4*x^6 - 4*x^5 + 2*x^4 + x^3 - 4*x + 7, -5*x^6 + 4*x^5 - 3*x^3 + 4*x^2 + x, 2*x^7 + 2*x^6 + 2*x^5 - x^4 - 2*x^2 + 3*x + 1]

    .. automethod:: __init__
    .. automethod:: __call__
    """
    def __init__(self, P, n, sigma):
        r"""
        Construct a sampler for univariate polynomials of degree ``n-1``
        where coefficients are drawn independently with standard deviation
        ``sigma``.

        INPUT:

        - ``P`` - a univariate polynomial ring over the Integers
        - ``n`` - number of coefficients to be sampled
        - ``sigma`` - coefficients `x` are accepted with probability
          proportional to `\exp(-x²/(2σ²))`. If an object of type
          :class:`sage.stats.distributions.discrete_gaussian_integer.DiscreteGaussianDistributionIntegerSampler`
          is passed, then this sampler is used to sample coefficients.

        EXAMPLE::

            sage: from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianDistributionPolynomialSampler
            sage: DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], 8, 3.0)()
            3*x^7 + 3*x^6 - 3*x^5 - x^4 - 5*x^2 + 3
            sage: gs = DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], 8, 3.0)
            sage: [gs() for _ in xrange(3)]
            [4*x^7 + 4*x^6 - 4*x^5 + 2*x^4 + x^3 - 4*x + 7, -5*x^6 + 4*x^5 - 3*x^3 + 4*x^2 + x, 2*x^7 + 2*x^6 + 2*x^5 - x^4 - 2*x^2 + 3*x + 1]
        """
        if isinstance(sigma, DiscreteGaussianDistributionIntegerSampler):
            self.D = sigma
        else:
            self.D = DiscreteGaussianDistributionIntegerSampler(RR(sigma))
        self.n = ZZ(n)
        self.P = P

    def __call__(self):
        """
        Return a new sample.

        EXAMPLE::

            sage: from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianDistributionPolynomialSampler
            sage: sampler = DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], 8, 12.0)
            sage: sampler()
            8*x^7 - 11*x^5 - 19*x^4 + 6*x^3 - 34*x^2 - 21*x + 9
        """
        coeffs = [self.D() for _ in range(self.n)]
        f = self.P(coeffs)
        return f

    def _repr_(self):
        """
        EXAMPLE::

            sage: from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianDistributionPolynomialSampler
            sage: DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], 8, 3.0)
            Discrete Gaussian sampler for polynomials of degree < 8 with σ=3.000000 in each component
        """
        return "Discrete Gaussian sampler for polynomials of degree < %d with σ=%f in each component"%(self.n, self.D.sigma)
