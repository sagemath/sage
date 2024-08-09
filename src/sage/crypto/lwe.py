# sage.doctest: needs scipy sage.symbolic
"""
(Ring-)LWE oracle generators

The Learning with Errors problem (LWE) is solving linear systems of equations
where the right hand side has been disturbed 'slightly' where 'slightly' is made
precise by a noise distribution - typically a discrete Gaussian
distribution. See [Reg09]_ for details.

The Ring Learning with Errors problem (LWE) is solving a set of univariate
polynomial equations - typically in a cyclotomic field - where the right hand
side was disturbed 'slightly'. See [LPR2010]_ for details.

This module implements generators of LWE samples where parameters are chosen
following proposals in the cryptographic literature.

EXAMPLES:

We get 30 samples from an LWE oracle parameterised by security parameter
``n=20`` and where the modulus and the standard deviation of the noise are
chosen as in [Reg09]_::

    sage: from sage.crypto.lwe import samples
    sage: S = samples(30, 20, 'Regev')
    sage: len(S)
    30
    sage: S[0][0].parent(), S[0][1].parent()
    (Vector space of dimension 20 over Ring of integers modulo 401,
     Ring of integers modulo 401)

We may also pass classes to the samples function, which is useful for users
implementing their own oracles::

    sage: from sage.crypto.lwe import samples, LindnerPeikert
    sage: S = samples(30, 20, LindnerPeikert)
    sage: len(S)
    30
    sage: S[0][0].parent(), S[0][1].parent()
    (Vector space of dimension 20 over Ring of integers modulo 2053,
     Ring of integers modulo 2053)

Finally, :func:`samples` also accepts instances of classes::

    sage: from sage.crypto.lwe import LindnerPeikert
    sage: lwe = LindnerPeikert(20)
    sage: S = samples(30, 20, lwe)
    sage: len(S)
    30
    sage: S[0][0].parent(), S[0][1].parent()
    (Vector space of dimension 20 over Ring of integers modulo 2053,
     Ring of integers modulo 2053)

Note that Ring-LWE samples are returned as vectors::

    sage: # needs sage.libs.pari
    sage: from sage.crypto.lwe import RingLWE
    sage: from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianDistributionPolynomialSampler
    sage: D = DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], euler_phi(16), 5)
    sage: ringlwe = RingLWE(16, 257, D, secret_dist='uniform')
    sage: p = samples(30, euler_phi(16), ringlwe)[0][0].parent(); p
    Vector space of dimension 8 over Ring of integers modulo 257
    sage: assert all(c.parent() is p for b in samples(30, euler_phi(16), ringlwe) for c in b)

One technical issue when working with these generators is that by default they
return vectors and scalars over/in rings modulo some `q`. These are represented
as elements in `(0,q-1)` by Sage. However, it usually is more natural to think
of these entries as integers in `(-q//2,q//2)`. To allow for this, this module
provides the option to balance the representation. In this case vectors and
scalars over/in the integers are returned::

    sage: from sage.crypto.lwe import samples
    sage: for s in samples(30, 20, 'Regev', balanced=True):
    ....:     s1 = list(s[0]) + [s[1]]
    ....:     assert all(-401//2 <= b <= 401//2 for b in s1)

AUTHORS:

- Martin Albrecht
- Robert Fitzpatrick
- Daniel Cabracas
- Florian Göpfert
- Michael Schneider

REFERENCES:

- [Reg09]_

- [LP2011]_

- [LPR2010]_

- [CGW2013]_
"""

from sage.arith.misc import euler_phi, next_prime
from sage.functions.log import log
from sage.functions.other import floor, ceil
from sage.misc.functional import cyclotomic_polynomial, round
from sage.misc.functional import sqrt
from sage.misc.prandom import randint
from sage.misc.randstate import set_random_seed
from sage.modules.free_module import FreeModule
from sage.modules.free_module_element import random_vector, vector
from sage.numerical.optimize import find_root
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianDistributionPolynomialSampler
from sage.structure.element import parent
from sage.structure.sage_object import SageObject
from sage.symbolic.constants import pi
from sage.symbolic.ring import SR


class UniformSampler(SageObject):
    """
    Uniform sampling in a range of integers.

    EXAMPLES::

        sage: from sage.crypto.lwe import UniformSampler
        sage: sampler = UniformSampler(-2, 2); sampler
        UniformSampler(-2, 2)
        sage: sampler() in range(-2, 3)
        True

    .. automethod:: __init__
    .. automethod:: __call__
    """
    def __init__(self, lower_bound, upper_bound):
        """
        Construct a uniform sampler with bounds ``lower_bound`` and
        ``upper_bound`` (both endpoints inclusive).

        INPUT:

        - ``lower_bound`` -- integer
        - ``upper_bound`` -- integer

        EXAMPLES::

            sage: from sage.crypto.lwe import UniformSampler
            sage: UniformSampler(-2, 2)
            UniformSampler(-2, 2)
        """
        if lower_bound > upper_bound:
            raise TypeError("lower bound must be <= upper bound.")
        self.lower_bound = ZZ(lower_bound)
        self.upper_bound = ZZ(upper_bound)

    def __call__(self):
        """
        Return a new sample.

        EXAMPLES::

            sage: from sage.crypto.lwe import UniformSampler
            sage: sampler = UniformSampler(-12, 12)
            sage: sampler() in range(-12, 13)
            True
        """
        return randint(self.lower_bound, self.upper_bound)

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.crypto.lwe import UniformSampler
            sage: UniformSampler(-2, 2)
            UniformSampler(-2, 2)
        """
        return "UniformSampler(%d, %d)" % (self.lower_bound, self.upper_bound)


class UniformPolynomialSampler(SageObject):
    """
    Uniform sampler for polynomials.

    EXAMPLES::

        sage: from sage.crypto.lwe import UniformPolynomialSampler
        sage: UniformPolynomialSampler(ZZ['x'], 8, -2, 2)().parent()
        Univariate Polynomial Ring in x over Integer Ring

    .. automethod:: __init__
    .. automethod:: __call__
    """
    def __init__(self, P, n, lower_bound, upper_bound):
        """
        Construct a sampler for univariate polynomials of degree ``n-1`` where
        coefficients are drawn uniformly at random between ``lower_bound`` and
        ``upper_bound`` (both endpoints inclusive).

        INPUT:

        - ``P`` -- a univariate polynomial ring over the Integers
        - ``n`` -- number of coefficients to be sampled
        - ``lower_bound`` -- integer
        - ``upper_bound`` -- integer

        EXAMPLES::

            sage: from sage.crypto.lwe import UniformPolynomialSampler
            sage: UniformPolynomialSampler(ZZ['x'], 10, -10, 10)
            UniformPolynomialSampler(10, -10, 10)
        """
        self.n = ZZ(n)
        self.P = P
        if lower_bound > upper_bound:
            raise TypeError("lower bound must be <= upper bound.")
        self.lower_bound = ZZ(lower_bound)
        self.upper_bound = ZZ(upper_bound)
        self.D = UniformSampler(self.lower_bound, self.upper_bound)

    def __call__(self):
        """
        Return a new sample.

        EXAMPLES::

            sage: from sage.crypto.lwe import UniformPolynomialSampler
            sage: sampler = UniformPolynomialSampler(ZZ['x'], 8, -12, 12)
            sage: sampler().parent()
            Univariate Polynomial Ring in x over Integer Ring
        """
        coeff = [self.D() for _ in range(self.n)]
        f = self.P(coeff)
        return f

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.crypto.lwe import UniformPolynomialSampler
            sage: UniformPolynomialSampler(ZZ['x'], 8, -3, 3)
            UniformPolynomialSampler(8, -3, 3)
        """
        return "UniformPolynomialSampler(%d, %d, %d)" % (self.n, self.lower_bound, self.upper_bound)


class LWE(SageObject):
    """
    Learning with Errors (LWE) oracle.

    .. automethod:: __init__
    .. automethod:: __call__
    """
    def __init__(self, n, q, D, secret_dist='uniform', m=None):
        r"""
        Construct an LWE oracle in dimension ``n`` over a ring of order
        ``q`` with noise distribution ``D``.

        INPUT:

        - ``n`` -- dimension (integer > 0)
        - ``q`` -- modulus typically > n (integer > 0)
        - ``D`` -- an error distribution such as an instance of
          :class:`DiscreteGaussianDistributionIntegerSampler` or :class:`UniformSampler`
        - ``secret_dist`` -- distribution of the secret (default: ``'uniform'``); one of

          - ``'uniform'`` -- secret follows the uniform distribution in `\Zmod{q}`
          - ``'noise'`` -- secret follows the noise distribution
          - ``(lb, ub)`` -- the secret is chosen uniformly from ``[lb,...,ub]``
            including both endpoints

        - ``m`` -- number of allowed samples or ``None`` if no such limit exists
          (default: ``None``)

        EXAMPLES:

        First, we construct a noise distribution with standard deviation 3.0::

            sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
            sage: D = DiscreteGaussianDistributionIntegerSampler(3.0)

        Next, we construct our oracle::

            sage: from sage.crypto.lwe import LWE
            sage: lwe = LWE(n=20, q=next_prime(400), D=D); lwe
            LWE(20, 401, Discrete Gaussian sampler over the Integers with sigma = 3.000000 and c = 0.000000, 'uniform', None)

        and sample 1000 samples::

            sage: L = []
            sage: def add_samples():
            ....:     global L
            ....:     L += [lwe() for _ in range(100)]
            sage: add_samples()

        To test the oracle, we use the internal secret to evaluate the samples
        in the secret::

            sage: S = lambda : [ZZ(a.dot_product(lwe._LWE__s) - c) for (a,c) in L]

        However, while Sage represents finite field elements between 0 and q-1
        we rely on a balanced representation of those elements here. Hence, we
        fix the representation and recover the correct standard deviation of the
        noise::

            sage: from numpy import std                                                             # needs numpy
            sage: while abs(std([e if e <= 200 else e-401 for e in S()]) - 3.0) > 0.01:             # needs numpy
            ....:     L = []  # reset L to avoid quadratic behaviour
            ....:     add_samples()

        If ``m`` is not ``None`` the number of available samples is restricted::

            sage: from sage.crypto.lwe import LWE
            sage: lwe = LWE(n=20, q=next_prime(400), D=D, m=30)
            sage: _ = [lwe() for _ in range(30)]
            sage: lwe()  # 31
            Traceback (most recent call last):
            ...
            IndexError: Number of available samples exhausted.
        """
        self.n = ZZ(n)
        self.m = m
        self.__i = 0
        self.K = IntegerModRing(q)
        self.FM = FreeModule(self.K, n)
        self.D = D

        self.secret_dist = secret_dist
        if secret_dist == 'uniform':
            self.__s = random_vector(self.K, self.n)
        elif secret_dist == 'noise':
            self.__s = vector(self.K, self.n, [self.D() for _ in range(n)])
        else:
            try:
                lb, ub = map(ZZ, secret_dist)
                self.__s = vector(self.K, self.n, [randint(lb,ub) for _ in range(n)])
            except (IndexError, TypeError):
                raise TypeError("Parameter secret_dist=%s not understood." % (secret_dist))

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
            sage: from sage.crypto.lwe import LWE
            sage: D = DiscreteGaussianDistributionIntegerSampler(3.0)
            sage: lwe = LWE(n=20, q=next_prime(400), D=D); lwe
            LWE(20, 401, Discrete Gaussian sampler over the Integers with sigma = 3.000000 and c = 0.000000, 'uniform', None)

            sage: lwe = LWE(n=20, q=next_prime(400), D=D, secret_dist=(-3, 3)); lwe
            LWE(20, 401, Discrete Gaussian sampler over the Integers with sigma = 3.000000 and c = 0.000000, (-3, 3), None)
        """
        if isinstance(self.secret_dist, str):
            return "LWE(%d, %d, %s, '%s', %s)" % (self.n,self.K.order(),self.D,self.secret_dist, self.m)
        else:
            return "LWE(%d, %d, %s, %s, %s)" % (self.n,self.K.order(),self.D,self.secret_dist, self.m)

    def __call__(self):
        """
        EXAMPLES::

            sage: from sage.crypto.lwe import DiscreteGaussianDistributionIntegerSampler, LWE
            sage: LWE(10, 401, DiscreteGaussianDistributionIntegerSampler(3))()[0].parent()
            Vector space of dimension 10 over Ring of integers modulo 401
            sage: LWE(10, 401, DiscreteGaussianDistributionIntegerSampler(3))()[1].parent()
            Ring of integers modulo 401
        """
        if self.m is not None:
            if self.__i >= self.m:
                raise IndexError("Number of available samples exhausted.")
        self.__i += 1
        a = self.FM.random_element()
        return a, a.dot_product(self.__s) + self.K(self.D())


class Regev(LWE):
    """
    LWE oracle with parameters as in [Reg09]_.

    .. automethod:: __init__
    """
    def __init__(self, n, secret_dist='uniform', m=None):
        """
        Construct LWE instance parameterised by security parameter ``n`` where
        the modulus ``q`` and the ``stddev`` of the noise are chosen as in
        [Reg09]_.

        INPUT:

        - ``n`` -- security parameter (integer > 0)
        - ``secret_dist`` -- distribution of the secret. See documentation of :class:`LWE`
          for details (default='uniform')
        - ``m`` -- number of allowed samples or ``None`` if no such limit exists
          (default: ``None``)

        EXAMPLES::

            sage: from sage.crypto.lwe import Regev
            sage: Regev(n=20)
            LWE(20, 401, Discrete Gaussian sampler over the Integers with sigma = 1.915069 and c = 401.000000, 'uniform', None)
        """
        q = ZZ(next_prime(n**2))
        s = RR(1/(RR(n).sqrt() * log(n, 2)**2) * q)
        D = DiscreteGaussianDistributionIntegerSampler(s/sqrt(2*pi.n()), q)
        LWE.__init__(self, n=n, q=q, D=D, secret_dist=secret_dist, m=m)

class LindnerPeikert(LWE):
    """
    LWE oracle with parameters as in [LP2011]_.

    .. automethod:: __init__
    """
    def __init__(self, n, delta=0.01, m=None):
        """
        Construct LWE instance parameterised by security parameter ``n`` where
        the modulus ``q`` and the ``stddev`` of the noise is chosen as in
        [LP2011]_.

        INPUT:

        - ``n`` -- security parameter (integer > 0)
        - ``delta`` -- error probability per symbol (default: 0.01)
        - ``m`` -- number of allowed samples or ``None`` in which case ``m=2*n +
          128`` as in [LP2011]_ (default: ``None``)

        EXAMPLES::

            sage: from sage.crypto.lwe import LindnerPeikert
            sage: LindnerPeikert(n=20)
            LWE(20, 2053, Discrete Gaussian sampler over the Integers with sigma = 3.600954 and c = 0.000000, 'noise', 168)
        """
        if m is None:
            m = 2*n + 128
        # Find c>=1 such that c*exp((1-c**2)/2))**(2*n) == 2**-40
        #         (c*exp((1-c**2)/2))**(2*n) == 2**-40
        #    log((c*exp((1-c**2)/2))**(2*n)) == -40*log(2)
        #       (2*n)*log(c*exp((1-c**2)/2)) == -40*log(2)
        #  2*n*(log(c)+log(exp((1-c**2)/2))) == -40*log(2)
        #            2*n*(log(c)+(1-c**2)/2) == -40*log(2)
        #              2*n*log(c)+n*(1-c**2) == -40*log(2)
        #  2*n*log(c)+n*(1-c**2) + 40*log(2) == 0
        c = SR.var('c')
        c = find_root(2*n*log(c)+n*(1-c**2) + 40*log(2) == 0, 1, 10)
        # Upper bound on s**2/t
        s_t_bound = (sqrt(2) * pi / c / sqrt(2*n*log(2/delta))).n()
        # Interpretation of "choose q just large enough to allow for a Gaussian parameter s>=8" in [LP2011]_
        q = next_prime(floor(2**round(log(256 / s_t_bound, 2))))
        # Gaussian parameter as defined in [LP2011]_
        s = sqrt(s_t_bound*floor(q/4))
        # Transform s into stddev
        stddev = s/sqrt(2*pi.n())
        D = DiscreteGaussianDistributionIntegerSampler(stddev)
        LWE.__init__(self, n=n, q=q, D=D, secret_dist='noise', m=m)


class UniformNoiseLWE(LWE):
    """
    LWE oracle with uniform secret with parameters as in [CGW2013]_.

    .. automethod:: __init__
    """
    def __init__(self, n, instance='key', m=None):
        """
        Construct LWE instance parameterised by security parameter ``n`` where
        all other parameters are chosen as in [CGW2013]_.

        INPUT:

        - ``n`` -- security parameter (integer >= 89)
        - ``instance`` -- one of

          - ``'key'`` -- the LWE-instance that hides the secret key is generated
          - ``'encrypt'`` -- the LWE-instance that hides the message is generated
            (default: ``'key'``)

        - ``m`` -- number of allowed samples or ``None`` in which case ``m`` is
          chosen as in [CGW2013]_.  (default: ``None``)

        EXAMPLES::

            sage: from sage.crypto.lwe import UniformNoiseLWE
            sage: UniformNoiseLWE(89)
            LWE(89, 64311834871, UniformSampler(0, 6577), 'noise', 131)

            sage: UniformNoiseLWE(89, instance='encrypt')
            LWE(131, 64311834871, UniformSampler(0, 11109), 'noise', 181)
        """

        if n < 89:
            raise TypeError("Parameter too small")

        n2 = n
        C = 4/sqrt(2*pi)
        kk = floor((n2-2*log(n2, 2)**2)/5)
        n1 = (3*n2-5*kk) // 2
        ke = floor((n1-2*log(n1, 2)**2)/5)
        l = (3*n1-5*ke) // 2 - n2
        sk = ceil((C*(n1+n2))**(ZZ(3)/2))
        se = ceil((C*(n1+n2+l))**(ZZ(3)/2))
        q = next_prime(max(ceil((4*sk)**(ZZ(n1+n2)/n1)),
                           ceil((4*se)**(ZZ(n1+n2+l)/(n2+l))),
                           ceil(4*(n1+n2)*se*sk+4*se+1)))

        if kk <= 0:
            raise TypeError("Parameter too small")

        if instance == 'key':
            D = UniformSampler(0, sk-1)
            if m is None:
                m = n1
            LWE.__init__(self, n=n2, q=q, D=D, secret_dist='noise', m=m)
        elif instance == 'encrypt':
            D = UniformSampler(0, se-1)
            if m is None:
                m = n2+l
            LWE.__init__(self, n=n1, q=q, D=D, secret_dist='noise', m=m)
        else:
            raise TypeError("Parameter instance=%s not understood." % (instance))

class RingLWE(SageObject):
    """
    Ring Learning with Errors oracle.

    .. automethod:: __init__
    .. automethod:: __call__
    """
    def __init__(self, N, q, D, poly=None, secret_dist='uniform', m=None):
        """
        Construct a Ring-LWE oracle in dimension ``n=phi(N)`` over a ring of order
        ``q`` with noise distribution ``D``.

        INPUT:

        - ``N`` -- index of cyclotomic polynomial (integer > 0, must be power of 2)
        - ``q`` -- modulus typically > N (integer > 0)
        - ``D`` -- an error distribution such as an instance of
          :class:`DiscreteGaussianDistributionPolynomialSampler` or :class:`UniformSampler`
        - ``poly`` -- a polynomial of degree ``phi(N)``. If ``None`` the
          cyclotomic polynomial used (default: ``None``).
        - ``secret_dist`` -- distribution of the secret. See documentation of
          :class:`LWE` for details (default='uniform')
        - ``m`` -- number of allowed samples or ``None`` if no such limit exists
          (default: ``None``)

        EXAMPLES::

            sage: from sage.crypto.lwe import RingLWE
            sage: from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianDistributionPolynomialSampler
            sage: D = DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], n=euler_phi(20), sigma=3.0)                # needs sage.libs.pari
            sage: RingLWE(N=20, q=next_prime(800), D=D)                                 # needs sage.libs.pari
            RingLWE(20, 809, Discrete Gaussian sampler for polynomials of degree < 8 with σ=3.000000 in each component, x^8 - x^6 + x^4 - x^2 + 1, 'uniform', None)
        """
        self.N = ZZ(N)
        self.n = euler_phi(N)
        self.m = m
        self.__i = 0
        self.K = IntegerModRing(q)

        if self.n != D.n:
            raise ValueError("Noise distribution has dimensions %d != %d" % (D.n, self.n))

        self.D = D
        self.q = q
        if poly is not None:
            self.poly = poly
        else:
            self.poly = cyclotomic_polynomial(self.N, 'x')

        self.R_q = self.K['x'].quotient(self.poly, 'x')

        self.secret_dist = secret_dist
        if secret_dist == 'uniform':
            self.__s = self.R_q.random_element()  # uniform sampling of secret
        elif secret_dist == 'noise':
            self.__s = self.D()
        else:
            raise TypeError("Parameter secret_dist=%s not understood." % (secret_dist))

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.crypto.lwe import DiscreteGaussianDistributionPolynomialSampler, RingLWE
            sage: D = DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], n=8, sigma=3.0)
            sage: RingLWE(N=16, q=next_prime(400), D=D)
            RingLWE(16, 401, Discrete Gaussian sampler for polynomials of degree < 8 with σ=3.000000 in each component, x^8 + 1, 'uniform', None)
        """
        if isinstance(self.secret_dist, str):
            return "RingLWE(%d, %d, %s, %s, '%s', %s)" % (self.N, self.K.order(), self.D, self.poly, self.secret_dist, self.m)
        else:
            return "RingLWE(%d, %d, %s, %s, %s, %s)" % (self.N, self.K.order(), self.D, self.poly, self.secret_dist, self.m)

    def __call__(self):
        """
        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: from sage.crypto.lwe import DiscreteGaussianDistributionPolynomialSampler, RingLWE
            sage: N = 16
            sage: n = euler_phi(N)
            sage: D = DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], n, 5)
            sage: ringlwe = RingLWE(N, 257, D, secret_dist='uniform')
            sage: ringlwe()[0].parent()
            Vector space of dimension 8 over Ring of integers modulo 257
            sage: ringlwe()[1].parent()
            Vector space of dimension 8 over Ring of integers modulo 257
        """
        if self.m is not None:
            if self.__i >= self.m:
                raise IndexError("Number of available samples exhausted.")
        self.__i += 1
        a = self.R_q.random_element()
        return vector(a), vector(a * (self.__s) + self.D())

class RingLindnerPeikert(RingLWE):
    """
    Ring-LWE oracle with parameters as in [LP2011]_.

    .. automethod:: __init__
    """
    def __init__(self, N, delta=0.01, m=None):
        """
        Construct a Ring-LWE oracle in dimension ``n=phi(N)`` where
        the modulus ``q`` and the ``stddev`` of the noise is chosen as in
        [LP2011]_.

        INPUT:

        - ``N`` -- index of cyclotomic polynomial (integer > 0, must be power of 2)
        - ``delta`` -- error probability per symbol (default: 0.01)
        - ``m`` -- number of allowed samples or ``None`` in which case ``3*n`` is
          used (default: ``None``)

        EXAMPLES::

            sage: from sage.crypto.lwe import RingLindnerPeikert
            sage: RingLindnerPeikert(N=16)
            RingLWE(16, 1031, Discrete Gaussian sampler for polynomials of degree < 8 with σ=2.803372 in each component, x^8 + 1, 'noise', 24)
        """
        n = euler_phi(N)
        if m is None:
            m = 3*n
        # Find c>=1 such that c*exp((1-c**2)/2))**(2*n) == 2**-40
        #  i.e c>=1 such that 2*n*log(c)+n*(1-c**2) + 40*log(2) == 0
        c = SR.var('c')
        c = find_root(2*n*log(c)+n*(1-c**2) + 40*log(2) == 0, 1, 10)
        # Upper bound on s**2/t
        s_t_bound = (sqrt(2) * pi / c / sqrt(2*n*log(2/delta))).n()
        # Interpretation of "choose q just large enough to allow for a Gaussian parameter s>=8" in [LP2011]_
        q = next_prime(floor(2**round(log(256 / s_t_bound, 2))))
        # Gaussian parameter as defined in [LP2011]_
        s = sqrt(s_t_bound*floor(q/4))
        # Transform s into stddev
        stddev = s/sqrt(2*pi.n())
        D = DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], n, stddev)
        RingLWE.__init__(self, N=N, q=q, D=D, poly=None, secret_dist='noise', m=m)

class RingLWEConverter(SageObject):
    """
    Wrapper callable to convert Ring-LWE oracles into LWE oracles by
    disregarding the additional structure.

    .. automethod:: __init__
    .. automethod:: __call__
    """
    def __init__(self, ringlwe):
        """
        INPUT:

        - ``ringlwe`` -- an instance of a :class:`RingLWE`

        EXAMPLES::

            sage: from sage.crypto.lwe import DiscreteGaussianDistributionPolynomialSampler, RingLWE, RingLWEConverter
            sage: D = DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], euler_phi(16), 5)
            sage: lwe = RingLWEConverter(RingLWE(16, 257, D, secret_dist='uniform'))
            sage: set_random_seed(1337)
            sage: lwe()
            ((171, 197, 58, 125, 3, 216, 32, 130), ...)
        """
        self.ringlwe = ringlwe
        self._i = 0
        self._ac = None
        self.n = self.ringlwe.n

    def __call__(self):
        """
        EXAMPLES::

            sage: from sage.crypto.lwe import DiscreteGaussianDistributionPolynomialSampler, RingLWE, RingLWEConverter
            sage: D = DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], euler_phi(16), 5)
            sage: lwe = RingLWEConverter(RingLWE(16, 257, D, secret_dist='uniform'))
            sage: set_random_seed(1337)
            sage: lwe()
            ((171, 197, 58, 125, 3, 216, 32, 130), ...)
        """
        R_q = self.ringlwe.R_q

        if (self._i % self.n) == 0:
            self._ac = self.ringlwe()
        a, c = self._ac
        x = R_q.gen()
        r = vector((x**(self._i % self.n) * R_q(a.list())).list()), c[self._i % self.n]
        self._i += 1
        return r

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.crypto.lwe import DiscreteGaussianDistributionPolynomialSampler, RingLWE, RingLWEConverter
            sage: D = DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], euler_phi(20), 5)
            sage: rlwe = RingLWE(20, 257, D)
            sage: lwe = RingLWEConverter(rlwe)
            sage: lwe
            RingLWEConverter(RingLWE(20, 257, Discrete Gaussian sampler for polynomials of degree < 8 with σ=5.000000 in each component, x^8 - x^6 + x^4 - x^2 + 1, 'uniform', None))
        """
        return "RingLWEConverter(%s)" % str(self.ringlwe)

def samples(m, n, lwe, seed=None, balanced=False, **kwds):
    """
    Return ``m`` LWE samples.

    INPUT:

    - ``m`` -- the number of samples (integer > 0)
    - ``n`` -- the security parameter (integer > 0)
    - ``lwe`` -- either

      - a subclass of :class:`LWE` such as :class:`Regev` or :class:`LindnerPeikert`
      - an instance of :class:`LWE` or any subclass
      - the name of any such class (e.g., "Regev", "LindnerPeikert")

    - ``seed`` -- seed to be used for generation or ``None`` if no specific seed
      shall be set (default: ``None``)
    - ``balanced`` -- use function :func:`balance_sample` to return balanced
      representations of finite field elements (default: ``False``)
    - ``**kwds`` -- passed through to LWE constructor

    EXAMPLES::

        sage: from sage.crypto.lwe import samples, Regev
        sage: samples(2, 20, Regev, seed=1337)
        [((199, 388, 337, 53, 200, 284, 336, 215, 75, 14, 274, 234, 97, 255, 246, 153, 268, 218, 396, 351), 15),
         ((365, 227, 333, 165, 76, 328, 288, 206, 286, 42, 175, 155, 190, 275, 114, 280, 45, 218, 304, 386), 143)]

        sage: from sage.crypto.lwe import samples, Regev
        sage: samples(2, 20, Regev, balanced=True, seed=1337)
        [((199, -13, -64, 53, 200, -117, -65, -186, 75, 14, -127, -167, 97, -146, -155, 153, -133, -183, -5, -50), 15),
         ((-36, -174, -68, 165, 76, -73, -113, -195, -115, 42, 175, 155, 190, -126, 114, -121, 45, -183, -97, -15), 143)]

        sage: from sage.crypto.lwe import samples
        sage: samples(2, 20, 'LindnerPeikert')
        [((506, 1205, 398, 0, 337, 106, 836, 75, 1242, 642, 840, 262, 1823, 1798, 1831, 1658, 1084, 915, 1994, 163), 1447),
         ((463, 250, 1226, 1906, 330, 933, 1014, 1061, 1322, 2035, 1849, 285, 1993, 1975, 864, 1341, 41, 1955, 1818, 1357), 312)]
    """
    if seed is not None:
        set_random_seed(seed)

    if isinstance(lwe, str):
        lwe = eval(lwe)

    if isinstance(lwe, type):
        lwe = lwe(n, m=m, **kwds)
    else:
        if lwe.n != n:
            raise ValueError("Passed LWE instance has n=%d, but n=%d was passed to this function." % (lwe.n, n))

    if balanced is False:
        f = lambda a_c: a_c
    else:
        f = balance_sample
    return [f(lwe()) for _ in range(m)]


def balance_sample(s, q=None):
    r"""
    Given ``(a,c) = s`` return a tuple ``(a',c')`` where ``a'`` is an integer
    vector with entries between -q//2 and q//2 and ``c`` is also within these
    bounds.

    If ``q`` is given ``(a,c) = s`` may live in the integers. If ``q`` is not
    given, then ``(a,c)`` are assumed to live in `\Zmod{q}`.

    INPUT:

    - ``s`` -- sample of the form (a,c) where a is a vector and c is a scalar
    - ``q`` -- modulus (default: ``None``)

    EXAMPLES::

        sage: from sage.crypto.lwe import balance_sample, samples, Regev
        sage: for s in samples(10, 5, Regev):
        ....:     b = balance_sample(s)
        ....:     assert all(-29//2 <= c <= 29//2 for c in b[0])
        ....:     assert -29//2 <= b[1] <= 29//2
        ....:     assert all(s[0][j] == b[0][j] % 29 for j in range(5))
        ....:     assert s[1] == b[1] % 29


        sage: from sage.crypto.lwe import balance_sample, DiscreteGaussianDistributionPolynomialSampler, RingLWE, samples
        sage: D = DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], 8, 5)
        sage: rlwe = RingLWE(20, 257, D)
        sage: for s in samples(10, 8, rlwe):
        ....:     b = balance_sample(s)
        ....:     assert all(-257//2 <= c <= 257//2 for bi in b for c in bi)
        ....:     assert all(s[i][j] == b[i][j] % 257 for i in range(2) for j in range(8))

    .. NOTE::

        This function is useful to convert between Sage's standard
        representation of elements in `\Zmod{q}` as integers between 0 and q-1
        and the usual representation of such elements in lattice cryptography as
        integers between -q//2 and q//2.
    """
    a, c = s

    try:
        c[0]
        scalar = False
    except TypeError:
        c = vector(c.parent(),[c])
        scalar = True

    if q is None:
        q = parent(c[0]).order()
        a = a.change_ring(ZZ)
        c = c.change_ring(ZZ)
    else:
        K = IntegerModRing(q)
        a = a.change_ring(K).change_ring(ZZ)
        c = c.change_ring(K).change_ring(ZZ)

    q2 = q//2

    if scalar:
        return vector(ZZ, len(a), [e if e <= q2 else e-q for e in a]), c[0] if c[0] <= q2 else c[0]-q
    else:
        return vector(ZZ, len(a), [e if e <= q2 else e-q for e in a]), vector(ZZ, len(c), [e if e <= q2 else e-q for e in c])
