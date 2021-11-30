"""
Markov Switching Multifractal model

REFERENCE:

*How to Forecast Long-Run Volatility: Regime Switching and
the Estimation of Multifractal Processes*, Calvet and Fisher, 2004.

AUTHOR:

- William Stein, 2008

TESTS::

    sage: msm = finance.MarkovSwitchingMultifractal(8,1.4,1.0,0.95,3)
    doctest:warning...
    DeprecationWarning: the package sage.finance is deprecated...
    sage: loads(dumps(msm)) == msm
    True
"""
import math

class MarkovSwitchingMultifractal:
    def __init__(self, kbar, m0, sigma, gamma_kbar, b):
        """
        INPUT:

        - ``kbar`` -- positive integer

        - ``m0`` -- float with ``0 <= m0 <= 2``

        - ``sigma`` -- positive float

        - ``gamma_kbar`` -- float with ``0 <= gamma_kbar < 1``

        - ``b`` -- float > 1

        EXAMPLES::

            sage: msm = finance.MarkovSwitchingMultifractal(8,1.4,0.5,0.95,3); msm
            Markov switching multifractal model with m0 = 1.4, sigma = 0.5, b = 3.0, and gamma_8 = 0.95
            sage: yen_usd = finance.MarkovSwitchingMultifractal(10,1.448,0.461,0.998,3.76)
            sage: cad_usd = finance.MarkovSwitchingMultifractal(10,1.278,0.262,0.644,2.11)
            sage: dm = finance.MarkovSwitchingMultifractal(10,1.326,0.643,0.959,2.7)
        """
        self.__m0 = float(m0)
        assert self.__m0 >= 0 and self.__m0 <= 2, "m0 must be between 0 and 2"
        self.__sigma = float(sigma)
        assert self.__sigma > 0, "sigma must be positive"
        self.__b = float(b)
        assert self.__b > 1, "b must be bigger than 1"
        self.__gamma_kbar = float(gamma_kbar)
        assert self.__gamma_kbar >= 0 and self.__gamma_kbar < 1, \
               "gamma_kbar must be between 0 and 1"
        self.__kbar = int(kbar)
        assert self.__kbar > 0, "kbar must be positive"

    def __eq__(self, other):
        """
        Test equality of ``self`` and ``other``.

        Comparison is done on the tuple ``(m0, sigma, b, gamma_kbar, kbar)``.

        EXAMPLES::

            sage: msm = finance.MarkovSwitchingMultifractal(8,1.4,1.0,0.95,3)

            sage: msm == msm
            True
            sage: cad_usd = finance.MarkovSwitchingMultifractal(10,1.278,0.262,0.644,2.11); cad_usd
            Markov switching multifractal model with m0 = 1.278, sigma = 0.262, b = 2.11, and gamma_10 = 0.644
            sage: msm == cad_usd
            False
        """
        if not isinstance(other, MarkovSwitchingMultifractal):
            return False
        return (self.__m0 == other.__m0 and
                self.__sigma == other.__sigma and
                self.__b == other.__b and
                self.__gamma_kbar == other.__gamma_kbar and
                self.__kbar == other.__kbar)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: msm = finance.MarkovSwitchingMultifractal(8,1.4,1.0,0.95,3)
            sage: H = hash(msm)
        """
        return hash((self.__m0, self.__sigma, self.__b, self.__gamma_kbar,
                     self.__kbar))

    def __ne__(self, other):
        """
        Test inequality of ``self`` and ``other``.

        EXAMPLES::

            sage: msm = finance.MarkovSwitchingMultifractal(8,1.4,1.0,0.95,3)

            sage: msm != msm
            False
            sage: cad_usd = finance.MarkovSwitchingMultifractal(10,1.278,0.262,0.644,2.11); cad_usd
            Markov switching multifractal model with m0 = 1.278, sigma = 0.262, b = 2.11, and gamma_10 = 0.644
            sage: msm != cad_usd
            True
        """
        return not (self == other)

    def __repr__(self):
        """
        Return string representation of Markov switching multifractal model.

        EXAMPLES::

            sage: msm = finance.MarkovSwitchingMultifractal(8,1.4,1,0.95,3)
            sage: msm.__repr__()
            'Markov switching multifractal model with m0 = 1.4, sigma = 1.0, b = 3.0, and gamma_8 = 0.95'
        """
        return "Markov switching multifractal model with m0 = %s, sigma = %s, b = %s, and gamma_%s = %s"%(self.m0(), self.sigma(), self.b(), self.kbar(), self.gamma_kbar())

    def m0(self):
        """
        Return parameter m0 of Markov switching multifractal model.

        EXAMPLES::

            sage: msm = finance.MarkovSwitchingMultifractal(8,1.4,1,0.95,3)
            sage: msm.m0()
            1.4
        """
        return self.__m0

    def sigma(self):
        """
        Return parameter sigma of Markov switching multifractal model.

        EXAMPLES::

            sage: msm = finance.MarkovSwitchingMultifractal(8,1.4,1,0.95,3)
            sage: msm.sigma()
            1.0
        """
        return self.__sigma

    def b(self):
        """
        Return parameter b of Markov switching multifractal model.

        EXAMPLES::

            sage: msm = finance.MarkovSwitchingMultifractal(8,1.4,1,0.95,3)
            sage: msm.b()
            3.0
        """
        return self.__b

    def gamma_kbar(self):
        """
        Return parameter ``gamma_kbar`` of Markov switching multifractal model.

        EXAMPLES::

            sage: msm = finance.MarkovSwitchingMultifractal(8,1.4,0.01,0.95,3)
            sage: msm.gamma_kbar()
            0.95
        """
        return self.__gamma_kbar

    def kbar(self):
        """
        Return parameter ``kbar`` of Markov switching multifractal model.

        EXAMPLES::

            sage: msm = finance.MarkovSwitchingMultifractal(8,1.4,0.01,0.95,3)
            sage: msm.kbar()
            8
        """
        return self.__kbar

    def gamma(self):
        """
        Return the vector of the kbar transitional probabilities.

        OUTPUT:

        - gamma -- a tuple of ``self.kbar()`` floats.

        EXAMPLES::

            sage: msm = finance.MarkovSwitchingMultifractal(8,1.4,1.0,0.95,3)
            sage: msm.gamma()
            (0.001368852970712986, 0.004100940201672509, 0.012252436441829..., 0.03630878209190..., 0.10501923017634..., 0.28312883556311..., 0.6315968501359..., 0.95000000000000...)
        """
        try:
            return self.__gamma
        except AttributeError:
            pass

        b          = self.__b
        gamma_kbar = self.__gamma_kbar
        kbar       = self.__kbar

        # We compute gamma1 from gamma_kbar by inverting the relation
        # that defines the gamma_k given on page 54 of Calvet-Fisher:
        gamma1 = 1 - math.exp(math.log(1-gamma_kbar)/(b**(kbar-1)))

        gamma  = tuple([1 - (1 - gamma1)**(b**k) for k in range(kbar)])
        self.__gamma = gamma
        return gamma

    def simulation(self, n):
        """
        Same as ``self.simulations``, but run only 1 time, and returns a time
        series instead of a list of time series.

        INPUT:

        - ``n`` -- a positive integer.

        EXAMPLES::

            sage: msm = finance.MarkovSwitchingMultifractal(8,1.4,1.0,0.95,3)
            sage: m = msm.simulation(5); m  # random
            [0.0059, -0.0097, -0.0101, -0.0110, -0.0067]
            sage: len(m)
            5
            sage: m = msm.simulation(3); m  # random
            [0.0055, -0.0084, 0.0141]
            sage: len(m)
            3
        """
        return self.simulations(n, 1)[0]

    def simulations(self, n, k=1):
        """
        Return ``k`` simulations of length ``n`` using this Markov switching
        multifractal model for ``n`` time steps.

        INPUT:

        - ``n`` -- positive integer; number of steps.

        - ``k`` -- positive integer (default: 1); number of simulations.

        OUTPUT:

        list -- a list of TimeSeries objects.

        EXAMPLES::

            sage: cad_usd = finance.MarkovSwitchingMultifractal(10,1.278,0.262,0.644,2.11); cad_usd
            Markov switching multifractal model with m0 = 1.278, sigma = 0.262, b = 2.11, and gamma_10 = 0.644
        """
        from . import markov_multifractal_cython
        return markov_multifractal_cython.simulations(n, k,
                   self.__m0, self.__sigma,
                   self.__kbar, self.gamma())



## def ml_estimation(v, kbar, M):
##     """
##     Compute parameters that model the time series v,

##     INPUT:
##         v -- series of returns; e.g., sequence of
##              differences of logs of price
##         kbar -- positive integer; model parameter
##         m -- finite list of the values that the multiplier
##              M takes on.

##     OUTPUT:
##         m0, sigma, gamma_kbar, b
##     """

