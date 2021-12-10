"""
`L`-series of modular abelian varieties

AUTHOR:

- William Stein (2007-03)

TESTS::

    sage: L = J0(37)[0].padic_lseries(5)
    sage: loads(dumps(L)) == L
    True
    sage: L = J0(37)[0].lseries()
    sage: loads(dumps(L)) == L
    True
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.integer import Integer
from sage.rings.infinity import infinity
from sage.rings.cc import CC
from sage.modules.free_module import span
from sage.misc.misc_c import prod


class Lseries(SageObject):
    """
    Base class for `L`-series attached to modular abelian varieties.

    This is a common base class for complex and `p`-adic `L`-series
    of modular abelian varieties.
    """
    def __init__(self, abvar):
        """
        Called when creating an L-series.

        INPUT:

        - ``abvar`` -- a modular abelian variety

        EXAMPLES::

            sage: J0(11).lseries()
            Complex L-series attached to Abelian variety J0(11) of dimension 1
            sage: J0(11).padic_lseries(7)
            7-adic L-series attached to Abelian variety J0(11) of dimension 1
        """
        self.__abvar = abvar

    def abelian_variety(self):
        """
        Return the abelian variety that this `L`-series is attached to.

        OUTPUT:

        a modular abelian variety

        EXAMPLES::

            sage: J0(11).padic_lseries(7).abelian_variety()
            Abelian variety J0(11) of dimension 1
        """
        return self.__abvar


class Lseries_complex(Lseries):
    """
    A complex `L`-series attached to a modular abelian variety.

    EXAMPLES::

        sage: A = J0(37)
        sage: A.lseries()
        Complex L-series attached to Abelian variety J0(37) of dimension 2
    """
    def __call__(self, s, prec=53):
        """
        Evaluate this complex `L`-series at `s`.

        INPUT:

        - ``s`` -- complex number

        - ``prec`` -- integer (default: 53) the number of bits of precision
          used in computing the lseries of the newforms.

        OUTPUT:

        a complex number L(A, s).

        EXAMPLES::

            sage: L = J0(23).lseries()
            sage: L(1)
            0.248431866590600
            sage: L(1, prec=100)
            0.24843186659059968120725033931

            sage: L = J0(389)[0].lseries()
            sage: L(1) # long time (2s) abstol 1e-10
            -1.33139759782370e-19
            sage: L(1, prec=100) # long time (2s) abstol 1e-20
            6.0129758648142797032650287762e-39
            sage: L.rational_part()
            0

            sage: L = J1(23)[0].lseries()
            sage: L(1)
            0.248431866590600

            sage: J = J0(11) * J1(11)
            sage: J.lseries()(1)
            0.0644356903227915

            sage: L = JH(17,[2]).lseries()
            sage: L(1)
            0.386769938387780

        """
        abelian_variety = self.abelian_variety()
        # Check for easy dimension zero case
        if abelian_variety.dimension() == 0:
            return CC(1)
        try:
            factors = self.__factors[prec]
            return prod(L(s) for L in factors)
        except AttributeError:
            self.__factors = {}
        except KeyError:
            pass
        abelian_variety = self.abelian_variety()
        newforms = abelian_variety.newform_decomposition('a')

        factors = [newform.lseries(embedding=i, prec=prec)
                for newform in newforms
                for i in range(newform.base_ring().degree())]
        self.__factors[prec] = factors

        return prod(L(s) for L in factors)

    def __eq__(self, other):
        """
        Compare this complex `L`-series to another one.

        INPUT:

        - ``other`` -- object

        OUTPUT:

        boolean

        EXAMPLES::

            sage: L = J0(37)[0].lseries()
            sage: M = J0(37)[1].lseries()
            sage: L == M
            False
            sage: L == L
            True
        """
        if not isinstance(other, Lseries_complex):
            return False
        return self.abelian_variety() == other.abelian_variety()

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        INPUT:

        - ``other`` -- object

        OUTPUT:

        boolean

        EXAMPLES::

            sage: L = J0(37)[0].lseries()
            sage: M = J0(37)[1].lseries()
            sage: L != M
            True
            sage: L != L
            False
        """
        return not (self == other)

    def _repr_(self):
        """
        String representation of `L`-series.

        OUTPUT:

        a string

        EXAMPLES::

            sage: L = J0(37).lseries()
            sage: L._repr_()
            'Complex L-series attached to Abelian variety J0(37) of dimension 2'
        """
        return "Complex L-series attached to %s" % self.abelian_variety()

    def vanishes_at_1(self):
        """
        Return True if `L(1)=0` and return False otherwise.

        OUTPUT:

        a boolean

        EXAMPLES:

        Numerically, the `L`-series for `J_0(389)` appears to vanish
        at 1.  This is confirmed by this algebraic computation::

            sage: L = J0(389)[0].lseries(); L
            Complex L-series attached to Simple abelian subvariety 389a(1,389) of dimension 1 of J0(389)
            sage: L(1) # long time (2s) abstol 1e-10
            -1.33139759782370e-19
            sage: L.vanishes_at_1()
            True

        Numerically, one might guess that the `L`-series for `J_1(23)`
        and `J_1(31)` vanish at 1.  This algebraic computation shows
        otherwise::

            sage: L = J1(23).lseries(); L
            Complex L-series attached to Abelian variety J1(23) of dimension 12
            sage: L(1)  # long time (about 3 s)
            0.0001295198...
            sage: L.vanishes_at_1()
            False
            sage: abs(L(1, prec=100)- 0.00012951986142702571478817757148) < 1e-32  # long time (about 3 s)
            True

            sage: L = J1(31).lseries(); L
            Complex L-series attached to Abelian variety J1(31) of dimension 26
            sage: abs(L(1) - 3.45014267547611e-7) < 1e-15  # long time (about 8 s)
            True
            sage: L.vanishes_at_1()  # long time (about 6 s)
            False
        """
        abelian_variety = self.abelian_variety()
        # Check for easy dimension zero case
        if abelian_variety.dimension() == 0:
            return False
        if not abelian_variety.is_simple():
            from .constructor import AbelianVariety
            decomp = (AbelianVariety(f) for f in
                      abelian_variety.newform_decomposition('a'))
            return any(S.lseries().vanishes_at_1() for S in decomp)
        modular_symbols = abelian_variety.modular_symbols()
        Phi = modular_symbols.rational_period_mapping()
        ambient_module = modular_symbols.ambient_module()

        e = ambient_module([0, infinity])
        return Phi(e).is_zero()

    def rational_part(self):
        """
        Return the rational part of this `L`-function at the central critical
        value 1.

        OUTPUT:

        a rational number

        EXAMPLES::

            sage: A, B = J0(43).decomposition()
            sage: A.lseries().rational_part()
            0
            sage: B.lseries().rational_part()
            2/7
        """
        abelian_variety = self.abelian_variety()
        modular_symbols = abelian_variety.modular_symbols()
        Phi = modular_symbols.rational_period_mapping()
        ambient_module = modular_symbols.ambient_module()

        if self.vanishes_at_1():
            return QQ(0)
        else:
            s = ambient_module.sturm_bound()
            I = ambient_module.hecke_images(0, range(1, s+1))
            PhiTe = span([Phi(ambient_module(I[n]))
                for n in range(I.nrows())], ZZ)

        ambient_plus = ambient_module.sign_submodule(1)
        ambient_plus_cusp = ambient_plus.cuspidal_submodule()
        PhiH1plus = span([Phi(x) for
            x in ambient_plus_cusp.integral_basis()], ZZ)

        return PhiTe.index_in(PhiH1plus)

    lratio = rational_part


class Lseries_padic(Lseries):
    """
    A `p`-adic `L`-series attached to a modular abelian variety.
    """
    def __init__(self, abvar, p):
        """
        Create a `p`-adic `L`-series.

        EXAMPLES::

            sage: J0(37)[0].padic_lseries(389)
            389-adic L-series attached to Simple abelian subvariety 37a(1,37) of dimension 1 of J0(37)
        """
        Lseries.__init__(self, abvar)
        p = Integer(p)
        if not p.is_prime():
            raise ValueError("p (=%s) must be prime"%p)
        self.__p = p

    def __eq__(self, other):
        """
        Compare this `p`-adic `L`-series to another one.

        First the abelian varieties are compared; if they are the same,
        then the primes are compared.

        INPUT:

        other -- object

        OUTPUT:

        boolean

        EXAMPLES::

            sage: L = J0(37)[0].padic_lseries(5)
            sage: M = J0(37)[1].padic_lseries(5)
            sage: K = J0(37)[0].padic_lseries(3)
            sage: L == K
            False
            sage: L == M
            False
            sage: L == L
            True
        """
        if not isinstance(other, Lseries_padic):
            return False
        return (self.abelian_variety() == other.abelian_variety() and
                self.__p == other.__p)

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        INPUT:

        other -- object

        OUTPUT:

        boolean

        EXAMPLES::

            sage: L = J0(37)[0].padic_lseries(5)
            sage: M = J0(37)[1].padic_lseries(5)
            sage: K = J0(37)[0].padic_lseries(3)
            sage: L != K
            True
            sage: L != M
            True
            sage: L != L
            False
        """
        return not (self == other)

    def prime(self):
        """
        Return the prime `p` of this `p`-adic `L`-series.

        EXAMPLES::

            sage: J0(11).padic_lseries(7).prime()
            7
        """
        return self.__p

    def power_series(self, n=2, prec=5):
        """
        Return the `n`-th approximation to this `p`-adic `L`-series as
        a power series in `T`.

        Each coefficient is a `p`-adic number
        whose precision is provably correct.

        NOTE: This is not yet implemented.

        EXAMPLES::

            sage: L = J0(37)[0].padic_lseries(5)
            sage: L.power_series()
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: L.power_series(3,7)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _repr_(self):
        """
        String representation of this `p`-adic `L`-series.

        EXAMPLES::

            sage: L = J0(37)[0].padic_lseries(5)
            sage: L._repr_()
            '5-adic L-series attached to Simple abelian subvariety 37a(1,37) of dimension 1 of J0(37)'
        """
        return "%s-adic L-series attached to %s" % (self.__p,
                                                    self.abelian_variety())
